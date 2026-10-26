"""Classify individual machines and their remote group from saved full-model samples."""
from pathlib import Path
import json
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1];D=P/'data';C=P/'checkpoints';B=[1,2,3,7];R=[2,3,7]
S=P/'provenance/source_snapshot/studies/04_clearing/08_state_space/output'
states=pd.read_csv(S/'modal_analysis/states.csv')
eq=np.loadtxt(S/'parameter_fine/condition_007_equilibria.csv',delimiter=',')[:,0]
def targets(vec):return np.array([vec[int(states.loc[states.symbol==f'inv{b}₊droop₊δ','index'].iloc[0])-1] for b in B])
target=targets(eq)
def load(tag):
    fs=[f for f in (C/tag).glob('*.csv') if f.stem in ['pre','fault'] or (f.stem.startswith('post') and f.stem[4:].isdigit())]
    if not fs:
        # Delivery archives can omit duplicate raw observations when the complete
        # enriched CSV retains every original column, saved point and stage.
        full=D/(tag+'.csv.gz')
        return pd.read_csv(full) if full.exists() else None
    d=pd.concat([pd.read_csv(f) for f in sorted(fs)],ignore_index=True).sort_values('t_absolute_s',kind='stable')
    # Exact event left/right samples retained, ordered by stage instead of arbitrary filename.
    d['_order']=d.stage.map(lambda s:0 if s=='pre' else 1 if s=='fault' else 2)
    d=d.sort_values(['t_absolute_s','_order'],kind='stable').drop(columns='_order').reset_index(drop=True)
    for b in B:
        d=d.copy()
        u=d[f'b{b}_U_r_pu']+1j*d[f'b{b}_U_i_pu'];valid=abs(u)>=.1
        d[f'b{b}_terminal_phase_unwrapped_rad']=np.unwrap(np.angle(u))
        d[f'b{b}_terminal_phase_reliable']=valid
        d[f'b{b}_phase_segment']=((valid!=valid.shift(fill_value=False))&valid).cumsum().where(valid,0)
        d[f'b{b}_df_Hz']=d[f'b{b}_omega_rad_s']/(2*np.pi)-60
        d[f'b{b}_V_pu']=abs(u)
    d['theta_mean_rad']=d[[f'b{b}_delta_rad' for b in R]].mean(axis=1)
    d['df_mean_Hz']=d[[f'b{b}_df_Hz' for b in R]].mean(axis=1)
    return d
def measure(tag,d,eqtarget=target):
    post=d[d.stage.str.startswith('post')].drop_duplicates('t_post_s',keep='last')
    end=post.t_post_s.iloc[-1];tail=post[post.t_post_s>=max(0,end-10)]
    machines=[]
    for j,b in enumerate(B):
        delta=post[f'b{b}_delta_rad'];df=tail[f'b{b}_df_Hz'];shift=delta-delta.iloc[0]
        er=abs(np.angle(np.exp(1j*(delta.iloc[-1]-eqtarget[j]))))
        phase=post[f'b{b}_terminal_phase_unwrapped_rad'];valid=post[f'b{b}_terminal_phase_reliable'].all()
        turns=shift.iloc[-1]/(2*np.pi)
        if abs(df).max()<1e-4 and er<.002:
            classification='transient_negative_recovery' if shift.min()<-.01 else 'synchronized_recovery'
        elif turns< -2 and df.max()<-.001:classification='sustained_negative_rotation'
        elif turns>2 and df.min()>.001:classification='positive_rotation'
        else:classification='unresolved'
        sat=d[f'b{b}_saturated'];time=d.t_absolute_s.to_numpy()
        inputmag=np.hypot(d[f'b{b}_ref_in_d_pu'],d[f'b{b}_ref_in_q_pu'])
        outputmag=np.hypot(d[f'b{b}_ref_out_d_pu'],d[f'b{b}_ref_out_q_pu'])
        machines.append(dict(case=tag,bus=b,horizon_s=end,classification=classification,internal_turns=turns,min_shift_rad=shift.min(),phase_error_rad=er,tail_df_mean_Hz=df.mean(),tail_df_absmax_Hz=abs(df).max(),terminal_turns=(phase.iloc[-1]-phase.iloc[0])/(2*np.pi) if valid else np.nan,terminal_phase_all_reliable=valid,voltage_min_pu=d[f'b{b}_V_pu'].min(),actual_current_max_pu=np.hypot(d[f'b{b}_If_r_pu'],d[f'b{b}_If_i_pu']).max(),limiter_input_max_pu=inputmag.max(),limiter_output_max_pu=outputmag.max(),saturation_sampled_duration_s=np.trapezoid(sat,time),tail_voltage_min_pu=tail[f'b{b}_V_pu'].min(),tail_voltage_max_pu=tail[f'b{b}_V_pu'].max()))
    classes=[m['classification'] for m in machines]
    if all('recovery' in s for s in classes):label='transient_negative_recovery' if 'transient_negative_recovery' in classes else 'synchronized_recovery'
    elif 'sustained_negative_rotation' in classes and 'positive_rotation' in classes:label='other_instability'
    elif 'sustained_negative_rotation' in classes:label='sustained_negative_rotation'
    elif 'positive_rotation' in classes:label='positive_rotation'
    else:label='unresolved'
    remote=post[[f'b{b}_delta_rad' for b in R]].to_numpy();freq=post[[f'b{b}_df_Hz' for b in R]].to_numpy()
    out=dict(case=tag,horizon_s=end,classification=label,remote_internal_turns=(post.theta_mean_rad.iloc[-1]-post.theta_mean_rad.iloc[0])/(2*np.pi),tail_remote_df_Hz=tail.df_mean_Hz.mean(),max_remote_phase_spread_rad=np.ptp(remote,axis=1).max(),max_remote_frequency_spread_Hz=np.ptp(freq,axis=1).max(),max_algebraic_residual=d.alg_residual.max(),max_phase_error_rad=max(m['phase_error_rad'] for m in machines),max_tail_df_Hz=max(m['tail_df_absmax_Hz'] for m in machines))
    return out,machines
def main():
    result=[];allmachines=[]
    cases=pd.read_csv(D/'cases.csv').drop_duplicates('case',keep='last')
    for _,r in cases.iterrows():
        d=load(r['case'])
        if d is None or r.retcode!='Success':result.append(dict(**r.to_dict(),classification='numerical_failure'));continue
        eqtarget=target
        if r.event in ['none','fault_only']:
            pre=np.loadtxt(C/r['case']/'pre_initial.csv',delimiter=',');eqtarget=targets(pre)
        z,ms=measure(r['case'],d,eqtarget);z.update({k:r[k] for k in ['model','tau_s','T_s','event','level','solver','retcode']});z['model']='line' if r['case'].startswith('line_') else 'bus';result.append(z);allmachines+=ms
        d.to_csv(D/(r['case']+'.csv.gz'),index=False,compression={'method':'gzip','mtime':0})
    pd.DataFrame(result).to_csv(D/'classification.csv',index=False);pd.DataFrame(allmachines).to_csv(D/'machine_metrics.csv',index=False)
    print(pd.DataFrame(result)[['case','horizon_s','classification','max_phase_error_rad','max_tail_df_Hz']].to_string(index=False))
    if (D/'basin_runs.csv').exists():
        out=[];bs=[]
        for _,r in pd.read_csv(D/'basin_runs.csv').drop_duplicates('case',keep='last').iterrows():
            d=load(r['case'])
            if d is None or r.retcode!='Success':out.append(dict(**r.to_dict(),classification='initialization_or_solver_failure'));continue
            z,ms=measure(r['case'],d);z.update({k:r[k] for k in ['theta_rad','df_Hz','kind','retcode']});z['winding_branch']=round((d.theta_mean_rad.iloc[-1]-target[1:].mean())/(2*np.pi));out.append(z);bs+=ms
            d.to_csv(D/(r['case']+'.csv.gz'),index=False,compression={'method':'gzip','mtime':0})
        pd.DataFrame(out).to_csv(D/'basin_classification.csv',index=False);pd.DataFrame(bs).to_csv(D/'basin_machine_metrics.csv',index=False)
        print('Basin:',pd.DataFrame(out).classification.value_counts().to_dict())
if __name__=='__main__':main()
