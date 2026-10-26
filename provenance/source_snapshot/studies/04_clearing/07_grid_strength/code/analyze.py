"""Reconstruct actual network powers, external-grid SCR proxies, and classifications."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent.parent
OUT=HERE/'output'
BUSES=[1,2,3,7]
LINES=[(1,4,2.14),(4,6,3.43),(5,6,.68),(3,5,2.21),(2,3,1.91),(2,5,2.14),(4,5,1.8),(2,8,2.17),(2,7,6.39),(7,8,1.69),(5,8,8.96)]
def impedances(r):
    return [(r.k/m*np.exp(1j*np.deg2rad(r.angle)) if r['mode']=='all' or j in [2,3]
             else 1/m*np.exp(1j*np.deg2rad(75))) for j,(_,_,m) in enumerate(LINES,1)]
def admittance(r,opened=False):
    y=np.zeros((8,8),complex)
    for j,((a,b,_),z) in enumerate(zip(LINES,impedances(r)),1):
        if j==7 and opened:continue
        a-=1;b-=1;v=1/z;y[a,a]+=v;y[b,b]+=v;y[a,b]-=v;y[b,a]-=v
    return y
def scr(r):
    ans={};free=[0,1,2,3,4,6,7]
    for opened in [False,True]:
        Z=np.linalg.inv(admittance(r,opened)[np.ix_(free,free)])
        for b in BUSES:
            z=Z[free.index(b-1),free.index(b-1)];stage='post' if opened else 'pre'
            ans[f'{stage}_scr{b}']=1/abs(z)
            ans[f'{stage}_zth_angle{b}']=np.rad2deg(np.angle(z))
    return ans
def read(r):
    d=pd.read_csv(OUT/'runs'/(r['case']+'.csv')).drop_duplicates('t',keep='last')
    t=d.t.to_numpy();V=np.column_stack([d[f'net_b{b}_u_r']+1j*d[f'net_b{b}_u_i'] for b in range(1,9)])
    into=np.column_stack([d[f'net_b{b}_i_r']+1j*d[f'net_b{b}_i_i'] for b in range(1,9)])
    net=np.zeros_like(V);loss=np.zeros(len(t));clear=1+r.duration
    for j,((a,b,_),zz) in enumerate(zip(LINES,impedances(r)),1):
        active=np.ones(len(t)) if j!=7 else np.where(t>=clear,0.,1.)
        y=active/zz;i=y*(V[:,a-1]-V[:,b-1]);net[:,a-1]+=i;net[:,b-1]-=i
        loss+=y.real*abs(V[:,a-1]-V[:,b-1])**2
    z=pd.DataFrame({'t':t,'V4':abs(V[:,3]),'V5':abs(V[:,4]),'Pslack':(V[:,5]*np.conj(net[:,5])).real,'line_loss':loss})
    for b in BUSES:
        for key in ['δ','Pfilt']:
            z[f'{"delta" if key=="δ" else key}{b}']=d[f'b{b}_droop₊{key}'].to_numpy()
        z[f'df{b}']=d[f'b{b}_droop₊ω'].to_numpy()/(2*np.pi)-60
        z[f'phase{b}']=np.unwrap(np.angle(V[:,b-1]))
        z[f'V{b}']=abs(V[:,b-1]);z[f'Pterminal{b}']=(V[:,b-1]*np.conj(net[:,b-1])).real
        z[f'Pmeas{b}']=(d[f'b{b}_vsrc₊filter₊V_C_r']*d[f'b{b}_vsrc₊filter₊i_g_r']+d[f'b{b}_vsrc₊filter₊V_C_i']*d[f'b{b}_vsrc₊filter₊i_g_i']).to_numpy()
        z[f'I{b}']=np.hypot(d[f'b{b}_vsrc₊filter₊i_f_r'],d[f'b{b}_vsrc₊filter₊i_f_i']).to_numpy()
        z[f'Iref{b}']=np.hypot(d[f'b{b}_vsrc₊cc1₊i_f_ref_d'],d[f'b{b}_vsrc₊cc1₊i_f_ref_q']).to_numpy()
    return z,float(np.max(abs(net+into)))
def main():
    rows=[]
    validation=pd.read_csv(OUT/'equilibrium_validation.csv') if (OUT/'equilibrium_validation.csv').exists() else pd.DataFrame()
    certificates=pd.read_csv(OUT/'static_checks.csv') if (OUT/'static_checks.csv').exists() else pd.DataFrame()
    for p in sorted(OUT.glob('*_status.csv')):
        for _,r in pd.read_csv(p).iterrows():
            a=r.to_dict();a.update(scr(r))
            a['post_equilibrium']='stable' if r.post_residual<1e-8 and r.post_eigen<0 and r.post_current<=1.200001 else ('unstable' if r.post_residual<1e-8 else 'unresolved')
            if len(validation):
                vr=validation[(validation['mode']==r['mode'])&(validation.k==r.k)&(validation.angle==r.angle)&(validation.tau==r.tau)&(validation.residual<1e-8)&(validation.max_real_eigenvalue<0)&(validation.max_current<=1.200001)]
                if len(vr):
                    q=vr.iloc[0];a.update(post_equilibrium='stable',validated_post_residual=q.residual,validated_post_eigen=q.max_real_eigenvalue,post_delta2=q.delta2,post_V2=q.V2)
            if len(certificates):
                cr=certificates[(certificates['mode']==r['mode'])&(certificates.k==r.k)&(certificates.angle==r.angle)]
                if len(cr) and cr.iloc[0].nonexistence_certified:a['post_equilibrium']='nonexistent_certified'
            a['classification']=r.status
            if r.status!='Success':rows.append(a);continue
            z,kcl=read(r);a['max_kcl']=kcl
            tail=z.t>=r.horizon-10;post=z.t>=1+r.duration+.02;fault=(z.t>=1)&(z.t<1+r.duration)
            neg=[];pos=[];rec=[]
            def avg(v,mask):
                return float(np.trapezoid(v[mask],z.t[mask])/(z.t[mask].iloc[-1]-z.t[mask].iloc[0]))
            for b in BUSES:
                turns=(z[f'phase{b}'].iloc[-1]-z[f'phase{b}'].iloc[0])/(2*np.pi)
                iturns=(z[f'delta{b}'].iloc[-1]-z[f'delta{b}'].iloc[0])/(2*np.pi)
                df=avg(z[f'df{b}'],tail)
                a.update({f'turns{b}':turns,f'internal_turns{b}':iturns,f'tail_df{b}':df,f'tail_max_df{b}':z.loc[tail,f'df{b}'].abs().max(),f'fault_min_V{b}':z.loc[fault,f'V{b}'].min(),f'fault_P{b}':avg(z[f'Pmeas{b}'],fault),f'fault_limit_fraction{b}':avg((z[f'Iref{b}']>=1.1999).astype(float),fault),f'post_limit_fraction{b}':avg((z[f'Iref{b}']>=1.1999).astype(float),post)})
                if turns<-1 and iturns<-1 and df<-.001:neg.append(b)
                if turns>1 and iturns>1 and df>.001:pos.append(b)
                if z.loc[tail,f'df{b}'].abs().max()<1e-5 and np.ptp(z.loc[tail,f'phase{b}'])<.001:rec.append(b)
            a['classification']='negative_slip' if neg else ('positive_slip' if pos else ('recovered' if len(rec)==4 else 'unresolved_tail'))
            if kcl>1e-5:a['classification']='audit_failed'
            a.update(negative_buses=';'.join(map(str,neg)),positive_buses=';'.join(map(str,pos)),recovered_buses=';'.join(map(str,rec)),post_min_V5=z.loc[post,'V5'].min(),fault_min_V4=z.loc[fault,'V4'].min(),tail_Pslack=avg(z.Pslack,tail),tail_loss=avg(z.line_loss,tail),max_actual_current=max(z[f'I{b}'].max() for b in BUSES),max_voltage=max(z[f'V{b}'].max() for b in BUSES))
            rows.append(a)
    m=pd.DataFrame(rows);m.to_csv(OUT/'metrics.csv',index=False)
    print(m[['case','mode','k','angle','tau','duration','post_scr2','post_equilibrium','classification','turns2','tail_df2']].to_string(index=False))
    return m
if __name__=='__main__':main()
