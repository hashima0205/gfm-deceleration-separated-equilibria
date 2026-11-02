# Finite-horizon criteria shared by trajectory-neighborhood and voltage probes.
function interp010(t,a,q)
    j=clamp(searchsortedlast(t,q),1,length(t)-1)
    a[j]+(a[j+1]-a[j])*(q-t[j])/(t[j+1]-t[j])
end

function classify_probe(t,u,target; phase_tolerance=fill(.002,4))
    finish=t[end];tail=findall(t .>= finish-min(50.,finish/3))
    machines=[]
    for j in eachindex(ALLG)
        angle=[x[DI010[j]] for x in u]
        freq=[KP010[j]*(PS010[j]-x[PI010[j]])/(2pi) for x in u]
        shift=angle .- angle[1];turns=shift[end]/(2pi)
        error=abs(atan(sin(angle[end]-target[DI010[j]]),cos(angle[end]-target[DI010[j]])))
        width=maximum(angle[tail])-minimum(angle[tail])
        windows=diff([interp010(t,angle,q) for q in range(.4*finish,finish;length=4)])/(2pi)
        dfmax=maximum(abs.(freq[tail]))
        label=if dfmax<1e-4 && error<phase_tolerance[j] && width<.005
            maximum(abs.(shift))>=2pi ? "recovery_after_slip" : "recovery"
        elseif turns< -2 && maximum(windows)< -.2 && sum(windows)/(.6*finish)< -1e-4
            "negative_running"
        elseif turns>2 && minimum(windows)>.2 && sum(windows)/(.6*finish)>1e-4
            "positive_running"
        else
            "unresolved"
        end
        push!(machines,(ALLG[j],label,turns,dfmax,error,width,minimum(windows),maximum(windows)))
    end
    labels=[m[2] for m in machines]
    label=if all(l in ["recovery","recovery_after_slip"] for l in labels)
        "recovery_after_slip" in labels ? "recovery_after_slip" : "recovery"
    elseif "negative_running" in labels && "positive_running" in labels
        "other"
    elseif "negative_running" in labels
        "negative_running"
    elseif "positive_running" in labels
        "positive_running"
    else
        "unresolved"
    end
    turns=sum(u[end][DELTA]-u[1][DELTA])/(6pi)
    label,turns,machines
end
