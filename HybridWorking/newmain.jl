using NLPModels, LinearAlgebra, CUTEst, Printf
include("hybrid_new.jl")
include("linesearch.jl")

output = open("info.dat","w")

#custom_filter = x->x["origin"]=="real" 
#problems = ["MNISTS5LS", "MGH10SLS", "ARGLINA", "SENSORS", "RAT43LS", "CLUSTERLS", "HELIX", "MEXHAT", "DENSCHNB", "MNISTS0LS", "MGH09LS", "MEYER3", "GULF", "DMN15333LS", "DENSCHND", "BROWNBS", "SSI", "HATFLDE", "LUKSAN11LS", "KIRBY2LS", "LUKSAN16LS", "COATING", "ERRINRSM", "BENNETT5LS", "DANWOODLS", "ARGTRIGLS", "DANIWOODLS", "RAT42LS", "DEVGLA1", "HATFLDFL", "STREG", "AKIVA", "QING", "LANCZOS2LS", "ARGLINB", "PARKCH", "ROSENBR", "HYDC20LS", "HATFLDFLS", "ERRINROS", "LUKSAN14LS", "HIELOW", "BROWNAL", "HUMPS", "BARD", "HATFLDGLS", "MGH10LS", "DEVGLA2", "CHWIRUT1LS", "LUKSAN15LS", "EGGCRATE", "HILBERTB", "BA-L1LS", "RECIPELS", "CHNROSNB", "TOINTPSP", "LANCZOS1LS", "PALMER3C", "ELATVIDU", "DMN15103LS", "TOINTQOR", "ROSZMAN1LS", "POWERSUM", "ROSENBRTU", "WATSON", "GAUSSIAN", "GENROSE", "PALMER4C", "MISRA1CLS", "DENSCHNE", "OSBORNEB", "CLIFF", "STRTCHDV", "TRIGON2", "HEART6LS", "POWELLSQLS", "EXP2", "METHANB8LS", "LANCZOS3LS", "DIAMON3DLS", "LOGHAIRY", "CERI651BLS", "CUBE"]
problems = CUTEst.select(contype="unc",min_var=2,max_var=100)
#problems = CUTEst.select(objtype="sum_of_squares", contype="unc", custom_filter=custom_filter)
#problems = CUTEst.select(objtype="linear", contype="unc", custom_filter=custom_filter)
#problems = CUTEst.select(objtype="quadratic", contype="unc", custom_filter=custom_filter)

length(problems)

DELTA_tests = [2.0, 1.5, 1.0, 1.e-2, 1.e-4, 1.e-5]

epsilon = 1.e-6
maxiter = 50000
gamma = 1.e-3
stpmin = 1.e-12

#LG = armijo
#LN = armijo
LineSearchs = [armijo, wolfe]

for p in problems
    nlp = CUTEstModel(p)

    global function f(x)
    return obj(nlp,x)
    end

    global function g(x)
    return grad(nlp,x)
    end

    global function h(x)
    return hess(nlp,x)
    end

    x_0 = nlp.meta.x0
    for LG in LineSearchs
        for LN in LineSearchs

            for delta in DELTA_tests
                x,error,iter,N_count,norm_gradf_x,et = hybrid(x_0,f,g,h,epsilon,maxiter,gamma,delta,LG,LN,stpmin)
                @printf(       "%15s %10.1e  %8d %5d  %20.15e %10s %10s %15.10e %3d\n",p,delta,iter,N_count,norm_gradf_x,string(LG),string(LN),et,error)
                @printf(output,"%15s %10.1e  %8d %5d  %20.15e %10s %10s %15.10e %3d\n",p,delta,iter,N_count,norm_gradf_x,string(LG),string(LN),et,error)
            end
        end
    end
    finalize(nlp)
end


close(output)