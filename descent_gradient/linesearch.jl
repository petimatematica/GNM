#
#   Linesearch library
#

#
#  Armijo
#
function armijo(x,f,gfx,sigma,stpmin)
    error = 0
    stp = 1.0
    fx = f(x)
    gtg = sigma * gfx' * gfx
    while true
        q = x - stp * gfx
        fq = f(q)
        stptest = fq - fx - stp * gtg
        if stptest > 0.0
            stp = stp / 2.0
            if stp < stpmin
                return(stp,q,error)
            end
        else
            return(stp,q,error)
        end
    end
end

#
# Goldstein
# 
function goldstein(x,f,gfx,sigma,alpha,beta,stpmin)   
    error = 0
    stp = 1.0
    fx = f(x)
    while true
        gtg = sigma * gfx' * gfx
        p = x + stp * gfx 
        fp = f(p) 
        steptestone = fp - fx - alpha * sigma * gtg
        steptestwo = fp -fx - beta * sigma * gtg
        if steptestone < 0.0 < steptestwo
            stp = stp / 2.0
        else 
            return (stp, p, error) 
        end    

    end
end