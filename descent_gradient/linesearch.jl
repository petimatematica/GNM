#
# Linesearch library
#

#
#  Armijo
#
# f(x^k+alpha f'(x^k))<=f(x^k)-sigma alpha||f'(x^k)||^2
#
function armijo(x,f,gfx,sigma,stpmin)
    error = 0
    stp = 1.0
    fx = f(x)
    gtg = sigma * gfx' * gfx
    while true
        q = x - stp * gfx
        fq = f(q)
        stptest = fq - fx + stp * gtg
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
#
function goldstein(x,f,gfx,alpha,beta,stpmin)   
    error = 0
    stp = 1.0
    fx = f(x)
    while true
        gtg = gfx' * gfx
        p = x - stp * gfx 
        fp = f(p) 
        steptestone = fp - fx - alpha * gtg * stp #armijo
        steptestwo = fp - fx - beta * gtg * stp 
        if steptestone > 0.0 || steptestwo < 0
            stp = stp / 2.0  
         else 
            if  steptestwo < 0
                stp = stp * alpha
                if stp<stpmin
                   error=1
                   return(stp,p,error) 
                end
            end
            if steptestone > 0.0
                stp = stp / beta
                if stp<stpmin
                    error=1
                    return(stp,p,error) 
                 end 
            end
            return(stp,p,error)
        end    
    end
end
