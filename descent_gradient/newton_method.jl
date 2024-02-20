using ForwardDiff
using LinearAlgebra

#f(x) = (x[1] - 1)^2 + 4 * (x[2] - x[1]^2)^2

#function newtonmethod(x, f)
    #G(x) = ForwardDiff.gradient(f, x)
    #H(x) = ForwardDiff.hessian(f, x)

    #for k = 1:10
        #g = G(x)
        #h = H(x)
        #d = -h \ g
        #x += d
    #end
    #return x
#end    

#resultado = newtonmethod([-1.2; 1.0], f)
#println(resultado)




#function newtonmethod(x, f; t=1.e-9, max_iter=10)
    #for iter = 1:max_iter
       # g = ForwardDiff.gradient(f, x)
       # h = ForwardDiff.hessian(f, x)
       # d = -h \ g
       # x_k = x + d
       # delta_f = abs(f(x_k) - f(x))

        # if delta_f < t
           # println("Foram realizadas $iter iterações.")
           # return x_k
        # end
        #x = x_k
    #end
    
    #println("O número máximo de iterações foi atingido.")
    #return x
#end    

#f(x) = (x[1] - 1)^2 + 4 * (x[2] - x[1]^2)^2
#resultado = newtonmethod([-1.2; 1.0], f)
#println(resultado)






function newtonmethod(x, f; stpmin=1.e-9, maxiter=15)
   for iter = 1:maxiter
        error=0
        g = ForwardDiff.gradient(f, x)
        h = ForwardDiff.hessian(f, x)
        d = -h \ g
        x_k = x + d
        grad_norm = norm(g)
        
        # Convergence checking
        if grad_norm < stpmin
            return(x_k,grad_norm,iter,error)
        end
        iter= iter+1
        x = x_k
        # Maximum iteration numbers checking
        if iter > maxiter
            error = 1
            return(x_k,grad_norm,iter,error)
        end
    end

end    

include("testfunctions.jl")
newtonmethod([1.0;-1.0], f)




