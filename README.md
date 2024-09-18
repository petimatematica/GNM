# Hybrid Globalization of Newton's Method


This repository contains the implementation of the Hybrid Method, which combines the global characteristic of the Gradient Method with the fast convergence rates of the Newton Method. By switching between the descent directions of these two methods, the Hybrid Method is used as a globalization strategy for Newton's method.

# hybrid.jl

This document contains the "hybrid" function which has the following inputs: 
- x (vector):The initial point.
- fun (function): Objective function to be minimized.
- grad (function): Gradient of the objective function.
- hess (function): Hessian of the objective function.
- epsilon (float64): Convergence tolerance. 
- maxiter (int): The maximum number of iterations allowed.
- delta (float64): Parameter responsible for changing directions
- stpmin (float64): The minimum step length for line search.
- linesearch (function): The line search function.

# linesearch.jl
For each direction of descent, we employ a linear search that aims to guarantee the decrease of the objective function. The Hybrid Method is evaluated considering the three monotonic searches of Armijo, Goldtein and Wolfe.  These will be tested both in the initial phase led by the MG and in the acceleration phase promoted by the MN. The objective is to investigate the performance of these combinations equipped in the method when determining the global minimum of the objective function.
## Armijo
Armijo's Rule consists of finding a step length that ensures a sufficient reduction in $f$, that is, determining a value $\alpha^k > 0$ that satisfies, for $\sigma\in(0,1)$ the condition


$f(x^k + \alpha^k d^k) \le f(x^k) + \sigma\alpha^k\langle \nabla f(x^k), d^k \rangle.$
## Goldstein 
In Goldstein's Rule, an inequality is added to Armijo's condition in order to rule out excessively small step lengths. This rule seeks to determine, for $0<\sigma_1<\sigma_2<1,$ a $\alpha^k > 0$ such that

$f(x^k) + \sigma_1\alpha^k\langle \nabla f(x^k), d^k \rangle \le f(x^k + \alpha^k d^k) \le f(x^k) + \sigma_2\alpha^k\langle \nabla f(x^k), d^k \rangle.$

## Wolfe
Wolfe's Rule complements Armijo's inequality with a curvature condition, stipulating that, in addition to the sufficient reduction in $f$, for $0<\sigma_1<\sigma_2<1,$ the step length $\alpha^k >0$ must comply

$f(x^k + \alpha^k d^k) \le f(x^k) + \sigma_1\alpha^k\langle \nabla f(x^k), d^k \rangle$ 


$\langle \nabla f(x^k + \alpha^k d^k), d^k \rangle \ge \sigma_2\langle \nabla f(x^k), d^k \rangle.$

### Parameters
- x_k (vector): Current point in the iteration.
- fx_k (function): Objective function to be minimized at point x_k.
- gradf_x (function): Gradient of the objective function.
- d_k (vector): Direction of descent.


