# Main ODE
function g(du,u,p,t)
  N = Nstar(u,p)
  print("\ru=",u)
  if any(N.<0)
    println(N)
    println("Negative abundances!")
  end
  S = length(u)
  for i in 1:S
    du[i] = (drdx(u[i],p) - sum([dadx(u[i],u[j],p) for j in 1:S].*N)) # array comprehension
  # du[i] = N[i]*(drdx(u[i],p) - sum([dadx(u[i],u[j],p) for j in 1:S].*N)) # array comprehension

    # println(du[i])
  end
end

# Checks stability
function d2w(u,f,N,p)
  d2rdx(u[f],p)- sum([d2adx(u[f],u[j],p) for j in 1:length(u)].*N)
end


# Callbacks ------------------------------

function check_ss(u,t,integrator)
  DiffEqCallbacks.allDerivPass(integrator, 1e-10, 1e-10)
  # anyDerivPass(integrator, 1e-08, 1e-10) # should match steady val below!
end

function eval_ss!(integrator)
  abs_du = abs.(get_du(integrator))
  println("\ndu=",  abs_du)
  p = integrator.p
  u = integrator.u
  s = length(u)
  if s > 8
    println("Terminating further branching, S=",s)
    terminate!(integrator)
    return
  end
  N = Nstar(u,p)
  println("Nstar=",N)
  if any(N.<0) error("Neg abunds!") end
  d2 = [d2w(u,f,N,p) for f in 1:s]
  println("d2=",d2," ")
    if any(d2 .> 0)
      # steady = findall(abs_du .< 1e-07) # idx morphs for which dx is approaching zero - should match abstol above!
      steady = 1:s
      unstable = findall(d2 .> 0) # indexes
      branch = intersect(steady,unstable)
      println("Branching:",branch)
      if length(branch) == 0
        println("?")
        terminate!(integrator)
      end
      #***
      # if length(branch) > 1
      #   branch=sample(branch, Weights(N[branch]), 1) # N[branch]
      #   println("Random mutant:",branch)
      # end
      #***
      if length(branch) > 1
        # branch = branch[findmax(d2[branch])[2]] #subset du -> find min -> new id -> old id. Branch the most unstable?
        branch = branch[findmax(N[branch])[2]] # branch most abundant?
        println("Mutant:",branch)
      end
      #***
      x_rs = u[branch]
      s = length(u)
      resize!(integrator,s+length(branch)) # initialize N_mutatnt
      for i in 1:length(branch)
        tmp = u[branch[i]]
        u[branch[i]] = tmp-0.01
        u[s+i] = tmp+0.01
      end
      println("\nN=",Nstar(u,p),"\nu=",u,"\nt=", integrator.t)
      if any(Nstar(u,p).<0) error("Neg abundances after branch!") end # not guaranteed to be feasible...
    else
      println("CSS")
      terminate!(integrator)
    end
end

mutation = DiscreteCallback(check_ss,eval_ss!)


# Continuous callback
# check abundance, return least abundant
# continuous callbacks trigger effect when f(u) = 0
function check_abundance(u,t,integrator)
  p = integrator.p
  return(minimum(Nstar(u),p)-0.001)
end

# always true - for DiscreteCallback
function check_abundance_d(u,t,integrator)
  true
end


function extinct!(integrator)
  u = integrator.u
  p = integrator.p
  N = Nstar(u,p)
  if any(N .< 0.001)
    idx = findall(N .< 0.001)
    println("N=",N," deleting at position ",idx)
    deleteat!(integrator, idx)
    println("New N=", Nstar(u,p))
  end
  nothing
end

extinction = ContinuousCallback(check_abundance,extinct!)
d_extinction = DiscreteCallback(check_abundance_d,extinct!)

# Checks that dx/dt below threshold
function anyDerivPass(integrator, abstol, reltol)
    if DiffEqBase.isinplace(integrator.sol.prob)
        testval = first(get_tmp_cache(integrator))
        DiffEqBase.get_du!(testval, integrator)
        if typeof(integrator.sol.prob) <: DiffEqBase.DiscreteProblem
            @. testval =  testval - integrator.u
        end
    else
        testval = get_du(integrator)
        if typeof(integrator.sol.prob) <: DiffEqBase.DiscreteProblem
            testval =  testval - integrator.u
        end
    end
    all(abs(d) > abstol && abs(d) > reltol*abs(u) for (d,abstol, reltol, u) =
           zip(testval, Iterators.cycle(abstol), Iterators.cycle(reltol), integrator.u)) && (return false)
    return true
end

# simulate LV models:
# from vector of x -> needs rfun, alphafun
# function ode_x(x, u0 = "Ks")
#   function f(du,u,p,t) # need to unpack p instead of using global variables
#     s_rs, s_alphas = p
#     temp = s_alphas*u
#     for i in 1:length(s_rs)
#       du[i] =  u[i]*(s_rs[i] - temp[i])
#     end
#   end
#   r_ = map(rfun, x)
#   a_ = make_alpha(x)
#   p = [r_,a_]

#   if u0 == "Ks"
#     u0 = r_./diag(a_)
#   elseif u0 == "Cst"
#     u0 = fill(mean(r_./diag(a_)),length(x))
#   elseif !isa(u0,String)
#     if length(u0) == length(x)
#       u0 = u0
#     elseif length(u0) == 1
#       u0 = fill(u0[1],length(x))
#     else
#       error("Invalid u0!")
#     end
#   else
#     error("Invalid u0")
#   end

#   clbck = TerminateSteadyState()
#   tspan = (0.0,10.0)
#   prob = ODEProblem(f,u0,tspan,p)
#   cb = CallbackSet(d_extinction,mutation)
#   out = solve(prob,Tsit5(), callback = clbck, abstol = 1e-8, reltol = 1e-6, maxiters = 1e9);
#   return out
# end


# # from vector of r and matrix of alphas
# function ode_r_alpha(r, alpha, u0 = "Ks")
#   function f(du,u,p,t) # need to unpack p instead of using global variables
#     s_rs, s_alphas = p
#     temp = s_alphas*u
#     for i in 1:length(s_rs)
#       du[i] =  u[i]*(s_rs[i] - temp[i])
#     end
#   end
#   p = [r,alpha]
#   if u0 == "Ks"
#     u0 = r./diag(alpha)
#   elseif u0 == "Cst"
#     u0 = fill(mean(r./diag(alpha)),length(r))
#   elseif !isa(u0,String)
#     if length(u0) == length(r)
#       u0 = u0
#     elseif length(u0) == 1
#       u0 = fill(u0[1],length(r))
#     else
#       error("Invalid u0!")
#     end
#   else
#     error("Invalid u0")
#   end
#   clbck = TerminateSteadyState()
#   tspan = (0.0,10.0)
#   prob = ODEProblem(f,u0,tspan,p)
#   cb = CallbackSet(d_extinction,mutation)
#   out = solve(prob,Tsit5(), callback = clbck, abstol = 1e-8, reltol = 1e-6, maxiters = 1e9);
#   return out
# end