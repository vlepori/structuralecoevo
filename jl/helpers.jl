# Generate alpha matrix given a vector of trait values
function make_alpha(v::Any,p)
  n = length(v);
  out = Array{Float64,2}(undef, n,n);
  for i in eachindex(v)
    for j in eachindex(v)
      out[i,j]=alphafun(v[i],v[j],p)
    end
  end
  return(out)
end

function make_r(v::Any,p)
  map(x->rfun(x,p),v)  
end

function Nstar(x,p)
  alpha = make_alpha(x,p)
  r = map(y -> rfun(y,p), x)
  out = alpha\r # Solves system, works for (S=1)  too (N=r/a)
  # if any(out .< 0) error("Negative abundances") end
  return(out)
end

# Make a PIP for given range of x
function PIP(ran,p)
  (length(ran)==2) && (ran=range(ran[1],stop=ran[2],length=200))
  rnge = collect(ran)
  nr = length(rnge)
  out = Array{Float64,2}(undef, nr,nr);
  for i in 1:nr
    for j in 1:nr
      out[j,i]= sign(rfun(rnge[j],p) - alphafun(rnge[j],rnge[i],p)*(rfun(rnge[i],p)/alphafun(rnge[i],rnge[i],p)))
    end
  end
  heatmap(rnge,rnge,out)
end

# reshapes an array from wide to long format
function gather(array,time=nothing)
    out = [0,0,0]';
    for i in 1:length(array)
        tmpx = array[i]
        tmpi = fill(i, length(tmpx))
        t = (time == nothing) ? 0 : time[i]
        tmpt = fill(t,length(tmpx))
        chunk = hcat(tmpi,tmpt,tmpx)
        out = vcat(out,chunk)
    end
    return out[2:end,:]
end

# gathers + plot a trajectory from array
function outplot(array,time=nothing;logtime=false)
  out = [0,0,0]';
  for i in 1:length(array)
      tmpx = array[i]
      tmpi = fill(i, length(tmpx))
      t = (time == nothing) ? 0 : time[i]
      tmpt = fill(t,length(tmpx))
      chunk = hcat(tmpi,tmpt,tmpx)
      out = vcat(out,chunk)
  end
  out=out[2:end,:]
  if time==nothing
    plot(out[:,1], out[:,3], seriestype = :scatter,
        markersize = 3, markeralpha = 0.5, markerstrokewidth = 0,
        markercolor = :black, legend = false)
    else
      if logtime
        plot(log.(out[:,2]),out[:,3], seriestype = :scatter,
            markersize = 3, markeralpha = 0.5, markerstrokewidth = 0,
            markercolor = :black, legend = false)
      else
        plot((out[:,2]),out[:,3], seriestype = :scatter,
            markersize = 3, markeralpha = 0.5, markerstrokewidth = 0,
            markercolor = :black, legend = false)
      end
    end
end

# returns matrix
function draw_mc(range::Array,S::Int, n::Int, p)
    out = Array{Float64}(undef, n, S)
    i = 1
    while i <= n
        tmp = range[1] .+ rand(S)*(range[2]-range[1])
        N = [-1,-1]
        try
            N = Nstar(tmp,p)
        catch
            println("Singularity")
        end
        if prod(N .> 0)
            out[i,:] = tmp
            i += 1
        end
    end
    return out
end

function distance(a,b)
  sqrt(sum((a-b).^2))  
end

# method for S array -> draw communities of varying size. 
# different ranges can be supplied for different S (array of array)
# returns array of array
function draw_mc(range::Array,S::Array, n::Array, p; min_d)
    length(n) == 1 || length(n) == length(S) || error("n is neither constant nor length S")
    length(n) == 1 && (n=fill(n[1],length(S)))
    out = []

    if (length(range)==2) & (length(range[1])==1)
        range = fill(range,length(S)) # length for array of arrays is the number of arrays
    end
    (length(range) == length(S)) || error("length of range should be 2 or length of S!")
    (length.(range) == fill(2,length(range))) || error("ranges should have a max and a min")

    for s in collect(1:length(S))
        println("Starting s=",S[s])
        j = 1
        i = 1
        while j <= n[s]
            # println(range[s][1])
            tmp = range[s][1] .+ rand(S[s])*(range[s][2]-range[s][1])
            N = [-1,-1]
            try
                N = Nstar(tmp,p)
            catch
                println("Singularity")
            end
            if prod(N .> 0) 
              # only compare distances for same size communities
              euclid = (min_d & (j>1)) ? minimum(map(x -> distance(x, tmp), out[(end-j+2):end] )) : Inf 

              if euclid > .001
                push!(out, tmp)
                j += 1
              end
            end
            i +=1
            if i > (n[s]*1000) #acceptance below .001
              println("reached maxiter") 
              break
            end
        end
        println("Acceptance rate = ", round(j/i, digits=4))
    end
    return out
end


function logrange(start;stop,length_out)
  st_ = log(start)
  sp_ = log(stop)
  o_ = collect(range(st_,stop=sp_, length= length_out))
  return exp.(o_)
end



# import Arrow
# Arrow.write("df.feather",df)

import RCall
function save_RDS(obj::DataFrame, path::String)
    RCall.reval("library(dplyr)")
    rdf = RCall.robject(obj) # necessary?
    RCall.@rput rdf
    RCall.reval("rdf=as_tibble(rdf)")
    RCall.reval("saveRDS(rdf, file = '$path')")
end

# extract from array of array the elements where size changes
function extract_branchings(a::Array; after_ = false, last_ = true, indices = false)
  branchings = findall(x->x==1,  diff(map(length,a)))
  after_ && (branchings = sort([branchings ; branchings .+ 1])) # add first points after branching
  last_ && push!(branchings, length(a)) # add last state
  indices && return branchings # return indices only?
  return a[branchings]
end



# Simulate once for a given initial condition and set of parameters
function sim_ad(init, p ; t = [0,Inf], tol = [1e-13,1e-11], kall = false)

    prob = ODEProblem(g,init,t, p)
    cb = CallbackSet(d_extinction,mutation)
    # sol = solve(prob, Tsit5(), callback=cb, abstol=tol[1], reltol=tol[2], maxiters=1e7)
    sol = solve(prob, Rosenbrock23(autodiff=false), callback=cb, abstol=tol[1], reltol=tol[2], maxiters=1e7)


    save_id = 
    [collect(range(1,stop = length(sol.t),length = 5000));
    logrange(1, stop=length(sol.t),length_out=5000)] |>  x -> Int.(unique(round.(sort(x))))



    if length(sol.t) > 10000 & kall == false
        solN = sol.u[save_id]
        solt = sol.t[save_id]
    else
        solN = sol.u
        solt = sol.t
    end

    return solt, solN
end

# do multiple simulations across a range of parameters - keep only branchings
function sim_loop(par_array, init)
    
    df = DataFrame(sigmar = Float64[], par = Dict[], ss = Array[], stability = String[])

    for i in par_array
        tmpt,tmpN = sim_ad(init, i, kall = true)

        tmpN = extract_branchings(tmpN)     # keep only branching points
        tmpP = fill(i[:sigma_r], length(tmpN))
        tmpPar = fill(i, length(tmpN))
        tmpS = [fill("branching", length(tmpN)-1);"css"]

        append!(df, DataFrame(sigmar = tmpP, par = tmpPar, ss = tmpN, stability = tmpS))
        println("finished.")
    end

    return df
end


function trajectory_ranges(x::Array)
    branchings = findall(y->y==1,  diff(map(length,x)))
    push!(branchings, length(x)) # add last state
    ranges = map(y->collect(extrema(y)),x[branchings])
    return ranges
end

function isbetween(x::Real,lo::Real,hi::Real)
    (x>lo) & (x<hi)    
end

isbetween(x::Real, v::NTuple{2,Real}) = isbetween(x,v[1],v[2])