# 
function draw_low(range::Array,S::Int,n::Int,p)
    out = []
    singularities = 0
    j = 1
    # i = 1
        while j <= n
            tmp = range[1] .+ rand(S)*(range[2]-range[1])
            N = [-1,-1]
            try
                N = Nstar(tmp,p)
            catch
                singularities += 1
            end
            if prod(N .> 0) 
                push!(out, tmp)
                j += 1
            end
            # i +=1
            # if i > (n[s]*1000) #acceptance below .001
            #   println("reached maxiter") 
            #   break
            # end
        end 
    println("Singularities: ", singularities)   
    return out
end



# 
function draw_low_pl(range::Array,S::Int,n::Int,p)
    outout = []
    singularities = 0
    tries = 0
    batch = 2000 # initial guess
    while(length(outout) < n)
        print("\rstatus: ", round(length(outout)/ n, digits=2))
        out = Vector{Vector{Vector{Float64}}}()
        for i in 1:Threads.nthreads()
            push!(out,Float64[])
        end   
        Threads.@threads for i in 1:batch # increase batch size if rejection rate high
            tmp = range[1] .+ rand(S)*(range[2]-range[1])
            N = [-1,-1]
            try
                N = Nstar(tmp,p)
            catch
                singularities += 1
            end
            if prod(N .> 0) 
                push!(out[Threads.threadid()], tmp)
            end
        end
        tries += batch
        append!(outout,vcat(out...))
        # estimate rejection rate, adapt batch size
        length(outout) > 5 && (batch=round(0.05*n/(length(outout)/tries)))  
    end
    println("...done. singularities: ", singularities)
    # should not waste too many calculations
    return outout[1:n]
end


# high-level wrapper
# do not parallelize at high-level
function draw_high(range::Array,S::Array, n, p; parallel = false)

    isa(range[1],Array) || (range = fill(range,length(S)))
    length(n)>1 || (n=fill(n[1],length(S)))

    parallel && println("using ",Threads.nthreads()," threads")
    out=[]
    for s in collect(1:length(S)) 
        println("starting s=",S[s])
        tmp = parallel ? draw_low_pl(range[s], S[s],n[s], p) : draw_low(range[s], S[s],n[s], p)
        out = append!(out,tmp)        
    end
    return out
end




