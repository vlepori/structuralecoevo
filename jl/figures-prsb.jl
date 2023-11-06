# This file generates the simulations for the figures.
# Use together with R/figures.r 

using Plots
using LinearAlgebra
using DifferentialEquations
using StatsBase
using Distributions
using DataFrames
import Serialization
# using Zygote # only for autodiff

include("helpers.jl")
include("ode.jl")
include("draw.jl")
include("r_norm.jl")
include("a_mixed.jl")

mkdir("out")
#*******************************************************

parameters = Dict{Symbol,Float64}(:a_r => 10.0, :mu_r => 0.0, :sigma_r => 1.6, # 3.5
    :a => 3.5, :sigma => 1.0,
    :c => 1, :v => 1, :k => 1,
    :m => 0.1,
    :beta => 0)

initial = [1.00]

t, xs = sim_ad(initial, parameters; tol=[1e-12, 1e-12], kall=true);
sim1 = DataFrame(x=xs, t=t)
sim1[!, "rs"] = map(x -> make_r(x.x, parameters), eachrow(sim1))
sim1[!, "alpha"] = map(x -> make_alpha(x.x, parameters), eachrow(sim1))

save_RDS(sim1, "out/fig1.rds")

## Randomization
# First : within ranges
branchings = findall(y -> y == 1, diff(map(length, xs)))
push!(branchings, length(xs)) # add last state
ranges = map(y -> collect(extrema(y)), xs[branchings])
comm_size = length(xs[end])

ns = 1000 * (2 .^ (1:comm_size))
mc = draw_high(ranges, collect(1:comm_size), ns, parameters; parallel=true)

mcdf = DataFrame(x=mc)
mcdf[!, "rs"] = map(y -> make_r(y.x, parameters), eachrow(mcdf))
mcdf[!, "alpha"] = map(y -> make_alpha(y.x, parameters), eachrow(mcdf))

save_RDS(mcdf, string("out/mc.rds"))

# Second
x = collect(-100:0.01:100)
y = map(z -> rfun(z, parameters), x)

xmin_rpos = minimum(x[y.>0])
xmax_rpos = maximum(x[y.>0])

ns = 1000 * (2 .^ (1:comm_size))

mc2 = draw_high([xmin_rpos, xmax_rpos], collect(1:comm_size), ns, parameters; parallel=true)


mcdf2 = DataFrame(x=mc2)
mcdf2[!, "rs"] = map(y -> make_r(y.x, parameters), eachrow(mcdf2))
mcdf2[!, "alpha"] = map(y -> make_alpha(y.x, parameters), eachrow(mcdf2))

save_RDS(mcdf2, string("out/mc2.rds"))
