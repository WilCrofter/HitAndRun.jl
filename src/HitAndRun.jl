module HitAndRun

using Distributions
using Mamba

export mcmc_tables, hitrun2chains, chisq_tables, snee, near_random, victorians

include("hrutils.jl")
include("contingency_tables.jl")
include("knapsack.jl")

end # module
