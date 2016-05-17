module HitAndRun

using Distributions
using Mamba

export mcmc_tables, hitrun2chains, chisq_tables, snee, near_random, victorians
export init_knapsack, mcmc_knapsack
export mcmc_mallow

include("hrutils.jl")
include("contingency_tables.jl")
include("knapsack.jl")
include("mallow.jl")

end # module
