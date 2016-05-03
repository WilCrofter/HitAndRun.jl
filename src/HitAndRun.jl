module HitAndRun

using Distributions
using Mamba

export mcmc_tables, tables2chains, chisq_tables, snee, near_random, victorians

include("hrutils.jl")
include("contingency_tables.jl")

end # module
