#=

Given a vector of positive weights, w=(w_1,...,w_n), and a maximum weight, W, this chain converges to a uniform sampler of binary vectors a=(a_1,...,a_n) such that <a,w) <= W, i.e., to a uniform sampler of feasible solutions (not optimal feasible solutions) to a knapsack problem. The algorithm is that of reference [1]. 

Originally I had thought to combine the algorithms of [1] and [2] in one unit of code, but [1] involves a discrete random walk along the vertices of a truncated cube, while [2] involves a continuous walk in the interior of a convex set (which I thought to specialize to a truncated unit cube.) The diversity of code required argues for separate treatment, so this file implements only the algorithm of [1].

[1] B Morris, A Sinclair [Random walks on truncated cubes and sampling 0-1 knapsack solutions](https://scholar.google.com/scholar?q=RANDOM+WALKS+ON+TRUNCATED+CUBES+AND+SAMPLING+0-1+KNAPSACK+SOLUTIONS&btnG=&hl=en&as_sdt=0%2C21) SIAM journal on computing, 2004 - SIAM

[2] L Lovász, S Vempala, [Hit and Run from a Corner](https://scholar.google.com/scholar?hl=en&q=Hit-and-Run+from+a+Corner&btnG=&as_sdt=1%2C21&as_sdtp=) SIAM Journal on Computing, 2006 - SIAM

=#


#= Define problem:
   A vector of weights (Float64).
   A maximum weight.
These are simple enough. There's little point in combining them in a special type.
 =#

#= Initialize
  One can always begin from an empty knapsack. Alternatively, one could begin by filling the knapsack randomly until the weight limit is reached. Is this significantly different from stepping the algorithm?
=#

#= Step
=#

#= MCMC
=#

#= Compute #items, total weight
=#

#= Convert to Mamba.Chains object for convergence tests
=#


