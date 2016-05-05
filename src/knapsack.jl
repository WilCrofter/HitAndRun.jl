#=

Given a vector of positive weights, w=(w_1,...,w_n), and a maximum weight, W, this chain converges to a uniform sampler of binary vectors a=(a_1,...,a_n) such that <a,w) <= W, i.e., to a uniform sampler of feasible solutions (not *optimal* feasible solutions) to a knapsack problem. The algorithm is that of reference [1]. 

Originally I had thought to combine the algorithms of [1] and [2] in one unit of code, but [1] involves a discrete random walk along the vertices of a truncated cube, while [2] involves a continuous walk in the interior of a convex set (which I thought to specialize to a truncated unit cube.) The diversity of code required argues for separate treatment, so this file implements only the algorithm of [1].

[1] B Morris, A Sinclair [Random walks on truncated cubes and sampling 0-1 knapsack solutions](https://scholar.google.com/scholar?q=RANDOM+WALKS+ON+TRUNCATED+CUBES+AND+SAMPLING+0-1+KNAPSACK+SOLUTIONS&btnG=&hl=en&as_sdt=0%2C21) SIAM journal on computing, 2004 - SIAM

[2] L Lovász, S Vempala, [Hit and Run from a Corner](https://scholar.google.com/scholar?hl=en&q=Hit-and-Run+from+a+Corner&btnG=&as_sdt=1%2C21&as_sdtp=) SIAM Journal on Computing, 2006 - SIAM

=#


#= A knapsack problem consists of:
   a vector of real, positive weights.
   a maximum weight.
These specifications are simple enough. There's little point in combining them in a special type.
 =#

""" init_knapsack{T<:Real}(item_weights::Array{T,1}, max_weight::T)

  Fill the knapsack with items in random order, skipping over items which would
cause the weight limit to be exceeded.
"""
function init_knapsack{T<:Real}(item_weights::Array{T,1}, max_weight::T)
    all(item_weights .> 0) || error("Item weights must all be positive.");
    max_weight > 0 || error("Maximum weight must be positive.");
    W = 0.0;
    n = length(item_weights);
    a = zeros(T,n);
    perm = randperm(n);
    for i in 1:n
        k = perm[i];
        w = item_weights[k];
        if W+w <= max_weight
            W += w;
            a[k] = 1;
        end
    end
    return a;
end


""" step_knapsack{T<:Real}(state::T, item_weights::Array{T,1}, max_weight::T)

Step the knapsack chain using the Morris-Sinclair algorithm:
  1. With probability 1/2 do nothing, to avoid a periodic chain.
  2. Pick an item, i, at random
  3. If i is in the knapsack take it out
  4. If i in not int knapsack put it in, provided doing so doesn't exceed the weight limit.

"""
function step_knapsack{T<:Real}(state::Array{T,1}, current_weight::T, item_weights::Array{T,1}, max_weight::T)
    # do nothing with probability 1/2 to avoid a periodic chain
    if rand(1)[1] < .5
        return state, current_weight;
    end
    # pick a random item
    item = rand(1:length(item_weights));
    if state[item] > 0
        # if it is in the knapsack take it out
        state[item] = 0;
        current_weight -= item_weights[item];
    elseif current_weight + item_weights[item] <= max_weight
        # else put it in, provided it doesn't cause weight limit to be exceeded
        state[item] = 1;
        current_weight += item_weights[item];
    end
    return state, current_weight
end

""" mcmc_knapsack{T<:Real}(item_weights::Array{T,1}, max_weight::T, burnin::Int, thin::Int, iterations::Int, initial_state::Array{T,1}=init_knapsack(item_weights,max_weight))

      
"""
function mcmc_knapsack{T<:Real}(item_weights::Array{T,1}, max_weight::T, burnin::Int, thin::Int, iterations::Int, initial_state::Array{T,1}=init_knapsack(item_weights,max_weight))
    I = length(item_weights);
    N = length((burnin+thin):thin:iterations);
    record = zeros(T,(N, I+1));
    state=initial_state;
    w = dot(item_weights, state);
    k = 0;
    for i in 1:iterations
        state, w = step_knapsack(state, w, item_weights, max_weight);
        if i > burnin && mod(i,thin)==0 && k < N
            k += 1;
            record[k,2:(I+1)] = state;
            record[k,1] = w;
        end
    end
    return record;
end

#= Convert to Mamba.Chains object for convergence tests
=#


