
function step_mallow!{T<:Real}(betarho::Array{T,1}, state::Array{Int,1})
    n = length(state)
    # u = rand(n) .* exp(betarho .* state)  =>
    # log(u) = log(rand(n)) + betarho .* state =>
    b = (log(rand(n))./ betarho) + state 
    # Note that state=target will satisfy the suppositions
    # necessary to complete the loop.
    for j in 1:n
        i = rand(find(b .<= j))
        state[j] = i
        b[i] = n+1
    end
end

""" mcmc_mallow{T<:Real, U<:Real}(beta::Array{T,1}, target::Array{U,1}, burnin::Int, thin::Int, iterations::Int; state::Array{Int,1}=target)

Return a record of MCMC samples of Mallow's model as described below.

(If beta is a number rather than an array, first expand it to an array.)

Reference: "Hit and Run as a Unifying Device." by P. Diaconis & H. Andersen, (2007), *Journal de la Société Française de Statistique*, 148(4):5-28, [[PDF](http://statweb.stanford.edu/~cgates/PERSI/papers/hitandrun062207.pdf)].

The reference describes the following MCMC step for Mallow's model in which p(σ)=exp(Σ β[i]*ρ[i]*σ[i] )/Z. (ρ fixed and Z a normalizing factor.) Equivalently, p(σ)=exp(Σ β[i]*(ρ[i]-σ[i])^2)/Z.

Choose u[j] uniformly on [0, exp(β*ρ[j]*σ[j])].
Given u, choose new σ uniformly on {σ | σ[j] ≤ b[j] = log(u[j])/(ρ[j]*β) } as follows:

Look at places j with b[j] ≤ 1 and place symbol 1 at a uniform choice among these. Look at places with b[j] ≤ 2 and place symbol 2 in a uniform choice among those (with the place where 1 was placed deleted). In general, look at places j with b[j] ≤ k. Of these k − 1 will be occupied. There will always be one or more available. Choose uniformly among these and place k there.

The last sentence is justified as follows

u[j] ≤ exp(β[j]*ρ[j]*σ[j]) =>
b[j]=log(u[j])/(ρ[j]*β[j]) ≤ σ[j],
Since σ is a permutation of {1,...,n} there are at least k values of i such that σ[i] ≤ k, hence such that b[j] ≤ k. 


The two fixed parameters, β and ρ, always appear as a product, hence can be replaced by a single fixed parameter representing the product. This is done internally, but the function signature distinguishes between them, where ρ is called the target.

"""
function mcmc_mallow{T<:Real, U<:Real}(beta::Array{T,1}, target::Array{U,1}, burnin::Int, thin::Int, iterations::Int; state::Array{Int,1}=target)
    I = length(state)
    I == length(beta) || error("beta and initial state must have the same length.")
    I == length(target) || error("beta and target must have the same length.")
    all( beta .> 0) || error("All elements of beta must be positive.")
    betarho = beta .* target
    N = length((burnin+thin):thin:iterations)
    record = zeros(T,(N, I+1))    
    k = 0
    for i in 1:iterations
        step_mallow!(betarho, state)
        if i > burnin && mod(i,thin)==0 && k < N
            k += 1;
            record[k,2:(I+1)] = state
            record[k,1] = dot(beta, state)
        end
    end
    return record
end
                             
function mcmc_mallow{T<:Real, U<:Real}(beta::T, target::Array{U,1}, burnin::Int, thin::Int, iterations::Int; state::Array{Int,1}=target)
     return mcmc_mallow(fill(beta, length(target)), target, burnin, thin, iterations; state=state)
end
