#= Mallow's model

Reference: "Hit and Run as a Unifying Device." by P. Diaconis & H. Andersen, (2007), *Journal de la Société Française de Statistique*, 148(4):5-28, [[PDF](http://statweb.stanford.edu/~cgates/PERSI/papers/hitandrun062207.pdf)].

The reference describes the following MCMC step for Mallow's model in which p(σ)=exp(Σ β[i]*ρ[i]*σ[i] )/Z. (ρ fixed and Z a normalizing factor.)


Choose u[j] uniformly on [0, exp(β*ρ[j]*σ[j])].
Given u, choose new σ uniformly on {σ | σ[j] ≤ b[j] = log(u[j])/(ρ[j]*β) } as follows:

Look at places j with b[j] ≤ 1 and place symbol 1 at a uniform choice among these. Look at places with b[j] ≤ 2 and place symbol 2 in a uniform choice among those (with the place where 1 was placed deleted). In general, look at places j with b[j] ≤ k. Of these k − 1 will be occupied. There will always be one or more available. Choose uniformly among these and place k there.

The last sentence is justified as follows

u[j] ≤ exp(β[j]*ρ[j]*σ[j]) =>
b[j]=log(u[j])/(ρ[j]*β[j]) ≤ σ[j], hence
b[j] ≤ j for j=1,...,n, hence
b[i] ≤ j for i=1,...,j.

=#
function step_mallow!{T<:Real}(beta::Array{T,1}, target::Array{Int,1}, state::Array{Int,1})
    n = length(target)
    # u = rand(n) .* exp(beta .* target .* state)  =>
    # log(u) = log(rand(n)) + beta .* target .* state =>
    b = (log(rand(n))./ (beta .* target)) + state 
    # Note that state=target will satisfy the suppositions
    # necessary to complete the loop.
    for j in 1:n
        i = rand(find(b .<= j))
        state[j] = i
        b[i] = n+1
    end
end

function mcmc_mallow{T<:Real, U<:Real}(beta::Array{T,1},burnin::Int, thin::Int, iterations::Int; target::Array{U,1}=sortperm(sortperm(beta)), state::Array{Int,1}=target)
    I = length(state)
    I == length(beta) || error("beta and initial state must have the same length.")
    I == length(target) || error("beta and target must have the same length.")
    all( beta .> 0) || error("All elements of beta must be positive.")
    N = length((burnin+thin):thin:iterations)
    record = zeros(T,(N, I+1))    
    k = 0
    for i in 1:iterations
        step_mallow!(beta, target, state)
        if i > burnin && mod(i,thin)==0 && k < N
            k += 1;
            record[k,2:(I+1)] = state
            record[k,1] = dot(beta, state)
        end
    end
    return record
end
                             
