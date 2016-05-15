#= Mallow's model

Reference: "Hit and Run as a Unifying Device." by P. Diaconis & H. Andersen, (2007), *Journal de la Société Française de Statistique*, 148(4):5-28, [[PDF](http://statweb.stanford.edu/~cgates/PERSI/papers/hitandrun062207.pdf)].

The reference describes the following MCMC step for Mallow's model in which p(σ)=exp(Σ β[i]*i*σ[i] )/Z. (Z a normalizing factor.) Given u, choose new σ uniformly on {σ | σ[j] ≤ b[j] = log(u[j])/β } as follows:

Look at places j with b[j] ≤ 1 and place symbol 1 at a uniform choice among these. Look at places with b[j] ≤ 2 and place symbol 2 in a uniform choice among those (with the place where 1 was placed deleted). In general, look at places j with b[j] ≤ k. Of these k − 1 will be occupied. There will always be one or more available. Choose uniformly among these and place k there.

The last sentence is justified as follows

u[j] ≤ exp(β[j]*j*σ[j]) <=> b[j]=log(u[j])/β[j] ≤ j*σ[j], hence
b[i] ≤ j for i=1,...,j.

But by the same argument we also have b[j]/β[j] ≤ σ[j]. Hence,

b[σ<sup>-1</sup>[j]/β[σ<sup>-1</sup>[j]] ≤ j

This allows a generalization in which p(σ)= exp(Σ β[i]*σ[i] )/Z for any positive  β[i].

=#
function step_mallow{T<:Real}(beta::Array{T,1}, sigma::Array{Int,1})
    n = length(sigma)
    u = rand(n) .* exp(beta .* sigma)
    b = log(u) ./ beta
    # Apply sigma inverse to b.
    b = b[sortperm(sigma)]
    # Form new permutation. Note that sigma itself will satisfy
    # the suppositions necessary to complete the for loop.
    s = zeros(Int, n)
    for i in 1:n
        s[i] = rand(find(b .<= i))
        b[s[i]] = n+1;
    end
    return s
end

function mcmc_mallow{T<:Real}(beta::Array{T,1},burnin::Int, thin::Int, iterations::Int; state=randperm(length(beta)))
    I = length(state)
    I == length(beta) || error("beta and initial state must have the same length.")
    all( beta .> 0) || error("All elements of beta must be positive.")
    N = length((burnin+thin):thin:iterations)
    record = zeros(T,(N, I+1))    
    k = 0
    for i in 1:iterations
        state = step_mallow(beta, state)
        if i > burnin && mod(i,thin)==0 && k < N
            k += 1;
            record[k,2:(I+1)] = state
            record[k,1] = dot(beta, state)
        end
    end
    return record
end
                             
