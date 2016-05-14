#= Mallow's model
=#

#= Do a straight Mallow's model first, then try to generalize. Use the identity as central wlog. (Just apply the real central to results to compensate.)
=#

#= Given sigma choose u[i] uniformly on [0, exp(beta*i*sigma[i])]
=#

#= Given u, choose new sigma uniformly on sigma[j] ≤ b[j] = log(u[j])/beta as follows:

Look at places j with b[j] ≤ 1 and place symbol 1 at a uniform choice among these. Look at places with b[j] ≤ 2 and place symbol 2 in a uniform choice among those (with the place where 1 was placed deleted). In general, look at places j with b[j] ≤ k. Of these k − 1 will be occupied. There will always be one or more available. Choose uniformly among these and place k there.
=#

#= Can't the algorithm above be applied to any set of numbers for which there is at least k numbers ≤ k for each k? But how is that guaranteed? It must be by the choice of the u[j]. That is, the u from the previous sigma must satisfy it:

u[j] ≤ exp(beta[j]*j*sigma[j])
<=> log(u[j])/beta[j] ≤ j*sigma[j] ≤ j, hence
   b[i] ≤ j for i=1,...,j.

But suppose the pmf is redefined: p(sigma) = exp(dot(b',sigma)) where b'[i] = i*b[i]. The same step works for any positive b', whether associated with a permutation or not. 
=#

function step_mallow{T<:Real}(beta::Array{T,1}, sigma::Array{Int,1})
    all( beta .> 0) || error("All elements of beta must be positive.")
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
