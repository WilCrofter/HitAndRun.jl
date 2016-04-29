
#= Generating contingency tables.

=#

""" choose2(m::Int)
Return a sorted, uniformly random pair of distinct integers from {1,...,m}.    
"""
function choose2(m::Int)
    m > 1 || error("m must exceed 1.");
    a = rand(1:m);
    b = rand(2:m);
    return sort([a, b == a ? 1 : b]);
end

""" function h_sub(T::Array{Int,2})
Given a 2x2 table, return another with the same row and column sums as chosen by a hypergeometric distribution on the 1,1 entry.
"""
function h_sub(T::Array{Int,2})
    ns = sum(T[1,:]); # "successes"
    nf = sum(T[2,:]); # "failures" 
    draws = sum(T[:,1]); # 1st column represents number of draws
    #  draw from the associated hypergeometric
    a11 = rand(Distributions.Hypergeometric(ns,nf,draws));
    a12 = ns - a11;
    a21 = draws-a11;
    a22 = sum(T[:,2])-a12;
    return reshape([a11,a21,a12,a22],(2,2));
end

"""
Given a 2x2 table, return another with the same row and column sums as chosen by a uniform distribution on the 1,1 entry.
"""
function u_sub(T::Array{Int, 2})
    s1, s2 = sum(T,2); # row sums
    t1, t2 = sum(T,1); # column sums
    lb = max(s1-t2,0); # lower bound for 1,1 entry
    ub = min(s1,t1); # upper bound for 1,1 entry
    a = rand(Distributions.DiscreteUniform(lb, ub));
    return reshape([a, t1-a, s1-a, t2-(s1-a)],(2,2));
end

""" step!(T::Array{Int,2}, uniform=true)
Given a contingency table, T, replace a random 2x2 subtable using the default hypergeometric distribution or the optional uniform.
"""
function step!(T::Array{Int,2}, uniform=false)
    n, m = size(T);
    # Generate index i,j
    i = choose2(n);
    j = choose2(m);
    # Generate new table given i,j
    T[i,j] = uniform ? u_sub(T[i,j]) : h_sub(T[i,j]);
end

""" mcmc_tables(T::Array{Int,2}, burnin::Int, thin::Int, iterations::Int; uniform::Bool=true, name::AbstractString="T")
Using `T` as the initial state, step `burnin` times without saving, then save every `thin` times thereafter until a total of `iterations` steps have been taken. Return an size(T) by N array of samples, where N is length((burnin+thin):thin:iterations).
"""
function mcmc_tables(T::Array{Int,2}, burnin::Int, thin::Int, iterations::Int; uniform::Bool=true, name::AbstractString="T")
    state = deepcopy(T);
    I, J = size(state);
    N = length((burnin+thin):thin:iterations);
    record = zeros(Int,(I,J,N));
    k = 0;
    for i in 1:iterations
        step!(state, uniform);
        if i > burnin && mod(i,thin)==0
            k += 1;
            record[:,:,k] = state;
        end
    end
    return record;
end


""" tables2chains(record::Array{Int,3}, burnin::Int, thin::Int; name::AbstractString="T")
    
Convert an I,J,K array of IxJ contingency tables to an object, C, of type Mamba.Chains in which C.value[K,:] contains the entries of table K, reshaped to a vector.
"""
function tables2chains(record::Array{Int,3}, burnin::Int, thin::Int; name::AbstractString="T")
    I,J,N = size(record);
    return Mamba.Chains(reshape(record, (I*J, N))',
                        start=burnin+thin,
                        thin=thin,
                        names=make_names(name,(I,J))
                        );
end

"""  chisq_tables(T::Array{Int,2})
Compute Pearson's chi-squared statistic for the given, IxJ, contingency table. Under the null hypothesis, that entries in each column represent independent draws from a distribution common to all columns, the statistic should be approximately chi-squared with (I-1)(J-1) degrees of freedom.
"""    
function chisq_tables(T::Array{Int,2})
    ex = sum(T,2)*sum(T,1)/sum(T); # expectations under H0, independence
    return sum((T-ex).^2 ./ ex); 
end

"""
If the argument is an IxJxK array, the statistic will be computed for each IxJ table in the array.    
"""    
function chisq_tables(record::Array{Int, 3})
    return [chisq_tables(record[:,:,n]) for n in 1:size(record,3)];
end


""" snee

Contingency table, eye color (row) vs hair color (column), from Snee, Graphical display of two-way contingency tables. Amer. Statist. 38 9᎐12. (1974)

	Black 	Brunette Red 	Blonde
Brown	68	119	26	7
Blue	20	84	17	94
Hazel	15	54	14	10
Green	5	29	14	16
"""
snee = [68 119 26 7;
        20 84 71 94;
        15 54 14 10;
        5  29 14 16];

""" near_random
A contingency table with the same row and column sums as snee, whose entries are close to their expectations.
"""
near_random = [38 96 43 43;
               45 119 52 53;
               16 41 18 18;
               11 28 12 13];

""" victorians

Contingency table, birth month (row) vs death month (column) for 82 descendants of Queen Victoria, from Diaconis and Sturmfels, Algebraic Algorithms for Sampling from Conditional Distributions, The Annals of Statistics 1998, Vol. 26, No. 1, 363 ᎐ 397

    J F M A M J J A S O N D
Jan 1 0 0 0 1 2 0 0 1 0 1 0
Feb 1 0 0 1 0 0 0 0 0 1 0 2
Mar 1 0 0 0 2 1 0 0 0 0 0 1
Apr 3 0 2 0 0 0 1 0 1 3 1 1
May 2 1 1 1 1 1 1 1 1 1 1 0
Jun 2 0 0 0 1 0 0 0 0 0 0 0
Jul 2 0 2 1 0 0 0 0 1 1 1 2
Aug 0 0 0 3 0 0 1 0 0 1 0 2
Sep 0 0 0 1 1 0 0 0 0 0 1 0
Oct 1 1 0 2 0 0 1 0 0 1 1 0
Nov 0 1 1 1 2 0 0 2 0 1 1 0
Dec 0 1 1 0 0 0 1 0 0 0 0 0
"""
victorians =
    [1 0 0 0 1 2 0 0 1 0 1 0;
     1 0 0 1 0 0 0 0 0 1 0 2;
     1 0 0 0 2 1 0 0 0 0 0 1;
     3 0 2 0 0 0 1 0 1 3 1 1;
     2 1 1 1 1 1 1 1 1 1 1 0;
     2 0 0 0 1 0 0 0 0 0 0 0;
     2 0 2 1 0 0 0 0 1 1 1 2;
     0 0 0 3 0 0 1 0 0 1 0 2;
     0 0 0 1 1 0 0 0 0 0 1 0;
     1 1 0 2 0 0 1 0 0 1 1 0;
     0 1 1 1 2 0 0 2 0 1 1 0;
     0 1 1 0 0 0 1 0 0 0 0 0];

