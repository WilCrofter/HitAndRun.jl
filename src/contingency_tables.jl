
#= Generating contingency tables.

=#

""" choose2(m::Int)
Return a  uniformly random pair of distinct integers from {1,...,m}.    
"""
function choose2(m::Int)
    m > 1 || error("m must exceed 1.");
    a = rand(1:m);
    b = rand(2:m);
    return [a, b == a ? 1 : b];
end

"""
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

"""
Given a contingency table, T, replace a random 2x2 subtable using the default uniform distribution or the optional hypergeometric.
"""
function step!(T::Array{Int,2}, uniform=true)
    n, m = size(T);
    i = choose2(n);
    j = choose2(m);
    T[i,j] = uniform ? u_sub(T[i,j]) : h_sub(T[i,j]);
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

