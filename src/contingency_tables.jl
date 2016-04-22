""" Generating contingency tables.
"""

""" Model

A Mamba model has 3 node types, Stochastic, Logical, and Input, which are natural to a Bayesian network, but somewhat unnatural to an ensemble of contingency tables with specified row and column sums. The shoe fits, however.

Input consists of everything fixed. In this case, input would be the initial table.

Logical nodes are deterministic functions of other nodes. In the present case these would be table entries.

Stochastic nodes would be 2x2 subtables, i.e., 2 row indices and 2 column indices, along with a random upper left entry. 

"""

"""
Return an initialized model, given an initial contingency table
"""
function contingency_table_model(T::Array{Int,2})
    rows, cols = size(T);
    rowsums = sum(T,2);
    colsums = sum(T,1);
    model = Model(
                  
                  );
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

