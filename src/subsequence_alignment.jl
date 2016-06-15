#=
References:

[1] "Hit and Run as a Unifying Device." by P. Diaconis & H. Andersen, (2007), *Journal de la Société Française de Statistique*, 148(4):5-28, [[PDF](http://statweb.stanford.edu/~cgates/PERSI/papers/hitandrun062207.pdf)].

[2] "Algebraic algorithms for sampling from conditional distributions" by P. Diaconis & B. Sturmfels, The Annals of Statistics 1998, Vol. 26, No. 1, 363-397, [[PDF](http://projecteuclid.org/download/pdf_1/euclid.aos/1030563990)]

An attempt to apply a hit and run technique to multiple sequence alignment.

A set of strings with a special character added to make them all the same length is given. They are aligned as (say) rows of a matrix. A column’s measure of misalignment is a pairwise distance between distinct pairs in the column. Permute the spaces in each row so as to minimize total misalignment (summed over all columns.)

We impose an exponential penalty on misalignment and use a hit and run technique somewhat similar to that used in [1] for the Mallows model. Given a putative alignment, A(t) with score s(t), we choose a row i and a special character in that row at random. We choose a random number, u(t), uniformly on [0, s(t)]. We then exhaustively score alignments with all alternate placements of the chosen special character within the chosen row. We choose the putative alignment A(t+1) uniformly among those which score at least as well as u(t).

=#
