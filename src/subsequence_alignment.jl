#=
References:

[1] "Hit and Run as a Unifying Device." by P. Diaconis & H. Andersen, (2007), *Journal de la Société Française de Statistique*, 148(4):5-28, [[PDF](http://statweb.stanford.edu/~cgates/PERSI/papers/hitandrun062207.pdf)].

[2] "Algebraic algorithms for sampling from conditional distributions" by P. Diaconis & B. Sturmfels, The Annals of Statistics 1998, Vol. 26, No. 1, 363-397, [[PDF](http://projecteuclid.org/download/pdf_1/euclid.aos/1030563990)]

An attempt to apply a hit and run technique to multiple sequence alignment.

A set of strings with a special character added to make them all the same length is given. They are aligned as (say) rows of a matrix. A column’s measure of misalignment is a pairwise distance between distinct pairs in the column. Permute the spaces in each row so as to minimize total misalignment (summed over all columns.)

We impose an exponential penalty on misalignment and use a hit and run technique somewhat similar to that used in [1] for the Mallows model. Given a putative alignment, A(t), with score, s(t), we choose a row, i, and a special character in that row at random. We choose a random number, u(t), uniformly on [0, s(t)]. We then exhaustively score alignments with all alternate placements of the chosen special character within the chosen row. We choose the putative alignment, A(t+1), uniformly among those which score at least as large as u(t).

If a row, i, has m(i) special characters, at most m(i) movements of individual special characters are needed to permute one arrangement of the row into any other. Thus, any alignment may be permuted into any other by at most m(1)+...+m(n) operations, where n is the number of rows. Thus the graph induced by such operations has diameter at most m(1)+...+m(n). Since this is small relative to the number of such arrangments (the product of L choose m(i) where L is the length of a row,) there is a good chance the procedure will mix quickly.

=#

"""
Single letter codes for the amino acids in DNA.
"""    
dna = [Symbol(i) for i in split("A C T G")]

"""
Single letter codes for the amino acids in RNA.
"""    
rna = [Symbol(i) for i in split("A C U G")]

"""
Single letter codes for the amino acids in proteins. Codes for ambiguous reads are omitted.
"""    
proteins = [Symbol(i) for i in split("G P A V L I M C F Y W H K R Q N E D S T")]

function random_alignment_data(alphabet::Array{Symbol,1}, sequence_length::Int, nr_sequences::Int, deletion_probability::Real; special_symbol=:_)
    A = fill(:_, (nr_sequences, sequence_length))
    # base sequence
    A[1,:] = alphabet[rand(1:length(alphabet), sequence_length)]
    # subsequences
    d = Distributions.Bernoulli(1.0-deletion_probability)
    for i in 2:nr_sequences
        subseq = A[1, find(rand(d, sequence_length))]
        A[i,1:length(subseq)] = subseq
    end
    return A
end
