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

function random_alignment_data(alphabet::Array{Symbol,1}, sequence_length::Int, nr_sequences::Int, deletion_probability::Real; special_symbol=:*)
    A = fill(special_symbol, (nr_sequences, sequence_length))
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

function count_symbols(A::Array{Symbol,2}, alphabet::Array{Symbol,1})
    n = size(A,1)
    M = zeros(Int,(length(alphabet)+1, size(A,2)))
    for i in 1:length(alphabet)
        M[i,:] = sum(A .== alphabet[i], 1)
    end
    M[end,:] = size(A,1)-sum(M,1)
    return M
end

function log_score_alignment(A::Array{Symbol, 2}, alphabet::Array{Symbol,1}; variance::Float64=1.0)
    counts = count_symbols(A, alphabet)
    return -sum(counts.^2)/(2*variance)
end

function shift_pos!{T<:Int}(A::Array{Symbol,2}, fixedrow::T, altrow::Array{Symbol,1} , oldpos::Int, newpos::Int)
    for i in 1:size(A,2)
        if i == newpos
            altrow[i] = A[fixedrow,oldpos]
        elseif oldpos <= i < newpos
            altrow[i] = A[fixedrow, i+1]
        elseif newpos < i <= oldpos
            altrow[i] = A[fixedrow, i-1]
        else # i < min(oldpos, newpos) or i > max(oldpos, newpos)
            altrow[i] = A[fixedrow, i]
        end
    end
end

"""
Will alter A in place and return log score.    
"""
function step_alignment!(current_log_score::Float64, A::Array{Symbol, 2}, alphabet::Array{Symbol,1}; special_symbol::Symbol=:*, variance::Float64=1.0)
    nrows, ncols = size(A)
    nsym = length(alphabet)
    # Pick a row (other than the first) at random
    rw = rand(2:nrows)
    # Pick a special symbol within the row at random. If there are none, return
    ss = find(special_symbol .== A[rw,:])
    if length(ss)==0
        return current_log_score
    end
    # Generate the log of uniform random variate on [0,1] and add the current
    # log score. The result is equivalent to the log of a uniform random variate
    # on [0, exp(current_log_score)]. Multiply by twice the variance for
    # convenience.
    logu = (log(rand(1))+current_log_score)*2*variance
    # Choose a special symbol at random in the chosen row.
    oldpos = rand(ss)
    # Score all distinct re-positionings of this symbol within the row
    #  For convenience, determine symbol counts exclusive of chosen row
    counts = count_symbols(A[ setdiff(1:nrows, rw), :], alphabet)
    base_score = -sum(counts.^2)
    #  Score all placements of the chosen symbol. Note that placements
    #  within a given block of special symbols all result in the same
    #  sequence. Thus, except for position 1, it suffices to place
    #  the chosen symbol at positions which follow a symbol in the
    #  alphabet.
    positions = A[rw,1] == special_symbol ? vcat(1, find(A[rw,2:ncols] .!= special_symbol)) : find(A[rw,2:ncols] .!= special_symbol)
    scores = fill(base_score,length(positions))
    altrow = fill(special_symbol, ncols)
    for i in 1:length(positions)
        newpos = positions[i]
        shift_pos!(A, rw, altrow, oldpos, newpos) # fills altrow in place
        # Compute score. The base score is the summed squares of the counts
        # of distinct symbols in each column, exclusive of those in the row which
        # is currently being altered. The altered row would add 1 to the count,
        # counts[j, k], of the kth symbol in each column, where k is the index
        # of the symbol in the altered row. Since
        #        -(counts[j, k]+1)^2 = -counts[j,k]^2 - 2*counts[j,k] - 1
        # and since counts[j,k]^2 is in the base score, the altered row
        # contributes increments, -1-2*counts[j,k], to the base score. 
        for j in 1:ncols
            k = altrow[j] == special_symbol ? nsym : find(altrow[j] .== alphabet)[1]
            scores[i] -= 1 + 2*counts[k,j]
        end
    end
    # Choose an altered row at random among those which score as well as logu
    if any(scores .>= logu) == false
        println(scores')
        println(logu)
    end
    k = rand(find(scores .>= logu))
    newpos = positions[k]
    newsc  = scores[k]
    # Generate the altered row
    shift_pos!(A, rw, altrow, oldpos, newpos)
    # Alter A in place accordingly
    for i in 1:ncols
        A[rw,i]=altrow[i]
    end
    # Return its score
    return newsc/(2*variance)
end
