#=
NOTES:

Mamba convergence diagnostics operatate on Mamba.Chains objects. A Mamba.Chains object is immutable, hence must be constructed after simulation. It has 4 fields as follows.

* :chains -- a tuple, (n,), where n is the number of chains
* :range  -- a range object: burnin+thin:thin:iterations
* :names  -- an array of variable names
* :values -- a 3D array of size (|T|, length(range), chains)
where |T| would be the number of elements in a contingengy table and :names would be e.g., `["T[1,1]", "T[2,1]", ..."T[m,n]"]`.

=#

"""  make_names(name::AbstractString, dim::Tuple{Int,Vararg{Int}})
Return a vector of names suitable for the `names` field of `Mamba.Chains`. E.g., make_names("T",(2,3)) should return the vector, `["T[1,1], "T[2,1]", "T[1,2]", "T[2,2]", "T[1,3]", "T[2,3]"]`.
"""
function make_names(name::AbstractString, dim::Tuple{Int,Vararg{Int}})
    [name*"["*i for i in make_indices(dim)]
end

"""  make_indices(dim::Tuple{Int,Vararg{Int}})
Return a vector of strings of a form `i]`, `i,j]`, etc. reflecting indices of an array whose dimensions are given in the argument. 
"""
function make_indices(dim::Tuple{Int,Vararg{Int}})
    n = length(dim);
    if n == 1
        return [string(k)*"]" for k in 1:dim[end]];
    end
    tail = make_indices(dim[2:n]);
    head = [string(k)*"," for k in 1:dim[1]];
    m = length(tail)*length(head);
    return reshape([h*t for h in head, t in tail], m);
end
