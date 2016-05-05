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



""" hitrun2chains(record::Array{Int,3}, burnin::Int, thin::Int; name::AbstractString="T")
    
Convert an I,J,K array of IxJ contingency tables to an object, C, of type Mamba.Chains in which C.value[K,:,1] contains the entries of table K, reshaped to a vector.
"""
function hitrun2chains{U <: Real}(record::Array{U,3}, burnin::Int, thin::Int; name::AbstractString="T")
    I,J,N = size(record);
    return Mamba.Chains(reshape(record, (I*J, N))',
                        start=burnin+thin,
                        thin=thin,
                        names=make_names(name,(I,J))
                        );
end

"""
If `record` is IxK, convert to an object, C,  of type Mamba.Chains where C.value is KxIx1.     
"""    
function  hitrun2chains{U <: Real}(record::Array{U,2}, burnin::Int, thin::Int; name::AbstractString="T")
    I, K = size(record);
    return Mamba.Chains(record',
                        start=burnin+thin,
                        thin=thin,
                        names=make_names(name,(I,K)));
end

"""
If `record` is a vector of length K, convert to an object, C, of type Mamba.Chains where C.value is Kx1x1.
"""    
function  hitrun2chains{U <: Real}(record::Array{U,1}, burnin::Int, thin::Int; name::AbstractString="T")
    K = size(record,1);
    return Mamba.Chains(reshape(record, (K,1)),
                        start=burnin+thin,
                        thin=thin,
                        names=[name]);
end                        
