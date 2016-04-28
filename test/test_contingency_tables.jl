# Is make_names working?
Base.Test.@test HitAndRun.make_names("A",(2,)) == ["A[1]","A[2]"];
Base.Test.@test HitAndRun.make_names("A",(2,2)) == ["A[1,1]", "A[2,1]", "A[1,2]", "A[2,2]"];

# initial contingency table
snee = HitAndRun.snee;
# mcmc samples starting from snee, using default option, uniform
record = HitAndRun.mcmc_tables(snee, 500, 50, 50000);
# different representation of same samples
chain = HitAndRun.tables2chains(record,500,50);

I,J,N = size(record);
N1, IJ, C = size(chain.value);
Base.Test.@test N == N1;   # same number of samples?
Base.Test.@test I*J == IJ; # same number of parameters?
Base.Test.@test C == 1;    # one chain?

n = rand(1:N); # random place in the chain
tbl1 = record[:,:,n];   # 4x4 table
tbl2 = chain.value[n,:]; # 16 long vector
Base.Test.@test tbl1 == reshape(tbl2,(4,4)); # equal?

Base.Test.@test sum(tbl1,1) == sum(snee,1); # same col sums?
Base.Test.@test sum(tbl1,2) == sum(snee,2); # same row sums?
