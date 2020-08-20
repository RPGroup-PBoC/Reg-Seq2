using Test, DataFrames, Pandas

include("../src/seq_utils.jl")

df = Pandas.read_csv("../data/ykgEAnaerodataset_alldone_with_large", delim_whitespace=true)
df = DataFrames.DataFrame(df)
ndf = select_region(df, 30, 50) 

@testset begin
    @test length(DataFrames.groupby(ndf, :seq)) == nrow(ndf)
    @test ~any((ndf[!, :ct_0] .+ ndf[!, :ct_1]) .!= ndf[!, :ct])
    @test ~any(map(x -> length(x) != 21, ndf[!, :seq]))
end

@testset begin
   @test transform_seq(transform_seq(ndf.seq[1])) == ndf.seq[1]
   @test transform_seq(transform_seq(ndf.seq[1])) == ndf.seq[1]
   @test_throws ArgumentError transform_seq(["X"])
   @test_throws ArgumentError transform_seq(["5"])
   @test typeof(transform_seq(ndf.seq[1])) == Array{Int64,1}
end