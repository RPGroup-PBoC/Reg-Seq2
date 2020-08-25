using Test, DataFrames, Pandas, CSV, jregseq


#df = Pandas.read_csv("../data/ykgEAnaerodataset_alldone_with_large", delim_whitespace=true)
#df = DataFrames.DataFrame(df)
#ndf = select_region(df, 30, 50) 

@testset begin
    @test_broken length(DataFrames.groupby(ndf, :seq)) == nrow(ndf)
    @test_broken ~any((ndf[!, :ct_0] .+ ndf[!, :ct_1]) .!= ndf[!, :ct])
    @test_broken ~any(map(x -> length(x) != 21, ndf[!, :seq]))
end

@testset begin
   @test_broken transform_seq(transform_seq(ndf.seq[1])) == ndf.seq[1]
   @test_broken transform_seq(transform_seq(ndf.seq[1])) == ndf.seq[1]
   #@test_throws_broken ArgumentError transform_seq(["X"])
   #@test_throws ArgumentError transform_seq(["5"])
   @test_broken typeof(transform_seq(ndf.seq[1])) == Array{Int64,1}
end

df = select!(CSV.read("../../data/RegSeq/wtsequences.csv"), Not(:Column1))
sequence = df[df.name.=="ykgE", :geneseq][1]

@testset "scrambles" begin
    @test_throws ArgumentError jregseq.scramble(sequence, 10, 3)
end
