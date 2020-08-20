using Pandas, DataFrames, CmdStan, RecursiveArrayTools, JLD, LinearAlgebra


df = Pandas.read_csv("data/ykgEAnaerodataset_alldone_with_large", delim_whitespace=true)
df = DataFrames.DataFrame(df)

function select_region(df::DataFrames.DataFrame, a, b)
    # Sort for sequences with mutations in region
    df = df[map(x-> x[a:b], df.seq) .!= df.seq[1][a:b], :]
    
    # Take selected region
    df.seq = map(x-> x[a:b], df.seq)
    
    # Join equal sequences
    gdf = DataFrames.groupby(df, :seq)
    df = combine(gdf, :ct => sum, :ct_0 => sum, :ct_1 => sum)
    
    # Rename colums
    rename!(df, :ct_sum => :ct, :ct_0_sum => :ct_0, :ct_1_sum => :ct_1)
    return df
end


function transform_seq(seq)
    nuc_dic = Dict('A'=>'1', 'C'=>'2', 'G'=>'3', 'T'=>'4')
    revnuc_dic = Dict('1'=>'A', '2'=>'C', '3'=>'G', '4'=>'T')
    
    if typeof(seq) == String
        if seq[1] in ['A', 'C', 'G', 'T']
            seq = map(x -> nuc_dic[x], seq)
            seq = parse.(Int, collect(seq))
        else 
            throw(ArgumentError("Sequence of unknown type."))
        end
    elseif typeof(seq) == Array{Int64,1}
        if seq[1] in [1, 2, 3, 4]
            seq = String(map(x -> revnuc_dic[Char('0' + x)], seq))
        else
            throw(ArgumentError("Sequence of unknown type."))
        end
    else 
        throw(ArgumentError("Sequence of unknown type."))
    end
    return seq
end

ndf = select_region(df, 30, 40)

stanfile = open("stan_files/model.stan") do file
    read(file, String)
end

stanmodel = Stanmodel(
    model=stanfile,
    name="model",
    nchains=4,
    thin=10,
    num_warmup=1000,
    num_samples=10000,
    #Sample(algorithm=CmdStan.Fixed_param())

)

data = Dict(
    "L_S" => length(ndf.seq[1]), 
    "N_S" => nrow(ndf), 
    "seqs" => arr_real, 
    "ct_0"=>convert(Array{Int64}, ndf.ct_0),
    "ct_1"=>convert(Array{Int64}, ndf.ct_1),
    "ct"=>convert(Array{Int64}, ndf.ct),
    "n"=>sum(convert(Array{Int64}, ndf.ct)),
    "n_1"=>sum(convert(Array{Int64}, ndf.ct_1)),
    "n_0"=>sum(convert(Array{Int64}, ndf.ct_0)),
    "theta_length"=>4*length(ndf.seq[1])
)

a, chains, b = stan(stanmodel, data, summary=false);

JLD.save("rc.jld", "data", a)
JLD.save("chains.jld", "data", chains)
JLD.save("names.jld", "data", b)
