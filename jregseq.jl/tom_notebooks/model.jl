using DataFrames, CmdStan, RecursiveArrayTools, JLD, LinearAlgebra
import Pandas.read_csv

df = read_csv("data/ykgEAnaerodataset_alldone_with_large", delim_whitespace=true)
df = DataFrames.DataFrame(df)


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


s = transform_seq.(df.seq)
VA = VectorOfArray(s)
arr = transpose(convert(Array,VA))


stanfile = open("stan_files/model.stan") do file
    read(file, String)
end
stanmodel = Stanmodel(
    model=stanfile,
    name="model",
    nchains=2,
    thin=10,
    num_warmup=2000,
    num_samples=1000,
    #Sample(algorithm=CmdStan.Fixed_param())

)

data = Dict(
    "L_S" => length(df.seq[1]), 
    "N_S" => nrow(df), 
    "seqs" => arr, 
    "ct_0"=>convert(Array{Int64}, df.ct_0),
    "ct_1"=>convert(Array{Int64}, df.ct_1),
    "n"=>sum(convert(Array{Int64}, df.ct)),
)

a, chains, b = stan(stanmodel, data, summary=false)
chains = JLD.load("chains.jld")["chains"]
b = JLD.load("b.jld")["b"]