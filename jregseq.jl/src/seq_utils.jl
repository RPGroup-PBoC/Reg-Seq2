"""
    function select_region(df::DataFrames.DataFrame, a, b)

Select region in sequence and find unique regions. Duplicates will be joined.

#Arguments
----------
- `df` : DataFrames.Dataframe
    Dataframe containing sequences `seq`, DNA counts `ct_0`, RNA_counts `ct_1` and total counts `ct`
- `a` : first base
- `b` : last base
"""
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


"""
    transform_seq!(seq)

Transforms nucleotides into integers or vice versa.
"""
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