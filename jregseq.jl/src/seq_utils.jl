import Random.shuffle
Random.shuffle(s::String) = isascii(s) ? s[randperm(end)] : join(shuffle!(collect(s)))


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


"""
    scramble(sequence, windowsize::Int, overlap::Int)
"""
function scramble(sequence, windowsize::Int, overlap::Int)
    l_sequence = length(sequence)
    n_scrambles = (l_sequence - windowsize) / (windowsize - overlap)
    
    # Test if scrambles can created throughout whole site
    if ~isinteger(n_scrambles)
        throw(ArgumentError("Cannot make scrambles with windowsize $windowsize and overlap $overlap."))
    end
    scrambled_sequences = [sequence]
    indices(i) = (windowsize-overlap)*(i-1)+1 : (windowsize-overlap)*(i-1) + windowsize
    for i in 1:Int(n_scrambles)+1
        temp_sequence = collect(sequence)
        temp_sequence[indices(i)] = shuffle(temp_sequence[indices(i)])
        push!(scrambled_sequences, string(temp_sequence...))
    end
    return scrambled_sequences
end
    

function site_single_mutations(sequence, site_start::Int, site_end::Int)
    scrambled_sequences = [sequence]
    letters = ["A", "C", "G", "T"]
    if site_end < site_start
        throw(ArgumentError("The end of the site should be larger than the start."))
    end
    for i in site_start:site_end
        mutations = filter(x-> x[1] != sequence[i], letters)
        for x in mutations
            temp_sequence = collect(sequence)
            temp_sequence[i] = x[1]
            push!(scrambled_sequences, string(temp_sequence...))
        end
    end
    return scrambled_sequences
end
