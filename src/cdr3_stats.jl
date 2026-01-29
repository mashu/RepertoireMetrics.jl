# CDR3 length statistics
# These operate on DataFrames since CDR3 length is a sequence-level property

using Statistics: mean, std, median

"""
    CDR3Stats

Container for CDR3 length statistics.

# Fields
- `mean_length::Float64`: Mean CDR3 length
- `median_length::Float64`: Median CDR3 length
- `std_length::Float64`: Standard deviation of CDR3 length
- `min_length::Int`: Minimum CDR3 length
- `max_length::Int`: Maximum CDR3 length
- `n_sequences::Int`: Number of sequences analyzed
- `weighted::Bool`: Whether statistics were weighted by count
"""
struct CDR3Stats
    mean_length::Float64
    median_length::Float64
    std_length::Float64
    min_length::Int
    max_length::Int
    n_sequences::Int
    weighted::Bool
end

function Base.show(io::IO, s::CDR3Stats)
    w = s.weighted ? " (count-weighted)" : ""
    println(io, "CDR3Stats$w:")
    println(io, "  Mean length:   ", round(s.mean_length, digits=2))
    println(io, "  Median length: ", round(s.median_length, digits=2))
    println(io, "  Std length:    ", round(s.std_length, digits=2))
    println(io, "  Range:         ", s.min_length, " - ", s.max_length)
    print(io, "  N sequences:   ", s.n_sequences)
end

"""
    cdr3_stats(df::DataFrame; cdr3_column=:cdr3, count_column=nothing, use_aa=false) -> CDR3Stats

Compute CDR3 length statistics from a DataFrame.

# Arguments
- `df`: DataFrame containing CDR3 sequences
- `cdr3_column`: Column name containing CDR3 sequences (default `:cdr3`)
- `count_column`: If provided, weight statistics by this count column. If `nothing`, 
  each row is weighted equally.
- `use_aa`: If `true`, expects amino acid sequences. If `false` (default), expects 
  nucleotide sequences and divides length by 3 for amino acid length.

# Returns
A `CDR3Stats` object with mean, median, std, min, max lengths and sequence count.

# Example
```julia
using RepertoireMetrics, CSV, DataFrames

df = CSV.read("sequences.tsv", DataFrame)

# Unweighted statistics
stats = cdr3_stats(df)
println("Mean CDR3 length: ", stats.mean_length)

# Weighted by count column
stats = cdr3_stats(df; count_column=:count)

# Using amino acid sequences directly
stats = cdr3_stats(df; cdr3_column=:cdr3_aa, use_aa=true)
```
"""
function cdr3_stats(
    df::DataFrame;
    cdr3_column::Symbol = :cdr3,
    count_column::Union{Symbol,Nothing} = nothing,
    use_aa::Bool = false
)
    hasproperty(df, cdr3_column) || 
        throw(ArgumentError("Column $cdr3_column not found in DataFrame"))
    
    # Extract CDR3 lengths
    lengths = Int[]
    weights = Int[]
    
    has_count = count_column !== nothing && hasproperty(df, count_column)
    
    for row in eachrow(df)
        cdr3 = row[cdr3_column]
        
        # Skip missing/empty
        if cdr3 === missing || cdr3 === nothing || isempty(string(cdr3))
            continue
        end
        
        len = length(string(cdr3))
        if !use_aa
            len = len ÷ 3  # Convert nucleotide to amino acid length
        end
        
        push!(lengths, len)
        
        if has_count
            c = row[count_column]
            push!(weights, (c === missing || c === nothing) ? 1 : round(Int, c))
        else
            push!(weights, 1)
        end
    end
    
    isempty(lengths) && return CDR3Stats(0.0, 0.0, 0.0, 0, 0, 0, has_count)
    
    # Compute statistics
    if has_count && sum(weights) > 0
        # Weighted statistics
        total_weight = sum(weights)
        weighted_mean = sum(lengths .* weights) / total_weight
        
        # Weighted standard deviation
        weighted_var = sum(weights .* (lengths .- weighted_mean).^2) / total_weight
        weighted_std = sqrt(weighted_var)
        
        # Weighted median (expand and compute)
        expanded = reduce(vcat, [fill(l, w) for (l, w) in zip(lengths, weights)])
        weighted_median = Float64(median(expanded))
        
        return CDR3Stats(
            weighted_mean,
            weighted_median,
            weighted_std,
            minimum(lengths),
            maximum(lengths),
            total_weight,
            true
        )
    else
        # Unweighted statistics
        return CDR3Stats(
            mean(lengths),
            Float64(median(lengths)),
            length(lengths) > 1 ? std(lengths) : 0.0,
            minimum(lengths),
            maximum(lengths),
            length(lengths),
            false
        )
    end
end

"""
    cdr3_stats(filepath::AbstractString, strategy::AbstractLineageDefinition; kwargs...) -> CDR3Stats

Convenience method to compute CDR3 statistics directly from a file.

# Example
```julia
stats = cdr3_stats("sequences.tsv", VJCdr3Definition(); count_column=:count)
```
"""
function cdr3_stats(
    filepath::AbstractString,
    strategy::AbstractLineageDefinition;
    cdr3_column::Symbol = :cdr3,
    count_column::Union{Symbol,Nothing} = :count,
    use_aa::Bool = false,
    kwargs...
)
    csv_kwargs = Dict{Symbol,Any}(
        :missingstring => ["", "NA", "na", "N/A", "n/a"],
        :delim => '\t',
    )
    merge!(csv_kwargs, Dict(kwargs))
    
    df = CSV.read(filepath, DataFrame; csv_kwargs...)
    
    # Use cdr3_column from strategy if it's VJCdr3Definition
    if strategy isa VJCdr3Definition
        cdr3_column = strategy.cdr3_column
    end
    
    return cdr3_stats(df; cdr3_column=cdr3_column, count_column=count_column, use_aa=use_aa)
end

"""
    cdr3_length_distribution(df::DataFrame; cdr3_column=:cdr3, count_column=nothing, use_aa=false) -> Dict{Int,Int}

Compute the distribution of CDR3 lengths.

# Returns
A `Dict{Int,Int}` mapping CDR3 length to count.

# Example
```julia
dist = cdr3_length_distribution(df)
for (length, cnt) in sort(collect(dist))
    println("Length \$length: \$cnt sequences")
end
```
"""
function cdr3_length_distribution(
    df::DataFrame;
    cdr3_column::Symbol = :cdr3,
    count_column::Union{Symbol,Nothing} = nothing,
    use_aa::Bool = false
)
    hasproperty(df, cdr3_column) || 
        throw(ArgumentError("Column $cdr3_column not found in DataFrame"))
    
    distribution = Dict{Int,Int}()
    has_count = count_column !== nothing && hasproperty(df, count_column)
    
    for row in eachrow(df)
        cdr3 = row[cdr3_column]
        
        if cdr3 === missing || cdr3 === nothing || isempty(string(cdr3))
            continue
        end
        
        len = length(string(cdr3))
        if !use_aa
            len = len ÷ 3
        end
        
        weight = 1
        if has_count
            c = row[count_column]
            weight = (c === missing || c === nothing) ? 1 : round(Int, c)
        end
        
        distribution[len] = get(distribution, len, 0) + weight
    end
    
    return distribution
end
