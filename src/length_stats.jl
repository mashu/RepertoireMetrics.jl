# Length statistics - generic sequence length metrics
# Can be used for CDR3, junction, or any string column

using Statistics: mean, std, median

"""
    LengthStats

Container for sequence length statistics, stored in Repertoire metadata.

# Fields
- `mean_length::Float64`: Mean sequence length
- `median_length::Float64`: Median sequence length  
- `std_length::Float64`: Standard deviation of sequence length
- `min_length::Int`: Minimum sequence length
- `max_length::Int`: Maximum sequence length
- `n_sequences::Int`: Number of sequences analyzed
- `column::Symbol`: Column used for length calculation
- `use_aa::Bool`: Whether lengths are in amino acids

Computed during repertoire construction when a `length_column` is specified.
"""
struct LengthStats
    mean_length::Float64
    median_length::Float64
    std_length::Float64
    min_length::Int
    max_length::Int
    n_sequences::Int
    column::Symbol
    use_aa::Bool
end

function Base.show(io::IO, s::LengthStats)
    unit = s.use_aa ? "aa" : "nt→aa"
    println(io, "LengthStats ($(s.column), $unit):")
    println(io, "  Mean:   ", round(s.mean_length, digits=2))
    println(io, "  Median: ", round(s.median_length, digits=2))
    println(io, "  Std:    ", round(s.std_length, digits=2))
    println(io, "  Range:  ", s.min_length, " - ", s.max_length)
    print(io, "  N:      ", s.n_sequences)
end

"""
    compute_length_stats(df::DataFrame; length_column=:cdr3, count_column=nothing, use_aa=false) -> LengthStats

Compute length statistics from a DataFrame column.

# Arguments
- `df`: DataFrame containing sequences
- `length_column`: Column to measure length from (default `:cdr3`)
- `count_column`: If provided, weight statistics by this count column
- `use_aa`: If `true`, column contains amino acid sequences. If `false` (default),
  assumes nucleotide sequences and divides length by 3.

# Returns
A `LengthStats` object with computed statistics.
"""
function compute_length_stats(
    df::DataFrame;
    length_column::Symbol = :cdr3,
    count_column::Union{Symbol,Nothing} = nothing,
    use_aa::Bool = false
)
    hasproperty(df, length_column) || 
        return LengthStats(0.0, 0.0, 0.0, 0, 0, 0, length_column, use_aa)
    
    lengths = Int[]
    weights = Int[]
    
    has_count = count_column !== nothing && hasproperty(df, count_column)
    
    for row in eachrow(df)
        val = row[length_column]
        
        # Skip missing/empty
        if val === missing || val === nothing
            continue
        end
        
        s = string(val)
        if isempty(s)
            continue
        end
        
        len = length(s)
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
    
    isempty(lengths) && return LengthStats(0.0, 0.0, 0.0, 0, 0, 0, length_column, use_aa)
    
    total_weight = sum(weights)
    
    if has_count && total_weight > 0
        # Weighted statistics
        weighted_mean = sum(lengths .* weights) / total_weight
        weighted_var = sum(weights .* (lengths .- weighted_mean).^2) / total_weight
        weighted_std = sqrt(weighted_var)
        
        # Weighted median
        expanded = reduce(vcat, [fill(l, w) for (l, w) in zip(lengths, weights)])
        weighted_median = Float64(median(expanded))
        
        return LengthStats(
            weighted_mean,
            weighted_median,
            weighted_std,
            minimum(lengths),
            maximum(lengths),
            total_weight,
            length_column,
            use_aa
        )
    else
        return LengthStats(
            mean(lengths),
            Float64(median(lengths)),
            length(lengths) > 1 ? std(lengths) : 0.0,
            minimum(lengths),
            maximum(lengths),
            length(lengths),
            length_column,
            use_aa
        )
    end
end

"""
    length_distribution(df::DataFrame; length_column=:cdr3, count_column=nothing, use_aa=false) -> Dict{Int,Int}

Compute the distribution of sequence lengths.

# Arguments
- `df`: DataFrame containing sequences
- `length_column`: Column to measure length from (default `:cdr3`)
- `count_column`: If provided, weight by this count column
- `use_aa`: If `true`, column contains amino acid sequences

# Returns
A `Dict{Int,Int}` mapping length to count.

# Example
```julia
dist = length_distribution(df; length_column=:cdr3)
for (len, cnt) in sort(collect(dist))
    println("Length \$len: \$cnt sequences")
end
```
"""
function length_distribution(
    df::DataFrame;
    length_column::Symbol = :cdr3,
    count_column::Union{Symbol,Nothing} = nothing,
    use_aa::Bool = false
)
    hasproperty(df, length_column) || 
        throw(ArgumentError("Column $length_column not found in DataFrame"))
    
    distribution = Dict{Int,Int}()
    has_count = count_column !== nothing && hasproperty(df, count_column)
    
    for row in eachrow(df)
        val = row[length_column]
        
        if val === missing || val === nothing
            continue
        end
        
        s = string(val)
        if isempty(s)
            continue
        end
        
        len = length(s)
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

# ============================================================================
# Composable length metric types
# ============================================================================

"""
    MeanLength <: AbstractMetric

Metric type for mean sequence length. Requires `length_column` to be specified
during repertoire construction.
"""
struct MeanLength <: AbstractMetric end

"""
    MedianLength <: AbstractMetric

Metric type for median sequence length.
"""
struct MedianLength <: AbstractMetric end

"""
    StdLength <: AbstractMetric

Metric type for standard deviation of sequence length.
"""
struct StdLength <: AbstractMetric end

"""
    MinLength <: AbstractMetric

Metric type for minimum sequence length.
"""
struct MinLength <: AbstractMetric end

"""
    MaxLength <: AbstractMetric

Metric type for maximum sequence length.
"""
struct MaxLength <: AbstractMetric end

# Accessor functions for Repertoire
"""
    mean_length(rep::Repertoire) -> Float64

Return the mean sequence length. Requires length statistics to be computed
during repertoire construction (via `length_column` parameter).
"""
function mean_length(rep::Repertoire)
    stats = get(rep.metadata, "length_stats", nothing)
    stats === nothing && error("Length statistics not computed. Use length_column parameter when creating repertoire.")
    return stats.mean_length
end

"""
    median_length(rep::Repertoire) -> Float64

Return the median sequence length.
"""
function median_length(rep::Repertoire)
    stats = get(rep.metadata, "length_stats", nothing)
    stats === nothing && error("Length statistics not computed. Use length_column parameter when creating repertoire.")
    return stats.median_length
end

"""
    std_length(rep::Repertoire) -> Float64

Return the standard deviation of sequence length.
"""
function std_length(rep::Repertoire)
    stats = get(rep.metadata, "length_stats", nothing)
    stats === nothing && error("Length statistics not computed. Use length_column parameter when creating repertoire.")
    return stats.std_length
end

"""
    min_length(rep::Repertoire) -> Int

Return the minimum sequence length.
"""
function min_length(rep::Repertoire)
    stats = get(rep.metadata, "length_stats", nothing)
    stats === nothing && error("Length statistics not computed. Use length_column parameter when creating repertoire.")
    return stats.min_length
end

"""
    max_length(rep::Repertoire) -> Int

Return the maximum sequence length.
"""
function max_length(rep::Repertoire)
    stats = get(rep.metadata, "length_stats", nothing)
    stats === nothing && error("Length statistics not computed. Use length_column parameter when creating repertoire.")
    return stats.max_length
end

"""
    has_length_stats(rep::Repertoire) -> Bool

Check if length statistics are available for this repertoire.
"""
has_length_stats(rep::Repertoire) = haskey(rep.metadata, "length_stats")

"""
    length_stats(rep::Repertoire) -> LengthStats

Return the LengthStats object for this repertoire.
"""
function length_stats(rep::Repertoire)
    stats = get(rep.metadata, "length_stats", nothing)
    stats === nothing && error("Length statistics not computed. Use length_column parameter when creating repertoire.")
    return stats
end

# compute_metric implementations for length metrics
compute_metric(rep::Repertoire, ::MeanLength) = mean_length(rep)
compute_metric(rep::Repertoire, ::MedianLength) = median_length(rep)
compute_metric(rep::Repertoire, ::StdLength) = std_length(rep)
compute_metric(rep::Repertoire, ::MinLength) = Float64(min_length(rep))
compute_metric(rep::Repertoire, ::MaxLength) = Float64(max_length(rep))

# Metric symbols for output (required by compute_metrics)
metric_symbol(::MeanLength) = :mean_length
metric_symbol(::MedianLength) = :median_length
metric_symbol(::StdLength) = :std_length
metric_symbol(::MinLength) = :min_length
metric_symbol(::MaxLength) = :max_length

# Predefined length metric set
"""
    LENGTH_METRICS

Predefined MetricSet containing all length metrics:
`MeanLength() + MedianLength() + StdLength() + MinLength() + MaxLength()`
"""
const LENGTH_METRICS = MeanLength() + MedianLength() + StdLength() + MinLength() + MaxLength()
