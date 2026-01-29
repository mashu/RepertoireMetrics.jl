# Composable metric selection system
# Allows users to select which metrics to compute using + operator

"""
    AbstractMetric

Abstract supertype for individual metric selectors.
Each concrete subtype represents a specific diversity/clonality metric.
"""
abstract type AbstractMetric end

# ============================================================================
# Individual Metric Types
# ============================================================================

"""Richness: Number of unique lineages (S)"""
struct Richness <: AbstractMetric end

"""TotalCount: Total sequence/cell count (N)"""
struct TotalCount <: AbstractMetric end

"""ShannonEntropy: H = -Σ(pᵢ log pᵢ)"""
struct ShannonEntropy <: AbstractMetric end

"""ShannonDiversity: exp(H) - effective number of lineages"""
struct ShannonDiversity <: AbstractMetric end

"""NormalizedShannon: H / log(S), normalized to [0,1]"""
struct NormalizedShannon <: AbstractMetric end

"""SimpsonIndex: D = Σpᵢ² - probability two random sequences are same lineage"""
struct SimpsonIndex <: AbstractMetric end

"""SimpsonDiversity: 1 - D (Gini-Simpson index)"""
struct SimpsonDiversity <: AbstractMetric end

"""InverseSimpson: 1/D - effective number of lineages (Hill q=2)"""
struct InverseSimpson <: AbstractMetric end

"""BergerParker: Proportion of most abundant lineage"""
struct BergerParker <: AbstractMetric end

"""Evenness: Pielou's J = H / log(S)"""
struct Evenness <: AbstractMetric end

"""Clonality: 1 - normalized Shannon entropy"""
struct Clonality <: AbstractMetric end

"""GiniCoefficient: Inequality measure [0,1]"""
struct GiniCoefficient <: AbstractMetric end

"""D50: Minimum lineages comprising 50% of repertoire"""
struct D50 <: AbstractMetric end

"""Chao1: Richness estimator including unobserved species"""
struct Chao1 <: AbstractMetric end

# ============================================================================
# MetricSet - Composable collection of metrics
# ============================================================================

"""
    MetricSet{M<:Tuple}

A composable set of metrics to compute. Use the `+` operator to combine metrics.

# Examples
```julia
# Select specific metrics
metrics = Richness() + ShannonEntropy() + Clonality()
result = compute_metrics(rep, metrics)

# Use predefined sets
result = compute_metrics(rep, DIVERSITY_METRICS)
result = compute_metrics(rep, ALL_METRICS)

# Default computes all metrics
result = compute_metrics(rep)
```
"""
struct MetricSet{M<:Tuple}
    metrics::M
end

MetricSet(m::AbstractMetric) = MetricSet((m,))
MetricSet(ms::AbstractMetric...) = MetricSet(ms)

# Combine metrics with + operator
Base.:+(a::AbstractMetric, b::AbstractMetric) = MetricSet((a, b))
Base.:+(a::MetricSet, b::AbstractMetric) = MetricSet((a.metrics..., b))
Base.:+(a::AbstractMetric, b::MetricSet) = MetricSet((a, b.metrics...))
Base.:+(a::MetricSet, b::MetricSet) = MetricSet((a.metrics..., b.metrics...))

# Iteration support
Base.iterate(ms::MetricSet) = iterate(ms.metrics)
Base.iterate(ms::MetricSet, state) = iterate(ms.metrics, state)
Base.length(ms::MetricSet) = length(ms.metrics)

function Base.show(io::IO, ms::MetricSet)
    print(io, "MetricSet(")
    for (i, m) in enumerate(ms.metrics)
        i > 1 && print(io, " + ")
        print(io, nameof(typeof(m)))
    end
    print(io, ")")
end

# ============================================================================
# Predefined MetricSets
# ============================================================================

"""
    ALL_METRICS

All available metrics. This is the default when calling `compute_metrics(rep)`.
"""
const ALL_METRICS = MetricSet((
    Richness(),
    TotalCount(),
    ShannonEntropy(),
    ShannonDiversity(),
    NormalizedShannon(),
    SimpsonIndex(),
    SimpsonDiversity(),
    InverseSimpson(),
    BergerParker(),
    Evenness(),
    Clonality(),
    GiniCoefficient(),
    D50(),
    Chao1()
))

"""
    DIVERSITY_METRICS

Common diversity metrics: Shannon entropy/diversity, Simpson diversity, inverse Simpson.
"""
const DIVERSITY_METRICS = MetricSet((
    Richness(),
    TotalCount(),
    ShannonEntropy(),
    ShannonDiversity(),
    SimpsonDiversity(),
    InverseSimpson(),
    Evenness()
))

"""
    CLONALITY_METRICS

Metrics focused on clonal expansion: clonality, Gini, Berger-Parker, D50.
"""
const CLONALITY_METRICS = MetricSet((
    Richness(),
    TotalCount(),
    Clonality(),
    GiniCoefficient(),
    BergerParker(),
    D50()
))

"""
    RICHNESS_METRICS

Richness-related metrics: observed richness and Chao1 estimator.
"""
const RICHNESS_METRICS = MetricSet((
    Richness(),
    TotalCount(),
    Chao1()
))

# ============================================================================
# Metrics result container
# ============================================================================

"""
    Metrics <: AbstractMetricResult

Result container for computed metrics. Access values by property name.
All values are stored as `Float64` for type stability.

# Example
```julia
result = compute_metrics(rep)
println(result.shannon_entropy)
println(result.clonality)
println(result.d50)
```
"""
struct Metrics <: AbstractMetricResult
    values::Dict{Symbol, Float64}
    metric_names::Vector{Symbol}
end

function Base.getproperty(m::Metrics, name::Symbol)::Union{Dict{Symbol,Float64}, Vector{Symbol}, Float64, Missing}
    if name === :values
        return getfield(m, :values)
    elseif name === :metric_names
        return getfield(m, :metric_names)
    end
    return get(getfield(m, :values), name, missing)
end

Base.propertynames(m::Metrics) = getfield(m, :metric_names)

function Base.show(io::IO, m::Metrics)
    println(io, "Metrics ($(length(m.metric_names)) computed):")
    for name in m.metric_names
        val = m.values[name]
        if val isa AbstractFloat
            println(io, "  $name: ", round(val, digits=4))
        else
            println(io, "  $name: ", val)
        end
    end
end

# Symbol name for each metric type
metric_symbol(::Richness) = :richness
metric_symbol(::TotalCount) = :total_count
metric_symbol(::ShannonEntropy) = :shannon_entropy
metric_symbol(::ShannonDiversity) = :shannon_diversity
metric_symbol(::NormalizedShannon) = :normalized_shannon
metric_symbol(::SimpsonIndex) = :simpson_index
metric_symbol(::SimpsonDiversity) = :simpson_diversity
metric_symbol(::InverseSimpson) = :inverse_simpson
metric_symbol(::BergerParker) = :berger_parker
metric_symbol(::Evenness) = :evenness
metric_symbol(::Clonality) = :clonality
metric_symbol(::GiniCoefficient) = :gini_coefficient
metric_symbol(::D50) = :d50
metric_symbol(::Chao1) = :chao1

# ============================================================================
# Individual metric computation (dispatch on metric type)
# ============================================================================

"""
    compute_metric(rep::Repertoire, metric::AbstractMetric) -> value

Compute a single metric for a repertoire.
"""
function compute_metric end

compute_metric(rep::Repertoire, ::Richness) = richness(rep)
compute_metric(rep::Repertoire, ::TotalCount) = total_count(rep)
compute_metric(rep::Repertoire, ::ShannonEntropy) = shannon_entropy(rep)
compute_metric(rep::Repertoire, ::ShannonDiversity) = exp(shannon_entropy(rep))
compute_metric(rep::Repertoire, ::SimpsonIndex) = simpson_index(rep)
compute_metric(rep::Repertoire, ::SimpsonDiversity) = simpson_diversity(rep)
compute_metric(rep::Repertoire, ::InverseSimpson) = inverse_simpson(rep)
compute_metric(rep::Repertoire, ::BergerParker) = berger_parker_index(rep)
compute_metric(rep::Repertoire, ::Evenness) = evenness(rep)
compute_metric(rep::Repertoire, ::Clonality) = clonality(rep)
compute_metric(rep::Repertoire, ::GiniCoefficient) = gini_coefficient(rep)
compute_metric(rep::Repertoire, ::D50) = d50(rep)
compute_metric(rep::Repertoire, ::Chao1) = chao1(rep)

function compute_metric(rep::Repertoire, ::NormalizedShannon)
    S = richness(rep)
    S <= 1 && return 0.0
    H = shannon_entropy(rep)
    return H / log(S)
end

# ============================================================================
# Main compute_metrics functions
# ============================================================================

"""
    compute_metrics(rep::Repertoire) -> Metrics
    compute_metrics(rep::Repertoire, metrics::MetricSet) -> Metrics

Compute metrics for a repertoire. Without a MetricSet argument, computes all metrics.

# Examples
```julia
# Compute all metrics (default)
result = compute_metrics(rep)

# Compute specific metrics
result = compute_metrics(rep, ShannonEntropy() + Clonality() + D50())

# Use predefined sets
result = compute_metrics(rep, DIVERSITY_METRICS)
```
"""
function compute_metrics(rep::Repertoire, ms::MetricSet)
    values = Dict{Symbol, Float64}()
    names = Symbol[]
    
    for m in ms.metrics
        name = metric_symbol(m)
        values[name] = Float64(compute_metric(rep, m))
        push!(names, name)
    end
    
    return Metrics(values, names)
end

# Default: compute all metrics
compute_metrics(rep::Repertoire) = compute_metrics(rep, ALL_METRICS)

# Single metric shortcut
compute_metrics(rep::Repertoire, m::AbstractMetric) = compute_metrics(rep, MetricSet(m))

"""
    compute_metrics(collection::RepertoireCollection) -> Vector{Metrics}
    compute_metrics(collection::RepertoireCollection, metrics::MetricSet) -> Vector{Metrics}

Compute metrics for all repertoires in a collection.
"""
function compute_metrics(collection::RepertoireCollection, ms::MetricSet)
    return [compute_metrics(rep, ms) for rep in collection]
end

compute_metrics(collection::RepertoireCollection) = compute_metrics(collection, ALL_METRICS)

compute_metrics(collection::RepertoireCollection, m::AbstractMetric) = 
    compute_metrics(collection, MetricSet(m))

# ============================================================================
# DataFrame conversion for Metrics
# ============================================================================

"""
    metrics_to_dataframe(m::Metrics, donor_id::String="") -> DataFrame

Convert Metrics to a single-row DataFrame.
"""
function metrics_to_dataframe(m::Metrics, donor_id::String="")
    df = DataFrame(donor_id = donor_id)
    for name in m.metric_names
        df[!, name] = [m.values[name]]
    end
    return df
end

"""
    metrics_to_dataframe(collection::RepertoireCollection, metrics::Vector{Metrics}) -> DataFrame

Convert a collection's metrics to a DataFrame with one row per donor.
"""
function metrics_to_dataframe(collection::RepertoireCollection, metrics::Vector{Metrics})
    dfs = [metrics_to_dataframe(m, collection.donor_ids[i]) 
           for (i, m) in enumerate(metrics)]
    return reduce(vcat, dfs)
end
