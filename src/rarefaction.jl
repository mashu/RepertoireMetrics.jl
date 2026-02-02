# Rarefaction curves and averaged rarefaction for depth-normalized comparisons
#
# Provides:
#   - RarefactionCurve: result type for rarefaction curves
#   - rarefaction_curve: compute a metric at multiple subsampled depths
#   - average_rarefaction: average metrics over repeated draws at one depth

import Random
using Statistics: mean, std

# ============================================================================
# Result types
# ============================================================================

"""
    RarefactionPoint

A single point on a rarefaction curve: metric value at a given depth.
"""
struct RarefactionPoint
    depth::Int
    mean::Float64
    std::Float64
    n_iter::Int
end

"""
    RarefactionCurve

Result of a rarefaction curve computation.

# Fields
- `points::Vector{RarefactionPoint}`: Metric values at each depth
- `metric::Symbol`: Name of the metric computed
- `donor_id::String`: Donor identifier
- `original_depth::Int`: Original (unrarefied) total count

# Example
```julia
curve = rarefaction_curve(rep, Richness())
for p in curve.points
    println("\$(p.depth) => \$(p.mean) ± \$(p.std)")
end
```
"""
struct RarefactionCurve
    points::Vector{RarefactionPoint}
    metric::Symbol
    donor_id::String
    original_depth::Int
end

function Base.show(io::IO, rc::RarefactionCurve)
    n = length(rc.points)
    dmin = n > 0 ? first(rc.points).depth : 0
    dmax = n > 0 ? last(rc.points).depth : 0
    print(io, "RarefactionCurve(:$(rc.metric), $(n) depths from $(dmin) to $(dmax)")
    !isempty(rc.donor_id) && print(io, ", donor=$(rc.donor_id)")
    print(io, ")")
end

"""
    AveragedMetrics

Result of averaging metrics over multiple rarefaction draws at a fixed depth.

# Fields
- `means::Dict{Symbol, Float64}`: Mean of each metric
- `stds::Dict{Symbol, Float64}`: Standard deviation of each metric
- `depth::Int`: Rarefaction depth used
- `n_iter::Int`: Number of iterations
- `metric_names::Vector{Symbol}`: Ordered metric names
"""
struct AveragedMetrics
    means::Dict{Symbol, Float64}
    stds::Dict{Symbol, Float64}
    depth::Int
    n_iter::Int
    metric_names::Vector{Symbol}
end

function Base.getproperty(m::AveragedMetrics, name::Symbol)
    name in (:means, :stds, :depth, :n_iter, :metric_names) && return getfield(m, name)
    return get(getfield(m, :means), name, missing)
end

function Base.show(io::IO, am::AveragedMetrics)
    println(io, "AveragedMetrics (depth=$(am.depth), n_iter=$(am.n_iter)):")
    for name in am.metric_names
        m = am.means[name]
        s = am.stds[name]
        println(io, "  $(name): $(round(m, digits=4)) ± $(round(s, digits=4))")
    end
end

# ============================================================================
# Optimized sampling kernel
# ============================================================================

"""
    _build_pool(rep::Repertoire) -> Vector{Int}

Expand counts into an index pool for shuffle-based sampling.
Cached outside the iteration loop for efficiency.
"""
function _build_pool(rep::Repertoire{T}) where T
    cnts = counts(rep)
    tc = total_count(rep)
    pool = Vector{Int}(undef, tc)
    pos = 1
    @inbounds for (i, c) in enumerate(cnts)
        for _ in one(T):c
            pool[pos] = i
            pos += 1
        end
    end
    pool
end

"""
    _sample_counts!(new_counts, pool, depth, rng) -> nothing

In-place shuffle-and-count. Mutates `new_counts` (must be pre-zeroed)
and partially shuffles `pool`.
"""
function _sample_counts!(new_counts::Vector{T}, pool::Vector{Int},
                         depth::Integer, rng) where T
    n = length(pool)
    # Fisher-Yates partial shuffle: only shuffle first `depth` elements
    @inbounds for i in 1:depth
        j = rand(rng, i:n)
        pool[i], pool[j] = pool[j], pool[i]
    end
    @inbounds for k in 1:depth
        new_counts[pool[k]] += one(T)
    end
    nothing
end

"""
    _rarefied_repertoire(rep, pool, new_counts, depth, rng) -> Repertoire

Build a rarefied Repertoire from a pre-allocated pool and count buffer.
Zeros `new_counts` before sampling.
"""
function _rarefied_repertoire(rep::Repertoire{T}, pool::Vector{Int},
                              new_counts::Vector{T}, depth::Integer,
                              rng) where T
    fill!(new_counts, zero(T))
    _sample_counts!(new_counts, pool, depth, rng)

    nonzero = findall(>(zero(T)), new_counts)
    ids = lineage_ids(rep)

    Repertoire(
        new_counts[nonzero],
        ids[nonzero];
        donor_id = donor_id(rep),
        metadata = copy(rep.metadata)
    )
end

# ============================================================================
# Default depth schedule
# ============================================================================

"""
    default_depths(rep::Repertoire; n_steps::Integer=20) -> Vector{Int}

Generate logarithmically-spaced depth steps from 1% of total count
to the full total count.
"""
function default_depths(rep::Repertoire; n_steps::Integer=20)
    tc = total_count(rep)
    tc <= 0 && return Int[]
    lo = max(1, tc ÷ 100)
    steps = unique(round.(Int, exp.(range(log(lo), log(tc), length=n_steps))))
    sort!(steps)
    steps
end

# ============================================================================
# Rarefaction curve
# ============================================================================

"""
    rarefaction_curve(rep::Repertoire, metric::AbstractMetric;
                      depths=default_depths(rep),
                      n_iter::Integer=50,
                      rng=nothing) -> RarefactionCurve

Compute a rarefaction curve: the given metric evaluated at multiple
subsampled depths, averaged over `n_iter` random draws per depth.

# Arguments
- `rep`: Input repertoire
- `metric`: Which metric to compute (e.g. `Richness()`, `Clonality()`, `InverseSimpson()`)
- `depths`: Vector of depths to evaluate (default: ~20 log-spaced steps)
- `n_iter`: Number of random draws per depth for averaging
- `rng`: Random number generator for reproducibility

# Example
```julia
using Random

# Richness rarefaction curve
curve = rarefaction_curve(rep, Richness(); n_iter=50, rng=MersenneTwister(42))

# Clonality curve
curve_c = rarefaction_curve(rep, Clonality(); depths=100:100:10000, n_iter=30)

# Compare donors
for rep in collection
    curve = rarefaction_curve(rep, Richness())
    println(donor_id(rep), ": ", last(curve.points).mean)
end
```
"""
function rarefaction_curve(rep::Repertoire{T}, metric::AbstractMetric;
                           depths::AbstractVector{<:Integer}=default_depths(rep),
                           n_iter::Integer=50,
                           rng=nothing) where T
    n_iter >= 1 || throw(ArgumentError("n_iter must be ≥ 1"))
    tc = total_count(rep)

    if rng === nothing
        rng = Random.default_rng()
    end

    pool = _build_pool(rep)
    new_counts = zeros(T, richness(rep))
    values = Vector{Float64}(undef, n_iter)
    points = Vector{RarefactionPoint}(undef, length(depths))

    sym = metric_symbol(metric)

    for (di, depth) in enumerate(depths)
        depth <= 0 && throw(ArgumentError("depth must be positive"))
        depth > tc && throw(ArgumentError("depth $depth exceeds total count $tc"))

        for iter in 1:n_iter
            rarefied = _rarefied_repertoire(rep, pool, new_counts, depth, rng)
            values[iter] = Float64(compute_metric(rarefied, metric))
        end

        points[di] = RarefactionPoint(depth, mean(values), std(values), n_iter)
    end

    RarefactionCurve(points, sym, donor_id(rep), tc)
end

# ============================================================================
# Averaged rarefaction at a single depth
# ============================================================================

"""
    average_rarefaction(rep::Repertoire, depth::Integer,
                        metrics::MetricSet=ALL_METRICS;
                        n_iter::Integer=50,
                        rng=nothing) -> AveragedMetrics

Rarefy a repertoire to `depth` repeatedly and return averaged metrics.

# Arguments
- `rep`: Input repertoire
- `depth`: Target subsampling depth
- `metrics`: Which metrics to compute (default: all)
- `n_iter`: Number of random draws for averaging
- `rng`: Random number generator for reproducibility

# Example
```julia
# Normalize all donors to the shallowest depth
target = minimum(total_count(rep) for rep in collection)
averaged = [average_rarefaction(rep, target, ROBUST_METRICS; n_iter=100) for rep in collection]

# Access results
averaged[1].clonality          # mean clonality
averaged[1].stds[:clonality]   # standard deviation
```
"""
function average_rarefaction(rep::Repertoire{T}, depth::Integer,
                             ms::MetricSet=ALL_METRICS;
                             n_iter::Integer=50,
                             rng=nothing) where T
    n_iter >= 1 || throw(ArgumentError("n_iter must be ≥ 1"))
    depth <= 0 && throw(ArgumentError("depth must be positive"))
    depth > total_count(rep) && throw(ArgumentError("depth exceeds total count"))

    if rng === nothing
        rng = Random.default_rng()
    end

    pool = _build_pool(rep)
    new_counts = zeros(T, richness(rep))

    # Collect metric names from the MetricSet
    names = Symbol[metric_symbol(m) for m in ms.metrics]
    accum = Dict(name => Vector{Float64}(undef, n_iter) for name in names)

    for iter in 1:n_iter
        rarefied = _rarefied_repertoire(rep, pool, new_counts, depth, rng)
        result = compute_metrics(rarefied, ms)
        for name in names
            accum[name][iter] = result.values[name]
        end
    end

    means = Dict(name => mean(accum[name]) for name in names)
    stds = Dict(name => std(accum[name]) for name in names)

    AveragedMetrics(means, stds, depth, n_iter, names)
end

# Single metric convenience
function average_rarefaction(rep::Repertoire, depth::Integer, m::AbstractMetric;
                             n_iter::Integer=50, rng=nothing)
    average_rarefaction(rep, depth, MetricSet(m); n_iter, rng)
end

# ============================================================================
# DataFrame conversion for curves
# ============================================================================

"""
    rarefaction_curve_to_dataframe(curve::RarefactionCurve) -> DataFrame

Convert a rarefaction curve to a DataFrame for plotting.

# Columns
- `depth::Int`
- `mean::Float64`
- `std::Float64`
- `donor_id::String`
- `metric::Symbol`
"""
function rarefaction_curve_to_dataframe(curve::RarefactionCurve)
    DataFrame(
        depth = [p.depth for p in curve.points],
        mean = [p.mean for p in curve.points],
        std = [p.std for p in curve.points],
        donor_id = fill(curve.donor_id, length(curve.points)),
        metric = fill(curve.metric, length(curve.points))
    )
end

"""
    rarefaction_curve_to_dataframe(curves::Vector{RarefactionCurve}) -> DataFrame

Combine multiple rarefaction curves (e.g. from different donors) into
a single DataFrame for faceted plotting.
"""
function rarefaction_curve_to_dataframe(curves::Vector{RarefactionCurve})
    reduce(vcat, rarefaction_curve_to_dataframe.(curves))
end

"""
    averaged_metrics_to_dataframe(collection::RepertoireCollection,
                                  averaged::Vector{AveragedMetrics}) -> DataFrame

Convert averaged rarefaction results for a collection to a DataFrame.
Includes mean and std columns for each metric.
"""
function averaged_metrics_to_dataframe(collection::RepertoireCollection,
                                       averaged::Vector{AveragedMetrics})
    length(collection) == length(averaged) ||
        throw(ArgumentError("collection and averaged must have same length"))

    rows = []
    for (i, am) in enumerate(averaged)
        row = Dict{Symbol, Any}(:donor_id => collection.donor_ids[i], :depth => am.depth)
        for name in am.metric_names
            row[name] = am.means[name]
            row[Symbol(name, :_std)] = am.stds[name]
        end
        push!(rows, row)
    end

    DataFrame(rows)
end
