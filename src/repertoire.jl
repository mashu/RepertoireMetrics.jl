# Core Repertoire data structure with type-stable design

"""
    Repertoire{T<:Real}

Immutable representation of a B cell repertoire for diversity analysis.
Stores lineage counts in a type-stable manner for efficient metric computation.

# Type parameter
- `T`: Numeric type for counts (typically `Int` or `Float64`)

# Fields
- `counts::Vector{T}`: Sorted vector of lineage counts (descending order)
- `lineage_ids::Vector{String}`: Lineage identifiers corresponding to counts
- `donor_id::String`: Identifier for the donor/sample
- `total_count::T`: Total count (cached for efficiency)
- `metadata::Dict{String,Any}`: Optional metadata

# Construction
Use the constructors or `read_repertoire` for type-safe creation.

# Example
```julia
# From counts directly
rep = Repertoire([100, 50, 25, 10, 5], donor_id="Donor1")

# From DataFrame
rep = read_repertoire("data.tsv", VJCdr3Definition())
```
"""
struct Repertoire{T<:Real}
    counts::Vector{T}
    lineage_ids::Vector{String}
    donor_id::String
    total_count::T
    metadata::Dict{String,Any}
    
    function Repertoire{T}(
        counts::Vector{T},
        lineage_ids::Vector{String},
        donor_id::String,
        metadata::Dict{String,Any}
    ) where T<:Real
        length(counts) == length(lineage_ids) || 
            throw(ArgumentError("counts and lineage_ids must have same length"))
        all(c -> c >= zero(T), counts) || 
            throw(ArgumentError("counts must be non-negative"))
        
        # Sort by count descending for consistent representation
        perm = sortperm(counts, rev=true)
        sorted_counts = counts[perm]
        sorted_ids = lineage_ids[perm]
        total = sum(sorted_counts)
        
        new{T}(sorted_counts, sorted_ids, donor_id, total, metadata)
    end
end

# Convenience constructors
function Repertoire(
    counts::Vector{T},
    lineage_ids::Vector{String};
    donor_id::String = "",
    metadata::Dict{String,Any} = Dict{String,Any}()
) where T<:Real
    Repertoire{T}(counts, lineage_ids, donor_id, metadata)
end

function Repertoire(
    counts::Vector{T};
    donor_id::String = "",
    metadata::Dict{String,Any} = Dict{String,Any}()
) where T<:Real
    lineage_ids = ["lineage_$i" for i in 1:length(counts)]
    Repertoire{T}(counts, lineage_ids, donor_id, metadata)
end

# Type-stable accessors
"""
    richness(rep::Repertoire) -> Int

Return the number of unique lineages (richness/species count).
"""
richness(rep::Repertoire) = length(rep.counts)

"""
    total_count(rep::Repertoire) -> T

Return the total count across all lineages.
"""
total_count(rep::Repertoire) = rep.total_count

"""
    frequencies(rep::Repertoire{T}) -> Vector{Float64}

Return normalized frequencies (proportions) for each lineage.
Frequencies sum to 1.0.
"""
function frequencies(rep::Repertoire{T}) where T
    rep.total_count == zero(T) && return Float64[]
    return rep.counts ./ rep.total_count
end

"""
    counts(rep::Repertoire) -> Vector{T}

Return the count vector (sorted descending).
"""
counts(rep::Repertoire) = rep.counts

"""
    lineage_ids(rep::Repertoire) -> Vector{String}

Return lineage identifiers (sorted by count descending).
"""
lineage_ids(rep::Repertoire) = rep.lineage_ids

"""
    donor_id(rep::Repertoire) -> String

Return the donor/sample identifier.
"""
donor_id(rep::Repertoire) = rep.donor_id

# Base interface implementations
Base.length(rep::Repertoire) = richness(rep)
Base.isempty(rep::Repertoire) = isempty(rep.counts)

function Base.show(io::IO, rep::Repertoire{T}) where T
    println(io, "Repertoire{$T}:")
    println(io, "  Donor:      ", isempty(rep.donor_id) ? "<unknown>" : rep.donor_id)
    println(io, "  Lineages:   ", richness(rep))
    print(io, "  Total count: ", rep.total_count)
end

function Base.show(io::IO, ::MIME"text/plain", rep::Repertoire{T}) where T
    show(io, rep)
    if !isempty(rep)
        println(io)
        n_show = min(5, richness(rep))
        println(io, "  Top $n_show lineages:")
        for i in 1:n_show
            freq = rep.counts[i] / rep.total_count * 100
            println(io, "    $(rep.lineage_ids[i]): $(rep.counts[i]) ($(round(freq, digits=2))%)")
        end
        if richness(rep) > 5
            print(io, "    ... and $(richness(rep) - 5) more")
        end
    end
end

# ============================================================================
# Repertoire collection for multi-donor analysis
# ============================================================================

"""
    RepertoireCollection{T<:Real}

Collection of repertoires from multiple donors for comparative analysis.

# Fields
- `repertoires::Vector{Repertoire{T}}`: Vector of repertoires
- `donor_ids::Vector{String}`: Donor identifiers (for quick lookup)

# Example
```julia
collection = RepertoireCollection([rep1, rep2, rep3])
metrics = compute_metrics(collection)  # Returns vector of DiversityMetrics
```
"""
struct RepertoireCollection{T<:Real}
    repertoires::Vector{Repertoire{T}}
    donor_ids::Vector{String}
    
    function RepertoireCollection{T}(repertoires::Vector{Repertoire{T}}) where T<:Real
        donor_ids = [donor_id(r) for r in repertoires]
        new{T}(repertoires, donor_ids)
    end
end

function RepertoireCollection(repertoires::Vector{Repertoire{T}}) where T<:Real
    RepertoireCollection{T}(repertoires)
end

Base.length(col::RepertoireCollection) = length(col.repertoires)
Base.getindex(col::RepertoireCollection, i::Int) = col.repertoires[i]
Base.iterate(col::RepertoireCollection) = iterate(col.repertoires)
Base.iterate(col::RepertoireCollection, state) = iterate(col.repertoires, state)

function Base.show(io::IO, col::RepertoireCollection{T}) where T
    print(io, "RepertoireCollection{$T} with $(length(col)) repertoires")
end
