# Abstract types and concrete implementations for lineage definition strategies
# Uses Julia's type system for compile-time dispatch and type stability

"""
    AbstractLineageDefinition

Abstract supertype for all lineage definition strategies.
Concrete subtypes define how sequences are grouped into clonal lineages.

# Implementing a new strategy
To implement a new lineage definition strategy:
1. Create a subtype of `AbstractLineageDefinition`
2. Implement `lineage_key(strategy, row)` returning a hashable key
"""
abstract type AbstractLineageDefinition end

"""
    LineageIDDefinition <: AbstractLineageDefinition

Define lineages using the `lineage_id` column directly.
This is the simplest strategy when lineage assignment has already been performed.

# Fields
- `column::Symbol`: Column name containing lineage IDs (default: `:lineage_id`)
"""
struct LineageIDDefinition <: AbstractLineageDefinition
    column::Symbol
end

LineageIDDefinition() = LineageIDDefinition(:lineage_id)

"""
    VJCdr3Definition <: AbstractLineageDefinition

Define lineages using the combination of V gene, J gene, and CDR3 sequence.
This is consistent with LineageCollapse.jl's default behavior.

# Fields
- `v_column::Symbol`: Column name for V gene call (default: `:v_call`)
- `j_column::Symbol`: Column name for J gene call (default: `:j_call`)
- `cdr3_column::Symbol`: Column name for CDR3 sequence (default: `:cdr3`)
- `use_first_allele::Bool`: Use only first allele from calls (default: `true`)
"""
struct VJCdr3Definition <: AbstractLineageDefinition
    v_column::Symbol
    j_column::Symbol
    cdr3_column::Symbol
    use_first_allele::Bool
end

function VJCdr3Definition(;
    v_column::Symbol = :v_call,
    j_column::Symbol = :j_call,
    cdr3_column::Symbol = :cdr3,
    use_first_allele::Bool = true
)
    VJCdr3Definition(v_column, j_column, cdr3_column, use_first_allele)
end

"""
    CustomDefinition{F} <: AbstractLineageDefinition

Define lineages using a custom function that extracts a key from each row.

# Fields
- `key_func::F`: Function `row -> key` that returns a hashable lineage key

# Example
```julia
# Group by V gene family only
strategy = CustomDefinition(row -> first(split(string(row.v_call), "-")))
```
"""
struct CustomDefinition{F} <: AbstractLineageDefinition
    key_func::F
end

# ============================================================================
# Lineage key extraction - multiple dispatch on strategy type
# ============================================================================

"""
    lineage_key(strategy::AbstractLineageDefinition, row) -> key

Extract a lineage key from a data row using the given strategy.
Returns a hashable value that identifies the lineage.
"""
function lineage_key end

function lineage_key(strategy::LineageIDDefinition, row)
    getproperty(row, strategy.column)
end

function lineage_key(strategy::VJCdr3Definition, row)
    v = getproperty(row, strategy.v_column)
    j = getproperty(row, strategy.j_column)
    cdr3 = getproperty(row, strategy.cdr3_column)
    
    if strategy.use_first_allele
        v = first_allele(v)
        j = first_allele(j)
    end
    
    return (v, j, cdr3)
end

function lineage_key(strategy::CustomDefinition, row)
    strategy.key_func(row)
end

"""
    first_allele(call::AbstractString) -> String

Extract the first allele from a gene call string.
Handles comma-separated multiple calls and allele notation.

# Examples
```julia
first_allele("IGHV1-2*01,IGHV1-2*02") # "IGHV1-2*01"
first_allele("IGHV1-2*01")             # "IGHV1-2*01"
```
"""
function first_allele(call::AbstractString)
    isempty(call) && return ""
    first_call = first(split(call, ','))
    return strip(String(first_call))
end

first_allele(call::Missing) = missing
first_allele(::Nothing) = nothing

# ============================================================================
# Abstract type for metric results
# ============================================================================

"""
    AbstractMetricResult

Abstract supertype for metric computation results.
"""
abstract type AbstractMetricResult end

# ============================================================================
# Hill numbers - unified diversity framework
# ============================================================================

"""
    HillNumber{Q} <: AbstractMetricResult

Hill number of order Q, providing a unified framework for diversity.

# Type parameter
- `Q`: Order of the Hill number (compile-time constant for type stability)

# Fields  
- `value::Float64`: The computed Hill number
- `richness::Int`: Number of unique lineages
- `total_count::Int`: Total count

# Special cases
- Q=0: Richness (number of species)
- Q=1: Exponential of Shannon entropy
- Q=2: Inverse Simpson index
- Q=∞: Reciprocal of Berger-Parker index
"""
struct HillNumber{Q} <: AbstractMetricResult
    value::Float64
    richness::Int
    total_count::Int
end

function Base.show(io::IO, h::HillNumber{Q}) where Q
    print(io, "HillNumber{$Q}(", round(h.value, digits=4), ")")
end
