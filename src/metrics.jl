# Core diversity and clonality metric computation functions
# These operate on frequency vectors and Repertoire objects

using Statistics: mean
import Random

# ============================================================================
# Core metric computation functions (operate on frequency vectors)
# ============================================================================

"""
    shannon_entropy(freqs::Vector{Float64}) -> Float64

Compute Shannon entropy H = -Σ(pᵢ log pᵢ) using natural logarithm.
Zero frequencies are handled correctly (0 * log(0) = 0).
"""
function shannon_entropy(freqs::Vector{Float64})
    isempty(freqs) && return 0.0
    H = 0.0
    @inbounds for p in freqs
        if p > 0
            H -= p * log(p)
        end
    end
    return H
end

"""
    simpson_index(freqs::Vector{Float64}) -> Float64

Compute Simpson's index D = Σpᵢ².
This is the probability that two randomly selected individuals belong to the same lineage.
"""
function simpson_index(freqs::Vector{Float64})
    isempty(freqs) && return 0.0
    D = 0.0
    @inbounds for p in freqs
        D += p * p
    end
    return D
end

"""
    berger_parker_index(freqs::Vector{Float64}) -> Float64

Compute Berger-Parker index: proportion of the most abundant lineage.
Assumes frequencies are sorted in descending order.
"""
function berger_parker_index(freqs::Vector{Float64})
    isempty(freqs) && return 0.0
    return first(freqs)  # Assumes sorted descending
end

"""
    gini_coefficient(freqs::Vector{Float64}) -> Float64

Compute Gini coefficient measuring inequality in lineage abundances.
Returns 0 for perfect equality, approaches 1 for maximum inequality.
"""
function gini_coefficient(freqs::Vector{Float64})
    n = length(freqs)
    n <= 1 && return 0.0
    
    sorted_freqs = sort(freqs)  # Sort ascending
    
    # Gini = (2Σᵢ i*xᵢ) / (n*Σxᵢ) - (n+1)/n
    cumsum_weighted = 0.0
    @inbounds for i in 1:n
        cumsum_weighted += i * sorted_freqs[i]
    end
    
    total = sum(sorted_freqs)
    total == 0 && return 0.0
    
    return (2 * cumsum_weighted) / (n * total) - (n + 1) / n
end

"""
    d50(counts::Vector{<:Real}) -> Int

Compute D50: minimum number of lineages comprising 50% of total repertoire.
Assumes counts are sorted in descending order.
"""
function d50(counts::Vector{T}) where T<:Real
    isempty(counts) && return 0
    
    total = sum(counts)
    total == zero(T) && return 0
    
    threshold = total / 2
    cumsum = zero(T)
    
    @inbounds for (i, c) in enumerate(counts)
        cumsum += c
        if cumsum >= threshold
            return i
        end
    end
    
    return length(counts)
end

"""
    chao1(counts::Vector{<:Integer}) -> Float64

Compute Chao1 richness estimator.
Estimates total species richness including unobserved species.

Formula: S_chao1 = S_obs + f₁²/(2f₂)
where f₁ = singletons, f₂ = doubletons

If f₂ = 0, uses bias-corrected formula: S_obs + f₁(f₁-1)/2
"""
function chao1(counts::Vector{T}) where T<:Integer
    isempty(counts) && return 0.0
    
    s_obs = length(counts)
    
    # Count singletons (f1) and doubletons (f2)
    f1 = count(==(one(T)), counts)
    f2 = count(==(T(2)), counts)
    
    if f2 > 0
        return s_obs + (f1 * f1) / (2 * f2)
    elseif f1 > 0
        # Bias-corrected formula when no doubletons
        return s_obs + (f1 * (f1 - 1)) / 2
    else
        return Float64(s_obs)
    end
end

# Float version - rounds to integers for singleton/doubleton counting
function chao1(counts::Vector{T}) where T<:AbstractFloat
    int_counts = round.(Int, counts)
    return chao1(int_counts)
end

# ============================================================================
# Hill numbers - unified diversity framework
# ============================================================================

"""
    hill_number(freqs::Vector{Float64}, q::Real) -> Float64

Compute Hill number of order q.
Hill numbers provide a unified framework for diversity metrics.

# Special cases
- q=0: Richness (number of species with non-zero frequency)
- q=1: exp(Shannon entropy) - uses L'Hôpital's rule limit
- q=2: Inverse Simpson (1/Σpᵢ²)
- q→∞: 1/max(pᵢ) (reciprocal of Berger-Parker)

# Formula
ᵍD = (Σpᵢᵍ)^(1/(1-q)) for q ≠ 1
"""
function hill_number(freqs::Vector{Float64}, q::Real)
    isempty(freqs) && return 0.0
    
    # Handle special cases
    if q ≈ 0
        return Float64(count(>(0), freqs))
    elseif q ≈ 1
        # L'Hôpital's rule: limit is exp(Shannon entropy)
        return exp(shannon_entropy(freqs))
    elseif isinf(q) && q > 0
        # q → ∞: reciprocal of max frequency
        max_p = maximum(freqs)
        return max_p > 0 ? 1.0 / max_p : 0.0
    else
        # General case
        sum_pq = 0.0
        @inbounds for p in freqs
            if p > 0
                sum_pq += p^q
            end
        end
        return sum_pq^(1 / (1 - q))
    end
end

"""
    hill_number(rep::Repertoire, q::Real) -> HillNumber{Q}

Compute Hill number of order q for a repertoire.
Returns a typed `HillNumber{Q}` for compile-time known orders.
"""
function hill_number(rep::Repertoire, q::Real)
    freqs = frequencies(rep)
    value = hill_number(freqs, q)
    return HillNumber{q}(value, richness(rep), total_count(rep))
end

# ============================================================================
# Repertoire-level metric functions
# ============================================================================

"""
    shannon_entropy(rep::Repertoire) -> Float64

Compute Shannon entropy for a repertoire.
"""
shannon_entropy(rep::Repertoire) = shannon_entropy(frequencies(rep))

"""
    simpson_index(rep::Repertoire) -> Float64

Compute Simpson's index for a repertoire.
"""
simpson_index(rep::Repertoire) = simpson_index(frequencies(rep))

"""
    simpson_diversity(rep::Repertoire) -> Float64

Compute Simpson's diversity (1 - D) for a repertoire.
"""
simpson_diversity(rep::Repertoire) = 1 - simpson_index(rep)

"""
    inverse_simpson(rep::Repertoire) -> Float64

Compute inverse Simpson index (1/D) for a repertoire.
"""
function inverse_simpson(rep::Repertoire)
    D = simpson_index(rep)
    return D > 0 ? 1 / D : Inf
end

"""
    clonality(rep::Repertoire) -> Float64

Compute clonality: 1 - normalized Shannon entropy.
High clonality indicates oligoclonal expansion.
Range: [0, 1] where 1 = maximally clonal (single dominant lineage).
"""
function clonality(rep::Repertoire)
    S = richness(rep)
    S <= 1 && return 1.0
    
    H = shannon_entropy(rep)
    H_max = log(S)
    
    return 1 - H / H_max
end

"""
    evenness(rep::Repertoire) -> Float64

Compute Pielou's evenness J = H / H_max.
Range: [0, 1] where 1 = perfectly even distribution.
"""
function evenness(rep::Repertoire)
    S = richness(rep)
    S <= 1 && return 0.0
    
    H = shannon_entropy(rep)
    H_max = log(S)
    
    return H / H_max
end

"""
    berger_parker_index(rep::Repertoire) -> Float64

Compute Berger-Parker index for a repertoire.
"""
berger_parker_index(rep::Repertoire) = berger_parker_index(frequencies(rep))

"""
    gini_coefficient(rep::Repertoire) -> Float64

Compute Gini coefficient for a repertoire.
"""
gini_coefficient(rep::Repertoire) = gini_coefficient(frequencies(rep))

"""
    d50(rep::Repertoire) -> Int

Compute D50 for a repertoire.
"""
d50(rep::Repertoire) = d50(counts(rep))

"""
    chao1(rep::Repertoire) -> Float64

Compute Chao1 richness estimator for a repertoire.
"""
chao1(rep::Repertoire) = chao1(counts(rep))

# ============================================================================
# Rarefaction
# ============================================================================

"""
    rarefaction(rep::Repertoire, depth::Integer; rng=nothing) -> Repertoire

Randomly subsample a repertoire to a specified depth (total count).

Many diversity metrics are sensitive to sample size—deeper sequencing captures more 
rare lineages. Rarefaction normalizes sample sizes by randomly subsampling all 
repertoires to the same depth, enabling fair comparisons.

# Arguments
- `rep`: Input repertoire
- `depth`: Target total count (must be ≤ current total count of `rep`)
- `rng`: Optional random number generator for reproducibility

# Returns
A new `Repertoire` with subsampled counts. Lineages that received zero counts 
in the subsample are removed.

# Example
```julia
# Compare two repertoires of different sizes
rep1 = read_repertoire("donor1.tsv", VJCdr3Definition())  # 50,000 sequences
rep2 = read_repertoire("donor2.tsv", VJCdr3Definition())  # 12,000 sequences

# Rarefy to the smaller depth
target = min(total_count(rep1), total_count(rep2))
rep1_rare = rarefaction(rep1, target)
rep2_rare = rarefaction(rep2, target)

# Now compare fairly
compute_metrics(rep1_rare)
compute_metrics(rep2_rare)

# For reproducibility, use a fixed RNG
using Random
rng = MersenneTwister(42)
rarefied = rarefaction(rep, 10000; rng=rng)
```

# Notes
- Rarefaction is stochastic; results vary between runs unless `rng` is fixed
- For robust estimates, consider averaging metrics over multiple rarefactions
- Metrics like Simpson index are naturally less sensitive to sample size
"""
function rarefaction(rep::Repertoire{T}, depth::Integer; rng=nothing) where T
    depth <= 0 && throw(ArgumentError("depth must be positive"))
    depth > total_count(rep) && throw(ArgumentError("depth exceeds total count"))
    
    if rng === nothing
        rng = Random.default_rng()
    end
    
    cnts = counts(rep)
    ids = lineage_ids(rep)
    
    # Create expanded pool of indices
    pool = Vector{Int}()
    sizehint!(pool, total_count(rep))
    for (i, c) in enumerate(cnts)
        for _ in 1:c
            push!(pool, i)
        end
    end
    
    # Sample without replacement
    Random.shuffle!(rng, pool)
    sampled = pool[1:depth]
    
    # Count occurrences
    new_counts = zeros(T, length(cnts))
    for i in sampled
        new_counts[i] += one(T)
    end
    
    # Filter out zeros
    nonzero_mask = new_counts .> 0
    
    return Repertoire(
        new_counts[nonzero_mask],
        ids[nonzero_mask];
        donor_id = donor_id(rep),
        metadata = copy(rep.metadata)
    )
end
