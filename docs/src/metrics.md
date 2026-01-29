# Diversity and Clonality Metrics

This page provides detailed explanations of all metrics implemented in RepertoireMetrics.jl, including mathematical formulas, interpretation, and usage guidelines.

## Quick Reference: Which Metric Should I Use?

| Your Question | Recommended Metrics | Why |
|--------------|---------------------|-----|
| How many unique clones? | `Richness`, `Chao1` | Chao1 estimates unseen clones |
| Is the repertoire clonally expanded? | `Clonality`, `Gini`, `D50` | Directly measure dominance |
| How diverse overall? | `Shannon Diversity`, `Inverse Simpson` | Effective number of clones |
| How dominant is the top clone? | `Berger-Parker`, `D50` | Simple, interpretable |
| Comparing samples of different sizes? | `Simpson`, `Berger-Parker`, `Gini` | Frequency-based, depth-robust |
| Need richness comparison across depths? | `rarefaction` + `Richness` | Rarefaction equalizes depth |
| Full diversity profile? | `Hill numbers` with varying q | Unified theoretical framework |

### Common Pitfalls

| Pitfall | Problem | Solution |
|---------|---------|----------|
| Comparing richness across sample sizes | Larger samples always show more clones | Use rarefaction or Chao1 |
| Ignoring rare clones | Simpson/Berger-Parker may miss important rare populations | Also compute Shannon or richness |
| Over-interpreting small differences | Sampling noise can be substantial | Use confidence intervals or bootstrapping |
| Assuming one metric tells the whole story | Different metrics capture different aspects | Report multiple complementary metrics |

### Metric Relationships

```
Clonality = 1 - Evenness = 1 - (Shannon Entropy / log(Richness))
Shannon Diversity = exp(Shannon Entropy) = Hill number q=1
Inverse Simpson = 1 / Simpson Index = Hill number q=2
```

---

## Notation

Throughout this document:

- ``S`` = number of unique lineages (richness)
- ``N`` = total count (sequences or cells)
- ``n_i`` = count of lineage ``i``
- ``p_i = n_i / N`` = relative frequency of lineage ``i``
- ``f_k`` = number of lineages with exactly ``k`` sequences (frequency of frequencies)

## Basic Counts

### Richness (S)

**Formula:**

```math
S = \sum_{i=1}^{S} \mathbf{1}[n_i > 0]
```

**Interpretation:** The total number of unique lineages observed in the repertoire. This is the simplest measure of diversity—more lineages means more diverse.

**Limitations:** Richness is highly sensitive to sample size. Larger samples almost always have higher richness, making direct comparisons difficult without rarefaction.

**Use when:** You want a basic count of unique clones, especially when comparing samples of similar size.

### Total Count (N)

**Formula:**

```math
N = \sum_{i=1}^{S} n_i
```

**Interpretation:** The total number of sequences or cells in the repertoire. Essential for normalization and understanding sampling depth.

---

## Shannon Entropy Family

### Shannon Entropy (H)

**Formula:**

```math
H = -\sum_{i=1}^{S} p_i \log(p_i)
```

where we use the convention that ``0 \log(0) = 0``.

**Interpretation:** Shannon entropy quantifies the uncertainty in predicting the lineage of a randomly selected sequence. Higher entropy means more uncertainty, indicating a more diverse repertoire.

**Range:** ``[0, \log(S)]``
- ``H = 0`` when one lineage dominates completely
- ``H = \log(S)`` when all lineages are equally abundant

**Properties:**
- Uses natural logarithm (nats), though log₂ (bits) is also common
- Sensitive to rare lineages but less so than richness
- Affected by both richness and evenness

### Shannon Diversity (Exponential of H)

**Formula:**

```math
{}^1D = e^H = \exp\left(-\sum_{i=1}^{S} p_i \log(p_i)\right)
```

**Interpretation:** The "effective number of lineages"—how many equally abundant lineages would produce the same entropy. Also known as the Hill number of order 1.

**Range:** ``[1, S]``
- ``{}^1D = 1`` for a single dominant lineage
- ``{}^1D = S`` for perfect evenness

**Why use this over raw entropy?** Shannon diversity is in interpretable units (number of lineages) rather than abstract entropy units.

### Normalized Shannon Entropy

**Formula:**

```math
H_{norm} = \frac{H}{\log(S)} = \frac{H}{H_{max}}
```

**Interpretation:** Shannon entropy scaled to [0, 1]. Represents how close the distribution is to maximum entropy (perfect evenness).

**Range:** ``[0, 1]``
- ``H_{norm} = 0`` for a single dominant lineage
- ``H_{norm} = 1`` for perfect evenness

### Clonality

**Formula:**

```math
\text{Clonality} = 1 - H_{norm} = 1 - \frac{H}{\log(S)}
```

**Interpretation:** The complement of normalized Shannon entropy. High clonality indicates oligoclonal expansion—a repertoire dominated by few lineages.

**Range:** ``[0, 1]``
- ``\text{Clonality} = 0`` for perfect evenness (maximally polyclonal)
- ``\text{Clonality} = 1`` for a single lineage (monoclonal)

**Use when:** Investigating antigen-driven clonal expansion, disease states, or immune responses.

### Evenness (Pielou's J)

**Formula:**

```math
J = \frac{H}{\log(S)} = H_{norm}
```

**Interpretation:** Identical to normalized Shannon entropy. Measures how evenly sequences are distributed among lineages, independent of richness.

**Range:** ``[0, 1]``
- ``J = 0`` for maximum unevenness
- ``J = 1`` for perfect evenness

**Note:** Evenness = 1 - Clonality.

---

## Simpson Family

### Simpson's Index (D)

**Formula:**

```math
D = \sum_{i=1}^{S} p_i^2
```

**Interpretation:** The probability that two randomly selected sequences belong to the same lineage. Higher D means less diversity (more concentrated).

**Range:** ``[1/S, 1]``
- ``D = 1/S`` for perfect evenness
- ``D = 1`` for a single lineage

**Note:** This is sometimes called Simpson's concentration index.

### Simpson's Diversity (Gini-Simpson Index)

**Formula:**

```math
1 - D = 1 - \sum_{i=1}^{S} p_i^2
```

**Interpretation:** The probability that two randomly selected sequences belong to *different* lineages. Higher values indicate more diversity.

**Range:** ``[0, 1 - 1/S]``
- ``1 - D \to 0`` as repertoire becomes monoclonal
- ``1 - D \to 1`` as repertoire becomes maximally diverse

### Inverse Simpson (Hill Number q=2)

**Formula:**

```math
{}^2D = \frac{1}{D} = \frac{1}{\sum_{i=1}^{S} p_i^2}
```

**Interpretation:** The "effective number of lineages" when weighted by frequency squared. Emphasizes dominant lineages more than Shannon diversity.

**Range:** ``[1, S]``
- ``{}^2D = 1`` for a single lineage
- ``{}^2D = S`` for perfect evenness

**Use when:** You want a diversity measure less sensitive to rare lineages than Shannon entropy.

---

## Dominance Metrics

### Berger-Parker Index

**Formula:**

```math
d = p_{max} = \max_i(p_i)
```

**Interpretation:** The relative frequency of the most abundant lineage. A simple measure of dominance.

**Range:** ``[1/S, 1]``
- ``d = 1/S`` for perfect evenness
- ``d = 1`` when one lineage has all sequences

**Use when:** You want a simple, interpretable measure of how dominant the top clone is.

### D50

**Formula:**

```math
D50 = \min\left\{k : \sum_{i=1}^{k} p_{(i)} \geq 0.5\right\}
```

where ``p_{(1)} \geq p_{(2)} \geq \cdots`` are frequencies sorted in descending order.

**Interpretation:** The minimum number of lineages needed to account for 50% of the repertoire. Low D50 indicates high clonality (few lineages dominate).

**Range:** ``[1, \lceil S/2 \rceil]``
- ``D50 = 1`` when the top lineage has ≥50% of sequences
- ``D50 = \lceil S/2 \rceil`` for perfect evenness

**Use when:** You want an intuitive measure of how "top-heavy" the distribution is.

---

## Inequality Metrics

### Gini Coefficient

**Formula:**

```math
G = \frac{\sum_{i=1}^{S} \sum_{j=1}^{S} |p_i - p_j|}{2S \sum_{i=1}^{S} p_i} = \frac{2\sum_{i=1}^{S} i \cdot p_{(i)}}{S} - \frac{S+1}{S}
```

where ``p_{(1)} \leq p_{(2)} \leq \cdots`` are frequencies sorted in ascending order.

**Interpretation:** Measures inequality in lineage abundances. Originally developed for economic inequality (wealth distribution).

**Range:** ``[0, 1]``
- ``G = 0`` for perfect equality (all lineages equal)
- ``G \to 1`` for maximum inequality

**Geometric interpretation:** Half the relative mean absolute difference between all pairs of lineages.

**Use when:** You want to quantify how unequally sequences are distributed, independent of richness.

---

## Richness Estimators

### Chao1

**Formula:**

```math
\hat{S}_{Chao1} = S_{obs} + \frac{f_1^2}{2 f_2}
```

where:
- ``S_{obs}`` = observed richness
- ``f_1`` = number of singletons (lineages with count = 1)
- ``f_2`` = number of doubletons (lineages with count = 2)

If ``f_2 = 0``, use the bias-corrected formula:

```math
\hat{S}_{Chao1} = S_{obs} + \frac{f_1(f_1 - 1)}{2}
```

**Interpretation:** Estimates the true number of lineages, including those not observed due to sampling limitations. Based on the idea that rare species inform us about unobserved species.

**Range:** ``[S_{obs}, \infty)``

**Use when:** You suspect many rare lineages exist that weren't sampled, and want to estimate true richness.

---

## Hill Numbers (Unified Framework)

### Definition

**Formula:**

```math
{}^qD = \left(\sum_{i=1}^{S} p_i^q\right)^{1/(1-q)}
```

**Special cases:**

| Order | Formula | Interpretation |
|-------|---------|----------------|
| ``q = 0`` | ``{}^0D = S`` | Richness |
| ``q \to 1`` | ``{}^1D = e^H`` | exp(Shannon entropy) |
| ``q = 2`` | ``{}^2D = 1/D`` | Inverse Simpson |
| ``q \to \infty`` | ``{}^\infty D = 1/p_{max}`` | Reciprocal of Berger-Parker |

**Interpretation:** Hill numbers provide a unified framework for diversity. The order ``q`` controls sensitivity to rare vs. common lineages:
- Low ``q``: More sensitive to rare lineages
- High ``q``: More sensitive to dominant lineages

**Use when:** You want to explore the full diversity profile of a repertoire, or need a theoretically grounded framework that unifies different diversity measures.

---

## Comparing Metrics

### Relationship Between Diversity and Clonality

| Aspect | Diversity Metrics | Clonality Metrics |
|--------|------------------|-------------------|
| **Focus** | Variability | Dominance |
| **High value** | Polyclonal | Oligoclonal |
| **Examples** | Shannon diversity, Inverse Simpson | Clonality, Gini, D50 (inverse) |

Key relationships:
- `Clonality = 1 - Evenness`
- `Shannon Diversity = exp(Shannon Entropy)`
- `Inverse Simpson = 1 / Simpson Index`

### Sensitivity to Different Features

| Metric | Sensitive to rare lineages | Sensitive to dominant lineages | Sensitive to sample size |
|--------|---------------------------|-------------------------------|-------------------------|
| Richness | Very high | Low | Very high |
| Shannon Entropy | Moderate | Moderate | Moderate |
| Simpson Index | Low | High | Low |
| Berger-Parker | Very low | Very high | Low |
| Chao1 | Very high | Low | High |

### Recommendations

1. **For detecting clonal expansion:** Use Clonality, Gini, D50, or Berger-Parker.

2. **For a complete picture:** Compute multiple metrics or use the Hill number profile across different orders.

3. **For estimating true diversity:** Use Chao1, especially when many singletons are present.

---

## Handling Different Sequencing Depths

When comparing repertoires with different total counts, you have several strategies:

### Strategy 1: Use Depth-Robust Metrics (Recommended)

Many metrics are computed from **frequencies** (proportions), not raw counts. Mathematically, these metrics depend only on the shape of the distribution, not the total N:

| Metric | Depth Sensitivity | Why |
|--------|-------------------|-----|
| Simpson Index | **Low** | Dominated by abundant clones |
| Inverse Simpson | **Low** | Same as Simpson |
| Berger-Parker | **Very Low** | Only looks at top clone |
| Gini Coefficient | **Low** | Measures inequality of frequencies |
| Shannon Entropy | **Moderate** | More sensitive to rare clones |
| Richness | **Very High** | More sequences = more rare clones observed |
| Chao1 | **High** | Based on singletons/doubletons |

**When to use:** When you have accurate count data and want to avoid information loss.

```julia
# Use the predefined ROBUST_METRICS set (includes Depth for reporting)
metrics = compute_metrics(rep, ROBUST_METRICS)

# Or select specific robust metrics
metrics = compute_metrics(rep, Depth() + SimpsonDiversity() + Clonality() + GiniCoefficient())
```

**Always report depth:** Include `Depth()` in your metrics so readers know the sequencing depth of each sample.

### Strategy 2: Rarefaction (Conservative)

Randomly subsample all repertoires to the same depth. See [Rarefaction](quickstart.md#Rarefaction-for-Comparing-Uneven-Samples).

**Pros:** Makes richness comparable; theoretically rigorous  
**Cons:** Discards information from deeper samples; introduces stochasticity

**When to use:** When comparing richness specifically, or when you want the most conservative comparison.

### Strategy 3: Normalize to Total Count

If you have the `count` column, the total count reflects the original sequencing depth:

```julia
rep = read_repertoire("data.tsv", VJCdr3Definition())
println("Total sequences: ", total_count(rep))  # Original depth preserved
println("Unique lineages: ", richness(rep))     # Observed clones

# Frequencies are already normalized
freqs = frequencies(rep)  # Sum to 1.0
```

The frequency-based metrics already use this normalization internally.

### Strategy 4: Report Depth Alongside Metrics

Always report sequencing depth so readers can assess comparability:

```julia
metrics = compute_metrics(rep)
df = metrics_to_dataframe(rep, metrics)
# DataFrame includes total_count for context
```

### The Sampling Reality

Even with "depth-robust" metrics, there's a subtlety: **observed frequencies are estimates of true population frequencies**. With shallow sequencing:
- Rare clones may not be observed at all
- Observed frequencies of rare clones have higher variance

This means Simpson/Shannon computed from shallow data may still differ from deep data, not because the metric is biased, but because the *input frequencies* are different estimates of the same underlying distribution.

**Bottom line:** For robust comparisons, prefer Simpson-family metrics and report total counts. Use rarefaction when richness comparison is essential.
