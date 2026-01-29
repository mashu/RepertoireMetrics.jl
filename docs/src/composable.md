# Composable Metrics

RepertoireMetrics.jl provides a flexible system for selecting which metrics to compute. This is useful when you:

- Only need a few specific metrics
- Want to optimize performance by avoiding unnecessary computations
- Are building analysis pipelines with specific requirements

## Basic Usage

### The `+` Operator

Combine metrics using the `+` operator:

```julia
using RepertoireMetrics

rep = Repertoire([100, 50, 25, 10, 5])

# Select specific metrics
metrics = ShannonEntropy() + Clonality() + D50()
result = compute_metrics(rep, metrics)

println(result)
# Output:
# Metrics (3 computed):
#   shannon_entropy: 1.2891
#   clonality: 0.1991
#   d50: 2
```

### Default: All Metrics

Without a second argument, `compute_metrics` computes all available metrics:

```julia
result = compute_metrics(rep)  # Same as compute_metrics(rep, ALL_METRICS)
```

### Available Metric Types

Each metric has a corresponding type:

| Type | Computed Metric |
|------|-----------------|
| `Richness()` | Number of unique lineages |
| `TotalCount()` | Total sequence count |
| `ShannonEntropy()` | Shannon entropy H |
| `ShannonDiversity()` | exp(H) |
| `NormalizedShannon()` | H / log(S) |
| `SimpsonIndex()` | Σpᵢ² |
| `SimpsonDiversity()` | 1 - Σpᵢ² |
| `InverseSimpson()` | 1 / Σpᵢ² |
| `BergerParker()` | max(pᵢ) |
| `Evenness()` | Pielou's J |
| `Clonality()` | 1 - normalized Shannon |
| `GiniCoefficient()` | Gini index |
| `D50()` | D50 |
| `Chao1()` | Chao1 estimator |
| `MeanLength()` | Mean sequence length* |
| `MedianLength()` | Median sequence length* |
| `StdLength()` | Std dev of sequence length* |
| `MinLength()` | Minimum sequence length* |
| `MaxLength()` | Maximum sequence length* |

\* Length metrics require `length_column` to be specified when creating the repertoire.

### Accessing Results

The result is a `Metrics` object with property access:

```julia
result = compute_metrics(rep, ShannonEntropy() + Clonality())

# Access by property name
println(result.shannon_entropy)
println(result.clonality)

# Missing metrics return `missing`
println(result.d50)  # missing (wasn't computed)
```

## Predefined MetricSets

For convenience, several predefined metric sets are available:

### `ALL_METRICS` (default)

All available metrics:

```julia
result = compute_metrics(rep)  # Uses ALL_METRICS by default
```

### `DIVERSITY_METRICS`

Metrics focused on repertoire diversity:

```julia
result = compute_metrics(rep, DIVERSITY_METRICS)
# Includes: richness, total_count, shannon_entropy, shannon_diversity,
#           simpson_diversity, inverse_simpson, evenness
```

### `CLONALITY_METRICS`

Metrics focused on clonal expansion:

```julia
result = compute_metrics(rep, CLONALITY_METRICS)
# Includes: richness, total_count, clonality, gini_coefficient,
#           berger_parker, d50
```

### `RICHNESS_METRICS`

Richness and its estimator:

```julia
result = compute_metrics(rep, RICHNESS_METRICS)
# Includes: richness, total_count, chao1
```

### `LENGTH_METRICS`

All sequence length statistics (requires `length_column` during repertoire creation):

```julia
rep = read_repertoire("data.tsv", VJCdr3Definition(); length_column=:cdr3)
result = compute_metrics(rep, LENGTH_METRICS)
# Includes: mean_length, median_length, std_length, min_length, max_length
```

## Combining Sets and Metrics

You can combine predefined sets with additional metrics:

```julia
# Add Chao1 to diversity metrics
extended = DIVERSITY_METRICS + Chao1()
result = compute_metrics(rep, extended)
```

## Multi-Donor Analysis

The composable system works with repertoire collections:

```julia
collection = read_repertoires_from_directory("data/", VJCdr3Definition())

# Compute only clonality metrics for all donors
all_results = compute_metrics(collection, CLONALITY_METRICS)

# Convert to DataFrame
df = metrics_to_dataframe(collection, all_results)
```

## Performance Considerations

Computing only the metrics you need can improve performance, especially when:

1. Processing many repertoires
2. Running in a loop or pipeline
3. Only needing basic counts (richness, total)

```julia
# Fast: computes only 2 metrics
result = compute_metrics(rep, Richness() + TotalCount())

# Slower: computes all 14 metrics
result = compute_metrics(rep)
```

## Type Stability

The composable system is designed for type stability:

```julia
# MetricSet is parameterized by the tuple of metric types
metrics = ShannonEntropy() + Clonality()
typeof(metrics)  # MetricSet{Tuple{ShannonEntropy, Clonality}}

# compute_metric dispatches on metric type
compute_metric(rep, ShannonEntropy())  # Calls shannon_entropy(rep)
compute_metric(rep, Clonality())       # Calls clonality(rep)
```

This means the compiler can specialize code for your specific metric selection, potentially enabling optimizations.
