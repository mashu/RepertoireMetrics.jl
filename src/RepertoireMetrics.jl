"""
    RepertoireMetrics

A Julia package for computing diversity and clonality metrics for B cell repertoire
sequencing data from MIAIRR-formatted files.

## Quick Start

```julia
using RepertoireMetrics

# Read a repertoire using V-J-CDR3 lineage definition
rep = read_repertoire("data.tsv", VJCdr3Definition())

# Compute all metrics (default)
metrics = compute_metrics(rep)
println(metrics.clonality)
println(metrics.shannon_entropy)

# Or select specific metrics
metrics = compute_metrics(rep, ShannonEntropy() + Clonality() + D50())
```

## Multi-donor Analysis

```julia
# Read multiple files
collection = read_repertoires_from_directory("data/", VJCdr3Definition())

# Compute metrics for all
all_metrics = compute_metrics(collection)

# Export to DataFrame
df = metrics_to_dataframe(collection, all_metrics)
write_metrics("results.tsv", df)
```
"""
module RepertoireMetrics

# External dependencies
using CSV
using DataFrames
using Statistics

# Include source files in dependency order
include("types.jl")
include("repertoire.jl")
include("metrics.jl")
include("metric_sets.jl")
include("io.jl")
include("cdr3_stats.jl")

# ============================================================================
# Exports
# ============================================================================

# Types - Lineage definition strategies
export AbstractLineageDefinition,
       LineageIDDefinition,
       VJCdr3Definition,
       CustomDefinition

# Types - Results
export AbstractMetricResult,
       Metrics,
       HillNumber

# Types - Data structures
export Repertoire,
       RepertoireCollection

# Core functions - Repertoire operations
export lineage_key,
       richness,
       total_count,
       frequencies,
       counts,
       lineage_ids,
       donor_id

# Core functions - Metrics computation
export compute_metrics,
       compute_metric,
       shannon_entropy,
       simpson_index,
       simpson_diversity,
       inverse_simpson,
       berger_parker_index,
       gini_coefficient,
       clonality,
       evenness,
       d50,
       chao1,
       hill_number

# Composable metric selection
export AbstractMetric,
       Richness,
       TotalCount,
       ShannonEntropy,
       ShannonDiversity,
       NormalizedShannon,
       SimpsonIndex,
       SimpsonDiversity,
       InverseSimpson,
       BergerParker,
       Evenness,
       Clonality,
       GiniCoefficient,
       D50,
       Chao1,
       MetricSet,
       # Predefined metric sets
       DIVERSITY_METRICS,
       CLONALITY_METRICS,
       RICHNESS_METRICS,
       ALL_METRICS

# Core functions - Sampling
export rarefaction

# I/O functions
export read_repertoire,
       read_repertoires,
       read_repertoires_from_directory,
       repertoire_from_dataframe,
       split_by_donor,
       metrics_to_dataframe,
       write_metrics

# Utilities
export first_allele

# CDR3 statistics
export CDR3Stats,
       cdr3_stats,
       cdr3_length_distribution

end # module RepertoireMetrics
