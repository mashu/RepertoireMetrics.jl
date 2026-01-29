# Quick Start

This guide covers the essential workflows for analyzing B cell repertoire diversity.

## Reading Repertoire Data

RepertoireMetrics.jl reads MIAIRR-formatted TSV files. You need to specify how lineages are defined.

### Using `lineage_id` Column

If your data has been pre-processed with lineage assignment (e.g., by LineageCollapse.jl):

```julia
using RepertoireMetrics

rep = read_repertoire("collapsed_data.tsv", LineageIDDefinition())
```

### Using V-J-CDR3 Combination

Define lineages by V gene, J gene, and CDR3 sequence (compatible with LineageCollapse.jl's default behavior):

```julia
rep = read_repertoire("sequences.tsv", VJCdr3Definition())
```

The `VJCdr3Definition` automatically uses the first allele from gene calls. You can customize this:

```julia
# Use full gene call (with all alleles)
rep = read_repertoire("data.tsv", VJCdr3Definition(use_first_allele=false))

# Use different column names
strategy = VJCdr3Definition(
    v_column = :v_gene,
    j_column = :j_gene,
    cdr3_column = :cdr3_aa
)
```

### Handling Counts

By default, the package looks for a `count` column. If not present, each row is counted as 1:

```julia
# Explicit: use count column
rep = read_repertoire("data.tsv", VJCdr3Definition(); count_column=:count)

# Each row = 1 sequence
rep = read_repertoire("data.tsv", VJCdr3Definition(); count_column=nothing)
```

## Computing Metrics

### All Metrics (Default)

```julia
metrics = compute_metrics(rep)
println(metrics)
```

Output:
```
Metrics (14 computed):
  richness: 1523
  total_count: 50000
  shannon_entropy: 5.2341
  ...
  clonality: 0.2854
  ...
```

### Accessing Individual Values

```julia
println("Clonality: ", metrics.clonality)
println("Shannon entropy: ", metrics.shannon_entropy)
println("D50: ", metrics.d50)
```

### Computing Specific Metrics

Select specific metrics using the `+` operator:

```julia
# Define which metrics to compute
selected = ShannonEntropy() + Clonality() + D50() + Chao1()

# Compute only those
result = compute_metrics(rep, selected)
println(result.shannon_entropy)
println(result.clonality)
```

Use predefined sets:

```julia
# Diversity-focused metrics
result = compute_metrics(rep, DIVERSITY_METRICS)

# Clonality-focused metrics  
result = compute_metrics(rep, CLONALITY_METRICS)
```

### Direct Metric Functions

For single metrics, you can also call functions directly:

```julia
H = shannon_entropy(rep)
C = clonality(rep)
d = d50(rep)
```

## Multi-Donor Analysis

### Reading Multiple Files

```julia
# From a list of files
files = ["donor1.tsv", "donor2.tsv", "donor3.tsv"]
collection = read_repertoires(files, VJCdr3Definition())

# From a directory
collection = read_repertoires_from_directory("data/", VJCdr3Definition())
```

### Splitting a Multi-Donor File

```julia
using CSV, DataFrames

df = CSV.read("all_donors.tsv", DataFrame)
collection = split_by_donor(df, :library_id, VJCdr3Definition())
```

### Computing Metrics for All Donors

```julia
# Compute all metrics
all_metrics = compute_metrics(collection)

# Or selected metrics
all_metrics = compute_metrics(collection, CLONALITY_METRICS)
```

### Exporting Results

```julia
# Convert to DataFrame
df_results = metrics_to_dataframe(collection, all_metrics)

# Write to TSV
write_metrics("results.tsv", df_results)
```

## Custom Lineage Definition

For advanced use cases, define your own grouping:

```julia
# Group by V gene family only
strategy = CustomDefinition(row -> begin
    v = string(row.v_call)
    # Extract family (e.g., "IGHV1" from "IGHV1-2*01")
    m = match(r"^(IGH?V\d+)", v)
    return m === nothing ? v : m.captures[1]
end)

rep = read_repertoire("data.tsv", strategy)
```

## Hill Numbers

For a unified view of diversity across different orders:

```julia
# q=0: Richness
h0 = hill_number(rep, 0)

# q=1: exp(Shannon entropy)
h1 = hill_number(rep, 1)

# q=2: Inverse Simpson
h2 = hill_number(rep, 2)

# Arbitrary order
h_half = hill_number(rep, 0.5)
```

## Comparing Samples of Different Depths

### The Problem with Uneven Sample Sizes

A repertoire with 100,000 sequences will almost always show more unique lineages than one with 10,000 sequences, even if the underlying distributions are identical. This is because deeper sequencing captures more rare clones.

### Solution 1: Use Depth-Robust Metrics (Recommended Default)

Many metrics work on **frequencies** (proportions) rather than counts. Since frequencies sum to 1.0 regardless of total count, these metrics are mathematically less sensitive to depth:

```julia
# Recommended: use the predefined ROBUST_METRICS set
metrics = compute_metrics(rep, ROBUST_METRICS)

println("Depth: ", metrics.depth)  # Always report sequencing depth!
println("Simpson: ", metrics.simpson_diversity)
println("Clonality: ", metrics.clonality)
```

Or select specific robust metrics:

```julia
metrics = compute_metrics(rep, Depth() + SimpsonDiversity() + Clonality() + GiniCoefficient())
```

| Metric | Depth Sensitivity |
|--------|-------------------|
| Simpson Index/Diversity | **Low** |
| Inverse Simpson | **Low** |
| Berger-Parker | **Very Low** |
| Gini Coefficient | **Low** |
| Shannon Entropy | Moderate |
| Richness | **Very High** |

**Key insight:** If you have accurate count data (e.g., from UMI deduplication), the frequencies computed from your counts already reflect the true clonal distribution. Simpson-family metrics on these frequencies are directly comparable across samples.

### Solution 2: Rarefaction (For Richness Comparisons)

When you specifically need to compare **richness** (number of unique clones), rarefaction is the standard approach. It randomly subsamples all repertoires to the same depth.

Metrics most affected by sample size (where rarefaction helps):
- **Richness** - heavily affected (more sampling = more lineages detected)
- **Chao1** - heavily affected (relies on singletons/doubletons)
- **Shannon entropy** - moderately affected

#### What is Rarefaction?

Rarefaction randomly subsamples a repertoire to a fixed depth (total count). By subsampling all repertoires to the same depth, you can fairly compare richness.

#### Basic Usage

```julia
# Check sample sizes
println("Rep1: ", total_count(rep1))  # e.g., 50,000
println("Rep2: ", total_count(rep2))  # e.g., 12,000

# Rarefy both to the smaller size
target_depth = min(total_count(rep1), total_count(rep2))
rep1_rarefied = rarefaction(rep1, target_depth)
rep2_rarefied = rarefaction(rep2, target_depth)

# Now compare fairly
metrics1 = compute_metrics(rep1_rarefied)
metrics2 = compute_metrics(rep2_rarefied)
```

#### Reproducibility with Random Seed

Rarefaction involves random sampling. For reproducible results, provide a random number generator:

```julia
using Random

rng = MersenneTwister(42)  # Fixed seed
rarefied = rarefaction(rep, 10000; rng=rng)
```

#### Rarefaction Curves

To understand how metrics change with sequencing depth, compute metrics at multiple depths:

```julia
depths = [1000, 5000, 10000, 25000, 50000]
results = []

for d in depths
    if d <= total_count(rep)
        rarefied = rarefaction(rep, d)
        m = compute_metrics(rarefied, Richness() + ShannonEntropy())
        push!(results, (depth=d, richness=m.richness, shannon=m.shannon_entropy))
    end
end
```

#### Best Practices for Rarefaction

1. **Choose a common depth** - Use the minimum depth across all samples you want to compare
2. **Don't rarefy too aggressively** - Very low depths lose too much information
3. **Consider multiple rarefactions** - Average results over several random subsamples for robustness
4. **Report the depth used** - Always document the rarefaction depth in your analysis
5. **Alternative: use robust metrics** - Simpson index and Berger-Parker are naturally less sensitive to sample size

#### When NOT to Use Rarefaction

- When comparing samples of similar depth (within ~2-fold)
- When using depth-robust metrics (Simpson, Berger-Parker, Gini)
- When absolute richness (not relative comparison) is the goal
- When you want to preserve information from deeply sequenced samples

### Summary: Depth Strategy Decision Tree

```
Do you need to compare RICHNESS specifically?
├── YES → Use rarefaction, then compare richness
└── NO → Use depth-robust metrics (Simpson, Berger-Parker, Gini)
         These work directly on frequencies and don't need rarefaction
```

## Length Metrics (CDR3 and more)

Sequence length is an important repertoire characteristic. CDR3 length in particular is associated with binding specificity. RepertoireMetrics provides **composable length metrics** that integrate seamlessly with diversity metrics.

### Enabling Length Statistics

Specify a `length_column` when creating your repertoire:

```julia
using RepertoireMetrics

# Read repertoire with CDR3 length statistics
rep = read_repertoire("sequences.tsv", VJCdr3Definition(); 
    length_column=:cdr3)  # Default: nucleotide → amino acid conversion

# Or from DataFrame
df = CSV.read("sequences.tsv", DataFrame)
rep = repertoire_from_dataframe(df, VJCdr3Definition();
    count_column=:count,
    length_column=:cdr3)
```

### Accessing Length Statistics

```julia
# Direct accessors
println("Mean CDR3 length: ", mean_length(rep))
println("Median: ", median_length(rep))
println("Std: ", std_length(rep))
println("Range: ", min_length(rep), " - ", max_length(rep))

# Get full stats object
stats = length_stats(rep)
println(stats)
```

### Composable Length Metrics

Length metrics work with the same `+` operator as diversity metrics:

```julia
# Just length metrics
metrics = compute_metrics(rep, MeanLength() + MedianLength())

# Mix with diversity metrics
metrics = compute_metrics(rep, 
    Richness() + ShannonEntropy() + MeanLength() + StdLength())

# All length metrics at once
metrics = compute_metrics(rep, LENGTH_METRICS)
```

### Flexible Column Selection

The `length_column` parameter accepts any column name, making it useful for:
- CDR3 nucleotide (`:cdr3`) - converted to AA length
- CDR3 amino acid (`:cdr3_aa`) - use with `length_aa=true`
- Junction (`:junction`)
- Any other sequence column

```julia
# Using amino acid column directly
rep = read_repertoire("data.tsv", VJCdr3Definition();
    length_column=:cdr3_aa,
    length_aa=true)

# Using junction length
rep = read_repertoire("data.tsv", VJCdr3Definition();
    length_column=:junction)
```

### Length Distribution

For detailed analysis, get the full distribution:

```julia
dist = length_distribution(df; length_column=:cdr3, count_column=:count)

for (len, cnt) in sort(collect(dist))
    println("Length $len: $cnt sequences")
end
```

### Multi-Donor Analysis with Length

```julia
# Read collection with length stats
files = ["donor1.tsv", "donor2.tsv"]
collection = read_repertoires(files, VJCdr3Definition();
    length_column=:cdr3)

# Compute diversity + length metrics for all donors
results = compute_metrics(collection, 
    Richness() + ShannonEntropy() + MeanLength())

# Convert to DataFrame for analysis
df = metrics_to_dataframe(collection, results)
```
