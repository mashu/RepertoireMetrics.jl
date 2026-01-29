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

## Rarefaction for Comparing Uneven Samples

### The Problem with Uneven Sample Sizes

Many diversity metrics are sensitive to sample size (sequencing depth). A repertoire with 100,000 sequences will almost always appear more diverse than one with 10,000 sequences, even if the underlying distributions are identical. This is because deeper sequencing captures more rare lineages.

Metrics most affected by sample size:
- **Richness** - heavily affected (more sampling = more lineages detected)
- **Chao1** - heavily affected (relies on singletons/doubletons)
- **Shannon entropy** - moderately affected
- **Simpson index** - less affected (dominated by abundant lineages)

### What is Rarefaction?

Rarefaction is the process of randomly subsampling a repertoire to a fixed depth (total count). By subsampling all repertoires to the same depth, you can make fair comparisons.

### Basic Usage

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

### Reproducibility with Random Seed

Rarefaction involves random sampling. For reproducible results, provide a random number generator:

```julia
using Random

rng = MersenneTwister(42)  # Fixed seed
rarefied = rarefaction(rep, 10000; rng=rng)
```

### Rarefaction Curves

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

### Best Practices

1. **Choose a common depth** - Use the minimum depth across all samples you want to compare
2. **Don't rarefy too aggressively** - Very low depths lose too much information
3. **Consider multiple rarefactions** - Average results over several random subsamples for robustness
4. **Report the depth used** - Always document the rarefaction depth in your analysis
5. **Alternative: use robust metrics** - Simpson index and Berger-Parker are naturally less sensitive to sample size

### When NOT to Use Rarefaction

- When comparing samples of similar depth (within ~2-fold)
- When using metrics insensitive to depth (Simpson, Berger-Parker)
- When absolute richness (not relative comparison) is the goal

## CDR3 Length Statistics

In addition to diversity metrics, CDR3 length is an important repertoire characteristic. Longer CDR3s are associated with higher specificity and are often enriched in certain immune responses.

### Basic CDR3 Statistics

```julia
using RepertoireMetrics, CSV, DataFrames

df = CSV.read("sequences.tsv", DataFrame)

# Compute CDR3 length statistics (nucleotide sequences converted to AA length)
stats = cdr3_stats(df; cdr3_column=:cdr3)
println(stats)
# Output:
# CDR3Stats:
#   Mean length:   14.5
#   Median length: 14.0
#   Std length:    3.2
#   Range:         8 - 24
#   N sequences:   50000
```

### Weighted by Count

If your data has a count column (e.g., from UMI deduplication), weight the statistics:

```julia
stats = cdr3_stats(df; cdr3_column=:cdr3, count_column=:count)
```

### Using Amino Acid Sequences

If your CDR3 column contains amino acid sequences:

```julia
stats = cdr3_stats(df; cdr3_column=:cdr3_aa, use_aa=true)
```

### CDR3 Length Distribution

Get the full distribution of CDR3 lengths:

```julia
dist = cdr3_length_distribution(df; cdr3_column=:cdr3, count_column=:count)

# Print distribution
for (len, count) in sort(collect(dist))
    println("Length $len: $count sequences")
end
```
