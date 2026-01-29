# API Reference

## Lineage Definition Strategies

```@docs
AbstractLineageDefinition
LineageIDDefinition
VJCdr3Definition
CustomDefinition
lineage_key
first_allele
```

## Data Structures

```@docs
Repertoire
RepertoireCollection
richness
total_count
frequencies
counts
lineage_ids
donor_id
```

## Reading Data

```@docs
read_repertoire
read_repertoires
read_repertoires_from_directory
repertoire_from_dataframe
split_by_donor
```

## Computing Metrics

### Main Functions

```@docs
compute_metrics
Metrics
compute_metric
```

### Individual Metric Functions

```@docs
shannon_entropy
simpson_index
simpson_diversity
inverse_simpson
berger_parker_index
gini_coefficient
clonality
evenness
d50
chao1
hill_number
```

### Hill Numbers

```@docs
HillNumber
```

### Composable Metric Selection

```@docs
AbstractMetric
MetricSet
```

### Predefined Metric Sets

```@docs
ALL_METRICS
DIVERSITY_METRICS
CLONALITY_METRICS
RICHNESS_METRICS
ROBUST_METRICS
```

### Metric Types

```@docs
Richness
TotalCount
Depth
ShannonEntropy
ShannonDiversity
NormalizedShannon
SimpsonIndex
SimpsonDiversity
InverseSimpson
BergerParker
Evenness
Clonality
GiniCoefficient
D50
Chao1
```

## Sampling

```@docs
rarefaction
```

## Exporting Results

```@docs
metrics_to_dataframe
write_metrics
```

## Length Statistics

### Types and Functions

```@docs
LengthStats
compute_length_stats
length_distribution
has_length_stats
length_stats
```

### Composable Length Metrics

```@docs
MeanLength
MedianLength
StdLength
MinLength
MaxLength
mean_length
median_length
std_length
min_length
max_length
LENGTH_METRICS
```

## Index

```@index
```
