![RepertoireMetrics_logo](https://github.com/user-attachments/assets/b9234102-f8c4-49a0-89b0-688b18bf383f)

# RepertoireMetrics.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mashu.github.io/RepertoireMetrics.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mashu.github.io/RepertoireMetrics.jl/dev/)
[![Build Status](https://github.com/mashu/RepertoireMetrics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mashu/RepertoireMetrics.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mashu/RepertoireMetrics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/RepertoireMetrics.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/mashu/RepertoireMetrics.jl/blob/main/LICENSE)

A Julia package for computing diversity and clonality metrics for B cell repertoire sequencing data from MIAIRR-formatted files.

## Features

- **Flexible lineage definition**: Use `lineage_id`, V-J-CDR3 combinations (LineageCollapse.jl compatible), or custom strategies
- **Comprehensive metrics**: Shannon, Simpson, Gini, Hill numbers, Chao1, D50, and more
- **Composable metric selection**: Choose which metrics to compute with `+` operator
- **Type-stable design**: Proper Julia abstractions for performance
- **Multi-donor support**: Process single files or entire directories

## Installation

```julia
using Pkg
Pkg.add("RepertoireMetrics")
```

## Quick Start

```julia
using RepertoireMetrics

# Read a repertoire file (supports .tsv and .tsv.gz)
rep = read_repertoire("sample-001.tsv.gz", VJCdr3Definition(); 
                      length_column=:cdr3)
```

```
Repertoire{Int64}:
  Donor:      sample-001
  Lineages:   15420
  Total count: 198350
  Top 5 lineages:
    IGHV3-23*01|IGHJ4*02|GCGAGAGATCTTGACTACTGGGGCCAGGGAACC: 856 (0.43%)
    IGHV4-39*01|IGHJ5*02|TGTGCGAGAGTCGATTACTATGATAGTAGTGGT: 724 (0.37%)
    IGHV1-69*01|IGHJ3*02|GCGAGAGATAGTGGCTACGATTTTGACTACTGG: 512 (0.26%)
    IGHV5-51*01|IGHJ4*02|TGTGCGAGACATATTGTGGTGGTAACTGCCCC: 398 (0.2%)
    IGHV3-48*01|IGHJ6*02|GCGAGAGGGGATAGCAGCAGCTGGTACTTTGAC: 287 (0.14%)
    ... and 15415 more
```

```julia
# Compute selected metrics using the composable + operator
metrics = compute_metrics(rep, Depth() + SimpsonDiversity() + Clonality() + GiniCoefficient())
```

```
Metrics (4 computed):
  depth: 198350.0
  simpson_diversity: 0.9997
  clonality: 0.0612
  gini_coefficient: 0.548
```

```julia
# Or compute all available metrics (default)
metrics = compute_metrics(rep)
```

