![RepertoireMetrics_logo](https://github.com/user-attachments/assets/b9234102-f8c4-49a0-89b0-688b18bf383f)

# RepertoireMetrics.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mashu.github.io/RepertoireMetrics.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mashu.github.io/RepertoireMetrics.jl/dev/)
[![Build Status](https://github.com/mashu/RepertoireMetrics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mashu/RepertoireMetrics.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mashu/RepertoireMetrics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/RepertoireMetrics.jl)

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

# Read and analyze a repertoire
rep = read_repertoire("data.tsv", VJCdr3Definition())

# Compute all metrics (default)
metrics = compute_metrics(rep)

# Or select specific metrics
metrics = compute_metrics(rep, ShannonEntropy() + Clonality() + D50())
```

## Documentation

See the [documentation](https://mashu.github.io/RepertoireMetrics.jl/stable/) for detailed usage, metric explanations with formulas, and API reference.

## License

MIT License
