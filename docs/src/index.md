# RepertoireMetrics.jl

A Julia package for computing diversity and clonality metrics for B cell repertoire sequencing data.

## Overview

RepertoireMetrics.jl provides a comprehensive toolkit for analyzing the clonal composition of adaptive immune repertoires. It reads MIAIRR-formatted TSV files and computes standard diversity and clonality metrics used in immunology research.

### What are Diversity and Clonality?

**Diversity** and **clonality** are two perspectives on the same underlying phenomenon: the distribution of lineages (clones) in a B cell repertoire.

- **Diversity** measures how varied a repertoire is. High diversity indicates many unique lineages with relatively even abundances (a *polyclonal* repertoire).
  
- **Clonality** measures the opposite: how dominated a repertoire is by a few lineages. High clonality indicates oligoclonal or monoclonal expansion.

Mathematically, clonality is typically defined as:

```math
\text{Clonality} = 1 - \frac{H}{\log(S)}
```

where ``H`` is Shannon entropy and ``S`` is richness. This means clonality and normalized diversity sum to 1.

### When to Use Which?

Both metrics describe the same distribution, so the choice is often a matter of interpretation:

| Metric | Interpretation | High value means |
|--------|----------------|------------------|
| Diversity | Variability | Polyclonal, healthy repertoire |
| Clonality | Dominance | Oligoclonal expansion, possible antigen-driven response |

In practice, researchers often report both, or choose based on their biological question:
- Studying immune reconstitution? Focus on **diversity**.
- Investigating antigen-specific responses? Focus on **clonality**.

## Features

- **Flexible lineage definition**: Use `lineage_id`, V-J-CDR3 combinations, or custom strategies
- **Comprehensive metrics**: Shannon, Simpson, Gini, Hill numbers, Chao1, D50, and more
- **Composable metric selection**: Choose which metrics to compute with the `+` operator
- **Type-stable design**: Proper Julia abstractions for performance
- **Multi-donor support**: Process single files or entire directories

## Installation

```julia
using Pkg
Pkg.add("RepertoireMetrics")
```

## Package Contents

```@contents
Pages = ["quickstart.md", "metrics.md", "composable.md", "api.md"]
Depth = 2
```
