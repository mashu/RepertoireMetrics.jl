# I/O functions for reading MIAIRR TSV files and constructing Repertoires

using CSV
using DataFrames

# ============================================================================
# DataFrame to Repertoire conversion
# ============================================================================

"""
    repertoire_from_dataframe(
        df::DataFrame,
        strategy::AbstractLineageDefinition;
        count_column::Union{Symbol,Nothing} = :count,
        donor_id::String = "",
        donor_column::Union{Symbol,Nothing} = nothing,
        length_column::Union{Symbol,Nothing} = nothing,
        length_aa::Bool = false
    ) -> Repertoire{Int}

Convert a DataFrame to a Repertoire using the specified lineage definition strategy.

# Arguments
- `df`: Input DataFrame in MIAIRR format
- `strategy`: Lineage definition strategy (e.g., `VJCdr3Definition()`, `LineageIDDefinition()`)
- `count_column`: Column containing sequence counts (default `:count`). If `nothing` or 
  column doesn't exist, each row counts as 1.
- `donor_id`: Donor/sample identifier. If empty and `donor_column` is specified, 
  extracts from first row.
- `donor_column`: Column containing donor IDs (e.g., `:library_id`)
- `length_column`: Column to compute length statistics from (e.g., `:cdr3`, `:junction`).
  If `nothing` (default), no length statistics are computed.
- `length_aa`: If `true`, `length_column` contains amino acid sequences. If `false` (default),
  assumes nucleotide sequences and converts to amino acid length (÷3).

# Returns
A `Repertoire{Int}` with aggregated lineage counts. If `length_column` is specified,
length statistics are stored in metadata and accessible via composable metrics like
`MeanLength()`, `MedianLength()`, etc.

# Example
```julia
df = CSV.read("sequences.tsv", DataFrame)

# Basic repertoire
rep = repertoire_from_dataframe(df, VJCdr3Definition())

# With CDR3 length statistics
rep = repertoire_from_dataframe(df, VJCdr3Definition(); length_column=:cdr3)
println(mean_length(rep))  # Access via function
metrics = compute_metrics(rep, MeanLength() + MedianLength())  # Or compose

# Using amino acid column directly
rep = repertoire_from_dataframe(df, VJCdr3Definition(); 
    length_column=:cdr3_aa, length_aa=true)
```
"""
function repertoire_from_dataframe(
    df::DataFrame,
    strategy::AbstractLineageDefinition;
    count_column::Union{Symbol,Nothing} = :count,
    donor_id::String = "",
    donor_column::Union{Symbol,Nothing} = nothing,
    length_column::Union{Symbol,Nothing} = nothing,
    length_aa::Bool = false
)
    isempty(df) && return Repertoire(Int[], String[]; donor_id=donor_id)
    
    # Determine count column
    has_count = count_column !== nothing && hasproperty(df, count_column)
    
    # Extract donor_id if not provided
    if isempty(donor_id) && donor_column !== nothing && hasproperty(df, donor_column)
        donor_id = string(df[1, donor_column])
    end
    
    # Aggregate counts by lineage
    lineage_counts = Dict{Any,Int}()
    
    for row in eachrow(df)
        key = lineage_key(strategy, row)
        
        # Skip rows with missing keys
        if _has_missing(key)
            continue
        end
        
        count_val = has_count ? _get_count(row, count_column) : 1
        
        lineage_counts[key] = get(lineage_counts, key, 0) + count_val
    end
    
    # Convert to vectors
    n = length(lineage_counts)
    counts_vec = Vector{Int}(undef, n)
    ids_vec = Vector{String}(undef, n)
    
    for (i, (key, count)) in enumerate(lineage_counts)
        counts_vec[i] = count
        ids_vec[i] = _key_to_string(key)
    end
    
    metadata = Dict{String,Any}(
        "strategy" => string(typeof(strategy)),
        "original_rows" => nrow(df)
    )
    
    # Compute length statistics if requested
    if length_column !== nothing && hasproperty(df, length_column)
        length_stats = compute_length_stats(
            df;
            length_column=length_column,
            count_column=has_count ? count_column : nothing,
            use_aa=length_aa
        )
        metadata["length_stats"] = length_stats
    end
    
    return Repertoire(counts_vec, ids_vec; donor_id=donor_id, metadata=metadata)
end

# Helper to safely get count value
function _get_count(row, count_column::Symbol)
    val = getproperty(row, count_column)
    if val === missing || val === nothing
        return 1
    elseif val isa Real
        return round(Int, val)
    else
        parsed = tryparse(Int, string(val))
        return parsed === nothing ? 1 : parsed
    end
end

# Helper to check for missing values in key
_has_missing(key::Missing) = true
_has_missing(key::Nothing) = true
_has_missing(key::Tuple) = any(_has_missing, key)
_has_missing(key) = key === missing

# Convert lineage key to string representation
_key_to_string(key::AbstractString) = String(key)
_key_to_string(key::Tuple) = join(string.(key), "|")
_key_to_string(key) = string(key)

# ============================================================================
# File reading functions
# ============================================================================

"""
    read_repertoire(
        filepath::AbstractString,
        strategy::AbstractLineageDefinition;
        count_column::Union{Symbol,Nothing} = :count,
        donor_id::String = "",
        donor_column::Union{Symbol,Nothing} = nothing,
        length_column::Union{Symbol,Nothing} = nothing,
        length_aa::Bool = false,
        kwargs...
    ) -> Repertoire{Int}

Read a MIAIRR TSV file and convert to a Repertoire.

# Arguments
- `filepath`: Path to TSV file
- `strategy`: Lineage definition strategy
- `count_column`: Column containing sequence counts (default `:count`)
- `donor_id`: Donor/sample identifier
- `donor_column`: Column to extract donor ID from (e.g., `:library_id`)
- `length_column`: Column to compute length statistics from (e.g., `:cdr3`)
- `length_aa`: If `true`, `length_column` contains amino acid sequences
- `kwargs...`: Additional arguments passed to `CSV.read`

# Example
```julia
# Using lineage_id column
rep = read_repertoire("collapsed_data.tsv", LineageIDDefinition())

# Using V-J-CDR3 definition with CDR3 length stats
rep = read_repertoire("sequences.tsv", VJCdr3Definition(); length_column=:cdr3)

# With donor ID from file
rep = read_repertoire("data.tsv", VJCdr3Definition(); donor_column=:library_id)

# Compute length from amino acid column
rep = read_repertoire("data.tsv", VJCdr3Definition(); 
    length_column=:cdr3_aa, length_aa=true)
```
"""
function read_repertoire(
    filepath::AbstractString,
    strategy::AbstractLineageDefinition;
    count_column::Union{Symbol,Nothing} = :count,
    donor_id::String = "",
    donor_column::Union{Symbol,Nothing} = nothing,
    length_column::Union{Symbol,Nothing} = nothing,
    length_aa::Bool = false,
    kwargs...
)
    # Determine file type and set appropriate CSV options
    csv_kwargs = Dict{Symbol,Any}(
        :missingstring => ["", "NA", "na", "N/A", "n/a"],
        :delim => '\t',
    )
    merge!(csv_kwargs, Dict(kwargs))
    
    # Handle gzipped files
    if endswith(lowercase(filepath), ".gz")
        error("Gzipped files require CodecZlib.jl. Please decompress first or add CodecZlib.jl support.")
    end
    
    df = CSV.read(filepath, DataFrame; csv_kwargs...)
    
    # Extract donor_id from filename if not provided
    if isempty(donor_id) && donor_column === nothing
        donor_id = _extract_donor_from_filename(filepath)
    end
    
    rep = repertoire_from_dataframe(
        df, strategy;
        count_column=count_column,
        donor_id=donor_id,
        donor_column=donor_column,
        length_column=length_column,
        length_aa=length_aa
    )
    
    # Add filepath to metadata
    rep.metadata["filepath"] = filepath
    
    return rep
end

function _extract_donor_from_filename(filepath::AbstractString)
    basename_str = basename(filepath)
    # Remove common extensions
    for ext in [".tsv.gz", ".tsv", ".csv.gz", ".csv", ".txt"]
        if endswith(lowercase(basename_str), ext)
            basename_str = basename_str[1:end-length(ext)]
            break
        end
    end
    return basename_str
end

"""
    read_repertoires(
        filepaths::Vector{<:AbstractString},
        strategy::AbstractLineageDefinition;
        kwargs...
    ) -> RepertoireCollection{Int}

Read multiple MIAIRR TSV files and return a collection.

# Example
```julia
files = ["donor1.tsv", "donor2.tsv", "donor3.tsv"]
collection = read_repertoires(files, VJCdr3Definition())
all_metrics = compute_metrics(collection)
```
"""
function read_repertoires(
    filepaths::Vector{<:AbstractString},
    strategy::AbstractLineageDefinition;
    kwargs...
)
    repertoires = [read_repertoire(fp, strategy; kwargs...) for fp in filepaths]
    return RepertoireCollection(repertoires)
end

"""
    read_repertoires_from_directory(
        dirpath::AbstractString,
        strategy::AbstractLineageDefinition;
        pattern::Regex = r"\\.tsv\$"i,
        kwargs...
    ) -> RepertoireCollection{Int}

Read all matching files from a directory and return a collection.

# Arguments
- `dirpath`: Path to directory containing TSV files
- `strategy`: Lineage definition strategy
- `pattern`: Regex pattern to match filenames (default: `.tsv` files)
- `kwargs...`: Additional arguments passed to `read_repertoire`

# Example
```julia
collection = read_repertoires_from_directory("data/", VJCdr3Definition())
```
"""
function read_repertoires_from_directory(
    dirpath::AbstractString,
    strategy::AbstractLineageDefinition;
    pattern::Regex = r"\.tsv$"i,
    kwargs...
)
    isdir(dirpath) || throw(ArgumentError("Not a directory: $dirpath"))
    
    files = filter(f -> occursin(pattern, f), readdir(dirpath))
    filepaths = [joinpath(dirpath, f) for f in files]
    
    isempty(filepaths) && @warn "No files matching pattern found in $dirpath"
    
    return read_repertoires(filepaths, strategy; kwargs...)
end

# ============================================================================
# Multi-donor DataFrame processing
# ============================================================================

"""
    split_by_donor(
        df::DataFrame,
        donor_column::Symbol,
        strategy::AbstractLineageDefinition;
        kwargs...
    ) -> RepertoireCollection{Int}

Split a multi-donor DataFrame into separate Repertoires by donor.

# Arguments
- `df`: DataFrame containing data from multiple donors
- `donor_column`: Column containing donor identifiers
- `strategy`: Lineage definition strategy
- `kwargs...`: Additional arguments passed to `repertoire_from_dataframe`

# Example
```julia
df = CSV.read("all_donors.tsv", DataFrame)
collection = split_by_donor(df, :library_id, VJCdr3Definition())
```
"""
function split_by_donor(
    df::DataFrame,
    donor_column::Symbol,
    strategy::AbstractLineageDefinition;
    kwargs...
)
    hasproperty(df, donor_column) || 
        throw(ArgumentError("Column $donor_column not found in DataFrame"))
    
    donor_groups = groupby(df, donor_column)
    
    repertoires = Repertoire{Int}[]
    for group in donor_groups
        donor_id = string(first(group[!, donor_column]))
        rep = repertoire_from_dataframe(
            DataFrame(group), strategy;
            donor_id=donor_id, kwargs...
        )
        push!(repertoires, rep)
    end
    
    return RepertoireCollection(repertoires)
end

# ============================================================================
# Export functions
# ============================================================================

"""
    write_metrics(filepath::AbstractString, df::DataFrame; kwargs...)

Write metrics DataFrame to a TSV file.
"""
function write_metrics(filepath::AbstractString, df::DataFrame; kwargs...)
    CSV.write(filepath, df; delim='\t', kwargs...)
end
