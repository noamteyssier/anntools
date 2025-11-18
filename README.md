# anntools

A command line utility for working with [`AnnData`](https://anndata.readthedocs.io/en/stable/) files.

`anntools` provides tools to manipulate, analyze, and transform AnnData (`.h5ad`) files commonly used in single-cell genomics. It offers efficient operations for quality control, downsampling, concatenation, pseudobulking, and more.

## Installation

`anntools` can be installed using `uv` or `pip`:

```bash
# Using uv (recommended)
uv tool install anntools-bio

# Using pip
pip install anntools-bio

# From source
git clone https://github.com/noamteyssier/anntools.git
cd anntools
uv tool install -U -e .
```

## Usage

```bash
# Get help information
anntools --help

# Get help for specific commands
anntools downsample --help
anntools concat --help
anntools info --help
```

## Commands

### `info`
Display comprehensive information about an h5ad file including shape, sparsity, layers, and metadata.

```bash
# Basic info
anntools info input.h5ad

# Detailed info with metadata summaries
anntools info input.h5ad -v
```

### `downsample`
Downsample UMIs or cells to a specified fraction.

```bash
# Downsample UMIs using binomial sampling
anntools downsample input.h5ad 0.5 -o output.h5ad

# Downsample cells
anntools downsample input.h5ad 0.8 --which cells

# Use multinomial sampling
anntools downsample input.h5ad 0.5 --method multinomial

# Set random seed for reproducibility
anntools downsample input.h5ad 0.5 --seed 42
```

### `concat`
Concatenate multiple h5ad files along the observation axis.

```bash
# Concatenate with inner join (only common variables)
anntools concat output.h5ad file1.h5ad file2.h5ad file3.h5ad

# Concatenate with outer join (all variables)
anntools concat output.h5ad file1.h5ad file2.h5ad --join outer

# Add batch labels
anntools concat output.h5ad file1.h5ad file2.h5ad --batch-key batch

# Specify custom batch labels
anntools concat output.h5ad file1.h5ad file2.h5ad --batch-key batch --batch-categories "ctrl,treat"
```

### `pseudobulk`
Aggregate single-cell data into pseudobulk profiles.
Pseudobulking is done using [`adpbulk`](https://github.com/noamteyssier/adpbulk).

```bash
# Pseudobulk by cell type
anntools pseudobulk input.h5ad cell_type -o output.h5ad

# Pseudobulk by multiple grouping variables
anntools pseudobulk input.h5ad cell_type sample_id

# Use different aggregation methods
anntools pseudobulk input.h5ad cell_type --method sum
anntools pseudobulk input.h5ad cell_type --method median

# Pseudobulk a specific layer
anntools pseudobulk input.h5ad cell_type --layer counts
```

### `qc`
Calculate quality control metrics.

```bash
# Add QC metrics to a new file
anntools qc input.h5ad -o output_qc.h5ad

# Replace existing file with QC metrics
anntools qc input.h5ad --replace
```

### `sparse`
Convert data to CSR sparse format.

```bash
# Convert to sparse format
anntools sparse input.h5ad -o output_sparse.h5ad

# Replace existing file
anntools sparse input.h5ad --replace
```

### `view-obs` / `view-var`
Export observation or variable metadata to TSV.

```bash
# View observation metadata
anntools view-obs input.h5ad > obs.tsv

# View variable metadata
anntools view-var input.h5ad > var.tsv
```

## License

MIT License - see [LICENSE](LICENSE) for details.
