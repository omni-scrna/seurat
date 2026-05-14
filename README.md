# Seurat

Seurat-backed module for omnibenchmark scRNA pipelines.

## Setup

```sh
pixi install
pixi run check
```

`pixi run check` loads all runtime libraries and prints `OK`. Run it after install to confirm the environment is healthy.

## Usage

### Gene selection (`select` entrypoint)

Selects highly variable genes from normalized expression data using Seurat's VST method.

```sh
pixi run Rscript select.R \
  --output_dir <dir> \
  --name <name> \
  --normalized.h5 <normalized.h5> \
  --rawdata.h5ad <rawdata.h5ad> \
  --filtered.cellids <cellids.txt.gz> \
  --selection_type <seurat_vst|seurat_vst_batch> \
  --number_selected <int> \
  --batch_variable <obs_column>   # required for seurat_vst_batch only
```

Output: `<output_dir>/<name>_normalized_selected.h5`

**Selection types:**
- `seurat_vst` — VST on all cells jointly (`FindVariableFeatures`)
- `seurat_vst_batch` — VST run per batch, features aggregated with `SelectIntegrationFeatures`

## Conda environment export

```sh
pixi run export-env
```

Exports the resolved environment to `envs/seurat.yml`. The environment is named after the repo root folder.
