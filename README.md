# Seurat

Seurat-backed module for omnibenchmark scRNA pipelines.

## Setup

```sh
pixi install
pixi run check
```

`pixi run check` loads all runtime libraries and prints `OK`. Run it after install to confirm the environment is healthy.

Run the unit tests with `pixi run -e dev test`.

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

### kNN (`knn` entrypoint)

Seurat's `FindNeighbors` on any module's embedding, written in the layout the other
kNN modules write, so any clustering module can read it.

```sh
pixi run Rscript knn.R --output_dir <dir> --name <name> \
  --embedding_tsv <embedding.tsv> --k 20 --random_seed 0
# -> <name>_neighbors.h5: distances (kNN graph, Euclidean) and connectivities (SNN graph)
#    FindNeighbors defaults: k 20, annoy (50 trees), prune.SNN 1/15; --n_dims 0 = every PC
# phases (obkit-events.jsonl): load, knn, write
```

## Conda environment export

```sh
pixi run export-env
```

Exports the resolved environment to `envs/seurat.yml`. The environment is named after the repo root folder.
