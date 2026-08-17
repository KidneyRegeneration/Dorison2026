# Dorison2026

## This repository provides code to do general quantification  and colocalisation within glomerular and podocyte regions.

## Repository layout

- `segmentation/` — nuclei, glomerulus, and podocyte segmentation scripts
- `quantify/` — per-podocyte / per-glomerulus marker quantification and result aggregation
- `plotting/` — result plotting scripts
- `qc/` — QC image/PowerPoint generation
- `colocalisation/` — Manders coefficient (M1/M2) 
- `nf_workflow/` — the main Nextflow pipeline, its config, and an example samplesheet
- `standardised_pipelines/` — example Nextflow pipelines (`pipeline.nf`, `nextflow.config`, SLURM submission script) showing how the scripts above are orchestrated for specific analysis types (single/multi marker in glomeruli and/or nuclei, dual-marker colocalisation)

## Getting started

```bash
uv sync
source .venv/bin/activate
```

See `nf_workflow/example_samplesheet.csv` for the samplesheet format expected by `nf_workflow/colocalisation.nf`.
