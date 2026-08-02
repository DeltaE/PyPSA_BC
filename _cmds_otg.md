uv sync --extra dev

```
uv run snakemake \
  --snakefile workflow/Snakefile \
  --cores 1 \
  fetch_base_network_validation
```

```
uv run snakemake \
  --snakefile workflow/Snakefile \
  --cores 1 \
  base_network_preparation
```

```
uv run snakemake \
  --snakefile workflow/Snakefile \
  --cores 1 \
  validate_base_network
```