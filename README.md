# lra-stats

```bash
Quickly compute long read stats for alignment files from stdin

Usage: lra-stats [OPTIONS]

Options:
  -o, --out-prefix <OUT_PREFIX>  out prefix string to use for naming output files [default: lra-stats]
  -h, --help                     Print help
  -V, --version                  Print version
```

## Usage

```bash
export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
samtools view --bam \
  --reference gs://nygc-comp-s-1856-resources/reference_genomes/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  gs://nygc-comp-s-1856-resources/test_small_datasets/hg002.pass.chr2:50000000-51000000.cram | lra-stats -o out_prefix
```
