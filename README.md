# soydbean-variantCalling-pipeline

## Usage

Clone this repo to local
``` bash
git clone https://github.com/yeah-zin/soydbean-variantCalling-pipeline.git
```

Create a new environment with conda

```bash
conda create -n vcf_snake snakemake bcftools samtools bwa sambamba
# atviate
conda activate vcf_snake
```

Change sample path, reference path, and etc. in Snakefile

```bash
BASE_DIR = "./"
REF_NAME = "glyma.Wm82.gnm1.FCtY.genome_main.fna"
REF = BASE_DIR + "ref/" + REF_NAME
THREAD=15
DEPTH=5
```

Run command

```bash
snakemake
```
