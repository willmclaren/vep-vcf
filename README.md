# vep-vcf

* **vepvcf** - a python library based on [cyvcf2](https://github.com/brentp/cyvcf2) for parsing VCF files annotated by [VEP](https://www.ensembl.org/vep)
* **filter_vcf** - a tool for filtering VEP-annotated VCF files
* **lof_rollup** - a tool for calculating per-sample gene-burden statistics from VEP-annotated VCF files


## Installation

### Pre-install requirements

Instructions for installing on a standard AWS Amazon Linux instance:

```
sudo yum update
sudo yum install -y \
  build-essential \
  libcurl4-openssl-dev \
  libssl-dev \
  libbz2-dev \
  liblzma-dev \
  zlib1g-dev
```

Install python3 with Anaconda (NB: not required if you already have python3 installed):

```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
export PATH=${PATH}:${HOME}/miniconda3/bin
conda init bash
source ~/.bashrc
```

Log out and back in or start a new shell for conda to properly initialise.


### Check out repository and install

You probably already have the repository checked out, but if not:

```
git clone git@bitbucket.org:ggcuk/vep-vcf.git
cd vep-vcf
```

Set up virtual environment for vep-vcf

```
virtualenv -p python3.7 venv
source venv/bin/activate
```

Install the package:

```
pip install -r requirements.txt
pip install -e .
```

Test installation:

```
pytest
```


## Basic usage

Basic tool instructions can be viewed by running:

```
filter_vcf --help
```

and

```
lof_rollup --help
```


## Writing filters

Both `filter_vcf` and `lof_rollup` use filters to narrow down variants of interest in a VCF file.

The main fields of the VCF along with any of the fields in the encoded VEP annotation found in the VCF INFO field can be used as filter parameters.

Filters are written as python code, with a few extras. To denote a field from the VCF, use `#[field_name]`. For example, to use `filter_vcf` to filter for variants on chromosome 1:

```
filter_vep -f '#CHROM == "1"' [vcf]
```

Filter strings are fully evaluated, so any combination of logic permissible in python can be used. For example, to filter for variants on chromosome 1 where the reference allele is "C":

```
filter_vep -f '#CHROM == "1" and #REF == "C"' [vcf]
```


### Using data in files as filters

Filters can load data from files using the `~` notation. For example, suppose you have a file with gene names, one name per line, you can filter for variants with an annotation in those genes:

```
filter_vep -f `#SYMBOL in ~gene_names.txt` [vcf]
```

Or the inverse:

```
filter_vep -f `#SYMBOL not in ~gene_names.txt` [vcf]
```


### Shortcuts

The vepvcf library comes with a few commonly used shortcuts baked in. These can be accessed by using `@[shortcut_name]`. For example, to filter for variants that are "pathogenic" or "likely pathogenic" in ClinVar:

```
filter_vep -f '@pathogenic' [vcf]
```

Shortcuts may be combined with other filters:

```
filter_vep -f '@lof and #AF < 0.01' [vcf]
```

Valid shortcuts:

* **@pathogenic** - filter for variants that are "pathogenic" or "likely pathogenic" in ClinVar. The expectation is that the VCF is annotated by VEP via the ClinVar VCF, with the annotation in the sub-field named `clinvar_CLNSIG`
* **@lof** - filter for variants that are classified as high-confidence (HC) loss-of-function by LOFTEE or have the functional effect `start_lost`


### [effect] or worse shortcuts

You may also filter for variants that are, for example, `missense_variant` or worse using the `+` notation, e.g.:

```
filter_vep -f `@missense_variant+` [vcf]
```

The following ranking is used to order the VEP effect types (worst first):

```
transcript_ablation
splice_acceptor_variant
splice_donor_variant
stop_gained
frameshift_variant
stop_lost
start_lost
transcript_amplification
inframe_insertion
inframe_deletion
protein_altering_variant
missense_variant
splice_region_variant
incomplete_terminal_codon_variant
stop_retained_variant
synonymous_variant
start_retained_variant
coding_sequence_variant
mature_miRNA_variant
5_prime_UTR_variant
3_prime_UTR_variant
non_coding_transcript_exon_variant
intron_variant
NMD_transcript_variant
non_coding_transcript_variant
upstream_gene_variant
downstream_gene_variant
TFBS_ablation
TFBS_amplification
TF_binding_site_variant
regulatory_region_ablation
regulatory_region_amplification
regulatory_region_variant
feature_elongation
feature_truncation
intergenic_variant
sequence_variant
```
