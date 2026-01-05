# PowerAnalysis: Tumor-Informed cfDNA Variant Validation

[![License](https://img.shields.io/badge/License-Custom-blue.svg)](LICENSE)
[![R Version](https://img.shields.io/badge/R-%3E%3D4.0.0-blue)](https://www.r-project.org/)

## Overview

PowerAnalysis is an R-based toolkit for validating tumor-informed variants in cell-free DNA (cfDNA) using binomial power analysis. The package determines minimum sequencing depth requirements to distinguish true biological negatives from under-sampled sites with sufficient statistical power.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Methodology](#methodology)
- [Workflow](#workflow)
- [Usage Examples](#usage-examples)
- [Configuration](#configuration)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

## Features

- **Binomial Power Analysis**: Statistical framework for cfDNA variant validation
- **Flexible SNV Categorization**: Analyze founder, shared, and private SNVs
- **Customizable Parameters**: Adjustable detection thresholds and probability cutoffs
- **Visualization Tools**: Generate publication-ready plots
- **Batch Processing**: Handle multiple samples efficiently

## Installation

### Prerequisites

- R (≥ 4.0.0)
- Required R packages:

```r
# Install required packages
install.packages(c("tidyverse", "data.table", "ggplot2", "patchwork"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "rtracklayer"))
```

### From GitHub

```r
# Install devtools if not already installed
install.packages("devtools")

# Install PowerAnalysis
devtools::install_github("manasvitavashisth/PowerAnalysis")
```

## Quick Start

```r
library(PowerAnalysis)

run_power_analysis(
    mutect_force_calling_path = "results/mutect_force_calling/",
    output_file_path = "output/detection_rates.tsv",
   sample_file_path = "data/sample_manifest.tsv",
   snv_list = "results/pyclone/",
   min_alt_reads = 3,
   prob_cutoff = 0.8,
   snv_type = "data/cluster_categories.tsv"
 )

# Visualize results
results <- compare_detection_rates(
   track = detection_results,
   metadata = sample_info,
   variables = c("Organ", "Phenotype","CCP score")
 )
```

## Methodology

### Statistical Framework

The probability of observing *x* or more alternate reads at a tumor mutant locus in cfDNA follows a binomial distribution:

$$
P(X \geq x) = 1 - F_{\text{binom}}(x-1; N, p)
$$

Where:
- **N**: Read depth at the cfDNA locus
- **p**: Expected probability of observing the mutant allele
- **F_binom**: Cumulative binomial distribution function

### Expected Probability Calculation

$$
p = \text{VAF}_{\text{tumor}} \times \frac{\text{cfDNA Purity}}{\text{Tumor Purity}}
$$

### Power Analysis Parameters

Different hyperparameters can be used for how many alt alleles you expect to observe and with what probability as shown in the table below, we went with there should be enough read depth to observe atleast 3 or more alt alleles with a probability of 0.8.

<img width="543" height="226" alt="image" src="https://github.com/user-attachments/assets/b252b714-4fe8-4a80-a513-8f53431b14c9" />

The following table shows the relationship between minimum alternate reads and detection probability:

| Min Alt Reads | Probability | Interpretation |
|---------------|-------------|----------------|
| 1 | 0.5 | Very sensitive, higher false positives |
| 2 | 0.7 | Balanced approach |
| **3** | **0.8** | **Recommended default** |
| 4 | 0.9 | More conservative |
| 5 | 0.95 | Highly conservative |


## Workflow

### Step 1: Prepare BED Files

Generate BED files of tumor-detected mutations:

```r
source("R/prepare_bed_files.R")

prepare_bed_file(tumor_file_path,output_file_path,patient_file_path,validate_input = TRUE)
```

**SNV Categories:**
- **Founder SNVs**: Detected in all metastatic sites
- **Shared SNVs**: Detected in >1 but not all sites
- **Private SNVs**: Detected in only one site

*Note: Categorization is limited by the number of sequenced metastatic sites*

### Step 2: Force Call cfDNA Alleles

Use GATK Mutect2 for force calling at tumor-informed loci:

```bash
# Using Mutect2 (recommended for consistency with tumor calling)
gatk Mutect2 \
  -R reference.fasta \
  -I cfdna_sample.bam \
  -L tumor_mutations.bed \
  --force-active \
  -O cfdna_forced_calls.vcf

# Alternative: GATK CollectAllelicCounts (faster but fewer filters)
gatk CollectAllelicCounts \
  -R reference.fasta \
  -I cfdna_sample.bam \
  -L tumor_mutations.bed \
  -O cfdna_allelic_counts.txt
```

Use Mutect force calling to get the cfDNA allele counts at each tumor informed loci. Alternately, GATK Collect Allelic Counts can also be used, though it applies fewer filters on reads mutect force calling and hence can have different read depth count. We suggest using Mutect force calling if you use Mutect2 for mutation calling in tumors to keep consistent workflows, though the trade-off is that Mutect2 needs more compute time compared to Collect Allelic Count.

### Step 3: Run Power Analysis

```r
source("R/power_analysis.R")

run_power_analysis(
    mutect_force_calling_path = "results/mutect_force_calling/",
    output_file_path = "output/detection_rates.tsv",
   sample_file_path = "data/sample_manifest.tsv",
   snv_list = "results/pyclone/",
   min_alt_reads = 3,
   prob_cutoff = 0.8,
   snv_type = "data/cluster_categories.tsv"
 )
```

### Step 4: Downstream Analysis

```r
# Analyze detection by metastatic site
plot_by_organ_site(results, metadata = clinical_data)

# Generate summary statistics
results <- compare_detection_rates(
   track = detection_results,
   metadata = sample_info,
   variables = c("Organ", "Phenotype","CCP score")
 )
```

### Development Setup

```bash
git clone https://github.com/manasvitavashisth/PowerAnalysis.git
cd PowerAnalysis
```

## Citation

If you use PowerAnalysis in your research, please cite:

```bibtex
@software{poweranalysis2024,
  author = {Vashisth, Manasvita},
  title = {PowerAnalysis: Tumor-Informed cfDNA Variant Validation},
  year = {2024},
  url = {https://github.com/manasvitavashisth/PowerAnalysis}
}
```

## License

Copyright 2025 Fred Hutchinson Cancer Center

Permission is hereby granted, free of charge, to any government or not-for-profit entity, or to any person employed at one of the foregoing (each, an "Academic Licensee") who obtains a copy of this software and associated documentation files (the “Software”), to deal in the Software purely for non-commercial research and educational purposes, including the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or share copies of the Software, and to permit other Academic Licensees to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. 
No Academic Licensee shall be permitted to sell or use the Software or derivatives thereof in any service for commercial benefit. For the avoidance of doubt, any use by or transfer to a commercial entity shall be considered a commercial use and will require a separate license with Fred Hutchinson Cancer Center.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Contact

- **Author**: Manasvita Vashisth
- **Institution**: Fred Hutchinson Cancer Center
- **Issues**: [GitHub Issues](https://github.com/manasvitavashisth/PowerAnalysis/issues)


