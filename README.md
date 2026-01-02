# Tumor-Informed cfDNA Variant Validation 🧬

This repository contains a workflow to validate tumor-informed variants in cell-free DNA (cfDNA). It applies a binomial power analysis to distinguish true biological negatives from under-sampled sites.


## 📊 Statistical Power Analysis

To establish the minimum sequencing depth required to achieve sufficient statistical power, we model the probability of observing alternate reads using a Binomial Distribution.

### Probability of Detection
The probability ($prob$) of observing $x$ or more alternate reads in cfDNA at a tumor mutant locus is:

$$
prob = 1 - pbinom(x - 1, N, p)
$$

Where:
* **pbinom**: Cumulative binomial probability distribution function.
* **N**: Read Depth of cfDNA at the mutant locus.
* **p**: Expected probability of observing the mutant allele.

### Expected Probability ($p$)
The expected probability of observing a mutant allele is calculated as:

$$
p = VAF_{\text{tumor}} \times \frac{\text{cfDNA Purity}}{\text{Tumor Purity}}
$$

## Workflow 
Use Prepare_bed_files.R to generate bed files of mutations detected in tumor sites. For our analysis, we chose to query all tumor SNVs (single nucleotide variants), and additionally categorized them as founder, shared and private SNVs. Private SNVs are SNVs which were detected only in one metastatic site in the same patient. Shared SNVs are SNVs detected in more than one but not all metastatic sites. Founder SNVs are those detected in all metastatic sites in the patient. Note: This categorization is limited by the number of metastatic sites in the patient sequenced, for example, a shared SNV can be categorized as a private SNV if the metastatic site that shares the SNV is not sequenced.
Use Mutect force calling to get the cfDNA allele counts at each tumor informed loci. Alternately, GATK Collect Allelic Counts can also be used, though it applies fewer filters on reads mutect force calling and hence can have different read depth count. We suggest using Mutect force calling if you use Mutect2 for mutation calling in tumors to keep consistent workflows, though the trade-off is that Mutect2 needs more compute time compared to Collect Allelic Count.
Then run PowerAnalysis.R to identify sites in cfDNA that have enough read depth to detect alt allele in the tumor. Different hyperparameters can be used for how many alt alleles you expect to observe and with what probability as shown in the table below, we went with there should be enough read depth to observe atleast 3 or more alt alleles with a probability of 0.8.

<img width="543" height="226" alt="image" src="https://github.com/user-attachments/assets/b252b714-4fe8-4a80-a513-8f53431b14c9" />

Further downstream custom analysis includes plotting differences in private vs shared vs founder SNV detection in cfDNA, plotting differences across metastatic organ sites etc.

## License

Copyright 2025 Fred Hutchinson Cancer Center

Permission is hereby granted, free of charge, to any government or not-for-profit entity, or to any person employed at one of the foregoing (each, an "Academic Licensee") who obtains a copy of this software and associated documentation files (the “Software”), to deal in the Software purely for non-commercial research and educational purposes, including the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or share copies of the Software, and to permit other Academic Licensees to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. 
No Academic Licensee shall be permitted to sell or use the Software or derivatives thereof in any service for commercial benefit. For the avoidance of doubt, any use by or transfer to a commercial entity shall be considered a commercial use and will require a separate license with Fred Hutchinson Cancer Center.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


THE SOFTWARE IS PROVIDED ÒAS ISÓ, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

