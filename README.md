# Challenges of estimating sex-biased admixture from X chromosomal and autosomal data 

### Table of Contents
- [Models for estimating sex ratios assuming single pulse or constant admixture](#models-for-estimating-sex-ratios-assuming-single-pulse-or-constant-admixture)
   - [Software requirements](#software-requirements)
   - [Input files](#input-files)
   - [Analysis](#analysis)
   - [Output files](#output-files)
   - [Example usage](#example-usage)
- [Recalculation of sex ratios based on summary statistics from Micheletti et al., calculation of possible sex ratios given the unassigend ancestry in their study, and sensitivity analysis of the applied models](#recalculation-of-sex-ratios-based-on-summary-statistics-from-micheletti-et-al-calculation-of-possible-sex-ratios-given-the-unassigend-ancestry-in-their-study-and-sensitivity-analysis-of-the-applied-models)
   - [Software requirements](#software-requirements-1)
   - [Input files](#input-files-1)
   - [Analysis](#analysis-1)
   - [Output files](#output-files-1)
   - [Example usage](#example-usage-1)
- [Assessing the robustness of single admixture models to violations of demographic assumptions using simulations](#assessing-the-robustness-of-single-admixture-models-to-violations-of-demographic-assumptions-using-simulations)
   - [Software requirements](#software-requirements-2)
   - [Analysis](#analysis-2)
   - [Output files](#output-files-2)
   - [Example usage](#example-usage-2)
- [References](#references)

### Models for estimating sex ratios assuming single pulse or constant admixture

#### Software requirements
The code was tested with the following python package version:

 - matplotlib v3.34
 - numpy v1.20.2
 - pandas v1.2.4
 - scipy v1.6.2
 
To set up a conda environment with the following versions and activate the environment, run the following commands:
```
conda env create -f sex_ratios.yml
conda activate sex_ratios
```
 
#### Input files

Either:

1. A comma-separated file, containing the ancestry proportions of various source populations and fraction of sampled females in different regions. The first row of the file is the header, and the file has the following structure: &nbsp;  
&nbsp;  
,H<sub>1</sub><sup>X</sup>,H<sub>1</sub><sup>A</sup>,H<sub>2</sub><sup>X</sup>,H<sub>2</sub><sup>A</sup>,...,p<sub>f</sub><br>
region1,x,a,xx,aa,...,f<br>
region2,x,a,xx,aa,...,f<br>
&nbsp;  
where H<sub>i</sub><sup>X</sup> and H<sub>i</sub><sup>A</sup> are the X chromosomal and autosomal ancestry proportions that trace to source population 1, and p<sub>f</sub> are the fraction of females in the sample in a given region. Comments are indicated by a "#".

Or (only accepted by <code>scripts/estimate_sex_ratios_single_pulse.py</code>):

2. <code>.Q</code> files for autosomes and X chromosomes, which were generated with ADMIXTURE.<sup>1</sup> Additionally, <code>.fam</code> file generated with plink2 needs to be provided to get individual IDs.<sup>2</sup>

If a comma-separated file with mean ancestry proprotions is provided as input, sex ratios are calculated based on the provided proportions. If .Q and .fam files are provided, individual ancestry proportions are bootstrapped and for each bootstrap resample mean ancestry proportions are computed and sex ratios are estimated, yielding confidence intervals. 

#### Analysis
If `scripts/estimate_sex_ratios_single_pulse.py` is used, sex ratios are calculated using a dynamic model proposed by Golderberg and Rosenberg (2015) and an equilibirium model.<sup>3</sup> Both models assume a single admixture event, a constant population size, no subsequent gene flow, random mating, no genetic drift, and no selection. The `scripts/estimate_sex_ratios_constant_admixture.py` implements another model proposed by Goldberg and Rosenberg that assumes constant, nonzero admixture.<sup>3</sup> This model attempts to find parameter combinations of constant migration rates that fit the data. For details on any of these models see either Goldberg and Rosenberg or our publication.<sup>3</sup>

#### Output files

* Single pulse: Generates an excel file with estimated sex ratios for 2-15 generations of random within the admixed population and the equilibrium sex ratio. If a comma-separated file is provided this is done for each ancestry and each region. The output file name can be specified with the `-o` flag.
* Constant admixture: Generates a box plot of possible sex ratios for each ancestry and region. The name of the output file can be specified with the `-o` flag.

#### Example usage
1. To infer sex ratios assuming a single admixture event and provided mean ancestry proportions in `input.txt`, do:
```
scripts/estimate_sex_ratios_single_pulse.py -i input.txt -o output.xlsx
```

2. To infer sex ratios with confidence intervals assuming a single admixture event and provided individual level ADMIXTURE results for *K* clusters, do:
```
scripts/estimate_sex_ratios_single_pulse.py -qa admixture_results_autosomes.K.Q -qx admixture_results_x_chromosome.K.Q -f plink.fam -b 10000 -o output.xlsx
```
where confidence intervals are obtained from 10,000 bootstrap resamples(`-b`).

3. To infer sex ratios under a demographic scenarios of constant, nonzero admixture and given mean ancestry proportions provided in `input.txt`, do:
```
scripts/estimate_sex_ratios_constant_admixture.py -i input.txt -n 0.02 -d 0.01 -g 15 -o output_figure.pdf --sf0 0.5 --sm0 0.5
```
where 0.02 increments are used for the gridsearch (`-n`), accept only constant admixture parameter that have an Eucleadean distance less than 0.01 to observed mean ancestry proportions (`-d`), admixture is assumed to have occured 15 generations ago (`-g`), the output figure is saved to output_figure.pdf (`-o`), and initial female and male contributions are assumed to have been 0.5 (`--sf0` and `--sm0`; the effect of initial contributions is erased if g > 10).

Figure S1 in our manuscript can be reproduced by replacing `input.txt` with `data/ancestry_proportions_micheletti.txt` in the above command.

### Recalculation of sex ratios based on summary statistics from Micheletti et al., calculation of possible sex ratios given the unassigend ancestry in their study, and sensitivity analysis of the applied models

#### Software requirements

See above.

#### Input files

African, European, and Native American mean ancestry proportions for different geographic regions were inferred from summary statistics of Micheletti et al. We weighted mean ancestry estimates for each granular region reported in Table S10 in Micheletti et al. by the corresponding sample size reported in Table S2 in Micheletti et al.<sup>4</sup> The ancestry proportions were stored in `data/ancestry_proportions_micheletti.txt`.

#### Analysis

First, the script `sensitivity_analysis_and_distribute_unassigned_ancestry.py` applies the single pulse admixture model from Goldberg and Rosenberg (2015) to the provided ancestry proportions.<sup>3</sup> If unassigned ancestry is present, such as in the study by Micheletti et al.,<sup>4</sup> the unassigned ancestry is distributed in all permissible ways such that ancestry proportions sum to 100%. For each possible distribution, sex ratios are computed, yielding a distribution of sex ratios that is plotted as a box plot for each ancestry and region (this is only done for regions specified with the `-p` flag). In addition, it conducts a sensitivity analysis of the equilibirium model, which also holds for Goldberg's and Rosenberg's model. For more details, see our publication.

#### Output files

1. An excel file with estimated sex ratios for 2-15 generations of random within the admixed population and the equilibrium sex ratio for each ancestry and region based on provided ancestry proportions. The name of the output file can be specified with the `-o` flag (Table 1 in our manuscript).
2. A figure with possible sex ratios under consideration of the unassigned ancestry. The name of the figure can be specified with `--figure` flag (Figure 1 in our manuscript). 
3. A figure illustrating the sensitivity of the applied models. The name of the figure can be specified with the `--figure_sensitivity`flag (Figure 2 in our manuscript).

#### Example usage
The following command takes mean ancestry proportions from Micheletti et al., distributes all unassigned ancestry for the geographic regions of the Latin Caribbean, central South America, northern South America, and Central America, and saves box plots of corresponding distributions of possible s<sup>f</sup>/s<sup>m</sup> values to figure1.pdf. Additionally, it conducts a sensitivity analysis, which will be plotted to figure2.pdf. Point estimates of s<sup>f</sup>/s<sup>m</sup> for each ancestry and georgraphic regions based on the reported mean ancestry proportions will be written to sex_ratios_micheletti.xlsx.

```
scripts/sensitivity_analysis_and_distribute_unassigned_ancestry.py -i data/ancestry_proportions_micheletti.txt --figure figure1.pdf -p "Latin Caribbean" "C South America" "N South America" "Central America" --figure_sensitivity figure2.pdf -o sex_ratios_micheletti.xlsx
```

### Assessing the robustness of single admixture models to violations of demographic assumptions using simulations

#### Software requirements

The following software is required:

- ADMIXTURE v1.3.0
- plink2 v2.00a2.3
- bcftools v
- tabix v1.11
- SliM3 v3.7.1

Additionally, the following Python packages are required:

- numpy v1.22.2
- pandas v1.4.1
- multiprocessing-loggin v0.3.1
- pyslim v0.700
- msprime v1.1.1
- matplotlib v3.5.1
- scipy v1.8.0
- scikit-learn v1.0.2

The simulations and subsequent analysis are implemented in a Snakemake workflow. Therefore, snakemake needs to be installed. Snakemake will automatically install all required software from the .yml files provided in the envs directory. Snakemake can be installed from the provided .yml file using the following commands:
```
conda env create -f snakemake.yml
conda activate snakemake
```

#### Analysis

To assess the effects of violations of demographic assumptions made by the single pulse admixture model, such as population growth, subsequent gene flow, and sex-specific assortative mating, we simulated these scenarios and evaluate the effects on inferred sex ratios.
Three continental ancestral population are simulated according to Gravel et al. (2011),<sup>5</sup> and admixture is simulated 15 generations ago. Then, random individuals are sampled from the ancestral and admixed population (we recommend 10,000), and LD pruning is performed using plink2.<sup>2</sup> On the pruned data, ADMIXTURE is run to separately infer X chromosomal and autosomal ancestry proportions using the individuals from the ancestral populations as training data.<sup>6</sup> Mitochondrial and Y chromosome haplogroup frequencies are inferred by clustering the ancestral haplogroups and assigning haplogroups of admixed individuals based proximity to ancestral haplogroup clusters. Finally, sex ratios are estimated using the single admixture models and bootstrapping individual ancestry proportions. Additionally, PCA is performed on the simulated populations to determine that the correct population structure was achievied, and site-frequency spectra (SFS) of the admixed population is plotted.
Simulations are done with SLiM3 using tree sequence recording.<sup>7,8</sup> Pyslim and msprime are then used for recapitation and superimposing mutations in Python.<sup>9,10</sup> 

#### Output files
1. Tree sequences will be stored in the directory specified with `tree_dir` in `config/config.yaml`
2. VCF (unpruned and pruned), plink files, haplogroup assignments, ADMIXTURE results will be stored in the directory specified with `data_dir` in `config/config.yaml`
3. Figures of the PCA and SFS will stored in the directory specified with `plot_dir` in `config/config.yaml`
4. Summary files of bootstrap results will stored in directory specified with `results_dir` in `config/config.yaml`

#### Example usage
The entire process is implemented in a snakemake workflow, which can be executed using the command below.<sup>11</sup>

```
snakemake --use-conda --conda-frontend --cores 16 
```

All simulations parameters, such as admixture proportions, sex biases, migration rates, sample sizes etc., can be set via the `config/config.yaml` file. To re-run all simulations that were presented in our paper (Table 2), run `./run_all_simulations.sh`. This scripts also generates Figure 3 in our manuscript, which assess the effect of sample size on the inferred sex ratios. Inferred sex ratios are compared for different sample sizes using a custom script:

```
scripts/assess_effect_of_sampling_size.py -qa <list_of_q_files_for_autosomes_for_different_sample_sizes> -qx <list_of_q_files_for_autosomes_for_different_sample_sizes> -f <fam_files_for_different_sample_sizes> -p ./ -o figure3.pdf
```


### References
1. Alexander, David H, John Novembre, and Kenneth Lange. 2009. “Fast Model-Based Estimation of Ancestry in Unrelated Individuals.” Genome Research 19 (9): 1655–64. https://doi.org/10.1101/gr.094052.109.
2. Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ. Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience. 2015;4:7. Published 2015 Feb 25. doi:10.1186/s13742-015-0047-8
3. Goldberg, A., and Rosenberg, N.A. (2015). Beyond 2/3 and 1/3: The Complex Signatures of Sex-Biased Admixture on the X Chromosome. Genetics 201, 263 LP – 279.
4. Micheletti, S.J., Bryc, K., Ancona Esselmann, S.G., Freyman, W.A., Moreno, M.E., Poznik, G.D., Shastri, A.J., Agee, M., Aslibekyan, S., Auton, A., et al. (2020). Genetic Consequences of the Transatlantic Slave Trade in the Americas. Am. J. Hum. Genet. 107, 265–277.<br>
5. Gravel, S., Henn, B.M., Gutenkunst, R.N., Indap, A.R., Marth, G.T., Clark, A.G., Yu, F., Gibbs, R.A., and Bustamante, C.D. (2011). Demographic history and rare allele sharing among human populations. Proc. Natl. Acad. Sci. U. S. A. 108, 11983–11988.
6. Alexander, D.H., Novembre, J., and Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. Genome Res. 19, 1655–1664.
7. Haller, B.C., and Messer, P.W. (2019). SLiM 3: Forward Genetic Simulations Beyond the Wright–Fisher Model. Mol. Biol. Evol. 36, 632–637. 
8. Haller, B.C., Galloway, J., Kelleher, J., Messer, P.W., and Ralph, P.L. (2019). Tree-sequence recording in SLiM opens new horizons for forward-time simulation of whole genomes. Mol. Ecol. Resour. 19, 552–566
9. Kelleher, J., Etheridge, A.M., and McVean, G. (2016). Efficient Coalescent Simulation and Genealogical Analysis for Large Sample Sizes. PLOS Comput. Biol. 12, 1–22
10. Nelson, D., Kelleher, J., Ragsdale, A.P., Moreau, C., McVean, G., and Gravel, S. (2020). Accounting for long-range correlations in genome-wide simulations of large cohorts. PLOS Genet. 16, 1–12.
11. Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33.
