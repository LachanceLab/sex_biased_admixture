# Estimating sex-biased admixture when genetic ancestry proportions do not add up 

### Abstract
<p>Sex-biased admixture can be inferred from ancestry-specific proportions of X chromosome and autosomes. This approach was used in a recent high-profile paper: Genetic Consequences of the Transatlantic Slave Trade in the Americas. In this paper, Micheletti et al. controversially claim that the contributions of African women and European men to American gene pools were much more biased than previously thought. Here we show that extreme sex-specific contributions can be attributed to unassigned genetic ancestries and the sensitivity of the applied mathematical models towards small differences in ancestry proportions. After accounting for unassigned ancestry proportions, we reanalyzed summary statistic data from Micheletti et al. Our reanalysis yielded modest sex ratios that are more consistent with observed haplogroup imbalances between mitochondrial DNA and Y chromosomes. These results underscore the importance of using adequate models, as well as accurate and complete ancestry information to infer sex-biased patterns of demography.</p>

<p>We reanalyze sex ratios for different ancestral populations for contemporary gene pools in the Americas using mean observed ancestry proportions inferred from the summary statistics by Micheletti et al. (2020)<sup>1</sup>. The observed proportions of African, European, and Native American ancestry are given in the <code>admixture.txt</code> file. We estimate the sex ratios for generation 2-15 after admixture using a dynamic model proposed by Goldberg and Rosenberg as well as an equilibrium model. The model by Goldberg and Rosenberg quickly converges to the equilibrium model.<sup>2</sup> For more details, see our publication.</p>

<p>The input file should be comma-separated, the first row should be the header, and have the following structure:
  
  ,H<sub>1</sub><sup>X</sup>,H<sub>1</sub><sup>A</sup>,H<sub>2</sub><sup>X</sup>,H<sub>2</sub><sup>A</sup>,...,p<sub>f</sub><br>
region1,x,a,xx,aa,...,f<br>
region2,x,a,xx,aa,...,f<br>

where H<sub>i</sub><sup>X</sup> and H<sub>i</sub><sup>A</sup> are the X chromosomal and autosomal ancestry proportions that trace to source population 1, and p<sub>f</sub> are the fraction of females in the sample in a given region. Comments are indicated by a "#". </p>

<p> To set up the conda environment run the following commands:<br>
  <code>conda env create -f environment.yml</code><br>
  <code>conda activate sex_ratios</code>

<p>Running the script:

<code> python estimate_sex_ratios_single_admixture_pulse.py -i admixture.txt -o <output.xlsx></code>

yields three output files:
1. output.xlsx: Sex ratios for generations 2-15 and the equilibrium value. This is a corrected version of Table S9B in Micheletti et al. (2020)<sup>1</sup>. This table corresponds to Table S3 in our publication.</p>
2. output_equally.xlsx: Sex ratios for generations 2-15 after admixture and the equilibrium after distributing unassigned ancestries equally to all ancestral populations. This table corresponds to Table S4 in our publication.
3. figure1.pdf: **A)** Illustration of the sensitivity of the applied models towards small changes in the observed ancestry proportions. **B)** Exploration of the parameter space leading to male-biased admixture, female-biased admixture, or model breakage. ![](figure1.pdf)

By default, we distribute the unassigned ancestry equally to all ancestral bins, but we also support distributing them proportionally by setting the <code>-p</code> flag. In this case, output file two is named output_proportionally.xlsx (Table S5 in our publication).

</p>

<p>Although not used in our publication, we also provide a Python implementation of the model for a constant admixture scenario proposed by Goldberg and Rosenberg.<sup>2</sup> The usage of the script is equivalent to the usage of the script for a single admixture event. The additional parameter <code>-n <float></code> allows adjusting the increment for the grid-search to find the optimal set of values for s<sub>1</sub><sup>f</sup>, s<sub>1</sub><sup>m</sup>, s<sub>2</sub><sup>f</sup>, and s<sub>2</sub><sup>m</sup>. Note that this script explicitly assumes that there are only two source populations.

<code>python estimate_sex_ratios_constant_admixture_pulse.py -i admixture.txt -o <output.xlsx> -n <float></code>
</p>
  
### References
<sup>1</sup> Micheletti, S.J., Bryc, K., Ancona Esselmann, S.G., Freyman, W.A., Moreno, M.E., Poznik, G.D., Shastri, A.J., Agee, M., Aslibekyan, S., Auton, A., et al. (2020). Genetic Consequences of the Transatlantic Slave Trade in the Americas. Am. J. Hum. Genet. 107, 265–277.<br>
<sup>2</sup> Goldberg, A., and Rosenberg, N.A. (2015). Beyond 2/3 and 1/3: The Complex Signatures of Sex-Biased Admixture on the X Chromosome. Genetics 201, 263 LP – 279. 
