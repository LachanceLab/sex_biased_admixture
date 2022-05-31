#!/usr/bin bash

# effect of sampling size
snakemake --config sample_size=100 results_dir=results_sex_bias_no_growth_sample_100 demographic_models=["single_pulse_admixture"] --use-conda --conda-frontend conda --latency-wait 20 --cores 16 --rerun-incomplete --keep-incomplete
snakemake --config sample_size=500 results_dir=results_sex_bias_no_growth_sample_500 demographic_models=["single_pulse_admixture"] --use-conda --conda-frontend conda --latency-wait 20 --cores 16 --rerun-incomplete --keep-incomplete
snakemake --config sample_size=1000 results_dir=results_sex_bias_no_growth_sample_1000 demographic_models=["single_pulse_admixture"] --use-conda --conda-frontend conda --latency-wait 20 --cores 16 --rerun-incomplete --keep-incomplete
snakemake --config sample_size=2500 results_dir=results_sex_bias_no_growth_sample_2500 demographic_models=["single_pulse_admixture"] --use-conda --conda-frontend conda --latency-wait 20 --cores 16 --rerun-incomplete --keep-incomplete
snakemake --config sample_size=5000 results_dir=results_sex_bias_no_growth_sample_5000 demographic_models=["single_pulse_admixture"] --use-conda --conda-frontend conda --latency-wait 20 --cores 16 --rerun-incomplete --keep-incomplete
snakemake --config sample_size=10000 results_dir=results_sex_bias_no_growth_sample_10000 demographic_models=["single_pulse_admixture"] --use-conda --conda-frontend conda --latency-wait 20 --cores 16 --rerun-incomplete --keep-incomplete
snakemake --config sample_size=20000 results_dir=results_sex_bias_no_growth_sample_20000 demographic_models=["single_pulse_admixture"] --use-conda --conda-frontend conda --latency-wait 20 --cores 16 --rerun-incomplete --keep-incomplete
snakemake --config sample_size=30000 results_dir=results_sex_bias_no_growth_sample_30000 demographic_models=["single_pulse_admixture"] --use-conda --conda-frontend conda --latency-wait 20 --cores 16 --rerun-incomplete --keep-incomplete
mv trees/single_pulse_admixture* results_sex_bias_no_growth_sample_10000
conda activate sex_ratios
scripts/assess_effect_of_sampling_size.py -qa results_sex_bias_no_growth_sample_100/single_pulse_admixture_g15_chr1.3.Q results_sex_bias_no_growth_sample_500/single_pulse_admixture_g15_chr1.3.Q results_sex_bias_no_growth_sample_1000/single_pulse_admixture_g15_chr1.3.Q results_sex_bias_no_growth_sample_2500/single_pulse_admixture_g15_chr1.3.Q results_sex_bias_no_growth_sample_5000/single_pulse_admixture_g15_chr1.3.Q results_sex_bias_no_growth_sample_10000/single_pulse_admixture_g15_chr1.3.Q results_sex_bias_no_growth_sample_20000/single_pulse_admixture_g15_chr1.3.Q results_sex_bias_no_growth_sample_30000/single_pulse_admixture_g15_chr1.3.Q -qx results_sex_bias_no_growth_sample_100/single_pulse_admixture_g15_chrX.3.Q results_sex_bias_no_growth_sample_500/single_pulse_admixture_g15_chrX.3.Q results_sex_bias_no_growth_sample_1000/single_pulse_admixture_g15_chrX.3.Q results_sex_bias_no_growth_sample_2500/single_pulse_admixture_g15_chrX.3.Q results_sex_bias_no_growth_sample_5000/single_pulse_admixture_g15_chrX.3.Q results_sex_bias_no_growth_sample_10000/single_pulse_admixture_g15_chrX.3.Q results_sex_bias_no_growth_sample_20000/single_pulse_admixture_g15_chrX.3.Q results_sex_bias_no_growth_sample_30000/single_pulse_admixture_g15_chrX.3.Q -f results_sex_bias_no_growth_sample_100/single_pulse_admixture_g15_chr1.fam results_sex_bias_no_growth_sample_500/single_pulse_admixture_g15_chr1.fam results_sex_bias_no_growth_sample_1000/single_pulse_admixture_g15_chr1.fam results_sex_bias_no_growth_sample_2500/single_pulse_admixture_g15_chr1.fam results_sex_bias_no_growth_sample_5000/single_pulse_admixture_g15_chr1.fam results_sex_bias_no_growth_sample_10000/single_pulse_admixture_g15_chr1.fam results_sex_bias_no_growth_sample_20000/single_pulse_admixture_g15_chr1.fam results_sex_bias_no_growth_sample_30000/single_pulse_admixture_g15_chr1.fam
conda deactivate


# constant admxiture
snakemake --config sample_size=10000 results_dir=results_sex_bias_no_growth_constant_1_sample_10000 demographic_models=["constant_admixture"] afr_const=0.05 eur_constant=0.025 ea_const=0.01 --use-conda --conda-frontend conda --latency-wait 20 --cores 16 --rerun-incomplete --keep-incomplete
mv trees/constant_admixture* results_sex_bias_no_growth_constant_1_sample_10000/

# population growth
snakemake --config sample_size=10000 results_dir=results_sex_bias_pop_growth_5_constant_5_sample_10000 demographic_models=["single_pulse_admixture"] amr_growth_rate=0.05 --use-conda --conda-frontend conda --latency-wait 20 --cores 16 --rerun-incomplete --keep-incomplete
mv trees/single_pulse_admixture* results_sex_bias_pop_growth_5_constant_5_sample_10000
snakemake --config sample_size=10000 results_dir=results_sex_bias_pop_growth_5_constant_5_sample_10000 demographic_models=["constant_admixture"] afr_const=0.05 eur_constant=0.025 ea_const=0.01 amr_growth_rate=0.05 --use-conda --conda-frontend conda --latency-wait 20 --cores 16 --rerun-incomplete --keep-incomplete
mv trees/constant_admixture* results_sex_bias_pop_growth_5_constant_5_sample_10000

# non-random mating
snakemake --config sample_size=10000 results_dir=results_sex_bias_non_random_mating_90_sample_10000 rejection_rate_female_eur_male_non_eur=0.9 demographic_models=["single_pulse_non_random_mating"] --use-conda --conda-frontend conda --latency-wait 20 --cores 16 --rerun-incomplete --keep-incomplete
mv trees/single_pulse_non_random* results_sex_bias_non_random_mating_sample_10000/

# non-random mating and constant admixture
snakemake --config afr_const=0.05 eur_constant=0.025 ea_const=0.01 sample_size=10000 results_dir=results_sex_bias_non_random_mating_sample_10000 rejection_rate_female_eur_male_non_eur=0.4 demographic_models=["constant_admixture_non_random_mating"] --use-conda --conda-frontend conda --latency-wait 20 --cores 16 --rerun-incomplete --keep-incomplete
mv trees/constant_admixture_non_random* results_sex_bias_non_random_mating_sample_10000/

# non-random mating, constant admixture and population growth
snakemake --config afr_const=0.05 eur_constant=0.025 ea_const=0.01 sample_size=10000 results_dir=results_sex_bias_growth_non_random_mating_sample_10000 amr_growth_rate=0.05 rejection_rate_female_eur_male_non_eur=0.4 demographic_models=["constant_admixture_non_random_mating"] --use-conda --conda-frontend conda --latency-wait 20 --cores 16 --rerun-incomplete --keep-incomplete
mv trees/constant_admixture_non_random* results_sex_bias_growth_non_random_mating_sample_10000/