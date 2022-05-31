import numpy as np

configfile: "config/config.yaml"

mutation_rate=config["mutation_rate"]
recombination_rate=config["recombination_rate"]
mtDNA_mutation_rate=config['mtDNA_mutation_rate']
chromosome_length=config['chrom_length']
mtDNA_length=config['mtDNA_length']
tree_dir=config["tree_dir"]
#data_dir=config['data_dir']
data_dir=config["results_dir"]
populations=['afr', 'eur', 'ea', 'amr']
afr=config['afr']
eur=config['eur']
ea=config['ea']
sm_afr=config['sm_afr']
sm_eur=config['sm_eur']
sm_ea=config['sm_ea']
N_amr=config['adx_n']
chromosomes=['chr1', 'chrX', 'chrY', 'mtDNA']
plot_dir=config['plot_dir']
results_dir=config['results_dir']
first_gen=config['first_gen']
last_gen=config['last_gen']



rule all:
    input:
        expand("{results_dir}/{demographic_model}_g{gen}_ancestry_proportions.txt", results_dir=results_dir,
               gen=np.arange(first_gen, last_gen +1), demographic_model=config["demographic_models"]),
        expand("{results_dir}/{demographic_model}_g{gen}_haplogroup_frequencies.txt", results_dir=results_dir,
               gen=np.arange(first_gen, last_gen +1), demographic_model=config["demographic_models"]),
       	expand("{results_dir}/{demographic_model}_g{gen}_bootstrapped.txt", results_dir=results_dir,
               gen=np.arange(first_gen, last_gen +1), demographic_model=config["demographic_models"]),
        expand("{plot_dir}/{demographic_model}_g{gen}_{chrom}_pca.png",
               chrom=['chr1', 'chrX'], demographic_model=config['demographic_models'], plot_dir=plot_dir,
               gen=np.arange(first_gen, last_gen +1))

rule simulate_afr_eur_ea:
    output:
        f"{tree_dir}/afr_eur_ea.trees"
    threads: 1
    conda: "envs/slim.yml"
    params:
        recipe=config["afr_eur_ea_recipe"]
    shell:
        "scripts/afr_eur_ea_populations.py -r {recombination_rate} --recipe {params.recipe} " \
        " -l {chromosome_length} --mtDNA_length {mtDNA_length} -o {output}"

use rule simulate_afr_eur_ea as simulate_afr_eur_ea_non_random_mating with:
    output:
        f"{tree_dir}/afr_eur_ea_non_random_mating.trees"
    params:
        recipe=config["afr_eur_ea_non_random_mating_recipe"]

rule single_pulse_admixture:
    input:
        rules.simulate_afr_eur_ea.output
    output:
        expand("{tree_dir}/single_pulse_admixture_g{gen}.trees", tree_dir=tree_dir, gen=np.arange(first_gen, last_gen + 1))
    threads: 1
    conda: "envs/slim.yml"
    shell:
        "scripts/american_admixture.py -r {recombination_rate} --recipe {config[single_pulse_admixture_recipe]} " \
        " -l {chromosome_length} --mtDNA_length {mtDNA_length} --afr {afr} --eur {eur} --ea {ea} --sm_afr {sm_afr} "\
        "--sm_eur {sm_eur} --sm_ea {sm_ea} --N_amr {N_amr} --tree_sequences {input} "
        "--amr_growth_rate {config[amr_growth_rate]}"

rule single_pulse_admixture_non_random_mating:
    input:
        rules.simulate_afr_eur_ea_non_random_mating.output
    output:
        expand("{tree_dir}/single_pulse_non_random_mating_g{gen}.trees", tree_dir={tree_dir}, gen=np.arange(first_gen, last_gen + 1))
    threads: 1
    conda: "envs/slim.yml"
    shell:
        "scripts/american_admixture.py -r {recombination_rate} --recipe {config[single_pulse_admixture_non_random_mating_recipe]} " \
        " -l {chromosome_length} --mtDNA_length {mtDNA_length} --afr {afr} --eur {eur} --ea {ea} --sm_afr {sm_afr} "\
        "--sm_eur {sm_eur} --sm_ea {sm_ea} --N_amr {N_amr} --tree_sequences {input} "
        "--amr_growth_rate {config[amr_growth_rate]} "
        "--rejection_rate_female_eur_male_non_eur {config[rejection_rate_female_eur_male_non_eur]}"

rule constant_admixture:
    input:
        rules.simulate_afr_eur_ea.output
    output:
        expand("{tree_dir}/constant_admixture_g{gen}.trees", tree_dir=tree_dir, gen=np.arange(first_gen, last_gen +1))
    threads: 1
    conda: "envs/slim.yml"
    shell:
        "scripts/american_admixture.py -r {recombination_rate} --recipe {config[constant_admixture_recipe]} " \
        " -l {chromosome_length} --mtDNA_length {mtDNA_length} --afr {afr} --eur {eur} --ea {ea} --sm_afr {sm_afr} "\
        "--sm_eur {sm_eur} --sm_ea {sm_ea} --N_amr {N_amr} --tree_sequences {input} "
        "--afr_const {config[afr_const]} --eur_const {config[eur_const]} "\
        "--ea_const {config[ea_const]} --sm_afr_const {config[sm_afr_const]} --sm_eur_const {config[sm_eur_const]} "\
        "--sm_ea_const {config[sm_ea_const]} --amr_growth_rate {config[amr_growth_rate]}"

rule constant_admixture_non_random_mating:
    input:
        rules.simulate_afr_eur_ea_non_random_mating.output
    output:
        expand("{tree_dir}/constant_admixture_non_random_mating_g{gen}.trees", tree_dir={tree_dir}, gen=np.arange(first_gen, last_gen + 1))
    threads: 1
    conda: "envs/slim.yml"
    shell:
        "scripts/american_admixture.py -r {recombination_rate} --recipe {config[constant_admixture_non_random_mating_recipe]} " \
        " -l {chromosome_length} --mtDNA_length {mtDNA_length} --afr {afr} --eur {eur} --ea {ea} --sm_afr {sm_afr} "\
        "--sm_eur {sm_eur} --sm_ea {sm_ea} --N_amr {N_amr} --tree_sequences {input} "\
        "--amr_growth_rate {config[amr_growth_rate]} --afr_const {config[afr_const]} --eur_const {config[eur_const]} "\
        "--ea_const {config[ea_const]} --sm_afr_const {config[sm_afr_const]} --sm_eur_const {config[sm_eur_const]} "\
        "--sm_ea_const {config[sm_ea_const]} "
        "--rejection_rate_female_eur_male_non_eur {config[rejection_rate_female_eur_male_non_eur]}"

rule generate_vcf_single_pulse:
    input:
        rules.single_pulse_admixture.output
    output:
        ts=expand("{tree_dir}/single_pulse_admixture_g{gen}_mutated.trees", tree_dir=tree_dir,
                  gen=np.arange(first_gen, last_gen +1)),
        pkl=expand("{tree_dir}/single_pulse_admixture_g{gen}_mutated.pkl",  tree_dir=tree_dir,
                   gen=np.arange(first_gen, last_gen +1)),
        vcf_chr1=expand("{data_dir}/single_pulse_admixture_g{gen}_chr1.vcf", gen=np.arange(first_gen, last_gen +1),
                        data_dir=data_dir),
        vcf_chrX=expand("{data_dir}/single_pulse_admixture_g{gen}_chrX.vcf",gen=np.arange(first_gen, last_gen +1),
                        data_dir=data_dir),
        vcf_chrY=expand("{data_dir}/single_pulse_admixture_g{gen}_chrY.vcf", gen=np.arange(first_gen, last_gen +1),
                        data_dir=data_dir),
        vcf_mtDNA=expand("{data_dir}/single_pulse_admixture_g{gen}_mtDNA.vcf", gen=np.arange(first_gen, last_gen +1),
                         data_dir=data_dir),
        sex=expand("{data_dir}/single_pulse_admixture_g{gen}_sex.tab", gen=np.arange(first_gen, last_gen +1),
                   data_dir=data_dir),
        pop_chr1=expand("{data_dir}/single_pulse_admixture_g{gen}_chr1.pop",gen=np.arange(first_gen,last_gen + 1),
                         data_dir=data_dir),
        pop_chrX=expand("{data_dir}/single_pulse_admixture_g{gen}_chrX.pop",gen=np.arange(first_gen,last_gen + 1),
                         data_dir=data_dir),
        pop_chrY=expand("{data_dir}/single_pulse_admixture_g{gen}_chrY.pop",gen=np.arange(first_gen,last_gen + 1),
                         data_dir=data_dir),
        pop_mtDNA=expand("{data_dir}/single_pulse_admixture_g{gen}_mtDNA.pop",gen=np.arange(first_gen,last_gen + 1),
                          data_dir=data_dir),

    threads: 4
    conda: "envs/slim.yml"
    log: "log/adding_mutations_single_pulse.log"
    shell:
        "scripts/add_mutations_to_tree_sequence.py -N {config[ancestral_Ne]} -m {mutation_rate} "\
        "--mtDNA_mutation_rate {mtDNA_mutation_rate} -o {output.ts} --output_vcf_dir {data_dir} "
        "-l {chromosome_length} --mtDNA_length {mtDNA_length} --populations {populations} --log {log} "\
        "--tree_sequences {input} --threads {threads} --sample_size {config[sample_size]} --plot_sfs "
        "--recapitate"

use rule generate_vcf_single_pulse as generate_vcf_constant with:
    input:
        rules.constant_admixture.output
    output:
        ts=expand("{tree_dir}/constant_admixture_g{gen}_mutated.trees", tree_dir=tree_dir,
                  gen=np.arange(first_gen, last_gen +1)),
        pkl=expand("{tree_dir}/constant_admixture_g{gen}_mutated.pkl",  tree_dir=tree_dir,
                   gen=np.arange(first_gen, last_gen +1)),
        vcf_chr1=expand("{data_dir}/constant_admixture_g{gen}_chr1.vcf", gen=np.arange(first_gen, last_gen +1),
                        data_dir=data_dir),
        vcf_chrX=expand("{data_dir}/constant_admixture_g{gen}_chrX.vcf",gen=np.arange(first_gen, last_gen +1),
                        data_dir=data_dir),
        vcf_chrY=expand("{data_dir}/constant_admixture_g{gen}_chrY.vcf", gen=np.arange(first_gen, last_gen +1),
                        data_dir=data_dir),
        vcf_mtDNA=expand("{data_dir}/constant_admixture_g{gen}_mtDNA.vcf", gen=np.arange(first_gen, last_gen +1),
                         data_dir=data_dir),
        sex=expand("{data_dir}/constant_admixture_g{gen}_sex.tab", gen=np.arange(first_gen, last_gen +1),
                   data_dir=data_dir),
        pop_chr1=expand("{data_dir}/constant_admixture_g{gen}_chr1.pop",gen=np.arange(first_gen,last_gen + 1),
                         data_dir=data_dir),
        pop_chrX=expand("{data_dir}/constant_admixture_g{gen}_chrX.pop",gen=np.arange(first_gen,last_gen + 1),
                         data_dir=data_dir),
        pop_chrY=expand("{data_dir}/constant_admixture_g{gen}_chrY.pop",gen=np.arange(first_gen,last_gen + 1),
                         data_dir=data_dir),
        pop_mtDNA=expand("{data_dir}/constant_admixture_g{gen}_mtDNA.pop",gen=np.arange(first_gen,last_gen + 1),
                         data_dir=data_dir),
    log: "log/adding_mutations_constant.log"


use rule generate_vcf_single_pulse as generate_vcf_non_random_mating with:
    input:
        rules.single_pulse_admixture_non_random_mating.output
    output:
        ts=expand("{tree_dir}/single_pulse_non_random_mating_g{gen}_mutated.trees", tree_dir=tree_dir,
                  gen=np.arange(first_gen, last_gen +1)),
        pkl=expand("{tree_dir}/single_pulse_non_random_mating_g{gen}_mutated.pkl",  tree_dir=tree_dir,
                   gen=np.arange(first_gen, last_gen +1)),
        vcf_chr1=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chr1.vcf", gen=np.arange(first_gen, last_gen +1),
                        data_dir=data_dir),
        vcf_chrX=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrX.vcf",gen=np.arange(first_gen, last_gen +1),
                        data_dir=data_dir),
        vcf_chrY=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrY.vcf", gen=np.arange(first_gen, last_gen +1),
                        data_dir=data_dir),
        vcf_mtDNA=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_mtDNA.vcf", gen=np.arange(first_gen, last_gen +1),
                         data_dir=data_dir),
        sex=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_sex.tab", gen=np.arange(first_gen, last_gen +1),
                   data_dir=data_dir),
        pop_chr1=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chr1.pop",gen=np.arange(first_gen,last_gen + 1),
                         data_dir=data_dir),
        pop_chrX=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrX.pop",gen=np.arange(first_gen,last_gen + 1),
                         data_dir=data_dir),
        pop_chrY=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrY.pop",gen=np.arange(first_gen,last_gen + 1),
                         data_dir=data_dir),
        pop_mtDNA=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_mtDNA.pop",gen=np.arange(first_gen,last_gen + 1),
                          data_dir=data_dir),
    log: "log/adding_mutations_non_random_mating.log"


use rule generate_vcf_single_pulse as generate_vcf_non_random_mating_constant with:
    input:
        rules.constant_admixture_non_random_mating.output
    output:
        ts=expand("{tree_dir}/constant_admixture_non_random_mating_g{gen}_mutated.trees", tree_dir=tree_dir,
                  gen=np.arange(first_gen, last_gen +1)),
        pkl=expand("{tree_dir}/constant_admixture_non_random_mating_g{gen}_mutated.pkl",  tree_dir=tree_dir,
                   gen=np.arange(first_gen, last_gen +1)),
        vcf_chr1=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chr1.vcf", gen=np.arange(first_gen, last_gen +1),
                        data_dir=data_dir),
        vcf_chrX=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrX.vcf",gen=np.arange(first_gen, last_gen +1),
                        data_dir=data_dir),
        vcf_chrY=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrY.vcf", gen=np.arange(first_gen, last_gen +1),
                        data_dir=data_dir),
        vcf_mtDNA=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_mtDNA.vcf", gen=np.arange(first_gen, last_gen +1),
                         data_dir=data_dir),
        sex=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_sex.tab", gen=np.arange(first_gen, last_gen +1),
                   data_dir=data_dir),
        pop_chr1=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chr1.pop",gen=np.arange(first_gen,last_gen + 1),
                         data_dir=data_dir),
        pop_chrX=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrX.pop",gen=np.arange(first_gen,last_gen + 1),
                         data_dir=data_dir),
        pop_chrY=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrY.pop",gen=np.arange(first_gen,last_gen + 1),
                         data_dir=data_dir),
        pop_mtDNA=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_mtDNA.pop",gen=np.arange(first_gen,last_gen + 1),
                          data_dir=data_dir),
    log: "log/adding_mutations_non_random_mating_constant.log"

rule compress_vcf_single_pulse:
    input:
        rules.generate_vcf_single_pulse.output.vcf_chr1,
        rules.generate_vcf_single_pulse.output.vcf_chrX,
        rules.generate_vcf_single_pulse.output.vcf_chrY,
        rules.generate_vcf_single_pulse.output.vcf_mtDNA
    output:
        vcf_chr1=expand("{data_dir}/single_pulse_admixture_g{gen}_chr1.vcf.gz", data_dir=data_dir,
                        gen=np.arange(first_gen, last_gen +1)),
        vcf_chrX=expand("{data_dir}/single_pulse_admixture_g{gen}_chrX.vcf.gz", data_dir=data_dir,
                        gen=np.arange(first_gen, last_gen +1)),
        vcf_chrY=expand("{data_dir}/single_pulse_admixture_g{gen}_chrY.vcf.gz", data_dir=data_dir,
                        gen=np.arange(first_gen, last_gen +1)),
        vcf_mtDNA=expand("{data_dir}/single_pulse_admixture_g{gen}_mtDNA.vcf.gz", data_dir=data_dir,
                         gen=np.arange(first_gen, last_gen +1)),
        tbi=expand("{data_dir}/single_pulse_admixture_g{gen}_{chrom}.vcf.gz.tbi", chrom=chromosomes,
                   data_dir=data_dir, gen=np.arange(first_gen, last_gen +1))
    threads: 8
    conda: "envs/tabix.yaml"
    shell:
        "for vcf_file in {input}; "
        "do "
        "bcftools norm -m - --threads {threads} ${{vcf_file}} | bcftools annotate --threads {threads} "
        "--set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o $( echo ${{vcf_file}} | sed 's/_filtered//' ).gz; "
        "tabix ${{vcf_file}}.gz; "
        "done"

use rule compress_vcf_single_pulse as compress_vcf_constant with:
    input:
        rules.generate_vcf_constant.output.vcf_chr1,
        rules.generate_vcf_constant.output.vcf_chrX,
        rules.generate_vcf_constant.output.vcf_chrY,
        rules.generate_vcf_constant.output.vcf_mtDNA
    output:
        vcf_chr1=expand("{data_dir}/constant_admixture_g{gen}_chr1.vcf.gz",data_dir=data_dir,
                        gen=np.arange(first_gen,last_gen + 1)),
        vcf_chrX=expand("{data_dir}/constant_admixture_g{gen}_chrX.vcf.gz",data_dir=data_dir,
                        gen=np.arange(first_gen,last_gen + 1)),
        vcf_chrY=expand("{data_dir}/constant_admixture_g{gen}_chrY.vcf.gz",data_dir=data_dir,
                        gen=np.arange(first_gen,last_gen + 1)),
        vcf_mtDNA=expand("{data_dir}/constant_admixture_g{gen}_mtDNA.vcf.gz",data_dir=data_dir,
                         gen=np.arange(first_gen,last_gen + 1)),
        tbi=expand("{data_dir}/constant_admixture_g{gen}_{chrom}.vcf.gz.tbi",
                   chrom=chromosomes, data_dir=data_dir, gen=np.arange(first_gen, last_gen +1))

use rule compress_vcf_single_pulse as compress_vcf_non_random_mating with:
    input:
        rules.generate_vcf_non_random_mating.output.vcf_chr1,
        rules.generate_vcf_non_random_mating.output.vcf_chrX,
        rules.generate_vcf_non_random_mating.output.vcf_chrY,
        rules.generate_vcf_non_random_mating.output.vcf_mtDNA
    output:
        vcf_chr1=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chr1.vcf.gz", data_dir=data_dir,
                        gen=np.arange(first_gen, last_gen +1)),
        vcf_chrX=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrX.vcf.gz", data_dir=data_dir,
                        gen=np.arange(first_gen, last_gen +1)),
        vcf_chrY=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrY.vcf.gz", data_dir=data_dir,
                        gen=np.arange(first_gen, last_gen +1)),
        vcf_mtDNA=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_mtDNA.vcf.gz", data_dir=data_dir,
                         gen=np.arange(first_gen, last_gen +1)),
        tbi=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_{chrom}.vcf.gz.tbi",
                   chrom=chromosomes, data_dir=data_dir, gen=np.arange(first_gen, last_gen +1))

use rule compress_vcf_single_pulse as compress_vcf_non_random_mating_constant with:
    input:
        rules.generate_vcf_non_random_mating_constant.output.vcf_chr1,
        rules.generate_vcf_non_random_mating_constant.output.vcf_chrX,
        rules.generate_vcf_non_random_mating_constant.output.vcf_chrY,
        rules.generate_vcf_non_random_mating_constant.output.vcf_mtDNA
    output:
        vcf_chr1=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chr1.vcf.gz", data_dir=data_dir,
                        gen=np.arange(first_gen, last_gen +1)),
        vcf_chrX=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrX.vcf.gz", data_dir=data_dir,
                        gen=np.arange(first_gen, last_gen +1)),
        vcf_chrY=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrY.vcf.gz", data_dir=data_dir,
                        gen=np.arange(first_gen, last_gen +1)),
        vcf_mtDNA=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_mtDNA.vcf.gz", data_dir=data_dir,
                         gen=np.arange(first_gen, last_gen +1)),
        tbi=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_{chrom}.vcf.gz.tbi",
                   chrom=chromosomes, data_dir=data_dir, gen=np.arange(first_gen, last_gen +1))

rule ld_prune_chr1_single_pulse:
    input:
        vcf=rules.compress_vcf_single_pulse.output.vcf_chr1,
        sex=rules.generate_vcf_single_pulse.output.sex
    output:
        in_eq=expand("{data_dir}/single_pulse_admixture_g{gen}_chr1.prune.in", data_dir=data_dir,
                     gen=np.arange(first_gen, last_gen +1)),
        out=expand("{data_dir}/single_pulse_admixture_g{gen}_chr1.prune.out", data_dir=data_dir,
                   gen=np.arange(first_gen, last_gen +1))
    conda: "envs/plink.yaml"
    params:
        out=expand("{data_dir}/single_pulse_admixture_g{gen}_chr1", data_dir=data_dir, gen=np.arange(first_gen, last_gen +1))
    threads: 1
    shell:
       "scripts/ld_pruning_plink2.py --threads {threads} --vcf {input.vcf} --output_base {params.out} --sex {input.sex}"

use rule ld_prune_chr1_single_pulse as ld_prune_chr1_constant with:
    input:
        vcf=rules.compress_vcf_constant.output.vcf_chr1,
        sex=rules.generate_vcf_constant.output.sex
    output:
        in_eq=expand("{data_dir}/constant_admixture_g{gen}_chr1.prune.in", data_dir=data_dir,
                     gen=np.arange(first_gen, last_gen +1)),
        out=expand("{data_dir}/constant_admixture_g{gen}_chr1.prune.out", data_dir=data_dir,
                   gen=np.arange(first_gen, last_gen +1))
    params:
        out=expand("{data_dir}/constant_admixture_g{gen}_chr1", data_dir=data_dir,
                   gen=np.arange(first_gen, last_gen +1))

use rule ld_prune_chr1_single_pulse as ld_prune_chr1_non_random_mating with:
    input:
        vcf=rules.compress_vcf_non_random_mating.output.vcf_chr1,
        sex=rules.generate_vcf_non_random_mating.output.sex
    output:
        in_eq=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chr1.prune.in", data_dir=data_dir,
                     gen=np.arange(first_gen, last_gen +1)),
        out=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chr1.prune.out", data_dir=data_dir,
                   gen=np.arange(first_gen, last_gen +1))
    params:
        out=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chr1", data_dir=data_dir,
                   gen=np.arange(first_gen, last_gen +1))

use rule ld_prune_chr1_single_pulse as ld_prune_chr1_non_random_mating_constant with:
    input:
        vcf=rules.compress_vcf_non_random_mating_constant.output.vcf_chr1,
        sex=rules.generate_vcf_non_random_mating_constant.output.sex
    output:
        in_eq=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chr1.prune.in", data_dir=data_dir,
                     gen=np.arange(first_gen, last_gen +1)),
        out=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chr1.prune.out", data_dir=data_dir,
                   gen=np.arange(first_gen, last_gen +1))
    params:
        out=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chr1", data_dir=data_dir,
                   gen=np.arange(first_gen, last_gen +1))

use rule ld_prune_chr1_single_pulse as ld_prune_chrX_single_pulse with:
    input:
        vcf=rules.compress_vcf_single_pulse.output.vcf_chrX,
        sex=rules.generate_vcf_single_pulse.output.sex
    output:
        in_eq=expand("{data_dir}/single_pulse_admixture_g{gen}_chrX.prune.in", data_dir=data_dir,
                     gen=np.arange(first_gen, last_gen +1)),
        out=expand("{data_dir}/single_pulse_admixture_g{gen}_chrX.prune.out", data_dir=data_dir,
                   gen=np.arange(first_gen, last_gen +1))
    params:
        out=expand("{data_dir}/single_pulse_admixture_g{gen}_chrX", data_dir=data_dir, gen=np.arange(first_gen, last_gen +1))

use rule ld_prune_chr1_single_pulse as ld_prune_chrX_constant with:
    input:
        vcf=rules.compress_vcf_constant.output.vcf_chrX,
        sex=rules.generate_vcf_constant.output.sex
    output:
        in_eq=expand("{data_dir}/constant_admixture_g{gen}_chrX.prune.in",data_dir=data_dir,gen=np.arange(first_gen,
                                                                                                          last_gen + 1)),
        out=expand("{data_dir}/constant_admixture_g{gen}_chrX.prune.out",data_dir=data_dir,gen=np.arange(first_gen,
                                                                                                         last_gen + 1))
    params:
        out=expand("{data_dir}/constant_admixture_g{gen}_chrX",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))

use rule ld_prune_chr1_single_pulse as ld_prune_chrX_non_random_mating with:
    input:
        vcf=rules.compress_vcf_non_random_mating.output.vcf_chrX,
        sex=rules.generate_vcf_non_random_mating.output.sex
    output:
        in_eq=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrX.prune.in",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        out=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrX.prune.out",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    params:
        out=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrX",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))

use rule ld_prune_chr1_single_pulse as ld_prune_chrX_non_random_mating_constant with:
    input:
        vcf=rules.compress_vcf_non_random_mating_constant.output.vcf_chrX,
        sex=rules.generate_vcf_non_random_mating_constant.output.sex
    output:
        in_eq=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrX.prune.in",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        out=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrX.prune.out",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    params:
        out=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrX",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))

rule convert_vcf_to_plink_chr1_single_pulse:
    input:
        pruned=rules.ld_prune_chr1_single_pulse.output.in_eq,
        vcf=rules.compress_vcf_single_pulse.output.vcf_chr1,
        sex=rules.generate_vcf_single_pulse.output.sex,
    output:
        fam=expand("{data_dir}/single_pulse_admixture_g{gen}_chr1.fam", data_dir=data_dir,
                    gen=np.arange(first_gen, last_gen +1)),
        bed=expand("{data_dir}/single_pulse_admixture_g{gen}_chr1.bed", data_dir=data_dir,
                   gen=np.arange(first_gen, last_gen +1)),
        eigenval=expand("{data_dir}/single_pulse_admixture_g{gen}_chr1.eigenval", data_dir=data_dir,
                        gen=np.arange(first_gen, last_gen +1)),
        eigenvec=expand("{data_dir}/single_pulse_admixture_g{gen}_chr1.eigenvec", data_dir=data_dir,
                        gen=np.arange(first_gen, last_gen +1))
    params:
        out=expand("{data_dir}/single_pulse_admixture_g{gen}_chr1", data_dir=data_dir,
                   gen=np.arange(first_gen, last_gen +1))
    conda: "envs/plink.yaml"
    threads: 4
    shell:
        "scripts/convert_vcf_to_plink.py --threads {threads} --vcf {input.vcf} --sex {input.sex} "
        "--output_base {params.out} --pruned {input.pruned}"


use rule convert_vcf_to_plink_chr1_single_pulse as convert_vcf_to_plink_chrX_single_pulse with:
    input:
        vcf=rules.compress_vcf_single_pulse.output.vcf_chrX,
        sex=rules.generate_vcf_single_pulse.output.sex,
        pruned=rules.ld_prune_chrX_single_pulse.output.in_eq
    output:
        fam=expand("{data_dir}/single_pulse_admixture_g{gen}_chrX.fam",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        bed=expand("{data_dir}/single_pulse_admixture_g{gen}_chrX.bed",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenval=expand("{data_dir}/single_pulse_admixture_g{gen}_chrX.eigenval",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenvec=expand("{data_dir}/single_pulse_admixture_g{gen}_chrX.eigenvec",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    params:
        out=expand("{data_dir}/single_pulse_admixture_g{gen}_chrX",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))

use rule convert_vcf_to_plink_chr1_single_pulse as convert_vcf_to_plink_chr1_constant with:
    input:
        vcf=rules.compress_vcf_constant.output.vcf_chr1,
        sex=rules.generate_vcf_constant.output.sex,
        pruned=rules.ld_prune_chr1_constant.output.in_eq
    output:
        fam=expand("{data_dir}/constant_admixture_g{gen}_chr1.fam",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        bed=expand("{data_dir}/constant_admixture_g{gen}_chr1.bed",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenval=expand("{data_dir}/constant_admixture_g{gen}_chr1.eigenval",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenvec=expand("{data_dir}/constant_admixture_g{gen}_chr1.eigenvec",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    params:
        out=expand("{data_dir}/constant_admixture_g{gen}_chr1",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))

use rule convert_vcf_to_plink_chr1_single_pulse as convert_vcf_to_plink_chrX_constant with:
    input:
        vcf=rules.compress_vcf_constant.output.vcf_chrX,
        sex=rules.generate_vcf_constant.output.sex,
        pruned=rules.ld_prune_chrX_constant.output.in_eq
    output:
        fam=expand("{data_dir}/constant_admixture_g{gen}_chrX.fam",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        bed=expand("{data_dir}/constant_admixture_g{gen}_chrX.bed",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenval=expand("{data_dir}/constant_admixture_g{gen}_chrX.eigenval",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenvec=expand("{data_dir}/constant_admixture_g{gen}_chrX.eigenvec",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    params:
        out=expand("{data_dir}/constant_admixture_g{gen}_chrX",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))

use rule convert_vcf_to_plink_chr1_single_pulse as convert_vcf_to_plink_chr1_non_random_mating with:
    input:
        vcf=rules.compress_vcf_non_random_mating.output.vcf_chr1,
        sex=rules.generate_vcf_non_random_mating.output.sex,
        pruned=rules.ld_prune_chr1_non_random_mating.output.in_eq
    output:
        fam=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chr1.fam",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        bed=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chr1.bed",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenval=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chr1.eigenval",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenvec=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chr1.eigenvec",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    params:
        out=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chr1",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))

use rule convert_vcf_to_plink_chr1_single_pulse as convert_vcf_to_plink_chrX_non_random_mating with:
    input:
        vcf=rules.compress_vcf_non_random_mating.output.vcf_chrX,
        sex=rules.generate_vcf_non_random_mating.output.sex,
        pruned=rules.ld_prune_chrX_non_random_mating.output.in_eq
    output:
        fam=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrX.fam",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        bed=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrX.bed",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenval=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrX.eigenval",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenvec=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrX.eigenvec",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    params:
        out=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrX",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))

use rule convert_vcf_to_plink_chr1_single_pulse as convert_vcf_to_plink_chr1_non_random_mating_constant with:
    input:
        vcf=rules.compress_vcf_non_random_mating_constant.output.vcf_chr1,
        sex=rules.generate_vcf_non_random_mating_constant.output.sex,
        pruned=rules.ld_prune_chr1_non_random_mating_constant.output.in_eq
    output:
        fam=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chr1.fam",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        bed=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chr1.bed",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenval=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chr1.eigenval",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenvec=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chr1.eigenvec",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    params:
        out=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chr1",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))

use rule convert_vcf_to_plink_chr1_single_pulse as convert_vcf_to_plink_chrX_non_random_mating_constant with:
    input:
        vcf=rules.compress_vcf_non_random_mating_constant.output.vcf_chrX,
        sex=rules.generate_vcf_non_random_mating_constant.output.sex,
        pruned=rules.ld_prune_chrX_non_random_mating_constant.output.in_eq
    output:
        fam=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrX.fam",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        bed=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrX.bed",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenval=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrX.eigenval",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenvec=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrX.eigenvec",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    params:
        out=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrX",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))

rule convert_vcf_to_plink_chrY_single_pulse:
    input:
        vcf=rules.compress_vcf_single_pulse.output.vcf_chrY,
        sex=rules.generate_vcf_single_pulse.output.sex
    output:
        fam=expand("{data_dir}/single_pulse_admixture_g{gen}_chrY.fam",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        bed=expand("{data_dir}/single_pulse_admixture_g{gen}_chrY.bed",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenval=expand("{data_dir}/single_pulse_admixture_g{gen}_chrY.eigenval",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenvec=expand("{data_dir}/single_pulse_admixture_g{gen}_chrY.eigenvec",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    params:
        out=expand("{data_dir}/single_pulse_admixture_g{gen}_chrY",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    conda: "envs/plink.yaml"
    threads: 4
    shell:
        "scripts/convert_vcf_to_plink.py --threads {threads} --vcf {input.vcf} --sex {input.sex} "
        "--output_base {params.out}"

use rule convert_vcf_to_plink_chrY_single_pulse  as convert_vcf_to_plink_mtDNA_single_pulse with:
    input:
        vcf=rules.compress_vcf_single_pulse.output.vcf_mtDNA,
        sex=rules.generate_vcf_single_pulse.output.sex
    output:
        fam=expand("{data_dir}/single_pulse_admixture_g{gen}_mtDNA.fam",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        bed=expand("{data_dir}/single_pulse_admixture_g{gen}_mtDNA.bed",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenval=expand("{data_dir}/single_pulse_admixture_g{gen}_mtDNA.eigenval",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenvec=expand("{data_dir}/single_pulse_admixture_g{gen}_mtDNA.eigenvec",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    params:
        out=expand("{data_dir}/single_pulse_admixture_g{gen}_mtDNA",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))


use rule convert_vcf_to_plink_chrY_single_pulse as convert_vcf_to_plink_chrY_constant with:
    input:
        vcf=rules.compress_vcf_constant.output.vcf_chrY,
        sex=rules.generate_vcf_constant.output.sex
    output:
        fam=expand("{data_dir}/constant_admixture_g{gen}_chrY.fam",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        bed=expand("{data_dir}/constant_admixture_g{gen}_chrY.bed",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenval=expand("{data_dir}/constant_admixture_g{gen}_chrY.eigenval",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenvec=expand("{data_dir}/constant_admixture_g{gen}_chrY.eigenvec",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    params:
        out=expand("{data_dir}/constant_admixture_g{gen}_chrY",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))

use rule convert_vcf_to_plink_chrY_single_pulse  as convert_vcf_to_plink_mtDNA_constant with:
    input:
        vcf=rules.compress_vcf_constant.output.vcf_mtDNA,
        sex=rules.generate_vcf_constant.output.sex
    output:
        fam=expand("{data_dir}/constant_admixture_g{gen}_mtDNA.fam",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        bed=expand("{data_dir}/constant_admixture_g{gen}_mtDNA.bed",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenval=expand("{data_dir}/constant_admixture_g{gen}_mtDNA.eigenval",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenvec=expand("{data_dir}/constant_admixture_g{gen}_mtDNA.eigenvec",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    params:
        out=expand("{data_dir}/constant_admixture_g{gen}_mtDNA",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))


use rule convert_vcf_to_plink_chrY_single_pulse as convert_vcf_to_plink_chrY_non_random_mating with:
    input:
        vcf=rules.compress_vcf_non_random_mating.output.vcf_chrY,
        sex=rules.generate_vcf_non_random_mating.output.sex
    output:
        fam=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrY.fam",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        bed=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrY.bed",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenval=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrY.eigenval",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenvec=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrY.eigenvec",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    params:
        out=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrY",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))

use rule convert_vcf_to_plink_chrY_single_pulse  as convert_vcf_to_plink_mtDNA_non_random_mating with:
    input:
        vcf=rules.compress_vcf_non_random_mating.output.vcf_mtDNA,
        sex=rules.generate_vcf_non_random_mating.output.sex
    output:
        fam=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_mtDNA.fam",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        bed=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_mtDNA.bed",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenval=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_mtDNA.eigenval",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenvec=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_mtDNA.eigenvec",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    params:
        out=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_mtDNA",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))


use rule convert_vcf_to_plink_chrY_single_pulse as convert_vcf_to_plink_chrY_non_random_mating_constant with:
    input:
        vcf=rules.compress_vcf_non_random_mating_constant.output.vcf_chrY,
        sex=rules.generate_vcf_non_random_mating_constant.output.sex
    output:
        fam=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrY.fam",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        bed=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrY.bed",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenval=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrY.eigenval",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenvec=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrY.eigenvec",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    params:
        out=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrY",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))

use rule convert_vcf_to_plink_chrY_single_pulse  as convert_vcf_to_plink_mtDNA_non_random_mating_constant with:
    input:
        vcf=rules.compress_vcf_non_random_mating_constant.output.vcf_mtDNA,
        sex=rules.generate_vcf_non_random_mating_constant.output.sex
    output:
        fam=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_mtDNA.fam",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        bed=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_mtDNA.bed",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenval=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_mtDNA.eigenval",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1)),
        eigenvec=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_mtDNA.eigenvec",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))
    params:
        out=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_mtDNA",data_dir=data_dir,gen=np.arange(first_gen, last_gen + 1))

rule assign_haplogroups_mtDNA_single_pulse:
    input:
        rules.convert_vcf_to_plink_mtDNA_single_pulse.output.eigenvec
    output:
        expand("{data_dir}/single_pulse_admixture_g{gen}_mtDNA_haplogroups.tab", data_dir=data_dir,
               gen=np.arange(first_gen, last_gen +1))
    log: 'log/classifying_mtDNA_haplogroups_single_pulse.log'
    threads: 1
    conda: 'envs/haplogroups.yml'
    shell:
        "scripts/classify_mtDNA_Y_haplogroups.py --eigenvec {input} --log {log} --output_file {output}"

use rule assign_haplogroups_mtDNA_single_pulse as assign_haplogroups_chrY_single_pulse with:
    input:
        rules.convert_vcf_to_plink_chrY_single_pulse.output.eigenvec
    output:
        expand("{data_dir}/single_pulse_admixture_g{gen}_chrY_haplogroups.tab", data_dir=data_dir,
               gen=np.arange(first_gen, last_gen +1))
    log: 'log/classifying_chrY_haplogroups_single_pulse.log'

use rule assign_haplogroups_mtDNA_single_pulse as assign_haplogroups_mtDNA_constant with:
    input:
        rules.convert_vcf_to_plink_mtDNA_constant.output.eigenvec
    output:
        expand("{data_dir}/constant_admixture_g{gen}_mtDNA_haplogroups.tab", data_dir=data_dir,
               gen=np.arange(first_gen, last_gen +1))
    log: 'log/classifying_mtDNA_haplogroups_constant.log'

use rule assign_haplogroups_mtDNA_single_pulse as assign_haplogroups_chrY_constant with:
    input:
        rules.convert_vcf_to_plink_chrY_constant.output.eigenvec
    output:
        expand("{data_dir}/constant_admixture_g{gen}_chrY_haplogroups.tab", data_dir=data_dir,
               gen=np.arange(first_gen, last_gen +1))
    log: 'log/classifying_chrY_haplogroups_constant.log'

use rule assign_haplogroups_mtDNA_single_pulse as assign_haplogroups_mtDNA_non_random_mating with:
    input:
        rules.convert_vcf_to_plink_mtDNA_non_random_mating.output.eigenvec
    output:
        expand("{data_dir}/single_pulse_non_random_mating_g{gen}_mtDNA_haplogroups.tab", data_dir=data_dir,
               gen=np.arange(first_gen, last_gen +1))
    log: 'log/classifying_mtDNA_haplogroups_non_random_mating.log'

use rule assign_haplogroups_mtDNA_single_pulse as assign_haplogroups_chrY_non_random_mating with:
    input:
        rules.convert_vcf_to_plink_chrY_non_random_mating.output.eigenvec
    output:
        expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrY_haplogroups.tab", data_dir=data_dir,
               gen=np.arange(first_gen, last_gen +1))
    log: 'log/classifying_chrY_haplogroups_non_random_mating.log'

use rule assign_haplogroups_mtDNA_single_pulse as assign_haplogroups_mtDNA_non_random_mating_constant with:
    input:
        rules.convert_vcf_to_plink_mtDNA_non_random_mating_constant.output.eigenvec
    output:
        expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_mtDNA_haplogroups.tab", data_dir=data_dir,
               gen=np.arange(first_gen, last_gen +1))
    log: 'log/classifying_mtDNA_haplogroups_non_random_mating_constant.log'

use rule assign_haplogroups_mtDNA_single_pulse as assign_haplogroups_chrY_non_random_mating_constant with:
    input:
        rules.convert_vcf_to_plink_chrY_non_random_mating_constant.output.eigenvec
    output:
        expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrY_haplogroups.tab", data_dir=data_dir,
               gen=np.arange(first_gen, last_gen +1))
    log: 'log/classifying_chrY_haplogroups_non_random_mating_constant.log'

rule plot_pca_single_pulse:
    input:
        eigenvec_1 = rules.convert_vcf_to_plink_chr1_single_pulse.output.eigenvec,
       	eigenvec_x = rules.convert_vcf_to_plink_chrX_single_pulse.output.eigenvec,
    output:
        expand("{plot_dir}/{admixture_type}_admixture_g{gen}_{chrom}_pca.png",
               chrom=['chr1', 'chrX'], admixture_type='single_pulse', plot_dir=plot_dir,
               gen=np.arange(first_gen, last_gen +1))
    shell:
        "scripts/plot_pca.py --eigenvec {input.eigenvec_1} {input.eigenvec_x} --plot_dir {plot_dir}"


use rule plot_pca_single_pulse as plot_pca_constant with:
    input:
        eigenvec_1 = rules.convert_vcf_to_plink_chr1_constant.output.eigenvec,
        eigenvec_x = rules.convert_vcf_to_plink_chrX_constant.output.eigenvec,
    output:
        expand("{plot_dir}/{admixture_type}_admixture_g{gen}_{chrom}_pca.png",
               chrom=['chr1', 'chrX'], admixture_type='constant', plot_dir=plot_dir,
               gen=np.arange(first_gen, last_gen +1))


use rule plot_pca_single_pulse as plot_pca_non_random_mating with:
    input:
        eigenvec_1=rules.convert_vcf_to_plink_chr1_non_random_mating.output.eigenvec,
        eigenvec_x=rules.convert_vcf_to_plink_chrX_non_random_mating.output.eigenvec,
    output:
        expand("{plot_dir}/{admixture_type}_g{gen}_{chrom}_pca.png",
               chrom=['chr1', 'chrX'], admixture_type='single_pulse_non_random_mating', plot_dir=plot_dir,
               gen=np.arange(first_gen, last_gen +1))

use rule plot_pca_single_pulse as plot_pca_non_random_mating_constant with:
    input:
        eigenvec_1=rules.convert_vcf_to_plink_chr1_non_random_mating_constant.output.eigenvec,
        eigenvec_x=rules.convert_vcf_to_plink_chrX_non_random_mating_constant.output.eigenvec,
    output:
        expand("{plot_dir}/{admixture_type}_g{gen}_{chrom}_pca.png",
               chrom=['chr1', 'chrX'], admixture_type='constant_admixture_non_random_mating', plot_dir=plot_dir,
               gen=np.arange(first_gen, last_gen +1))

rule run_admixture_chr1_single_pulse:
    input:
        bed=rules.convert_vcf_to_plink_chr1_single_pulse.output.bed,
        population_file=rules.generate_vcf_single_pulse.output.pop_chr1
    output:
        p_file_single_pulse=expand("{data_dir}/single_pulse_admixture_g{gen}_chr1.3.P", data_dir=data_dir,
                                   gen=np.arange(first_gen, last_gen +1)),
        q_file_single_pulse=expand("{data_dir}/single_pulse_admixture_g{gen}_chr1.3.Q", data_dir=data_dir,
                                   gen=np.arange(first_gen, last_gen +1)),
    conda: "envs/admixture.yaml"
    threads: 8
    shell:
        "admixture --supervised -j{threads} {input.bed} 3; mv *3.Q {data_dir}; mv *3.P {data_dir}"

use rule run_admixture_chr1_single_pulse as run_admixture_chr1_constant with:
    input:
        bed=rules.convert_vcf_to_plink_chr1_constant.output.bed,
        population_file=rules.generate_vcf_constant.output.pop_chr1
    output:
        p_file_constant=expand("{data_dir}/constant_admixture_g{gen}_chr1.3.P", data_dir=data_dir,
                               gen=np.arange(first_gen, last_gen +1)),
        q_file_constant=expand("{data_dir}/constant_admixture_g{gen}_chr1.3.Q", data_dir=data_dir,
                               gen=np.arange(first_gen, last_gen +1))

use rule run_admixture_chr1_single_pulse as run_admixture_chr1_non_random_mating with:
    input:
        bed=rules.convert_vcf_to_plink_chr1_non_random_mating.output.bed,
        population_file=rules.generate_vcf_non_random_mating.output.pop_chr1
    output:
        p_file_non_random_mating=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chr1.3.P", data_dir=data_dir,
                                         gen=np.arange(first_gen, last_gen +1)),
        q_file_non_random_mating=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chr1.3.Q", data_dir=data_dir,
                                         gen=np.arange(first_gen, last_gen +1))

use rule run_admixture_chr1_single_pulse as run_admixture_chr1_non_random_mating_constant with:
    input:
        bed=rules.convert_vcf_to_plink_chr1_non_random_mating_constant.output.bed,
        population_file=rules.generate_vcf_non_random_mating_constant.output.pop_chr1
    output:
        p_file_non_random_mating=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chr1.3.P", data_dir=data_dir,
                                         gen=np.arange(first_gen, last_gen +1)),
        q_file_non_random_mating=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chr1.3.Q", data_dir=data_dir,
                                         gen=np.arange(first_gen, last_gen +1))

rule run_admixture_chrX_single_pulse:
    input:
        bed=rules.convert_vcf_to_plink_chrX_single_pulse.output.bed,
        population_file=rules.generate_vcf_single_pulse.output.pop_chrX
    output:
        p_file_single_pulse=expand("{data_dir}/single_pulse_admixture_g{gen}_chrX.3.P", data_dir=data_dir,
                                   gen=np.arange(first_gen, last_gen +1)),
        q_file_single_pulse=expand("{data_dir}/single_pulse_admixture_g{gen}_chrX.3.Q", data_dir=data_dir,
                                   gen=np.arange(first_gen, last_gen +1)),
    conda: "envs/admixture.yaml"
    threads: 1
    shell:
        "admixture --supervised --haploid=\"male:23\" {input.bed} 3; "
        "mv *3.Q {data_dir}; mv *3.P {data_dir}"

use rule run_admixture_chrX_single_pulse as run_admixture_chrX_constant with:
    input:
        bed=rules.convert_vcf_to_plink_chrX_constant.output.bed,
        population_file=rules.generate_vcf_constant.output.pop_chrX
    output:
        p_file_constant=expand("{data_dir}/constant_admixture_g{gen}_chrX.3.P", data_dir=data_dir,
                               gen=np.arange(first_gen, last_gen +1)),
        q_file_constant=expand("{data_dir}/constant_admixture_g{gen}_chrX.3.Q", data_dir=data_dir,
                               gen=np.arange(first_gen, last_gen +1))

use rule run_admixture_chrX_single_pulse as run_admixture_chrX_non_random_mating with:
    input:
        bed=rules.convert_vcf_to_plink_chrX_non_random_mating.output.bed,
        population_file=rules.generate_vcf_non_random_mating.output.pop_chrX
    output:
        p_file_non_random_mating=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrX.3.P", data_dir=data_dir,
                               gen=np.arange(first_gen, last_gen +1)),
        q_file_non_random_mating=expand("{data_dir}/single_pulse_non_random_mating_g{gen}_chrX.3.Q", data_dir=data_dir,
                               gen=np.arange(first_gen, last_gen +1))

use rule run_admixture_chrX_single_pulse as run_admixture_chrX_non_random_mating_constant with:
    input:
        bed=rules.convert_vcf_to_plink_chrX_non_random_mating_constant.output.bed,
        population_file=rules.generate_vcf_non_random_mating_constant.output.pop_chrX
    output:
        p_file_non_random_mating=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrX.3.P", data_dir=data_dir,
                               gen=np.arange(first_gen, last_gen +1)),
        q_file_non_random_mating=expand("{data_dir}/constant_admixture_non_random_mating_g{gen}_chrX.3.Q", data_dir=data_dir,
                               gen=np.arange(first_gen, last_gen +1))

rule summarize_results_single_pulse:
    input:
        autosome_q=rules.run_admixture_chr1_single_pulse.output.q_file_single_pulse,
        x_q=rules.run_admixture_chrX_single_pulse.output.q_file_single_pulse,
        mtDNA=rules.assign_haplogroups_mtDNA_single_pulse.output,
        chrY=rules.assign_haplogroups_chrY_single_pulse.output,
        fam=rules.convert_vcf_to_plink_chr1_single_pulse.output.fam
    output:
        autosomes_x=expand("{results_dir}/single_pulse_admixture_g{gen}_ancestry_proportions.txt", results_dir=results_dir,
                           gen=np.arange(first_gen, last_gen +1)),
        haplogroups=expand("{results_dir}/single_pulse_admixture_g{gen}_haplogroup_frequencies.txt",
                           results_dir=results_dir, gen=np.arange(first_gen, last_gen +1)),
        bootstrapped=expand("{results_dir}/single_pulse_admixture_g{gen}_bootstrapped.txt",
                           results_dir=results_dir, gen=np.arange(first_gen, last_gen +1))
    log: "log/summarizing_results.log"
    conda: "envs/haplogroups.yml"
    threads: 1
    shell:
        "scripts/summarize_results.py -hy {input.chrY} -hm {input.mtDNA} "
        "-qa {input.autosome_q} -qx {input.x_q} -f {input.fam} -rax {output.autosomes_x} "
        "-rmy {output.haplogroups} -rb {output.bootstrapped} -p {plot_dir} --log {log}"

use rule summarize_results_single_pulse as summarize_results_constant with:
    input:
        autosome_q=rules.run_admixture_chr1_constant.output.q_file_constant,
        x_q=rules.run_admixture_chrX_constant.output.q_file_constant,
        mtDNA=rules.assign_haplogroups_mtDNA_constant.output,
        chrY=rules.assign_haplogroups_chrY_constant.output,
        fam=rules.convert_vcf_to_plink_chr1_constant.output.fam
    output:
        autosomes_x=expand("{results_dir}/constant_admixture_g{gen}_ancestry_proportions.txt",
                           results_dir=results_dir, gen=np.arange(first_gen, last_gen +1)),
        haplogroups=expand("{results_dir}/constant_admixture_g{gen}_haplogroup_frequencies.txt", results_dir=results_dir,
                           gen=np.arange(first_gen, last_gen +1)),
       	bootstrapped=expand("{results_dir}/constant_admixture_g{gen}_bootstrapped.txt",
                           results_dir=results_dir, gen=np.arange(first_gen, last_gen +1))

use rule summarize_results_single_pulse as summarize_results_non_random_mating with:
    input:
        autosome_q=rules.run_admixture_chr1_non_random_mating.output.q_file_non_random_mating,
        x_q=rules.run_admixture_chrX_non_random_mating.output.q_file_non_random_mating,
        mtDNA=rules.assign_haplogroups_mtDNA_non_random_mating.output,
        chrY=rules.assign_haplogroups_chrY_non_random_mating.output,
        fam=rules.convert_vcf_to_plink_chr1_non_random_mating.output.fam
    output:
        autosomes_x=expand("{results_dir}/single_pulse_non_random_mating_g{gen}_ancestry_proportions.txt",
                           results_dir=results_dir, gen=np.arange(first_gen, last_gen +1)),
        haplogroups=expand("{results_dir}/single_pulse_non_random_mating_g{gen}_haplogroup_frequencies.txt", results_dir=results_dir,
                           gen=np.arange(first_gen, last_gen +1)),
       	bootstrapped=expand("{results_dir}/single_pulse_non_random_mating_g{gen}_bootstrapped.txt",
                           results_dir=results_dir, gen=np.arange(first_gen, last_gen +1))

use rule summarize_results_single_pulse as summarize_results_non_random_mating_constant with:
    input:
        autosome_q=rules.run_admixture_chr1_non_random_mating_constant.output.q_file_non_random_mating,
        x_q=rules.run_admixture_chrX_non_random_mating_constant.output.q_file_non_random_mating,
        mtDNA=rules.assign_haplogroups_mtDNA_non_random_mating_constant.output,
        chrY=rules.assign_haplogroups_chrY_non_random_mating_constant.output,
        fam=rules.convert_vcf_to_plink_chr1_non_random_mating_constant.output.fam
    output:
        autosomes_x=expand("{results_dir}/constant_admixture_non_random_mating_g{gen}_ancestry_proportions.txt",
                           results_dir=results_dir, gen=np.arange(first_gen, last_gen +1)),
        haplogroups=expand("{results_dir}/constant_admixture_non_random_mating_g{gen}_haplogroup_frequencies.txt", results_dir=results_dir,
                           gen=np.arange(first_gen, last_gen +1)),
       	bootstrapped=expand("{results_dir}/constant_admixture_non_random_mating_g{gen}_bootstrapped.txt",
                           results_dir=results_dir, gen=np.arange(first_gen, last_gen +1))

