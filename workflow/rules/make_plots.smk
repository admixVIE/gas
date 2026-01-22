# Copyright 2025 Xin Huang and Simon Chen
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


ruleorder: plot_fitted_dfe > plot_fitted_dm
ruleorder: plot_fitted_2pop_dfe > plot_fitted_2pop_dm
ruleorder: plot_mutation_proportion > plot_fitted_dm
ruleorder: plot_mutation_proportion > plot_fitted_dfe


rule plot_selscan_xp:
    input:
        scores=rules.merge_selscan_xp_scores.output.merged_scores,
    output:
        candidates1="results/selscan/{species}/2pop/{pair}/{method}_{maf}/candidates/{pair}.normalized.{method}.top.0.00005.candidates.scores",
        candidates2="results/selscan/{species}/2pop/{pair}/{method}_{maf}/candidates/{pair}.normalized.{method}.top.0.0005.candidates.scores",
        plot="results/plots/selscan/{species}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.scores.png",
    params:
        score_column="normalized_{method}",
        top_fraction1=0.00005,
        top_fraction2=0.0005,
        use_absolute="TRUE",
        width=640,
        height=240,
        title="{pair}",
    resources:
        mem_gb=16,
    shell:
        """
        Rscript workflow/scripts/manhattan.R {input.scores} {params.score_column} {params.top_fraction1} {params.top_fraction2} {params.use_absolute} {params.width} {params.height} {params.title} {output.candidates1} {output.candidates2} {output.plot}
        """


rule plot_betascan:
    input:
        scores=rules.merge_b1_scores.output.merged_scores,
    output:
        candidates1="results/betascan/{species}/{ppl}/m_{core_freq}/candidates/{ppl}.b1.top.0.00005.candidates.scores",
        candidates2="results/betascan/{species}/{ppl}/m_{core_freq}/candidates/{ppl}.b1.top.0.0005.candidates.scores",
        plot="results/plots/betascan/{species}/{ppl}/m_{core_freq}/{ppl}.b1.scores.png",
    params:
        score_column="B1",
        top_fraction1=0.00005,
        top_fraction2=0.0005,
        use_absolute="FALSE",
        width=640,
        height=240,
        title="{ppl}",
    resources:
        mem_gb=16,
    shell:
        """
        Rscript workflow/scripts/manhattan.R {input.scores} {params.score_column} {params.top_fraction1} {params.top_fraction2} {params.use_absolute} {params.width} {params.height} {params.title} {output.candidates1} {output.candidates2} {output.plot}
        """


rule plot_fitted_dm:
    input:
        fs=rules.generate_1pop_fs.output.syn_fs,
        dm_popt=rules.infer_1pop_dm_fine_tune.output.bestfit,
        pop_info=rules.create_pop_info.output.pop_info,
    output:
        fs_plot="results/plots/dadi/{species}/dfe/{ppl}/{ppl}.{demog}.png",
    shell:
        """
        dadi-cli Plot --fs {input.fs} --demo-popt {input.dm_popt} --model {wildcards.demog} --projections $(awk 'END {{print NR * 2}}' {input.pop_info}) --output {output.fs_plot}
        """


rule plot_fitted_2pop_dm:
    input:
        fs=rules.generate_2pop_fs.output.syn_fs,
        dm_popt=rules.infer_2pop_dm_fine_tune.output.bestfit,
        pop_info=rules.create_pair_info.output.pair_info,
    output:
        fs_plot="results/plots/dadi/{species}/jdfe/{pair}/{pair}.{demog}.png",
    shell:
        """
        pair="{wildcards.pair}"
        pop1=$(echo $pair | cut -d'_' -f1)
        pop2=$(echo $pair | cut -d'_' -f2)

        dadi-cli Plot --fs {input.fs} --demo-popt {input.dm_popt} --model {wildcards.demog} --projections $(grep -w $pop1 {input.pop_info} | awk 'END {{print NR * 2}}') $(grep -w $pop2 {input.pop_info} | awk 'END {{print NR * 2}}') --output {output.fs_plot}
        """


rule plot_fitted_dfe:
    input:
        fs=rules.generate_1pop_fs.output.nonsyn_fs,
        dfe_popt=rules.infer_dfe_fine_tune.output.bestfit,
        cache=rules.generate_1d_cache.output.cache,
        pop_info=rules.create_pop_info.output.pop_info,
    output:
        fs_plot="results/plots/dadi/{species}/dfe/{ppl}/{ppl}.{demog}.{dfe}.png",
    shell:
        """
        dadi-cli Plot --fs {input.fs} --cache1d {input.cache} --dfe-popt {input.dfe_popt} --pdf1d {wildcards.dfe} --projections $(awk 'END {{print NR * 2}}' {input.pop_info}) --output {output.fs_plot}
        """


rule plot_fitted_2pop_dfe:
    input:
        fs=rules.generate_2pop_fs.output.nonsyn_fs,
        dfe_popt=rules.infer_jdfe_fine_tune.output.bestfit,
        cache=rules.generate_2d_cache.output.cache,
        pop_info=rules.create_pair_info.output.pair_info,
    output:
        fs_plot="results/plots/dadi/{species}/jdfe/{pair}/{pair}.{demog}.biv_{dfe}.png",
    shell:
        """
        pair="{wildcards.pair}"
        pop1=$(echo $pair | cut -d'_' -f1)
        pop2=$(echo $pair | cut -d'_' -f2)

        dadi-cli Plot --fs {input.fs} --cache2d {input.cache} --dfe-popt {input.dfe_popt} --pdf2d biv_{wildcards.dfe} --projections $(grep -w $pop1 {input.pop_info} | awk 'END {{print NR * 2}}') $(grep -w $pop2 {input.pop_info} | awk 'END {{print NR * 2}}') --output {output.fs_plot}
        """


rule plot_mutation_proportion:
    input:
        dfe_popt=rules.infer_dfe_fine_tune.output.bestfit,
    output:
        plot="results/plots/dadi/{species}/dfe/{ppl}/{ppl}.{demog}.{dfe}.mut.prop.png",
    script:
        "../scripts/plot_mutation_prop.py"


rule plot_selscan_functional_enrichment:
    input:
        enrichment=rules.selscan_functional_enrichment.output.enrichment,
    output:
        plot="results/plots/enrichment/selscan/{species}/lineages/{method}_{maf}/{lineage}.{method}_{maf}.top.{cutoff}.gowinda.enrichment.png",
    resources:
        mem_gb=8,
    script:
        "../scripts/plot_gowinda_enrichment.py"


rule plot_selscan_overlap_within_species:
    input:
        genes=expand(
            "results/selscan/{species}/lineages/{method}_{maf}/candidates/{lineage}.{method}_{maf}.top.{cutoff}.candidate.genes",
            species=species,
            lineage=lineages_for_species,
            allow_missing=True,
        ),
    output:
        plot="results/plots/selscan/lineages/{species}.{method}_{maf}.top.{cutoff}.candidate.genes.png",
    script:
        "../scripts/plot_upset_diagram.py"
