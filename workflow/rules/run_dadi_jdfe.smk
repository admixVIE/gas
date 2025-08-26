# Copyright 2025 Xin Huang
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


rule generate_2pop_fs:
    input:
        syn_vcf="results/processed_data/{species}/2pop/{pair}/{pair}.biallelic.syn.snps.vcf.gz",
        nonsyn_vcf="results/processed_data/{species}/2pop/{pair}/{pair}.biallelic.nonsyn.snps.vcf.gz",
        pop_info=rules.create_pair_info.output.pair_info,
    output:
        syn_fs="results/dadi/{species}/jdfe/{pair}/fs/{pair}.syn.unfolded.dadi.fs",
        nonsyn_fs="results/dadi/{species}/jdfe/{pair}/fs/{pair}.nonsyn.unfolded.dadi.fs",
    shell:
        """
        pair="{wildcards.pair}"
        pop1=$(echo $pair | cut -d'_' -f1)
        pop2=$(echo $pair | cut -d'_' -f2)

        dadi-cli GenerateFs --vcf {input.syn_vcf} --pop-ids $pop1 $pop2 --pop-info {input.pop_info} --projections $(grep -w $pop1 {input.pop_info} | awk 'END {{print NR * 2}}') $(grep -w $pop2 {input.pop_info} | awk 'END {{print NR * 2}}') --output {output.syn_fs} --polarized --mask-singletons
        dadi-cli GenerateFs --vcf {input.nonsyn_vcf} --pop-ids $pop1 $pop2 --pop-info {input.pop_info} --projections $(grep -w $pop1 {input.pop_info} | awk 'END {{print NR * 2}}') $(grep -w $pop2 {input.pop_info} | awk 'END {{print NR * 2}}') --output {output.nonsyn_fs} --polarized --mask-singletons
        """


rule infer_2pop_dm_warm_up:
    input:
        fs=rules.generate_2pop_fs.output.syn_fs,
    output:
        p0="results/dadi/{species}/jdfe/{pair}/InferDM/{pair}.{demog}.InferDM.opts.0",
    params:
        p0=lambda wildcards: demog_params[wildcards.demog]["p0"],
        ubounds=lambda wildcards: demog_params[wildcards.demog]["ubounds"],
        lbounds=lambda wildcards: demog_params[wildcards.demog]["lbounds"],
        output_prefix="results/dadi/{species}/jdfe/{pair}/InferDM/{pair}.{demog}",
        grid_size="300 400 500",
        optimizations=100,
    resources:
        time=4320,
        cpus=368,
    shell:
        """
        dadi-cli InferDM --fs {input.fs} --model {wildcards.demog} --p0 {params.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} --grids {params.grid_size} --output-prefix {params.output_prefix} --cpus {resources.cpus} --optimizations {params.optimizations}
        """


rule infer_2pop_dm_fine_tune:
    input:
        fs=rules.generate_2pop_fs.output.syn_fs,
        p0=rules.infer_2pop_dm_warm_up.output.p0,
    output:
        bestfit="results/dadi/{species}/jdfe/{pair}/InferDM/{pair}.{demog}.InferDM.bestfits",
    params:
        ubounds=lambda wildcards: demog_params[wildcards.demog]["ubounds"],
        lbounds=lambda wildcards: demog_params[wildcards.demog]["lbounds"],
        output_prefix="results/dadi/{species}/jdfe/{pair}/InferDM/{pair}.{demog}",
        grid_size="300 400 500",
        optimizations=100,
    resources:
        time=4320,
        cpus=368,
    shell:
        """
        dadi-cli InferDM --fs {input.fs} --model {wildcards.demog} --bestfit-p0-file {input.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} --grids {params.grid_size} --output-prefix {params.output_prefix} --cpus {resources.cpus} --force-convergence {params.optimizations}
        """


rule generate_2d_cache:
    input:
        bestfit=rules.infer_2pop_dm_fine_tune.output.bestfit,
        pop_info=rules.create_pair_info.output.pair_info,
    output:
        cache="results/dadi/{species}/jdfe/{pair}/InferDFE/{pair}.{demog}.spectra.bpkl",
    params:
        grid_size="300 400 500",
        gamma_pts=100,
    resources:
        time=4320,
        cpus=368,
        mem_gb=512,
    shell:
        """
        pair="{wildcards.pair}"
        pop1=$(echo $pair | cut -d'_' -f1)
        pop2=$(echo $pair | cut -d'_' -f2)

        dadi-cli GenerateCache --model {wildcards.demog}_sel --demo-popt {input.bestfit} --sample-size $(grep -w $pop1 {input.pop_info} | awk 'END {{print NR * 2}}') $(grep -w $pop2 {input.pop_info} | awk 'END {{print NR * 2}}') --grids {params.grid_size} --gamma-pts {params.gamma_pts} --output {output.cache} --cpus {resources.cpus} --dimensionality 2 
        """


rule get_dfe_bestfits:
    input:
        cache=rules.generate_2d_cache.output.cache,
    output:
        constants="results/dadi/{species}/jdfe/{pair}/InferDFE/{pair}.{demog}.constants",
    shell:
        """
        pair="{wildcards.pair}"
        pop1=$(echo $pair | cut -d'_' -f1)
        pop2=$(echo $pair | cut -d'_' -f2)
        
        read p1_v1 p1_v2 <<< $(grep Converged -A 2 results/dadi/{wildcards.species}/dfe/$pop1/InferDFE/$pop1.two_epoch.lognormal.InferDFE.bestfits | tail -1 | awk '{{printf "%.2f %.2f\\n", $2, $3}}')
        read p2_v1 p2_v2 <<< $(grep Converged -A 2 results/dadi/{wildcards.species}/dfe/$pop2/InferDFE/$pop2.two_epoch.lognormal.InferDFE.bestfits | tail -1 | awk '{{printf "%.2f %.2f\\n", $2, $3}}')
        echo "$p1_v1 $p2_v1 $p1_v2 $p2_v2" > {output.constants}
        """


rule infer_jdfe_warm_up:
    input:
        fs=rules.generate_2pop_fs.output.nonsyn_fs,
        cache2d=rules.generate_2d_cache.output.cache,
        dm_bestfit=rules.infer_2pop_dm_fine_tune.output.bestfit,
        constants=rules.get_dfe_bestfits.output.constants,
    output:
        p0="results/dadi/{species}/jdfe/{pair}/InferDFE/{pair}.{demog}.biv_{dfe}.InferDFE.opts.0",
    params:
        p0="1 1 1 1 .5 .5",
        ubounds="100 100 100 100 0.99 1",
        lbounds="0.01 0.01 0.01 0.01 -0.99 0.01",
        ratio=2.31,
        output_prefix="results/dadi/{species}/jdfe/{pair}/InferDFE/{pair}.{demog}.biv_{dfe}",
        optimizations=100,
    resources:
        time=4320,
        cpus=368,
        mem_gb=512,
    shell:
        """
        constants=$(cat {input.constants})
        dadi-cli InferDFE --fs {input.fs} --cache2d {input.cache2d} --demo-popt {input.dm_bestfit} --output-prefix {params.output_prefix} --pdf2d biv_{wildcards.dfe} --p0 {params.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} --constants $constants -1 -1 --ratio {params.ratio} --cpus {resources.cpus} --optimizations {params.optimizations}
        """


rule infer_jdfe_fine_tune:
    input:
        fs=rules.generate_2pop_fs.output.nonsyn_fs,
        cache2d=rules.generate_2d_cache.output.cache,
        dm_bestfit=rules.infer_2pop_dm_fine_tune.output.bestfit,
        constants=rules.get_dfe_bestfits.output.constants,
        p0=rules.infer_jdfe_warm_up.output.p0,
    output:
        bestfit="results/dadi/{species}/jdfe/{pair}/InferDFE/{pair}.{demog}.biv_{dfe}.InferDFE.bestfits",
    params:
        ubounds="100 100 100 100 0.99 1",
        lbounds="0.01 0.01 0.01 0.01 -0.99 0.01",
        ratio=2.31,
        output_prefix="results/dadi/{species}/jdfe/{pair}/InferDFE/{pair}.{demog}.biv_{dfe}",
        optimizations=100,
    resources:
        time=4320,
        cpus=368,
        mem_gb=512,
    shell:
        """
        constants=$(cat {input.constants})
        dadi-cli InferDFE --fs {input.fs} --cache2d {input.cache2d} --demo-popt {input.dm_bestfit} --output-prefix {params.output_prefix} --pdf2d biv_{wildcards.dfe} --bestfit-p0-file {input.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} --constants $constants -1 -1 --ratio {params.ratio} --cpus {resources.cpus} --force-convergence {params.optimizations}
        """


rule jdfe_godambe_ci:
    input:
        syn_vcf="results/processed_data/{species}/2pop/{pair}/{pair}.biallelic.syn.snps.vcf.gz",
        nonsyn_vcf="results/processed_data/{species}/2pop/{pair}/{pair}.biallelic.nonsyn.snps.vcf.gz",
        pop_info=rules.create_pair_info.output.pair_info,
        nonsyn_fs=rules.generate_2pop_fs.output.nonsyn_fs,
        cache=rules.generate_2d_cache.output.cache,
        dfe_bestfit=rules.infer_jdfe_fine_tune.output.bestfit,
        constants=rules.get_dfe_bestfits.output.constants,
    output:
        dfe_godambe_ci="results/dadi/{species}/jdfe/{pair}/StatDFE/{pair}.{demog}.biv_{dfe}.godambe.ci",
    params:
        syn_dir="results/dadi/{species}/jdfe/{pair}/StatDFE/{pair}_bootstrapping_syn",
        nonsyn_dir="results/dadi/{species}/jdfe/{pair}/StatDFE/{pair}_bootstrapping_non",
        syn_output_prefix="results/dadi/{species}/jdfe/{pair}/StatDFE/{pair}_bootstrapping_syn/{pair}.syn.unfolded",
        nonsyn_output_prefix="results/dadi/{species}/jdfe/{pair}/StatDFE/{pair}_bootstrapping_non/{pair}.nonsyn.unfolded",
    shell:
        """
        pair="{wildcards.pair}"
        pop1=$(echo $pair | cut -d'_' -f1)
        pop2=$(echo $pair | cut -d'_' -f2)
        constants=$(cat {input.constants})

        [ -d {params.syn_dir} ] || mkdir {params.syn_dir}
        [ -d {params.nonsyn_dir} ] || mkdir {params.nonsyn_dir}
        dadi-cli GenerateFs --vcf {input.syn_vcf} --pop-info {input.pop_info} --pop-ids $pop1 $pop2 --projections $(grep -w $pop1 {input.pop_info} | awk 'END {{print NR * 2}}') $(grep -w $pop2 {input.pop_info} | awk 'END {{print NR * 2}}') --polarized --bootstrap 100 --chunk-size 10000000 --output {params.syn_output_prefix} --mask-singletons
        dadi-cli GenerateFs --vcf {input.nonsyn_vcf} --pop-info {input.pop_info} --pop-ids $pop1 $pop2 --projections $(grep -w $pop1 {input.pop_info} | awk 'END {{print NR * 2}}') $(grep -w $pop2 {input.pop_info} | awk 'END {{print NR * 2}}') --polarized --bootstrap 100 --chunk-size 10000000 --output {params.nonsyn_output_prefix} --mask-singletons
        dadi-cli StatDFE --fs {input.nonsyn_fs} --dfe-popt {input.dfe_bestfit} --cache2d {input.cache} --pdf2d biv_{wildcards.dfe} --bootstrapping-nonsynonymous-dir {params.nonsyn_dir} --bootstrapping-synonymous-dir {params.syn_dir} --output {output.dfe_godambe_ci} --constant $constants -1 -1
        """
