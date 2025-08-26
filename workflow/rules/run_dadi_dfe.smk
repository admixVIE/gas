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


demog_params = {
    "two_epoch": {
        "p0": "0.5 0.5 0.5",
        "ubounds": "10 10 1",
        "lbounds": "0.01 0.01 0.01",
    },
    "split_mig": {
        "p0": "0.5 0.5 0.5 0.5 0.5",
        "ubounds": "10 10 10 10 1",
        "lbounds": "0.01 0.01 0.01 0.01 0.01",
    },
}

dfe_params = {
    "lognormal": {
        "p0": "5 5 0.5",
        "ubounds": "100 100 1",
        "lbounds": "0.01 0.01 0.01",
    },
}


rule create_pop_info:
    input:
        metadata=rules.download_greatapes_metadata.output.txt,
    output:
        pop_info="results/dadi/{species}/dfe/{ppl}/pop.list",
    shell:
        """
        grep -v captive {input.metadata} | awk -v pop={wildcards.ppl} '$2 ~ ("^" pop "$"){{print $4"\\t"$2}}' > {output.pop_info}
        """


rule generate_1pop_fs:
    input:
        syn_vcf="results/processed_data/{species}/1pop/{ppl}/{ppl}.biallelic.syn.snps.vcf.gz",
        nonsyn_vcf="results/processed_data/{species}/1pop/{ppl}/{ppl}.biallelic.nonsyn.snps.vcf.gz",
        pop_info=rules.create_pop_info.output.pop_info,
    output:
        syn_fs="results/dadi/{species}/dfe/{ppl}/fs/{ppl}.syn.unfolded.dadi.fs",
        nonsyn_fs="results/dadi/{species}/dfe/{ppl}/fs/{ppl}.nonsyn.unfolded.dadi.fs",
    shell:
        """
        dadi-cli GenerateFs --vcf {input.syn_vcf} --pop-ids {wildcards.ppl} --pop-info {input.pop_info} --projections $(awk 'END {{print NR * 2}}' {input.pop_info}) --output {output.syn_fs} --polarized --mask-singletons
        dadi-cli GenerateFs --vcf {input.nonsyn_vcf} --pop-ids {wildcards.ppl} --pop-info {input.pop_info} --projections $(awk 'END {{print NR * 2}}' {input.pop_info}) --output {output.nonsyn_fs} --polarized --mask-singletons
        """


rule infer_1pop_dm_warm_up:
    input:
        fs=rules.generate_1pop_fs.output.syn_fs,
    output:
        p0="results/dadi/{species}/dfe/{ppl}/InferDM/{ppl}.{demog}.InferDM.opts.0",
    params:
        p0=lambda wildcards: demog_params[wildcards.demog]["p0"],
        ubounds=lambda wildcards: demog_params[wildcards.demog]["ubounds"],
        lbounds=lambda wildcards: demog_params[wildcards.demog]["lbounds"],
        output_prefix="results/dadi/{species}/dfe/{ppl}/InferDM/{ppl}.{demog}",
        grid_size="300 400 500",
        optimizations=100,
    resources:
        time=4320,
        cpus=32,
        mem_gb=32,
    shell:
        """
        dadi-cli InferDM --fs {input.fs} --model {wildcards.demog} --p0 {params.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} --grids {params.grid_size} --output-prefix {params.output_prefix} --cpus {resources.cpus} --optimizations {params.optimizations}
        """


rule infer_1pop_dm_fine_tune:
    input:
        fs=rules.generate_1pop_fs.output.syn_fs,
        p0=rules.infer_1pop_dm_warm_up.output.p0,
    output:
        bestfit="results/dadi/{species}/dfe/{ppl}/InferDM/{ppl}.{demog}.InferDM.bestfits",
    params:
        ubounds=lambda wildcards: demog_params[wildcards.demog]["ubounds"],
        lbounds=lambda wildcards: demog_params[wildcards.demog]["lbounds"],
        output_prefix="results/dadi/{species}/dfe/{ppl}/InferDM/{ppl}.{demog}",
        grid_size="300 400 500",
        optimizations=100,
    resources:
        time=4320,
        cpus=32,
        mem_gb=32,
    shell:
        """
        dadi-cli InferDM --fs {input.fs} --model {wildcards.demog} --bestfit-p0-file {input.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} --grids {params.grid_size} --output-prefix {params.output_prefix} --cpus {resources.cpus} --force-convergence {params.optimizations}
        """


rule generate_1d_cache:
    input:
        bestfit=rules.infer_1pop_dm_fine_tune.output.bestfit,
        pop_info=rules.create_pop_info.output.pop_info,
    output:
        cache="results/dadi/{species}/dfe/{ppl}/InferDFE/{ppl}.{demog}.spectra.bpkl",
    params:
        grid_size="300 400 500",
        gamma_pts=50,
    resources:
        time=4320,
        cpus=32,
        mem_gb=32,
    shell:
        """
        dadi-cli GenerateCache --model {wildcards.demog}_sel --demo-popt {input.bestfit} --sample-size $(awk 'END {{print NR * 2}}' {input.pop_info}) --grids {params.grid_size} --gamma-pts {params.gamma_pts} --output {output.cache} --cpus {resources.cpus}
        """


rule infer_dfe_warm_up:
    input:
        fs=rules.generate_1pop_fs.output.nonsyn_fs,
        cache=rules.generate_1d_cache.output.cache,
        dm_bestfit=rules.infer_1pop_dm_fine_tune.output.bestfit,
    output:
        p0="results/dadi/{species}/dfe/{ppl}/InferDFE/{ppl}.{demog}.{dfe}.InferDFE.opts.0",
    params:
        p0=lambda wildcards: dfe_params[wildcards.dfe]["p0"],
        ubounds=lambda wildcards: dfe_params[wildcards.dfe]["ubounds"],
        lbounds=lambda wildcards: dfe_params[wildcards.dfe]["lbounds"],
        ratio=2.31,
        output_prefix="results/dadi/{species}/dfe/{ppl}/InferDFE/{ppl}.{demog}.{dfe}",
        optimizations=100,
    resources:
        time=4320,
        cpus=32,
        mem_gb=32,
    shell:
        """
        dadi-cli InferDFE --fs {input.fs} --cache1d {input.cache} --demo-popt {input.dm_bestfit} --output-prefix {params.output_prefix} --pdf1d {wildcards.dfe} --p0 {params.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} --ratio {params.ratio} --cpus {resources.cpus} --optimizations {params.optimizations}
        """


rule infer_dfe_fine_tune:
    input:
        fs=rules.generate_1pop_fs.output.nonsyn_fs,
        cache=rules.generate_1d_cache.output.cache,
        dm_bestfit=rules.infer_1pop_dm_fine_tune.output.bestfit,
        p0=rules.infer_dfe_warm_up.output.p0,
    output:
        bestfit="results/dadi/{species}/dfe/{ppl}/InferDFE/{ppl}.{demog}.{dfe}.InferDFE.bestfits",
    params:
        ubounds=lambda wildcards: dfe_params[wildcards.dfe]["ubounds"],
        lbounds=lambda wildcards: dfe_params[wildcards.dfe]["lbounds"],
        ratio=2.31,
        output_prefix="results/dadi/{species}/dfe/{ppl}/InferDFE/{ppl}.{demog}.{dfe}",
        optimizations=100,
    resources:
        time=4320,
        cpus=32,
        mem_gb=32,
    shell:
        """
        dadi-cli InferDFE --fs {input.fs} --cache1d {input.cache} --demo-popt {input.dm_bestfit} --output-prefix {params.output_prefix} --pdf1d {wildcards.dfe} --bestfit-p0-file {input.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} --ratio {params.ratio} --cpus {resources.cpus} --force-convergence {params.optimizations}
        """


rule dfe_godambe_ci:
    input:
        syn_vcf="results/processed_data/{species}/1pop/{ppl}/{ppl}.biallelic.syn.snps.vcf.gz",
        nonsyn_vcf="results/processed_data/{species}/1pop/{ppl}/{ppl}.biallelic.nonsyn.snps.vcf.gz",
        pop_info=rules.create_pop_info.output.pop_info,
        nonsyn_fs=rules.generate_1pop_fs.output.nonsyn_fs,
        cache=rules.generate_1d_cache.output.cache,
        dfe_bestfit=rules.infer_dfe_fine_tune.output.bestfit,
    output:
        dfe_godambe_ci="results/dadi/{species}/dfe/{ppl}/StatDFE/{ppl}.{demog}.{dfe}.godambe.ci",
    params:
        syn_dir="results/dadi/{species}/dfe/{ppl}/StatDFE/{ppl}_bootstrapping_syn",
        nonsyn_dir="results/dadi/{species}/dfe/{ppl}/StatDFE/{ppl}_bootstrapping_non",
        syn_output_prefix="results/dadi/{species}/dfe/{ppl}/StatDFE/{ppl}_bootstrapping_syn/{ppl}.syn.unfolded",
        nonsyn_output_prefix="results/dadi/{species}/dfe/{ppl}/StatDFE/{ppl}_bootstrapping_non/{ppl}.nonsyn.unfolded",
    shell:
        """
        [ -d {params.syn_dir} ] || mkdir {params.syn_dir}
        [ -d {params.nonsyn_dir} ] || mkdir {params.nonsyn_dir}
        dadi-cli GenerateFs --vcf {input.syn_vcf} --pop-info {input.pop_info} --pop-ids {wildcards.ppl} --projections $(awk 'END {{print NR * 2}}' {input.pop_info}) --polarized --bootstrap 100 --chunk-size 10000000 --output {params.syn_output_prefix} --mask-singletons
        dadi-cli GenerateFs --vcf {input.nonsyn_vcf} --pop-info {input.pop_info} --pop-ids {wildcards.ppl} --projections $(awk 'END {{print NR * 2}}' {input.pop_info}) --polarized --bootstrap 100 --chunk-size 10000000 --output {params.nonsyn_output_prefix} --mask-singletons
        dadi-cli StatDFE --fs {input.nonsyn_fs} --dfe-popt {input.dfe_bestfit} --cache1d {input.cache} --pdf1d {wildcards.dfe} --bootstrapping-nonsynonymous-dir {params.nonsyn_dir} --bootstrapping-synonymous-dir {params.syn_dir} --output {output.dfe_godambe_ci}
        """
