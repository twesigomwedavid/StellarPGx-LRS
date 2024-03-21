#!/usr/bin/env nextflow

nextflow.enable.dsl=2

gene_name = params.gene
up_gene_symbol = gene_name.toUpperCase()
d_base = params.db_init
res_base = params.res_init
caller_base = params.caller_init
output_folder = params.out_dir

params.build='hg38'

db = "${d_base}/${gene_name}/hg38"
res_dir = "${res_base}/${gene_name}/cyp_hg38"
caller_dir = "${caller_base}/${gene_name}/hg38/bin"

chrom = "chr22"
region_a1 = "chr22:42126000-42137500"
region_a2 = "042126000-042137500"
region_b1 = "chr22:42126300-42132400"
region_b2 = "042126300-042132400" 
transcript = "ENST00000645361"

params.outdir_vcfs = "$PWD/results/$gene_name/vcfs"

params.format='binary'

if (params.format=='compressed') {
    ext = "cram"
    ind = "crai"
    cram_options = "--force_use_input_ref_for_cram_reading" 

} else { 
    ext = "bam"
    ind = "bai"
    cram_options = ""

}

ref_dir_val = new File("${params.ref_file}").getParent()
ref_genome = new File("${params.ref_file}").getName()


workflow {
    ref_ch = Channel.value("${ref_dir_val}")
    reads_ch = Channel.fromFilePairs(params.in_bam, type: 'file') {  file -> file.name.replaceAll(/.${ext}|.${ind}$/,'') }
//    reads_ch.view()
    call_snvs(ref_ch, reads_ch)
    call_snvs.out.join(reads_ch).set {phase_ch}
    phase_snvs(ref_ch, phase_ch)
    get_core_var(phase_snvs.out, ref_ch, res_dir, res_base, caller_base)
    get_hap_snvs(get_core_var.out)
    get_depth(reads_ch, ref_ch, res_dir)
    query_dup_profile(phase_snvs.out)
    get_hap_snvs.out.join(get_depth.out).set{fin_files1}
    fin_files1.join(query_dup_profile.out).set{fin_files}
    call_stars(caller_dir, db, fin_files)
}

// reads_ch.view()


process call_snvs {
    tag "$name"

    // publishDir params.outdir_vcfs

 
    input:
    path ref_dir
    tuple val(name), path (reads)
 
    output:
    tuple val(name), path("*")
 
    script:

    """
    graphtyper genotype_lr ${ref_dir}/${ref_genome} --sam=${reads[0]} --region=${region_a1} --output=${name}_var_1 

    """
}


process phase_snvs {
    tag "$name"

    publishDir params.outdir_vcfs


    input:
    path ref_dir
    tuple val(name), path ("${name}_var_1"), path(bam)

    output:
    tuple val(name), path("*")

    script:

    """
    whatshap phase -o ${name}.${gene_name}.phased.vcf --reference ${ref_dir}/${ref_genome} --chromosome ${chrom} ${name}_var_1/${chrom}/${region_a2}.vcf.gz ${bam[0]} --indels
    """

}


process get_core_var {
//   maxForks 10
   
    errorStrategy 'ignore'
    tag "${name}"   
    label 'test2'

    input:
    tuple val(name), path("${name}.${gene_name}.phased.vcf") 
    path ref_dir
    path res_dir
    path res_base
    path caller_base

    output:
    tuple val(name), path("${name}_int") 

    script:
 
    """
    bgzip ${name}.${gene_name}.phased.vcf
    tabix ${name}.${gene_name}.phased.vcf.gz
    bcftools isec ${name}.${gene_name}.phased.vcf.gz ${res_dir}/allele_def_var.vcf.gz -p ${name}_int -Oz
    bcftools norm -m - ${name}_int/0002.vcf.gz > ${name}_int/${name}_core_int1.vcf
    bcftools csq -p m -v 0 -f ${ref_dir}/${ref_genome} -g ${res_base}/annotation/Homo_sapiens.GRCh38.110.gff3.gz ${name}_int/0000.vcf.gz -o ${name}_int/0000_annot.vcf
    python3 ${caller_base}/novel/core_var.py ${name}_int/0000_annot.vcf ${up_gene_symbol} ${transcript} >> ${name}_int/${name}_core_int1.vcf
    bcftools sort ${name}_int/${name}_core_int1.vcf -T ${name}_int | bgzip -c > ${name}_int/${name}_${gene_name}_core.vcf.gz
    tabix ${name}_int/${name}_${gene_name}_core.vcf.gz

    """

}


process get_hap_snvs {
//   maxForks 10

    errorStrategy 'ignore'
    tag "${name}"
    label 'test2'

    input:
    tuple val(name), path("${name}_int") 

    output:
    tuple val(name), path("*")

    script:
    """
    bcftools query -f'[%POS~%REF>%ALT\t%GT\n]' ${name}_int/${name}_${gene_name}_core.vcf.gz | grep -v '1|0' | cut -f1 > ${name}_${gene_name}_hap1.dip
    bcftools query -f'[%POS~%REF>%ALT\t%GT\n]' ${name}_int/${name}_${gene_name}_core.vcf.gz | grep -v '0|1' | cut -f1 > ${name}_${gene_name}_hap2.dip
    """

}



process get_depth {
    tag "$name"
    label 'test2'

    // publishDir params.outdir, mode: 'copy', overwrite: 'true' 

    input:
    tuple val(name), path (reads)
    path ref_dir 
    path res_dir

    output:
    tuple val(name), path ("*")

    script:

    """
    samtools bedcov --reference ${ref_dir}/${ref_genome} ${res_dir}/test3.bed ${name}.${ext} > ${name}_${gene_name}_ctrl.depth  

    """

}


process query_dup_profile {
//   maxForks 10

    errorStrategy 'ignore'
    tag "${name}"
    label 'test2'

    input:
    tuple val(name), path("${name}.${gene_name}.phased.vcf")

    output:
    tuple val(name), path("${name}_${gene_name}_dup_phased_summary.txt")

    script:

    """
    bcftools query -f'%POS~%REF>%ALT\t%QUAL\t[%GT\t%DP\t%AD]\n' -i'GT="alt"' ${name}.${gene_name}.phased.vcf > ${name}_${gene_name}_dup_phased_summary.txt

    """

}


process call_stars {
//   maxForks 10

    publishDir "$output_folder/$gene_name/alleles", mode: 'copy', overwrite: 'true'

    errorStrategy 'ignore'
    tag "${name}"
    label 'test2'

    input:
    path caller_dir
    path db
    tuple val(name), path(haps), path("${name}_${gene_name}_dp"), path("${name}_gene_dup_phased_summary.txt")

    output:
    tuple val(name), file("${name}_${gene_name}.alleles")

    script:
   
    """
    python3 ${caller_dir}/stellarpgx.py ${db}/alleles.dbs ${haps[0]} ${haps[1]} ${name}_${gene_name}_dp ${name}_gene_dup_phased_summary.txt ${db}/a_scores.dbs ${name} > ${name}_${gene_name}.alleles  

    """

}
