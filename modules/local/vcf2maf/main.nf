process VCF2MAF {
    tag "${meta.id}"
    label "process_small"
    
    container "awsecr111/bioinfo:vcf2maf"

    input:
    tuple val(meta), path(input_vcf)
    val species
    val assembly
    path ref_fasta

    output:
    path "*.maf.gz",        emit: maf
    path "versions.yml",    emit: versions

    script:
    def input_vcf_decomp = "${input_vcf}".replaceAll(/vcf.gz/, "vcf")
    def out_maf = "${input_vcf}".replaceAll(/vcf.gz/, "maf")
    def out_maf_gz = "${input_vcf}".replaceAll(/vcf.gz/, "maf.gz")
    println "Running VCF2MAF module for ${meta.id}"
    println "input vcf: ${input_vcf}"
    println "input vcf decomp: ${input_vcf_decomp}"
    println "output maf   :${out_maf}"
    println "output maf gz:${out_maf_gz}"
    """
    gzip -cd  ${input_vcf} > ${input_vcf_decomp}
    echo 'Input VCF file ${input_vcf} de-compressed successfully.'
    vcf2maf.pl \\
     --species ${species} \\
     --ncbi-build ${assembly} \\
     --ref-fasta ${ref_fasta} \\
     --inhibit-vep \\
     --verbose \\
     --input-vcf ${input_vcf_decomp} \\
     --output-maf ${out_maf}
    echo 'Output MAF ${out_maf} compressing...' 
    gzip -c ${out_maf} > ${out_maf_gz}
    echo 'Output MAF ${out_maf_gz} compressed successfully'
    """
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: 1.6.21
    END_VERSIONS
    """
}
