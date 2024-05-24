process VCF2MAF_UNGZ {
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
    println "Running VCF2MAF UNGZ module for ${meta.id}"
    def out_maf = "${input_vcf}".replaceAll(/vcf/, "maf")
    def out_maf_gz = "${input_vcf}".replaceAll(/vcf/, "maf.gz")
    println "input vcf: ${input_vcf}"
    println "output maf   :${out_maf}"
    println "output maf gz:${out_maf_gz}"
    """
    echo 'Input VCF file ${input_vcf} is already de-compressed'
    
    num_samples=\$(bcftools query -l ${input_vcf} | wc -l)
    tumor_id=\$(bcftools query -l ${input_vcf} | head -1)
    suffix=""

    if [ "\${num_samples}" -eq 2 ]; then
        echo "VCF with tumor-normal pair found.."
        normal_id=\$(bcftools query -l ${input_vcf} | head -2 | tail -1)
        suffix="--vcf-tumor-id \${tumor_id} --vcf-normal-id \${normal_id}"
    else
        suffix="--vcf-tumor-id \${tumor_id}"
    fi

    vcf2maf_cmd="vcf2maf.pl \\
     --species ${species} \\
     --ncbi-build ${assembly} \\
     --ref-fasta ${ref_fasta} \\
     --inhibit-vep \\
     --verbose \\
     --input-vcf ${input_vcf} \\
     --output-maf ${out_maf} \\
     \${suffix}"


    echo "----------------------------"
    echo "\${vcf2maf_cmd}"
    echo "----------------------------"

    # Execute the vcf2maf command
    eval \${vcf2maf_cmd}
    
    echo 'Output MAF ${out_maf} compressing...' 
    gzip -c ${out_maf} > ${out_maf_gz}
    echo 'Output MAF ${out_maf_gz} compressed successfully'
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: 1.6.21
        bcftools: \$(echo \$(bcftools --version|head -1|cut -d' ' -f2))
    END_VERSIONS
    """
}
