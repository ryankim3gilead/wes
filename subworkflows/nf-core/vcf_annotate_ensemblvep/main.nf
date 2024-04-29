//
// Run VEP to annotate VCF files
//

include { ENSEMBLVEP_VEP } from '../../../modules/nf-core/ensemblvep/vep/main'
include { VCF2MAF    } from '../../../modules/local/vcf2maf/main.nf'
include { TABIX_TABIX    } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_ANNOTATE_ENSEMBLVEP {
    take:
    ch_vcf                      // channel: [ val(meta), path(vcf), [path(custom_file1), path(custom_file2)... (optionnal)]]
    ch_fasta                    // channel: [ val(meta2), path(fasta) ] (optional)
    val_genome                  //   value: genome to use
    val_species                 //   value: species to use
    val_cache_version           //   value: cache version to use
    ch_cache                    // channel: [ val(meta3), path(cache) ] (optional)
    ch_extra_files              // channel: [ path(file1), path(file2)... ] (optional)
    ch_fasta_path
    
    main:
    ch_versions = Channel.empty()

    ENSEMBLVEP_VEP(
        ch_vcf,
        val_genome,
        val_species,
        val_cache_version,
        ch_cache,
        ch_fasta,
        ch_extra_files
    )
    VCF2MAF(ENSEMBLVEP_VEP.out.vcf, val_species, val_genome, ch_fasta_path)

    TABIX_TABIX(ENSEMBLVEP_VEP.out.vcf)

    ch_vcf_tbi = ENSEMBLVEP_VEP.out.vcf.join(TABIX_TABIX.out.tbi, failOnDuplicate: true, failOnMismatch: true)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)
    ch_versions = ch_versions.mix(VCF2MAF.out.versions)

    emit:
    vcf_tbi  = ch_vcf_tbi                  // channel: [ val(meta), path(vcf), path(tbi) ]
    maf      = VCF2MAF.out.maf // channel: [  path(maf) ]
    json     = ENSEMBLVEP_VEP.out.json     // channel: [ val(meta), path(json) ]
    tab      = ENSEMBLVEP_VEP.out.tab      // channel: [ val(meta), path(tab) ]
    reports  = ENSEMBLVEP_VEP.out.report   // channel: [ path(html) ]
    versions = ch_versions                 // channel: [ versions.yml ]
}
