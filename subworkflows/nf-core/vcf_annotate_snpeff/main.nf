//
// Run SNPEFF to annotate VCF files
//

include { SNPEFF_SNPEFF    } from '../../../modules/nf-core/snpeff/snpeff/main.nf'
include { VCF2MAF_UNGZ    } from '../../../modules/local/vcf2maf_ungz/main.nf'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix/main.nf'

workflow VCF_ANNOTATE_SNPEFF {
    take:
    ch_vcf          // channel: [ val(meta), path(vcf) ]
    val_snpeff_db   // string:  db version to use
    ch_snpeff_cache // channel: [ path(cache) ] (optional)
    ch_species
    ch_genome
    ch_fasta

    main:
    ch_versions = Channel.empty()

    SNPEFF_SNPEFF(ch_vcf, val_snpeff_db, ch_snpeff_cache)
    VCF2MAF_UNGZ(SNPEFF_SNPEFF.out.vcf, ch_species, ch_genome, ch_fasta)
    TABIX_BGZIPTABIX(SNPEFF_SNPEFF.out.vcf)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(SNPEFF_SNPEFF.out.versions)
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)
    ch_versions = ch_versions.mix(VCF2MAF_UNGZ.out.versions)

    emit:
    vcf_tbi  = TABIX_BGZIPTABIX.out.gz_tbi // channel: [ val(meta), path(vcf), path(tbi) ]
    maf = VCF2MAF_UNGZ.out.maf // channel: [  path(maf) ]
    reports  = SNPEFF_SNPEFF.out.report    // channel: [ path(html) ]
    versions = ch_versions                 // channel: [ path(versions.yml) ]
}
