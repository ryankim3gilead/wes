//
// MARKDUPLICATES AND/OR QC after mapping
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { CRAM_QC_MOSDEPTH_SAMTOOLS              } from '../cram_qc_mosdepth_samtools/main'
include { GATK4_ESTIMATELIBRARYCOMPLEXITY        } from '../../../modules/nf-core/gatk4/estimatelibrarycomplexity/main'
include { GATK4_MARKDUPLICATES_SPARK             } from '../../../modules/nf-core/gatk4/markduplicatesspark/main'
include { SAMTOOLS_INDEX as INDEX_MARKDUPLICATES } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_MARKDUPLICATES_SPARK {
    take:
    bam                           // channel: [mandatory] meta, bam
    dict                          // channel: [mandatory] dict
    fasta                         // channel: [mandatory] fasta
    fasta_fai                     // channel: [mandatory] fasta_fai
    intervals_bed_combined        // channel: [optional]  intervals_bed

    main:
    versions = Channel.empty()
    reports = Channel.empty()

    // RUN MARKUPDUPLICATES SPARK
    GATK4_MARKDUPLICATES_SPARK(bam, fasta, fasta_fai, dict)

    // Index cram
    INDEX_MARKDUPLICATES(GATK4_MARKDUPLICATES_SPARK.out.output)

    // Join with the crai file
    cram = GATK4_MARKDUPLICATES_SPARK.out.output.join(INDEX_MARKDUPLICATES.out.crai, failOnDuplicate: true, failOnMismatch: true)

    // QC on CRAM
    CRAM_QC_MOSDEPTH_SAMTOOLS(cram, fasta, intervals_bed_combined)

    // When running Marduplicates spark, and saving reports
    GATK4_ESTIMATELIBRARYCOMPLEXITY(bam, fasta, fasta_fai, dict)

    // Gather all reports generated
    reports = reports.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.metrics)
    reports = reports.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.reports)

    // Gather versions of all tools used
    versions = versions.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.versions)
    versions = versions.mix(GATK4_MARKDUPLICATES_SPARK.out.versions)
    versions = versions.mix(INDEX_MARKDUPLICATES.out.versions)
    versions = versions.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.versions)

    emit:
    cram
    reports

    versions // channel: [ versions.yml ]
}
