version 1.0

## Copyright Broad Institute, 2018
##
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF
## generation) according to the GATK Best Practices (June 2016) for germline SNP and
## Indel discovery in human whole-genome data.
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
## - GVCF output names must end in ".g.vcf.gz"
## - Reference genome must be Hg38 with ALT contigs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

import "../../../../../../pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSampleFastq.wdl" as WholeGenomeGermlineSingleSampleFastq
import "../../../../../../structs/dna_seq/DNASeqStructs.wdl"

# WORKFLOW DEFINITION
workflow WholeGenomeGermlineMultipleSamplesFastq {

  String pipeline_version = "1.0.0"

  input {
    File sample_map
    DNASeqSingleSampleReferences references
    DragmapReference? dragmap_reference
    VariantCallingScatterSettings scatter_settings
    PapiSettings papi_settings
    QCSettings qc_settings

    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index

    File wgs_coverage_interval_list

    String zones

    Boolean provide_bam_output = false
    Boolean use_gatk3_haplotype_caller = true

    Boolean dragen_functional_equivalence_mode = false
    Boolean dragen_maximum_quality_mode = false

    Boolean run_dragen_mode_variant_calling = false
    Boolean use_spanning_event_genotyping = true
    Boolean unmap_contaminant_reads = true
    Boolean perform_bqsr = true
    Boolean use_bwa_mem = true
    Boolean use_dragen_hard_filtering = false
  }

  Array[Array[String]] sample_map_lines = read_tsv(sample_map)
  Int num_samples = length(sample_map_lines)

  scatter (idx in range(num_samples)) {
    String sample_name = sample_map_lines[idx][0]
    File R1 = sample_map_lines[idx][1]
    File R2 = sample_map_lines[idx][2]
    String RGID = select_first([sample_map_lines[idx][3], "None"])
    String RGPL = select_first([sample_map_lines[idx][4], "ILLUMINA"])
    String RGPU = select_first([sample_map_lines[idx][5], "None"])
    String RGLB = select_first([sample_map_lines[idx][6], "None"])
    String RGCN = select_first([sample_map_lines[idx][7], "None"])

    SampleInput sample_and_fastqs = object {
      sample_name: sample_name,
      base_file_name: sample_name,
      input_R1: R1,
      input_R2: R2,
      final_gvcf_base_name: sample_name,
      RGID: RGID,
      RGPL: RGPL,
      RGPU: RGPU,
      RGLB: RGLB,
      RGCN: RGCN
    }

    call WholeGenomeGermlineSingleSampleFastq.WholeGenomeGermlineSingleSampleFastq {
      input:
        sample_and_fastqs = sample_and_fastqs,
        references = references,
        dragmap_reference = dragmap_reference,
        scatter_settings = scatter_settings,
        papi_settings = papi_settings,
        qc_settings = qc_settings,

        fingerprint_genotypes_file = fingerprint_genotypes_file,
        fingerprint_genotypes_index = fingerprint_genotypes_index,

        wgs_coverage_interval_list = wgs_coverage_interval_list,

        provide_bam_output = provide_bam_output,
        use_gatk3_haplotype_caller = use_gatk3_haplotype_caller,

        dragen_functional_equivalence_mode = dragen_functional_equivalence_mode,
        dragen_maximum_quality_mode = dragen_maximum_quality_mode,

        run_dragen_mode_variant_calling = run_dragen_mode_variant_calling,
        use_spanning_event_genotyping = use_spanning_event_genotyping,
        unmap_contaminant_reads = unmap_contaminant_reads,
        perform_bqsr = perform_bqsr,
        use_bwa_mem = use_bwa_mem,
        use_dragen_hard_filtering = use_dragen_hard_filtering,
        zones = zones
    }
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File?] read_group_alignment_summary_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.read_group_alignment_summary_metrics)
    Array[File?] read_group_gc_bias_detail_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.read_group_gc_bias_detail_metrics)
    Array[File?] read_group_gc_bias_pdf = select_all(WholeGenomeGermlineSingleSampleFastq.read_group_gc_bias_pdf)
    Array[File?] read_group_gc_bias_summary_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.read_group_gc_bias_summary_metrics)

    Array[File?] cross_check_fingerprints_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.cross_check_fingerprints_metrics)

    Array[File] selfSM = select_all(WholeGenomeGermlineSingleSampleFastq.selfSM)
    Array[Float] contamination = select_all(WholeGenomeGermlineSingleSampleFastq.contamination)

    Array[File?] calculate_read_group_checksum_md5 = select_all(WholeGenomeGermlineSingleSampleFastq.calculate_read_group_checksum_md5)

    Array[File?] agg_alignment_summary_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.agg_alignment_summary_metrics)
    Array[File?] agg_bait_bias_detail_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.agg_bait_bias_detail_metrics)
    Array[File?] agg_bait_bias_summary_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.agg_bait_bias_summary_metrics)
    Array[File?] agg_gc_bias_detail_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.agg_gc_bias_detail_metrics)
    Array[File?] agg_gc_bias_pdf = select_all(WholeGenomeGermlineSingleSampleFastq.agg_gc_bias_pdf)
    Array[File?] agg_gc_bias_summary_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.agg_gc_bias_summary_metrics)
    Array[File?] agg_insert_size_histogram_pdf = select_all(WholeGenomeGermlineSingleSampleFastq.agg_insert_size_histogram_pdf)
    Array[File?] agg_insert_size_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.agg_insert_size_metrics)
    Array[File?] agg_pre_adapter_detail_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.agg_pre_adapter_detail_metrics)
    Array[File?] agg_pre_adapter_summary_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.agg_pre_adapter_summary_metrics)
    Array[File?] agg_quality_distribution_pdf = select_all(WholeGenomeGermlineSingleSampleFastq.agg_quality_distribution_pdf)
    Array[File?] agg_quality_distribution_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.agg_quality_distribution_metrics)
    Array[File?] agg_error_summary_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.agg_error_summary_metrics)

    Array[File?] fingerprint_summary_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.fingerprint_summary_metrics)
    Array[File?] fingerprint_detail_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.fingerprint_detail_metrics)

    Array[File?] wgs_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.wgs_metrics)
    Array[File?] raw_wgs_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.raw_wgs_metrics)

    Array[File] duplicate_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.duplicate_metrics)
    Array[File?] output_bqsr_reports = select_all(WholeGenomeGermlineSingleSampleFastq.output_bqsr_reports)

    Array[File?] gvcf_summary_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.gvcf_summary_metrics)
    Array[File?] gvcf_detail_metrics = select_all(WholeGenomeGermlineSingleSampleFastq.gvcf_detail_metrics)

    Array[File?] output_bam = select_all(WholeGenomeGermlineSingleSampleFastq.output_bam)
    Array[File?] output_bam_index = select_all(WholeGenomeGermlineSingleSampleFastq.output_bam_index)

    Array[File] output_cram = select_all(WholeGenomeGermlineSingleSampleFastq.output_cram)
    Array[File] output_cram_index = select_all(WholeGenomeGermlineSingleSampleFastq.output_cram_index)
    Array[File] output_cram_md5 = select_all(WholeGenomeGermlineSingleSampleFastq.output_cram_md5)

    Array[File?] validate_cram_file_report = select_all(WholeGenomeGermlineSingleSampleFastq.validate_cram_file_report)

    Array[File] output_vcf = select_all(WholeGenomeGermlineSingleSampleFastq.output_vcf)
    Array[File] output_vcf_index = select_all(WholeGenomeGermlineSingleSampleFastq.output_vcf_index)
  }
  meta {
    allowNestedInputs: true
  }
}
