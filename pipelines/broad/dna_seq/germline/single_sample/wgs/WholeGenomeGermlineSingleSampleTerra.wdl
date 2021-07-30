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

import "WholeGenomeGermlineSingleSample.wdl" as WholeGenomeGermlineSingleSample
import "../../../../../../structs/dna_seq/DNASeqStructs.wdl"

# WORKFLOW DEFINITION
workflow WholeGenomeGermlineSingleSampleTerra {

  String pipeline_version = "3.0.0"

  input {
    String sample_and_unmapped_bams_sample_name
    String sample_and_unmapped_bams_base_file_name
    File sample_and_unmapped_bams_flowcell_unmapped_bam
    String sample_and_unmapped_bams_final_gvcf_base_name
    String sample_and_unmapped_bams_unmapped_bam_suffix

    File reference_contamination_sites_ud
    File reference_contamination_sites_bed
    File reference_contamination_sites_mu
    File reference_calling_interval_list
    File reference_reference_fasta_ref_dict
    File reference_reference_fasta_ref_fasta
    File reference_reference_fasta_ref_fasta_index
    File reference_reference_fasta_ref_alt
    File reference_reference_fasta_ref_sa
    File reference_reference_fasta_ref_amb
    File reference_reference_fasta_ref_bwt
    File reference_reference_fasta_ref_ann
    File reference_reference_fasta_ref_pac
    File? reference_reference_fasta_ref_str
    Array[File] reference_known_indels_sites_vcfs
    Array[File] reference_known_indels_sites_indices
    File reference_dbsnp_vcf
    File reference_dbsnp_vcf_index
    File reference_evaluation_interval_list
    File reference_haplotype_database_file

    File? dragmap_reference_reference_bin
    File? dragmap_reference_reference_index_bin
    File? dragmap_reference_hash_table_cmp

    Int scatter_settings_haplotype_scatter_count
    Int scatter_settings_break_bands_at_multiples_of

    Int papi_settings_preemptible_tries
    Int papi_settings_agg_preemptible_tries

    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index

    File wgs_coverage_interval_list

    Boolean provide_bam_output = false
    Boolean use_gatk3_haplotype_caller = false
    Boolean run_dragen_mode = true
    Boolean use_spanning_event_genotyping = true
    Boolean perform_bqsr = true
    Boolean use_bwa_mem = true
  }

  SampleAndUnmappedBams sample_and_unmapped_bams = object {
    sample_name: sample_and_unmapped_bams_sample_name,
    base_file_name: sample_and_unmapped_bams_base_file_name,
    flowcell_unmapped_bams: [sample_and_unmapped_bams_flowcell_unmapped_bam],
    final_gvcf_base_name: sample_and_unmapped_bams_final_gvcf_base_name,
    unmapped_bam_suffix: sample_and_unmapped_bams_unmapped_bam_suffix
  }

  VariantCallingScatterSettings scatter_settings = object {
    haplotype_scatter_count: scatter_settings_haplotype_scatter_count,
    break_bands_at_multiples_of: scatter_settings_break_bands_at_multiples_of
  }

  DNASeqSingleSampleReferences references = object {
    contamination_sites_ud: reference_contamination_sites_ud,
    contamination_sites_bed: reference_contamination_sites_bed,
    contamination_sites_mu: reference_contamination_sites_mu,
    calling_interval_list: reference_calling_interval_list,
    reference_fasta: object {
      ref_dict: reference_reference_fasta_ref_dict,
      ref_fasta: reference_reference_fasta_ref_fasta,
      ref_fasta_index: reference_reference_fasta_ref_fasta_index,
      ref_alt: reference_reference_fasta_ref_alt,
      ref_sa: reference_reference_fasta_ref_sa,
      ref_amb: reference_reference_fasta_ref_amb,
      ref_bwt: reference_reference_fasta_ref_bwt,
      ref_ann: reference_reference_fasta_ref_ann,
      ref_pac: reference_reference_fasta_ref_pac,
      ref_str: reference_reference_fasta_ref_str
    },
    known_indels_sites_vcfs: reference_known_indels_sites_vcfs,
    known_indels_sites_indices: reference_known_indels_sites_indices,
    dbsnp_vcf: reference_dbsnp_vcf,
    dbsnp_vcf_index: reference_dbsnp_vcf_index,
    evaluation_interval_list: reference_evaluation_interval_list,
    haplotype_database_file: reference_haplotype_database_file
  }

  if (defined(dragmap_reference_reference_bin) && defined(dragmap_reference_reference_index_bin) && defined(dragmap_reference_hash_table_cmp)) {
    DragmapReference dragmap_reference = object {
      reference_bin: dragmap_reference_reference_bin,
      reference_index_bin: dragmap_reference_reference_index_bin,
      hash_table_cmp: dragmap_reference_hash_table_cmp
    }
  }

  PapiSettings papi_settings = object {
    preemptible_tries: papi_settings_preemptible_tries,
    agg_preemptible_tries: papi_settings_agg_preemptible_tries
  }
  
  call WholeGenomeGermlineSingleSample.WholeGenomeGermlineSingleSample as WholeGenomeGermlineSingleSampleTask {
    input:
      sample_and_unmapped_bams = sample_and_unmapped_bams,
      references = references,
      dragmap_reference = dragmap_reference,
      scatter_settings = scatter_settings,
      papi_settings = papi_settings,

      fingerprint_genotypes_file = fingerprint_genotypes_file,
      fingerprint_genotypes_index = fingerprint_genotypes_index,

      wgs_coverage_interval_list = wgs_coverage_interval_list,

      provide_bam_output = provide_bam_output,
      use_gatk3_haplotype_caller = use_gatk3_haplotype_caller,
      run_dragen_mode = run_dragen_mode,
      use_spanning_event_genotyping = use_spanning_event_genotyping,
      perform_bqsr = perform_bqsr,
      use_bwa_mem = use_bwa_mem
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] quality_yield_metrics = WholeGenomeGermlineSingleSampleTask.quality_yield_metrics

    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = WholeGenomeGermlineSingleSampleTask.unsorted_read_group_base_distribution_by_cycle_pdf
    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = WholeGenomeGermlineSingleSampleTask.unsorted_read_group_base_distribution_by_cycle_metrics
    Array[File] unsorted_read_group_insert_size_histogram_pdf = WholeGenomeGermlineSingleSampleTask.unsorted_read_group_insert_size_histogram_pdf
    Array[File] unsorted_read_group_insert_size_metrics = WholeGenomeGermlineSingleSampleTask.unsorted_read_group_insert_size_metrics
    Array[File] unsorted_read_group_quality_by_cycle_pdf = WholeGenomeGermlineSingleSampleTask.unsorted_read_group_quality_by_cycle_pdf
    Array[File] unsorted_read_group_quality_by_cycle_metrics = WholeGenomeGermlineSingleSampleTask.unsorted_read_group_quality_by_cycle_metrics
    Array[File] unsorted_read_group_quality_distribution_pdf = WholeGenomeGermlineSingleSampleTask.unsorted_read_group_quality_distribution_pdf
    Array[File] unsorted_read_group_quality_distribution_metrics = WholeGenomeGermlineSingleSampleTask.unsorted_read_group_quality_distribution_metrics

    File read_group_alignment_summary_metrics = WholeGenomeGermlineSingleSampleTask.read_group_alignment_summary_metrics
    File read_group_gc_bias_detail_metrics = WholeGenomeGermlineSingleSampleTask.read_group_gc_bias_detail_metrics
    File read_group_gc_bias_pdf = WholeGenomeGermlineSingleSampleTask.read_group_gc_bias_pdf
    File read_group_gc_bias_summary_metrics = WholeGenomeGermlineSingleSampleTask.read_group_gc_bias_summary_metrics

    File? cross_check_fingerprints_metrics = WholeGenomeGermlineSingleSampleTask.cross_check_fingerprints_metrics

    File selfSM = WholeGenomeGermlineSingleSampleTask.selfSM
    Float contamination = WholeGenomeGermlineSingleSampleTask.contamination

    File calculate_read_group_checksum_md5 = WholeGenomeGermlineSingleSampleTask.calculate_read_group_checksum_md5

    File agg_alignment_summary_metrics = WholeGenomeGermlineSingleSampleTask.agg_alignment_summary_metrics
    File agg_bait_bias_detail_metrics = WholeGenomeGermlineSingleSampleTask.agg_bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = WholeGenomeGermlineSingleSampleTask.agg_bait_bias_summary_metrics
    File agg_gc_bias_detail_metrics = WholeGenomeGermlineSingleSampleTask.agg_gc_bias_detail_metrics
    File agg_gc_bias_pdf = WholeGenomeGermlineSingleSampleTask.agg_gc_bias_pdf
    File agg_gc_bias_summary_metrics = WholeGenomeGermlineSingleSampleTask.agg_gc_bias_summary_metrics
    File agg_insert_size_histogram_pdf = WholeGenomeGermlineSingleSampleTask.agg_insert_size_histogram_pdf
    File agg_insert_size_metrics = WholeGenomeGermlineSingleSampleTask.agg_insert_size_metrics
    File agg_pre_adapter_detail_metrics = WholeGenomeGermlineSingleSampleTask.agg_pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = WholeGenomeGermlineSingleSampleTask.agg_pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = WholeGenomeGermlineSingleSampleTask.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = WholeGenomeGermlineSingleSampleTask.agg_quality_distribution_metrics
    File agg_error_summary_metrics = WholeGenomeGermlineSingleSampleTask.agg_error_summary_metrics

    File? fingerprint_summary_metrics = WholeGenomeGermlineSingleSampleTask.fingerprint_summary_metrics
    File? fingerprint_detail_metrics = WholeGenomeGermlineSingleSampleTask.fingerprint_detail_metrics

    File wgs_metrics = WholeGenomeGermlineSingleSampleTask.wgs_metrics
    File raw_wgs_metrics = WholeGenomeGermlineSingleSampleTask.raw_wgs_metrics

    File duplicate_metrics = WholeGenomeGermlineSingleSampleTask.duplicate_metrics
    File? output_bqsr_reports = WholeGenomeGermlineSingleSampleTask.output_bqsr_reports

    File gvcf_summary_metrics = WholeGenomeGermlineSingleSampleTask.gvcf_summary_metrics
    File gvcf_detail_metrics = WholeGenomeGermlineSingleSampleTask.gvcf_detail_metrics

    File? output_bam = WholeGenomeGermlineSingleSampleTask.output_bam
    File? output_bam_index = WholeGenomeGermlineSingleSampleTask.output_bam_index

    File output_cram = WholeGenomeGermlineSingleSampleTask.output_cram
    File output_cram_index = WholeGenomeGermlineSingleSampleTask.output_cram_index
    File output_cram_md5 = WholeGenomeGermlineSingleSampleTask.output_cram_md5

    File validate_cram_file_report = WholeGenomeGermlineSingleSampleTask.validate_cram_file_report

    File output_vcf = WholeGenomeGermlineSingleSampleTask.output_vcf
    File output_vcf_index = WholeGenomeGermlineSingleSampleTask.output_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
}
