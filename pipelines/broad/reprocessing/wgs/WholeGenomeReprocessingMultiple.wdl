version 1.0

import "../../../../pipelines/broad/reprocessing/wgs/WholeGenomeReprocessing.wdl" as WholeGenomeReprocessing
import "../../../../structs/dna_seq/DNASeqStructs.wdl"

workflow WholeGenomeReprocessingMultiple {

  String pipeline_version = "1.0.0"

  input {
    File sample_map

    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings
    PapiSettings papi_settings
    QCSettings qc_settings

    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index

    File wgs_coverage_interval_list
  }

  Array[Array[String]] sample_map_lines = read_tsv(sample_map)
  Int num_samples = length(sample_map_lines)

  scatter (idx in range(num_samples)) {
    String sample_name = sample_map_lines[idx][0]
    File input_cram = sample_map_lines[idx][1]
    File cram_ref_fasta = sample_map_lines[idx][2]
    File cram_ref_fasta_index = sample_map_lines[idx][3]

    call WholeGenomeReprocessing.WholeGenomeReprocessing {
      input:
        input_cram = input_cram,
        sample_name = sample_name,
        base_file_name = sample_name,
        final_gvcf_base_name = sample_name,
        unmapped_bam_suffix = ".unmapped.bam",
        cram_ref_fasta = cram_ref_fasta,
        cram_ref_fasta_index = cram_ref_fasta_index,
        references = references,
        scatter_settings = scatter_settings,
        papi_settings = papi_settings,
        qc_settings = qc_settings,
        fingerprint_genotypes_file = fingerprint_genotypes_file,
        fingerprint_genotypes_index = fingerprint_genotypes_index,
        wgs_coverage_interval_list = wgs_coverage_interval_list
    }

  }

  output {

    Array[File?] read_group_alignment_summary_metrics = select_all(WholeGenomeReprocessing.read_group_alignment_summary_metrics)
    Array[File?] read_group_gc_bias_detail_metrics = select_all(WholeGenomeReprocessing.read_group_gc_bias_detail_metrics)
    Array[File?] read_group_gc_bias_pdf = select_all(WholeGenomeReprocessing.read_group_gc_bias_pdf)
    Array[File?] read_group_gc_bias_summary_metrics = select_all(WholeGenomeReprocessing.read_group_gc_bias_summary_metrics)

    Array[File?] cross_check_fingerprints_metrics = select_all(WholeGenomeReprocessing.cross_check_fingerprints_metrics)

    Array[File] selfSM = select_all(WholeGenomeReprocessing.selfSM)
    Array[Float] contamination = select_all(WholeGenomeReprocessing.contamination)

    Array[File?] calculate_read_group_checksum_md5 = select_all(WholeGenomeReprocessing.calculate_read_group_checksum_md5)

    Array[File?] agg_alignment_summary_metrics = select_all(WholeGenomeReprocessing.agg_alignment_summary_metrics)
    Array[File?] agg_bait_bias_detail_metrics = select_all(WholeGenomeReprocessing.agg_bait_bias_detail_metrics)
    Array[File?] agg_bait_bias_summary_metrics = select_all(WholeGenomeReprocessing.agg_bait_bias_summary_metrics)
    Array[File?] agg_gc_bias_detail_metrics = select_all(WholeGenomeReprocessing.agg_gc_bias_detail_metrics)
    Array[File?] agg_gc_bias_pdf = select_all(WholeGenomeReprocessing.agg_gc_bias_pdf)
    Array[File?] agg_gc_bias_summary_metrics = select_all(WholeGenomeReprocessing.agg_gc_bias_summary_metrics)
    Array[File?] agg_insert_size_histogram_pdf = select_all(WholeGenomeReprocessing.agg_insert_size_histogram_pdf)
    Array[File?] agg_insert_size_metrics = select_all(WholeGenomeReprocessing.agg_insert_size_metrics)
    Array[File?] agg_pre_adapter_detail_metrics = select_all(WholeGenomeReprocessing.agg_pre_adapter_detail_metrics)
    Array[File?] agg_pre_adapter_summary_metrics = select_all(WholeGenomeReprocessing.agg_pre_adapter_summary_metrics)
    Array[File?] agg_quality_distribution_pdf = select_all(WholeGenomeReprocessing.agg_quality_distribution_pdf)
    Array[File?] agg_quality_distribution_metrics = select_all(WholeGenomeReprocessing.agg_quality_distribution_metrics)

    Array[File?] fingerprint_summary_metrics = select_all(WholeGenomeReprocessing.fingerprint_summary_metrics)
    Array[File?] fingerprint_detail_metrics = select_all(WholeGenomeReprocessing.fingerprint_detail_metrics)

    Array[File?] wgs_metrics = select_all(WholeGenomeReprocessing.wgs_metrics)
    Array[File?] raw_wgs_metrics = select_all(WholeGenomeReprocessing.raw_wgs_metrics)

    Array[File] duplicate_metrics = select_all(WholeGenomeReprocessing.duplicate_metrics)
    Array[File?] output_bqsr_reports = select_all(WholeGenomeReprocessing.output_bqsr_reports)

    Array[File?] gvcf_summary_metrics = select_all(WholeGenomeReprocessing.gvcf_summary_metrics)
    Array[File?] gvcf_detail_metrics = select_all(WholeGenomeReprocessing.gvcf_detail_metrics)

    Array[File] output_cram = select_all(WholeGenomeReprocessing.output_cram)
    Array[File] output_cram_index = select_all(WholeGenomeReprocessing.output_cram_index)
    Array[File] output_cram_md5 = select_all(WholeGenomeReprocessing.output_cram_md5)

    Array[File?] validate_cram_file_report = select_all(WholeGenomeReprocessing.validate_cram_file_report)

    Array[File] output_vcf = select_all(WholeGenomeReprocessing.output_vcf)
    Array[File] output_vcf_index = select_all(WholeGenomeReprocessing.output_vcf_index)
  }
  meta {
    allowNestedInputs: true
  }
}
