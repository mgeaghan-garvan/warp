version 1.0

import "../../../../pipelines/broad/reprocessing/wgs/WholeGenomeFromFastq.wdl" as WholeGenomeFromFastq
import "../../../../structs/dna_seq/DNASeqStructs.wdl"

workflow WholeGenomeFromMultipleFastq {

  String pipeline_version = "1.0.0"

  input {
    File sample_map

    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings
    PapiSettings papi_settings
    QCSettings qc_settings

    File wgs_coverage_interval_list

    String zones
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

    call WholeGenomeFromFastq.WholeGenomeFromFastq {
      input:
        input_R1 = R1,
        input_R2 = R2,
        RGID = RGID,
        RGPL = RGPL,
        RGLB = RGLB,
        RGPU = RGPU,
        RGCN = RGCN,
        sample_name = sample_name,
        base_file_name = sample_name,
        final_gvcf_base_name = sample_name,
        unmapped_bam_suffix = ".unmapped.bam",

        references = references,
        scatter_settings = scatter_settings,
        papi_settings = papi_settings,
        qc_settings = qc_settings,

        wgs_coverage_interval_list = wgs_coverage_interval_list,
        zones = zones
    }
  }

  output {
    Array[File?] read_group_alignment_summary_metrics = select_all(WholeGenomeFromFastq.read_group_alignment_summary_metrics)
    Array[File?] read_group_gc_bias_detail_metrics = select_all(WholeGenomeFromFastq.read_group_gc_bias_detail_metrics)
    Array[File?] read_group_gc_bias_pdf = select_all(WholeGenomeFromFastq.read_group_gc_bias_pdf)
    Array[File?] read_group_gc_bias_summary_metrics = select_all(WholeGenomeFromFastq.read_group_gc_bias_summary_metrics)

    Array[File?] cross_check_fingerprints_metrics = select_all(WholeGenomeFromFastq.cross_check_fingerprints_metrics)

    Array[File] selfSM = select_all(WholeGenomeFromFastq.selfSM)
    Array[Float] contamination = select_all(WholeGenomeFromFastq.contamination)

    Array[File?] calculate_read_group_checksum_md5 = select_all(WholeGenomeFromFastq.calculate_read_group_checksum_md5)

    Array[File?] agg_alignment_summary_metrics = select_all(WholeGenomeFromFastq.agg_alignment_summary_metrics)
    Array[File?] agg_bait_bias_detail_metrics = select_all(WholeGenomeFromFastq.agg_bait_bias_detail_metrics)
    Array[File?] agg_bait_bias_summary_metrics = select_all(WholeGenomeFromFastq.agg_bait_bias_summary_metrics)
    Array[File?] agg_gc_bias_detail_metrics = select_all(WholeGenomeFromFastq.agg_gc_bias_detail_metrics)
    Array[File?] agg_gc_bias_pdf = select_all(WholeGenomeFromFastq.agg_gc_bias_pdf)
    Array[File?] agg_gc_bias_summary_metrics = select_all(WholeGenomeFromFastq.agg_gc_bias_summary_metrics)
    Array[File?] agg_insert_size_histogram_pdf = select_all(WholeGenomeFromFastq.agg_insert_size_histogram_pdf)
    Array[File?] agg_insert_size_metrics = select_all(WholeGenomeFromFastq.agg_insert_size_metrics)
    Array[File?] agg_pre_adapter_detail_metrics = select_all(WholeGenomeFromFastq.agg_pre_adapter_detail_metrics)
    Array[File?] agg_pre_adapter_summary_metrics = select_all(WholeGenomeFromFastq.agg_pre_adapter_summary_metrics)
    Array[File?] agg_quality_distribution_pdf = select_all(WholeGenomeFromFastq.agg_quality_distribution_pdf)
    Array[File?] agg_quality_distribution_metrics = select_all(WholeGenomeFromFastq.agg_quality_distribution_metrics)

    Array[File?] fingerprint_summary_metrics = select_all(WholeGenomeFromFastq.fingerprint_summary_metrics)
    Array[File?] fingerprint_detail_metrics = select_all(WholeGenomeFromFastq.fingerprint_detail_metrics)

    Array[File?] wgs_metrics = select_all(WholeGenomeFromFastq.wgs_metrics)
    Array[File?] raw_wgs_metrics = select_all(WholeGenomeFromFastq.raw_wgs_metrics)

    Array[File] duplicate_metrics = select_all(WholeGenomeFromFastq.duplicate_metrics)
    Array[File?] output_bqsr_reports = select_all(WholeGenomeFromFastq.output_bqsr_reports)

    Array[File?] gvcf_summary_metrics = select_all(WholeGenomeFromFastq.gvcf_summary_metrics)
    Array[File?] gvcf_detail_metrics = select_all(WholeGenomeFromFastq.gvcf_detail_metrics)

    Array[File] output_cram = select_all(WholeGenomeFromFastq.output_cram)
    Array[File] output_cram_index = select_all(WholeGenomeFromFastq.output_cram_index)
    Array[File] output_cram_md5 = select_all(WholeGenomeFromFastq.output_cram_md5)

    Array[File?] validate_cram_file_report = select_all(WholeGenomeFromFastq.validate_cram_file_report)

    Array[File] output_vcf = select_all(WholeGenomeFromFastq.output_vcf)
    Array[File] output_vcf_index = select_all(WholeGenomeFromFastq.output_vcf_index)
  }
  meta {
    allowNestedInputs: true
  }
}
