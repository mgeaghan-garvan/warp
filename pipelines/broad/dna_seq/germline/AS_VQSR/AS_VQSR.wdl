version 1.0

import "../../../../../tasks/broad/JointGenotypingTasks.wdl" as Tasks

workflow AS_VQSR {
  input {
    Int cohort_size
    File unpadded_intervals_file
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File input_vcf
    String callset_name

    Array[String] snp_recalibration_tranche_values
    Array[String] snp_recalibration_annotation_values
    Array[String] indel_recalibration_tranche_values
    Array[String] indel_recalibration_annotation_values

    Float snp_filter_level
    Float indel_filter_level

    Int small_disk
    Int medium_disk

    File hapmap_resource_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf
    File one_thousand_genomes_resource_vcf_index
    File mills_resource_vcf
    File mills_resource_vcf_index
    File axiomPoly_resource_vcf
    File axiomPoly_resource_vcf_index
    File dbsnp_resource_vcf
    File dbsnp_resource_vcf_index
  }

  Boolean is_small_callset = if cohort_size < 1000 then true else false
  Boolean is_huge_callset = if cohort_size >= 10000 then true else false

  Int snp_max_gaussians = if is_small_callset then 4 else if is_huge_callset then 8 else 6


  call Tasks.SplitIntervalList as SplitIntervalList {
    input:
      interval_list = unpadded_intervals_file,
      scatter_count = 10,
      scatter_mode = "INTERVAL_SUBDIVISION",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      disk_size = small_disk,
      sample_names_unique_done = true
  }

  call Tasks.CompressAndTabix as CompressAndTabix {
    input:
      input_vcf = input_vcf
  }

  call Tasks.IndelsVariantRecalibrator as IndelsVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf = CompressAndTabix.vcfgz,
      sites_only_variant_filtered_vcf_index = CompressAndTabix.tbi,
      recalibration_filename = callset_name + ".indels.recal",
      tranches_filename = callset_name + ".indels.tranches",
      recalibration_tranche_values = indel_recalibration_tranche_values,
      recalibration_annotation_values = indel_recalibration_annotation_values,
      mills_resource_vcf = mills_resource_vcf,
      mills_resource_vcf_index = mills_resource_vcf_index,
      axiomPoly_resource_vcf = axiomPoly_resource_vcf,
      axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
      use_allele_specific_annotations = true,
      disk_size = small_disk
  }


  call Tasks.SNPsVariantRecalibrator as SNPsVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf = CompressAndTabix.vcfgz,
      sites_only_variant_filtered_vcf_index = CompressAndTabix.tbi,
      recalibration_filename = callset_name + ".snps.recal",
      tranches_filename = callset_name + ".snps.tranches",
      recalibration_tranche_values = snp_recalibration_tranche_values,
      recalibration_annotation_values = snp_recalibration_annotation_values,
      hapmap_resource_vcf = hapmap_resource_vcf,
      hapmap_resource_vcf_index = hapmap_resource_vcf_index,
      omni_resource_vcf = omni_resource_vcf,
      omni_resource_vcf_index = omni_resource_vcf_index,
      one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
      one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
      use_allele_specific_annotations = true,
      max_gaussians = snp_max_gaussians,
      disk_size = small_disk
    }


  call Tasks.ApplyRecalibration as ApplyRecalibration {
    input:
      recalibrated_vcf_filename = callset_name + ".filtered.vcf.gz",
      input_vcf = CompressAndTabix.vcfgz,
      input_vcf_index = CompressAndTabix.tbi,
      indels_recalibration = IndelsVariantRecalibrator.recalibration,
      indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
      indels_tranches = IndelsVariantRecalibrator.tranches,
      snps_recalibration = SNPsVariantRecalibrator.recalibration,
      snps_recalibration_index = SNPsVariantRecalibrator.recalibration_index,
      snps_tranches = SNPsVariantRecalibrator.tranches,
      indel_filter_level = indel_filter_level,
      snp_filter_level = snp_filter_level,
      use_allele_specific_annotations = true,
      disk_size = medium_disk
  }

  call Tasks.AnnotateSB as AnnotateSB {
    input:
      input_vcf = ApplyRecalibration.recalibrated_vcf,
      input_vcf_index = ApplyRecalibration.recalibrated_vcf_index,
      annotated_vcf_filename = callset_name + ".filtered.sb.vcf.gz"
  }

  output {
    File as_vqsr_vcf = AnnotateSB.annotated_vcf
    File as_vqsr_vcf_index = AnnotateSB.annotated_vcf_index
  }
}