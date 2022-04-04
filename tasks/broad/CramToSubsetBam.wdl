version 1.0

import "../../tasks/broad/Utilities.wdl" as Utils

workflow CramToSubsetBam {
  input {
    File input_cram
    File input_cram_index
    File ref_fasta
    File ref_fasta_index
    String mt_string

    String base_file_name
    Int agg_preemptible_tries
    String zones
  }

  call Utils.CramToSubsetBam as CramToSubsetBam {
    input:
      input_cram = input_cram,
      input_cram_index = input_cram_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      mt_string = mt_string,
      output_prefix = base_file_name,
      preemptible_tries = agg_preemptible_tries,
      zones = zones
  }

  output {
    File mt_bam = CramToSubsetBam.mt_bam
    File mt_bam_index = CramToSubsetBam.mt_bam_index
  }
  meta {
    allowNestedInputs: true
  }
}