version 1.0

import "../../tasks/broad/JointGenotypingTasks.wdl" as JointGenotypingTasks

workflow ReblockGVCFs {
  input {
    File input_gvcf
    File input_tbi
    File ref_fasta
    File ref_fastq_index
    File ref_dict
    String output_prefix
  }

  call JointGenotypingTasks.ReblockGVCFs as ReblockGVCFs {
    input:
      input_gvcf = input_gvcf,
      input_tbi = input_tbi,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fastq_index = ref_fastq_index,
      output_prefix = output_prefix
  }

  meta {
    allowNestedInputs: true
  }

  output {
    File reblock_gvcf = ReblockGVCFs.reblock_gvcf
    File reblock_tbi = ReblockGVCFs.reblock_tbi
  }
}
