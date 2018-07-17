#!/usr/bin/env nextflow

//This pipeline uses a reference genome (of high quality ideally)
//to create a repeat library which it then uses to mask
//repeats in the genomes we feed it through the input directory

//+++++++++++++++++ SETUP++++++++++++++++++++++++
params.workdir = "/home/johannes/rdrive/PPG_SEQ_DATA-LICHTJ-SE00182/johannes/notebook/bcinerea_test"

//Reference genome fasta file
params.reference = "${params.workdir}/reference/*.fasta"

//Folder containing fasta files of genomes to be masked
params.input = "${params.workdir}/input/*.fasta"

params.outdir = "${params.workdir}/output"

//Minimum repeat length
params.minRepLength = 250
//+++++++++++++++++++++++++++++++++++++++++++++++

reference = file(params.reference)

//Creating a channel containing all the sequences we want masked
sequences = Channel
.fromPath(params.input)
.map{[it.getBaseName(), it]}
.set{genomes}

log.info "====================================================================="
log.info "Repeat annotation using a reference genome"
log.info "Reference : ${reference}"
log.info "Output    : ${params.outdir}"
log.info "====================================================================="


//Creating a repeat library from the reference genome
process RepeatModeler {
  container 'robsyme/nf-repeatmasking'
  publishDir "${params.outdir}/output/repeatlibrary", mode: 'copy'

  cpus 8

  input:
  file 'genome.fasta' from reference

  output:
  file 'RM*/consensi.fa.classified' into repeatlibrary

  """
BuildDatabase \
 -name reference \
 genome.fasta

RepeatModeler \
 -database reference \
 -pa ${task.cpus}
  """
}

//Masking repeats using the library created above
process RepeatMasker {
  container 'robsyme/nf-repeatmasking'
  tag {sampleID}
  cpus 8

  input:
  set sampleID, "genome.fasta", "library.fasta" from genomes.combine(repeatlibrary)

  output:
  set sampleID, '*.fasta.out', 'genome.fasta' into repeatmaskerout

  """
RepeatMasker \
 -no_is \
 -nolow \
 -gff \
 -pa ${task.cpus} \
 -lib library.fasta \
 genome.fasta
  """
}

process removeShortMatches {
  tag {sampleID}
  container 'robsyme/nf-repeatmasking'
  publishDir "${params.outdir}/${sampleID}", mode: 'copy'
  
  input:
  set sampleID, "rm.out", 'genome.fasta' from repeatmaskerout
  
  output:
  set sampleID, "${sampleID}.out", "${sampleID}*fasta" into repeatMaskerKnowns


  """
head -n 3 rm.out > ${sampleID}.out
tail -n +4 rm.out | awk '\$7 - \$6 > ${params.minRepLength}' >> ${sampleID}.out
tail -n +4 rm.out | awk 'BEGIN{OFS="\\t"} \$7 - \$6 > ${params.minRepLength} {print \$5, \$6, \$7}' >> tmp.bed
maskFastaFromBed -fi genome.fasta -bed tmp.bed -fo ${sampleID}.masked.soft.fasta -soft
maskFastaFromBed -fi genome.fasta -bed tmp.bed -fo ${sampleID}.masked.hard.fasta
  """
}