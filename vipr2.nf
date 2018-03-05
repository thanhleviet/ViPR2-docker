#!/usr/bin/env nextflow


// General configuration variables
params.help = false
params.pwd = "$PWD"
params.input = "$baseDir/" + "data"
params.output = "den-hcm"
params.reads = "*_L001_R{1,2}_001.fastq.gz"
params.read_pairs = params.input + "/" + params.reads
params.out_dir = "$baseDir/" + params.output
params.refdb = "$baseDir/" + "ref"
params.primerdb = "$baseDir/" + "primer"
params.threads = 8
params.trim = true
params.normalize = false
threads = params.threads

// Trimmomatic configuration variables
params.leading = 30
params.trailing = 30
params.slidingwindow = "4:30"
params.minlen = 15
params.adapters = "NexteraPE-PE.fa"

leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen
adapters = params.adapters

LB = "Nextera"

// Reference
ref_db = file(params.refdb)
ref_denv = file(ref_db + "/denv.fasta")
ref_denv1 = "DENV1_FJ687432.1.fasta"
ref_denv2 = "DENV2_FJ639718.1.fasta"
ref_denv3 = "DENV3_AB189127.1.fasta"
ref_denv4 = "DENV4_AY618993.1.fasta"

// Primer
primerdb = file(params.primerdb)
primer_denv1 = "primer_d1.fasta"
primer_denv2 = "primer_d2.fasta"
primer_denv3 = "primer_d3.fasta"
primer_denv4 = "primer_d4.fasta"


Channel
    .fromFilePairs(params.read_pairs, flat: true)
    .ifEmpty { exit 1, "Read pairs could not be found: ${params.read_pairs}" }
    .set { trimmomatic_read_pairs }


if (params.trim){
/*
 * Remove adapter sequences and low quality base pairs with Trimmomatic
 */
process RunQC {
    // publishDir "${params.out_dir}/${dataset_id}/", mode: "copy", pattern: "fastqc"

    tag { dataset_id }

    input:
    set dataset_id, file(forward), file(reverse) from trimmomatic_read_pairs

    output:
    set dataset_id, file("${dataset_id}_R1.fastq.gz"), file("${dataset_id}_R2.fastq.gz") into QC
    file("fastqc") into ch_fastqc

    """
    mkdir fastqc

    trimmomatic PE -threads ${threads} \
    $forward $reverse -baseout ${dataset_id} \
    ILLUMINACLIP:${TRIMMOMATIC}/${adapters}:2:30:10:3:TRUE \
    LEADING:${leading} \
    TRAILING:${trailing} \
    SLIDINGWINDOW:${slidingwindow} \
    MINLEN:${minlen}

    cat ${dataset_id}_1P | gzip -1 > ${dataset_id}_R1.fastq.gz
    cat ${dataset_id}_2P | gzip -1 > ${dataset_id}_R2.fastq.gz
    rm ${dataset_id}_1P
    rm ${dataset_id}_2P

    fastqc $forward $reverse *_R*.fastq.gz -o fastqc
    multiqc -o fastqc fastqc
    """
    }
} else {
    trimmomatic_read_pairs.into {QC}
}

process PreMapping {
    // publishDir "${params.out_dir}/${dataset_id}", mode: "copy", pattern: "*.{bam,bai,fasta,fai,dict,tsv}"

    tag {dataset_id}

    input:
    set dataset_id, file(forward), file(reverse) from QC
    file(ref_db)
    file(primerdb)

    output:
    set dataset_id, file(forward), file(reverse), file("new_ref.fasta"), file("primer.fasta") into ch_pre_check

    shell:
    '''
    bwa mem -t !{threads} !{ref_db}/denv.fasta !{forward} !{reverse} | \
    samtools view -@ !{threads} -Sb - | samtools sort  -@ !{threads} -o tmp.bam -

    samtools index -@ !{threads} tmp.bam

    st=`samtools idxstats tmp.bam | sort -k3 -nr | awk 'NR==1 {print $1}' | tr -d '\n'`

    get_ref.py $st \
    !{ref_db}/!{ref_denv1} \
    !{ref_db}/!{ref_denv2} \
    !{ref_db}/!{ref_denv3} \
    !{ref_db}/!{ref_denv4} \
    !{primerdb}/!{primer_denv1} \
    !{primerdb}/!{primer_denv2} \
    !{primerdb}/!{primer_denv3} \
    !{primerdb}/!{primer_denv4}

    rm tmp.bam
    '''

}

process ViPR2 {
    publishDir "${params.out_dir}/", mode: "copy"

    tag {dataset_id}

    input:
    set dataset_id, file(forward), file(reverse), file("new_ref.fasta"), file("primer.fasta") from ch_pre_check
    file("fastqc") from ch_fastqc
    output:
    file("${dataset_id}")

    """
    vipr2t.py -1 ${forward} -2 ${reverse} -o ${dataset_id} -r new_ref.fasta -p primer.fasta -n ${dataset_id}

    cp new_ref.fasta ${dataset_id}/ref.fa

    mv fastqc ${dataset_id}/fastqc

    mv ${dataset_id}/1_R1.fastq.gz ${dataset_id}/${dataset_id}_R1.fastq.gz
    mv ${dataset_id}/1_R2.fastq.gz ${dataset_id}/${dataset_id}_R2.fastq.gz
    """
}

// Display information about the completed run
// See https://www.nextflow.io/docs/latest/metadata.html for more
// information about available onComplete options
workflow.onComplete {
    log.info "Nextflow Version: $workflow.nextflow.version"
    log.info "Command Line:     $workflow.commandLine"
    log.info "Container:        $workflow.container"
    log.info "Duration:     $workflow.duration"
    log.info "Output Directory: $params.out_dir"
}
