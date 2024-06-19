#!/usr/bin/env nextflow

#Usage: nextflow run main.nf --input_dir input_data --metadata_file metadata_file.tsv --classifier classifier_file.qza --results_output_dir results-output

nextflow.enable.dsl=2

// Create the results output directory before the workflow starts
def createResultsOutputDir() {
    def results_output_dir = file(params.results_output_dir).toString()
    def dir = new File(results_output_dir)
    if (!dir.exists()) {
        dir.mkdirs()
    }
    return results_output_dir
}

params.results_output_dir = createResultsOutputDir()

process import_data {
    input:
    path fastq_dir

    output:
    path("q1-importdemux/q1-demux-pe.qza"), emit: demux_qza
    path("q1-importdemux/q1-demux-pe.qzv"), emit: demux_qzv

    publishDir "${params.results_output_dir}", mode: 'copy'

    script:
    """
    mkdir -p q1-importdemux

    qiime tools import \\
        --type 'SampleData[PairedEndSequencesWithQuality]' \\
        --input-path ${fastq_dir} \\
        --input-format CasavaOneEightSingleLanePerSampleDirFmt \\
        --output-path q1-importdemux/q1-demux-pe.qza

    qiime demux summarize \\
        --i-data q1-importdemux/q1-demux-pe.qza \\
        --o-visualization q1-importdemux/q1-demux-pe.qzv
    """
}

process denoise {
    input:
    path demux_file
    path metadata_file

    output:
    path("q2-denois-ftb/q2-rep-seqs-dada2.qza"), emit: rep_seqs_qza
    path("q2-denois-ftb/q2-table-dada2.qza"), emit: table_qza
    path("q2-denois-ftb/q2-stats-dada2.qza"), emit: stats_qza
    path("q2-denois-ftb/q2-stats-dada2.qzv"), emit: stats_qzv
    path("q2-denois-ftb/q2-table-dada2.qzv"), emit: table_qzv
    path("q2-denois-ftb/q2-rep-seqs-dada2.qzv"), emit: rep_seqs_qzv

    publishDir "${params.results_output_dir}", mode: 'copy'

    script:
    """
    mkdir -p q2-denois-ftb

    qiime dada2 denoise-paired \\
        --i-demultiplexed-seqs ${demux_file} \\
        --p-trim-left-f 0 \\
        --p-trunc-len-f 0 \\
        --p-trim-left-r 0 \\
        --p-trunc-len-r 0 \\
        --o-representative-sequences q2-denois-ftb/q2-rep-seqs-dada2.qza \\
        --o-table q2-denois-ftb/q2-table-dada2.qza \\
        --o-denoising-stats q2-denois-ftb/q2-stats-dada2.qza

    qiime metadata tabulate \\
        --m-input-file q2-denois-ftb/q2-stats-dada2.qza \\
        --o-visualization q2-denois-ftb/q2-stats-dada2.qzv

    qiime feature-table summarize \\
        --i-table q2-denois-ftb/q2-table-dada2.qza \\
        --o-visualization q2-denois-ftb/q2-table-dada2.qzv

    qiime feature-table tabulate-seqs \\
        --i-data q2-denois-ftb/q2-rep-seqs-dada2.qza \\
        --o-visualization q2-denois-ftb/q2-rep-seqs-dada2.qzv
    """
}

process taxonomy {
    input:
    path rep_seqs_file
    path table_file
    path classifier
    path metadata_file

    output:
    path("q3-taxoasn/q3-taxoasn-table-s138.qza"), emit: taxoasn_qza
    path("q3-taxoasn/q3-taxoasn-table-s138.qzv"), emit: taxoasn_qzv
    path("q3-taxoasn/q3-taxoasn-table-s138-barplots.qzv"), emit: taxoasn_barplots_qzv

    publishDir "${params.results_output_dir}", mode: 'copy'

    script:
    """
    mkdir -p q3-taxoasn

    qiime feature-classifier classify-sklearn \\
        --i-classifier ${classifier} \\
        --i-reads ${rep_seqs_file} \\
        --o-classification q3-taxoasn/q3-taxoasn-table-s138.qza

    qiime metadata tabulate \\
        --m-input-file q3-taxoasn/q3-taxoasn-table-s138.qza \\
        --o-visualization q3-taxoasn/q3-taxoasn-table-s138.qzv

    qiime taxa barplot \\
        --i-table ${table_file} \\
        --i-taxonomy q3-taxoasn/q3-taxoasn-table-s138.qza \\
        --m-metadata-file ${metadata_file} \\
        --o-visualization q3-taxoasn/q3-taxoasn-table-s138-barplots.qzv
    """
}

workflow {
    // Create channels
    input_channel = Channel.fromPath(params.input_dir)
    metadata_channel = Channel.fromPath(params.metadata_file)
    classifier_channel = Channel.fromPath(params.classifier)
    
    // Run the processes and manage dependencies
    import_output = import_data(input_channel)
    denoise_output = denoise(import_output.demux_qza, metadata_channel)
    taxonomy_output = taxonomy(denoise_output.rep_seqs_qza, denoise_output.table_qza, classifier_channel, metadata_channel)
}
