nextflow.enable.dsl=2

params.str = 'Hello world!'
params.fastq = "/mnt/disks/sec/David/10x/fastq"
params.in_fn = "/mnt/disks/sec/David/10x/prefixes.txt"
params.out = "/mnt/disks/sec/David/10x/nextflow_output"
params.ncells=10000
params.dev = false
params.number_of_inputs = 2
params.umi_len = 12

def osn_conda_env = file(params.home + "/miniconda3/envs/$params.osn_env")
def vel_conda_env = file(params.home + "/miniconda3/envs/$params.vel_env")

log.info("Using str: ${params.str}")
log.info("Using home: ${params.home}")
log.info("found output directory param: ${params.out}")
log.info("Using conda_paths: $osn_conda_env and $vel_conda_env")


allLines = file(params.in_fn).readLines()
println(allLines)

Channel
    .fromList(allLines)
    .take( params.dev ? params.number_of_inputs : -1 )
    .set { prefixes }

log.info("Found these prefixes....")
prefixes.view()

process cellranger_count {
    label "big"

    publishDir path: "${params.out}/${id}/cell_ranger", mode: "copy"

    input:
        val(id)

    output:
        tuple val(id), path("${id}/outs"), emit: out_folder
        path("${id}/_*")
        path("${id}/*.tgz")

    // also specify the sample?
    script:
    """
    echo $id
    ls ${params.fastq} | grep ${id}
    cellranger count --id=${id} --fastqs=${params.fastq} --sample=${id} \
    --transcriptome=${params.cell_reference} --expect-cells=${params.ncells}
    """
}


process sort_cells {
    label "small"
    conda osn_conda_env

    input:
        tuple val(id), path(cell)

    output:
        tuple val(id), path('cell_sorted_possorted_genome_bam.unique.bam'), path('barcodes.tsv')

    script:
    """
    echo "id: $id, cell: $cell"
    which sort_cells
    sort_cells $cell/possorted_genome_bam.bam $cell/filtered_feature_bc_matrix/barcodes.tsv.gz
    """
}


process velocyto {
    label "small"
    conda vel_conda_env
    publishDir path: "${params.out}/${id}", mode: "copy"

    input:
        tuple val(id), path(cell)

    output:
        tuple val(id), path('velocyto/*.loom')

    script:
    """
    echo "id: $id, cell: $cell"
    which velocyto
    # subset bam file
    zcat $cell/filtered_feature_bc_matrix/barcodes.tsv.gz > barcodes.tsv
    samtools view -b -F 0x4 -e '([NH] && [NH]==1) && [UB]' \
    -D CB:barcodes.tsv \
    -o possorted_genome_bam.unique.bam -@ ${task.cpus} $cell/possorted_genome_bam.bam
    samtools index -@ ${task.cpus} possorted_genome_bam.unique.bam
    velocyto run \
    -b barcodes.tsv \
    -e $id \
    -m $params.vel_mask \
    --samtools-threads 15 \
    possorted_genome_bam.unique.bam $params.gtf
    """
}


process make_adata {
    label "small"
    conda osn_conda_env
    publishDir path: "${params.out}/$id", mode: "copy", pattern: "${id}_counts.*"

    input:
        tuple val(id), path(cell_bam), path(barcodes)

    output:
        tuple val(id), path("${id}_counts.*")

    script:
    """
    echo $id
    which filter_bam
    filter_bam -i $cell_bam -b $barcodes -o "${id}_counts.tsv"
    make_adata -f "${id}_counts.tsv" -o "${id}_counts.h5ad" -p ${id}
    pigz "${id}_counts.tsv"
    """
}


workflow {
    cellranger_count(prefixes)
    cellranger_count.out.out_folder.view()
    velocyto(cellranger_count.out.out_folder)
    velocyto.out.view()
    sort_cells(cellranger_count.out.out_folder)
    sort_cells.out.view()
    make_adata(sort_cells.out)
    make_adata.out.view()
}

workflow.onComplete {
log.info("##Pipeline execution summary##")
log.info("---------------------------")
log.info("##Completed at: $workflow.complete")
log.info("##Duration: ${workflow.duration}")
log.info("##Success: ${workflow.success ? 'OK' : 'failed' }")
log.info("##Exit status: ${workflow.exitStatus}")
}