#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.cpus = 63
params.outdir = './results/'
params.tmpdir = '/data/tmp'

params.fastq = ''

params.star_index_dir_human = '/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/'
params.fasta_human = '/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta'
params.gtf_human = '/data/local/reference/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf'

params.dbsnp = '/data/local/reference/GATK_resource_bundle/hg38/hg38/dbsnp_146.hg38.vcf.gz'
params.known_indels = '/data/local/reference/GATK_resource_bundle/hg38/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'

params.transcriptome_fasta = '/data/bin/bcbio/genomes/Hsapiens/hg38/rnaseq/ref-transcripts.fa'
params.transcriptome_gtf = '/data/bin/bcbio/genomes/Hsapiens/hg38/rnaseq/ref-transcripts.gtf'

// For PDX
params.star_index_dir_mouse = '/data/proj/skcm_perkins/Pipelines/rna/preprocessing/rnaSeqPreprocessing/genome_mouse'
params.fasta_mouse = '/data/proj/skcm_perkins/Pipelines/rna/preprocessing/rnaSeqPreprocessing/genome_mouse/GRCm38.primary_assembly.genome.fa'
params.gtf_mouse = '/data/proj/skcm_perkins/Pipelines/rna/preprocessing/rnaSeqPreprocessing/genome_mouse/gencode.vM22.annotation.gtf'

input_reads_ch_1 = Channel.fromPath( params.fastq )
input_reads_ch_2 = Channel.fromFilePairs( params.fastq, flat:true )

// For classification with SCOPE and CUP_AI_DX
params.ext = '.htseq.counts'
params.genome = 'hg38'
params.tx2gene_fname = "/data/bin/bcbio/genomes/Hsapiens/hg38/rnaseq/tx2gene.csv"
params.map_hg19_fname = "/data/proj/cup/Investigations/rna/map.ensembl_entrez.hg19.rda"
params.map_hg38_fname = "/data/proj/cup/Investigations/rna/map.ensembl_entrez.hg38.rda"
params.features_fname = "/data/proj/cup/Investigations/rna/CUP-AI-Dx/Features_817.csv"

process FASTQC {

    maxForks 1

    publishDir "${params.outdir}/fastqc", mode: 'symlink'

    input:
    val(reads)

    output:
        path("*.zip")
        path("*.html")

    script:
    """
    mkdir ./tmp
    fastqc -q --threads ${params.cpus} --outdir "./" --dir "./tmp" ${reads}
    """
}

process STAR_ALIGN_HUMAN {

    maxForks 2

    publishDir params.is_pdx ? "${params.outdir}/STAR" : null, pattern: "*.Aligned.out.bam", mode: 'symlink'
    publishDir params.is_pdx ? "${params.outdir}/STAR" : null, pattern: "*.out", mode: 'symlink'
    publishDir params.is_pdx ? "${params.outdir}/STAR" : null, pattern: "*.sam", mode: 'symlink'
    publishDir params.is_pdx ? "${params.outdir}/STAR" : null, pattern: "*.tab", mode: 'symlink'
    publishDir params.is_pdx ? "${params.outdir}/STAR" : null, pattern: "*.Aligned.toTranscriptome.out.bam", mode: 'symlink'
    publishDir params.is_pdx ? "${params.outdir}/STAR" : null, pattern: "*.junction", mode: 'symlink'

    input:
        tuple val(sample), val(fastq_1), val(fastq_2)

	output:
    tuple val(sample), path('*.Aligned.out.bam'), emit: star_aligned_sample_bam

    shell:
    '''
	ID=$(zcat "!{fastq_1}" | head -n 1 | awk -F":" '{print $3}').$(zcat "!{fastq_1}" | head -n 1 | awk -F":" '{print $4}')
	SM="!{sample}"
	PL="ILLUMINA"
	LB="$SM"
	#read_group="ID:${ID} PL:${PL} LB:${LB} SM:${SM}"
    #echo $read_group
    #--outSAMattrRGline "$read_group" \

    STAR \
    --readFilesIn "!{fastq_1}" "!{fastq_2}" \
    --outSAMattrRGline ID:${ID} PL:${PL} LB:${LB} SM:${SM} \
    --alignIntronMax 1000000 \
    --alignIntronMin 20 \
    --alignMatesGapMax 1000000 \
    --alignSJDBoverhangMin 1 \
    --alignSJoverhangMin 8 \
    --alignSoftClipAtReferenceEnds Yes \
    --chimJunctionOverhangMin 15 \
    --chimMainSegmentMultNmax 1 \
    --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
    --chimSegmentMin 15 \
    --genomeDir "!{params.star_index_dir_human}" \
    --genomeLoad NoSharedMemory \
    --limitSjdbInsertNsj 1200000 \
    --outFileNamePrefix "./!{sample}.$ID.human." \
    --outFilterIntronMotifs None \
    --outFilterMatchNminOverLread 0.33 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.1 \
    --outFilterMultimapNmax 20 \
    --outFilterScoreMinOverLread 0.33 \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS nM NM ch \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped Within \
    --quantMode TranscriptomeSAM GeneCounts \
    --readFilesCommand zcat \
    --runThreadN 20 \
    --twopassMode Basic
    '''
}

process STAR_ALIGN_MOUSE {

    maxForks 2

    input:
        tuple val(sample), val(fastq_1), val(fastq_2)

	output:
    tuple val(sample), path('*.Aligned.out.bam'), emit: star_aligned_sample_bam

    shell:
    '''
	ID=$(zcat "!{fastq_1}" | head -n 1 | awk -F":" '{print $3}').$(zcat "!{fastq_1}" | head -n 1 | awk -F":" '{print $4}')
	SM="!{sample}"
	PL="ILLUMINA"
	LB="$SM"
	#read_group="ID:${ID} PL:${PL} LB:${LB} SM:${SM}"
    #echo $read_group
    #--outSAMattrRGline "$read_group" \

    STAR \
    --readFilesIn "!{fastq_1}" "!{fastq_2}" \
    --outSAMattrRGline ID:${ID} PL:${PL} LB:${LB} SM:${SM} \
    --alignIntronMax 1000000 \
    --alignIntronMin 20 \
    --alignMatesGapMax 1000000 \
    --alignSJDBoverhangMin 1 \
    --alignSJoverhangMin 8 \
    --alignSoftClipAtReferenceEnds Yes \
    --chimJunctionOverhangMin 15 \
    --chimMainSegmentMultNmax 1 \
    --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
    --chimSegmentMin 15 \
    --genomeDir "!{params.star_index_dir_mouse}" \
    --genomeLoad NoSharedMemory \
    --limitSjdbInsertNsj 1200000 \
    --outFileNamePrefix "./!{sample}.$ID.mouse." \
    --outFilterIntronMotifs None \
    --outFilterMatchNminOverLread 0.33 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.1 \
    --outFilterMultimapNmax 20 \
    --outFilterScoreMinOverLread 0.33 \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS nM NM ch \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped Within \
    --readFilesCommand zcat \
    --runThreadN 20 \
    --twopassMode Basic
    '''
}

process SORT_NAME_HUMAN {

    maxForks 1

    input:
        tuple val(sname), path(bam)

    output:
        tuple val(sname), path("*.name_sort.bam"), emit: sort_name_human_ch

    script:
    """
    ulimit -n 65535
    /opt/conda/envs/samtools_env/bin/samtools sort -n -@ ${params.cpus} -m 2G -O BAM -o \$(basename "${bam}" ".bam").name_sort.bam "${bam}"
    """
}

process SORT_NAME_MOUSE {

    maxForks 1

    input:
        tuple val(sname), path(bam)

    output:
        tuple val(sname), path("*.name_sort.bam"), emit: sort_name_mouse_ch

    script:
    """
    ulimit -n 65535
    /opt/conda/envs/samtools_env/bin/samtools sort -n -@ ${params.cpus} -m 2G -O BAM -o \$(basename "${bam}" ".bam").name_sort.bam "${bam}"
    """
}

process DISAMBIGUATE {

    publishDir "${params.outdir}/disambiguated", pattern: "*.*", mode: 'symlink'

    input:
        tuple val(sname), path(bam_human), path(bam_mouse)

    output:
        tuple val(sname), path("*.disambiguatedSpeciesA.bam"), emit: bam_human_disambiguated_ch

    script:
    """
    echo "${sname}"
    ngs_disambiguate -s "${sname}" -o "./" -a star "${bam_human}" "${bam_mouse}"
    """
}

process MERGE_BAMS {

    publishDir "${params.outdir}/MergeSamFiles", pattern: "*.bam", mode: 'symlink'

    input:
    tuple val(sample), val(bams)

    output:
    path("*.bam"), emit: merged_bam
    
    // In the case of one input bam file only, essentially just copies the file.
    script:
    """
    files="${[bams].flatten().join(',')}"

    result=""
    IFS=',' read -ra file_array <<< "\$files"
    for file in "\${file_array[@]}"; do
        result+=" -I=\$file"
    done
    
    result=\$(echo "\$result" | tr -d '[]')

    gatk MergeSamFiles \
      \$result \
      --SORT_ORDER unsorted \
      -O=${sample}.merged.bam
    """
}

process HTSEQ_COUNT {

    maxForks 4

	publishDir "${params.outdir}/HTSeq", mode: 'symlink'

	input:
		val(bam)
	output:
		path("*.htseq.counts"), emit: htseq_ch

    script:
    sample = bam.getSimpleName()
    """
    /opt/conda/envs/samtools_env/bin/samtools sort -n -@ ${params.cpus} -m 250MB -O BAM -o ./"${sample}".name_sort.bam "${bam}"
    #HTSeq-0.6.1p1
    /opt/conda/envs/htseq_env/bin/htseq-count \
    -f bam \
    -r name \
    -s no \
    -a 10 \
    -t exon \
    -i gene_id \
    -m intersection-nonempty \
    "${sample}".name_sort.bam \
    "${params.gtf_human}" > ./${sample}.htseq.counts

    rm ./"${sample}".name_sort.bam
    """
}

process CLASSIFICATION {

    maxForks 6

    publishDir "${params.outdir}/classification", mode: 'symlink'

    input:
        val(htseq_counts)
    output:
        path("*")

    script:
    """
    Rscript --vanilla -e 'library(rnaSeqPanCanClassifier); classify(
            input_test="${htseq_counts}", input_train = "tcga", outdir = "./", k = 6,
            subset_to=c(), genome = "hg38", labels_train = NULL,
            extension_test = ".htseq.counts", extension_train = ".gene_counts",
            override_labels = F, is_matrix = F
    )'
    """
}

process PREPARE_COUNT_MATRIX_SCOPE {

    input:
        val(htseq_counts)

    output:
        path("*.txt"), emit: matrix_path

    script:
    """
    Rscript --vanilla -e 'library(tidyverse); library(Biobase);

    sname <- str_split_fixed(basename("${htseq_counts}"),"\\\\.",2)[,1]
    indir <- paste0("./tmp/", sname)
    dir.create(indir, recursive = T, showWarnings = F)
    system(paste0("ln -s ${htseq_counts} ", indir, "/"))

    inhouse <- rnaSeqPanCanClassifier:::.setup_expression(indir, extension = "${params.ext}")
    g_info <- rnaSeqPanCanClassifier:::.get_gene_info(genes = rownames(inhouse), genome = "${params.genome}", symbols = T)
    inhouse <- inhouse[rownames(g_info),]
    inhouse <- rnaSeqPanCanClassifier:::.rpkm_eset(inhouse, gene_info = g_info)

    inhouse.mat <- as.data.frame(exprs(inhouse), stringsAsFactors=F)
    inhouse.mat <- cbind(ENSEMBL = rownames(inhouse.mat), inhouse.mat)
    write.table(inhouse.mat, file = paste0("./",sname,".rpkm.txt"), col.names = T,
                row.names = F, sep = "\t", quote = F)'
    """
}

process PREPARE_COUNT_MATRIX_CUP_AI_DX {

    input:
        path(transcript_counts)

    output:
        path("*.csv"), emit: matrix_path

    script:
    """
    Rscript --vanilla -e '
    library(tidyverse)
    library(tximport)

    sname <- basename(dirname(normalizePath("${transcript_counts}")))
    fnames <- setNames("${transcript_counts}",sname)

    tx2gene <- read.table("${params.tx2gene_fname}", header=F, sep=",")
    colnames(tx2gene) <- c("tx","gene")

    txi.kallisto <- tximport(fnames, type = "kallisto", txOut = F,
                             countsFromAbundance = "lengthScaledTPM",
                             tx2gene = tx2gene,
                             ignoreTxVersion = F,
                             ignoreAfterBar = F,dropInfReps=T)

    log2_tpm_p1 <- log2(txi.kallisto[["counts"]]+1)

    extract_features <- function(log2_tpm_p1){
      databases <- list()
      databases[["hg19"]][["map"]] <- readRDS(
        file = "${params.map_hg19_fname}")
      databases[["hg38"]][["map"]] <- readRDS(
        file = "${params.map_hg38_fname}")

      features <- read.csv("${params.features_fname}", header=T)
      features <- features[,2]
      features <- gsub(features, pattern = "^X", replacement = "")

      ens_keep <- union(databases[["hg38"]][["map"]][["ensembl_gene_id"]][
          databases[["hg38"]][["map"]][["entrezgene_id"]] %in% features],
            databases[["hg19"]][["map"]][["ensembl_gene_id"]][
                databases[["hg19"]][["map"]][["entrezgene_id"]] %in% 
                setdiff(features,as.character(databases[["hg38"]][["map"]][["entrezgene_id"]]))])

      log2_tpm_p1.minimal <- log2_tpm_p1[rownames(log2_tpm_p1) %in% ens_keep,]
      log2_tpm_p1.minimal <- as.data.frame(log2_tpm_p1.minimal,stringsAsFactors=F)

      map <- rbind(databases[["hg38"]][["map"]][databases[["hg38"]][["map"]][["ensembl_gene_id"]] %in% ens_keep,],
        databases[["hg19"]][["map"]][databases[["hg19"]][["map"]][["ensembl_gene_id"]] %in% 
        databases[["hg19"]][["map"]][["ensembl_gene_id"]][databases[["hg19"]][["map"]][["entrezgene_id"]] %in% 
        setdiff(features,as.character(databases[["hg38"]][["map"]][["entrezgene_id"]]))],])

      stopifnot(all(rownames(log2_tpm_p1.minimal) %in% map[["ensembl_gene_id"]]))

      tmp <- merge(log2_tpm_p1.minimal, map, by.x="row.names", by.y="ensembl_gene_id", all.x=T, all.y=F)

      tmp <- tmp[tmp[["entrezgene_id"]] %in% features,]
      rownames(tmp) <- NULL
      tmp[["Row.names"]] <- NULL
      dup_ez <- names(which(table(tmp[["entrezgene_id"]])>1))
      tmp.keep <- list()
      for (dup in dup_ez){
        idx <- tmp[["entrezgene_id"]] == dup
        tmp.keep[[dup]] <- tmp[idx,][which.max(apply(tmp[idx,],1,mean)),]
      }
      tmp.keep <- do.call("rbind",tmp.keep)

      tmp <- rbind(tmp[! tmp[["entrezgene_id"]] %in% names(which(table(tmp[["entrezgene_id"]])>1)),],tmp.keep)
      stopifnot(nrow(tmp)==length(features))
      rownames(tmp) <- tmp[["entrezgene_id"]]
      tmp[["entrezgene_id"]] <- NULL
      stopifnot(all(rownames(tmp) %in% features))
      stopifnot(all(features %in% rownames(tmp)))
      tmp[["dummy"]] <- ""
      tmp <- tmp[features,]
      tmp[["dummy"]] <- NULL
      stopifnot(identical(features,rownames(tmp)))

      log2_tpm_p1.minimal <- tmp

      return(log2_tpm_p1.minimal)
    }

    format_output <- function(log2_tpm_p1.minimal){
      log2_tpm_p1.minimal <- as.data.frame(t(log2_tpm_p1.minimal))
      colnames(log2_tpm_p1.minimal) <- paste0("X",colnames(log2_tpm_p1.minimal))
      log2_tpm_p1.minimal[["X"]] <- rownames(log2_tpm_p1.minimal)
      rownames(log2_tpm_p1.minimal) <- NULL
      log2_tpm_p1.minimal <- log2_tpm_p1.minimal[,c(which(colnames(log2_tpm_p1.minimal)=="X"),which(colnames(log2_tpm_p1.minimal)!="X"))]
      colnames(log2_tpm_p1.minimal)[colnames(log2_tpm_p1.minimal)=="X"] <- ""

      return(log2_tpm_p1.minimal)
    }

    log2_tpm_p1.minimal <- extract_features(log2_tpm_p1)
    log2_tpm_p1.minimal <- format_output(log2_tpm_p1.minimal)

    write.table(log2_tpm_p1.minimal,
                file = paste0("./",sname,".log2_tpm_p1.csv"),
                sep = ",", col.names = T, row.names = F)'
    """
}

process CLASSIFY_SCOPE {
    
    maxForks 10

    publishDir "${params.outdir}/scope", mode: 'copy'

    input:
        path(matrix_path)

    output:
        path("*")

    script:
    """
    data_matrix="\$(readlink -f ${matrix_path})"
    sname="\$(basename \$data_matrix | awk -F".rpkm" '{print \$1}')"
    outdir="\$sname"
    if [ ! -d \$outdir ];
    then
        mkdir \$outdir
    fi

    singularity exec \
        -B "\$outdir":/output \
        -B "\$data_matrix":/rpkm.txt \
        docker://jowkar/scope_container:latest /usr/bin/python3.7 -c \
        'import cancerscope as cs; \
        scope_obj = cs.scope(); \
        predictions_from_file = scope_obj.get_predictions_from_file("/rpkm.txt",outdir = "/output")'
    """
}

process CLASSIFY_CUP_AI_DX {

    maxForks 10

    publishDir "${params.outdir}/cup_ai_dx", mode: 'copy'

    input:
        path matrix_path

    output:
        path("*")

    script:
    """
    data_csv="\$(readlink -f ${matrix_path})" # absolute path
    sname="\$(basename \$data_csv | awk -F".log2" '{print \$1}')"
    outdir="\$sname"
    if [ ! -d \$outdir ];
    then
        mkdir \$outdir
    fi

    singularity exec \
        -B "\$outdir":/scripts/output \
        -B "\$data_csv":/scripts/data/\$sname.csv \
        --no-home \
        --fakeroot \
        docker://yuz12012/ai4cancer:product \
        /bin/bash -rcfile /root/.bashrc -c "source /root/.bashrc && conda activate tf-cpu && cd /scripts/ && python RunPrediction.py --data_set=data/\$sname.csv"
    """
}

process SORT_INDEX {

    maxForks 1

    input:
        val(bam)
    output:
        path("*.sorted.bam"), emit: sort_index_ch

    script:
    """
    ulimit -n 65535
    /opt/conda/envs/samtools_env/bin/samtools sort -@ ${params.cpus} -m 250MB -O BAM -o \$(basename "${bam}" ".bam").sorted.bam "${bam}"
    /opt/conda/envs/samtools_env/bin/samtools index \$(basename "${bam}" ".bam").sorted.bam
    """
}

process MARK_DUPLICATES {

    maxForks 4

    input:
        val(bam)
    output:
        path("*.MarkDuplicates.bam"), emit: markduplicates_ch

    script:
    """
    gatk MarkDuplicates \
        --INPUT="${bam}" \
        --OUTPUT=\$(basename "${bam}" ".bam").MarkDuplicates.bam \
        --METRICS_FILE=\$(basename "${bam}" ".bam").MarkDuplicates.metrics.txt \
        --CREATE_INDEX true
    """
}

process SPLIT_N_TRIM {

    maxForks 4

    input:
        val(bam)
    output:
        path("*.split.bam"), emit: split_ch

    script:
    """
    gatk SplitNCigarReads --tmp-dir="${params.tmpdir}" -R="${params.fasta_human}" -I="${bam}" -O=\$(basename "${bam}" ".bam").split.bam
    """
}

process RECALIBRATION {

    maxForks 4

	publishDir "${params.outdir}/Recalibrated", mode: 'symlink'

    input:
        val(bam)
    output:
        path("*.pass2.bam"), emit: recal_ch

    script:
    """
    gatk BaseRecalibrator \
        -I="${bam}" \
        -R="${params.fasta_human}" \
        --known-sites="${params.dbsnp}" \
        --known-sites="${params.known_indels}" \
        -O="${bam}".recal_pass1.table

    gatk ApplyBQSR \
        -I="${bam}" \
        -R="${params.fasta_human}" \
        --bqsr-recal-file="${bam}".recal_pass1.table \
        -O=\$(basename "${bam}" ".bam").recal.pass1.bam

    gatk BaseRecalibrator \
        -I=\$(basename "${bam}" ".bam").recal.pass1.bam \
        -R="${params.fasta_human}" \
        --known-sites="${params.dbsnp}" \
        --known-sites="${params.known_indels}" \
        -O=\$(basename "${bam}" ".bam").recal.pass1.bam.recal_pass2.table

    gatk ApplyBQSR \
        -I=\$(basename "${bam}" ".bam").recal.pass1.bam \
        -R="${params.fasta_human}" \
        --bqsr-recal-file=\$(basename "${bam}" ".bam").recal.pass1.bam.recal_pass2.table \
        -O=\$(basename "${bam}" ".bam").recal.pass2.bam
    """
}

process VARIANT_CALLING {

	publishDir "${params.outdir}/variant_calling", mode: 'symlink'

    input:
        val(bam)
    output:
        path("*.vcf.gz"), emit: variant_calling_ch

    script:
    """
    gatk HaplotypeCaller \
        -R="${params.fasta_human}" \
        -I="${bam}" \
        --dont-use-soft-clipped-bases=true \
        --standard-min-confidence-threshold-for-calling=20.0 \
        -O=\$(basename "${bam}").haplotypecaller.vcf.gz
    """
}

process FILTER_VARIANTS {
    
    input:
        val(vcf)
    output:
        path("*.filtered.vcf.gz"), emit: filter_variants_ch

    script:
    """
    gatk VariantFiltration \
        -R="${params.fasta_human}" \
        -V="${vcf}" \
        --cluster-window-size 35 \
        --cluster-size 3 \
        --filter-name FS \
        --filter-expression "FS > 30.0" \
        --filter-name QD \
        --filter-expression "QD < 2.0" \
        -O=\$(basename "${vcf}" ".vcf.gz").filtered.vcf.gz
    """
}

process REMOVE_FAILING {

    input:
        val(vcf)
    output:
        path("*.pass.vcf.gz"), emit: remove_failing_ch

    script:
    """
    outname=\$(basename "${vcf}" ".vcf.gz").pass.vcf
    cat <(bgzip -d -c "${vcf}" | grep "#") <(bgzip -d -c "${vcf}" | grep -v "#" | grep -w PASS) > "\$outname"
    bgzip -f "\$outname"
    tabix -f "\$outname".gz
    """
}

process VEP {

	publishDir "${params.outdir}/variant_calling", mode: 'symlink'

    input:
        val(vcf_gz)
    output:
        path("output_vep_updated/*.ann.vcf"), emit: vep_ch

    script:
    """
    vcf=\$(basename "${vcf_gz}" ".vcf.gz").vcf
    bgzip -d -c -f "${vcf_gz}" > "\$vcf"

    if [ ! -d output_vep_updated ]
    then
        mkdir output_vep_updated
    fi
    singularity exec \
        -B \$(pwd)/output_vep_updated:/output_vep_updated \
        -B /data/bin/vep_cache:/.vep \
        -B "\$vcf":/\$(basename "\$vcf") \
        /data/proj/um_perkins/Pipelines/rna/preprocessing/rnaSeqPreprocessing/vep.sif /opt/vep/src/ensembl-vep/vep \
        --species homo_sapiens \
        --assembly GRCh38 \
        --offline \
        --cache \
        --dir /.vep \
        --input_file /\$(basename "\$vcf") \
        --output_file /output_vep_updated/\$(basename "\$vcf" ".vcf").ann.vcf \
        --everything \
        --vcf \
        --fasta /.vep/homo_sapiens/105_GRCh38/Homo_sapiens_assembly38.fasta \
        --force_overwrite
    """
}

process VCF2MAF {

	publishDir "${params.outdir}/variant_calling", mode: 'symlink'

    input:
        val(vcf)
    output:
        path("*.maf"), emit: vcf2maf_ch

    script:
    """
    tumor_id=\$(basename "${vcf}" ".ann.vcf" | awk -F"_L00" '{print \$1}' | awk -F"." '{print \$1}')
    vcf2maf.pl \
        --inhibit-vep \
        --input-vcf "${vcf}" \
        --output-maf \$(basename "${vcf}" ".vcf").maf \
        --tumor-id "\$tumor_id" \
        --ref-fasta /data/bin/vep_cache/homo_sapiens/105_GRCh38/Homo_sapiens_assembly38.fasta \
        --ncbi-build GRCh38
    """
}

process BAM_TO_FASTQ {

    input:
        tuple val(sname), path(bam)

    output:
        tuple val(sname), path("${sname}_human_R1.fastq.gz"), path("${sname}_human_R2.fastq.gz"), emit:bam_to_fastq_ch

    script:
    """
    TMP_DIR=./\$(basename "${bam}" ".bam")_tmp

    # Sort by read name
    gatk SortSam \
        --TMP_DIR=\$TMP_DIR \
        -I ${bam} \
        -O queryname_sorted.bam \
        -SO queryname

    TMP_DIR_2=./\$(basename "${bam}" ".bam")_tmp_2

    # Convert to FASTQ while keeping the full reads (including soft trimmed parts)
    gatk SamToFastq \
        --TMP_DIR=\$TMP_DIR_2 \
        -I queryname_sorted.bam \
        -F ${sname}_human_R1.fastq.gz \
        -F2 ${sname}_human_R2.fastq.gz \
        --VALIDATION_STRINGENCY LENIENT \
        --INCLUDE_NON_PF_READS true

    rm queryname_sorted.bam
    """
}

process KALLISTO_INDEX {

    publishDir "${params.outdir}/kallisto", mode: 'symlink'

    output:
        path("hg38_index"), emit: kallisto_index_ch

    script:
    """
    kallisto index -i hg38_index ${params.transcriptome_fasta}
    """
}

process KALLISTO {

    maxForks 2

    publishDir "${params.outdir}/kallisto", mode: 'symlink'

    input:
        tuple val(sample), val(fastq_1), val(fastq_2), path(hg38_index)

    output:
        path("${sample}/fusion.txt"), emit: kallisto_ch
        path("${sample}/abundance.h5"), emit: kallisto_h5_ch

    script:
    """
    mkdir "./${sample}"
    kallisto quant -i "${hg38_index}" -o "./${sample}/" -b 100 --fusion --gtf ${params.transcriptome_gtf} --threads 31 "${fastq_1}" "${fastq_2}"
    """
}

process PIZZLY {

    publishDir "${params.outdir}/pizzly", mode: 'symlink'

    input:
        val(fusion_txt)

    output:
        path("*"), emit: pizzly_ch

    script:
    sample = fusion_txt.getParent().getSimpleName()
    """
    pizzly \
        -k 31 \
        --gtf "${params.transcriptome_gtf}" \
        --cache /data/proj/um_perkins/Pipelines/rna/preprocessing/index.cache.fixed.txt \
        --align-score 2 \
        --insert-size 400 \
        --fasta "${params.transcriptome_fasta}" \
        --output "${sample}" \
        "${fusion_txt}"
    """
}

workflow {
    FASTQC(
		input_reads_ch_1
	)

    starOutputHuman = STAR_ALIGN_HUMAN(
		input_reads_ch_2
	)

    if (params.is_pdx){
        starOutputMouse = STAR_ALIGN_MOUSE(
		    input_reads_ch_2
	    )
    }

    sortNameOutputHuman = SORT_NAME_HUMAN(
        starOutputHuman
    )

    if (params.is_pdx){
        sortNameOutputMouse = SORT_NAME_MOUSE(
            starOutputMouse
        )
    }

    if (params.is_pdx){
        sortNameHumanMouse = sortNameOutputHuman.join(
            sortNameOutputMouse, by: 0)

        disambiguate_ch = DISAMBIGUATE(
            sortNameHumanMouse
        )
    }

    if (params.is_pdx){
        starOutputMerged = disambiguate_ch.view()
    } else {
        starOutputMerged = starOutputHuman
    }

    .collate(1).map{ it -> it.flatten() }
    .map { sample, file -> 
        def sname = sample.replaceAll('_L00.*', '')
        return tuple(sname, file)
    }.groupTuple().view()
    
    merged_bam_ch = MERGE_BAMS( starOutputMerged )

    htseq_ch = HTSEQ_COUNT(
		merged_bam_ch
	)

    CLASSIFICATION(
        htseq_ch
    )

    sort_index_ch = SORT_INDEX(
        merged_bam_ch
    )

    markduplicates_ch = MARK_DUPLICATES(
        sort_index_ch
    )

    split_ch = SPLIT_N_TRIM(
        markduplicates_ch
    )

    recal_ch = RECALIBRATION(
        split_ch
    )

    variant_calling_ch = VARIANT_CALLING(
        recal_ch
    )

    filter_variants_ch = FILTER_VARIANTS(
        variant_calling_ch
    )

    remove_failing_ch = REMOVE_FAILING(
        filter_variants_ch
    )

    vep_ch = VEP(
        remove_failing_ch
    )

    VCF2MAF(
        vep_ch
    )

    kallisto_index_ch = KALLISTO_INDEX()

    if (params.is_pdx){
        pdx_merged_ch = merged_bam_ch.map { file ->
            def basename = file.baseName.tokenize('.').first()
            return tuple(basename, file)
        }
        bam_to_fastq_from_merged_ch = BAM_TO_FASTQ(pdx_merged_ch)
        input_reads_ch_with_index = bam_to_fastq_from_merged_ch.combine(kallisto_index_ch)
    } else {
        input_reads_ch_with_index = input_reads_ch_2.combine(kallisto_index_ch)
    }

    KALLISTO(
        input_reads_ch_with_index
    )

    PIZZLY (
        KALLISTO.out.kallisto_ch
    )

    scope_matrix_path_ch = PREPARE_COUNT_MATRIX_SCOPE(
        htseq_ch
    )

    cup_ai_dx_matrix_path_ch = PREPARE_COUNT_MATRIX_CUP_AI_DX(
        KALLISTO.out.kallisto_h5_ch
    )

    CLASSIFY_SCOPE(
        scope_matrix_path_ch
    )

    CLASSIFY_CUP_AI_DX(
        cup_ai_dx_matrix_path_ch
    )
}
