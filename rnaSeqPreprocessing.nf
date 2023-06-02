#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.cpus = 63
params.fastq = '/data/proj/um_perkins/Pipelines/rna/preprocessing/fastqdir/*/*_R{1,2}_*.fastq.gz'
params.star_index_dir = '/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/'
params.fasta = '/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta'
params.gtf = '/data/local/reference/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf'
params.dbsnp = '/data/local/reference/GATK_resource_bundle/hg38/hg38/dbsnp_146.hg38.vcf.gz'
params.known_indels = '/data/local/reference/GATK_resource_bundle/hg38/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
params.outdir = './results/'
params.tmpdir = '/data/tmp'
params.transcriptome_fasta = '/data/bin/bcbio/genomes/Hsapiens/hg38/rnaseq/ref-transcripts.fa'
params.transcriptome_gtf = '/data/bin/bcbio/genomes/Hsapiens/hg38/rnaseq/ref-transcripts.gtf'

input_reads_ch_1 = Channel.fromPath( params.fastq )
input_reads_ch_2 = Channel.fromFilePairs( params.fastq, flat:true )

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

process STAR_ALIGN {

    maxForks 2

	publishDir "${params.outdir}/STAR", pattern: "*.Aligned.out.bam", mode: 'symlink'
	publishDir "${params.outdir}/STAR", pattern: "*.out", mode: 'symlink'
	publishDir "${params.outdir}/STAR", pattern: "*.sam", mode: 'symlink'
	publishDir "${params.outdir}/STAR", pattern: "*.tab", mode: 'symlink'
	publishDir "${params.outdir}/STAR", pattern: "*.Aligned.toTranscriptome.out.bam", mode: 'symlink'
    publishDir "${params.outdir}/STAR", pattern: "*.junction", mode: 'symlink'

    input:
        tuple val(sample), val(fastq_1), val(fastq_2)

	output:
	path('*.Aligned.out.bam'), emit: star_aligned
	path "*.out", emit: alignment_logs
    path "*.sam"
    path "*.tab"
    path "*.Aligned.toTranscriptome.out.bam"
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
    --genomeDir "!{params.star_index_dir}" \
    --genomeLoad NoSharedMemory \
    --limitSjdbInsertNsj 1200000 \
    --outFileNamePrefix "./!{sample}.$ID." \
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

process MERGE_BAMS {

    publishDir "${params.outdir}/MergeSamFiles", pattern: "*.bam", mode: 'symlink'

    input:
    tuple val(sample), val(bams)
    // tuple val(sample), val(r1_files), val(r2_files)

    output:
    path("*.bam"), emit: merged_bam
    
    script:
    """
    files="${bams.join(',')}"

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
    /data/miniconda3/envs/cup_star_htseq/bin/samtools sort -n -@ ${params.cpus} -m 250MB -O BAM -o ./"${sample}".name_sort.bam "${bam}"
    #HTSeq-0.6.1p1
    htseq-count \
    -f bam \
    -r name \
    -s no \
    -a 10 \
    -t exon \
    -i gene_id \
    -m intersection-nonempty \
    "${sample}".name_sort.bam \
    "${params.gtf}" > ./${sample}.htseq.counts

    rm ./"${sample}".name_sort.bam
    #"${bam}" \
    #rm \$(readlink -f "${bam}")
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
    #!/usr/bin/env Rscript --vanilla
    library(classification)
    library(stringr)
    outdir <- paste0("./",str_split_fixed(basename("${htseq_counts}"),"\\\\.",2)[,1])
    res <- classify(input_test="${htseq_counts}",genome="hg38",outdir=outdir)
    s <- sessionInfo()
    saveRDS(s,paste0(outdir,"sessionInfo.rda"))
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
    /data/miniconda3/envs/cup_star_htseq/bin/samtools sort -@ ${params.cpus} -m 250MB -O BAM -o \$(basename "${bam}" ".bam").sorted.bam "${bam}"
    /data/miniconda3/envs/cup_star_htseq/bin/samtools index \$(basename "${bam}" ".bam").sorted.bam
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
        --METRICS_FILE=\$(basename "${bam}" ".bam").MarkDuplicates.metrics.txt

    /data/miniconda3/envs/cup_star_htseq/bin/samtools index -@ ${params.cpus} \$(basename "${bam}" ".bam").MarkDuplicates.bam
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
    gatk SplitNCigarReads --tmp-dir="${params.tmpdir}" -R="${params.fasta}" -I="${bam}" -O=\$(basename "${bam}" ".bam").split.bam
    """
}

process RECALIBRATION {

    maxForks 4

	publishDir "${params.outdir}/Recalibrated", mode: 'symlink'

    input:
        val(bam)
    output:
        path("*.pass2.bam"), emit: recal_ch
        path("*.pass2.bai")

    script:
    """
    gatk BaseRecalibrator \
        -I="${bam}" \
        -R="${params.fasta}" \
        --known-sites="${params.dbsnp}" \
        --known-sites="${params.known_indels}" \
        -O="${bam}".recal_pass1.table

    gatk ApplyBQSR \
        -I="${bam}" \
        -R="${params.fasta}" \
        --bqsr-recal-file="${bam}".recal_pass1.table \
        -O=\$(basename "${bam}" ".bam").recal.pass1.bam

    gatk BaseRecalibrator \
        -I=\$(basename "${bam}" ".bam").recal.pass1.bam \
        -R="${params.fasta}" \
        --known-sites="${params.dbsnp}" \
        --known-sites="${params.known_indels}" \
        -O=\$(basename "${bam}" ".bam").recal.pass1.bam.recal_pass2.table

    gatk ApplyBQSR \
        -I=\$(basename "${bam}" ".bam").recal.pass1.bam \
        -R="${params.fasta}" \
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
        path("*.vcf.gz.tbi")

    script:
    """
    gatk HaplotypeCaller \
        -R="${params.fasta}" \
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
        -R="${params.fasta}" \
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
        path("*.pass.vcf.gz.tbi")

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
        path("*")

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

    STAR_ALIGN(
		input_reads_ch_2
	)

    starOutputMerged = STAR_ALIGN.out.star_aligned_sample_bam
    .collate(1).map{ it -> it.flatten() }
    .map { sample, file -> 
        def sname = sample.replaceAll('_L00.*', '')
        return tuple(sname, file)
    }.groupTuple().view()
    
    MERGE_BAMS( starOutputMerged )

    HTSEQ_COUNT(
		MERGE_BAMS.out.merged_bam
	)

    CLASSIFICATION(
        HTSEQ_COUNT.out.htseq_ch
    )

    SORT_INDEX(
        MERGE_BAMS.out.merged_bam
    )

    MARK_DUPLICATES(
        SORT_INDEX.out.sort_index_ch
    )

    SPLIT_N_TRIM(
        MARK_DUPLICATES.out.markduplicates_ch
    )

    RECALIBRATION(
        SPLIT_N_TRIM.out.split_ch
    )

    VARIANT_CALLING(
        RECALIBRATION.out.recal_ch
    )

    FILTER_VARIANTS(
        VARIANT_CALLING.out.variant_calling_ch
    )

    REMOVE_FAILING(
        FILTER_VARIANTS.out.filter_variants_ch
    )

    VEP(
        REMOVE_FAILING.out.remove_failing_ch
    )

    VCF2MAF(
        VEP.out.vep_ch
    )

    KALLISTO_INDEX()

    input_reads_ch_with_index=input_reads_ch_2.combine(KALLISTO_INDEX.out.kallisto_index_ch)

    KALLISTO(
        input_reads_ch_with_index
    )

    PIZZLY (
        KALLISTO.out.kallisto_ch
    )
}
