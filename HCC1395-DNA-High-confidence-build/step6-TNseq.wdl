task sentieon_TNseq{
    String sample_id
    File tumor_bam
    File tumor_bam_bai
    File? normal_bam
    File? normal_bam_bai
    String tumor_name
    String normal_name

    File ref_dir
    String ref_fasta
    File germline_resource
    File germline_resource_tbi
    File? regions
    Int? interval_padding

    # excute env
    String docker
    String cluster_config
    String disk_size
    String SENTIEON_LICENSE


    command <<<
        set -o pipefail
        set -exo
        export SENTIEON_LICENSE=${SENTIEON_LICENSE}
        nt=$(nproc)

        if [ ${regions} ]; then
            INTERVAL="--interval ${regions} --interval_padding ${interval_padding}"
        else
            INTERVAL=""
        fi

        #if [ ${normal_bam} ]; then
            INPUT="-i ${tumor_bam} -i ${normal_bam}"
            SAMPLE="--tumor_sample ${tumor_name} --normal_sample ${normal_name}"
        #else
        #   INPUT="-i ${tumor_bam}"
        #   SAMPLE="--tumor_sample ${tumor_name}"
        #fi

        sentieon driver -t $nt -r ${ref_dir}/${ref_fasta} \
        $INPUT $INTERVAL \
        --algo TNhaplotyper2 $SAMPLE \
        --germline_vcf ${germline_resource} \
        ${sample_id}.TNseq.raw.vcf \
        --algo OrientationBias --tumor_sample ${tumor_name} \
        ${sample_id}.orientation \
        --algo ContaminationModel $SAMPLE \
        --vcf ${germline_resource} \
        --tumor_segments ${sample_id}.contamination.segments \
        ${sample_id}.contamination

        sentieon driver -t $nt \
        -r ${ref_dir}/${ref_fasta} \
        --algo TNfilter $SAMPLE \
        -v ${sample_id}.TNseq.raw.vcf \
        --contamination ${sample_id}.contamination \
        --tumor_segments ${sample_id}.contamination.segments \
        --orientation_priors ${sample_id}.orientation \
        ${sample_id}.bwa_TNseq.vcf

    >>>

    runtime{
        docker:docker
        cluster:cluster_config
        systemDisk:"cloud_ssd 40"
        dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"

    }

    output{
        File raw_vcf = "${sample_id}.TNseq.raw.vcf"
        File raw_vcf_index = "${sample_id}.TNseq.raw.vcf.idx"
        File vcf = "${sample_id}.bwa_TNseq.vcf"
        File vcf_index = "${sample_id}.bwa_TNseq.vcf.idx"
        File contamination = "${sample_id}.contamination"
        File contamination_segments = "${sample_id}.contamination.segments"
        File orientation = "${sample_id}.orientation"


    }
}