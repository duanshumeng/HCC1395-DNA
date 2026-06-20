task sentieon_TNscope{
    String sample_id
    File tumor_bam
    File tumor_bam_bai
    File? normal_bam
    File? normal_bam_bai
    String tumor_name
    String normal_name
    File tumor_recall_data
    File normal_recall_data

    File ref_dir
    String ref_fasta
    File dbsnp_dir
    String dbsnp

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

        sentieon driver -t $nt -r ${ref_dir}/${ref_fasta} \
        -i ${tumor_bam} -q ${tumor_recall_data} \
        -i ${normal_bam} -q ${normal_recall_data} \
        --algo TNscope --tumor_sample ${tumor_name} --normal_sample ${normal_name} \
        --trim_soft_clip \
        --dbsnp ${dbsnp_dir}/${dbsnp} ${sample_id}.TNscope.vcf || { echo "TNscope failed"; exit 1; }
        
        # --disable_detector sv 
        ls ./

    >>>

    runtime{
        docker:docker
        cluster:cluster_config
        systemDisk:"cloud_ssd 40"
        dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
    }

    output{
        File vcf = "${sample_id}.TNscope.vcf"
        File vcf_index = "${sample_id}.TNscope.vcf.idx"

    }
}