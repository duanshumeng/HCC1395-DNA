task SentieonFastqToBam {
    # 工具输入文件和参数
    File fastq1
    File fastq2
    String sample_id
    String Seq_platform
    String ref_fasta
    File ref_fasta_dir
    String SENTIEON_LICENSE

    String docker
    String cluster_config
    String disk_size

    ## Extra driver parameters
    String qc_driver_args = ""
    String lc_driver_args = "--traverse_param=200000/10000"
    String dedup_driver_args = "--traverse_param=200000/10000"
    ## Extra algo parameters
    String bwa_args = "-Y -M"
    String bwa_chunk_size = "100000000"
    String lc_args = ""
    String bam_option = "--bam_compression 1"



    
    String out_bam = sample_id + ".dedup.bam"
    String out_bai = sample_id + ".dedup.bam.bai"

    # 工具运行命令
    command <<<
        set -exo pipefail
        export SENTIEON_LICENSE=${SENTIEON_LICENSE}
        nt=$(nproc)
        
        sentieon bwa mem -R "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:${Seq_platform}" ${bwa_args} -K ${bwa_chunk_size} -t $nt ${ref_fasta_dir}/${ref_fasta} ${fastq1} ${fastq2} \
        | sentieon util sort ${bam_option} -i - -r ${ref_fasta_dir}/${ref_fasta} -t $nt -o ${sample_id}.sorted.bam --sam2bam

        ls ./
                  
        sentieon driver -r ${ref_fasta_dir}/${ref_fasta} -t $nt -i ${sample_id}.sorted.bam ${qc_driver_args} \
        --algo MeanQualityByCycle ${sample_id}.mq_metrics.txt \
        --algo QualDistribution ${sample_id}.qd_metrics.txt \
        --algo GCBias --summary ${sample_id}.gc_summary_metrics.txt ${sample_id}.gc_metrics.txt \
        --algo AlignmentStat ${sample_id}.aln_metrics.txt \
        --algo InsertSizeMetricAlgo ${sample_id}.is_metrics.txt
                  
        ls ./

        sentieon driver -r ${ref_fasta_dir}/${ref_fasta} -t $nt -i ${sample_id}.sorted.bam ${lc_driver_args} \
         --algo LocusCollector \
         ${lc_args} \
         ${sample_id}.score.txt.gz
                  
        ls ./

        sentieon driver -r ${ref_fasta_dir}/${ref_fasta} -t $nt -i ${sample_id}.sorted.bam ${dedup_driver_args} \
         --algo Dedup \
         --score_info ${sample_id}.score.txt.gz \
         --metrics ${sample_id}.dedup_metrics.txt \
         ${bam_option} ${out_bam} 
         ls ./

    >>>


    runtime {
        docker:docker
        cluster:cluster_config
        systemDisk:"cloud_ssd 40"
        dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
        timeout:259200

    }
    

    # 工具运行输出结果
    output {
        File deduped_bam = out_bam
        File deduped_bam_bai = out_bai
        Array[File] qc_metrics = glob("*_metrics.txt")
    }

}