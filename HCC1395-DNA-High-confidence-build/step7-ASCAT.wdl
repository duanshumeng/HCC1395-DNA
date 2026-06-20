task ASCAT{
    String sample
    File tumor_bam
    File tumor_bam_index
    File? normal_bam
    File? normal_bam_index
    File? bed_file
    String gender
    String ref_genome
    String docker
    String cluster_config
    String disk_size

    command <<<
        set -o pipefail
        set -e
        
        out_dir='./'
        
        if [ ${bed_file} ]; then
           
           /root/miniconda3/envs/hrd/bin/Rscript /home/ASCAT/ascat.r ${tumor_bam} ${normal_bam} ${sample}'_T' ${sample}'_N' ${sample} WES $out_dir ${gender} ${ref_genome} ${bed_file}
        else     
         	/root/miniconda3/envs/hrd/bin/Rscript /home/ASCAT/ascat.r ${tumor_bam} ${normal_bam} ${sample}'_T' ${sample}'_N' ${sample} WGS $out_dir ${gender} ${ref_genome}

        fi
        ls ./

        tar cvf ${sample}.tar ./

    >>>

    runtime{
        docker:docker
        cluster:cluster_config
        systemDisk:"cloud_ssd 40"
        dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
    }

    output{
        File out_file = "${sample}.tar"
    }

}