task strelka_calling{
    File tumor_bam
    File tumor_bam_bai
    File normal_bam
    File normal_bam_bai
    String ref_fasta
    File ref_dir
    String sample_id
    String ref_bed
    File? manta_indel_vcf
    File? manta_indel_vcf_index
    
    # excute env
    String docker
    String cluster_config
    String disk_size
    
    String out_dir = "${sample_id}_result"
    command <<<
    set -exo pipefail
    nt=$(nproc)
    /home/biosoft/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam ${normal_bam} \
    --tumorBam ${tumor_bam} \
    --referenceFasta ${ref_dir}/${ref_fasta} \
    --callRegions ${ref_dir}/${ref_bed} \
    --runDir ${out_dir}
                  
    #--callRegions ${ref_dir}/${ref_bed} \
    
    ls ${out_dir}

    python2.7 ${out_dir}/runWorkflow.py -m local -j $nt

    ls ${out_dir}

    tar cvf ${out_dir}.tar ${out_dir}
    >>>

    runtime{
		docker:docker
		cluster:cluster_config
		systemDisk:"cloud_ssd 40"
		dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
    }

    output{
        File out_file = "${out_dir}.tar"
        File indel_vcf = "${out_dir}/results/variants/somatic.indels.vcf.gz"
        File indel_vcf_index = "${out_dir}/results/variants/somatic.indels.vcf.gz.tbi"
        File snv_vcf = "${out_dir}/results/variants/somatic.snvs.vcf.gz"
        File snv_vcf_index = "${out_dir}/results/variants/somatic.snvs.vcf.gz.tbi"
        
    }
}