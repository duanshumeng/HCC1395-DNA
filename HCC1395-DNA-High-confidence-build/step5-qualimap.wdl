task qualimap{
    String sample_id
    File bam_file
    File bam_bai
    File annot_gff

    String docker
    String cluster_config
    String disk_size

    String out_dir = sample_id+'_BamQC'

    command <<<
        set -o pipefail
        set -exo
        nt=$(nproc)
        /opt/qualimap/qualimap bamqc -bam ${bam_file} -gff ${annot_gff} -outformat PDF:HTML -nt $nt -nr 500 -nw 1500 -outdir ${out_dir} --java-mem-size=64G
        tar -zcvf ${out_dir}.tar ${out_dir}
    >>>

    runtime{
		docker:docker
		cluster:cluster_config
		systemDisk:"cloud_ssd 40"
		dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
    }

    output{
        File out_file = "${out_dir}.tar"
    }
}