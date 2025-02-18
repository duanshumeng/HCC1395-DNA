task fastp {
        
    # I/O options
    File in1
    File in2
    String sample_id

    Boolean? phred64 = false 
    Boolean? fix_mgi_id = false

    String? adapter_sequence
    String? adapter_sequence_r2

    Int? reads_to_process # specify how many reads/pairs to be processed. Default 0 means process all reads.

    # reporting options
    String json = sample_id+"fastp.json"
    String html = sample_id+"fastp.html"
    String report_title = "\'fastp report\'"

    # excute env
    String docker
    String cluster_config
    String disk_size


    String out1_name = sample_id+'clean_1.fastq.gz'
    String out2_name = sample_id+'clean_2.fastq.gz'

    command <<<
        set -o pipefail
        set -e
        nt=$(nproc)
                  
        # basic command
        /opt/conda/bin/fastp \
        --in1 ${in1} \
        --in2 ${in2} \
        --out1 ${out1_name} \
        --out2 ${out2_name} \
        --json ${json} \
        --html ${html} \
        --report_title ${report_title} \
        --thread $nt \
        --length_required 30
    >>>

    runtime {
		docker:docker
		cluster:cluster_config
		systemDisk:"cloud_ssd 40"
		dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
    }

    output {
        File out1 = out1_name
        File out2 = out2_name
        File json_report = json
        File html_report = html
    }

}