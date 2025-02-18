
task fastq_screen {
  
  String sample
  File read1
  File read2
  File screen_ref_dir
  File fastq_screen_conf

  String docker
  String cluster_config
  String disk_size

  command <<<
    set -o pipefail
    set -e
    nt=$(nproc)
    mkdir -p /cromwell_root/tmp
    cp -r ${screen_ref_dir} /cromwell_root/tmp/
    ln -s ${read1} ${sample}_R1.fastq.gz
    ln -s ${read2} ${sample}_R2.fastq.gz
    fastq_screen --aligner bowtie2 --conf ${fastq_screen_conf} --top 100000 --threads $nt ${sample}_R1.fastq.gz
    fastq_screen --aligner bowtie2 --conf ${fastq_screen_conf} --top 100000 --threads $nt ${sample}_R2.fastq.gz
  >>>

  runtime {
    docker:docker
    cluster: cluster_config
    systemDisk: "cloud_ssd 40"
    dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
  }

  output {
    File png1 = "${sample}_R1_screen.png"
    File txt1 = "${sample}_R1_screen.txt"
    File html1 = "${sample}_R1_screen.html"
    File png2 = "${sample}_R2_screen.png"
    File txt2 = "${sample}_R2_screen.txt"
    File html2 = "${sample}_R2_screen.html"
  }
}