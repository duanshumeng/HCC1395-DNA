task sequenza {
  
  String sample
  File ref_dir
  String fasta
  File gc
  File tumor_bam
  File tumor_bam_index
  File? normal_bam
  File? normal_bam_index
  File? bed_file
  String docker
  String cluster_config
  String disk_size
  
  command <<<
    set -o pipefail
    set -e
    nt=$(nproc)
    
    seqz=${sample}'.seqz.gz'
    small=${sample}'.small.seqz.gz'
    
    # bam2seqz
                
    sequenza-utils bam2seqz -gc ${gc} --fasta ${ref_dir}/${fasta} -n ${normal_bam} -t ${tumor_bam} -o $seqz -C chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY --parallel 24

    # merge and remove
    zcat ${sample}_*.seqz.gz | awk '{if (NR == 1 || (NR != 1 && $1 != "chromosome")) {print $0}}' | bgzip > $seqz
    tabix -f -s 1 -b 2 -e 2 -S 1 $seqz
    rm ${sample}_*.seqz.gz; rm ${sample}_*.seqz.gz.tbi
    
    # seqz_binning: WES: 50; WGS: 200
    if [ ${bed_file} ]; then
    	sequenza-utils seqz_binning --seqz $seqz -w 50 -o $small
    else
        sequenza-utils seqz_binning --seqz $seqz -w 200 -o $small
    fi
    # analysis in r
    Rscript /home/sequenza/sequenza.r '.' ${sample} 'XY'
  >>>
  
  runtime {
    docker: docker
    cluster: cluster_config
    systemDisk: "cloud_ssd 40"
    dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
  }
  
  output {
    File hrd="${sample}.HRD.txt"
    File alternative_fit="${sample}_alternative_fit.pdf"
    File alternative_solutions="${sample}_alternative_solutions.txt"
    File chromosome_depths="${sample}_chromosome_depths.pdf"
    File chromosome_view="${sample}_chromosome_view.pdf"
    File CN_bars="${sample}_CN_bars.pdf"
    File confints_CP="${sample}_confints_CP.txt"
    File contours_CP="${sample}_contours_CP.pdf"
    File CP_contours="${sample}_CP_contours.pdf"
    File gc_plots="${sample}_gc_plots.pdf"
    File genome_view="${sample}_genome_view.pdf"
    File model_fit="${sample}_model_fit.pdf"
    File mutations="${sample}_mutations.txt"
    File scarHRD_input="${sample}_scarHRD_input.txt"
    File segments="${sample}_segments.txt"
    File sequenza_cp_table="${sample}_sequenza_cp_table.RData"
    File sequenza_extract="${sample}_sequenza_extract.RData"
    File sequenza_log="${sample}_sequenza_log.txt"
    File small_seqz="${sample}.small.seqz.gz"
    File small_seqz_index="${sample}.small.seqz.gz.tbi"
  }
}