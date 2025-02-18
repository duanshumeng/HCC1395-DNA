task cnvkit{
    String sample_id
    File ref_dir
    String? fasta
    File ref_flat
    File? bed_file
    File hrd
    File tumor_bam
    File tumor_bam_index
    File? normal_bam
    File? normal_bam_index
    String docker
    String cluster_config
    String disk_size
    String sample=basename(tumor_bam,".bam")

    command <<<

       set -o pipefail
       set -e
       nt=$(nproc)
        
       echo ${sample}
                  
                  
       Ploidy=`awk -F'\t' '{print $7}' ${hrd} | sed -n '2p'`
                  
       echo $Ploidy
       
       center=`awk -v vv=$Ploidy 'BEGIN {print log(2/vv)/log(2)}'`
       
       echo $center 
                  
       #WES
       if [ ${bed_file} ]; then
           echo "WES"
           if [ ${normal_bam} ]; then
                   echo "WES with normal"
                   /root/miniconda2/bin/cnvkit.py target ${bed_file} --annotate ${ref_flat} --split --short-names -o my_baits.bed
                   /root/miniconda2/bin/cnvkit.py batch ${tumor_bam} \
                      --normal ${normal_bam} \
                      --targets my_baits.bed \
                      --fasta ${ref_dir} \
                      --annotate ${ref_flat} -p $nt \
                      --drop-low-coverage \
                      --output-dir ${sample}.reference.cnn
           else
                   echo "WES no normal"
                   /root/miniconda2/bin/cnvkit.py access ${ref_dir} -o access.bed
                   # Prepare the target bed
                  /root/miniconda2/bin/cnvkit.py target ${bed_file} --annotate ${ref_flat} --split --short-names -o my_baits.bed

                   /root/miniconda2/bin/cnvkit.py autobin ${tumor_bam} -t my_baits.bed -g access.bed

                   /root/miniconda2/bin/cnvkit.py coverage ${tumor_bam} my_baits.target.bed -o ${sample}.T.targetcoverage.cnn
                   /root/miniconda2/bin/cnvkit.py coverage ${tumor_bam} my_baits.antitarget.bed -o ${sample}.T.antitargetcoverage.cnn

                   /root/miniconda2/bin/cnvkit.py reference -o ${sample}.reference.cnn/reference.cnn -f ${ref_dir} -t my_baits.target.bed -a my_baits.antitarget.bed
           fi
    
       #WGS         
       else
           echo "WGS"
           if [ ${normal_bam} ]; then
                   echo "WGS with normal"
                   /root/miniconda2/bin/cnvkit.py batch ${tumor_bam} \
                      --normal ${normal_bam} \
                      --method wgs \
                      --fasta ${ref_dir} \
                      --annotate ${ref_flat} -p $nt \
                      --drop-low-coverage \
                      --output-dir ${sample}.reference.cnn
           else
                   echo "WGS no normal"
                   /root/miniconda2/bin/cnvkit.py access ${ref_dir} -o access.bed
                   # Prepare the target bed
                   #/root/miniconda2/bin/cnvkit.py --annotate ${ref_flat} --split --short-names -o my_baits.bed

                   /root/miniconda2/bin/cnvkit.py autobin ${tumor_bam} --method wgs -g access.bed

                   /root/miniconda2/bin/cnvkit.py coverage ${tumor_bam} ${sample}.target.bed -o ${sample}.T.targetcoverage.cnn
                   /root/miniconda2/bin/cnvkit.py coverage ${tumor_bam} ${sample}.antitarget.bed -o ${sample}.T.antitargetcoverage.cnn

                   /root/miniconda2/bin/cnvkit.py reference -o ${sample}.reference.cnn/reference.cnn -f ${ref_dir} -t ${sample}.target.bed -a ${sample}.antitarget.bed

           fi
       fi
                  
       ls ./
       
        
       /root/miniconda2/bin/cnvkit.py batch  ${tumor_bam} \
                  -r ${sample}.reference.cnn/reference.cnn \
                  --output-dir ${sample}.cns \
                  -p $nt
                  
       ls ./
       /root/miniconda2/bin/cnvkit.py call ${sample}.cns/${sample}.cns --center-at $center \
                  -o ${sample}.call.cns
                  
                

        ls ./

        tar cvf ${sample}.tar ${sample}*

    >>>

    runtime{
        docker:docker
        cluster:cluster_config
        systemDisk:"cloud_ssd 40"
        dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
        timeout:259200
    }

    output{
        File out_file = "${sample}.tar"
        #File cnv_bed = "${sample}.ratio_cnv.call.filter.bed"
    }

}