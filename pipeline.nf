#!/usr/bin/env nextflow

// Script parameters
params.resources = "/tmp/ref"
params.sample_path = "/tmp/in"
params.samplename = "sample"
params.run = "run"
params.unit = "0"
params.pipeline_name = "nips_viral"
params.pipeline_version = "v1.0.0"
params.threads = "16"
params.mail = "your_mail@mail.com"
params.config_name = "delta_viral.v1.0.0.torch.conf"

params.human_path = "Homo_sapiens/hg19"
params.adapter_path = "adapters/v1"

output_path = "$params.sample_path" + "/" + "$params.pipeline_name" + "/" + "$params.pipeline_version"
def read_numbers = [1, 2]

def viral_list = []
new File(params.resources + "/" + params.config_name).eachLine { line ->
  if (!line.contains("#")) {
      viral_list.add(line.split("\t")[0])
  }
}

log.info """\
Directories
-----------
Resources           : $params.resources
Sample path         : $params.sample_path
Sample Information
------------------
  * Sample          : $params.samplename
  * run             : $params.run
  * unit            : $params.unit
Pipeline
--------
  * name            : $params.pipeline_name
  * version         : $params.pipeline_version
Pipeline Specifics
------------------
  * config name     : $params.config_name
  * Human Ref       : $params.human_path
  * Adapters        : $params.adapter_path
  * # viral genomes : $viral_list.size
"""


process create_output_dir {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  maxRetries 3

  input:
    val sample from params.samplename
    val run from params.run
    val unit from params.unit
    val pipeline_name from params.pipeline_name
    val pipeline_version from params.pipeline_version
    val output_path from output_path

  script:
  """
  echo $output_path
  mkdir -p $output_path
  rm -rf $output_path/*
  """
}

process get_sample_code {

  input:
    val sample from params.samplename
    path sample_path from params.sample_path

  output:
    env sample_code into sample_code

  script:
  """
  if [ `ls -d ${sample_path}/fastq/*${sample}_S*_L0??_R1_001.fastq.gz | wc -l` -ne 0 ];
  then
    sample_code=`ls -d ${sample_path}/fastq/*${sample}_S*_L0??_R1_001.fastq.gz | head -n 1 | sed 's|${sample_path}/fastq/||g' | sed 's|_S.*_L0.._R1_001.fastq.gz||g'`
  elif [ `ls -d ${sample_path}/fastq/*${sample}.R1.fastq.gz | wc -l` -ne 0 ];
  then
    sample_code=`ls -d ${sample_path}/fastq/*${sample}.R1.fastq.gz | head -n 1 | sed 's|${sample_path}/fastq/||g' | sed 's|.R1.fastq.gz||g'`
  else
    echo "fastq file not found"
    exit 1
  fi
  """

}

process combine_fastq {

  input:
    path sample_path from params.sample_path
    val name from sample_code
    
  output:
    path "${name}.R1.fastq.gz" into fastq1
    path "${name}.R2.fastq.gz" into fastq2

  script:
  """
  if [ `ls -d ${sample_path}/fastq/*${name}_S*_L0??_R1_001.fastq.gz | wc -l` -ne 0 ];
  then
    zcat $sample_path/fastq/${name}*_R1_001.fastq.gz >> ${name}.R1.fastq
    zcat $sample_path/fastq/${name}*_R2_001.fastq.gz >> ${name}.R2.fastq
    gzip ${name}.R1.fastq
    gzip ${name}.R2.fastq
  else
    cp $sample_path/fastq/${name}.R1.fastq.gz .
    cp $sample_path/fastq/${name}.R2.fastq.gz .
  fi
  """
}

process get_lib_size_factor {

  input:
    path fastq1 from fastq1

  output:
    env all_reads_factor into all_reads_factor

  script:
  """
  all_reads_factor=`zcat $fastq1 | wc -l | awk '{print \$1/4}' | awk '{print \$1/10000000}'`
  """
}

process read_trimming {
    
  conda '/conda/envs/trimmomatic'
  afterScript 'rm *unpaired*'
  // publishDir "$output_path", mode: 'copy', pattern: '*{_summary.txt}'
  
  input:
    path sample_path from params.sample_path
    val name from sample_code
    val threads from params.threads
    path resources from params.resources
    val adapter_path from params.adapter_path
    path fastq1, stageAs: 'data.R1.fastq.gz' from fastq1
    path fastq2, stageAs: 'data.R2.fastq.gz' from fastq2

  output:
    path "${name}.R1.fastq.gz" into trimmomatic_fastq1
    path "${name}.R2.fastq.gz" into trimmomatic_fastq2
    path "${name}_summary.txt" into trimmomatic_summary

  script:
  """
  trimmomatic PE $fastq1 \
        $fastq2 \
        ${name}.R1.fastq.gz \
        ${name}.unpaired.R1.fastq.gz \
        ${name}.R2.fastq.gz \
        ${name}.unpaired.R2.fastq.gz \
        ILLUMINACLIP:$resources/$adapter_path/adapters.fa:3:6:6:5:keepBothReads \
        SLIDINGWINDOW:4:20 MINLEN:30 \
        -summary ${name}_summary.txt -threads $threads
  """
}

process cleanup1_remove_merged_fastq {
  input:
    path fastq1, stageAs: 'data.R1.fastq.gz' from fastq1
    path fastq2, stageAs: 'data.R2.fastq.gz' from fastq2
    val all_reads_factor from all_reads_factor
    path trim_sum from trimmomatic_summary

  script:
  """
  rm `readlink $fastq1`;
  rm `readlink $fastq2`;
  """

}

process human_mapping {

  conda '/conda/envs/mapping'
  cpus params.threads

  input:
    path resources from params.resources
    val human_path from params.human_path
    path fastq1 from trimmomatic_fastq1
    path fastq2 from trimmomatic_fastq2
    val name from sample_code
    val threads from params.threads
    
  output:
    path "filter_reads.txt" into map_human_filtered_read_names

  script:
  """
  # mapping the data to the Human genome, extract all unmapping fragments, and print the name
  bowtie2 --threads $threads \
    --very-fast-local -L 15 \
		-x $resources/$human_path/genome \
    -1 $fastq1 -2 $fastq2 | samtools view -f 0x4 -f 0x8 | awk '{print \$1}' | uniq > filter_reads.txt
  """
}

process filter_human_reads {

  input:
    path fastq1, stageAs: 'data.R1.fastq.gz' from trimmomatic_fastq1
    path fastq2, stageAs: 'data.R2.fastq.gz' from trimmomatic_fastq2
    each read from read_numbers
    val name from sample_code
    path filter_reads from map_human_filtered_read_names
    
  output:
    path "final.R1.fastq.gz" optional true into rest_fastq1
    path "final.R2.fastq.gz" optional true into rest_fastq2

  script:
  """
  # get all fragments where both reads did not mapp
  zcat "data.R"$read".fastq.gz" | awk '{if(NR%4==1){printf \$1"\\t"}if(NR%4==2){printf \$1"\\t"}if(NR%4==0){print \$1}}' | \
    grep -w -f $filter_reads | \
    awk '{print \$1; print \$2; print "+"; print \$3}' > "final.R"$read".fastq"
  gzip "final.R"$read".fastq"
  rm `readlink data.R"$read".fastq.gz`
  """


}


process fastq_remove_duplicates_and_low_complexity {

  publishDir "$output_path", mode: 'copy', pattern: 'viral*'

  input:
    path rest_fastq1 from rest_fastq1
    path rest_fastq2 from rest_fastq2
    
  output:
    path "viral.R1.fastq.gz" into viral_fastq1
    path "viral.R2.fastq.gz" into viral_fastq2

  script:
  """
  my_dir=`pwd`
  cd /code
  python3.9 -m pipenv run python src/low_complexity_removal.py \$my_dir/$rest_fastq1 \$my_dir/$rest_fastq2 \
        \$my_dir/viral.R1.fastq.gz \$my_dir/viral.R2.fastq.gz -v 
  cd -
  rm `readlink $rest_fastq1`
  rm `readlink $rest_fastq2`
  """
}

process viral_mapping {

  maxForks 4
  cpus params.threads
  conda '/conda/envs/mapping'

  publishDir "$output_path", mode: 'copy', pattern: '*.highq_fragments.txt'

  input:
    path resources from params.resources
    path viral_fastq1 from viral_fastq1
    path viral_fastq2 from viral_fastq2
    val threads from params.threads
    val config_name from params.config_name
    each viral_name from viral_list
    
  output:
    path "*.highq_fragments.txt" into viral_highq_fragments_files_txt

  script:
  """
  # Get the correct path
  viral_path=""
  while read -r line; 
  do
    name=`echo \$line | awk '{print \$1}'`
    if [ \$name == $viral_name ];
    then
      viral_path=`echo \$line | awk '{print \$2}'`
    fi
  done<$resources/$config_name

  # map the data, adjusted the parameters for higher matches
  # end-to-end
  # only select the reads with correct mapping, high quality and expected insert size
  bowtie2 --end-to-end -D 20 -R 3 -N 0 -i S,1,0 -L 10 -p $threads \
      -x $resources/\$viral_path/genome \
      -1 $viral_fastq1 -2 $viral_fastq2 | \
    elprep filter /dev/stdin /dev/stdout --sorting-order coordinate --nr-of-threads $threads \
      --mark-duplicates --mark-optical-duplicates $viral_name".fastq.metrics" | \
    samtools view -F 0x4 -F 0x8 -F 0x100 -F 0x200 -F 0x400 -F 0x800 -f 0x1 -f 0x2 -q 10 /dev/stdin \
      | awk 'function abs(v) {return v < 0 ? -v : v}{if(abs(\$9) < 300){print \$1}}' | sort -k1,1 | uniq > $viral_name".highq_fragments.txt"
  """
}

process viral_combine {

  publishDir "$output_path", mode: 'copy', pattern: '*.viral.json'
  
  input:
    val sample from sample_code
    val run from params.run
    val unit from params.unit
    path "*" from viral_highq_fragments_files_txt.toList()
    path resources from params.resources
    val config_name from params.config_name
    val all_reads_factor from all_reads_factor

  output:
    path "${sample}.viral.json" into viral_information_json

  script:
  """
  my_dir=`pwd`
  cd /code
  python3.9 -m pipenv run python src/combine_viral.py $run $unit $sample \$my_dir/$resources/$config_name $all_reads_factor \$my_dir \$my_dir/$sample".viral.json" -v
  cd -
  """
}

workflow.onComplete {

  myDir = file(output_path)
  myDir.eachFileRecurse { item ->
    item.setPermissions(7,7,7)
  }

}

workflow.onError {

  def date = new Date()
  String datePart = date.format("yyyy-MM-dd")
  String timePart = date.format("HH:mm:ss")
  
  
  def msg = """\
      Pipeline Failed
      ---------------
      Run               : $params.run
      Sample            : $params.samplename
      SampleCode        : $sample_code
      Unit              : $params.unit
      Panel             : $panel_name
      Panel version     : $panel_version


      Pipeline execution summary
      ---------------------------
      Completed at  : ${workflow.complete}
      Duration      : ${workflow.duration}
      Started at    : ${workflow.start}
      Current time  : $datePart $timePart
      Success       : ${workflow.success}
      workDir       : ${workflow.workDir}
      exit status   : ${workflow.exitStatus}


      Pipeline command
      ----------------
      command line  : ${workflow.commandLine}


      Error message and report
      ------------------------
      error message : ${workflow.errorMessage}
      error report  : ${workflow.errorReport}
      

      """
      .stripIndent()

  sendMail(to: params.mail, from: params.mail, subject: 'FAILED: ' + params.pipeline_name + ' of ' + params.samplename, body: msg)
}