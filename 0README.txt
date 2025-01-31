------------------------------------------------------------------------------------------------------------------------------------------------------
NGS pre-processing pipelines written by Allison
Note: the pipelines align to the mouse genome unless otherwise specified.
PAIRED END
- ATAC-seq pre-processing mouse (PAIRED END, updated June 2024)
- ChIP-seq pre-processing mouse (PAIRED END, updated June 2024)
- RNA-seq pre-processing mouse (PAIRED END, updated June 2024)
- RNA-seq pre-processing human (PAIRED END, updated June 2024)
- CUT&Tag pre-processing mouse (updated June 2024)
- CUT&Tag pre-processing human (updated June 2024)
SINGLE END
- ChIP-seq pre-processing (updated 2022/11/16 by ABS)
- ChIP-seq pre-processing w/ human genome (updated 2023/02/14 by ABS)
- ATAC-seq pre-processing (updated 2022/11/16 by ABS)
- CUT&Tag pre-processing mouse (updated Dec 2022)
- CUT&Tag pre-processing human (updated June 2024)

----------------------------------------------------------------------------------------------------------------
Before running any pipelines...
You will need to initialize conda. You can do that by adding the following block of code to your ~/.bashrc file.
----------------------------------------------------------------------------------------------------------------
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/usr/local/bin/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/usr/local/bin/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/usr/local/bin/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/usr/local/bin/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
----------------------------------------------------------------------------------------------------------------

Instructions for running the pipelines are below.

---------------------------------------------------------------------------------------------------------
PAIRED-END
---------------------------------------------------------------------------------------------------------
-------------------------
RNA-seq pre-processing PAIRED END
-------------------------
Steps:

1. fastqc on fastq files
2. adapter trimming with adapter AGATCGGAAGAGC using cutadapt
    1. fastqc
3. Alignment to GRCh38 genome using STAR
    1. qualimap
4. Filter
    1. -F 2828 remove read unmapped or with mate unmapped, not primary alignment reads, reads failed platform/vendor quality check, supplementary alignment,
    2. -f 2 keep read mapped in proper pair
    3. -q 30 remove low quality alignment
5. Infer strandedness if not set (infer_experiment.py)
6. Make wig IGV tracks (note: you may want to further process these by converting to bw)
7. multiqc

How to run:
-------------------------
/home/shared/pipelines/RNAseq_preprocessing_PE_human.sh takes 1-4 arguments:
1.fastq_list.txt: List of full paths to all R1 fastq files. 
Example command for making this list: ls -1 /home/shared/pipelines/example_PE/data/*R1.fastq > fastq_list.txt

2. example_directory: parent directory where your 01fastqc, 02trim, etc. directories will be stored.
/home/shared/pipelines/example_directory/

3. n: optional third argument for # of threads used for EACH sample. Default=1. The script already runs samples in parallel. If you want multiple threads for each alignment, etc. you can set this number to 2+.

4. strandedness: optional fourth argument for strandedness (0=unstranded, 1=forward, 2=reverse). The default value is unset. If this value is unset or anything other than 0, 1, or 2 is set, strandedness will be inferred.

Example command line for running script:
/home/shared/pipelines/RNAseq_preprocessing_PE_human.sh /path/to/fastq_list.txt [path/to/output_directory/] [num_cores_per_sample] [strandedness (0=unstranded, 1=forward, 2=reverse, default=unset)]

-------------------------
ATAC-seq pre-processing PAIRED END
-------------------------
Steps:
1. fastqc on fastq files
2. adapter trimming with adapter AGATCGGAAGAGC using cutadapt
    + fastqc
3. Alignment to mm10 genome using bowtie2
    + qualimap
4. Filter for reads with quality score > 30
    * new step not in previous script
5. Remove read duplicates with Picard
6. Make bigwig IGV tracks with reads NO EXTENSION
7. Quality control with spp
    + qualimap
Other: Saves stats file with # of reads passing each threshold. Also in stats directory, makes fingerprint plot.
NRF plot is not currently working.

How to run:
-------------------------
/home/shared/pipelines/ATACseq_preprocessing.sh takes two arguments:
1. ATAC_fastq_list_example.list: file with list of the full paths to all R1 fastq files. Make sure that all files in this list include "_R1" as part of the file name, and that the R2 file has the same name except that you replace "_R1" with "_R2". If you get any files back from TCGB that don't match this pattern, let Allison know ASAP.
Example command for making this list: ls -1 /home/shared/pipelines/example_PE/data/*R1.fastq > ATAC_fastq_list_example.list

2. example_directory: parent directory where your 01fastqc, 02trim, etc. directories will be stored.
/home/shared/pipelines/example_directory/

3. n: optional third argument for # of threads used for EACH sample. Default=1. The script already runs samples in parallel. If you want multiple threads for each alignment, etc. you can set this number to 2+.

Example command line for running script:
/home/shared/pipelines/ATACseq_preprocessing.sh /home/shared/pipelines/example_directory/ATAC_fastq_list_example.list /home/shared/pipelines/example_directory/ 2

-------------------------
ChIP-seq pre-processing PAIRED END
-------------------------
Steps:
1. fastqc on fastq files
2. adapter trimming with adapter AGATCGGAAGAGC using cutadapt
    + fastqc
3. Alignment to mm10 genome using bowtie2
    + qualimap
4. Filter for reads with quality score > 30
    * new step not in previous script
5. Remove read duplicates with Picard
6. Make bigwig IGV tracks with reads extended 200 bp
7. Quality control with spp
    + qualimap
Other: Saves stats file with # of reads passing each threshold. Also in stats directory, makes fingerprint plot.
NRF plot is not currently working.

How to run:
-------------------------
/home/shared/pipelines/ChIPseq_preprocessing.sh takes two arguments:
1. ChIP_fastq_list_example.list: file with list of the full paths to all R1 fastq files. file with list of the full paths to all R1 fastq files. Make sure that all files in this list include "_R1" as part of the file name, and that the R2 file has the same name except that you replace "_R1" with "_R2". If you get any files back from TCGB that don't match this pattern, let Allison know ASAP.
Example command for making this list: ls -1 /home/shared/pipelines/example_PE/data/*R1.fastq > ChIP_fastq_list_example.list

2. example_directory: parent directory where your 01fastqc, 02trim, etc. directories will be stored.
/home/shared/pipelines/example_directory/

3. n: optional third argument for # of threads used for EACH sample. Default=1. The script already runs samples in parallel. If you want multiple threads for each alignment, etc. you can set this number to 2+.

Example command line for running script:
/home/shared/pipelines/ChIPseq_preprocessing.sh /home/shared/pipelines/example_directory/ChIP_fastq_list_example.list /home/shared/pipelines/example_directory/ 2

-------------------------
CUT&Tag pre-processing PAIRED-end
-------------------------
Steps:
1. fastqc on fastq files
2. adapter trimming with adapter AGATCGGAAGAGC using cutadapt
    + fastqc
3. Alignment to mm10 genome using bowtie2
    + qualimap
4. Filter for reads with quality score > 30
    * new step not in previous script
5. Remove read duplicates with Picard
6. Make bigwig IGV tracks with reads extended 200 bp
7. Quality control with spp
    + qualimap
Other: Saves stats file with # of reads passing each threshold. Also in stats directory, makes fingerprint plot.
NRF plot is not currently working.

How to run:
-------------------------
/home/shared/pipelines/CUTandTag_preprocessing_PE_mouse.sh takes two arguments:
1. CT_fastq_list_example.list: file with list of the full paths to all R1 fastq files. file with list of the full paths to all R1 fastq files. Make sure that all files in this list include "_R1" as part of the file name, and that the R2 file has the same name except that you replace "_R1" with "_R2". If you get any files back from TCGB that don't match this pattern, let Allison know ASAP.
Example command for making this list: ls -1 /home/shared/pipelines/example_PE/data/*R1.fastq > CT_fastq_list_example.list

2. example_directory: parent directory where your 01fastqc, 02trim, etc. directories will be stored.
/home/shared/pipelines/example_directory/

3. n: optional third argument for # of threads used for EACH sample. Default=1. The script already runs samples in parallel. If you want multiple threads for each alignment, etc. you can set this number to 2+.

Example command line for running script:
/home/shared/pipelines/CUTandTag_preprocessing_PE_mouse.sh /home/shared/pipelines/example_directory/CT_fastq_list_example.list /home/shared/pipelines/example_directory/ 2

---------------------------------------------------------------------------------------------------------
SINGLE-END
---------------------------------------------------------------------------------------------------------
-------------------------
ChIP-seq pre-processing SINGLE-END
-------------------------
Steps:
1. fastqc on fastq files
2. adapter trimming with adapter AGATCGGAAGAGC using cutadapt
    + fastqc
3. Alignment to mm10 genome using bowtie2
    + qualimap
4. Filter for reads with quality score > 30
    * new step not in previous script
5. Remove read duplicates with Picard
6. Make bigwig IGV tracks with reads extended 200 bp
7. Quality control with spp
    + qualimap
Other: Saves stats file with # of reads passing each threshold. Also in stats directory, makes fingerprint plot.
NRF plot is not currently working.

How to run:
-------------------------
/home/shared/pipelines/ChIPseq_preprocessing.sh takes two arguments:
1. ChIP_fastq_list_example.list: file with list of the full paths to all fastq files. They don't have to be in the same directory, and can be on biggie, etc.
ChIP_fastq_list_example.list
    /home/allison/storage/exp_RelA_TAD/rep2/data/ChIP_1_S1_R1_001.fastq.gz
    /home/allison/storage/exp_RelA_TAD/rep2/data/ChIP_2_S2_R1_001.fastq.gz
2. example_directory: parent directory where your 01fastqc, 02trim, etc. directories will be stored.
/home/shared/pipelines/example_directory/
3. n: optional third argument for # of threads used for EACH sample. Default=1. The script already runs samples in parallel. If you want multiple threads for each alignment, etc. you can set this number to 2+.

Example command line for running script:
/home/shared/pipelines/ChIPseq_preprocessing.sh /home/shared/pipelines/example_directory/ChIP_fastq_list_example.list /home/shared/pipelines/example_directory/ 2

/home/shared/pipelines/ChIPseq_preprocessing_hg38.sh is exactly the same, but aligns to the human genome instead of the mouse genome.


-------------------------
ATAC-seq pre-processing SINGLE-END
-------------------------
Steps:
1. fastqc on fastq files
2. adapter trimming with adapter AGATCGGAAGAGC using cutadapt
    + fastqc
3. Alignment to mm10 genome using bowtie2
    + qualimap
4. Filter for reads with quality score > 30
    * new step not in previous script
5. Remove read duplicates with Picard
6. Make bigwig IGV tracks with reads NO EXTENSION
7. Quality control with spp
    + qualimap
Other: Saves stats file with # of reads passing each threshold. Also in stats directory, makes fingerprint plot.
NRF plot is not currently working.

How to run:
-------------------------
/home/shared/pipelines/ATACseq_preprocessing.sh takes two arguments:
1. ATAC_fastq_list_example.list: file with list of the full paths to all fastq files. They don't have to be in the same directory, and can be on biggie, etc.
ATAC_fastq_list_example.list
2. example_directory: parent directory where your 01fastqc, 02trim, etc. directories will be stored.

/home/shared/pipelines/example_directory/
3. n: optional third argument for # of threads used for EACH sample. Default=1. The script already runs samples in parallel. If you want multiple threads for each alignment, etc. you can set this number to 2+.

Example command line for running script:
/home/shared/pipelines/ATACseq_preprocessing.sh /home/shared/pipelines/example_directory/ATAC_fastq_list_example.list /home/shared/pipelines/example_directory/ 2





