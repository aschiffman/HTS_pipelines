#!/bin/bash

# RNAseq_preprocessing.sh - A script to preprocess RNAseq data. Last edited by Allison on 2024/05/06
# Usage: RNAseq_preprocessing.sh list_of_fastq_paths.txt [/path/to/parent_dir] [nproc_per_sample] [strandness]
# list_of_fastq_paths.txt: A file with full path to all R1 fastq.gz or fq.gz samples containing _R1. or _1.
# /path/to/parent_dir: Directory for analysis, default=current directory
# nproc_per_sample: Number of cores per sample, default=1
# strandness: Strandness, default=unset (0=unstranded, 1=forward, 2=reverse)

FASTQ_FILE=$1 # FIRST ARGUMENT, file with list of full path to all samples ending in _R1/fastq or .fastq.gz
PARENT_DIR=${2:-$(pwd)} # SECOND ARGUMENT, directory for analysis, default=current directory
FASTQ_R1=$(<$FASTQ_FILE) # read in list
SAMPLE_NO=$(sed -n '$=' $FASTQ_FILE) # get number of samples from fastq list
NPROC_PER_SAMPLE=${3:-"1"} # THIRD ARGUMENT, number of cores per sample, default=1
STRAND=${4:-"100"} # FOURTH ARGUMENT, strandness, default=unset (0=unstranded, 1=forward, 2=reverse)
# unset STRAND if not 0, 1, or 2
if [[ $STRAND != 0 ]] && [[ $STRAND != 1 ]] && [[ $STRAND != 2 ]]; then
    unset STRAND
fi

# log files
LOG_FILE=$PARENT_DIR"/run_log.txt";LOG_ERR_FILE=$PARENT_DIR"/run_err.txt"
STATS_DIR=$PARENT_DIR/stats
mkdir -p $STATS_DIR
# echo -e "\n\n\n\n\n" | tee -a $LOG_FILE $LOG_ERR_FILE
echo -e "### -------------------------------------------------------------------------------------------- ###" | tee $LOG_FILE $LOG_ERR_FILE
echo -e "(`date`) Starting RNAseq preprocessing... ------------------------------------------- " | tee -a $LOG_FILE $LOG_ERR_FILE

# # For each element in FASTQ_R1, replace instances of _R1 with _R2
FASTQ_R2=$(echo $FASTQ_R1 | sed 's@_R1@_R2@g')
# For each element in FASTQ_R1, replace instances of _1. with _2.
FASTQ_R2=$(echo $FASTQ_R2 | sed 's@_1.@_2.@g')

TOTAL_PROC_NO=$((SAMPLE_NO*NPROC_PER_SAMPLE < $(nproc) ? SAMPLE_NO*NPROC_PER_SAMPLE:$(nproc)-1))
echo -e "Output directory: $PARENT_DIR \n" | tee -a $LOG_FILE
echo -e "FASTQ_FILE=$FASTQ_FILE \n" | tee -a $LOG_FILE
echo -e "Number of samples: $SAMPLE_NO \nNumber of cores per sample: $NPROC_PER_SAMPLE \nNumber of cores to be used: $TOTAL_PROC_NO" | tee -a $LOG_FILE
echo -e "total: $TOTAL_PROC_NO threads wil be used for this analysis\n" | tee -a $LOG_FILE
echo -e "Log file: $LOG_FILE \nError file: $LOG_ERR_FILE\n" | tee -a $LOG_FILE

echo -e "\n\nChecking that R1 and R2 files exist..." | tee -a $LOG_FILE

#  Check that every file in FASTQ_R1 and FASTQ_R2 exists
for file in $FASTQ_R1 $FASTQ_R2; do
    if [ ! -f $file ]; then
        echo -e "File $file does not exist. Remember to specify FULL PATH to file. Exiting..." | tee -a $LOG_FILE $LOG_ERR_FILE
        exit 1
    else
        echo -e "Will process $file" | tee -a $LOG_FILE
    fi

done

# get main genome parameter
export GENOME_DIR_STAR="/home/shared/GRCm39/50bp"
export GENOME_GTF="/home/shared/GRCm39/gencode.vM31.primary_assembly.annotation.gtf"
export GENOME_BED="/home/shared/GRCm39/gencode.vM31_primary_assembly.bed"

cd $PARENT_DIR
echo -e "###########################################" | tee -a $LOG_FILE
echo -e "(`date`) Start trimming and QC process" | tee -a $LOG_FILE
echo -e "###########################################" | tee -a $LOG_FILE
echo -e "(`date`) 1 QC of FASTQ files" | tee -a $LOG_FILE


WORKING_DIR=$PARENT_DIR'/01fastqc'
mkdir -p $WORKING_DIR
parallel -j $SAMPLE_NO --no-notice fastqc -t $NPROC_PER_SAMPLE -outdir $WORKING_DIR {} ::: $FASTQ_R1 $FASTQ_R2 1>>$LOG_FILE 2>>$LOG_ERR_FILE
wait; echo -e "(`date`) Step 1 QC of FASTQ files finished." | tee -a $LOG_FILE

echo -e "###########################################" | tee -a $LOG_FILE
echo -e "(`date`) 2.1 Trimming" | tee -a $LOG_FILE
WORKING_DIR=$PARENT_DIR'/02trim'
mkdir -p $WORKING_DIR

trimfun() {
    FILE1=$1
    FILE2=$2
    echo -e "Trimming $FILE1 and $FILE2" | tee -a $LOG_FILE
    BASE=$(basename $FILE1 | sed 's@\(_R1.*\.fastq\|_R1.*\.fq\.gz\|_R1.*\.fq\|_R1.*\.fastq\.gz\|_1.*\.fastq\|_1.*\.fq\.gz\|_1.*\.fq\|_1.*\.fastq\.gz\)@@g')
    

    # echo -e "BASE: $BASE" | tee -a $LOG_FILE
    # echo -e "Output file: $WORKING_DIR/$BASE.R1.trim.fastq and $WORKING_DIR/$BASE.R2.trim.fastq" | tee -a $LOG_FILE
    WORKING_DIR=$3
    NPROC_PER_SAMPLE=$4
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -j $NPROC_PER_SAMPLE -e 0.1 -O 6 -q 20 -m 35  \
     -o $WORKING_DIR/$BASE.R1.trim.fastq.gz -p $WORKING_DIR/$BASE.R2.trim.fastq.gz $FILE1 $FILE2

}


export -f trimfun
parallel -j $TOTAL_PROC_NO --no-notice trimfun {} $WORKING_DIR $NPROC_PER_SAMPLE ::: $FASTQ_R1 :::+ $FASTQ_R2 1>>$LOG_FILE 2>>$LOG_ERR_FILE
unset -f trimfun
wait; echo -e "(`date`) 2.1 Trimming finished." | tee -a $LOG_FILE

echo -e "###########################################" | tee -a $LOG_FILE
echo -e "(`date`) 2.2 QC of trimmed files" | tee -a $LOG_FILE
cd $WORKING_DIR
WORKING_DIR=$WORKING_DIR'/fastqc'
mkdir -p $WORKING_DIR
ls -1 *.fastq.gz | xargs -n 1 -P $SAMPLE_NO fastqc -t $NPROC_PER_SAMPLE -outdir $WORKING_DIR 1>>$LOG_FILE 2>>$LOG_ERR_FILE
wait;echo -e "(`date`) 2.2 QC trimmed fastqc finished." | tee -a $LOG_FILE

echo -e "(`date`) Step 2 Trimming finished." | tee -a $LOG_FILE

echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 3.1: alignment" | tee -a $LOG_FILE

WORKING_DIR=$PARENT_DIR'/03alignment/'
mkdir -p $WORKING_DIR
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

date | tee -a $LOG_FILE_STEP

alignfun (){
    SAMPLE=$1
    samplename=${SAMPLE%%.*}
    LOG_FILE=$4
    echo -e "(`date`) align $samplename" | tee -a $LOG_FILE
    READ1=$samplename".R1.trim.fastq.gz"
    READ2=$samplename".R2.trim.fastq.gz"
    WORKING_DIR=$2
    NPROC_PER_SAMPLE=$3

    echo -e "If alignment is taking longer than expected, check $WORKING_DIR/${samplename}_Log.out. \n If it is stuck on \"Another job is still loading the genome, sleeping for 1 min\" then contact Allison or Helen to run the following command: \n sudo STAR --genomeDir $GENOME_DIR_STAR --genomeLoad Remove"

    if [ ! -f $READ1 ]; then
        echo -e "Alignment error. Read 1 $READ1 does not exist."  | tee -a $LOG_FILE $LOG_ERR_FILE
        return 1
    else
        echo -e "Read 1 file: $READ1" | tee -a $LOG_FILE
    fi

    if [ ! -f $READ2 ]; then
        echo -e "Alignment error. Read 2 $READ2 does not exist." | tee -a $LOG_FILE $LOG_ERR_FILE
        return 1
    else
        echo -e "Read 2 file: $READ2"  | tee -a $LOG_FILE
    fi

    STAR --genomeDir $GENOME_DIR_STAR \
         --runThreadN $NPROC_PER_SAMPLE \
         --genomeLoad LoadAndRemove \
         --readFilesIn $READ1 $READ2 \
         --outSAMunmapped Within \
         --outSAMtype BAM SortedByCoordinate \
         --limitBAMsortRAM 17179869184 \
         --outFilterType BySJout \
         --outFilterMultimapNmax 20 \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --outFilterMismatchNoverLmax 0.04 \
         --seedSearchStartLmax 30 \
         --readFilesCommand zcat \
         --outFileNamePrefix $WORKING_DIR"/"$samplename"_"
}

ulimit -n 4000

export -f alignfun
echo -e "(`date`) Performing alignment on files in here $PWD" | tee -a $LOG_FILE
TRIM_LIST=$(ls -1 *.R1.trim.fastq.gz)
# Run alignment on MAX 10 samples at once to avoid using too much RAM
parallel --halt soon,fail=1 -j 10 --no-notice alignfun {} $WORKING_DIR $NPROC_PER_SAMPLE $LOG_FILE ::: $TRIM_LIST 1>>$LOG_FILE_STEP 2>>$LOG_ERR_FILE
unset -f alignfun

if [ $? -eq 1 ]; then
    echo "alignfun failed, exiting script."
    exit 1
fi

wait

ulimit -n 1024

echo -e "This version of STAR was used for the analysis" | tee -a $LOG_FILE_STEP

STAR --version | tee -a $LOG_FILE_STEP

echo -e "This version of samtools was used for the analysis" | tee -a $LOG_FILE_STEP

samtools --version | tee -a $LOG_FILE_STEP
wait;echo -e "`date`: Step 3.1 alignment finished!" | tee -a $LOG_FILE

echo -e "############################################################" | tee -a $LOG_FILE

# ------------------------------------------------------------
cd $WORKING_DIR
echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 3.2: index bam file" | tee -a $LOG_FILE

date | tee -a $LOG_FILE_STEP
ls -1 *.bam | xargs -n 1 -P $SAMPLE_NO -I {} \
                    samtools index -@ $NPROC_PER_SAMPLE {} \
                    1>>$LOG_FILE_STEP 2>>$LOG_FILE_STEP
date | tee -a $LOG_FILE_STEP
wait;echo -e "(`date`) Step 3.2 index bam file finished!" | tee -a $LOG_FILE

echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 3.3: Calculate stats" | tee -a $LOG_FILE

alignstatfun (){
    while read data;do
    samplename=${data%.*}
    echo -e "(`date`) calc stats for $samplename" | tee -a $LOG_FILE
    samtools flagstat -@ $NPROC_PER_SAMPLE  $WORKING_DIR"/"$samplename".bam" > $WORKING_DIR"/"$samplename".stat"  &
done
}


ls -1 *.bam | alignstatfun | tee -a $LOG_FILE_STEP 2>> $LOG_FILE_STEP

wait;echo -e "(`date`) Step 3.3 calculate stats finished!" | tee -a $LOG_FILE

# ###################################" | tee -a $LOG_FILE
# echo -e "(`date`) Starting Step 3.4: Quality mapping" | tee -a $LOG_FILE
# STEP='qualimap'
# WORKING_DIR=$WORKING_DIR'/'$STEP
# mkdir -p $WORKING_DIR

# qualmapfun (){
#     int1=1;int2=1;
#     while read data; do
# 	nameStr=$(echo "$data"| cut -f1 -d".")
# 	mkdir -p $WORKING_DIR'/'$nameStr

# 	if [ `echo $int1" % 4" | bc` -eq 0 ]
# 	then
# 	    echo -e "caculating $int1/$SAMPLE_NO samples \n"
# 	    int1=$((int1+int2))
# 	    JAVA_OPTS="-Djava.awt.headless=true" qualimap bamqc -nt 2 -c  -bam "$data"  -outdir $WORKING_DIR/$nameStr --java-mem-size=4G
# 	else
# 	    echo -e "caculating $int1/$SAMPLE_NO samples \n"
# 	    int1=$((int1+int2))
# 	    JAVA_OPTS="-Djava.awt.headless=true" qualimap bamqc -nt 2 -c  -bam "$data"  -outdir $WORKING_DIR/$nameStr --java-mem-size=4G &
# 	fi
#     done
# }

# ls *.bam | qualmapfun 1>>$LOG_ERR_FILE 2>>$LOG_FILE_STEP
# echo -e "This version of qualimap was used for the analysis" | tee -a $LOG_FILE_STEP
# qualimap --version >> $LOG_FILE_STEP 2>> $LOG_FILE_STEP

# wait;echo -e "(`date`) Step 3.4 quality mapping finished!" | tee -a $LOG_FILE

echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 4.1: filtering the aligned bam files" | tee -a $LOG_FILE
STEP='04filter'
WORKING_DIR=$PARENT_DIR'/'$STEP
mkdir -p $WORKING_DIR
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

filter (){
    while read data;do
    samplename=${data%.*.*.*}
    echo -e "(`date`) Filter $samplename for quality >30" | tee -a $LOG_FILE

    samtools view -b -F 2820 -q 30 -@ $NPROC_PER_SAMPLE ${data} > $WORKING_DIR"/"$samplename".filter.bam" &
done
}


date | tee -a $LOG_FILE_STEP
ls -1 *.bam | filter 2>&1 | tee -a $LOG_FILE_STEP
wait
date | tee -a $LOG_FILE_STEP

wait; echo -e "(`date`)  Step 4.1 filter bam finished" | tee -a $LOG_FILE $LOG_FILE_STEP

echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 4.2: indexing the filtered bam" | tee -a $LOG_FILE

cd $WORKING_DIR
date | tee -a $LOG_FILE_STEP
ls *.bam | parallel --no-notice --progress -j $SAMPLE_NO samtools index -@ $NPROC_PER_SAMPLE {}
                    1>>$LOG_FILE_STEP 2>>$LOG_ERR_FILE


date | tee -a $LOG_FILE_STEP

wait; echo -e "(`date`)  finished Step 4.2: indexing the filtered bam" | tee -a $LOG_FILE $LOG_FILE_STEP

echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 4.3: Calculate stats" | tee -a $LOG_FILE

filterstatfun (){
    while read data;do
    samplename=${data%.*}
    echo -e "(`date`) calc stats for $samplename" | tee -a $LOG_FILE
    samtools flagstat -@ $NPROC_PER_SAMPLE  $WORKING_DIR"/"$samplename".bam" > $WORKING_DIR"/"$samplename".stat"  &
done
}

ls -1 *.bam | filterstatfun | tee -a $LOG_FILE_STEP 2>> $LOG_FILE_STEP

wait;echo -e "(`date`) Step 4.3 calculate stats finished!" | tee -a $LOG_FILE

echo -e "\n############################################################" | tee -a $LOG_FILE
echo "$(date) Starting step 4.4: Library strandness" | tee -a $LOG_FILE
echo "############################################################" | tee -a $LOG_FILE

if [[ ! $STRAND ]] # if strand not set; auto detect
then
    echo "Strandness is not set, it will be auto-detected" | tee -a $LOG_FILE
    STEP='Infer_strandness'
    INFER_OUTDIR=${WORKING_DIR}/${STEP}

    mkdir -p ${INFER_OUTDIR}

    infer_strandness (){
        FILE=$1
        OUT_DIR=$2

        BASENAME=$( basename -s ".filtered.bam" $FILE )

        echo "Infer strandness of file: $FILE"

        infer_experiment.py -r $GENOME_BED -i $FILE -q 30 > ${OUT_DIR}/${BASENAME}_infer_experiment.log
    }

    export -f infer_strandness
    parallel -j $(echo `expr $TOTAL_PROC_NO \* 2`) --no-notice infer_strandness {} $INFER_OUTDIR ::: $(ls *.bam ) 1>>$LOG_FILE 2>>$LOG_ERR_FILE
    unset -f infer_strandness

    wait

    # check result for all samples
    UNSTRANDED=0
    FORWARD=0
    REVERSE=0
    UNDETERMINE=0

    for INFER_FILE in $(ls ${INFER_OUTDIR}/*_infer_experiment.log )
    do
        sf=$(cat $INFER_FILE | grep '1++,1--,2+-,2-+' | cut -d":" -f 2) # get forward percentage
        sr=$(cat $INFER_FILE | grep '1+-,1-+,2++,2--' | cut -d":" -f 2) # get reverse percentage

        if (( $(echo "$sf >= 0.8" | bc -l) )); # if forward percentage >= 0.8, sample is forwardly stranded
        then
            FORWARD=$(( $FORWARD + 1))
        elif (( $(echo "$sr >= 0.8" | bc -l) )); # if reverse percentage >= 0.8, sample is reversely stranded
        then
            REVERSE=$(( $REVERSE + 1 ))
        elif (( $(echo " 0.4 <= $sf && $sf <= 0.6 && 0.4 <= $sr && $sr <= 0.6" | bc -l) )) # if forward and reverse percentage between 0.4-0.6 , sample is not stranded
        then
            UNSTRANDED=$(( $UNSTRANDED + 1 ))
        else
            UNDETERMINE=$(( $UNDETERMINE + 1 ))
        fi
    done

    echo "Infered strandness" | tee -a $LOG_FILE
    echo "$UNSTRANDED sample(s) are not strand-specific" | tee -a $LOG_FILE
    echo "$FORWARD sample(s) are forwardly stranded" | tee -a $LOG_FILE
    echo "$REVERSE sample(s) are reversely stranded" | tee -a $LOG_FILE
    echo "$UNDETERMINE sample(s) strandness could not be determined" | tee -a $LOG_FILE

    # determine strandness
    if [[ $FORWARD == 0 ]] && [[ $REVERSE == 0 ]] && [[ $UNDETERMINE == 0 ]] && [[ $UNSTRANDED == $(ls ${INFER_OUTDIR}/*_infer_experiment.log | wc -l) ]] # all non strand specific
    then
        STRAND=0 # unstranded
        echo "LIBRARY IS INFERED TO BE UNSTRANDED" | tee -a $LOG_FILE
    elif [[ $FORWARD == 0 ]] && [[ $UNSTRANDED == 0 ]] && [[ $UNDETERMINE == 0 ]] && [[ $REVERSE == $(ls ${INFER_OUTDIR}/*_infer_experiment.log | wc -l) ]] # all reverse strand specific
    then
        STRAND=2 # reversely stranded
        echo "LIBRARY IS INFERED TO BE REVERSELY STRANDED" | tee -a $LOG_FILE
    elif [[ $REVERSE == 0 ]] && [[ $UNSTRANDED == 0 ]] && [[ $UNDETERMINE == 0 ]] && [[ $FORWARD == $(ls ${INFER_OUTDIR}/*infer_experiment.log | wc -l) ]] # all forward strand specific
    then
        STRAND=1 # fowardly stranded
        echo "LIBRARY IS INFERED TO BE FORWARDLY STRANDED" | tee -a $LOG_FILE
    elif [[ $UNDETERMINE != 0 ]] # some sample not determine
    then
        echo "SOME SAMPLES STRANDNESS COULD NOT BE DETERMINED!!" | tee -a $LOG_FILE
        echo "ARBORTING PIPELINE" | tee -a $LOG_FILE
        echo "Check what type of library was used and rerun by manually giving the correct strandness information:" | tee -a $LOG_FILE
        echo "Where STRAND is replaced by 0 for unstranded library, 1 for forward stranded library, 2 for reversely stranded library" | tee -a $LOG_FILE
        exit
    else
        echo "SAMPLES HAVE INCONSISTENT LIBRARY STRANDNESS!!" | tee -a $LOG_FILE
        echo "ARBORTING PIPELINE" | tee -a $LOG_FILE
        echo "Check what type of library was used and rerun by manually giving the correct strandness information:" | tee -a $LOG_FILE
        echo "Where STRAND is replaced by 0 for unstranded library, 1 for forward stranded library, 2 for reversely stranded library" | tee -a $LOG_FILE
        exit
    fi
    unset UNSTRANDED FORWARD REVERSE UNDETERMINE
else
    echo "Strandness was set manually to $STRAND" | tee -a $LOG_FILE
fi

export STRAND

echo -e "############################################################" | tee -a $LOG_FILE
echo -e " (`date`) Starting Step 4.4: qualimap QC the mapping" | tee -a $LOG_FILE

STEP='qualimap'
WORKING_DIR=$WORKING_DIR'/'$STEP
mkdir -p $WORKING_DIR

qualmapfun (){
    int1=1;int2=1;
    while read data; do
	nameStr=$(echo "$data"| cut -f1 -d".")
	mkdir -p $WORKING_DIR'/'$nameStr

	if [ `echo $int1" % 4" | bc` -eq 0 ]
	then
	    echo -e "caculating $int1/$SAMPLE_NO samples \n"
	    int1=$((int1+int2))
	    JAVA_OPTS="-Djava.awt.headless=true" qualimap bamqc -nt 2 -c  -bam "$data"  -outdir $WORKING_DIR/$nameStr --java-mem-size=4G
	else
	    echo -e "caculating $int1/$SAMPLE_NO samples \n"
	    int1=$((int1+int2))
	    JAVA_OPTS="-Djava.awt.headless=true" qualimap bamqc -nt 2 -c  -bam "$data"  -outdir $WORKING_DIR/$nameStr --java-mem-size=4G &
	fi
    done
}

# ls *.bam | qualmapfun 1>>$LOG_ERR_FILE 2>>$LOG_FILE_STEP
export -f qualmapfun
parallel -j $SAMPLE_NO --no-notice qualmapfun {} $WORKING_DIR ::: $(ls *.bam) 1>>$LOG_ERR_FILE 2>>$LOG_FILE
unset -f qualmapfun
echo -e "This version of qualimap was used for the analysis" | tee -a $LOG_FILE_STEP
qualimap -h | grep "QualiMap" >> $LOG_FILE_STEP 2>> $LOG_FILE_STEP

wait;echo -e "(`date`) Step 4.4 quality mapping finished!" | tee -a $LOG_FILE

echo -e "(`date`)Starting Step 5: generate the counts file " | tee -a $LOG_FILE


STEP='05counts'
WORKING_DIR=$PARENT_DIR'/'$STEP; mkdir -p $WORKING_DIR;

featureCounts -p -T $TOTAL_PROC_NO -s $STRAND -t exon -g gene_id -a $GENOME_GTF -o $WORKING_DIR/gene_counts.tsv *.bam 1>>$LOG_ERR_FILE 2>>$LOG_FILE

wait;echo -e "(`date`) Step 5 finished!" | tee -a $LOG_FILE


echo -e "############################################################" | tee -a $LOG_FILE

echo -e "(`date`)Starting Step 6: make tracks using bamCoverage" | tee -a $LOG_FILE

STEP='06tracks'
WORKING_DIR=$PARENT_DIR'/'$STEP; mkdir -p $WORKING_DIR;

bcfun (){
	BW=$(basename $1 | sed "s|\.bam|.bw|")
	OUT_BW=$2/$BW	
	echo -e "(`date`) generate tracks for $1"
	bamCoverage -p $3 -b $1 -o $OUT_BW
	echo -e "(`date`) created $BW"

}

export -f bcfun
BAM_LIST=$(ls -1 *.bam)
parallel -j $TOTAL_PROC_NO --no-notice bcfun {} $WORKING_DIR $NPROC_PER_SAMPLE ::: $BAM_LIST 1>>$LOG_FILE_STEP 2>>$LOG_ERR_FILE
unset -f bcfun

wait

wait;echo -e "(`date`) Step 6 finished!" | tee -a $LOG_FILE

echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`)Starting Step 7: mutliqc" | tee -a $LOG_FILE


cd $PARENT_DIR
rm -rf 07multiqc
mkdir -p 07multiqc

multiqc . -o ./07multiqc 1>>$LOG_ERR_FILE 2>>$LOG_FILE
wait

# Save stats
echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Saving stats" | tee -a $LOG_FILE

# Save mapping stats
savemapfun(){
    while read data;do
        samplename=${data%.*}
        shortname=$(echo $samplename | rev | cut -d'_' -f2- | rev)
        PREMAP=$(cat $samplename.stat | grep "total" | cut -d' ' -f1)
        echo -e "Pre-mapping reads:\t$PREMAP" > $STATS_DIR/$shortname".stat"

        POSTMAP=$(cat $samplename.stat | grep "mapped (" | cut -d' ' -f1)
        PCTMAP=$(cat $samplename.stat | grep "mapped (" | cut -d' ' -f5)
        echo -e "Mapped reads:\t$POSTMAP $PCTMAP)" >> $STATS_DIR/$shortname".stat"

        # short sample name is everything before _Aligned
        UNIQMAP=$(cat $shortname"_Log.final.out" | grep "Uniquely mapped reads number" | cut -f2)
        UNIQPCT=$(cat $shortname"_Log.final.out" | grep "Uniquely mapped reads %" | cut -f2)
        echo -e "Uniquely mapped reads:\t$UNIQMAP ($UNIQPCT)" >> $STATS_DIR/$shortname".stat"
done
}

cd $PARENT_DIR/03alignment/
echo -e "(`date`) Saving mapping stats" | tee -a $LOG_FILE
ls -1 *.bam | savemapfun 1>> $LOG_FILE 2>> $LOG_ERR_FILE | tee -a $LOG_FILE

# Save library size
cd $PARENT_DIR"/04filter/"
echo -e "Sample\t library" > $STATS_DIR/library_size.tsv   
savefun (){
    while read data;do
        # sample name is everything before _Aligned
        samplename=$(echo $data | rev | cut -d'_' -f2- | rev)
        fullname=$(echo $data | rev | cut -d'.' -f2- | rev)

        # Save library size from stats file
        total=$(cat $fullname.stat | grep "in total" | cut -d' ' -f1)
        echo -e "$samplename\t$total" >> $STATS_DIR/library_size.tsv

        # Save library size to stats file
        echo -e "Mapped and filtered:\t $total" >> $STATS_DIR/$samplename".stat"
done
}
echo -e "(`date`) Saving library size" | tee -a $LOG_FILE
ls -1 *.bam | savefun 1>> $LOG_FILE 2>> $LOG_ERR_FILE | tee -a $LOG_FILE


echo -e "############################################################" | tee -a $LOG_FILE

unset GENOME_DIR_STAR GENOME_GTF GENOME_BED NPROC_PER_SAMPLE SAMPLE_NOR

echo -e "############################################################"
echo -e "This script has finished running"
