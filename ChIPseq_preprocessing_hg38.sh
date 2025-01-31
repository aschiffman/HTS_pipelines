#!/bin/bash

FASTQ_FILE=$1 # first argument, file with list of full path to all samples
PARENT_DIR=$2 # second argument, directory for analysis
FASTQ_LIST=$(<$FASTQ_FILE) # read in list
SAMPLE_NO=$(sed -n '$=' $FASTQ_FILE) # get number of samples from fastq list
NPROC_PER_SAMPLE=${3:-"1"} # third argument, number of cores per sample, default=1

# log files
LOG_FILE=$PARENT_DIR"/run_log.txt";LOG_ERR_FILE=$PARENT_DIR"/run_err.txt"
STATS_DIR=$PARENT_DIR/stats
mkdir -p $STATS_DIR
echo -e "\n\n\n\n\n" | tee -a $LOG_FILE $LOG_ERR_FILE
echo -e "### -------------------------------------------------------------------------------------------- ###" | tee -a $LOG_FILE $LOG_ERR_FILE
echo -e "(`date`) Starting preprocessing... ------------------------------------------- " | tee -a $LOG_FILE $LOG_ERR_FILE

TOTAL_PROC_NO=$((SAMPLE_NO*NPROC_PER_SAMPLE < $(nproc) ? SAMPLE_NO*NPROC_PER_SAMPLE:$(nproc)-1))
echo -e "FASTQ_FILE=$FASTQ_FILE \n\nsample fastq list:\n$FASTQ_LIST \n\nPARENT_DIR=$PARENT_DIR \nnumber of samples=$SAMPLE_NO, processors per sample = $NPROC_PER_SAMPLE" | tee -a $LOG_FILE
echo -e "total: $TOTAL_PROC_NO threads wil be used for this analysis\n" | tee -a $LOG_FILE
echo -e "Log file: $LOG_FILE \t Error file: $LOG_ERR_FILE\n" | tee -a $LOG_FILE


#debug
#WORKING_DIR=$PARENT_DIR'/05deDupped/'
#LOG_FILE_STEP=$WORKING_DIR"/log.txt"
#cd $PARENT_DIR/05deDupped

cd $PARENT_DIR
echo -e "###########################################" | tee -a $LOG_FILE
echo -e "(`date`) Start trimming and QC process" | tee -a $LOG_FILE
echo -e "###########################################" | tee -a $LOG_FILE
echo -e "(`date`) 1 QC of FASTQ files" | tee -a $LOG_FILE


WORKING_DIR=$PARENT_DIR'/01fastqc'
mkdir -p $WORKING_DIR
parallel -j $SAMPLE_NO --no-notice fastqc -t $NPROC_PER_SAMPLE -outdir $WORKING_DIR {} ::: $FASTQ_LIST 1>>$LOG_FILE 2>>$LOG_ERR_FILE
wait; echo -e "(`date`) 1 QC of FASTQ files finished." | tee -a $LOG_FILE
 

echo -e "###########################################" | tee -a $LOG_FILE
echo -e "(`date`) 2 Trimming" | tee -a $LOG_FILE
WORKING_DIR=$PARENT_DIR'/02trim'
mkdir -p $WORKING_DIR

trimfun() {
    FILE=$1
    echo -e "Trimming file $FILE" | tee -a $LOG_FILE
    BASE=$(basename $FILE | sed 's@\(\.fq\|\.fastq\|\.fastq.gz\)@@')
    WORKING_DIR=$2
    # cutadapt -f fastq -e 0.1 -O 6 -q 20 -m 35 -a AGATCGGAAGAGC -o $WORKING_DIR"/"$BASE".trim.fastq" $FILE 
    cutadapt -e 0.1 -O 6 -q 20 -m 35 -a AGATCGGAAGAGC $FILE -o $WORKING_DIR/$BASE.trim.fastq
}


export -f trimfun
parallel -j $TOTAL_PROC_NO --no-notice trimfun {} $WORKING_DIR  ::: $FASTQ_LIST 1>>$LOG_FILE 2>>$LOG_ERR_FILE
unset -f trimfun
wait; echo -e "(`date`) 2 Trimming finished." | tee -a $LOG_FILE

echo -e "###########################################" | tee -a $LOG_FILE
echo -e "(`date`) 2 QC of trimmed files" | tee -a $LOG_FILE
cd $WORKING_DIR
WORKING_DIR=$WORKING_DIR'/fastqc'
mkdir -p $WORKING_DIR
# ls -1 *.f*q 2>&1 | tee $LOG_FILE # debug
# set -x # debug
ls -1 *.f*q | xargs -n 1 -P $SAMPLE_NO fastqc -t $NPROC_PER_SAMPLE -outdir $WORKING_DIR 1>>$LOG_FILE 2>>$LOG_ERR_FILE
# set +x # debug
wait;echo -e "(`date`) 2 QC trimmed fastqc finished." | tee -a $LOG_FILE

echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 3.1: alignment" | tee -a $LOG_FILE

WORKING_DIR=$PARENT_DIR'/03alignment/'
mkdir -p $WORKING_DIR
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

date | tee -a $LOG_FILE_STEP

alignfun (){
    while read data;do
    samplename=${data%.*.*}
    echo -e "(`date`) align $samplename" | tee -a $LOG_FILE
    # bowtie2 -t --non-deterministic --mm --phred33 -p $NPROC_PER_SAMPLE --very-sensitive -x mm10 -U $data | samtools sort -@ $NPROC_PER_SAMPLE -T $samplename -o $WORKING_DIR"/"$samplename".bam" &
    bowtie2 -t --non-deterministic --mm --phred33 -p $NPROC_PER_SAMPLE --very-sensitive -x hg38 -U $data | samtools sort -@ $NPROC_PER_SAMPLE  -o $WORKING_DIR"/"$samplename".align.bam" &
done
}
echo -e "(`date`) Performing alignment on files in here $PWD" | tee -a $LOG_FILE
ls -1 *.fastq | alignfun | tee -a $LOG_FILE_STEP 2>> $LOG_FILE_STEP


echo -e "This version of bowtie2 was used for the analysis" | tee -a $LOG_FILE_STEP

bowtie2 --version | tee -a $LOG_FILE_STEP

echo -e "This version of samtools was used for the analysis" | tee -a $LOG_FILE_STEP

samtools --version | tee -a $LOG_FILE_STEP
wait;echo -e "`date`: Step 3.1 alignment finished!" | tee -a $LOG_FILE


#------------------------------------------------------------
cd $WORKING_DIR
echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 3.2: index bam file" | tee -a $LOG_FILE

date | tee -a $LOG_FILE_STEP
ls -1 *.bam | xargs -n 1 -P $SAMPLE_NO -I {} \
                    samtools index {} \
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

#alignstatfun2 (){
#    while read data;do
#    samplename=${data%.*}_align.bam
#    echo -e "(`date`) save stats for$samplename" | tee -a $LOG_FILE
#
#    echo "QC-passed reads:" >> $STATS_DIR/${samplename}.stats
#    grep total -m 1 $WORKING_DIR"/"$samplename"_align.stat" | cut -f1  -d " " >> $STATS_DIR/${samplename}.stats
#    echo "Mapped reads:" >> $STATS_DIR/${samplename}.stats
#    grep mapped -m 1 $WORKING_DIR"/"$samplename"_align.stat" | cut -f1  -d " " >> $STATS_DIR/${samplename}.stats
#done
#}
ls -1 *.bam | alignstatfun | tee -a $LOG_FILE_STEP 2>> $LOG_FILE_STEP

wait;echo -e "(`date`) Step 3.3 calculate stats finished!" | tee -a $LOG_FILE
#------------------------------------------------------------
echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 3.4: Quality mapping" | tee -a $LOG_FILE
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

ls *.bam | qualmapfun 1>>$LOG_ERR_FILE 2>>$LOG_FILE_STEP
echo -e "This version of qualimap was used for the analysis" | tee -a $LOG_FILE_STEP
qualimap --version >> $LOG_FILE_STEP 2>> $LOG_FILE_STEP

wait;echo -e "(`date`) Step 3.4 quality mapping finished!" | tee -a $LOG_FILE

echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 4.1: filtering the aligned bam files" | tee -a $LOG_FILE
STEP='04filter'
WORKING_DIR=$PARENT_DIR'/'$STEP
mkdir -p $WORKING_DIR
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

filter (){
    while read data;do
    samplename=${data%.*.*}
    echo -e "(`date`) Filter $samplename for quality >30" | tee -a $LOG_FILE

    samtools view -b -q 30 -@ $NPROC_PER_SAMPLE ${data} > $WORKING_DIR"/"$samplename".filter.bam"
done
}


date | tee -a $LOG_FILE_STEP
ls -1 *.bam | filter | tee -a $LOG_FILE_STEP 2>> $LOG_FILE_STEP
date | tee -a $LOG_FILE_STEP

wait; echo -e "(`date`)  Step 4.1 filter bam finished" | tee -a $LOG_FILE $LOG_FILE_STEP


echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 4.2: indexing the filtered bam" | tee -a $LOG_FILE

cd $WORKING_DIR
date | tee -a $LOG_FILE_STEP
ls *.bam | parallel --no-notice --progress -j $SAMPLE_NO samtools index {}
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



echo -e "############################################################" | tee -a $LOG_FILE

 echo -e "(`date`) Starting Step 4.4: Calculating NRF, PCR bottleneck coefficient" | tee -a $LOG_FILE


 date | tee -a $LOG_FILE_STEP

 # PBC File output
 # TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair

 NRFqc () {
   while read data
   do
 	samplename=${data%.*}
 	echo -e "(`date`) calculating NRF, PCR bottleneck coefficients for $samplename" | tee -a $LOG_FILE_STEP
 	bedtools bamtobed -i $data | awk `BEGIN{OFS="\t"}{print $1,$2,$3,$6}` | grep -v `chrM` | sort | uniq -c | awk `BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}` > ${samplename}.txt 2>>$LOG_ERR_FILE
   done
 }


 ls *.bam | NRFqc | tee -a 1>>$LOG_FILE_STEP 2>>$LOG_ERR_FILE

 date | tee -a $LOG_FILE_STEP

 echo -e "(`date`)  finished Step 4.4: Calculating NRF, PCR bottleneck coefficient" | tee -a $LOG_FILE


echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 5: deduplicate" | tee -a $LOG_FILE

STEP='05deDupped'
WORKING_DIR=$PARENT_DIR'/'$STEP
mkdir -p $WORKING_DIR;
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

date | tee -a $LOG_FILE_STEP
mkdir -p $WORKING_DIR"/tmp"
#JAVA_OPTS=-Xmx64g
dedupfun (){
    while read data; do
    namestr=$(echo "$data"|cut -f1 -d".")
    picard MarkDuplicates I=$data O=$WORKING_DIR"/"$namestr".bam" REMOVE_DUPLICATES=true ASSUME_SORTED=true M=$WORKING_DIR"/"$namestr".txt" TMP_DIR=$WORKING_DIR"/tmp" QUIET=F 
done
}
ls -1 *.bam | dedupfun 1>>$LOG_FILE_STEP 2>>$LOG_FILE_STEP
rm -r $WORKING_DIR"/tmp"
date | tee -a $LOG_FILE_STEP
echo -e "This version of MarkDuplicates from picard was used for the analysis" | tee -a $LOG_FILE_STEP
picard MarkDuplicates --version >> $LOG_FILE_STEP 2>> $LOG_FILE_STEP


cd $WORKING_DIR
date | tee -a $LOG_FILE_STEP
ls *.bam | parallel --no-notice --progress -j $SAMPLE_NO samtools index {}
1>>$LOG_ERR_FILE 2>>$LOG_FILE_STEP

#prefix=".stat"
#ls *.bam | while read data; do samtools flagstat "$data" > "$data"${prefix}  & done | tee -a $LOG_FILE_STEP

echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 5.2: Calculate stats" | tee -a $LOG_FILE

dstatfun (){
    while read data;do
    samplename=${data%.*}
    echo -e "(`date`) calc stats for $samplename" | tee -a $LOG_FILE
    samtools flagstat -@ $NPROC_PER_SAMPLE ${data}  > $WORKING_DIR"/"$samplename".stat"  &
done
}

ls -1 *.bam | dstatfun | tee -a $LOG_FILE_STEP 2>> $LOG_FILE_STEP

wait;echo -e "(`date`) Step 5.2 calculate stats finished!" | tee -a $LOG_FILE


echo -e "(`date`)  Step 5 deduplication finished" | tee -a $LOG_FILE


echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 6: generate genome browser tracks " | tee -a $LOG_FILE
STEP='06tracks'
WORKING_DIR=$PARENT_DIR'/'$STEP
mkdir -p $WORKING_DIR;
LOG_FILE_STEP=$WORKING_DIR"/log.txt"
LOG_ERR_FILE=$WORKING_DIR"/err.txt"

bcfun (){
	BW=$(basename $1 | sed "s|\.bam|.bw|")
	OUT_BW=$2/$BW	
	echo -e "(`date`) generate tracks for $1"
	bamCoverage --binSize 10 --smoothLength 30 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --extendReads 200 -p $3 -b $1 -o $OUT_BW
	echo -e "(`date`) created $BW"

}

export -f bcfun
BAM_LIST=$(ls -1 *.bam)
parallel -j $TOTAL_PROC_NO --no-notice bcfun {} $WORKING_DIR $NPROC_PER_SAMPLE ::: $BAM_LIST 1>>$LOG_FILE_STEP 2>>$LOG_ERR_FILE
unset -f bcfun

wait; echo -e "(`date`) Step 6: generate genome browser tracks finished" | tee -a $LOG_FILE
echo -e "############################################################" | tee -a $LOG_FILE
echo -e "This script has finished running" | tee -a $LOG_FILE


echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`)Starting Step 7: spp cross-correlation " | tee -a $LOG_FILE
STEP='07spp'
WORKING_DIR=$PARENT_DIR'/'$STEP
mkdir -p $WORKING_DIR;
LOG_FILE_STEP=$WORKING_DIR"/log.txt"
LOG_ERR_FILE=$WORKING_DIR"/err.txt"

sppfun (){
    while read data; do
	nameStr=$(echo "$data"| cut -f1 -d".")
	Rscript /home/shared/phantompeakqualtools/run_spp.R -c="$data" -savp=$WORKING_DIR"/"$nameStr".pdf" -p=$NPROC_PER_SAMPLE &
    done
}

ls -1 *.bam | sppfun 1>> $LOG_FILE_STEP 2>> $LOG_ERR_FILE



STEP='qulimap'
WORKING_DIR=$WORKING_DIR'/'$STEP
mkdir -p $WORKING_DIR;

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

ls *.bam | qualmapfun 1>>$LOG_ERR_FILE 2>>$LOG_FILE_STEP
echo -e "This version of qualimap was used for the analysis" | tee -a $LOG_FILE_STEP
qualimap --version >> $LOG_FILE_STEP 2>> $LOG_FILE_STEP


##plotFingerprint from deeptools is run
echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Plotting fingerprint" | tee -a $LOG_FILE
plotFingerprint --bamfiles *.bam --numberOfProcessors 2 -o $STATS_DIR/fingerprint.png --outQualityMetrics $STATS_DIR/fingerprintQC.txt | tee -a $LOG_FILE_STEP

echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Saving stats" | tee -a $LOG_FILE

echo -e "Sample\t library_size_dedup" > $STATS_DIR/library_size.tsv   

savefun (){
    while read data;do
    samplename=${data%.*}


    QC=$(grep total -m 1 $PARENT_DIR"/03alignment/"$samplename".align.stat" | cut -f1  -d " ")
    echo -e "QC-passed reads:\t $QC" > $STATS_DIR/${samplename}_stats.tsv

    MAP=$(grep mapped -m 1 $PARENT_DIR"/03alignment/"$samplename".align.stat" | cut -f1  -d " ")
    echo -e "Mapped reads:\t $MAP" >> $STATS_DIR/${samplename}_stats.tsv

    HQ=$(grep mapped -m 1 $PARENT_DIR"/04filter/"$samplename".filter.stat" | cut -f1  -d " ")
    echo -e "High quality (q>=30) reads:\t $HQ" >> $STATS_DIR/${samplename}_stats.tsv

    FIN=$(grep mapped -m 1 $PARENT_DIR"/05deDupped/"$samplename".stat" | cut -f1  -d " ")
    echo -e "De-duplicated reads:\t $FIN" >> $STATS_DIR/${samplename}_stats.tsv

    echo -e "${samplename}\t $FIN" >> $STATS_DIR/library_size.tsv   

done
}
ls -1 *.bam | savefun 1>> $LOG_FILE_STEP 2>> $LOG_ERR_FILE | tee -a $LOG_FILE

## multiqc
cd $PARENT_DIR
multiqc ./



wait; echo -e "(`date`) Step 7: QC generate spp plots finished" | tee -a $LOG_FILE


echo -e "############################################################"
echo -e "This script has finished running"
