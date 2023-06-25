#* this pipe is integrated pipeline for Macrohaplotype de-UMI and mapping
#* author: Xuewen Wang
#* version: 0.1, April 2023
#* tested working in Ubuntu 20.04

###Configure paths of tools, may be exported to $PATH
scriptdir="/mnt/data0/xuewen/macrohaplotype/scripts"
toolMinimap2=/home/xuewen/tools/minimap2-2.24_x64-linux
toolSamtools=/usr/bin/samtools #1.10, 1.11 tested
toolUmidir=/usr/local/bin/umi_tools #UMI-tools version: 1.1.2
toolcutadpt=/usr/bin/cutadapt # cutadapt 2.8 with Python 3.8.10

### input args
refGenome=/mnt/data0/xuewen/hg38giab/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta 
refn="GRCh38" #short Name for genome
b="hg002rep123.bam"
readfastq=hg002rep2.fastq.gz #output read file: suffix .fastq.gz
outdir=/mnt/data0/xuewen/macrohaplotype/hg002_trio_PBhifi
threadn=40
#prefix="18bp" #prefix of output files
prefix="trio" #prefix of output files
umiQualityThreshold=25
umiProcess="No" #Yes for umi process, No for not
#results in an outbam file with automated file name 

###prepare
date
if [ ! -d $outdir ];then    
	mkdir $outdir
fi
cd $outdir 
refdir=$(dirname $refGenome)
refseq=$(basename $refGenome)
qrydir=$(dirname $readfastq)
qseq=$(basename $readfastq)
fq=$(basename -s .fastq.gz $qseq)

if [ ! -f ${b}.bai ]; then
    samtools index -b $b -@ $threadn
fi

#### hifi bam reads to fastq reads
 #get fastq
#$ samtools fastq $b -@ $threadn > ${b/bam/fastq}
#$ gzip ${b/bam/fastq} ## output: ${b/bam/fastq.gz}
#$ samtools stats $b >${b}.stat.txt 
 python3 $scriptdir/Fastq_stat.py -i ${b/bam/fastq.gz} >${b/bam/fastq}.stat.xls
#$ echo fastq file: ${b/bam/fastq.gz}


#### read quality control : optional
## quality filter for reads: Q40 preferred
minlen=6000
maxlen=50000
infq=$readfastq
QualityMean=30 #same as -e 0.0001
bn=`basename $infq`  #m64254e_220923_234542.hifi_reads.fastq.gz
outfq=${bn/.fastq.gz}.${QualityMean}.fastq.gz
echo $outfq
#/usr/local/bin/fastq-filter from https://github.com/LUMC/fastq-filter
 fastq-filter  -l $minlen -L $maxlen -o $outfq -q $QualityMean $infq
readfastq=$outfq

####remove adaptors: 24bp, linked^ adptors left+right, no --revcomp
#cutadapt -a "CAAGCAGAAGACGGCATACGAGAT...GATCTCGGTGGTCGCCGTATCATT"  --overlap 20  --action trim --error-rate 0.1 --quality-base 33 --cores $threadn -o $fq.noAdpt.fastq.gz $qseq


#### umi precessing
	#fq="m64254e_221106_121349.hifi_reads.Q30"
if [ $umiProcess == "Yes" ] ; then
	umi_tools extract --extract-method=string --stdin=$readfastq --quality-encoding=phred33 --quality-filter-threshold=$umiQualityThreshold --bc-pattern=NNNNNNNNNNNNNNNNNN --log=processed.log --stdout ${prefix}.${fq}.umiprocessed.fastq.gz 

	
	#use regex for input fastq which is after cutting adaptors
	#umi_tools extract --extract-method=regex --stdin=$fq.noAdpt.fastq.gz --quality-encoding=phred33 --quality-filter-threshold=$umiQualityThreshold --bc-pattern='(?P<umi_1>.{16,18}).+(?P<umi_2>.{16,18})'--log=processed.log --stdout ${prefix}.${fq}.umiprocessed.regex.fastq.gz 
	# zcat m64254e_221106_121349.hifi_reads.Q30.processed.fastq.gz|grep "@m" > 42bp.precessed.read.id #bad
	echo "umi precessed: " ${prefix}.${fq}.umiprocessed.fastq.gz
	readfastq=${prefix}.${fq}.umiprocessed.fastq.gz
	pwd
	
fi	

###mapping to reference genome
    #qseq=${prefix}.${fq}.umiprocessed
	qseq=$(basename $readfastq)
	$toolMinimap2/minimap2  -x map-hifi -a -t $threadn --cs="long" $refdir/$refseq $outdir/$qseq >${qseq}_${refn}.sam
   ## samtools v1.11
	samtools view -b -@ $threadn ${qseq}_${refn}.sam >${qseq}_${refn}.bam
	samtools sort  --output-fmt BAM -@ $threadn ${qseq}_${refn}.bam >${qseq}_${refn}srt.bam
	mv ${qseq}_${refn}srt.bam ${qseq}_${refn}.bam
	samtools index -b -@ $threadn ${qseq}_${refn}.bam
	rm ${qseq}_${refn}.sam
	echo "Mapping done " ${qseq}_${refn}.bam
	pwd
	
###dedup: keep one best read per umi
if [ $umiProcess == "Yes" ] ; then
	inbam=${qseq}_${refn}.bam
	outdedupBam=${qseq}_${refn}.umi.dedup.bam 
	umi_tools dedup --stdin=$inbam --log=dedup.log  --stdout=$outdedupBam 
	samtools index -b $outdedupBam
	echo "umi-deduped bam:" $outdedupBam
	pwd
fi
	
###clean
if [ $umiProcess == "Yes" ] ; then
	#rm ${qseq}_${refn}.bam*
fi
	date
	

