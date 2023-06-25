# MacroHapCaller

* __phased multiple alleles calling for long sequencing reads, e.g. PacBio HiFi__

*  __call variants of STR/SSR, SNP, InDel at the same time__

* __run faster with parallel computing__

Macrohaplotype caller for STR, SNP, and InDel from NGS long-read sequencing data, especially for PacBio HiFi reads
The MacroHapCaller calls targeted STRs/SSRs, SNPs, and InDels sites simultaneously in a row from each single NGS read and clusters the variants into phased haplotype strings. MacroHapCaller is the best tool to analyze all haplotypic variants from genetically inherited DNA. It suits diploid, polyploid, and complex DNA mixtures from many individuals, e.g. DNA forensics. MacroHapCaller is programmed in Java with parallel computing enabled so it can run on any computing platform. 

The MacroHapCaller is also integrated into a pipeline with a shell script or Python3. The pipeline takes the BAM or fastq reads as input, align reads to the genome, sort, and index alignment, and then user can use the output to run MacroHapCaller. Additional tools are also available as util tools. For util usage : `java -jar appName.jar` where appName should be replaced with a detailed tool name. 


MacroHapCaller runs fast. It takes ~ 2 mins for merged data from four full runs of PacBio SMRTcell HiFi reads on a Linux machine with 12 threads, and a maximum allowed 120G RAM. 

## Latest version
V0.3

## Usage
For help: 

`java -jar MacroHapCaller##.jar`

where ## is the version number, e.g.: 0.3.

## Command

usage: `java -jar -Xmx100G MacroHapCaller0.3.jar [options]`

 **-a,--strAnchorFile** <arg>  &nbsp; configure file of anchors for STRs in
                                 tabular plain text
 
**-c,--MinReadCount** &nbsp; integer, minimum read count to report a
                                 macrohaplotype, default [2]
 
 **-d,--inDelBedFile** <arg> &nbsp; configure file of InDels in BED format in
                                 tabular plain text
 
 **-i,--input** <arg> &nbsp; input BAM file with a path
 
 **-l,--homopolymerLen** &nbsp; integer, homopolymer length threshold to
                                 skip call variant,default [4]
 
 **-m,--MaxMismatches**  &nbsp; integer, maximum distance of edit
                                 distance in anchor match, default [1]
 
 **-n,--snpPanelPosFile** <arg> &nbsp; configure file of SNPs in tabular plain
                                 text
 
 **-o,--output** <arg> &nbsp; output file name for macrohaplotype
                                 result
 
 **-p,--MinReadProportion** &nbsp; double, minimum read proportion of all
                                 reads for each locus to report a macrohaplotype, default
                                 [0.001]
 
 **-q,--LowQualThreshold** &nbsp; integer, phred-scale quality threshold
                                 for variant bases, default [15]
 
 **-r,--refGenomeFastaFile** <arg> &nbsp; Genome reference sequence file in .fasta
                                 format
 
 **-t,--ThreadNumber** &nbsp; integer, the number of computing threads,
                                 default [12]

 
 ## Example 
 
    java -jar -Xmx120G MacroHapCaller0.3.jar -i /mnt/data0/xuewen/macrohaplotype/hg002/Q40hifi/hg002.8kampl.Q40.fastq.gz_GRCh38.bam -o hg002.8kampl.Q40_GRCh38.tsv -a /mnt/data0/xuewen/macrohaplotype/scripts/MHvar_v0.3/CODISSTR_anchor.XW.config_v0.2.txt -d /mnt/data0/xuewen/macrohaplotype/scripts/MHvar_v0.3/MHindels_v0.2.bed -n /mnt/data0/xuewen/macrohaplotype/scripts/MHvar_v0.3/MHsnps.pos.txt -r /mnt/data0/xuewen/hg38giab/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta -l 2 -m 1 -q 15 -c 100 -p 0.01 -t 12

## Preprocess pipeline for MacroHapCaller
 This pipeline will preprocess the originall PacBio HiFi reads in BAM to fastq, read quality control, umi-analysis (optional), statistical summary, read-reference alignment, sort, and index. The major results are a read-reference alignment in the BAM format and BAM index, as well as statistical report files.
 
 `MH_umi_dedup_map_pipe.sh`
 
 

## Input of MacroHapCaller
 1. a read-reference alignment file in [BAM format](https://en.wikipedia.org/wiki/SAM_(file_format)). This alignment file could be generated from  `MH_umi_dedup_map_pipe.sh` or user-generated alignment file in .Bam format

 2. Configure files of targeted variant positional information for each of STRs, SNPs, and InDels

## Configure file format for STR
The configure file listed the STR loci information in a tabular separated text file. An example file for forensic CODIS STR of human genome is given in "CODISSTR_anchor.XW.config_v0.2.txt"
Each line for one STR locus. There is one fixed headline before starting the detailed STR locus. Every line has 11 columns. The integer value for the last column is the estimated maximum length of TR, where 500 is enough for most cases. The coordinate starts from 1.

    Chrom	ChromStartPos_Str	ChromEndPos_Str	Name	Repeat_unit_length	Motif_must_present	Inner_offset	Anchor_left	Anchor_right	LeftAnchorStartPos	RightAnchorEndPos	MaxRefLength
    chr1	230765214	230765259	D1S1656@STR_118267	2	CACACA	0	CAGAAAATGAGAACACATG	GGTTATGCCAAAAGGGC	230765190	230765281	500
    chr1	230769616	230769683	D1S1656	4	TCTATCTATCTA	0	TCAGAGAAATAGAATCACTA	TGAGCAACACAGGCTTGAAC	230769556	230769721	500
    chr2	1489653	1489684	TPOX	4	AATGAATG	0	CAGAACAGGCACTTAGGGAA	AACGCTGACAAGGACAGAAG	1489624	1489716	500

    

## Configure file format for InDel
The config file for InDel is in a BED format tabular separated plain text file, which coordinates start from 0. The REF and ALT are the bases for InDel in the reference genome and variant. An example file is given for 8-kb fragments containing human CODIS STR. For Details, please read our paper. There is one fixed headline and 5 columns in each line. Then each subsequent line is for one InDel site. If multiple InDel for the same coordinate, list one of them is enough. 

    Chrom	ChromStart	ChromEnd	REF	ALT
    chr1	230764747	230764748	GA	G
    chr1	230764749	230764750	T	TC
    chr2	1485704	1485705	C	CT
    chr2	1486795	1486803	TAGTGGGG	T
    chr2	1487522	1487526	CAGG	C
    chr2	1489924	1489975	GCACACAGGAGGAGTCACGACAGAGCAGTGTAAGAGCCGCCACGTGGGTCC	G

## Configure file format for SNP
The config file for SNP is in a tabular separated plain text file, which coordinates start from 1.  There is one fixed headline and 5 columns in each line. Then each subsequent line is for one SNP site. An example file is given for 8-kb fragments containing 20 human CODIS STRs. the "ID" column is the name of a DNA fragment, which uses the same name as one of the STR in this 8-kb fragment as STR configure. All SNPs with the same "ID" will be put together with STR alleles in the same fragment once the ID of STR is shared between SNP and STR configure file. 

    #CHROM	 POS 	ID	REF	ALT
    chr1	230765010	D1S1656	A	G
    chr1	230765280	D1S1656	G	A
    chr1	230765336	D1S1656	C	T
    chr1	230765352	D1S1656	G	C
 
 ## Output of MacroHapCaller

Macrohaplotypes and supported read count in tabular text file .tsv. 
E.g., the two macrohaplotypes around FBI's CODIS loci D2S441 in 8 kb PacBio HiFi reads of benchmark reference hg002. Locus details are at [FBI](https://www.fbi.gov/how-we-can-help-you/dna-fingerprint-act-of-2005-expungement-policy/codis-and-ndis-fact-sheet)
     
The first line shows the total reads for this locus D2441. The position lines present the coordinate of each targeted variant site on chr2 in the human genome Hg38. Then the macrohaplotypes, including three parts of SNPs, InDels and STRs separated by ";".
     
    #Total hapVar: 	D2S441	27694	
    #Markername	Counts	HapVarLen	hapVar(s)
    #			Position:68011922,68011990,68012066,68012109,68012198,68012335,68012415,68012463,68012596,68012616,68012754,68012776,68012877,68012888,68013086,68013124,68013388,68013405,68013406,68013664,68013692,68013770,68013988,68014065,68014213,68014245,68014581,68014700,68015087,68015100,68015130,68015131,68015330,68015388,68015498,68015787,68015976,68016033,68016045,68016108,68016576,68016583,68017126,68017147,68017242,68017364,68017504,68017573,68017732,68017739,68018169,68018176,68018322,68018330,68012435,68012478,68012608,68012694,68012722,68012762,68012781,68011917
    
    D2S441	8646	136	
    G,A,A,A,G,C,T,A,T,C,C,T,G,T,T,G,C,C,G,C,A,A,G,T,A,T,A,G,A,A,C,G,T,G,C,G,G,G,T,C,C,G,C,G,G,A,A,C,C,C,T,A,A,G;TA,CACATATATA,CA,CA,CACATATATA,CATATATA,ATGT;TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTA

    D2S441	7467	152	
    G,A,A,A,G,C,T,A,T,T,C,T,A,T,C,G,C,C,T,C,A,G,G,T,G,T,A,G,A,A,C,G,C,G,C,G,G,G,T,C,C,G,C,G,A,A,A,C,C,C,T,G,A,G;TAC,C,CA,CA,CACATATATA,CATATAT,T;TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATTTATCTATCTA
    
   ## Visualization 
   The macrohaplotypes can be visualized or edited in a text processor, Microsoft excels, and spreadsheett.
   
   The graphically view and compare the macrohaplotypes, please use our graphic tool: [USAT](https://github.com/XuewenWangUGA/USAT)
   
  ## Fund
  
  This work was supported in part by award 15PNIJ-21-GG-04159-RESS, awarded by the National Institute of Justice, Office of Justice Programs, U.S. Department of Justice.
  
   ## Citation
   coming soon





