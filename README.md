# WXS_VariantCalling
Variant Calling by GATK, Samtools, Varscan

Tools 
BWA=/export/tools/bwa
SAMTOOLS=/export/tools/samtools/bin/samtools
GATK=/export/tools/GATK/GenomeAnalysisTK.jar
PICARD=/export/tools/picard.jar
Reference
genomeFasta=$WORKDIR/GATK_resource/ucsc.hg19.fasta

Command Line Parameters
Step 1: Converting BAMs to FASTQs with Biobambam - biobambam2 2.0.54
bamtofastq \
collate=1 \
exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \
filename= <input.bam> \
gz=1 \
inputformat=bam
level=5 \
outputdir= <output_path> \

#Quality control
fastqc #FastQC for FASTQ/BAM file quality assessment (1 enable,0 disable).

#Initial alignment
#BWA-MEM, for read length 70bp-1Mbp, default algorithm for initial alignment
matchScore=1 #score for a sequence match. -A for bwa mem.
mismatchPenalty=4 #penalty for a mismatch. -B for bwa mem.
gapOpenPenalty=6 #penalty for gap open. -O for bwa mem.
gapExtensionPenalty=1 #penalty for gap extension; a gap of size k cost: gapOpenPenalty + gapExtensionPenalty*k. -E for bwa mem
clipPenalty=5 #penalty for clipping. -L for bwa mem
readSingletonPenalty=17 #penalty for unpaired read pair. -U for bwa mem.

Step 2: BWA Alignment - bwa 0.7.15 - samtools 1.3.1
If mean read length is greater than or equal to 70bp:
bwa mem \
-t 8 \
-T 0 \
-R <read_group> \
<reference> \
<fastq_1.fq.gz> \
<fastq_2.fq.gz> |
samtools view \
-Shb
-o <output.bam> -

If mean read length is less than 70bp:
bwa aln -t 8 <reference> <fastq_1.fq.gz> > <sai_1.sai> &&
bwa aln -t 8 <reference> <fastq_2.fq.gz> > <sai_2.sai> &&
bwa sampe -r <read_group> <reference> <sai_1.sai> <sai_2.sai> <fastq_1.fq.gz> <fastq_2.fq.gz> | samtools
view -Shb -o <output.bam> -
If the quality scores are encoded as Illumina 1.3 or 1.5, use BWA aln with the “-l” flag.

Step 3: BAM Sort - picard 2.6.0
java -jar picard.jar SortSam \
CREATE_INDEX=true \
INPUT=<input.bam> \
OUTPUT=<output.bam> \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=STRICT

Step 4: BAM Merge - picard 2.6.0
java -jar picard.jar MergeSamFiles \
ASSUME_SORTED=false \
CREATE_INDEX=true \
[INPUT= <input.bam>] \
MERGE_SEQUENCE_DICTIONARIES=false \
OUTPUT= <output_path> \
SORT_ORDER=coordinate \
USE_THREADING=true \
VALIDATION_STRINGENCY=STRICT

Step 5: Mark Duplicates - picard 2.6.0
java -jar picard.jar MarkDuplicates \
CREATE_INDEX=true \
INPUT=<input.bam> \
VALIDATION_STRINGENCY=STRICT

Step 6:  RealignTargetCreator
java -jar GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R <reference>
-known <known_indels.vcf>
[ -I <input.bam> ]
-o <realign_target.intervals>

Step 7: IndelRealigner
java -jar GenomeAnalysisTK.jar \
-T IndelRealigner \
-R <reference> \
-known <known_indels.vcf> \
-targetIntervals <realign_target.intervals> \
--noOriginalAlignmentTags \
[ -I <input.bam> ] \
-nWayOut <output.map>

Step 8: BaseRecalibrator
java -jar GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R <reference> \
-I <input.bam> \
-knownSites <dbsnp.vcf>
-o <bqsr.grp>

Step 9: PrintReads
java -jar GenomeAnalysisTK.jar \
-T PrintReads \
-R <reference> \
-I <input.bam> \
--BQSR <bqsr.grp> \
-o <output.bam>

MuTect2
GATK
java -jar GenomeAnalysisTK.jar \
-T MuTect2 \
-R <reference> \
-L <region> \
-I:tumor <tumor.bam> \
-I:normal <normal.bam> \
--normal_panel <pon.vcf> \
--cosmic <cosmic.vcf> \
--dbsnp <dbsnp.vcf> \
--contamination_fraction_to_filter 0.02 \
-o <mutect_variants.vcf> \
--output_mode EMIT_VARIANTS_ONLY \
--disable_auto_index_creation_and_locking_when_reading_rods

OR

SomaticSniper
Somatic-sniper v1.0.5.0
1 bam-somaticsniper \
2 -q 0 \
3 -Q 15 \
4 -s 0.01 \
5 -T 0.85 \
6 -N 2 \
7 -r 0.001 \
8 -n NORMAL \
9 -t TUMOR \
10 -F vcf \
11 -f ref.fa \
12 <tumor.bam> \
13 <normal.bam> \
14 <somaticsniper_variants.vcf>

VarScan
Mpileup; Samtools 1.1
samtools mpileup \
-f <reference> \
-q 1 \
-B \
<normal.bam> \
<tumor.bam> >
<intermediate_mpileup.pileup>

Varscan Somatic; Varscan.v2.3.9
java -jar VarScan.jar somatic \
<intermediate_mpileup.pileup> \
<output_path> \
--mpileup 1 \
--min-coverage 8 \
--min-coverage-normal 8 \
--min-coverage-tumor 6 \
--min-var-freq 0.10 \
--min-freq-for-hom 0.75 \
--normal-purity 1.0 \
--tumor-purity 1.00 \
--p-value 0.99 \
--somatic-p-value 0.05 \
--strand-filter 0 \
--output-vcf

Varscan ProcessSomatic; Varscan.v2.3.9
java -jar VarScan.jar processSomatic \
<intermediate_varscan_somatic.vcf> \
--min-tumor-freq 0.10 \
--max-normal-freq 0.05 \
--p-value 0.07

Variant Effect Predictor (VEP) v84
Variants in the VCF files are also matched to known variants from external mutation databases. The following databases are used for VCF annotation:
• GENCODE v.22
• sift v.5.2.2
• ESP v.20141103
• polyphen v.2.2.2
• dbSNP v.146
• Ensembl genebuild v.2014-07
• Ensembl regbuild v.13.0
• HGMD public v.20154
• ClinVar v.201601
