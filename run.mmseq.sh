
#Download software
cd '/bubo/home/h26/edsgard/opt'
wget 'http://bgx.org.uk/software/mmseq_0.9.18.zip'
wget 'http://bgx.org.uk/software/mmseq_0.10.0.zip'
wget 'http://bgx.org.uk/software/polyHap2_20110607.zip'

#Download annotation
#Transcriptome, Genome, GFF gene annot, Bowtie indexes. All Ensembl v64
cd '/bubo/home/h26/edsgard/glob/annotation/mmseq'
wget 'http://bgx.org.uk/software/Homo_sapiens.GRCh37.64.ref_transcripts.fa.gz' &
wget 'http://bgx.org.uk/software/Homo_sapiens.GRCh37.64.dna.toplevel.ref.fa.gz'
wget 'http://bgx.org.uk/software/Homo_sapiens.GRCh37.64.ref.gff.gz' #(TBD: If you sequenced polyadenylated transcripts only, include only the cDNA (not the ncRNA) transcript sequences in the FASTA file.)
wget 'http://bgx.org.uk/software/Homo_sapiens.GRCh37.64.dna.toplevel.ref.ebwt.tgz'
gunzip Homo_sapiens.GRCh37.64.ref_transcripts.fa.gz
gunzip Homo_sapiens.GRCh37.64.dna.toplevel.ref.fa.gz
gunzip Homo_sapiens.GRCh37.64.ref.gff.gz
tar -xvzf Homo_sapiens.GRCh37.64.dna.toplevel.ref.ebwt.tgz

#Add binaries to PATH
module load bioinfo-tools
module load bowtie/0.12.7
module load tophat/1.3.3
#module load samtools/0.1.18 Conflicts with above module load (either bowtie or tophat)
module load java/sun_jdk1.6.0_18
export CLASSPATH=$CLASSPATH:/bubo/home/h26/edsgard/opt/polyHap2_20110607/polyHap2.jar
export PATH=$PATH:/bubo/home/h26/edsgard/opt/mmseq_0.9.18


#Set variables
export BOWTIE_INDEXES='/bubo/home/h26/edsgard/glob/annotation/mmseq/' #trailing slash needed
GENOME_FASTA='Homo_sapiens.GRCh37.64.dna.toplevel.ref.fa'
GENOME_EBWT=`basename $GENOME_FASTA .fa`
GFF_FILE='Homo_sapiens.GRCh37.64.ref.gff'
annot=${BOWTIE_INDEXES}${GFF_FILE}
ref=$GENOME_EBWT
email='daniel.edsgard@scilifelab.se'
projid='b2010035'
execdir='/bubo/home/h26/edsgard/glob/code/ngs_pipes/hs_pipe'
sbatchdir='/proj/b2011075/analysis/mmseq/outscripts/run1'
outdir='/proj/b2011075/analysis/mmseq/outdata/run1'
PH2TEMP=${outdir}/phasing


###
#Tophat Post-QC
###
execdir='/bubo/home/h26/edsgard/glob/code/ase'
fastqdir='/proj/b2012046/edsgard/ase/sim/data/synt/fastqfilt'
fastqsuffix='.filter.fastq'
outdir='/proj/b2012046/edsgard/ase/sim/data/mmseq/tophat'
sbatchdir='/proj/b2012046/edsgard/ase/sim/scripts/mmseq/tophat'

#Vars
threads='8'
annotdir='/bubo/home/h26/edsgard/glob/annotation/mmseq/'
gtf=${annotdir}'Homo_sapiens.GRCh37.64.ref.gff'
ref='Homo_sapiens.GRCh37.64.dna.toplevel.ref'
time='35:00:00' #Longest run: 7h10min: Input: ~7.3G per readpair fastq file
readlen='100'
isizefile='/proj/b2011075/data/rnaseq/bioanalyzer/isizes.only.tab'
isizedev='75' #very rough approx from the bioanalyzer plots. Also, we use trimming which further adds to the deviation.

#TBD: Check if -g1 flag was used when ran the real data
cd $sbatchdir
rm cmds.sh
samples=(`cat samples.list`)
for sample in ${samples[@]}
do
    find $fastqdir -maxdepth 1 -name ${sample}'*' | sed 's/[12].filter.fastq//' | grep -v '.S.filter.fastq' | sort -u >fastqfiles.prefix.list
    cat fastqfiles.prefix.list | xargs -I% echo perl ${execdir}/'tophat.gensbatch.pl' $sample % $fastqsuffix $outdir $threads $gtf $ref $isizefile $sbatchdir $projid $email $time $readlen $isizedev $annotdir >>cmds.sh
done
sh cmds.sh
find . -name '*.tophat.sbatch' | xargs -I% sbatch %


###
#RMdups
###
find $outdir -name 'accepted_hits.bam' | awk -F'/' -v sbatchdir=${sbatchdir} '{print "perl /bubo/home/h26/edsgard/glob/code/ngs_pipes/hs_pipe/rmdups.pl", $0, $10, sbatchdir, "b2010035 daniel.edsgard@scilifelab.se";}' >rmdups.sh
sh rmdups.sh
find $sbatchdir -name '*.rmdups.sbatch' | xargs -I% sbatch %


###
#Merge BAMs
###
#Merge
#TBD: As to avoid dep on rnaseq_wf_master.pl: possibly check out: merge/mergebams.samtools.pl, merge/mergebams.pl
perl rnaseq_wf_master.pl -a $projid -id $dirs -s $sampleids -em $email -sbatchdir $sbatchdir -outdir $outdir -ref $ref -pMerge 1
find $sbatchdir -name 'merge.sbatch' | xargs -I% sbatch %

#Create indexes
find $outdir -name 'merged.bam' | xargs -I% echo samtools index % > cmds.sh
#Manually add sbatch header to the cmds.sh script
sbatch cmds.sh


###
#VarCall (SAMtools)
###
#create file listing the bam files
allbamsfile=${outdir}/varcalls/allbams.list
find $outdir -name 'merged.bam' >$allbamsfile

#create regions files
srun -p devel -A b2010035 -t 1:00:00 samtools faidx ${BOWTIE_INDEXES}$GENOME_FASTA &
vcfutils.pl splitchr -l 6600000 ${BOWTIE_INDEXES}${GENOME_FASTA}.fai >${outdir}/varcalls/genome.480.regs

#Run mpileup
cat ${outdir}/varcalls/genome.480.regs | xargs -I% echo 'perl /bubo/home/h26/edsgard/glob/code/ngs_pipes/hs_pipe/mpileup_multisample.pl' $allbamsfile % ${sbatchdir}/varcalls 'b2010035 daniel.edsgard@scilifelab.se' ${BOWTIE_INDEXES}$GENOME_FASTA >cmds.sh
sh cmds.sh
find $sbatchdir -name '*.mpileup.sbatch' | xargs -I% sbatch %
#Check status:
lt *.stderr >allerrfiles.tmp
squeue -u edsgard | grep -v JOBID | awk '{print $1;}' | xargs -I% grep % ${outdir}/varcalls/info/allerrfiles.tmp

#Catenate all region-specific vcfs
find $outdir -name 'allbams.*.vcf' | xargs cat | grep -v '^#' >allbams.vcf.tmp
cat allbams.list.MT:1-16569.vcf | grep '^#' >header.tmp
cat header.tmp allbams.vcf.tmp >allbams.vcf
rm header.tmp allbams.vcf.tmp
wc -l allbams.vcf #989918


###
#TBD
###
#*Should filter variants harder before phasing: DP, Qual (see log.sh in epiase)
#*Add SNP-array calls!!!


###
#Phasing
###

##
#Rm slashes in vcf sample header
##
outdir='/proj/b2011075/analysis/mmseq/outdata/run1'
cd ${outdir}/varcalls
cat allbams.vcf | grep '^#' >header.tmp
cat header.tmp | sed 's/\/proj\/b2011075\/analysis\/mmseq\/outdata\/run1\///g' | sed 's/\/merge\/merged.bam//g' >file.tmp
mv file.tmp header.tmp
grep -v '^#' allbams.vcf >allbams.vcf.noheader
cat header.tmp allbams.vcf.noheader >allbams.vcf

##
#Import SNPs
##
export BOWTIE_INDEXES='/bubo/home/h26/edsgard/glob/annotation/mmseq/' #trailing slash needed
GFF_FILE='Homo_sapiens.GRCh37.64.ref.gff'
outdir='/proj/b2011075/analysis/mmseq/outdata/run1'
(cat <<EOF
#!/bin/bash -l
#SBATCH -A b2010035
#SBATCH -t 1:00:00
#SBATCH -J importsnps
#SBATCH -p devel -n 1
#SBATCH -e /proj/b2011075/analysis/mmseq/outdata/run1/phasing/info/importsnps.jid_%j.stderr
#SBATCH -o /proj/b2011075/analysis/mmseq/outdata/run1/phasing/info/importsnps.jid_%j.stdout
#SBATCH --mail-type=All
#SBATCH --mail-user=$email
module load java/sun_jdk1.6.0_18
export CLASSPATH=$CLASSPATH:/bubo/home/h26/edsgard/opt/polyHap2_20110607/polyHap2.jar
export PATH=$PATH:/bubo/home/h26/edsgard/opt/mmseq_0.9.18
java dataFormat.ImportSNPs ${BOWTIE_INDEXES}$GFF_FILE ${outdir}/varcalls/allbams.vcf ${outdir}/phasing/all
EOF
) >importsnps.sbatch
sbatch importsnps.sbatch


###
#Polyhap
###
execdir='/bubo/home/h26/edsgard/glob/code/ase'
sbatchdir='/proj/b2011075/analysis/mmseq/outscripts/run1'
outdir='/proj/b2011075/analysis/mmseq/outdata/run1'

#Split phase.sh into sbatch jobs
cp ${outdir}/phasing/all/phase.sh ${sbatchdir}/phasing/.
cd ${sbatchdir}/phasing

#Get header
cat phase.sh | grep -v '^java' >phase.sh.header

#Manually edit phase.sh.header into a sbatch header: phase.sh.sbatchheader
cat phase.sh | grep '^java' > phase.sh.noheader
#1 sample 480 jobs: max time was 12min. So 5h should be more than enough when split on 30 jobs.

#Add "--dir $OUTPATH" to every polyhap call.
cat phase.sh.noheader | sed 's/--parentDir/--dir $OUTPATH --parentDir/' | tr '\t' ' ' >file.tmp
mv file.tmp phase.sh.noheader

#Split polyhap calls into njob files
perl ${execdir}/splitfile.pl 480 phase.sh.noheader

#Add header to every split file
find . -name 'phase.sh.noheader.*' | grep -v 'test' | grep -v 'ab2atgc' | grep -v 'sbatch' |  xargs -I% echo cat phase.sh.sbatchheader % '>' %.sbatch >cmds.sh
sh cmds.sh

#Write a specific outdir to each sbatch file
rm cmds.sh
\ls -1 phase.sh.noheader.*.sbatch | grep -v 'ab2atgc' | xargs -I% echo cat % "| sed 's/phasing\/all$/phasing\/all\/"%"/'" '>' %.new >>cmds.sh
sh cmds.sh
 
#Create the specific outdir for every sbatch file
rm cmds.sh
\ls -1 phase.sh.noheader.*.sbatch | grep -v 'ab2atgc' | xargs -I% echo 'mkdir -p' ${outdir}/phasing/all/% >>cmds.sh
sh cmds.sh

#Symlink zip and build files to every specific outdir
rm cmds.sh
find ${outdir}/phasing/all -type d -name '*.sbatch' | xargs -I% echo ln -s ${outdir}/phasing/all/*.zip %/. >>cmds.sh
find ${outdir}/phasing/all -type d -name '*.sbatch' | xargs -I% echo ln -s ${outdir}/phasing/all/build.txt %/. >>cmds.sh
sh cmds.sh

#Submit
find ${sbatchdir}/phasing -name 'phase.sh.noheader.*.sbatch.new' | xargs -I% echo sed 's/40:00/5:00:00/' % '>'%.new >cmds.sh
find ${sbatchdir}/phasing -name 'phase.sh.noheader.*.sbatch.new.new' | xargs -I% sbatch %

#Check:
lt phasing.jid_*.stderr | grep 'Mar 30' >files.tmp
cat files.tmp | awk '{print $9;}' | xargs -I% grep Exc % | sort -u
#Same three errors thrown as before
#java.lang.NullPointerException
#java.lang.RuntimeException: is infinite!
#java.lang.RuntimeException: no snps!


###
#AB to ATGC conversion
###

#Set vars
outdir='/proj/b2011075/analysis/mmseq/outdata/run1'
sbatchdir='/proj/b2011075/analysis/mmseq/outscripts/run1'
execdir='/bubo/home/h26/edsgard/glob/code/ase'
export BOWTIE_INDEXES='/bubo/home/h26/edsgard/glob/annotation/mmseq/' #trailing slash needed
GFF_FILE='Homo_sapiens.GRCh37.64.ref.gff'

#Create specific dirs for each and every ENSG*.out
rm -r ${outdir}/phasing/all/polyhap
find ${outdir}/phasing/all -name 'ENSG*.out' >ensg.tmp
cat ensg.tmp | xargs -I% basename  % | xargs -I% echo mkdir -p ${outdir}/phasing/all/polyhap/% >cmds.sh
sh cmds.sh

#symlink every ensg.out to its own dir
cat ensg.tmp | awk -F'/' -v outdir=${outdir} '{print "ln -s " $0 " " outdir"/phasing/all/polyhap/"$11"/.";}' >cmds.sh
sh cmds.sh

#symlink build.txt to every ensg dir
find ${outdir}/phasing/all/polyhap -name 'ENSG*.out' -type d >ensgdirs.list
cat ensgdirs.list | xargs -I% echo ln -s ${outdir}/phasing/all/build.txt %/. >cmds.sh
sh cmds.sh

#create list of java cmds
cat ensgdirs.list | xargs -I% echo java dataFormat.ConvertABtoATGC ${BOWTIE_INDEXES}$GFF_FILE % %/atgc >ab2atgc.sh

#create sbatch jobs
perl ${execdir}/splitfile.pl 480 ab2atgc.sh
(cat <<EOF
#!/bin/bash -l
#SBATCH -A b2010035
#SBATCH -t 2:00:00
#SBATCH -J ab2atgc
#SBATCH -p core -n 1
#SBATCH -e /proj/b2011075/analysis/mmseq/outdata/run1/phasing/all/polyhap/info/ab2atgc.jid_%j.stderr
#SBATCH -o /proj/b2011075/analysis/mmseq/outdata/run1/phasing/all/polyhap/info/ab2atgc.jid_%j.stdout
#SBATCH --mail-type=All
#SBATCH --mail-user=$email
module load java/sun_jdk1.6.0_18
export CLASSPATH=$CLASSPATH:/bubo/home/h26/edsgard/opt/polyHap2_20110607/polyHap2.jar
export PATH=$PATH:/bubo/home/h26/edsgard/opt/mmseq_0.9.18
EOF
) >sbatch.template
ls -1 ab2atgc.sh.* | xargs -I% echo cat sbatch.template % '>' %.sbatch >cmds.sh
sh cmds.sh
find . -name 'ab2atgc.sh.*.sbatch' | xargs -I% sbatch %

#Catenate .hap files
\ls -1 ${outdir} | grep '12h' >samples.list
find ${outdir}/phasing/all/polyhap -name '*.hap' >hap.files
cat samples.list | xargs -I% echo grep % 'hap.files | xargs -I@ cat @ >>' ${outdir}/phasing/%.hap >cmds.sh
sh cmds.sh

#Catenate SNPinfo files
find ${outdir}/phasing/all/polyhap -name 'SNPinfo' >snpinfo.files
cat snpinfo.files | xargs -I% cat % >${outdir}/phasing/SNPinfo &

#CHECK
cat ${outdir}/phasing/SNPinfo | wc -l #16931
wc -l ${outdir}/phasing/*.hap #33862
cat ${outdir}/phasing/4_unstim_12h.hap | awk -F'\t' '{print $1;}' | grep '_A' | sed 's/_A//' >hap.enst
cat ${outdir}/phasing/SNPinfo | awk -F'\t' '{print $1;}' >snpinfo.enst
comm -3 hap.enst snpinfo.enst


###
#Create custom fasta transcriptomes
###
#Set vars
outdir='/proj/b2011075/analysis/mmseq/outdata/run1'
sbatchdir='/proj/b2011075/analysis/mmseq/outscripts/run1'
execdir='/bubo/home/h26/edsgard/glob/code/ase/mmseq'
export BOWTIE_INDEXES='/bubo/home/h26/edsgard/glob/annotation/mmseq/' #trailing slash needed
GFF_FILE='Homo_sapiens.GRCh37.64.ref.gff'
TRANSCRIPT_FASTA='Homo_sapiens.GRCh37.64.ref_transcripts.fa'

#Ensure that SNPinfo and .hap is compatible
#find ${outdir}/phasing -maxdepth 1 -name '*12h.hap' -type f | xargs -I% basename % | sed 's/.hap//' >sample.list
#cat sample.list | xargs -I% echo sh ${execdir}/clean_haplorefinput.sh % ${outdir}/phasing >cmds.sh
#srun -p devel -A b2010035 -t 1:00:00 sh cmds.sh &
#haploref.rb ${BOWTIE_INDEXES}$TRANSCRIPT_FASTA ${BOWTIE_INDEXES}$GFF_FILE sample.SNPinfo.shared sample.hap.shared > sample.fa

#Generate haploref sbatch scripts
(cat <<EOF
#!/bin/bash -l
#SBATCH -A b2010035
#SBATCH -t 55:00
#SBATCH -J haploref
#SBATCH -p core -n 1
#SBATCH -e /proj/b2011075/analysis/mmseq/outdata/run1/phasing/info/haploref.jid_%j.stderr
#SBATCH -o /proj/b2011075/analysis/mmseq/outdata/run1/phasing/info/haploref.jid_%j.stdout
#SBATCH --mail-type=All
#SBATCH --mail-user=$email
export PATH=$PATH:/bubo/home/h26/edsgard/opt/mmseq_0.10.0
cd ${outdir}/phasing
haploref.rb ${BOWTIE_INDEXES}$TRANSCRIPT_FASTA ${BOWTIE_INDEXES}$GFF_FILE SNPinfo sample.hap > sample.fa
EOF
) >sbatch.template
cat sample.list | xargs -I% echo cat sbatch.template "| sed 's/sample/"%"/g' >" %.haploref.sbatch >cmds.sh
sh cmds.sh
find . -name '*.haploref.sbatch' | xargs -I% sbatch %

#CHECK
cat ${outdir}/phasing/info/haploref* | sort -u
cat ${outdir}/phasing/2_unstim_12h.fa | grep '^>' | grep '_A' | sort -u | wc -l #7622


###
#Create custom transcriptome bowtie indexes
###
mkdir -p ${outdir}/bowtie/info
find ${outdir}/phasing -maxdepth 1 -name '*.fa' | awk -F'/' -v outdir=${outdir} '{print "ln -s", $0, outdir"/bowtie/"$9;}' >cmds.sh
sh cmds.sh

(cat <<EOF
#!/bin/bash -l
#SBATCH -A b2010035
#SBATCH -t 5:00:00
#SBATCH -J bowtie-build
#SBATCH -p core -n 1
#SBATCH -e /proj/b2011075/analysis/mmseq/outdata/run1/bowtie/info/bowtie-build.jid_%j.stderr
#SBATCH -o /proj/b2011075/analysis/mmseq/outdata/run1/bowtie/info/bowtie-build.jid_%j.stdout
#SBATCH --mail-type=All
#SBATCH --mail-user=$email
module load bioinfo-tools
module load bowtie/0.12.7
cd ${outdir}/bowtie
bowtie-build -f sample.fa sample
EOF
) >sbatch.template
cat sample.list | xargs -I% echo cat sbatch.template "| sed 's/sample/"%"/g' >" %.bowtiebuild.sbatch >cmds.sh
sh cmds.sh
find . -name '*.bowtiebuild.sbatch' | xargs -I% sbatch %


###
#Mapping to custom transcriptomes
###

#Memory requirements: Bowtie uses approximately as much memory as the size of the bowtie indices.
#du -h *.ebwt | awk '{print $1;}' | sed 's/M//' | awk 'BEGIN{a=0}{a=a+$1;}END{print a;}' #5.8G / 16 => <1GB => Can run on core.

#Set vars
sbatchdir='/proj/b2011075/analysis/mmseq/outscripts/run1'
outdir='/proj/b2011075/analysis/mmseq/outdata/run1'
RESULTS='/proj/b2011075/analysis/mmseq/outdata/run1/bowtie'
CHUNKMBS=512
THREADS=1
execdir='/bubo/home/h26/edsgard/glob/code/ase/mmseq'

#Put all fastq files in one dir
find '/proj/b2011075/analysis/hs_pipe/outdata/run1' -type f -name '*.fastq.[12].filter.fastq' >fastq.files
find '/proj/b2011075/analysis/hs_pipe/outdata/run1' -type l -name '*.fastq.[12].filter.fastq' >>fastq.files
cat fastq.files | xargs -I% echo ln -s % /proj/b2011075/analysis/mmseq/outdata/run1/fastq/. >cmds.sh
sh cmds.sh

#For each fastq file write a sbatch file
find ${outdir}/fastq -name '*1.filter.fastq' | xargs -I% echo perl ${execdir}/bowtie.gen.sbatch.pl % $THREADS $CHUNKMBS $RESULTS >cmds.sh
sh cmds.sh

#Set insert sizes
cat /proj/b2011075/data/rnaseq/bioanalyzer/isizes.tab | awk -F'\t' '{print $1, $4;}'
#Round these to nearest ceiling-hundred, but at least 400: That is, 400 for all apart from: 7_LPS: 600, 9_LPS: 500
find '/proj/b2011075/analysis/mmseq/outscripts/run1/bowtie' -name '*.fastq.sbatch' | grep 7_LPS | xargs -I% echo cat % " | sed s/400/600/ >" %.new >cmds.sh
sh cmds.sh
find '/proj/b2011075/analysis/mmseq/outscripts/run1/bowtie' -name '*.fastq.sbatch' | grep 9_LPS | xargs -I% echo cat % " | sed s/400/500/ >" %.new >cmds.sh
sh cmds.sh
rename fastq.sbatch.new fastq.sbatch *.fastq.sbatch.new

#Submit
find '/proj/b2011075/analysis/mmseq/outscripts/run1/bowtie' -name '*fastq.sbatch' | xargs -I% sbatch %


###
#RMdups
###
#12G bam takes ~4.5h
sbatchdir='/proj/b2011075/analysis/mmseq/outscripts/run1/bowtie'
execdir='/bubo/home/h26/edsgard/glob/code/ase/mmseq'
bamdir='/proj/b2011075/analysis/mmseq/outdata/run1/bowtie'
outdir=$bamdir
cd $sbatchdir
find $bamdir -name '*.bam' > bamfiles.list
cat bamfiles.list | xargs -I% echo perl ${execdir}'/rmdups.pl' % $outdir $sbatchdir 'b2010035 daniel.edsgard@scilifelab.se' >cmds.sh
sh cmds.sh
find . -name '*.rmdups.sbatch' | xargs -I% sbatch %

#Check:
grep 'Exc' ${outdir}/info/*rmdups.stderr


###
#Sort with samtools (NB: by read-name)
###
sbatchdir='/proj/b2011075/analysis/mmseq/outscripts/run1/bowtie'
execdir='/bubo/home/h26/edsgard/glob/code/ase/mmseq'
bamdir='/proj/b2011075/analysis/mmseq/outdata/run1/bowtie'
outdir='/proj/b2011075/analysis/mmseq/outdata/run1/bowtie'
cd $sbatchdir
find $bamdir -name '*bam.nodup' >bamfiles.nodup.list
cat bamfiles.nodup.list | xargs -I% echo perl ${execdir}/'sortbams.samtools.pl' % $outdir $sbatchdir 'b2010035 daniel.edsgard@scilifelab.se' >cmds.sh
sh cmds.sh
find . -name '*.samtoolssort.sbatch' | xargs -I% sbatch %

#Check:
cat ${bamdir}/info/*.samtoolssort.stderr | sort -u


###
#Merge BAMs. NB: Sorted by read names rather than coordinate.
###
sbatchdir='/proj/b2011075/analysis/mmseq/outscripts/run1/bowtie'
execdir='/bubo/home/h26/edsgard/glob/code/ase/mmseq'
bamdir='/proj/b2011075/analysis/mmseq/outdata/run1/bowtie'
outdir='/proj/b2011075/analysis/mmseq/outdata/run1/bowtie'
cd $sbatchdir
find $bamdir -name '*nodup.sorted.bam' | sed 's/\..*//' | sort -u >samples.list
cat samples.list | xargs -I% basename % | xargs -I% echo perl ${execdir}'/mergebams.samtools.pl' % $bamdir $outdir $sbatchdir 'b2010035 daniel.edsgard@scilifelab.se' >cmds.sh
sh cmds.sh
find $sbatchdir -name '*.mergesams.sbatch' | xargs -I% sbatch %


#Check
#total size of bams:
lt -h ${bamdir}/*.nodup.sorted.bam | awk '{print $5;}' | sed 's/G//' | awk 'BEGIN{a=0;}{a=a+$1;}END{print a;}'


###
#BAM2HITS
###
sbatchdir='/proj/b2011075/analysis/mmseq/outscripts/run1/mmseq'
execdir='/bubo/home/h26/edsgard/glob/code/ase/mmseq'
bamdir='/proj/b2011075/analysis/mmseq/outdata/run1/bowtie'
fastadir='/proj/b2011075/analysis/mmseq/outdata/run1/bowtie'
bam2hitsdir='/proj/b2011075/analysis/mmseq/outdata/run1/mmseq'
isizedir='/proj/b2011075/data/rnaseq/bioanalyzer'
isizefile=${isizedir}/'isizes.tmp'

cd $sbatchdir
cat ${isizedir}/isizes.tab | awk -F'\t' '{print $1, $3;}' | grep -v isize > $isizefile
find $bamdir -name '*_12h.bam' | sort -u >bamfiles.sorted.list
#ls -1 $bam2hitsdir/*.bam.hits | xargs -I% basename % | xargs -I % echo '/proj/b2011075/analysis/mmseq/outdata/run1/bowtie/'% | sed 's/.hits//' | sort -u >submitted.bamhits
#comm -3 bamfiles.sorted.list submitted.bamhits >bams2submit.list
cat bams2submit.list | xargs -I% echo perl ${execdir}/'bam2hits.gensbatch.pl' % $fastadir $bam2hitsdir $sbatchdir $isizefile >cmds.sh
sh cmds.sh
find $sbatchdir -name '*.bam2hits.sbatch' | grep 8_LPS | xargs -I% sbatch %

#Check
cat $bam2hitsdir/info/bam2hits.*.stderr | sort -u
#TBD: Getting read lengths is wrong. The bam2hits script apparently doesnt handle varying read length...


###
#MMSEQ
###
#~7h10min was the longest run.

bam2hitsdir='/proj/b2011075/analysis/mmseq/outdata/run1/mmseq'
find $bam2hitsdir -name '*_12h.bam.hits' | sort >bam2hitsfiles.list
(cat <<EOF
#!/bin/bash -l
#SBATCH -A b2010035
#SBATCH -t 11:00:00
#SBATCH -J MMSEQ
#SBATCH -p core -n 1
#SBATCH -e ${bam2hitsdir}/info/mmseq.jid_%j.stderr
#SBATCH -o ${bam2hitsdir}/info/mmseq.jid_%j.stdout
#SBATCH --mail-type=All
#SBATCH --mail-user=$email
export PATH=$PATH:/bubo/home/h26/edsgard/opt/mmseq_0.10.0
cd ${bam2hitsdir}
mmseq sample sample.out
EOF
) >sbatch.template
cat bam2hitsfiles.list | xargs -I% basename % | xargs -I% echo cat sbatch.template "| sed 's/sample/"%"/g' >" %.mmseq.sbatch >cmds.sh
sh cmds.sh
find . -name '*.mmseq.sbatch' | xargs -I% sbatch %


#Check:
cat 7_unstim_12h.bam.hits.out.mmseq | awk '{print $1;}' | sed 's/_[AB]$//' | sort -u | wc -l #169259
cat 7_unstim_12h.bam.hits.out.mmseq | awk '{print $1;}' | grep '_A$' | sort -u | wc -l #3998
cat 7_unstim_12h.bam.hits.out.gene.mmseq | awk -F'\t' '{print $1;}' | sort -u | wc -l #54030


###
#Test for ASE
###
#Binomial test... See mmseq.stats.R, basecount.sh, snv.ase.R, skelly.sh, skelly.R.

#Get annot
ensg2enst2hugo='/Users/danieledsgard/Dropbox/postdoc/db/ensg2enst2hugo.v66'
wget -O $ensg2enst2hugo 'http://www.biomart.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "hgnc_symbol" /><Attribute name = "description" /></Dataset></Query>'

#Map chr pos to ensg pos
execdir='/bubo/home/h26/edsgard/glob/code/ase/mmseq'
datadir='/proj/b2011075/analysis/mmseq/outdata/run1/phasing'
annotdir='/bubo/home/h26/edsgard/glob/annotation/mmseq'
REF=${annotdir}'/Homo_sapiens.GRCh37.64.ref_transcripts.fa'
GFF=${annotdir}'/Homo_sapiens.GRCh37.64.ref.gff'
SNPS=${datadir}'/SNPinfo'
cd $datadir
ruby ${execdir}/map.chrpos2featpos.rb $REF $GFF $SNPS >chrpos2relpos.tab




