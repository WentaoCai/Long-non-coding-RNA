# Long-non-coding-RNA

1.Hisat2 build genome index：

First, using the python scripts included in the HISAT2 package, extract splice-site and exon information from the gene
annotation file:

$ extract_splice_sites.py gemome.gtf >genome.ss# Get splicing information
$ extract_exons.py genome.gtf >genome.exon#Get exon information

Second, build a HISAT2 index:

$ hisat2-build --ss genome.ss --exon genome.exon genome.fa genome

you can directly use the following command:
$ hisat2-build  genome.fa genome

2. Mapped to genome：

hisat2 -p 8 --dta -x indexes/genome -1 file1_1.fastq.gz -2 file1_2.fastq.gz -S file1.sam
hisat2 -p 8 --dta -x indexes/genome -1 file2_1.fastq.gz -2 file2_2.fastq.gz -S file2.sam



3. Sort sam to bam:

$ samtools sort -@ 8 -o file1.bam file1.sam
$ samtools sort -@ 8 -o file2.bam file2.sam

4. Assembly：

$ stringtie -p 8 -G genome.gtf -o file1.gtf –l file1 file1.bam
$ stringtie -p 8 -G genome.gtf -o file2.gtf –l file2 file2.bam
lncRNA (-f 0.01 -a 10 -j 1 -c 0.01)

5. Merge gtf

$ stringtie --merge -p 8 -G genome.gtf -o stringtie_merged.gtf mergelist.txt

6. get novel transcripts for lncRNA

gffcompare –r genomegtf –G –o merged stringtie_merged.gtf

gffcompare download from http://ccb.jhu.edu/software/stringtie/gff.shtml

7. Qualification：

$ stringtie –e –B -p 8 -G stringtie_merged.gtf -o ballgown/file1/file1.gtf file1.bam
$ stringtie –e –B -p 8 -G stringtie_merged.gtf -o ballgown/file2/file2.gtf file2.bam

8. Ballgown DGEs analysis：

>library(ballgown)
>library(RSkittleBrewer)
>library(genefilter)
>library(dplyr)
>library(devtools)
>pheno_data = read.csv("geuvadis_phenodata.csv")#Read phenotype
>bg_chrX = ballgown(dataDir = "ballgown", samplePattern = "file", pData=pheno_data)#Read expression
>bg_chrX_filt = subset(bg_chrX,"rowVars(texpr(bg_chrX)) >1",genomesubset=TRUE)#filtering low expression genes
>results_transcripts = stattest(bg_chrX_filt,feature="transcript",covariate="sex",adjustvars =c("population"), getFC=TRUE, meas="FPKM")#compare group:sex，factor：population
>results_genes = stattest(bg_chrX_filt, feature="gene",covariate="sex", adjustvars = c("population"), getFC=TRUE,meas="FPKM")
>results_transcripts=data.frame(geneNames=ballgown::geneNames(bg_chrX_filt),geneIDs=ballgown::geneIDs(bg_chrX_filt), results_transcripts)#Add gene names and id
>results_transcripts = arrange(results_transcripts,pval)#按pval sort
>results_genes = arrange(results_genes,pval)
>write.csv(results_transcripts, "chrX_transcript_results.csv",
row.names=FALSE)
>write.csv(results_genes, "chrX_gene_results.csv",
row.names=FALSE)
>subset(results_transcripts,results_transcripts$qval<0.05)
>subset(results_genes,results_genes$qval<0.05)

9. Visualization：

>tropical= c('darkorange', 'dodgerblue',
'hotpink', 'limegreen', 'yellow')
>palette(tropical)
>fpkm = texpr(bg_chrX,meas="FPKM")
>fpkm = log2(fpkm+1)
>boxplot(fpkm,col=as.numeric(pheno_data$sex),las=2,ylab='log2(FPKM+1)')
>ballgown::transcriptNames(bg_chrX)[12]
## 12
## "NM_012227"
>ballgown::geneNames(bg_chrX)[12]
## 12
## "GTPBP6"
>plot(fpkm[12,] ~ pheno_data$sex, border=c(1,2),
main=paste(ballgown::geneNames(bg_chrX)[12],' : ',
ballgown::transcriptNames(bg_chrX)[12]),pch=19, xlab="Sex",
ylab='log2(FPKM+1)')
>points(fpkm[12,] ~ jitter(as.numeric(pheno_data$sex)),
col=as.numeric(pheno_data$sex))
>plotTranscripts(ballgown::geneIDs(bg_chrX)[1729], bg_chrX, main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))
>plotMeans('MSTRG.56', bg_chrX_filt,groupvar="sex",legend=FALSE)
