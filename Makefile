WHERE = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/SC-SLX-18019
REF = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/scRNAseq/ref/refdata-cellranger-GRCh38-1.2.0
INDEX = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/hg38index
WHOLEGENOME = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
GTF = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf
FC = /home/sa941/software/subread-1.6.3-Linux-i386/bin
ANNOTATIONS = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations/Homo_sapiens/NCBI/GRCh38/Annotation/Genes



Homo_sapiens.GRCh38.93.gtf.gz: 
	wget ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz
	wget ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz


#Dataset SLX-18019.H7HVMDRXX

#---------------------------------UMI Tools-----------------------------------------------------------------
SAMPLES = D E F G 

G.whitelist.txt:
	 $(foreach i, $(SAMPLES), umi_tools whitelist --stdin $(i).R1.L002.fastq.gz --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --plot-prefix $(i)plot --log2stderr > $(i).whitelist.txt;) 
G.R2.extracted.fastq.gz: 
	$(foreach i, $(SAMPLES),  umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN  --stdin $(i).R1.L002.fastq.gz --stdout $(i).R1.extracted.fastq.gz --read2-in $(i).R2.L002.fastq.gz --read2-out $(i).R2.extracted.fastq.gz --filter-cell-barcode --whitelist=$(i).whitelist.txt;)

G.Aligned.sortedByCoord.out.bam:
	 $(foreach i, $(SAMPLES), STAR --runThreadN 4 --genomeDir ${INDEX} --readFilesIn $(i).R2.extracted.fastq.gz --readFilesCommand zcat --outFilterIntronMotifs None --runMode alignReads --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $(i). --outFilterMultimapNmax 1 ;) 
 
G.gene_assigned: 
	$(foreach i, $(SAMPLES), ${FC}/featureCounts -a  ${ANNOTATIONS}/genes.gtf -o $(i).gene_assigned -R BAM $(i).Aligned.sortedByCoord.out.bam -T 4;)
	$(foreach i, $(SAMPLES), samtools sort $(i).Aligned.sortedByCoord.out.bam.featureCounts.bam -o $(i).assigned_sorted.bam;)
	$(foreach i, $(SAMPLES), samtools index $(i).assigned_sorted.bam;)
             
G.counts.tsv.gz: 
	$(foreach i, $(SAMPLES), umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --paired --wide-format-cell-counts -I $(i).assigned_sorted.bam -S $(i).counts.tsv.gz;) 

G.counts.csv: 
	$(foreach i, $(SAMPLES), gunzip -c $(i).counts.tsv.gz > $(i).counts.tsv;) 
	$(foreach i, $(SAMPLES), Rscript seurat.R ${WHERE} $(i).counts.tsv;)












 
