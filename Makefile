WHERE = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/scRNAseq
REF = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/scRNAseq/ref/refdata-cellranger-GRCh38-1.2.0
INDEX = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/hg38index
WHOLEGENOME = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
GTF = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf
FC = /home/sa941/software/subread-1.6.3-Linux-i386/bin
ANNOTATIONS = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations/Homo_sapiens/NCBI/GRCh38/Annotation/Genes



Homo_sapiens.GRCh38.93.gtf.gz: 
	wget ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz
	wget ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz



#---------------------------------UMI Tools-----------------------------------------------------------------
#-----------Sample D
d.whitelist.txt: D.R1.L002.fastq.gz 
	umi_tools whitelist --stdin ${WHERE}/D.R1.L002.fastq.gz \
	--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
	--plot-prefix Dplot \
	--log2stderr > d.whitelist.txt
                                                                                 
d.R2.extracted.fastq.gz: d.whitelist.txt 
	umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
	--stdin ${WHERE}/D.R1.L002.fastq.gz \
	--stdout ${WHERE}/D.R1.extracted.fastq.gz \
	--read2-in ${WHERE}/D.R2.L002.fastq.gz \
	--read2-out ${WHERE}/D.R2.extracted.fastq.gz \
	--filter-cell-barcode \
	--whitelist=d.whitelist.txt 

d.Aligned.sortedByCoord.out.bam: D.R2.extracted.fastq.gz 
	STAR --runThreadN 4 --genomeDir ${INDEX} --readFilesIn  ${WHERE}/D.R2.extracted.fastq.gz --readFilesCommand zcat --outFilterIntronMotifs None --runMode alignReads --outSAMtype BAM SortedByCoordinate --outFileNamePrefix d. --outFilterMultimapNmax 1 
 
d.gene_assigned: d.Aligned.sortedByCoord.out.bam
	${FC}/featureCounts -a  ${ANNOTATIONS}/genes.gtf \
	-o d.gene_assigned \
	-R BAM d.Aligned.sortedByCoord.out.bam \
	-T 4
	samtools sort d.Aligned.sortedByCoord.out.bam.featureCounts.bam -o d.assigned_sorted.bam
	samtools index d.assigned_sorted.bam
              
d.counts.tsv.gz: d.assigned_sorted.bam
	umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --paired --wide-format-cell-counts -I d.assigned_sorted.bam -S d.counts.tsv.gz



#--------------Sample E, we can loop if we have larger numbers of samples
e.whitelist.txt: E.R1.L002.fastq.gz
	umi_tools whitelist --stdin ${WHERE}/E.R1.L002.fastq.gz \
	--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
	--plot-prefix Eplot \
	--log2stderr > e.whitelist.txt

e.R2.extracted.fastq.gz: e.whitelist.txt
	umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
	--stdin ${WHERE}/E.R1.L002.fastq.gz \
	--stdout ${WHERE}/E.R1.extracted.fastq.gz \
	--read2-in ${WHERE}/E.R2.L002.fastq.gz \
	--read2-out ${WHERE}/E.R2.extracted.fastq.gz \
	--filter-cell-barcode \
	--whitelist=e.whitelist.txt 

e.Aligned.sortedByCoord.out.bam: E.R2.extracted.fastq.gz
	STAR --runThreadN 4 --genomeDir ${INDEX} --readFilesIn ${WHERE}/E.R2.extracted.fastq.gz --readFilesCommand zcat --outFilterIntronMotifs None --runMode alignReads --outSAMtype BAM SortedByCoordinate --outFileNamePrefix e. --outFilterMultimapNmax 1 

e.gene_assigned: e.Aligned.sortedByCoord.out.bam
	${FC}/featureCounts -a  ${ANNOTATIONS}/genes.gtf \
	-o e.gene_assigned \
	-R BAM e.Aligned.sortedByCoord.out.bam \
	-T 4
	samtools sort e.Aligned.sortedByCoord.out.bam.featureCounts.bam -o e.assigned_sorted.bam
	samtools index e.assigned_sorted.bam

e.counts.tsv.gz: e.assigned_sorted.bam 
	umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --paired --wide-format-cell-counts -I e.assigned_sorted.bam -S e.counts.tsv.gz


#--------------Sample F 
f.whitelist.txt: F.R1.L002.fastq.gz
	umi_tools whitelist --stdin ${WHERE}/F.R1.L002.fastq.gz \
	--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
	--plot-prefix Fplot \
	--log2stderr > f.whitelist.txt

f.R2.extracted.fastq.gz: f.whitelist.txt
	umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
	--stdin ${WHERE}/F.R1.L002.fastq.gz \
	--stdout ${WHERE}/F.R1.extracted.fastq.gz \
	--read2-in ${WHERE}/F.R2.L002.fastq.gz \
	--read2-out ${WHERE}/F.R2.extracted.fastq.gz \
	--filter-cell-barcode \
	--whitelist=f.whitelist.txt 
 
f.Aligned.sortedByCoord.out.bam: F.R2.extracted.fastq.gz
	STAR --runThreadN 4 --genomeDir ${INDEX} --readFilesIn ${WHERE}/F.R2.extracted.fastq.gz --readFilesCommand zcat --outFilterIntronMotifs None --runMode alignReads --outSAMtype BAM SortedByCoordinate --outFileNamePrefix f. --outFilterMultimapNmax 1 
f.gene_assigned: f.Aligned.sortedByCoord.out.bam
	${FC}/featureCounts -a  ${ANNOTATIONS}/genes.gtf \
	-o f.gene_assigned \
	-R BAM f.Aligned.sortedByCoord.out.bam \
	-T 4
	samtools sort f.Aligned.sortedByCoord.out.bam.featureCounts.bam -o f.assigned_sorted.bam
	samtools index f.assigned_sorted.bam

f.counts.tsv.gz: f.assigned_sorted.bam
	umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --paired --wide-format-cell-counts -I f.assigned_sorted.bam -S f.counts.tsv.gz



#--------------Sample G
g.whitelist.txt: G.R1.L002.fastq.gz
	umi_tools whitelist --stdin ${WHERE}/G.R1.L002.fastq.gz \
	--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
	--plot-prefix Gplot \
	--log2stderr > g.whitelist.txt

g.R2.extracted.fastq.gz: g.whitelist.txt 
	umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
	--stdin ${WHERE}/G.R1.L002.fastq.gz \
	--stdout ${WHERE}/G.R1.extracted.fastq.gz \
	--read2-in ${WHERE}/G.R2.L002.fastq.gz \
	--read2-out ${WHERE}/G.R2.extracted.fastq.gz \
	--filter-cell-barcode \
	--whitelist=g.whitelist.txt 

g.Aligned.sortedByCoord.out.bam: G.R2.extracted.fastq.gz
	STAR --runThreadN 4 --genomeDir ${INDEX} --readFilesIn ${WHERE}/G.R2.extracted.fastq.gz --readFilesCommand zcat --outFilterIntronMotifs None --runMode alignReads --outSAMtype BAM SortedByCoordinate --outFileNamePrefix g. --outFilterMultimapNmax 1 
g.gene_assigned: g.Aligned.sortedByCoord.out.bam
	${FC}/featureCounts -a  ${ANNOTATIONS}/genes.gtf \
	-o g.gene_assigned \
	-R BAM g.Aligned.sortedByCoord.out.bam \
	-T 4
	samtools sort g.Aligned.sortedByCoord.out.bam.featureCounts.bam -o g.assigned_sorted.bam
	samtools index g.assigned_sorted.bam

g.counts.tsv.gz: g.assigned_sorted.bam 
	umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --paired --wide-format-cell-counts  -I g.assigned_sorted.bam -S g.counts.tsv.gz


g.counts.csv: 
	gunzip -c d.counts.tsv.gz > d.counts.tsv 
	gunzip -c e.counts.tsv.gz > e.counts.tsv 
	gunzip -c f.counts.tsv.gz > f.counts.tsv 
	gunzip -c g.counts.tsv.gz > g.counts.tsv 













 
