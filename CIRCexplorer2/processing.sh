###1### align reads that span to mulitple exons with TopHat-Fusion (because we have paired end data, CIRCexplorer reccomends running TopHat-Fusion, instead of tophat first and then using the unmapped reads (multi-spanning exons) with TopHat-fusion for single-end reads

#32gb memory
#15 threads
#24hr run time

module load tophat2
module load bowtie


bowtie1_index="/cluster/projects/bhklab/Data/ncRNA_detect/ref/Homo_sapiens/UCSC/hg38/Sequence/BowtieIndex/genome"
f1="sample_1_1.rnaseq.fastq.gz"
f2="sample_1_2.rnaseq.fastq.gz"
output="/cluster/projects/bhklab/Data/ncRNA_detect/circRNA/CIRCexplorer2/HCT-116/gCSI"

tophat2 -o $output -p 15 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search $bowtie1_index $f1 $f2

###2### parse back-spliced junctions together into organized format (BED file), before annotation
#10gb memory
#6 threads
#5hr run time
module load CIRCexplorer2

hits="/cluster/projects/bhklab/Data/ncRNA_detect/circRNA/CIRCexplorer2/HCT-116/gCSI/accepted_hits.bam"
outdir="/cluster/projects/bhklab/Data/ncRNA_detect/circRNA/CIRCexplorer2/HCT-116/gCSI/tophat2_only_hct116/parse"
mkdir -p $outdir

CIRCexplorer2 parse --pe -t TopHat-Fusion $hits > $outdir/CIRCexplorer2_parse.log

###3### annotate back-spliced junctions & circRNA/ciRNA's expressed

module load CIRCexplorer2


bed="back_spliced_junction.bed"
outdir="/cluster/projects/bhklab/Data/ncRNA_detect/circRNA/CIRCexplorer2/HCT-116/gCSI/annotate"
gtf="/cluster/projects/bhklab/Data/ncRNA_detect/circRNA/CIRCexplorer2/HCT-116/gCSI/gencode_v33_refFlat.txt"
genome="/cluster/projects/bhklab/Data/ncRNA_detect/ref/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"

mkdir -p $outdir
CIRCexplorer2 annotate --low-confidence -r $gtf -g $genome -b $bed -o circularRNA_known.txt > $outdir/CIRCexplorer2_annotate_parse.log


