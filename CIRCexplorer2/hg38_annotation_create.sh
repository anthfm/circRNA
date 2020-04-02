###1. Download human KnownGenes (GENCODE) gene annotation file using fetch_ucsc.py script from authors
fetch_ucsc.py hg38 kg hg38_kg.txt


###2. Convert gene annotation file to GTF format (require genePredToGtf)
#to get genePredtoGtf do the following:
rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/genePredToGtf ./

cut -f2-11 hg38_kg.txt|genePredToGtf file stdin hg38_kg.gtf


###3. Download human reference genome (hg38) - UCSC
fetch_ucsc.py hg38 fa hg38.fa
