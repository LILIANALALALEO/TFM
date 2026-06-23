cd ~/TFM

wget -P Data/Reference_genome \
https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-53/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz

gunzip -f Data/Reference_genome/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz

salmon index \
-t Data/Reference_genome/Arabidopsis_thaliana.TAIR10.cdna.all.fa \
-i Data/Reference_genome/salmon_index \
-k 31
