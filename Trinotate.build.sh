Trinotate_home=$(dirname $(readlink -f Trinotate))
cd $Trinotate_home
$Trinotate_home/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
makeblastdb -in uniprot_sprot.pep -dbtype prot
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
