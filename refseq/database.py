import subprocess
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq

# Rutas a archivos
fasta_genoma = "D:/ADNhumano/ncbi-genomes-2023-06-30/GCF_000001405.40_GRCh38.p14_genomic.fna"
fasta_transcriptoma =  "gencode.v43.transcripts.fa"
db1_name = "db1_genoma"  # Nombre de la base de datos
db2_name = "db2_transcriptoma"
# Crear la base de datos local
makeblastdb_cmd = ['makeblastdb',
                   '-in', fasta_genoma,
                   '-dbtype', 'nucl',
                   '-out', db1_name]
subprocess.run(makeblastdb_cmd, shell=True)
