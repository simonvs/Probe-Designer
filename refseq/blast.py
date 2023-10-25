import subprocess
from Bio import SearchIO

db1_name = "db1_genoma"
db2_name = "db2_transcriptoma"
query = "refseq/query.txt"

blastn_cmd = ['blastn',
              '-task', 'blastn-short',
              '-query', 'refseq/query.txt',
              '-db', db2_name,
              '-out', 'resultados.xml',
              '-outfmt', '5',
              '-word_size', '11',
              '-dust', 'no']

subprocess.run(blastn_cmd, shell=True)

blast_qresult = SearchIO.read("resultados.xml", "blast-xml")[:3]

print(blast_qresult[0].hsps[0].evalue)

