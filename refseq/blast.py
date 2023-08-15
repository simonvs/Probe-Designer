import subprocess
from Bio import SearchIO

db1_name = "db1_genoma"
db2_name = "db2_transcriptoma"
query = "query.txt"

blastn_cmd = ['blastn',
              '-task', 'blastn-short',
              '-query', 'query.txt',
              '-db', db1_name,
              '-out', 'resultados.xml',
              '-outfmt', '5',
              '-word_size', '11',
              '-dust', 'no']

subprocess.run(blastn_cmd, shell=True)

blast_qresult = SearchIO.read("resultados.xml", "blast-xml")[:3]

print(blast_qresult[0].hsps[0].evalue)
print(blast_qresult[0].hsps[1].evalue)
print(blast_qresult[0].hsps[2].evalue)
print(blast_qresult[1].hsps[0].evalue)
print(blast_qresult[1].hsps[1].evalue)
print(blast_qresult[1].hsps[2].evalue)
