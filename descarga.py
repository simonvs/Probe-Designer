import os
import argparse
from Bio import Entrez, SeqIO
Entrez.mail = "simonvergaraswett@hotmail.com"

def accnum_to_seqrecord(accesion_number):
    # Función que recibe un accession number de la base de datos 'nucleotide' de NCBI
    # https://www.ncbi.nlm.nih.gov/nuccore.
    # Genera el archivo .gbk (GenBank) en la carpeta 'files' y retorna la secuencia en formato SeqRecord.
    path = "files/"+accesion_number+".gbk"
    if not os.path.isfile(path):
        print("Descargando secuencia...")
        with Entrez.efetch(db="nucleotide", id=accesion_number, rettype="gb", retmode="text") as handle:
            with open(path, "w") as out_handle:
                out_handle.write(handle.read())
    record = SeqIO.read(path, "genbank")
    return record

def gbk_to_seqrecord(path_to_gbk):
    # Función que recibe una ruta a un archivo .gbk o .gb (GenBank) y lo retorna como SeqRecord.
    record = SeqIO.read(path_to_gbk, "genbank")
    return record

def fasta_to_seqrecord(path_to_fasta):
    # Función que recibe una ruta a un archivo FASTA y lo retorna como SeqRecord.
    record = SeqIO.read(path_to_fasta, "fasta")
    return record

def main():    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-an',
                        '--accessionnumber',
                        help='Accesion number de la secuencia a descargar.',
                        required=True,
                        default=argparse.SUPPRESS)

    args = parser.parse_args()
    dargs = vars(args)

    record = accnum_to_seqrecord(accesion_number=dargs["accessionnumber"])
    print(record)

if __name__ == "__main__":
    main()