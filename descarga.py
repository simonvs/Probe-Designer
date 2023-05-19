import os
import argparse
import shutil
from Bio import Entrez, SeqIO

def accnum_to_seqrecord(accesion_number):
    # Función que recibe un accession number de la base de datos 'nucleotide' de NCBI
    # https://www.ncbi.nlm.nih.gov/nuccore.
    # Genera el archivo .gbk (GenBank) en la carpeta 'files' y retorna la secuencia en formato SeqRecord.
    path = os.path.join("files", accesion_number+".gbk")
    if not os.path.isfile(path):
        print("Descargando secuencia...")
        with Entrez.efetch(db="nucleotide", id=accesion_number, rettype="gb", retmode="text") as handle:
            with open(path, "w") as out_handle:
                out_handle.write(handle.read())
    record = SeqIO.read(path, "genbank")
    return record

def parse_file_to_seqrecord(filepath, type):
    # Función que recibe el handler de un archivo, su nombre, su tipo y lo retorna como SeqRecord.
    # Genera una copia del achivo en la carpeta files
    # 'type' puede ser "genbank", "fasta", "swiss" y otros soportados por SeqIO.
    currentpath = os.path.abspath(filepath)
    newpath = os.path.join(os.getcwd(), "files", os.path.split(currentpath)[1])
    if not os.path.isfile(newpath):
        shutil.copy(currentpath, newpath)

    record = SeqIO.read(newpath, type)
    return record

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-an',
                        '--accessionnumber',
                        help='Accesion number de la secuencia a descargar.',
                        required=True,
                        default=argparse.SUPPRESS)
    parser.add_argument('-db', '--database',help='Base de datos de la que se quiere descargar (nucleotide para secuencias).', type=str, default='nucleotide')

    args = parser.parse_args()
    dargs = vars(args)
    Entrez.mail = "simonvergaraswett@hotmail.com"

    if (dargs["database"] == 'nucleotide'):
        record = accnum_to_seqrecord(accesion_number=dargs["accessionnumber"])
        print(len(record))

if __name__ == "__main__":
    main()