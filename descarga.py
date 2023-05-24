import os
import argparse
import shutil
from Bio import Entrez, SeqIO
from BCBio import GFF

def accnum_to_seqrecord(accesion_number):
    """
    Función que recibe un accession number de la base de datos 'nucleotide' de NCBI
    https://www.ncbi.nlm.nih.gov/nuccore.
    Genera el archivo .gbk (GenBank) en la carpeta 'files' y retorna la secuencia en formato SeqRecord.

    :param accesion_number: El accesion number del archivo en GenBank, por ejemplo NG_008617.1
    :return: La secuencia en formato SeqRecord con todas sus anotaciones.
    """
    folder_path = os.path.join(os.getcwd(), "files")
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)

    filepath = os.path.join(folder_path, accesion_number+".gbk")
    if not os.path.isfile(filepath):
        try:
            print("Descargando secuencia...")
            with Entrez.efetch(db="nucleotide", id=accesion_number, rettype="gb", retmode="text") as handle:
                with open(filepath, "w") as out_handle:
                    out_handle.write(handle.read())
        except:
            print("Error al descargar el archivo.")
    record = SeqIO.read(filepath, "genbank")
    return record

def parse_file_to_seqrecord(filepath):
    """
    Función que recibe la ruta de un archivo genético anotado y lo retorna como SeqRecord. Genera una copia del achivo en la carpeta files.
    Solo se aceptan archivos en formato "genbank" y "gff", o sea '.gbk', '.gb' y '.gff'.

    :param filepath: La ruta donde se encuentra el archivo que se quiere utilizar
    :return: La secuencia en formato SeqRecord con todas sus anotaciones.
    """
    currentpath = os.path.abspath(filepath)
    folder_path = os.path.join(os.getcwd(), "files")
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)

    newpath = os.path.join(folder_path, os.path.split(currentpath)[1])
    if not os.path.isfile(newpath):
        shutil.copy(currentpath, newpath)
    try:
        if os.path.splitext(newpath)[1] == '.gbk' or os.path.splitext(newpath)[1] == '.gb':
            record = SeqIO.read(newpath, 'genbank')
        elif os.path.splitext(newpath)[1] == '.gff':
            with open(newpath) as in_handle:
                records = GFF.parse(in_handle)
                record = records[0]
        else:
            raise TypeError
        return record
    except:
        print("Error al parsear el archivo.")


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-an',
                        '--accessionnumber',
                        help='Accesion number de la secuencia a descargar.',
                        required=True,
                        default=argparse.SUPPRESS)

    args = parser.parse_args()
    dargs = vars(args)
    Entrez.mail = "simonvergaraswett@hotmail.com"

    if (dargs["database"] == 'nucleotide'):
        record = accnum_to_seqrecord(accesion_number=dargs["accessionnumber"])
        print(len(record))

if __name__ == "__main__":
    main()