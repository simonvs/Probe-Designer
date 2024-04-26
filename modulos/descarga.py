import os
import argparse
import shutil
import gzip
#import filecmp
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
    #Si la carpeta de secuencias no existe, la crea
    folder_path = os.path.join(os.getcwd(), "files")
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)

    filepath = os.path.join(folder_path, accesion_number+".gbk")
    print(filepath)
    #Verifica si el archivo ya está descargado
    if not os.path.isfile(filepath):
        try:
            #Descargar la secuencia y escribirla en un archivo
            #Entrez.mail = "simonvergaraswett@hotmail.com"
            with Entrez.efetch(db="nucleotide", id=accesion_number, rettype="gb", retmode="text") as handle:
                with open(filepath, "w") as out_handle:
                    out_handle.write(handle.read())
        except:
            print("Error al descargar el archivo.")
    
    #Lee el archivo descargado y lo parsea a SeqRecord
    record = SeqIO.read(filepath, "genbank")
    return record

def parse_file_to_seqrecord(filepath):
    """
    Función que recibe la ruta de un archivo genético anotado y lo retorna como SeqRecord. Genera una copia
    del achivo en la carpeta files. Solo se aceptan archivos en formato "GenBank", o sea '.gbk' y
    '.gb'. Estos archivos pueden estar comprimidos en formato '.gz'.

    :param filepath: La ruta donde se encuentra el archivo que se quiere utilizar
    :return: La secuencia en formato SeqRecord con todas sus anotaciones.
    """
    currentpath = os.path.abspath(filepath)
    folder_path = os.path.join(os.getcwd(), "files")

    #Si la carpeta de secuencias no existe, la crea
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)

    #Se genera una copia del archivo selecionado en la carpeta de secuencias (si es que no existe)
    newpath = os.path.join(folder_path, os.path.split(currentpath)[1])
    if not os.path.isfile(newpath):# or not filecmp.cmp(currentpath, newpath):
        shutil.copy(currentpath, newpath)

    try:
        #Verificar si está comprimido
        if os.path.splitext(newpath)[1] == '.gz':
            with gzip.open(newpath, 'rb') as f_in:
                with open(os.path.splitext(newpath)[0], 'wb') as f_out:
                    f_out.write(f_in.read())

        #Parseo para GenBank
        if os.path.splitext(newpath)[1] == '.gbk' or os.path.splitext(newpath)[1] == '.gb':
            record = SeqIO.read(newpath, 'genbank')

        #Parseo para GFF3
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
    """
    Ejecutar directamente el archivo 'descarga' requiere entregar el parámetro --accessionnumber para
    descargar una secuencia mediante la función accnum_to_seqrecord.

    :param --accesionnumber: El accesion number del archivo en GenBank, por ejemplo NG_008617.1
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-an',
                        '--accessionnumber',
                        help='Accesion number de la secuencia a descargar.',
                        required=True,
                        default=argparse.SUPPRESS)

    args = parser.parse_args()
    dargs = vars(args)

    record = accnum_to_seqrecord(accesion_number=dargs["accessionnumber"])
    print(len(record))

if __name__ == "__main__":
    main()