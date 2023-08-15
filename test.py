# from dna_features_viewer import BiopythonTranslator

# graphic_record = BiopythonTranslator().translate_record("C:/Users/simon/Downloads/tp53.gb")
# ax, _ = graphic_record.plot(figure_width=15, strand_in_label_threshold=7)
# ax.figure.tight_layout()
# ax.figure.savefig('sequence_and_translation.png')
#
#####################################################
#
# from Bio.Align import PairwiseAligner

# # Crear un objeto PairwiseAligner
# aligner = PairwiseAligner()

# # Definir las dos secuencias que deseas alinear
# secuencia1 = "ACGT"
# secuencia2 = "GGGGGGACGGTGAAAATGTAACGGGT"

# # Configurar los parámetros del alineador
# aligner.mode = 'local'  # Alineamiento global
# aligner.match_score = 1  # Puntaje por coincidencia
# aligner.mismatch_score = -1  # Puntaje por desajuste
# aligner.open_gap_score = -1  # Puntaje por abrir un gap
# aligner.extend_gap_score = -0.5  # Puntaje por extender un gap

# # Realizar el alineamiento entre las dos secuencias
# alineamiento = aligner.align(secuencia1, secuencia2)

# print(sorted(alineamiento)[0])
# #for alignment in sorted(alineamiento):
# #    print("Score = %.3f:" % alignment.score)
# #    print(alignment)

#import re
#print([match.group() for match in re.finditer(r"([A-Z]+)\1{1,}", "ACTGTAGGGCTGTACTACTACTATG")])
import subprocess
from Bio import SearchIO


def verify_specificity(seq, iscentral):
    """
    Función que verifica la especificidad o unicidad de una secuencia usando la herramienta de alineamiento BLAST de manera local.
    Es necesario tener instalado Blast+ y las bases de datos en la carpeta 'refseq'.
    :param seq: Secuencia de nuecleótidos que se quiere validar.
    :param iscentral: Booleano que indica si se trata de una sonda central. Se utilizan diferentes criterios para chequear la especificidad de una sonda central.
    :return: Bolleano que indica si la secuencia se considera específica.
    """
    db1_name = "refseq/db1_genoma"
    db2_name = "refseq/db2_transcriptoma"
    query = "refseq/query.txt"
    outxml = 'refseq/resultados.xml'

    with open(query, "w") as archivo:
        archivo.write(">1\n")
        archivo.write(str(seq))

    blastn_cmd = ['blastn',
                '-task', 'blastn-short',
                '-query', query,
                '-db', db1_name,
                '-out', outxml,
                '-outfmt', '5',
                '-word_size', '11',
                '-dust', 'no']

    subprocess.run(blastn_cmd, shell=True)
    blast_qresult = SearchIO.read(outxml, "blast-xml")

    print(blast_qresult[0].hsps[0].evalue)
    print(blast_qresult[0].hsps[1].evalue)
    print(blast_qresult[0].hsps[2].evalue)
    print(blast_qresult[1].hsps[0].evalue)
    print(blast_qresult[1].hsps[1].evalue)
    print(blast_qresult[1].hsps[2].evalue)

    if iscentral:
        if blast_qresult[0].hsps[2].evalue / blast_qresult[0].hsps[1].evalue < 100:
            print('1')
            return False
        if blast_qresult[1].hsps[0].evalue / blast_qresult[0].hsps[1].evalue < 100:
            print('2')
            return False
        
    else:
        if blast_qresult[0].hsps[1].evalue / blast_qresult[0].hsps[0].evalue < 100:
            print('3')
            return False
        if blast_qresult[1].hsps[0].evalue / blast_qresult[0].hsps[0].evalue < 100:
            print('4')
            return False
    return True

print(verify_specificity('TGTCCTTACCAGAACGTTGTTTTCAGGAAGTAGTTTCCATAGGTCTGAAAATGTTTCCTGAC', False))
# def validar_repeticiones(seq):
#     secuencia_adn = str(seq)
#     n = len(secuencia_adn)

#     for i in range(n - 3):
#         subsecuencia = secuencia_adn[i:i+4]
        
#         if all(subseq):

#         if all(subsecuencia == secuencia_adn[j:j+4] for j in range(i+1, i+4)):
#             return True

#     return False

# # Ejemplo de uso
# secuencia_adn = Seq("ACGTAGTAGTAGTAGTTTACCC")
# print(secuencia_adn)
# resultado = validar_repeticiones(secuencia_adn)

# if resultado:
#     print("La secuencia contiene una subsecuencia de 1, 2 o 3 nucleótidos que se repite 4 veces seguidas.")
# else:
#     print("La secuencia no contiene una subsecuencia de 1, 2 o 3 nucleótidos que se repite 4 veces seguidas.")



#import primer3

#print(primer3.calc_hairpin("AGTTTCCATAGGTCTGAAAATGTTTCCTGCGTCGAGGGGGGCTCGACGCTAGGATCTGAC").dg)