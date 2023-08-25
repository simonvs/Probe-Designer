import matplotlib.pyplot as plt
import descarga
import pandas as pd
# Ejemplo de datos de secuencia en formato SeqRecord de Biopython
# Aquí debes cargar tu propio archivo o crear tus propios registros SeqRecord

sequence_data = descarga.parse_file_to_seqrecord('C:\Users\simon\Documents\GitHub\Probe-Designer\files\NM_000186.4.gbk')

def get_exon_positions(locations):
   
    positions = []
    for subloc in locations.parts:
        positions.append((subloc.start, subloc.end))

    return positions

def plot_isoforms(record, name):

    locations = []
    dict_cds = {
        'info':[],
        'exones':[]
    }

    for f in record.features:
        if f.type == 'CDS':
            loc = f.location
            if str(loc) not in locations:
                locations.append(str(loc))
                dict_cds['info'].append(f.qualifiers.get('product', ['Transctipto no identificado'])[0])
                dict_cds['exones'].append(get_exon_positions(loc))

    df_cds = pd.DataFrame(dict_cds)

    fig, ax = plt.subplots(figsize=(10, len(dict_cds) * 1.5))

    yticks = []
    yticklabels = []

    xmin = float('inf')
    xmax = float('-inf')

    for idx, row in df_cds.iterrows():
        yticks.append(idx)
        yticklabels.append(row['info'])

        exons = row['exones']
        for start, end in exons:
            ax.add_patch(plt.Rectangle((start, idx), end - start, 0.3, color='gray'))
            xmin = min(xmin, start)
            xmax = max(xmax, end)
        
    

    ax.set_xlim(xmin-500, xmax+500)
    ax.set_ylim(-0.5, len(dict_cds) - 0.5)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    ax.set_xlabel('Posición')
    ax.set_title('Transcripciones de '+name)

    plt.show()

# Llama a la función con tus datos de secuencia
plot_isoforms(sequence_data, str(sequence_data.id))

#########
# from dna_features_viewer import BiopythonTranslator

# graphic_record = BiopythonTranslator().translate_record("C:/Users/simon/Downloads/tp53.gb")
# ax, _ = graphic_record.plot(figure_width=15, strand_in_label_threshold=7)
# ax.figure.tight_layout()
# ax.figure.savefig('sequence_and_translation.png')
#
#####################################################


# import subprocess
# from Bio import SearchIO


# def verify_specificity(seq, iscentral):
#     """
#     Función que verifica la especificidad o unicidad de una secuencia usando la herramienta de alineamiento BLAST de manera local.
#     Es necesario tener instalado Blast+ y las bases de datos en la carpeta 'refseq'.
#     :param seq: Secuencia de nuecleótidos que se quiere validar.
#     :param iscentral: Booleano que indica si se trata de una sonda central. Se utilizan diferentes criterios para chequear la especificidad de una sonda central.
#     :return: Bolleano que indica si la secuencia se considera específica.
#     """
#     db1_name = "refseq/db1_genoma"
#     db2_name = "refseq/db2_transcriptoma"
#     query = "refseq/query.txt"
#     outxml = 'refseq/resultados.xml'

#     with open(query, "w") as archivo:
#         archivo.write(">1\n")
#         archivo.write(str(seq))

#     blastn_cmd = ['blastn',
#                 '-task', 'blastn-short',
#                 '-query', query,
#                 '-db', db1_name,
#                 '-out', outxml,
#                 '-outfmt', '5',
#                 '-word_size', '11',
#                 '-dust', 'no']

#     subprocess.run(blastn_cmd, shell=True)
#     blast_qresult = SearchIO.read(outxml, "blast-xml")

#     print(blast_qresult[0].hsps[0].evalue)
#     print(blast_qresult[0].hsps[1].evalue)
#     print(blast_qresult[0].hsps[2].evalue)
#     print(blast_qresult[1].hsps[0].evalue)
#     print(blast_qresult[1].hsps[1].evalue)
#     print(blast_qresult[1].hsps[2].evalue)

#     if iscentral:
#         if blast_qresult[0].hsps[2].evalue / blast_qresult[0].hsps[1].evalue < 100:
#             print('1')
#             return False
#         if blast_qresult[1].hsps[0].evalue / blast_qresult[0].hsps[1].evalue < 100:
#             print('2')
#             return False
        
#     else:
#         if blast_qresult[0].hsps[1].evalue / blast_qresult[0].hsps[0].evalue < 100:
#             print('3')
#             return False
#         if blast_qresult[1].hsps[0].evalue / blast_qresult[0].hsps[0].evalue < 100:
#             print('4')
#             return False
#     return True

# print(verify_specificity('TGTCCTTACCAGAACGTTGTTTTCAGGAAGTAGTTTCCATAGGTCTGAAAATGTTTCCTGAC', False))


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