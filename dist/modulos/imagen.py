import matplotlib
import matplotlib.pyplot as plt
import descarga
import pandas as pd
import os
matplotlib.use('agg')

def get_all_transcripts(seqrecord):
    transcripts = []
    for f in seqrecord.features:
        if f.type == 'CDS':
            info = f.qualifiers.get('product')[0]
            if str(info) not in transcripts:
                transcripts.append(str(info))
    return transcripts

def get_splicings(seqrecord, transcripts):
    empalmes_array = []
    for f in seqrecord.features:
        if f.type == 'CDS':
            if f.qualifiers.get('product')[0] in transcripts:
                loc = f.location
                pares_empalme = get_splicing_pairs(loc)
                for par in pares_empalme:
                    if par not in empalmes_array:
                        empalmes_array.append(par)
    return empalmes_array

def get_exon_positions(locations):
   
    positions = []
    for subloc in locations.parts:
        positions.append((subloc.start, subloc.end))

    return positions

def get_splicing_pairs(locations):
    positions = []
    for subloc in locations.parts:
        positions.append((subloc.start, subloc.end))

    empalmes = []
    for i in range(len(positions)-1):
        #empalmes.append((positions[i][1].position, positions[i+1][0].position))
        empalmes.append((int(positions[i][1]), int(positions[i+1][0])))


    return empalmes

def plot_isoforms(record, transcripts, empalmes_array, foldername):

    locations = []
    dict_cds = {
        'info':[],
        'exones':[],
        'empalmes':[]
    }

    for f in record.features:
        if f.type == 'CDS':
            loc = f.location
            info = f.qualifiers.get('product')[0]
            if str(loc) not in locations and info in transcripts:
                locations.append(str(loc))
                dict_cds['info'].append(info)
                dict_cds['exones'].append(get_exon_positions(loc))
                dict_cds['empalmes'].append(get_splicing_pairs(loc))
                

    df_cds = pd.DataFrame(dict_cds)

    fig, ax = plt.subplots(figsize=(23, len(transcripts) * 0.5+5))

    yticks = []
    yticklabels = []


    xmin = float('inf')
    xmax = float('-inf')

    for idx, row in df_cds.iterrows():
        yticks.append(idx)
        yticklabels.append(row['info'])

        exones = row['exones']
        empalmes = row['empalmes']
        for start, end in exones:
            ax.add_patch(plt.Rectangle((start, idx-0.15), end - start, 0.2, color='gray'))
            xmin = min(xmin, start)
            xmax = max(xmax, end)

        for line_start, line_end in empalmes:
            id_empalme = str(empalmes_array.index((line_start, line_end)) + 1)
            mitad = (line_end+line_start)/2
            ax.plot([line_start, mitad], [idx, idx+0.1], color='red')
            ax.plot([mitad, line_end], [idx+0.1, idx], color='red')
            ax.text(mitad, idx + 0.15, id_empalme, fontsize=9 ,ha='center', va='bottom', color='red')
        
        

    ax.set_xlim(xmin-500, xmax+500)
    ax.set_ylim(-0.6, len(transcripts)-0.5)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    ax.set_xlabel('Posición')
    ax.set_title('Transcripciones de '+record.id)

    plt.savefig(os.path.join('sondas',foldername,record.id+'.png'))
    #plt.savefig(os.path.join('sondas',record.id+'.png'))
    #plt.show()

def plot_show():
    plt.show()

if __name__ == "__main__":
    sequence_data = descarga.parse_file_to_seqrecord('files/NG_008604.gbk')
    plot_isoforms(sequence_data, get_all_transcripts(sequence_data), get_splicings(sequence_data, get_all_transcripts(sequence_data)),"images")



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