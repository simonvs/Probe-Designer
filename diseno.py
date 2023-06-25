import descarga
import argparse
import pandas as pd
import primer3
from Bio.SeqUtils import GC, MeltingTemp
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows


def check_probe(seq, query, minlen, maxlen, tmmin, tmmax, gcmin, gcmax):
    """
    Función que recibe una sonda (secuencia) y las restricciones para verificar que la
    sonda cumpla los parámetros establecidos.
    Retorna un valor booleano que determina si cumple las restricciones o no.
    :param seq: La secuencia que se quiere verificar.
    :param record:
    :param minlen: El largo mínimo de la sonda.
    :param maxlen: El largo máximo de la sonda.
    :return: La secuencia en formato SeqRecord con todas sus anotaciones.
    """

    #Chequear largo
    if len(seq) < minlen or len(seq) > maxlen:
        #print('Largo no valido')
        return False
    
    #Chequear Tm
    if MeltingTemp.Tm_NN(seq) < tmmin or MeltingTemp.Tm_NN(seq) > tmmax:
        #print('Tm fuera de rango')
        return False
    
    #Chequear %GC
    if GC(seq) < gcmin or GC(seq) > gcmax:
        #print('GC fuera de rango')
        return False
    
    #Chequear homopolímeros
    #dg_hairpin = primer3.calcHairpin(str(seq)).dg
    #dg_homodimero = primer3.calcHomodimer(str(seq)).dg

    #Chequear especificidad usando alineamiento
    #aligner = PairwiseAligner()
    #aligner.mode = 'local'
    #aligner.match_score = 1  # Puntaje por coincidencia
    #aligner.mismatch_score = -1  # Puntaje por desajuste
    #aligner.open_gap_score = -1  # Puntaje por abrir un gap
    #aligner.extend_gap_score = -0.5  # Puntaje por extender un gap

    #alineamiento = aligner.align(seq.complement(), str(query))

    #PENDIENTE: ¿Qué criterio uso para considerarlo específico?
    #if(sorted(alineamiento)[0].score >= float(len(seq)/4)):
    #    return False
    return True


def get_probe_from_pos(empalme, record, site, minlen, maxlen, tmmin, tmmax, gcmin, gcmax):
    i = minlen
    b = False
    while i <= maxlen and not b:
        if site == -1:
            probe = record.seq[empalme[0]-i:empalme[0]].reverse_complement()
            query = record.seq[:empalme[0]-i-1] + record.seq[empalme[0]+1:]
        elif site == 1:
            probe = record.seq[empalme[1]:empalme[1]+i].reverse_complement()
            query = record.seq[:empalme[1]-1] + record.seq[empalme[1]+i+1:]
        else:
            probe = (record.seq[empalme[0]-i//2:empalme[0]] + record.seq[empalme[1]:empalme[1]+i//2]).reverse_complement()
            query = record.seq
        b = check_probe(probe, query, minlen, maxlen, tmmin, tmmax, gcmin, gcmax)
        i += 1
    if not b:
        return Seq('AAA')
    return probe


def get_splicing_pairs(locations):  
    positions = []
    for subloc in locations.parts:
        positions.append((subloc.start, subloc.end))

    empalmes = []
    for i in range(len(positions)-1):
        empalmes.append((positions[i][1].position, positions[i+1][0].position))
    return empalmes


def probe_designer(record, minlen=60, maxlen=120, tmmin=65, tmmax=75, gcmin=30, gcmax=70):

    locations = []
    dict_cds = {
        'gen':[],
        'info':[],
        'empalmes':[]
    }

    for f in record.features:
        if f.type == 'CDS':
            loc = f.location
            if str(loc) not in locations:
                locations.append(str(loc))
                dict_cds['gen'].append(f.qualifiers.get('gene', ['Gen no identificado'])[0])
                dict_cds['info'].append(f.qualifiers.get('product', ['Transctipto no identificado'])[0])
                dict_cds['empalmes'].append(get_splicing_pairs(loc))

    df_cds = pd.DataFrame(dict_cds)

    dict_probes = {
        'gen':[],
        'transcrito':[],
        'empalme':[],
        'sonda':[],
        'lado':[],
        'tm':[],
        'gc':[],
        'largo':[]
    }

    for index, row in df_cds.iterrows():
        for empalme in row['empalmes']:
            donor = get_probe_from_pos(empalme, record, -1, minlen, maxlen, tmmin, tmmax, gcmin, gcmax)
            dict_probes['gen'].append(row['gen'])
            dict_probes['transcrito'].append(row['info'])
            dict_probes['empalme'].append(str(empalme))
            dict_probes['sonda'].append(str(donor))
            dict_probes['lado'].append('donor site')
            dict_probes['tm'].append(int(MeltingTemp.Tm_NN(donor)))
            dict_probes['gc'].append(int(GC(donor)))
            dict_probes['largo'].append(len(donor))

            acceptor = get_probe_from_pos(empalme, record, 1, minlen, maxlen, tmmin, tmmax, gcmin, gcmax)
            dict_probes['gen'].append(row['gen'])
            dict_probes['transcrito'].append(row['info'])
            dict_probes['empalme'].append(str(empalme))
            dict_probes['sonda'].append(str(acceptor))
            dict_probes['lado'].append('acceptor site')
            dict_probes['tm'].append(int(MeltingTemp.Tm_NN(acceptor)))
            dict_probes['gc'].append(int(GC(acceptor)))
            dict_probes['largo'].append(len(acceptor))

            central = get_probe_from_pos(empalme, record, 0, minlen, maxlen, tmmin, tmmax, gcmin, gcmax)
            dict_probes['gen'].append(row['gen'])
            dict_probes['transcrito'].append(row['info'])
            dict_probes['empalme'].append(str(empalme))
            dict_probes['sonda'].append(str(central))
            dict_probes['lado'].append('central')
            dict_probes['tm'].append(int(MeltingTemp.Tm_NN(central)))
            dict_probes['gc'].append(int(GC(central)))
            dict_probes['largo'].append(len(central))

    df_probes = pd.DataFrame(dict_probes)
    return df_probes

def generate_xlsx(df, minlen, maxlen, tmmin, tmmax, gcmin, gcmax):
    wb = Workbook()

    sheet1 = wb.active
    sheet1.title = 'Parámetros'

    sheet1.cell(row=1, column=1).value = 'PARÁMETROS DE DISEÑO'

    sheet1.cell(row=3, column=1).value = 'Largo de sonda mínimo'
    sheet1.cell(row=3, column=2).value = minlen

    sheet1.cell(row=4, column=1).value = 'Largo de sonda máximo'
    sheet1.cell(row=4, column=2).value = maxlen

    sheet1.cell(row=6, column=1).value = 'Temp Melting mínima'
    sheet1.cell(row=6, column=2).value = tmmin

    sheet1.cell(row=7, column=1).value = 'Temp Melting máxima'
    sheet1.cell(row=7, column=2).value = tmmax

    sheet1.cell(row=9, column=1).value = '%GC mínimo'
    sheet1.cell(row=9, column=2).value = gcmin

    sheet1.cell(row=10, column=1).value = '%GC máximo'
    sheet1.cell(row=10, column=2).value = gcmax

    for column_cells in sheet1.columns:
        max_length = 2
        column = column_cells[0].column_letter
        for cell in column_cells:
            try:
                if len(str(cell.value)) > max_length:
                    max_length = len(cell.value)
            except:
                pass
        adjusted_width = (max_length + 2) * 1.2
        sheet1.column_dimensions[column].width = adjusted_width



    sheet2 = wb.create_sheet("Sondas")
    for r in dataframe_to_rows(df, index=True, header=True):
        sheet2.append(r)

    for cell in sheet2['A'] + sheet2[1]:
        cell.style = 'Pandas'

    #Ajustar ancho de columnas
    for column_cells in sheet2.columns:
        max_length = 2
        column = column_cells[0].column_letter
        for cell in column_cells:
            try:
                if len(str(cell.value)) > max_length:
                    max_length = len(cell.value)
            except:
                pass
        adjusted_width = (max_length + 2) * 1.2
        sheet2.column_dimensions[column].width = adjusted_width

    

    wb.save('archivo_excel.xlsx') # nombre: datetime?
    return 0

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-an',
                        '--accessionnumber',
                        help='Accesion number de la secuencia a analizar.',
                        required=True,
                        default=argparse.SUPPRESS)
    parser.add_argument('-minlen', '--minimumlength',help='Largo mínimo de sonda.' ,type=int,default=60)
    parser.add_argument('-maxlen', '--maximumlength',help='Largo mínimo de sonda.' ,type=int,default=120)
    parser.add_argument('-tmmin', '--mintempmelting',help='Temperatura de melting mínima.' ,type=int,default=65)
    parser.add_argument('-tmmax', '--maxtempmelting',help='Temperatura de melting máxima.' ,type=int,default=75)
    parser.add_argument('-gcmin', '--minumumgc',help='Porcentaje de GC mínimo.' ,type=int,default=30)
    parser.add_argument('-gcmax', '--maximumgc',help='Porcentaje de GC máximo.' ,type=int,default=70)

    args = parser.parse_args()
    dargs = vars(args)
    rec = descarga.parse_file_to_seqrecord('C:/Users/simon/Downloads/tp53.gb')
    #rec = descarga.accnum_to_seqrecord(dargs["accessionnumber"])
    df = probe_designer(record=rec,
                   minlen=dargs["minimumlength"],
                   maxlen=dargs["maximumlength"],
                   tmmin=dargs["mintempmelting"],
                   tmmax=dargs["maxtempmelting"],
                   gcmin=dargs["minumumgc"],
                   gcmax=dargs["maximumgc"])
    
    generate_xlsx(df=df,
                minlen=dargs["minimumlength"],
                maxlen=dargs["maximumlength"],
                tmmin=dargs["mintempmelting"],
                tmmax=dargs["maxtempmelting"],
                gcmin=dargs["minumumgc"],
                gcmax=dargs["maximumgc"])

if __name__ == "__main__":
    main()