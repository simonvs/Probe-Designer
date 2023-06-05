import descarga

import argparse
import pandas as pd
from Bio.SeqUtils import GC, MeltingTemp
from Bio.Seq import Seq


def check_probe(seq, minlen, maxlen, tmmin, tmmax, gcmin, gcmax):
    if len(seq) < minlen or len(seq) > maxlen:
        #print('Largo no valido')
        return False
    if MeltingTemp.Tm_NN(seq) < tmmin or MeltingTemp.Tm_NN(seq) > tmmax:
        #print('Tm fuera de rango')
        return False
    if GC(seq) < gcmin or GC(seq) > gcmax:
        #print('GC fuera de rango')
        return False
    return True


def get_probe_from_pos(empalme, record, site, minlen, maxlen, tmmin, tmmax, gcmin, gcmax):
    i = minlen
    b = False
    while i <= maxlen and not b:
        if site == -1:
            probe = record.seq[empalme[0]-i:empalme[0]].reverse_complement()
        elif site == 1:
            probe = record.seq[empalme[1]:empalme[1]+i].reverse_complement()
        else:
            probe = record.seq[empalme[0]-i//2:empalme[0]] + record.seq[empalme[1]:empalme[1]+i//2].reverse_complement()
        b = check_probe(probe, minlen, maxlen, tmmin, tmmax, gcmin, gcmax)
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


def probe_designer(record, minlen=10, maxlen=50, tmmin=50, tmmax=80, gcmin=30, gcmax=70):

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
            dict_probes['empalme'].append(empalme)
            dict_probes['sonda'].append(str(donor))
            dict_probes['lado'].append('donor site')
            dict_probes['tm'].append(MeltingTemp.Tm_NN(donor))
            dict_probes['gc'].append(GC(donor))
            dict_probes['largo'].append(len(donor))

            acceptor = get_probe_from_pos(empalme, record, 1, minlen, maxlen, tmmin, tmmax, gcmin, gcmax)
            dict_probes['gen'].append(row['gen'])
            dict_probes['transcrito'].append(row['info'])
            dict_probes['empalme'].append(empalme)
            dict_probes['sonda'].append(str(acceptor))
            dict_probes['lado'].append('acceptor site')
            dict_probes['tm'].append(MeltingTemp.Tm_NN(acceptor))
            dict_probes['gc'].append(GC(acceptor))
            dict_probes['largo'].append(len(acceptor))

            central = get_probe_from_pos(empalme, record, 0, minlen, maxlen, tmmin, tmmax, gcmin, gcmax)
            dict_probes['gen'].append(row['gen'])
            dict_probes['transcrito'].append(row['info'])
            dict_probes['empalme'].append(empalme)
            dict_probes['sonda'].append(str(central))
            dict_probes['lado'].append('central')
            dict_probes['tm'].append(MeltingTemp.Tm_NN(central))
            dict_probes['gc'].append(GC(central))
            dict_probes['largo'].append(len(central))

    df_probes = pd.DataFrame(dict_probes)
    return df_probes


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-an',
                        '--accessionnumber',
                        help='Accesion number de la secuencia a analizar.',
                        required=True,
                        default=argparse.SUPPRESS)
    parser.add_argument('-minlen', '--minimumlength',help='Largo mínimo de sonda.' ,type=int,default=10)
    parser.add_argument('-maxlen', '--maximumlength',help='Largo mínimo de sonda.' ,type=int,default=50)
    parser.add_argument('-tmmin', '--mintempmelting',help='Temperatura de melting mínima.' ,type=int,default=40)
    parser.add_argument('-tmmax', '--maxtempmelting',help='Temperatura de melting máxima.' ,type=int,default=70)
    parser.add_argument('-gcmin', '--minumumgc',help='Porcentaje de GC mínimo.' ,type=int,default=35)
    parser.add_argument('-gcmax', '--maximumgc',help='Porcentaje de GC máximo.' ,type=int,default=65)

    args = parser.parse_args()
    dargs = vars(args)
    #rec = descarga.parse_file_to_seqrecord('C:/Users/simon/Downloads/tp53.gb')
    rec = descarga.accnum_to_seqrecord(dargs["accessionnumber"])
    probe_designer(record=rec,
                   minlen=dargs["minimumlength"],
                   maxlen=dargs["maximumlength"],
                   tmmin=dargs["mintempmelting"],
                   tmmax=dargs["maxtempmelting"],
                   gcmin=dargs["minumumgc"],
                   gcmax=dargs["maximumgc"]).to_csv('files/probes.csv')

if __name__ == "__main__":
    main()