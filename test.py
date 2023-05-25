import descarga
#NG_008617.1
#AC005600.1
#NM_000186.4 
#NC_000016.1 chr16
#NM_000297.4 pkd2

record = descarga.parse_file_to_seqrecord("C:/Users/simon/Documents/NG_008617.1.gbk")
exones = 0
for f in record.features:
    if f.type == 'exon':
        exones += 1
        print("Inicio del exon {}: {}".format(exones, record.seq[f.location.start:f.location.start+20]))


record2 = descarga.accnum_to_seqrecord('NM_001077690')
features = 0
exones2 = 0
for f in record2.features:
    if f.type == 'exon':
        exones2 += 1
        print("Inicio del exon {}: {}".format(exones2, record.seq[f.location.start:f.location.start+20]))
