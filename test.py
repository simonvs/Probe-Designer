import descarga
#NG_008617.1
#AC005600.1
#NM_000186.4
record = descarga.accnum_to_seqrecord("NM_000186.4")
for f in record.features:
    print(f.type, len(f.location))