seq = 'AGAGAGTGACGTA'
print(seq[0:2]*3)
checked = []
for i in range(len(seq)-2):
    subseq = str(seq[i:i+2])
    if subseq not in checked:
        if subseq * 3 in seq:
            print('Homopolimero de 2 nucleotidos')
        checked.append(subseq)