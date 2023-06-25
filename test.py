# from dna_features_viewer import BiopythonTranslator

# graphic_record = BiopythonTranslator().translate_record("C:/Users/simon/Downloads/tp53.gb")
# ax, _ = graphic_record.plot(figure_width=10, strand_in_label_threshold=7)
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

import re
print([match.group() for match in re.finditer(r"([A-Z]+)\1{1,}", "ACTGTAGGGCTGTACTACTACTATG")])