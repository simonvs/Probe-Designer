from dna_features_viewer import BiopythonTranslator

graphic_record = BiopythonTranslator().translate_record("C:/Users/simon/Downloads/tp53.gb")
ax, _ = graphic_record.plot(figure_width=10, strand_in_label_threshold=7)
ax.figure.tight_layout()
ax.figure.savefig('sequence_and_translation.png')