import networkx as nx
from primer3 import calc_heterodimer
import queue
from Bio.SeqUtils import MeltingTemp

def are_sequences_compatible(seq1, seq2, mindg, maxdt):
    """
    Esta función determina si dos secuencias de ADN son compatibles. 
    Debes personalizar las condiciones de compatibilidad según tus necesidades.

    :param seq1: Primera secuencia de ADN.
    :param seq2: Segunda secuencia de ADN.
    :return: True si son compatibles, False en caso contrario.
    """
    
    
    temp1 = MeltingTemp.Tm_NN(seq1)
    temp2 = MeltingTemp.Tm_NN(seq2)
    if abs(temp2-temp1) > maxdt:
        #print("No compatibles, diferencia tm: "+str(abs(temp2-temp1)))
        return False
    
    if len(seq1)<=60 or len(seq2)<=60:
        if calc_heterodimer(seq1,seq2).dg < mindg:
            #print('dg')
            return False
    else:
        if calc_heterodimer(seq1[:60],seq2[:60]).dg < mindg:
            #print('dg')
            return False
        if calc_heterodimer(seq1[:60],seq2[60:]).dg < mindg:
            #print('dg')
            return False
        if calc_heterodimer(seq1[60:],seq2[:60]).dg < mindg:
            #print('dg')
            return False
        if calc_heterodimer(seq1[60:],seq2[60:]).dg < mindg:
            #print('dg')
            return False
        
    #print("Compatibles!")
    return True
    

def multiplex_sequences(sequences, progress_queue, mindg=-13627, maxdt=5):
    """
    Función principal de multiplexación de sondas. Recibe las secuencias, los parámetros de multiplexación y la cola para informar progreso.

    :param sequences: Lista con todas las diferentes secuencias.
    :param progress_queue: Cola para informar el progreso a la interfaz.
    :param mindg: Mínimo delta G heterodimerización para que dos secuencias sean compatibles.
    :param maxdt: Máxima diferencia de Tm para que dos secuencias sean compatibles.
    :return: True si son compatibles, False en caso contrario.
    """
    # Crea un grafo no dirigido.
    G = nx.Graph()

    # Agrega nodos (secuencias) al grafo.
    G.add_nodes_from(sequences)

    progreso_actual = 100
    #print("num seqs:",len(sequences))
    # Comprueba la compatibilidad y agrega aristas entre secuencias no compatibles.
    for i, seq1 in enumerate(sequences):
        for j, seq2 in enumerate(sequences):
            if i < j and not are_sequences_compatible(seq1, seq2, mindg, maxdt):
                G.add_edge(seq1, seq2)        
        progreso_actual += 95 / (len(sequences))
        if (progress_queue):
            progress_queue.put(progreso_actual)
        #print(progreso_actual)


    # Colorea el grafo para encontrar grupos.
    color_map = nx.coloring.greedy_color(G, strategy="largest_first")
    if (progress_queue):
            progress_queue.put(200)

    # Crea un diccionario que asocie cada secuencia con su grupo.
    sequence_to_group = {sequence: color + 1 for sequence, color in color_map.items()}

    #for seq, group in sequence_to_group.items():
    #    print(f"Secuencia: {seq}, Grupo: {group}, tm: {MeltingTemp.Tm_NN(seq)}")
    return sequence_to_group

# Ejemplo de uso:
if __name__ == "__main__":
    sequences = [
        "ATCGATCGATCG",
        "ATCGAAGCTAGC",
        "TACGATTAGCTA",
        "ATCGTAGCTAGC",
        "GCTAGCTAGCTA",
        "TACGCTAGCTAG"
    ]

    result = multiplex_sequences(sequences)

    for seq, group in result.items():
        print(f"Secuencia: {seq}, Grupo: {group}, tm: {MeltingTemp.Tm_NN(seq)}")
