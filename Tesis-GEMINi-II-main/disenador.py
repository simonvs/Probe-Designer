import docx
import os
from os import path
import pdb
import math
import primer3
import json
from datetime import datetime
from PIL import Image, ImageDraw, ImageFont

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

import time

exon_position = []
full_gene = ""
primer_list = []

PKD1_exon_desde = 0
contador_amplicones_PKD1 = 0
#"PKD1 exon 2-21.docx"         
#"PKD2.docx"
#"Secuencia CFH hg38 RefSeq NM000186.4.docx"
#"CFI NM000204.4.docx"
#"TTN.docx"
#PKD1 exon 22-34
#"ALG9.docx"
#"LR Ex27-34.docx"
#"GANAB.docx"

var = {
    'gene_filename' : "PKD2.docx",
    'exones_seleccionados' : [],
    'margin' : 15,
    'margin_between_exons' : 100,
    'primer_min_length' : 21,
    'primer_max_length' : 25,
    'max_gc_percent' : 0.6,
    'min_gc_percent' : 0.4,
    'mv_conc' : 50,
    'tm_min' : 59,
    'delta_tm_min' : 0,
    'delta_tm_max' : 4,
    'largo_amplicon': 400,
    'traslape' : 25,
    'max_homopol' : 4,
    'margen_busqueda' : 100,
    'deltag' : -6000,
    'intervalo' : 50,
    'aplicar_intervalo' : True,
    'overhangF' : "",
    'overhangR' : ""
}

def calcuTm(seq):
    met=1
    if met==0:
        tm=seq.count('a')*2+seq.count('A')*2+seq.count('t')*2+seq.count('T')*2+seq.count('c')*4+seq.count('C')*4+seq.count('g')*4+seq.count('G')*4
        return tm
    if met==1:
        tm=primer3.calcTm(seq,mv_conc=var["mv_conc"])
        return round(tm)

def areas_of_interest_set(margin,margin_between_exons):
    """!
    En base a las posiciones de cada exón en el gen (secuencias de caracteres en azul) genera
    las zonas de interés, a las que se les incluye el margen (margin) en los extremos de estas
    y si los exones estan a una distancia minima entre ellos (margin_between_exons) se incluirán
    dentro de una misma zona.
    \param <areas_of_interest_pos> {lista, en esta se almacena una lista con el inicio y fin de cada zona}
    """    
    global exon_position
    areas_of_interest_pos = []
    areas_of_interest_pos.append([((exon_position[0][0])-(margin)),exon_position[0][1],True])
    i = 0
    j = 0

    while(i< len(exon_position)):

        if (i == (len(exon_position)-1) ):
            areas_of_interest_pos[j][1] = exon_position[i][1]+margin
            i+=1
        else:
            if( (exon_position[i+1][0]-exon_position[i][1]) < margin_between_exons):
                i+=1
            else:
                areas_of_interest_pos[j][1] = exon_position[i][1]+margin
                i+=1
                j+=1
                areas_of_interest_pos.append([((exon_position[i][0])-(margin)),((exon_position[i][1])+(margin)),True])
    
    return (areas_of_interest_pos)

def load_gene(gene_filename):

    """!
    Realiza la carga del gen, para lo cual es necesario que el usuario ingrese el nombre
    del archivo donde está contenido el gen. Esta función arma un string con todas las letras
    del archivo, para asi poder trabajar con un string "limpio" (sin numeros ni espacios).
    Al mismo tiempo va almacenanando la posicion de inicio y fin de cada uno de los exones
    en una lista (exon_position), estas posiciones estan en relación al string que contiene
    el gen completo (full_gene).
    \param <genes_dir> {String, tiene la ruta hacia la carpeta de genes, con el nombre del archivo que contiene el gen}
    \param <full_gene> {String, almacena las letras contenidas en el archivo del gen}
    \param <exon_position> {lista, contiene la posicion inicial y final de cada exon}
    """ 
    script_dir = os.path.abspath(__file__)
    last_backslash_index = script_dir.rfind('/')
    #*****************************************************************
    #genes_dir = script_dir[0:last_backslash_index]+"/Genes/"+gene_filename 
    #genes_dir = script_dir[0:37]+"\\Genes\\"+gene_filename
    genes_dir = script_dir[0:last_backslash_index]+"/Genes/"+gene_filename
    doc_gene = docx.Document(genes_dir)

    global exon_position
    global full_gene
    aux_position = 0
    inside_exon = False
    largo_gen = 0
    last_run_color = ""
    count_nums_spaces = 0

    for paragraph in doc_gene.paragraphs:

        full_gene = full_gene + paragraph.text
        
        for run in paragraph.runs:
                                                                                        #3300FF PURPLE
                                                                                        #22CCEE LIGHT BLUE
                                                                                        #FF0033 RED
                                                                                        
            #The letter found is light blue, which means that it can be the start or the end of an exon
            if ( (str(run.font.color.rgb) == "22CCEE") or ( (str(run.font.color.rgb) == "3300FF") and last_run_color == "FF0033" ) or ( (str(run.font.color.rgb) == "FF0033") and last_run_color == "3300FF" ) ):

                if (inside_exon):
                    exon_position.append([aux_position,largo_gen])
                else:
                    aux_position = largo_gen

                inside_exon = not inside_exon

            if (str(run.text).isdigit() or ( " " in str(run.text))):
                count_nums_spaces+=1

            else:
                last_run_color = str(run.font.color.rgb)

            largo_gen = largo_gen + len(run.text.replace(" ","").replace("0","").replace("1","").replace("2","").replace("3","").replace("4","").replace("5","").replace("6","").replace("7","").replace("8","").replace("9",""))
    
    full_gene = full_gene.replace(" ","").replace("0","").replace("1","").replace("2","").replace("3","").replace("4","").replace("5","").replace("6","").replace("7","").replace("8","").replace("9","")
    #print("Exons positions")
    #print(exon_position)

    return ()

def primers_search(primer_min_length, primer_max_length, min_gc_percent, max_gc_percent, areas_of_interest):
    """!
    Funcion que se encarga de encontrar todos los primers posibles y almacenarlos en una lista (primer_list)
    que será posteriormente utilizada en la selección de los primers para cada sección.
    
    \param <primer_min_length> {Numero, indica el largo mínimo que deben tener los primers}
    \param <primer_max_length> {Numero, indica el largo máximo que deben tener los primers}
    \param <min_gc_percent> {Numero, indica el porcentaje minimo de GC que debe tener un primer}
    \param <max_gc_percent> {Numero, indica el porcentaje maximo de GC que debe tener un primer}
    \param <areas_of_interest> {Lista, contiene listas con el inicio y fin de cada zona de interes}

    """
    gene_segment = ""
    search_margin = var['largo_amplicon']
    j = 0
    tm_primer = 0
    gc_percent = 0
    num_homopol_simple = 1
    max_homopol_simple = 1
    num_double_homopol = 0
    max_homopol = 0
    dg_hairpin = 0
    dg_homodimero = 0
    mvconc=var["mv_conc"]

    for i in range(len(full_gene)-primer_max_length+2):

        if(i >= (areas_of_interest[j][0]-search_margin) and i<=(areas_of_interest[j][1]+search_margin) ):
            gene_segment = full_gene[i:i+primer_max_length]
        
            for k in range ((primer_min_length), primer_max_length+1):
                max_homopol_simple = 1
                num_homopol_simple = 1
                num_double_homopol = 0
                max_homopol_double = 0
                max_homopol = 0
                tm_primer=calcuTm(gene_segment[:k])
                #tm_primer = gene_segment[:k].count('a')*2+gene_segment[:k].count('A')*2+gene_segment[:k].count('t')*2+gene_segment[:k].count('T')*2+gene_segment[:k].count('c')*4+gene_segment[:k].count('C')*4+gene_segment[:k].count('g')*4+gene_segment[:k].count('G')*4
                gc_percent = (gene_segment[:k].count('c') + gene_segment[:k].count('C') + gene_segment[:k].count('g') + gene_segment[:k].count('G')) /k
                dg_hairpin = primer3.calcHairpin(gene_segment[:k],mv_conc=mvconc).dg
                dg_homodimero = primer3.calcHomodimer(gene_segment[:k],mv_conc=mvconc).dg
                for m in range(len(gene_segment[:k])-1):

                    if ( (m <= (len(gene_segment[:k])-3)) and (gene_segment[m] == gene_segment[m+2]) and (gene_segment[m] != gene_segment[m+1]) ):
                        num_double_homopol +=1
                    else:
                        num_double_homopol = 1

                    if (gene_segment[m] == gene_segment[m+1]):
                        num_homopol_simple +=1
                    else:
                        num_homopol_simple = 1

                    if(num_homopol_simple > max_homopol_simple):
                        max_homopol_simple = num_homopol_simple
                    
                    if(num_double_homopol > max_homopol_double):
                        max_homopol_double = num_double_homopol

                if(max_homopol_simple > math.floor(max_homopol_double/2)):
                    max_homopol = max_homopol_simple
                else:
                    max_homopol = math.floor(max_homopol_double/2)

                if ( (gc_percent >= min_gc_percent) and (gc_percent<= max_gc_percent) ):
                    primer_list.append([gene_segment[:k], i, len(gene_segment[:k]), max_homopol, gc_percent, tm_primer, True, dg_hairpin, dg_homodimero])

        if (i>(areas_of_interest[j][1]+search_margin) and j < len(areas_of_interest)-1):
            j+=1

    return (primer_list)

def seleccion_partidores(zonas_de_interes,rel,recomendador_largo):
    """!
    Función que realiza la selección de los partidores para las zonas previamente especificadas.
    
    \param <zonas_de_interes> {Lista, contiene todas las zonas para las cuales se seleccionarán partidores.}
    \param <primer_list> {Lista, variable global que contiene todos los posibles primers candidatos para ser seleccionados.}
    \param <resultados_primers> {Lista, contiene los resultados obtenidos posterior a la seleccion de partidores.}
    \param <rel> {Objeto, inicializa la clase Relajo.}
    \param <en_zona_de_interes> {Booleano, indica si en la seleccion aún se encuentra dentro de una zona de interes especifica.}
    \param <compatibilidad_FR> {Booleano, indica si los partidores forward y reverse son compatibles entre sí.}
    \param <relajo_valido> {Booleano, indica si el relajo es valido o no, el cual se vuelve falso en caso de no encontrar partidores usando el relajo máximo.}
    \param <restart_R> {Booleano, reinicio de la busqueda para el partidor reverse.}
    \param <heterodim_par_primers> {Numero, contiene el delta G asociado a la pareja de partidores.}
    \param <factibilidad_f> {Booleano, indica si el partidor forward es factible para su posterior selección en relacion a los valores del diccionario (var).}
    \param <factibilidad_r> {Booleano, indica si el partidor reverse es factible para su posterior selección en relacion a los valores del diccionario (var).}
    \param <pos_busqueda_f> {Numero, indice a la lista de partidores que muestra la posición desde donde se iniciará la busqueda de partidor forward.}
    \param <pos_busqueda_f_final> {Numero, indice a la lista de partidores que muestra la posicion donde finalizará la busqueda de partidor forward.}
    \param <pos_busqueda_r> {Numero, indice a la lista de partidores que muestra la posición desde donde se iniciará la busqueda de partidor reverse.}
    \param <pos_busqueda_r_final> {Numero, indice a la lista de partidores que muestra la posicion donde finalizará la busqueda de partidor reverse.}
    \param <puntero_forward> {Numero, indice a la lista de partidores que hace referencia a donde esta en ese momento la busqueda del partidor forward.}
    \param <puntero_reverse> {Numero, indice a la lista de partidores que hace referencia a donde esta en ese momento la busqueda del partidor reverse.}
    \param <rev_en_zona> {Booleano, indica si la busqueda del partidor reverse esta aun dentro de la zona que se desea cubrir.}
    \param <extrayendo_porcion_final> {Booleano, indica si se esta extrayendo la ultima porción dentro del archivo.}

    """
    numero_amplicon = 1
    zona_actual = 0
    json_para_grafico = {}
    json_amplicones_final = ""
    global primer_list
    largo_amplicon_actual = 0
    resultados = []
    resultados_primers = []
    resultados_txt = open("Resultados.txt","w")
    #rel = Relajo()
    en_zona_de_interes = True
    compatibilidad_FR = False
    relajo_valido = True
    restart_R = False
    heterodim_par_primers = 0
    factibilidad_f = False
    factibilidad_r = False
    pos_busqueda_f = 0
    pos_busqueda_f_final = 0
    pos_busqueda_r = 0
    pos_busqueda_r_final = 0
    puntero_forward = 0
    puntero_reverse = 0
    rev_en_zona = True
    relajo_valido = True
    extrayendo_porcion_final = False
    ultimo_inicio_de_zona = 0
    zona_sin_diseño_de_partidores = True
    n = 0
    m = 0
    mvconc=var["mv_conc"]

    for i in range(len(zonas_de_interes)):
        zona_actual = i
        restart_R = False
        en_zona_de_interes = True
        compatibilidad_FR = False
        relajo_valido = True
        zona_sin_diseño_de_partidores = zonas_de_interes[i][2]
        inicio_zona_a_extraer = zonas_de_interes[i][0]-var['primer_max_length']
        fin_zona_a_extraer = zonas_de_interes[i][1]
        #pdb.set_trace()
        while (en_zona_de_interes and relajo_valido and zona_sin_diseño_de_partidores):
            print("inicio de zona: "+str(inicio_zona_a_extraer)+" ultimo inicio zona: "+str(ultimo_inicio_de_zona))
            if (ultimo_inicio_de_zona != inicio_zona_a_extraer):
                rel.reset_relajo()
                print("+++++++++++++++++++++++SE RESETEO EL RELAJO+++++++++++++++++++++++")

            relajo_valido = True
            factibilidad_f = False
            factibilidad_r = False
            compatibilidad_FR = False
            restart_R = False
            #m = 0
            #while (primer_list[m][1] < primer_list[puntero_forward][1]+var['largo_amplicon']):
            #    m += 1
            #pos_busqueda_f = m-1
            #n = 0
            #while (primer_list[n][1] < primer_list[puntero_forward][1]+var['largo_amplicon']-margen_busqueda):
            #    n +=1
            #pos_busqueda_f_final = n-1
            #m = 0
            #n = 0
            if ( ((inicio_zona_a_extraer+(var['largo_amplicon'])) >= len(full_gene)) or extrayendo_porcion_final ):
                for i in range(len(primer_list)):
                    if (primer_list[i][1] == fin_zona_a_extraer-var['largo_amplicon']+var['primer_max_length']):
                        pos_busqueda_f_final = i
                        print("POS_BUSQUEDA_F_FINAL: "+str(primer_list[pos_busqueda_f_final][1])+" FIN ZONA A EXTRAER: "+str(fin_zona_a_extraer))
                        break

                if ( (inicio_zona_a_extraer-primer_list[pos_busqueda_f_final][1]) > var['margen_busqueda'] ):
                    for i in range(len(primer_list)):
                        if ( primer_list[i][1] == (primer_list[pos_busqueda_f_final][1]+var['margen_busqueda']) ):
                            pos_busqueda_f = i
                else:
                    for i in range(len(primer_list)):
                        if ( primer_list[i][1] == inicio_zona_a_extraer ):
                            pos_busqueda_f = i
                
                extrayendo_porcion_final = True
                print("DEFINIENDO ZONA PARA PORCION FINAL DE FORWARD")
            else:
                for i in range(len(primer_list)):
                    if (primer_list[i][1] == inicio_zona_a_extraer):
                        pos_busqueda_f = i
                
                for i in range(len(primer_list)):
                    if (primer_list[i][1] == primer_list[pos_busqueda_f][1]-var['margen_busqueda']):
                        pos_busqueda_f_final = i
                        print("la posicion de busqueda f final es: "+str(primer_list[pos_busqueda_f_final][1]))
                        break


            puntero_forward = pos_busqueda_f
            factibilidad_f = validacion_primer(puntero_forward)

            while (compatibilidad_FR == False and relajo_valido):

                while ( ((factibilidad_f == False) and (puntero_forward > pos_busqueda_f_final) and relajo_valido ) or restart_R):

                    restart_R = False
                    puntero_forward = puntero_forward-1
                    factibilidad_f = validacion_primer(puntero_forward)
                    ###print("PUNTERO FORWARD: "+str(puntero_forward)+" Posicion_busqueda_f_final: "+str(pos_busqueda_f_final))
                    
                    if (puntero_forward <= pos_busqueda_f_final):
                        relajo_valido = rel.relajar()
                        print("----------------------------RELAJE EN EL FORWARD--------------------------------")
                        print("EL INDICE DE RELAJO ES: "+str(rel.indice))
                        print("RELAJO VALIDO ES: "+str(relajo_valido))
                        puntero_forward = pos_busqueda_f
                        factibilidad_f = validacion_primer(puntero_forward)
                
                print("El relajo post forward es: "+str(rel.indice)+" Y puntero forward es: "+str(puntero_forward)+" Pos forward: "+str(primer_list[puntero_forward][1]))
                print("POS_BUSQUEDA_F_FINAL: "+str(pos_busqueda_f_final)+" POS_BUSQUEDA_F: "+str(pos_busqueda_f)+" DESDE: "+str(primer_list[pos_busqueda_f])+" HASTA: "+str(primer_list[pos_busqueda_f_final]))
                factibilidad_r = False
                compatibilidad_FR = False
                rev_en_zona = True
                #m = 0
                #while (primer_list[m][1] < primer_list[puntero_forward][1]+var['largo_amplicon']):
                #    m += 1
                #pos_busqueda_r = m-1
                #n = 0
                #while (primer_list[n][1] < primer_list[puntero_forward][1]+var['largo_amplicon']-margen_busqueda):
                #    n +=1
                #pos_busqueda_r_final = n-1
                #m = 0
                #n = 0
                if (extrayendo_porcion_final):
                    for i in range(len(primer_list)):
                        if ( primer_list[i][1] == fin_zona_a_extraer ):
                            pos_busqueda_r_final = i
                            break

                    for i in range(len(primer_list)):
                        #print("EL MIN ES: "+str(min( (len(full_gene)-var['primer_max_length']),( (primer_list[puntero_forward][1])+var['largo_amplicon'] ) ) ))
                        if ( primer_list[i][1] == min( (len(full_gene)-var['primer_max_length']),( (primer_list[puntero_forward][1])+var['largo_amplicon']-var['primer_max_length'] ) ) ):
                            pos_busqueda_r = i
                            print("EL MINIMO ENTRE FULL GENE Y FORWARD ES: "+str(min( (len(full_gene)-var['primer_max_length']),( (primer_list[puntero_forward][1])+var['largo_amplicon'] ) ) ))
                            print("DONDE FULL GENE ES: "+str((len(full_gene)-var['primer_max_length']))+" Y FORWARD + LARGO ES: "+str(( (primer_list[puntero_forward][1])+var['largo_amplicon'] )))

                    print("DEFINIENDO ZONA PARA PORCION FINAL DE REVERSE CON R FINAL: "+str(primer_list[pos_busqueda_r_final][1])+" Y R: "+str(primer_list[pos_busqueda_r][1]))
                    print("EL INICIO DE ZONA A EXTRAER ES: "+str(inicio_zona_a_extraer)+" Y EL FIN DE ZONA: "+str(fin_zona_a_extraer))
                else:       
                    for i in range(len(primer_list)):
                        if (primer_list[i][1] == primer_list[puntero_forward][1]+var['largo_amplicon']-var['primer_max_length']):
                            pos_busqueda_r = i
                            print("0000000000000000000000----- DEFINIENDO POS BUSQUEDA R ------000000000000000000")

                    for i in range(len(primer_list)):        
                        if (primer_list[i][1] == primer_list[puntero_forward][1]+var['largo_amplicon']-var['margen_busqueda']):
                            pos_busqueda_r_final = i
                            break

                print("El forward pos: "+str(primer_list[puntero_forward][1])+" pos_busqueda_r: "+str(primer_list[pos_busqueda_r][1])+" pos_busqueda_r_final: "+str(primer_list[pos_busqueda_r_final][1]))
                puntero_reverse = pos_busqueda_r
                factibilidad_r = validacion_primer(puntero_reverse)
                compatibilidad_FR = compatibilidad(primer_list[puntero_forward][0],primer_list[puntero_reverse][0])

                while ( (factibilidad_r == False or compatibilidad_FR == False) and rev_en_zona and relajo_valido ):

                    puntero_reverse = puntero_reverse-1
                    factibilidad_r = validacion_primer(puntero_reverse)
                    compatibilidad_FR = compatibilidad(primer_list[puntero_forward][0],primer_list[puntero_reverse][0])
                    if (puntero_reverse <= pos_busqueda_r_final):

                        rev_en_zona = False
                        factibilidad_r = False
                        factibilidad_f = False
                        compatibilidad_FR = False
                        puntero_forward = puntero_forward-1
                        print("0000000000000000000000----- PUNTERO FORWARD EN SELECCION REVERSE: "+str(puntero_forward))
                        restart_R = True

            if ( (primer_list[puntero_forward][1] < primer_list[puntero_reverse][1]) and (factibilidad_f == True) and (compatibilidad_FR == True) and (factibilidad_r == True) and ( ((primer_list[puntero_reverse][1]) - (primer_list[puntero_forward][1])) >= (var['largo_amplicon'] - var['margen_busqueda'])) and relajo_valido ):
                
                resultados.append("Forward: "+primer_list[puntero_forward][0]+" Posicion: "+str(primer_list[puntero_forward][1])+" Reverse: "+primer_list[puntero_reverse][0]+" Posicion: "+str(primer_list[puntero_reverse][1]))
                print(resultados)
                                
                resultados_txt.write("Forward: "+primer_list[puntero_forward][0]+" Tm: "+str(primer_list[puntero_forward][5])+" Porcentaje_GC: "+str(primer_list[puntero_forward][4])+" Reverse: "+reverse_complement(primer_list[puntero_reverse][0])+" Tm: "+str(primer_list[puntero_reverse][5])+" Porcentaje_GC: "+str(primer_list[puntero_reverse][4])+"\n")
                amplicon = full_gene[(primer_list[puntero_forward][1]):((primer_list[puntero_reverse][1])+(primer_list[puntero_reverse][2]))]
                resultados_txt.write("Largo amplicón: "+str((primer_list[puntero_reverse][1])+(primer_list[puntero_reverse][2])-(primer_list[puntero_forward][1]))+"\n")
                
                if (len(amplicon)>0):
                    amplicon_gc = (amplicon.count('c')+amplicon.count('C')+amplicon.count('g')+amplicon.count('G'))/len(amplicon)
                else:
                    amplicon_gc = 0

                heterodim_par_primers = primer3.calcHeterodimer(primer_list[puntero_forward][0], primer_list[puntero_reverse][0],mv_conc=mvconc).dg
                resultados_primers.append([primer_list[puntero_forward][0],primer_list[puntero_forward][1],primer_list[puntero_forward][5],primer_list[puntero_forward][4],primer_list[puntero_reverse][0],primer_list[puntero_reverse][1],primer_list[puntero_reverse][5],primer_list[puntero_reverse][4],amplicon_gc,rel.indice,primer_list[puntero_forward][7],primer_list[puntero_forward][8],primer_list[puntero_reverse][7],primer_list[puntero_reverse][8],heterodim_par_primers,((primer_list[puntero_reverse][1])+(primer_list[puntero_reverse][2])-(primer_list[puntero_forward][1]))])
                resultados_txt.write("Porcentaje GC amplicón: "+str(amplicon_gc)+"\n")
                largo_amplicon_actual = ((primer_list[puntero_reverse][1])+(primer_list[puntero_reverse][2])-(primer_list[puntero_forward][1]))
                
                if(var['aplicar_intervalo']):
                    recomendador_largo.reportar_largo_encontrado(largo_amplicon_actual)
                    recomendador_largo.proximo_largo_objetivo()

                json_para_grafico['type'] = 'amplicon'
                json_para_grafico['name'] = str(numero_amplicon)
                json_para_grafico['begin'] = primer_list[puntero_forward][1]+1
                json_para_grafico['end'] = primer_list[puntero_forward][1]+largo_amplicon_actual
                #for v in range(len(zonas_de_interes)):
                #    if ( ( (primer_list[puntero_forward][1] < zonas_de_interes[v][1] ) and ( primer_list[puntero_reverse][1] > zonas_de_interes[v][0]) ) ):
                #        json_para_grafico['zone'] = v
                json_para_grafico['zone'] = zona_actual
                if (numero_amplicon != 1):
                    json_amplicones_final = json_amplicones_final+ "," + json.dumps(json_para_grafico)
                else:
                    json_amplicones_final = json.dumps(json_para_grafico)
                
                numero_amplicon+=1

            if( (factibilidad_f == True) and (compatibilidad_FR == True) and (factibilidad_r == True) and ( ((primer_list[puntero_reverse][1]) - (primer_list[puntero_forward][1])) >= (var['largo_amplicon'] - var['margen_busqueda']) ) ):
                
                print("forward: "+str(primer_list[puntero_forward][6])+" reverse: "+str(primer_list[puntero_reverse][6]))
                
                for i in range(len(primer_list)):
                    
                    if ( primer_list[i][1] > ((primer_list[puntero_forward][1])-var['traslape']) and primer_list[i][1] < (primer_list[puntero_forward][1])+var['traslape']+primer_list[puntero_forward][2] ):
                        primer_list[i][6] = False
                        ###print("Se borro el primer: "+primer_list[i][0]+" posición: "+str(primer_list[i][1])+" y su borrado esta en: "+str(primer_list[i][6]))
                    
                    if ( primer_list[i][1] > (primer_list[puntero_reverse][1])-var['traslape']-(primer_list[puntero_reverse][2]) and primer_list[i][1] < (primer_list[puntero_reverse][1])+var['traslape']+primer_list[puntero_reverse][2] ):
                        primer_list[i][6] = False
                
                print("EL RELAJO ES: "+str(rel.indice))
                #pdb.set_trace()
                print("forward: "+str(primer_list[puntero_forward][6])+" reverse: "+str(primer_list[puntero_reverse][6]))
                
            ultimo_inicio_de_zona = inicio_zona_a_extraer

            if( (primer_list[puntero_reverse][1]> fin_zona_a_extraer) and (factibilidad_f == True) and (compatibilidad_FR == True) and (factibilidad_r == True) ):
                en_zona_de_interes = False
                print("SALI DE LA ZONA DE INTERES")
            elif((primer_list[puntero_forward][1] < primer_list[puntero_reverse][1]) and factibilidad_f and factibilidad_r and compatibilidad_FR):
                print("F SELECCIONADO EN POS: "+str(primer_list[puntero_forward][1])+" R SELECCIONADO EN POS: "+str(primer_list[puntero_reverse][1]))
                inicio_zona_a_extraer = primer_list[puntero_reverse][1]- len(primer_list[puntero_reverse][0]) - var['traslape']

    resultados_txt.close()
    print(resultados)
    for i in range(len(resultados)):
        print(resultados[i])
    
    json_amplicones_txt = open("JSON_amplicones.txt","w")
    json_amplicones_txt.write(json_amplicones_final)
    json_amplicones_txt.close()
    rel.reset_relajo()
    print("JSON AMPLICONES FINAL")
    print(json_amplicones_final)
    print("RESULTADO PRIMERS")
    print(resultados_primers)
    return resultados_primers

def reverse_complement(reverse):
    """!
    Funcion que recibe un partidor cualquiera y le genera su complemento. Para el caso
    de este proyecto, es utilizado para obtener el complemento de los partidores reverse.
    
    \param <reverse> {String, es el partidor que recibe la función}
    \param <reverse_result> {String, que contiene el complemento del partidor ingresado}

    """

    reverse_result = ""
    for i in range(len(reverse)):
        if ( reverse[(len(reverse)-1-i)] == "A" or reverse[(len(reverse)-1-i)] == "a" ):
            reverse_result = reverse_result+"T"
        if ( reverse[(len(reverse)-1-i)] == "T" or reverse[(len(reverse)-1-i)] == "t" ):
            reverse_result = reverse_result+"A"
        if ( reverse[(len(reverse)-1-i)] == "G" or reverse[(len(reverse)-1-i)] == "g" ):
            reverse_result = reverse_result+"C"
        if ( reverse[(len(reverse)-1-i)] == "C" or reverse[(len(reverse)-1-i)] == "c" ):
            reverse_result = reverse_result+"G"
    return (reverse_result)

def compatibilidad(forward, reverse):
    """!
    Función que se encarga de verificar si dos partidores son compatibles entre sí.
    Esta verificación se realiza accediendo a los valores de ambos partidores, para
    luego ser comparados con los valores presentes en el diccionario de variables.
    
    \param <forward> {String, partidor forward que recibe la función}
    \param <reverse> {String, partidor reverse que recibe la función}
    \param <tm_forward> {Numero, temperatura de melting asociada al partidor forward ingresado}
    \param <tm_reverse> {Numero, temperatura de melting asociada al partidor reverse ingresado}
    \param <heterodimero_forward_reverse> {Numero, delta G asociado a la posibilidad de que ambos partidores generen un heterodimero}

    """
    mvconc=var["mv_conc"]
    #tm_forward = forward.count('a')*2+forward.count('A')*2+forward.count('t')*2+forward.count('T')*2+forward.count('c')*4+forward.count('C')*4+forward.count('g')*4+forward.count('G')*4
    #tm_reverse = reverse.count('a')*2+reverse.count('A')*2+reverse.count('t')*2+reverse.count('T')*2+reverse.count('c')*4+reverse.count('C')*4+reverse.count('g')*4+reverse.count('G')*4
    tm_forward=tm_primer=calcuTm(forward)
    tm_reverse=tm_primer=calcuTm(reverse)
    heterodimero_forward_reverse = primer3.calcHeterodimer(forward, reverse,mv_conc=mvconc).dg
    if ( abs(tm_forward-tm_reverse) <= var['delta_tm_max'] and abs(tm_forward-tm_reverse) >= var['delta_tm_min']  and (heterodimero_forward_reverse > var['deltag']) ):
        return True
    else:
        return False

def validacion_primer(primer):
    """!
    Funcion que se encarga de verificar si un partidor es valido para su posterior selección.
    Para realizar esta validación, se utilizan los valores presentes en el diccionario de variables (var).
    
    \param <primer_list> {Lista, variable global que contiene todos los partidores encontrados}
    \param <primer> {Numero, es el índice asociado a un partidor especifico dentro la lista de partidores(primer_list)}

    """
    global primer_list
    
    if ( (primer_list[primer][3] > var['max_homopol']) or (primer_list[primer][4] > var['max_gc_percent']) or (primer_list[primer][4] < var['min_gc_percent']) or (primer_list[primer][5] < var['tm_min']) or (primer_list[primer][6] == False) or (primer_list[primer][7] < var['deltag']) or (primer_list[primer][8] < var['deltag'])):
        return False
    elif var["focus_tm"]>=0 and abs(primer_list[primer][5]-var["focus_tm"]) > var["focus_d"]:
        return False
    else:
        return True

class RecomendadorDeLargo:

    largo_amplicon = None
    largo_arreglo = None
    num_tubos = None

    def __init__(self, largo_exoma, num_exones, margen_intronico, largo_min_amplicon, largo_max_amplicon):
        self.largo_amplicon = []
        self.largo_arreglo = math.floor((largo_max_amplicon - largo_min_amplicon) / (var['intervalo']))
        print("Largo max amplicon: "+str(largo_max_amplicon)+" Largo min amplicon: "+str(largo_min_amplicon)+" Intervalo: "+str(var['intervalo'])+" largo arreglo: "+str(self.largo_arreglo))
        largo_dna_en_tubo = 0
        print("Largo Arreglo: "+str(self.largo_amplicon))
        for i in range(self.largo_arreglo):
            self.largo_amplicon.append([(largo_max_amplicon-(i* (var['intervalo']) )),0])
            largo_dna_en_tubo += self.largo_amplicon[i][0]
        
        self.num_tubos = math.ceil( (largo_exoma + (num_exones*2*margen_intronico))/largo_dna_en_tubo )
        print("Largo arreglo: "+str(self.largo_arreglo)+" Largo amplicon len: "+str(len(self.largo_amplicon)))
    def proximo_largo_objetivo(self):
        i = 0
        num_tubos_superado = False
        suma_largos_totales = 0

        for k in range(len(self.largo_amplicon)):
            suma_largos_totales = suma_largos_totales + self.largo_amplicon[k][1]

        if (suma_largos_totales == (self.num_tubos*len(self.largo_amplicon) ) ):
            print("AUMENTAR NUMERO DE TUBOS: "+str(self.num_tubos))
            self.num_tubos += 1
            print("AUMENTADO A: "+str(self.num_tubos))

        while (i < self.largo_arreglo):
            
            if (self.largo_amplicon[i][1] < self.num_tubos):
                break

            print("i antes: "+str(i))
            i+=1
            print("i despues: "+str(i))

        print("Largo arreglo: "+str(self.largo_arreglo)+" Largo amplicon len: "+str(len(self.largo_amplicon))+" i: "+str(i))
        print(self.largo_amplicon)
        print(str(self.num_tubos))
        var['margen_busqueda'] = (self.largo_arreglo * var['intervalo']) - (i * var['intervalo'])
        var['largo_amplicon'] = self.largo_amplicon[i][0]
        #return (self.largo_amplicon[i][0])
    
    def reportar_largo_encontrado(self, largo_amplicon_encontrado):
        min_distancia = 9999999
        min_indice = 0
        for i in range(self.largo_arreglo):
            if ( abs(largo_amplicon_encontrado - self.largo_amplicon[i][0]) < min_distancia):
                min_distancia = abs(largo_amplicon_encontrado - self.largo_amplicon[i][0])
                min_indice = i
        
        self.largo_amplicon[min_indice][1] +=1


class Relajo:
    
    indice=None
    valores_de_relajo = []
    lista_variables = []

    def __init__(self):
        """!
        Funcion que inicializa la clase relajo con todos sus valores ingresados.
    
        \param <valores_de_relajo> {Lista, contiene los valores que seran utilizados para disminuir la rigurosidad con la que son seleccionados los partidores}
        \param <indice> {Numero, indice asociado a la lista valores_de_relajo el cual permite ir relajando los valores a medida que se va aumentado su valor}

        """
        self.indice=0
        #self.valores_de_relajo = [[0, 4, 0.40, 0.6, 4],
        #                          [0, 4, 0.35, 0.7, 4],
        #                          [0, 4, 0.30, 0.8, 5],
        #                          [0, 4, 0.30, 0.8, 6],
        #                          [0, 100, 0.1, 1.0, 10]]
        # self.valores_de_relajo = [[0, 0, 0.4, 0.6, 4],
        #                           [0, 2, 0.4, 0.6, 4],
        #                           [0, 4, 0.4, 0.6, 4],
        #                           [0, 0, 0.4, 0.7, 4],
        #                           [0, 2, 0.4, 0.7, 4],
        #                           [0, 4, 0.4, 0.7, 4],
        #                           [0, 0, 0.35, 0.7, 4],
        #                           [0, 2, 0.35, 0.7, 4],
        #                           [0, 4, 0.35, 0.7, 4],
        #                           [0, 0, 0.35, 0.75, 4],
        #                           [0, 2, 0.35, 0.75, 4],
        #                           [0, 4, 0.35, 0.75, 4],
        #                           [0, 0, 0.3, 0.75, 4],
        #                           [0, 2, 0.3, 0.75, 4],
        #                           [0, 4, 0.3, 0.75, 4],
        #                           [0, 0, 0.30, 0.80, 4],
        #                           [0, 2, 0.30, 0.80, 4],
        #                           [0, 4, 0.30, 0.80, 4],
        #                           [0, 0, 0.4, 0.6, 5],
        #                           [0, 2, 0.4, 0.6, 5],
        #                           [0, 4, 0.4, 0.6, 5],
        #                           [0, 0, 0.4, 0.7, 5],
        #                           [0, 2, 0.4, 0.7, 5],
        #                           [0, 4, 0.4, 0.7, 5],
        #                           [0, 0, 0.35, 0.7, 5],
        #                           [0, 2, 0.35, 0.7, 5],
        #                           [0, 4, 0.35, 0.7, 5],
        #                           [0, 0, 0.35, 0.75, 5],
        #                           [0, 2, 0.35, 0.75, 5],
        #                           [0, 4, 0.35, 0.75, 5],
        #                           [0, 0, 0.3, 0.75, 5],
        #                           [0, 2, 0.3, 0.75, 5],
        #                           [0, 4, 0.3, 0.75, 5],
        #                           [0, 0, 0.30, 0.80, 5],
        #                           [0, 2, 0.30, 0.80, 5],
        #                           [0, 4, 0.30, 0.80, 5],
        #                           [0, 0, 0.4, 0.6, 6],
        #                           [0, 2, 0.4, 0.6, 6],
        #                           [0, 4, 0.4, 0.6, 6],
        #                           [0, 0, 0.4, 0.7, 6],
        #                           [0, 2, 0.4, 0.7, 6],
        #                           [0, 4, 0.4, 0.7, 6],
        #                           [0, 0, 0.35, 0.7, 6],
        #                           [0, 2, 0.35, 0.7, 6],
        #                           [0, 4, 0.35, 0.7, 6],
        #                           [0, 0, 0.35, 0.75, 6],
        #                           [0, 2, 0.35, 0.75, 6],
        #                           [0, 4, 0.35, 0.75, 6],
        #                           [0, 0, 0.3, 0.75, 6],
        #                           [0, 2, 0.3, 0.75, 6],
        #                           [0, 4, 0.3, 0.75, 6],
        #                           [0, 0, 0.30, 0.80, 6],
        #                           [0, 2, 0.30, 0.80, 6],
        #                           [0, 4, 0.30, 0.80, 6],
        #                           [0, 0, 0.4, 0.6, 7],
        #                           [0, 2, 0.4, 0.6, 7],
        #                           [0, 4, 0.4, 0.6, 7],
        #                           [0, 0, 0.4, 0.7, 7],
        #                           [0, 2, 0.4, 0.7, 7],
        #                           [0, 4, 0.4, 0.7, 7],
        #                           [0, 0, 0.35, 0.7, 7],
        #                           [0, 2, 0.35, 0.7, 7],
        #                           [0, 4, 0.35, 0.7, 7],
        #                           [0, 0, 0.35, 0.75, 7],
        #                           [0, 2, 0.35, 0.75, 7],
        #                           [0, 4, 0.35, 0.75, 7],
        #                           [0, 0, 0.3, 0.75, 7],
        #                           [0, 2, 0.3, 0.75, 7],
        #                           [0, 4, 0.3, 0.75, 7],
        #                           [0, 0, 0.30, 0.80, 7],
        #                           [0, 2, 0.30, 0.80, 7],
        #                           [0, 4, 0.30, 0.80, 7],
        #                           [0, 0, 0.4, 0.6, 8],
        #                           [0, 2, 0.4, 0.6, 8],
        #                           [0, 4, 0.4, 0.6, 8],
        #                           [0, 0, 0.4, 0.7, 8],
        #                           [0, 2, 0.4, 0.7, 8],
        #                           [0, 4, 0.4, 0.7, 8],
        #                           [0, 0, 0.35, 0.7, 8],
        #                           [0, 2, 0.35, 0.7, 8],
        #                           [0, 4, 0.35, 0.7, 8],
        #                           [0, 0, 0.35, 0.75, 8],
        #                           [0, 2, 0.35, 0.75, 8],
        #                           [0, 4, 0.35, 0.75, 8],
        #                           [0, 0, 0.3, 0.75, 8],
        #                           [0, 2, 0.3, 0.75, 8],
        #                           [0, 4, 0.3, 0.75, 8],
        #                           [0, 0, 0.30, 0.80, 8],
        #                           [0, 2, 0.30, 0.80, 8],
        #                           [0, 4, 0.30, 0.80, 8],
        #                           [0, 0, 0.4, 0.6, 9],
        #                           [0, 2, 0.4, 0.6, 9],
        #                           [0, 4, 0.4, 0.6, 9],
        #                           [0, 0, 0.4, 0.7, 9],
        #                           [0, 2, 0.4, 0.7, 9],
        #                           [0, 4, 0.4, 0.7, 9],
        #                           [0, 0, 0.35, 0.7, 9],
        #                           [0, 2, 0.35, 0.7, 9],
        #                           [0, 4, 0.35, 0.7, 9],
        #                           [0, 0, 0.35, 0.75, 9],
        #                           [0, 2, 0.35, 0.75, 9],
        #                           [0, 4, 0.35, 0.75, 9],
        #                           [0, 0, 0.3, 0.75, 9],
        #                           [0, 2, 0.3, 0.75, 9],
        #                           [0, 4, 0.3, 0.75, 9],
        #                           [0, 0, 0.30, 0.80, 9],
        #                           [0, 2, 0.30, 0.80, 9],
        #                           [0, 4, 0.30, 0.80, 9],
        #                           [0, 0, 0.4, 0.6, 10],
        #                           [0, 2, 0.4, 0.6, 10],
        #                           [0, 4, 0.4, 0.6, 10],
        #                           [0, 0, 0.4, 0.7, 10],
        #                           [0, 2, 0.4, 0.7, 10],
        #                           [0, 4, 0.4, 0.7, 10],
        #                           [0, 0, 0.35, 0.7, 10],
        #                           [0, 2, 0.35, 0.7, 10],
        #                           [0, 4, 0.35, 0.7, 10],
        #                           [0, 0, 0.35, 0.75, 10],
        #                           [0, 2, 0.35, 0.75, 10],
        #                           [0, 4, 0.35, 0.75, 10],
        #                           [0, 0, 0.3, 0.75, 10],
        #                           [0, 2, 0.3, 0.75, 10],
        #                           [0, 4, 0.3, 0.75, 10],
        #                           [0, 0, 0.30, 0.80, 10],
        #                           [0, 2, 0.30, 0.80, 10],
        #                           [0, 4, 0.30, 0.80, 10],
        #                           [0, 100, 0.1, 1.0, 10]]
        # self.valores_de_relajo = [[0, 0, 0.4, 0.6],
        #                          [0, 2, 0.4, 0.6],
        #                          [0, 4, 0.4, 0.6],
        #                          [0, 0, 0.4, 0.7],
        #                          [0, 2, 0.4, 0.7],
        #                          [0, 4, 0.4, 0.7],
        #                          [0, 0, 0.35, 0.7],
        #                          [0, 2, 0.35, 0.7],
        #                          [0, 4, 0.35, 0.7],
        #                          [0, 0, 0.35, 0.75],
        #                          [0, 2, 0.35, 0.75],
        #                          [0, 4, 0.35, 0.75],
        #                          [0, 0, 0.3, 0.75],
        #                          [0, 2, 0.3, 0.75],
        #                          [0, 4, 0.3, 0.75],
        #                          [0, 0, 0.30, 0.80],
        #                          [0, 2, 0.30, 0.80],
        #                          [0, 4, 0.30, 0.80],
        #                           [0, 100, 0.1, 1.0]]
        self.valores_de_relajo = [[0, 0.4, 0.6],
                                 [1, 0.35, 0.7],
                                 [2, 0.30, 0.8],
                                 [4, 0.2, 0.9],
                                 ]
        #self.lista_variables = ['delta_tm_min', 'delta_tm_max', 'min_gc_percent', 'max_gc_percent']
        self.lista_variables = ['delta_tm_max','min_gc_percent','max_gc_percent']

    def asigna_valores(self):
        """!
        Funcion que asigna los nuevos valores al diccionario (var) en relación al indice que se tiene en esa iteración.
    
        \param <lista_variables> {Lista, contiene los nombres de las variables a ser relajadas}
        \param <valores_de_relajo> {Lista, contiene los valores de relajo asociados a las variables}
        
        """
        for i in range(len(self.lista_variables)):
            var[self.lista_variables[i]] = self.valores_de_relajo[self.indice][i]

    def relajar(self):
        """!
        Funcion que aumenta el indice asociado a la clase relajo.
    
        \param <indice> {Numero, indice que hace referencia al nivel de relajo de la lista 'valores_de_relajo'}
        
        """
		
        if self.indice < len(self.valores_de_relajo)-1:
            self.indice +=1
            self.asigna_valores()
            return True
        else:
            return False
    
    def reset_relajo(self):
        """!
        Funcion que reinicia el indice para reiniciar los valores de relajo.
    
        \param <indice> {Numero, indice que hace referencia al nivel de relajo de la lista 'valores_de_relajo'}
        
        """
        self.indice = 0
        self.asigna_valores()

def generar_csv(resultados_amplicon,areas_a_extraer,nombre_archivo,rel):
    """!
    Funcion que genera el archivo .csv que entrega los resultados. En estos resultados se encuentran
    los valores asociados a los partidores de forma independiente como tambien con su respectiva pareja,
    y por otro lado estan presentes los valores asociados a los amplicones que estos deberían generar.
    
    \param <resultados_amplicon> {Lista, contiene todos los valores asociados a los resultados obtenidos luego de ejecutar el software.}
    \param <resultados_csv> {Archivo, variable con el archivo .csv.}

    """
    nombre_archivo = nombre_archivo+'.csv'
    reverse_complementado = ""
    #rel = Relajo()
    resultados_csv = open(nombre_archivo,'w')
    
    for key in var:
        if(key!="exones_seleccionados"):
            resultados_csv.write(key+" : "+str(var[key])+"\n")
    
    resultados_csv.write("\n")

    for i in range(len(rel.lista_variables)):
        if (i != (len(rel.lista_variables)-1) ):
            resultados_csv.write(rel.lista_variables[i]+",")
        else:
            resultados_csv.write(rel.lista_variables[i]+","+"Nivel Relajo"+"\n")

    for i in range(len(rel.valores_de_relajo)):
        for j in range(len(rel.valores_de_relajo[i])):
            if (j != (len(rel.valores_de_relajo[i])-1) ):
                resultados_csv.write(str(rel.valores_de_relajo[i][j])+",")
            else:
                resultados_csv.write(str(rel.valores_de_relajo[i][j])+","+str(i)+"\n")

    resultados_csv.write("\n")

    contador_exones_por_amplicon = 0
    contador_letra_por_exon = 65
    codigo_forward = ""
    codigo_reverse = ""
    delta_tm = 0

    resultados_csv.write("Código primer,Descripcion,Posicion de inicio,Largo,Tm,Porcentaje GC,Hairpin,Homodimero,DeltaTm,Porcentaje GC amplicón,Largo amplicon,Heterodimero,Nivel de relajo,Exones cubiertos \n")
    exones_por_amplicon = ""
    for j in range(len(resultados_amplicon)):
        exones_por_amplicon = ""
        contador_exones_por_amplicon = 0
        for i in range(len(exon_position)):
            #seccion para codigo forward
            if (i==0 and j == 0 ):
              
                codigo_forward = str(PKD1_exon_desde+i+1)

            elif ( resultados_amplicon[j][1] > exon_position[i][0] and resultados_amplicon[j][1] <= exon_position[i][1]):
            
                contador_letra_por_exon +=1
                codigo_forward = str(PKD1_exon_desde+i+1)+chr(contador_letra_por_exon)
            
            elif (i>0 and resultados_amplicon[j][1] < exon_position[i][0] and resultados_amplicon[j][1] > exon_position[i-1][1]):

                contador_letra_por_exon = 65
                codigo_forward = str(PKD1_exon_desde+i+1)
            
            #Seccion para codigo reverse
            #print("largo amplicones: "+str(len(resultados_amplicon))+" len de exon pos: "+str(len(exon_position))+" el i es: "+str(i)+" y el j: "+str(j))
            #i == (len(exon_position)-1) and
            if ( j == (len(resultados_amplicon)-1) and resultados_amplicon[j][5] > exon_position[(len(exon_position)-1)][1]   ):
                codigo_reverse = str(PKD1_exon_desde+i+1)

            elif (i>0 and resultados_amplicon[j][5] > exon_position[i][0] and resultados_amplicon[j][5] < exon_position[i][1] and resultados_amplicon[j][1] < exon_position[i-1][1]):
                
                contador_letra_por_exon = 65
                codigo_reverse = str(PKD1_exon_desde+i+1)+chr(contador_letra_por_exon)

            elif ( resultados_amplicon[j][5] > exon_position[i][0] and resultados_amplicon[j][5] < exon_position[i][1] ):

                codigo_reverse = str(PKD1_exon_desde+i+1)+chr(contador_letra_por_exon)

            elif ( i < (len(exon_position)-1) and resultados_amplicon[j][5] < exon_position[i+1][0] and resultados_amplicon[j][5] > exon_position[i][1] ):

                codigo_reverse = str(PKD1_exon_desde+i+1)
                contador_letra_por_exon = 65

            if(resultados_amplicon[j][1] < exon_position[i][1] and resultados_amplicon[j][5] > exon_position[i][0]):

                contador_exones_por_amplicon +=1
                exones_por_amplicon = exones_por_amplicon+"  "+(var['exones_seleccionados'][i]['label'])

        codigo_forward= str(j+1+contador_amplicones_PKD1)+" F"
        codigo_reverse= str(j+1+contador_amplicones_PKD1)+" R"
        print("codigo forward: "+codigo_forward+" Forward: "+resultados_amplicon[j][0]+" codigo reverse: "+codigo_reverse+" Reverse: "+resultados_amplicon[j][4])

        delta_tm = abs(resultados_amplicon[j][2]-resultados_amplicon[j][6])
        reverse_complementado = reverse_complement(resultados_amplicon[j][4])
        resultados_csv.write(codigo_forward+","+var['overhangF']+resultados_amplicon[j][0]+","+str(resultados_amplicon[j][1]+1)+","+str(len(resultados_amplicon[j][0]))+","+str(resultados_amplicon[j][2])+","+str(resultados_amplicon[j][3])+","+str(resultados_amplicon[j][10]/1000)+","+str(resultados_amplicon[j][11]/1000)+","+str(delta_tm)+","+str(resultados_amplicon[j][8])+","+str(resultados_amplicon[j][15]+len(var['overhangF']+var['overhangR']))+","+str(resultados_amplicon[j][14]/1000)+","+str(resultados_amplicon[j][9])+","+exones_por_amplicon[2:]+"\n")
        resultados_csv.write(codigo_reverse+","+var['overhangR']+reverse_complementado+","+str(resultados_amplicon[j][5]+1)+","+str(len(resultados_amplicon[j][4]))+","+str(resultados_amplicon[j][6])+","+str(resultados_amplicon[j][7])+","+str(resultados_amplicon[j][12]/1000)+","+str(resultados_amplicon[j][13]/1000)+","+str(delta_tm)+","+str(resultados_amplicon[j][8])+","+str(resultados_amplicon[j][15]+len(var['overhangF']+var['overhangR']))+","+str(resultados_amplicon[j][14]/1000)+","+str(resultados_amplicon[j][9])+","+exones_por_amplicon[2:]+"\n")


    resultados_csv.close    
    return

def verificacion_unicidad(primers_seleccionados, zonas_de_interes,rel,recomendador_largo):
    """!
    Funcion que realiza la verificación de si los partidores ingresados son únicos o no.
    
    \param <primers_seleccionados> {Lista, contiene los partidres y sus valores asociados luego de realizar la seleccion.}
    \param <zonas_de_interes> {Lista, contiene las zonas para las cuales fueron diseñados los partidores.}
    \param <blastn_cline> {String, contiene la linea de comando que será posteriormente ejecutada con 'blastn_cline().}
    \param <primers_verificados> {Lista, contiene los primers posterior a la verificación de unicidad.}

    """

    primers_verificados = []
    archivo_con_primers = open("listado_de_primers.txt", "w")

    for i in primers_seleccionados:

        archivo_con_primers.write("> \n")
        archivo_con_primers.write(i[0]+"\n")
        archivo_con_primers.write("> \n")
        archivo_con_primers.write(i[4]+"\n")

    archivo_con_primers.close()

    blastn_cline = NcbiblastnCommandline(query = "listado_de_primers.txt", num_descriptions=5, num_alignments=5, num_threads=4, out="resultadosXML.xml", word_size=11, outfmt=5, db = "PKD1")
    stdout, stderr = blastn_cline()

    resultados_xml = open("resultadosXML.xml")
    blast_records = NCBIXML.parse(resultados_xml)

    for b in blast_records:
        for alignment in b.alignments:
            for i in range (len(alignment.hsps)):
                if (i==1):
                    if ( alignment.hsps[i].expect/alignment.hsps[i-1].expect > 100):
                        print("e.value CORRECTO")
                    else:
                        print("Volviendo a ejecutar Seleccion")
                        primers_verificados = seleccion_partidores(zonas_de_interes,rel,recomendador_largo)
                
                print("***Alignment***")
                print(alignment.hsps[i].query)
                print("length: ", alignment.length)
                print("e value: ", alignment.hsps[i].expect)



    return primers_verificados

def key_to_sort(elem):
    
    return elem[0]

class Exon:
    def __init__(self, name, begin, end, zone):
        self.name = name
        self.begin = begin
        self.end = end
        self.zone = zone
    def Display(self):
        print ("[Exon] name:" + self.name + "  begin:" + str(self.begin) + "  end:" + str(self.end) + "  zone:" + str(self.zone))

class Amplicon:
    def __init__(self, name, begin, end, zone):
        self.name = name
        self.begin = begin
        self.end = end
        self.zone = zone
    def Display(self):
        print ("[Amplicon] name:" + self.name + "  begin:" + str(self.begin) + "  end:" + str(self.end) + "  zone:" + str(self.zone))
    
class Zone:
    def __init__(self, begin, end, delta):
        self.begin = begin
        self.end = end
        self.delta = delta
    def Display(self):
        print ("[Zone]   begin:" + str(self.begin) + "  end:" + str(self.end) + "  delta:" + str(self.delta))

def generar_grafico(nombre_archivo):

    #with open(nombre_archivo) as f:
    #    values = json.load(f)
    
    #f.close()
    m = 30
    m2 = 20
    w = 2100
    s = 0
    values = json.loads(nombre_archivo)
    exons = []
    amplicons  = []
    zones  = []

    l = len(values)


    # https://www.programiz.com/python-programming/json
    #https://realpython.com/python-json/
    # obtiene nombre de gen y de archivo de imagen desde el JSON
    # https://linuxhint.com/search_json_python/
    #valores por defecto
    genname = "gen"
    imagefile = "salida"
    # recore JSON para buscar nombres de archivo, gen y crear arreglos de exones y amplicones
    numexons = 0
    numamplicons = 0
    numzones = 0
    numzonesaux = 0


    zonebeginaux = 2147483647
    zoneendaux = -2147483648

    for keyval in values:

        if ( keyval["type"] == "metadata" ):
            genname = keyval['genname']
            imagefile = keyval['imagefile']

        if ( keyval["type"] == "amplicon" ):
            numamplicons +=1
            amplicons.append(Amplicon(keyval["name"],keyval["begin"],keyval["end"],keyval["zone"]))

        if ( keyval["type"] == "exon" ):
            numexons +=1
            exons.append(Exon(keyval["name"],keyval["begin"],keyval["end"],keyval["zone"]))

        if ( keyval["type"] == "amplicon" or keyval["type"] == "exon" ) :
            numzonesaux = keyval["zone"]
            
            if ( numzonesaux > numzones ):
                numzones = numzonesaux
    
    numzones += 1
    
    for i in range(0,numzones):
        zones.append(Zone(2147483647,-2147483647,0))

    for i in range(numzones):

        for keyval in values:

            if ( keyval["type"] == "amplicon" or keyval["type"] == "exon" ) :
                
                if ( int(keyval["zone"]) == i ):

                    if( int(keyval["begin"]) < zones[i].begin ):
                        zones[i].begin = int(keyval["begin"])

                    if( int(keyval["end"]) > zones[i].end ):
                        zones[i].end = int(keyval["end"])
        

    print("exones: "+str(numexons)+" amplicones: "+str(numamplicons)+" zonas: "+str(numzones))
    for i in range(numzones):
        print("begin: "+str(zones[i].begin)+" end: "+str(zones[i].end))

    for i in range(numzones-1):

        if( zones[i].end > zones[i+1].begin ):
            zones[i+1].begin = zones[i].begin
            zones[i].end = zones[i+1].end

    for i in range(numzones-1,0,-1):

        if( zones[i-1].begin == zones[i].begin ):
            zones[i-1].end = zones[i].end

    zones[0].delta = m - zones[0].begin
    for i in range(1,numzones):
        
        if( zones[i].begin == zones[i-1].begin ):
            zones[i].delta = zones[i-1].delta
        else:
            zones[i].delta = 2 * m + zones[i-1].end + zones[i-1].delta - zones[i].begin

    #valores cualqueira de tamaño de imagen
    width = 500 + numamplicons*60
    height = 300
    
    showDNAlen = zones[numzones-1].end + zones[numzones-1].delta + m
    s = width/showDNAlen
    #define nombres de colores para no usar valores RGB más adelante
    white = (255, 255, 255)
    black = (0,0,0)
    blue = (0,0,255)
    red_gemini = (190, 39, 47)
    red = (255,0,0)
    gray = (200,200,200)
    gray2 = (85, 85, 85)
    gray_breaks = (64, 64, 64)

    print("nombre gen: "+genname+" imagefile: "+imagefile)
    #referencia para dibujo: https://code-maven.com/create-images-with-python-pil-pillow
    # crea imagen para archivo png
    imagen = Image.new("RGB", (width+2*m2, height), white)
    draw = ImageDraw.Draw(imagen)

    hEx = 150
    hAm = 45

    sw = 0
    sw2 = 0
    d = 40
    d2 = 30
    z = 0
    p1 = 0
    p2 = 0
    center = 0
    centroActual = 0
    centroProx = 0
    fnt = ImageFont.truetype('OrkneyRegular.ttf', 11)

    for i in range(numamplicons):
        z = amplicons[i].zone
        p1 = amplicons[i].begin + zones[z].delta
        p2 = amplicons[i].end + zones[z].delta
        
        #lineas grises??
        draw.line([s*p1+m2, hAm+sw*d, s*p1+m2, hEx],gray)
        draw.line([s*p2+m2, hAm+sw*d, s*p2+m2, hEx],gray)

        #linea central
        draw.line([s*p1+m2, hAm + sw*d, s*p2+m2, hAm + sw*d],red_gemini)

        #topes de extremo
        draw.line([s*p1+m2, hAm + sw*d , s*p1+m2, hAm + sw*d + 6],red_gemini)
        draw.line([s*p2+m2, hAm + sw*d , s*p2+m2, hAm + sw*d + 6],red_gemini)

        #texto
        center = int( ((amplicons[i].begin+amplicons[i].end)/2+zones[z].delta) )
        draw.text([s*center+m2,hAm + sw*d + -10],str(amplicons[i].begin) + "-" + str(amplicons[i].end),fill=red_gemini,font=fnt,anchor="mm")
        draw.text([s*center+m2,hAm + sw*d - 25],amplicons[i].name+" ("+str(amplicons[i].end-amplicons[i].begin +1)+")",fill=red_gemini,font=fnt,anchor="mm")

        sw+=1
        if( sw == 3 ):
            sw = 0

    #dibuja linea DNA
    draw.line([0,hEx,width+2*m2,hEx],black)

    #dibuja exones
    x1 = 0
    y1 = 0
    x2 = 0
    y2 = 0
    for i in range(numexons):
        z = exons[i].zone

        x1 = s*(exons[i].begin+zones[z].delta)+m2
        y1 = (hEx-8)
        x2 = x1 + s*(exons[i].end-exons[i].begin)
        y2 = y1 + 16
        draw.rectangle([(x1,y1),(x2,y2)],fill=gray)
        center=int( ((exons[i].begin+exons[i].end)/2+zones[z].delta) )

        draw.text([s*center+m2,hEx+25 + d2*sw2],exons[i].name+" ("+str(int(exons[i].end-exons[i].begin))+")",fill=gray2,font=fnt,anchor="mm")

        draw.text([s*center+m2, hEx+40 + d2*sw2],str(int(exons[i].begin)) + "-" + str(int(exons[i].end)),fill=gray2,font=fnt,anchor="mm")

        if( sw2 > 0 ):
            print("linea para textos lejanos?: "+str(sw2))
            draw.line([s*center+m2, hEx+12, s*center+m2, hEx + d2*sw2+15],fill=gray2)

        if( i < numexons-1 ):
            z2 = exons[i+1].zone
            centroActual = s*((exons[i].end   + exons[i].begin)/2   + zones[z].delta)
            centroProx   = s*((exons[i+1].end + exons[i+1].begin)/2 + zones[z2].delta)

            if(centroProx - centroActual < 100):
                sw2 += 1
            else:
                sw2 = 0

        if( sw2 == 4):
            sw2 = 0    

    #dibuja breaks
    shl = 5
    incl = 1
    cutPos = 0
    sep = 0
    x1 = 0
    y1 = 0
    x2 = 0
    y2 = 0
    for i in range(numzones-1):
        cutPos = zones[i].end + zones[i].delta + m
        sep = 2
        x1 = s*(cutPos-sep)
        y1 = hEx-2+m2

        x2 = x1+2*s*sep
        y2 = y1 + 4
        draw.rectangle([(x1,y1),(x2,y2)],fill=white)

        draw.line([s*(cutPos-sep)+incl+m2, hEx-shl, s*cutPos - s * sep-incl+m2, hEx+shl],gray_breaks)
        draw.line([s*(cutPos+sep)+incl+m2, hEx-shl, s*cutPos + s * sep-incl+m2, hEx+shl],gray_breaks)

    fnt_nombregen = ImageFont.truetype('OrkneyRegular.ttf', 20)
    draw.text([10,275],genname,fill=black,font=fnt_nombregen)
    #draw.line([10,10,50,50], red)
    #draw.text((30,30), genname, fill=black)
    
    #draw.text((30,110), str(numamplicons)+" amplicones", fill=red, font=fnt)
    #draw.rectangle((70,70,100,100), fill=gray, outline=blue)
    #graba imagen al archivo especificado en el JSON
    #***********************************************************************************
    imagen.save(imagefile+".png")
    #imagen.save("/var/www/html/imagenes/"+imagefile+".png")

def generatePrimers(gene,mim=15,mex=100,pmil=20,pmal=27,magcp=0.6,migcp=0.4,tmmi=0,dtmmi=0,dtmma=4,lam=450,
        translape=25,maxhomopol=4,margenbusqueda=300,dg=-6000,interval=50,mvconc=50,fotm=-1,fod=5, lamin=-1):
    global var
    var = {
    'gene_filename' : gene,
    'exones_seleccionados' : [],
    'margin' : mim,
    'margin_between_exons' : mex,
    'primer_min_length' : pmil,
    'primer_max_length' : pmal,
    'max_gc_percent' : magcp,
    'min_gc_percent' : migcp,
    'tm_min' : tmmi,
    'mv_conc' : mvconc,
    'focus_tm': fotm,
    'focus_d': fod,
    'delta_tm_min' : dtmmi,
    'delta_tm_max' : dtmma,
    'largo_amplicon': lam,
    'traslape' : translape,
    'max_homopol' : maxhomopol,
    'margen_busqueda' : margenbusqueda,
    'deltag' : dg,
    'intervalo' : interval,
    'aplicar_intervalo' : False,
    'overhangF' : "",
    'overhangR' : ""
    }


    load_gene(var['gene_filename'])
    
    ####################################SECCION EDITADA#####################################
    #exones_a_ejecutar = [0,1,2,3,4,5,6,7,8,9,10,11,12]
    exones_a_ejecutar = [15]
    ejecutar_todo = True
    aux_exon_position = []

    if(ejecutar_todo == False):  

        for i in range(len(exones_a_ejecutar)):
            aux_exon_position.append(exon_position[exones_a_ejecutar[i]]) 
        
        largo_exon_pos = len(exon_position)
        exon_position.clear()

        for k in range(len(aux_exon_position)):
            exon_position.append(aux_exon_position[k])

    areas_to_extract = areas_of_interest_set(var['margin'], var['margin_between_exons'])

    ####################################SECCION EDITADA#####################################
    print(areas_to_extract)
    print(exon_position)
    #test_list = [(1,3),(4,5),(5,2),(6,2),(2,7),(3,5),(5,5),(4,4)]
    #test_list.sort(key=key_to_sort)
    #print(test_list)
    #first_primer_set = primers_search(var['primer_min_length'],var['primer_max_length'],var['min_gc_percent'],var['max_gc_percent'],areas_to_extract)
    


    first_primer_set = primers_search(var['primer_min_length'],var['primer_max_length'],0,100,areas_to_extract)
    rel = Relajo()
    #promedio_tamaño =( (areas_to_extract[i][1] - areas_to_extract[i][0]) / (var['largo_amplicon'] - (var['margen_busqueda'] / 2) ) )
    exoma_suma = 0
    for i in range(len(areas_to_extract)):
        if( (areas_to_extract[i][1] - areas_to_extract[i][0]) > (var['largo_amplicon'] - ( var['margen_busqueda']/2 ) ) ):
            exoma_suma = exoma_suma + (areas_to_extract[i][1] - areas_to_extract[i][0]) + ( ( var['traslape'] +  var['margen_busqueda'] ) * ( (areas_to_extract[i][1] - areas_to_extract[i][0]) / (var['largo_amplicon'] - (var['margen_busqueda'] / 2) ) )  )
        else:
            exoma_suma = exoma_suma + ( areas_to_extract[i][1] - areas_to_extract[i][0] )+ var['margen_busqueda']
        #exoma_suma = exoma_suma + ( areas_to_extract[i][1] - areas_to_extract[i][0] )
    exoma_suma = exoma_suma * 1.7
    print("largo exoma: "+str(exoma_suma))
    if lamin==-1:
        largo_min_amplicon = var['largo_amplicon'] - var['margen_busqueda']
    else:
        largo_min_amplicon = lamin
    recomendador_largo = RecomendadorDeLargo(exoma_suma, len(exon_position), var['margin'], largo_min_amplicon, var['largo_amplicon'])
    ##print(first_primer_set)
    primers_seleccionados = seleccion_partidores(areas_to_extract,rel,recomendador_largo)
    print(areas_to_extract)
    print (primers_seleccionados)

    ####CORRECCION PARA GENERAR CSV DE FORMA LOCAL
    aux_dict = {}
    for i in range(len(exon_position)):
        aux_dict = {}
        aux_dict['label'] = str(i)
        var['exones_seleccionados'].append(aux_dict)

    csvname=var['gene_filename'].split(".")[0]
    #generar_csv(primers_seleccionados,areas_to_extract,"CSV_Primers",rel)
    generar_csv(primers_seleccionados,areas_to_extract,csvname,rel)
    print(exon_position)
    
    json_exones = {}
    json_exones_txt = open("JSON_exones_y_amplicones.txt","w")
    nombre_exon = ""
    numero_exon = 0
    json_exones_final = ""
    for i in range(len(exon_position)):
        for k in range(len(areas_to_extract)):
            ###print("Area a extraer 0: "+str(areas_to_extract[k][0])+" Area a extraer 1: "+str(areas_to_extract[k][1])+" Exon pos 0: "+str(exon_position[i][0])+" Exon pos 1: "+str(exon_position[i][1])+"")
            if ( (exon_position[i][0] > areas_to_extract[k][0]) and (exon_position[i][1] < areas_to_extract[k][1]) ):
                json_exones['type'] = "exon"
                numero_exon = i+1
                nombre_exon = "Exon "+str(numero_exon)
                json_exones['name'] = nombre_exon
                json_exones['begin'] = exon_position[i][0]
                json_exones['end'] = exon_position[i][1]
                json_exones['zone'] = k
                print(json_exones)
                if (numero_exon != 1):
                    json_exones_final = json_exones_final + "," + json.dumps(json_exones)
                else:
                    json_exones_final = json.dumps(json_exones)

    json_amplicones_txt = open("JSON_amplicones.txt","r")
    json_amplicones = json_amplicones_txt.readline()
    fecha = datetime.today().strftime('%d-%m-%Y-%H:%M:%S')
    fecha = fecha.replace(":","_")
    nombre_imagen = var['gene_filename'][:-5]+ "_" + fecha
    metadata = {}
    metadata['type'] = 'metadata'
    metadata['genname'] = var['gene_filename'][:-5]
    metadata['imagefile'] = nombre_imagen
    json_exones_final = "[" + json_amplicones + "," + json_exones_final + "," + json.dumps(metadata) + "]"
    #exones_object = json.dumps(json_exones_final)
    #json_exones_final = "["+json_exones_final+"]"
    json_exones_txt.write(json_exones_final)
    json_exones_txt.close()
    print("json_exones_final")
    print(json_exones_final)
    generar_grafico(json_exones_final)        
    #blast_string = ""
    #blastn_cline = NcbiblastnCommandline(query = "test.txt", num_descriptions=5, num_alignments=5, num_threads=4, out="resultadosXML.xml", word_size=11, outfmt=5, db = "PKD1")
    #print(blastn_cline)
    

    #primers_verificados = verificacion_unicidad(primers_seleccionados,areas_to_extract,rel)
    
    #print(areas_to_extract)
    #result_handle = NCBIWWW.qblast("blastn", "nt", Seq("ggctgagatgtttatctttctccat", generic_dna))
    #result_handle = NCBIWWW.qblast(program="blastn", database="nt", sequence= Seq("ggctgagatgtttatctttactccat", generic_dna))
    #result_handle = NCBIWWW.qblast(program='blastn', database='refseq_genomic',sequence=(Seq("agggtgggaggccatagcttggggatgctg", generic_dna)))

    
    #blast_records = NCBIXML.parse(result_handle)

    #for b in blast_records:
    #    for alignment in b.alignments:
    #        for hsp in alignment.hsps:
    #            print("***Alignment***")
    #            print("Sequence: ", alignment.title)
    #            print("length: ", alignment.length)
    #            print("e value: ", hsp.expect)
    #            print(hsp.query)
    #            print(hsp.match)
    #            print(hsp.sbjct)

    print("Se generó archivo <{}> en directorio.".format(csvname))

if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-gfn', '--genfilename', help='Nombre del archivo de secuencia del gen.', required=True, default=argparse.SUPPRESS)

    parser.add_argument('-marint', '--margenintron',help='Margen intronico.' ,type=int,default=15)
    parser.add_argument('-marex', '--margenexon',help='Margen entre exones.' ,type=int,default=100)
    parser.add_argument('-mb', '--marbus',help='Margen de busqueda.' ,type=int,default=300)
    
    parser.add_argument('-mvc', '--mvconc',help='Concentración de iones monovalentes en medio de reacción.',type=int ,default=50)
    parser.add_argument('-dgd', '--dgdesign',help='Minimo deltaG para diseño.' ,type=int,default=-6000)
    
    parser.add_argument('-pmax', '--primermaxlen',help='Largo maximo de primer.' ,type=int,default=27)
    parser.add_argument('-pmin', '--primerminlen',help='Largo minimo de primer.' ,type=int,default=20)
    
    parser.add_argument('-maxgcp', '--maxgcpercent',help='Maximo porcentaje GC.' ,type=float,default=0.6)
    parser.add_argument('-mingcp', '--mingcpercent',help='Minimo porcentaje GC.' ,type=float,default=0.4)
    
    parser.add_argument('-tmmi', '--mintemp',help='Tm minima.',type=int ,default=0)
    parser.add_argument('-dtmmi', '--mindeltatemp',help='delta Tm minima.' ,type=int,default=0)
    parser.add_argument('-dtmma', '--maxdeltatemp',help='delta Tm maxima.' ,type=int,default=4)

    parser.add_argument('-lamp', '--largoamplicon',help='Largo amplicon.' ,type=int,default=450)
    parser.add_argument('-lmin', '--largominamp',help='Largo mínimo amplicon.' ,type=int,default=-1)

    parser.add_argument('-tlp', '--traslape',help='Traslape amplicones minimos.' ,type=int,default=25)
    parser.add_argument('-mhp', '--maxhomopol',help='Maximo de homopolimeros.' ,type=int,default=4)
    
    parser.add_argument('-fotm', '--focustm',help='Temperatura objetivo para partidores.' ,type=int,default=-1)
    parser.add_argument('-fod', '--focusdeltatm',help='Maximo delta para temperatura objetivo para partidores.' ,type=int,default=5)

    
    args = parser.parse_args()
    dargs = vars(args)
    #dargs[""]
    
    generatePrimers(gene=dargs["genfilename"],mim=dargs["margenintron"],mex=dargs["margenexon"],pmil=dargs["primerminlen"], pmal=dargs["primermaxlen"],magcp=dargs["maxgcpercent"],migcp=dargs["mingcpercent"],tmmi=dargs["mintemp"],dtmmi=dargs["mindeltatemp"], dtmma=dargs["maxdeltatemp"],lam=dargs["largoamplicon"],translape=dargs["traslape"],maxhomopol=dargs["maxhomopol"], margenbusqueda=dargs["marbus"],dg=dargs["dgdesign"],mvconc=dargs["mvconc"],fotm=dargs["focustm"],fod=dargs["focusdeltatm"], lamin=dargs["largominamp"])
