#!/usr/bin/env python
# coding: utf-8

# In[1]:


import json
from multiplexador import multiplexPrimers,write2cvs
from disenador import generatePrimers
import sys
import importlib
import shutil


# In[2]:


def checkExons(jsonfile):
    """!
    Función para revisar si paneles creados por diseñador cubren todos los amplicones.
    \param <jsonfile> {String, nombre del json generado por diseñador.}
    """ 
    ampz=[]
    exz=[]
    with open(jsonfile) as json_file:
        data = json.load(json_file)
        #print(data)
    for d in data:
        if d["type"]=="exon":
            n=d["zone"]
            if n not in exz:
                exz.append(n)
        if d["type"]=="amplicon":
            n=d["zone"]
            if n not in ampz:
                ampz.append(n)
    print(ampz)
    print(exz)
    for x in exz:
        if x not in ampz:
            raise TypeError("No cubre todos los exones.") 
            return False
    return True


# In[3]:


def bestColor(dlist,fname):
    """!
    Función que encuentra panel con menos grupos y lo escribe a archivo.
    \param <dlist> {Lista, paneles de partidores para distintas temperaturas objetivo.}
    \param <fname> {String, nombre del archivo a escribir.}
    """ 
    colorn=9999
    for d in dlist:
        if d["cnumber"]<colorn:
            colorn=d["cnumber"]
            bf=d["ft"]
            bb=d["balance"]
            bc=d["colors"]
    print("Mejor focus es "+str(bf))
    fname2=fname.replace(".csv","({}).csv".format(str(bf)))
    write2cvs(bc,bb,fname2)
    


# In[4]:


#[-1,64,66,68,70]
def controllerProgram(seqfile,multi=False,mvc=50,tmmin=0,marex=100,lamp=400,larmin=-1,marbus=300 ,deltag=-6000,multidg=-13627,ftemps=[-1,64,66,68,70]):
    """!
    Función para realizar el diseño y multiplexación de primers segun parámetros dados.
    \param <seqfile> {String, nombre del archivo de secuencia del gen.}
    \param <multi> {Booleano, si corresponde a un archivo multigen.}
    \param <mvc> {Numero, concentración de iones monovalentes.}
    \param <tmmin> {Numero, temperatura mínima a considerar.}
    \param <marex> {Entero, margen entre exones.}
    \param <lamp> {Entero, largo amplicon.}
    \param <larmin> {Entero, largo mínimo amplicon.}
    \param <marbus> {Entero, margen de busqueda.}
    \param <deltag> {Numero, minima deltaG para homodimeros y horquillas (diseñador).}
    \param <multidg> {Numero, minima deltaG para heterodimeros (multiplexador).}
    \param <ftemps> {Lista, temperaturas a probar para diseño y multiplexación.}
    """ 
    results=[]
    for ft in ftemps:
        try:
            if multi:
                csvs=[]
                for f in seqfile:
                    importlib.reload(sys.modules['disenador'])
                    from disenador import generatePrimers
                    generatePrimers(f,fotm=ft,mvconc=mvc,mex=marex,tmmi=tmmin,lam=lamp,lamin=larmin,
                                    margenbusqueda=marbus,dg=deltag)
                    checkExons("JSON_exones_y_amplicones.txt")
                    name=f.split(".")[0]+".csv"
                    shutil.copyfile(name, "f({})_".format(ft)+name)
                    csvs.append(name)
                    csvname=csvs
            else:
                importlib.reload(sys.modules['disenador'])
                from disenador import generatePrimers
                generatePrimers(seqfile,fotm=ft,mvconc=mvc,mex=marex,tmmi=tmmin,lam=lamp,lamin=larmin,
                                    margenbusqueda=marbus,dg=deltag)
                checkExons("JSON_exones_y_amplicones.txt")
                csvname=seqfile.split(".")[0]+".csv"
                shutil.copyfile(csvname,"f({})_".format(ft)+csvname)
        except Exception as e:
            print(e)
            if ft==-1:
                print("NO SE PUDO ENCONTRAR SOLUCIONES CON PARAMETROS ESPECIFICADOS")
            else:
                continue
            
        #csvname=seqfile.split(".")[0]+".csv"
        colors,bald=multiplexPrimers(seqfile,"word",csvname,"GEMIN1",mv_conc=mvc,mindg=multidg,isMultiple=multi)
        #print("Se tienen "+str(len(colors)))
        td={"ft":ft,"colors":colors,"cnumber":len(colors),"balance":bald}
        results.append(td)
    
    try:
        if multi:
            if len(results)==0:
                raise TypeError("No hay soluciones.")
            else:
                csvs=[]
                for f in seqfile:
                    name=f.split(".")[0]+".csv"
                    csvs.append(name)
                mcsvname="_".join(f.replace(".csv","") for f in csvs)+".csv"
        else:
            mcsvname="multiplexed_"+csvname
            
        bestColor(results,mcsvname)
    except Exception as e:
        print(e)
        print("NO SE PUDO ENCONTRAR SOLUCIONES CON PARAMETROS ESPECIFICADOS")
        exit()


if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-gfn', '--genfilename', action='append', help='Nombre del archivo de secuencia del gen (Repetir para multiple).', required=True, default=argparse.SUPPRESS)
    parser.add_argument('-tm', '--tmin',help='Temperatura mínima a considerar.',type=int,default=0)
    parser.add_argument('-marex', '--margenexon',help='Margen entre exones.',type=int,default=100)
    parser.add_argument('-lamp', '--largoamplicon',help='Largo amplicon.' ,type=int,default=400)
    parser.add_argument('-lmin', '--larminamp',help='Largo mínimo amplicon.' ,type=int,default=-1)
    parser.add_argument('-mb', '--marbus',help='Margen de busqueda.' ,type=int,default=300)
    parser.add_argument('-dgd', '--dgdesign',help='Minimo deltaG para diseño.' ,type=int,default=-9000)
    parser.add_argument('-dgm', '--dgmulti',help='Minimo deltaG para multiplexar.' ,type=int,default=-7000)
    parser.add_argument('-mvc', '--mvconc',help='Concentración de iones monovalentes en medio de reacción.' ,type=float,default=50)

    parser.add_argument('-tmps', '--temperatures', action='append', help='Temperaturas a probar (Repetir para multiples). (default: [-1,64,66,68,70])',type=int, default=argparse.SUPPRESS)

    args = parser.parse_args()
    dargs = vars(args)
    
    #tm2try=dargs["temperatures"]
    #if len(tm2try)==0:
    #    tm2try=[-1,64,66,68,70]
    if "temperatures" in dargs.keys():
        tm2try=dargs["temperatures"]
    else:
        tm2try=[-1,64,66,68,70]

    #print(tm2try)

    argseq=dargs["genfilename"]
    multivar=len(dargs["genfilename"])>1
    if not(multivar):
        argseq=argseq[0]
    
    controllerProgram(argseq,tmmin=dargs["tmin"],marex=dargs["margenexon"],lamp=dargs["largoamplicon"],marbus=dargs["marbus"], deltag=dargs["dgdesign"],multidg=dargs["dgmulti"],larmin=dargs["larminamp"],multi=multivar,mvc=dargs["mvconc"], ftemps=tm2try)




