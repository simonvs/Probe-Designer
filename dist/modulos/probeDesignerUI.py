#import modulos.descarga  as descarga
#import modulos.diseno as diseno
#import modulos.imagen as imagen
import tkinter as tk
import pandas as pd
import customtkinter as ctk
import os
import sys
import sqlite3
import datetime
import platform
import shutil
import threading
import queue
import time
import openpyxl
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from PIL import Image, ImageTk
import importlib.util

print(sys.version)

def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS2
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

spec1 = importlib.util.spec_from_file_location("descarga", resource_path(os.path.join("modulos", "descarga.py")))
spec2 = importlib.util.spec_from_file_location("diseno", resource_path(os.path.join("modulos", "diseno.py")))
spec3 = importlib.util.spec_from_file_location("imagen", resource_path(os.path.join("modulos", "imagen.py")))

descarga = importlib.util.module_from_spec(spec1)
diseno = importlib.util.module_from_spec(spec2)
imagen = importlib.util.module_from_spec(spec3)
sys.modules["descarga"] = descarga
sys.modules["diseno"] = descarga
sys.modules["imagen"] = descarga

spec1.loader.exec_module(descarga)
spec2.loader.exec_module(diseno)
spec3.loader.exec_module(imagen)

class PantallaInicial:
    def __init__(self, root):
        self.root = root
        self.root.wm_title("GEMINi - Diseñador de sondas de hibridación para ARN mensajero")
        self.root.wm_iconbitmap(resource_path(os.path.join("images", "gemini2.ico")))
        #self.frame = tk.Frame(root)
        self.frame = ctk.CTkFrame(master=root)
        self.frame.place(in_=root, relx=0.5, rely=0.5, anchor='c')
        for widget in self.frame.winfo_children():
                widget.destroy()
        self.root.state("zoomed")
        #self.frame.pack(padx=100, pady=200)
        self.crear_interfaz()
        #self.root.minsize(800,450)
        #ancho_pantalla = self.root.winfo_screenwidth()
        #alto_pantalla = self.root.winfo_screenheight()
        #self.root.minsize(ancho_pantalla,alto_pantalla)
        #print(ancho_pantalla, alto_pantalla)
        
        #self.root.geometry("800x450")
        #self.root.attributes("-alpha", True)
        #self.root.geometry(f"{ancho_pantalla}x{alto_pantalla}")
    

    def crear_interfaz(self):

        if not os.path.exists(resource_path('databases')):
            os.makedirs(resource_path('databases'))
        if not os.path.exists(resource_path('sondas')):
            os.makedirs(resource_path('sondas'))

        # if not os.path.exists(os.path.join(os.getcwd(),'databases')):
        #     os.makedirs(os.path.join(os.getcwd(),'databases'))
        # if not os.path.exists(os.path.join(os.getcwd(),'sondas')):
        #     os.makedirs(os.path.join(os.getcwd(),'sondas'))

        # Crear widgets de la pantalla de inicio
        #fuente_titulo = font.Font(weight="bold", size=16)
        #titulo_archivo = tk.Label(self.frame, text="Seleccione un archivo GenBank", font=fuente_titulo)
        fuente_titulo = ctk.CTkFont(family='Helvetica', size=20, weight='bold')
        titulo_archivo = ctk.CTkLabel(self.frame, text="Seleccione un archivo GenBank (.gb o .gbk)\no introduzca su accession number", text_color="#000000", font=fuente_titulo)
        titulo_archivo.grid(row=0, column=0, columnspan=3, padx=100, pady=20)

        # Pantalla de selección de archivo
        #seleccion_button = tk.Button(self.frame, text="Seleccionar Archivo", command=self.seleccionar_archivo)
        seleccion_button = ctk.CTkButton(self.frame, text="Seleccionar Archivo", corner_radius=30, fg_color="#404040", command=self.seleccionar_archivo)
        seleccion_button.grid(row=1, column=0, rowspan=2, padx=15, pady=30)

        obien_label = ctk.CTkLabel(self.frame, text="o bien...", text_color="#000000")
        obien_label.grid(row=1, column=1, rowspan=2, pady=30)

        accession_entry = ctk.CTkEntry(self.frame, placeholder_text="Accession Number", corner_radius=30)
        accession_entry.grid(row=1, column=2, padx=15, pady=30)


        def descargar_secuencia():
            accnum = accession_entry.get()
            if len(accnum) > 5:
                try:                    
                    self.seqrecord = descarga.accnum_to_seqrecord(accnum)
                    self.continuar_button.configure(state='normal')
                    self.file_label.configure(text='Archivo seleccionado: '+accnum+".gbk")
                    #self.archivo = os.path.join(os.getcwd(), 'files', accnum+".gbk")
                    self.archivo = resource_path(os.path.join('files', accnum+".gbk"))
                    #app.mostrar_pantalla("transcripciones", seqrecord, filepath=self.archivo)
                except:
                    messagebox.showinfo("Error", "Error al descargar la secuencia")

        descargar_button = ctk.CTkButton(self.frame, text="Descargar Secuencia", corner_radius=30, fg_color="#404040", command=descargar_secuencia)
        descargar_button.grid(row=2, column=2, pady=10)

        self.frame.grid_rowconfigure(3, minsize=25)

        #historial_button = tk.Button(self.frame, text="Ver historial de sondas", command=self.ver_historial)
        historial_button = ctk.CTkButton(self.frame, text="Ver historial de sondas", corner_radius=30, fg_color="#7a7a7a",command=self.ver_historial)
        historial_button.grid(row=4, column=0, rowspan=2, pady=25)

        def continuar():
            if self.seqrecord:
                app.mostrar_pantalla("transcripciones", self.seqrecord, filepath=self.archivo)

        self.file_label = ctk.CTkLabel(self.frame, text="", text_color="#000000")
        self.file_label.grid(row=4, column=2, pady=5)

        self.continuar_button = ctk.CTkButton(self.frame, text="Continuar", corner_radius=30, fg_color="#404040", state='disabled', command=continuar)
        self.continuar_button.grid(row=5, column=2, pady=10)
        #self.root.geometry("800x450")

        
    def seleccionar_archivo(self):
        self.archivo = filedialog.askopenfilename(filetypes=[("Archivos GenBank", "*.gbk *.gb")])
        if self.archivo:            
            self.file = os.path.split(self.archivo)[1]
            self.continuar_button.configure(state='normal')
            self.file_label.configure(text='Archivo seleccionado: '+self.file)
            self.seqrecord = descarga.parse_file_to_seqrecord(self.archivo)
            #app.mostrar_pantalla("transcripciones", self.seqrecord, filepath=archivo)

    
    def ver_historial(self):
        app.mostrar_pantalla('historial')


class PantallaTranscripciones(PantallaInicial):
    def __init__(self, root, seqrecord, filepath):
        self.seqrecord = seqrecord
        self.filepath = filepath
        super().__init__(root)

    def crear_interfaz(self):
        # Crear widgets de la pantalla de inicio
        #transcripciones_label = tk.Label(self.frame, text=f"Seleccione las transcripciones para el diseño de sondas", font=fuente_titulo)
        fuente_titulo = ctk.CTkFont(family='Helvetica', size=20, weight='bold')
        transcripciones_label = ctk.CTkLabel(self.frame, text="Seleccione los transcritos para el diseño de sondas", text_color="#000000", font=fuente_titulo)
        transcripciones_label.grid(row=0, column=0, columnspan=3, padx=10, pady=20)

        label_izq = ctk.CTkLabel(self.frame, text="Transcritos no seleccionados", text_color="#000000")
        label_izq.grid(row=1, column=0, padx=10, pady=10)

        label_der = ctk.CTkLabel(self.frame, text="Transcritos seleccionados", text_color="#000000")
        label_der.grid(row=1, column=2, padx=10, pady=10)

        #Obtener lista con transcripciones
        
        self.opciones_disponibles = []
        self.opciones_seleccionadas = []

        for f in self.seqrecord.features:
            if f.type == 'CDS':
                info = f.qualifiers.get('product')[0]
                if str(info) not in self.opciones_disponibles:
                    self.opciones_disponibles.append(str(info))
        
        def maxlen(a):
            maximo = 0
            for i in a:
                length = len(i)
                if length > maximo:
                    maximo = length
            return maximo
        
        self.lista_disponibles = tk.Listbox(self.frame, selectmode=tk.MULTIPLE, selectbackground="#404040", width=maxlen(self.opciones_disponibles))
        self.lista_disponibles.grid(row=2, column= 0, padx=10, pady=10)

        for opcion in self.opciones_disponibles:
            self.lista_disponibles.insert(tk.END, opcion)

        self.lista_seleccionadas = tk.Listbox(self.frame, selectmode=tk.MULTIPLE, selectbackground="#404040", width=maxlen(self.opciones_disponibles))
        self.lista_seleccionadas.grid(row=2, column= 2, padx=10, pady=10)

        botones_frame = ctk.CTkFrame(master=self.frame, bg_color="transparent")
        botones_frame.grid(row=2, column=1)

        #self.boton_seleccionar = tk.Button(botones_frame, text="Seleccionar", command=self.seleccionar)
        self.boton_seleccionar = ctk.CTkButton(botones_frame, text=">", corner_radius=30, fg_color="#7a7a7a", command=self.seleccionar)
        self.boton_seleccionar.pack(pady=10)
        
        #self.boton_eliminar = tk.Button(botones_frame, text="Eliminar", command=self.eliminar)
        self.boton_eliminar = ctk.CTkButton(botones_frame, text="<", corner_radius=30, fg_color="#7a7a7a", command=self.eliminar)
        self.boton_eliminar.pack(pady=10)

        #self.boton_seleccionar_todo = tk.Button(botones_frame, text="Seleccionar Todo", command=self.seleccionar_todo)
        self.boton_seleccionar_todo = ctk.CTkButton(botones_frame, text=">>>>", corner_radius=30, fg_color="#404040", command=self.seleccionar_todo)
        self.boton_seleccionar_todo.pack(pady=10)

        #self.boton_eliminar_todo = tk.Button(botones_frame, text="Eliminar Todo", command=self.eliminar_todo)
        self.boton_eliminar_todo = ctk.CTkButton(botones_frame, text="<<<<", corner_radius=30, fg_color="#404040", command=self.eliminar_todo)
        self.boton_eliminar_todo.pack(pady=10)

        def obtener_seleccion():
            seleccionados = self.lista_seleccionadas.get(0, tk.END)
            #pantalla parametreos
            if len(seleccionados) >= 1:
                app.mostrar_pantalla("parametros", self.seqrecord, [seleccionado for seleccionado in seleccionados], filepath=self.filepath)
            else:
                messagebox.showinfo("Error", "Se debe seleccionar al menos una transcripción para el diseño de las sondas")

        # Botón para obtener las opciones seleccionadas
        #boton_volver = tk.Button(self.frame, text="Volver", command=self.volver_a_seleccion)
        boton_volver = ctk.CTkButton(self.frame, text="Volver", corner_radius=30, fg_color="#7a7a7a", command=self.volver_a_seleccion)
        boton_volver.grid(row=3, column=0, pady=20)

        # Botón para obtener las opciones seleccionadas
        #boton_obtener_seleccionadas = tk.Button(self.frame, text="Continuar", command=obtener_seleccion)
        boton_obtener_seleccionadas = ctk.CTkButton(self.frame, text="Continuar", corner_radius=30, fg_color="#404040", command=obtener_seleccion)
        boton_obtener_seleccionadas.grid(row=3, column=2, pady=20)

        self.root.geometry("800x500")

    def seleccionar(self):
        #seleccion = self.lista_disponibles.get(tk.ACTIVE)
        selecciones = self.lista_disponibles.curselection()
        selecciones = sorted(selecciones, reverse=True)
        print(selecciones)
        for indice in selecciones:
            seleccion = self.lista_disponibles.get(indice)
            if seleccion not in self.opciones_seleccionadas:
                self.opciones_seleccionadas.append(seleccion)
                self.lista_seleccionadas.insert(tk.END, seleccion)
                self.lista_disponibles.delete(indice)
                self.opciones_disponibles.remove(seleccion)

    def eliminar(self):
        #seleccion = self.lista_seleccionadas.get(tk.ACTIVE)
        selecciones = self.lista_seleccionadas.curselection()
        selecciones = sorted(selecciones, reverse=True)
        print(selecciones)
        for indice in selecciones:
            seleccion = self.lista_seleccionadas.get(indice)
            if seleccion in self.opciones_seleccionadas:
                self.opciones_disponibles.append(seleccion)
                self.lista_disponibles.insert(tk.END, seleccion)
                self.lista_seleccionadas.delete(indice)
                self.opciones_seleccionadas.remove(seleccion)

    def seleccionar_todo(self):
        for opcion in self.opciones_disponibles:
            if opcion not in self.opciones_seleccionadas:
                self.opciones_seleccionadas.append(opcion)
                self.lista_seleccionadas.insert(tk.END, opcion)
        self.lista_disponibles.delete(0, tk.END)
        self.opciones_disponibles.clear()

    def eliminar_todo(self):
        for opcion in self.opciones_seleccionadas:            
            self.lista_disponibles.insert(tk.END, opcion)
            self.opciones_disponibles.append(opcion)
        self.lista_seleccionadas.delete(0, tk.END)
        self.opciones_seleccionadas.clear()

    def volver_a_seleccion(self):
        app.mostrar_pantalla("inicial")


class PantallaParametros(PantallaInicial):
    def __init__(self, root, seqrecord, transcripciones, filepath):
        self.seqrecord = seqrecord
        self.transcripciones = transcripciones
        self.filepath = filepath
        super().__init__(root)
        #self.root.geometry('1000x750')
    
    def crear_interfaz(self):
        #Espaciadores
        self.frame.grid_rowconfigure(4, minsize=70)
        self.frame.grid_rowconfigure(7, minsize=20)
        self.frame.grid_rowconfigure(10, minsize=20)
        self.frame.grid_rowconfigure(13, minsize=30)
        self.frame.grid_rowconfigure(17, minsize=20)
        self.frame.grid_columnconfigure(0, minsize=260)
        self.frame.grid_columnconfigure(1, minsize=200)
        self.frame.grid_columnconfigure(2, minsize=200)
        self.frame.grid_columnconfigure(3, minsize=200)
        self.frame.grid_columnconfigure(4, minsize=20)

        #nombre_label = tk.Label(self.frame, text=f"Nombre de la secuencia: {self.seqrecord.id}")
        nombre_label = ctk.CTkLabel(self.frame, text=f"{self.seqrecord.id} | {self.seqrecord.description}", text_color="#000000")
        nombre_label.grid(row=0, column=0, columnspan=4)

        genes = diseno.get_all_genes(self.seqrecord)
        #descripcion_label = tk.Label(self.frame, text=f"Descripción: {self.seqrecord.description}")
        descripcion_label = ctk.CTkLabel(self.frame, text=f"Genes: {str(genes)[1:-1]}", text_color="#000000")
        descripcion_label.grid(row=1, column=0, columnspan=4)

        #largo_label = tk.Label(self.frame, text=f"Largo de la secuencia: {len(self.seqrecord)}")
        largo_label = ctk.CTkLabel(self.frame, text=f"Largo de la secuencia: {len(self.seqrecord)}", text_color="#000000")
        largo_label.grid(row=2, column=0, columnspan=4)

        locations = []
        for f in self.seqrecord.features:
            if f.type == 'CDS':
                loc = f.location
                if str(loc) not in locations:
                    locations.append(str(loc))

        #transcripciones_label = tk.Label(self.frame, text=f"Se considerará(n) {len(self.transcripciones)} de las {len(locations)} transcripciones.")
        transcripciones_label = ctk.CTkLabel(self.frame, text=f"Se considerará(n) {len(self.transcripciones)} de las {len(locations)} transcripciones. ({len(diseno.get_splicings(self.seqrecord, self.transcripciones))} puntos de empalme)", text_color="#000000")
        transcripciones_label.grid(row=3, column=0, columnspan=4)

        #fuente_negrita = font.Font(weight="bold")
        #parametros_label = tk.Label(self.frame, text="Ingrese los parámetros para el diseño de las sondas:", font=fuente_negrita)
        fuente_titulo = ctk.CTkFont(family='Helvetica', size=20, weight='bold')
        parametros_label = ctk.CTkLabel(self.frame, text="Ingrese los parámetros para el diseño de las sondas:", text_color="#000000", font=fuente_titulo)
        parametros_label.grid(row=4, column=0, columnspan=3, padx=20, pady=10)

        def mostrar_ayuda():
            ayuda_texto = "Usa el cursor para obtener más información de los parámetros"
            messagebox.showinfo("Ayuda", ayuda_texto)

        #boton_ayuda = tk.Button(self.frame, text="?", command=mostrar_ayuda)
        boton_ayuda = ctk.CTkButton(master=self.frame, text="?", fg_color="#7a7a7a", corner_radius=20, width=10, height=10,command=mostrar_ayuda)
        boton_ayuda.grid(row=4, column=3)

        #Largo mínimo
        #minlen_label = tk.Label(self.frame, text="Largo mínimo (nt)")
        minlen_label = ctk.CTkLabel(self.frame, text="Largo mínimo (nt):", text_color="#000000")
        minlen_label.grid(row=5, column=0, padx=5, pady=10, sticky="e")
        
        #minlen_spinbox = tk.Spinbox(self.frame, from_=10, to=300, increment=1, width=5)
        minlen_spinbox = IntegerSelector(self.frame, default_value=60, min_value=10, max_value=500, increment=1)
        minlen_spinbox.grid(row=5, column=1, pady=10)
        #minlen_spinbox.delete(0, "end")
        #minlen_spinbox.insert(0, "60")


        #Largo máximo
        #maxlen_label = tk.Label(self.frame, text="Largo máximo (nt)")
        maxlen_label = ctk.CTkLabel(self.frame, text="Largo máximo (nt):", text_color="#000000")
        maxlen_label.grid(row=6, column=0, padx=5, sticky="e")

        #maxlen_spinbox = tk.Spinbox(self.frame, from_=10, to=300, increment=1, width=5)
        maxlen_spinbox = IntegerSelector(self.frame, default_value=120, min_value=10, max_value=500, increment=1)
        maxlen_spinbox.grid(row=6, column=1, padx= 5)
        #maxlen_spinbox.delete(0, "end")
        #maxlen_spinbox.insert(0, "120")

        texto_largo = "Largo de las sondas: Corresponde a la longitud en\nnucleótidos de la sonda a diseñar. Se prioriza\nel diseño de sondas de menor tamaño."
        tooltip_minlen = ToolTip(minlen_label, texto_largo, resource_path(os.path.join("images","largo.png")))
        tooltip_maxlen = ToolTip(maxlen_label, texto_largo, resource_path(os.path.join("images","largo.png")))

        #TM mínima
        #tmmin_label = tk.Label(self.frame, text="Temp melting mínima (°C)")
        tmmin_label = ctk.CTkLabel(self.frame, text="Temp melting mínima (°C):", text_color="#000000")
        tmmin_label.grid(row=5, column=2, padx=5, pady=10, sticky="e")
        
        #tmmin_spinbox = tk.Spinbox(self.frame, from_=10, to=300, increment=1, width=5)
        tmmin_spinbox = IntegerSelector(self.frame, default_value=65, min_value=10, max_value=300, increment=1)
        tmmin_spinbox.grid(row=5, column=3, padx= 5, pady=10)
        #tmmin_spinbox.delete(0, "end")
        #tmmin_spinbox.insert(0, "65")


        #TM máxima
        #tmmax_label = tk.Label(self.frame, text="Temp melting máxima (°C)")
        tmmax_label = ctk.CTkLabel(self.frame, text="Temp melting máxima (°C):", text_color="#000000")
        tmmax_label.grid(row=6, column=2, padx=5, sticky="e")
        
        #tmmax_spinbox = tk.Spinbox(self.frame, from_=10, to=300, increment=1, width=5)
        tmmax_spinbox = IntegerSelector(self.frame, default_value=80, min_value=10, max_value=300, increment=1)
        tmmax_spinbox.grid(row=6, column=3, padx= 5)
        #tmmax_spinbox.delete(0, "end")
        #tmmax_spinbox.insert(0, "80")

        texto_tm = "La temperatura melting (Tm) o temperatura de fusión de\nuna secuencia se refiere a la temperatura en la cual\nse desnaturaliza la doble hebra. Más en concreto\nes la temperatura en la cual 50 % de las copias de esa secuencia\npresentes en una reacción se encuentra en forma monocatenaria\ny 50 % en forma bicatenaria, que interactúan con su secuencia\ncomplementaria."
        tooltip_tmmin = ToolTip(tmmin_label, texto_tm)
        tooltip_tmmax = ToolTip(tmmax_label, texto_tm)

        #%GC mínimo
        #gcmin_label = tk.Label(self.frame, text="%GC mínimo")
        gcmin_label = ctk.CTkLabel(self.frame, text="%GC mínimo:", text_color="#000000")
        gcmin_label.grid(row=8, column=0, padx=5, pady=10, sticky="e")
        
        #gcmin_spinbox = tk.Spinbox(self.frame, from_=0, to=100, increment=1, width=5)
        gcmin_spinbox = IntegerSelector(self.frame, default_value=30, min_value=0, max_value=100, increment=1)
        gcmin_spinbox.grid(row=8, column=1, padx= 5, pady=10)
        #gcmin_spinbox.delete(0, "end")
        #gcmin_spinbox.insert(0, "30")


        #%GC máximo
        #gcmax_label = tk.Label(self.frame, text="%GC máximo")
        gcmax_label = ctk.CTkLabel(self.frame, text="%GC máximo:", text_color="#000000")
        gcmax_label.grid(row=9, column=0, padx=5, sticky="e")

        #gcmax_spinbox = tk.Spinbox(self.frame, from_=0, to=100, increment=1, width=5)
        gcmax_spinbox = IntegerSelector(self.frame, default_value=70, min_value=0, max_value=100, increment=1)
        gcmax_spinbox.grid(row=9, column=1, padx= 5)
        #gcmax_spinbox.delete(0, "end")
        #gcmax_spinbox.insert(0, "70")

        texto_gc = "Porcentaje de GC: Se refiere a la fracción de bases\nnitrogenadas que son citosinas o guaninas dentro de\nla secuencia de nucleótidos."
        tooltip_gcmin = ToolTip(gcmin_label, texto_gc)
        tooltip_gcmax = ToolTip(gcmax_label, texto_gc)

        #Distancia mínima al borde del exón
        #mindist_label = tk.Label(self.frame, text="Distancia mínima al borde del exón (nt)")
        mindist_label = ctk.CTkLabel(self.frame, text="Distancia mínima al borde del exón (nt):", text_color="#000000")
        mindist_label.grid(row=8, column=2, padx=5, pady=10, sticky="e")
        
        #mindist_spinbox = tk.Spinbox(self.frame, from_=0, to=500, increment=1, width=5)
        mindist_spinbox = IntegerSelector(self.frame, default_value=0, min_value=0, max_value=500, increment=1)
        mindist_spinbox.grid(row=8, column=3, padx= 5, pady=10)
        #mindist_spinbox.delete(0, "end")
        #mindist_spinbox.insert(0, "0")


        #Distancia máxima al borde del exón
        #maxdist_label = tk.Label(self.frame, text="Distancia máxima al borde del exón (nt)")
        maxdist_label = ctk.CTkLabel(self.frame, text="Distancia máxima al borde del exón (nt):", text_color="#000000")
        maxdist_label.grid(row=9, column=2, padx=5, sticky="e")
        
        #maxdist_spinbox = tk.Spinbox(self.frame, from_=0, to=500, increment=1, width=5)
        maxdist_spinbox = IntegerSelector(self.frame, default_value=50, min_value=0, max_value=500, increment=1)
        maxdist_spinbox.grid(row=9, column=3, padx= 5)
        #maxdist_spinbox.delete(0, "end")
        #maxdist_spinbox.insert(0, "100")

        texto_dist = "La distancia al borde del exón se refiere a la cantidad de nucleótidos que hay entre\nel punto de empalme y las sondas internas donor y acceptor.\n(No aplica para sondas centrales)"
        tooltip_mindist = ToolTip(mindist_label, texto_dist, resource_path(os.path.join("images","distancia.png")))
        tooltip_maxdist = ToolTip(maxdist_label, texto_dist, resource_path(os.path.join("images","distancia.png")))


        #Sobrelape mínimo
        #minoverlap_label = tk.Label(self.frame, text="Sobrelape mínimo (%)")
        minoverlap_label = ctk.CTkLabel(self.frame, text="Sobrelape mínimo (%):", text_color="#000000")
        minoverlap_label.grid(row=11, column=0, padx=5, pady=10, sticky="e")
        
        #minoverlap_spinbox = tk.Spinbox(self.frame, from_=0, to=100, increment=1, width=5)
        minoverlap_spinbox = IntegerSelector(self.frame, default_value=25, min_value=0, max_value=100, increment=1)
        minoverlap_spinbox.grid(row=11, column=1, padx= 5, pady=10)
        #minoverlap_spinbox.delete(0, "end")
        #minoverlap_spinbox.insert(0, "25")


        #Sobrelape máximo
        #maxoverlap_label = tk.Label(self.frame, text="Sobrelape máximo (%)")
        maxoverlap_label = ctk.CTkLabel(self.frame, text="Sobrelape máximo (%):", text_color="#000000")
        maxoverlap_label.grid(row=12, column=0, padx=5, sticky="e")

        maxoverlap_spinbox = tk.Spinbox(self.frame, from_=0, to=100, increment=1, width=5)
        maxoverlap_spinbox = IntegerSelector(self.frame, default_value=50, min_value=0, max_value=100, increment=1)
        maxoverlap_spinbox.grid(row=12, column=1, padx= 5)
        #maxoverlap_spinbox.delete(0, "end")
        #maxoverlap_spinbox.insert(0, "50")

        texto_overlap = "El sobrelape u overlap es el porcentaje de la\nsonda interior que es compartido con las sondas\nexteriores, solo aplica para donor y acceptor."
        tooltip_minoverlap = ToolTip(minoverlap_label, texto_overlap, resource_path(os.path.join("images","sobrelape.png")))
        tooltip_maxoverlap = ToolTip(maxoverlap_label, texto_overlap, resource_path(os.path.join("images","sobrelape.png")))

        #Delta G mínimo homodimerización
        #dgmin_homodim_label = tk.Label(self.frame, text="Delta G mínimo homodimerización")
        dgmin_homodim_label = ctk.CTkLabel(self.frame, text="Delta G mínimo homodimerización:", text_color="#000000")
        dgmin_homodim_label.grid(row=11, column=2, padx=5, pady=10, sticky="e")
        
        #dgmin_homodim_spinbox = tk.Spinbox(self.frame, from_=-50000, to=10000, increment=1000, width=5)
        dgmin_homodim_spinbox = IntegerSelector(self.frame, default_value=-15000, min_value=-999999, max_value=1000000, increment=1000)
        dgmin_homodim_spinbox.grid(row=11, column=3, padx= 5, pady=10)
        #dgmin_homodim_spinbox.delete(0, "end")
        #dgmin_homodim_spinbox.insert(0, "-10000")


        #Delta G mínimo hairpin
        #dgmin_hairpin_label = tk.Label(self.frame, text="Delta G mínimo hairpin")
        dgmin_hairpin_label = ctk.CTkLabel(self.frame, text="Delta G mínimo hairpin:", text_color="#000000")
        dgmin_hairpin_label.grid(row=12, column=2, padx=5, sticky="e")
        
        #dgmin_hairpin_spinbox = tk.Spinbox(self.frame, from_=-50000, to=10000, increment=1000, width=6)
        dgmin_hairpin_spinbox = IntegerSelector(self.frame, default_value=-15000, min_value=-999999, max_value=1000000, increment=1000)
        dgmin_hairpin_spinbox.grid(row=12, column=3, padx= 5)
        #dgmin_hairpin_spinbox.delete(0, "end")
        #dgmin_hairpin_spinbox.insert(0, "-10000")

        texto_dgmin_homodim = "Se refiere a la evaluación de la capacidad\nde la secuencia para formar estructuras de dímeros\nhomólogos. Los dímeros homólogos son la unión de\ndos secuencias de ADN idénticas o muy similares\nentre sí. Básicamente que la sonda no se complemente\ncon sí misma. Esto se verifica con la librería\nprimer3, que calcula la energía máxima entre las\ndos secuencias iguales."
        tooltip_dgmin_homodim = ToolTip(dgmin_homodim_label, texto_dgmin_homodim)
        texto_hairpin_label = "Capacidad de una secuencia para formar estructuras\nde horquilla (hairpin) en sí misma. Una estructura\nde hairpin es una conformación en la que una región\nde la secuencia se pliega hacia atrás y se empareja\ncon su complementaria, formando una estructura en forma de\nhorquilla."
        tooltip_hairpin_label = ToolTip(dgmin_hairpin_label, texto_hairpin_label)

        #Máximo homopolímeros simples
        #maxhomopol_simple_label = tk.Label(self.frame, text="Máximo homopolímeros simples")
        maxhomopol_simple_label = ctk.CTkLabel(self.frame, text="Máximo homopolímeros simples:", text_color="#000000")
        maxhomopol_simple_label.grid(row=14, column=0, padx=5, sticky="e")
        
        #maxhomopol_simple_spinbox = tk.Spinbox(self.frame, from_=2, to=100, increment=1, width=5)
        maxhomopol_simple_spinbox = IntegerSelector(self.frame, default_value=6, min_value=2, max_value=1000, increment=1)
        maxhomopol_simple_spinbox.grid(row=14, column=1, padx= 5)
        #maxhomopol_simple_spinbox.delete(0, "end")
        #maxhomopol_simple_spinbox.insert(0, "6")


        #Máximo homopolímeros dobles
        #maxhomopol_double_label = tk.Label(self.frame, text="Máximo homopolímeros dobles")
        maxhomopol_double_label = ctk.CTkLabel(self.frame, text="Máximo homopolímeros dobles:", text_color="#000000")
        maxhomopol_double_label.grid(row=15, column=0, padx=5, pady=10, sticky="e")

        #maxhomopol_double_spinbox = tk.Spinbox(self.frame, from_=2, to=100, increment=1, width=5)
        maxhomopol_double_spinbox = IntegerSelector(self.frame, default_value=5, min_value=2, max_value=1000, increment=1)
        maxhomopol_double_spinbox.grid(row=15, column=1, padx= 5, pady=10)
        #maxhomopol_double_spinbox.delete(0, "end")
        #maxhomopol_double_spinbox.insert(0, "5")

        #Máximo homopolímeros triples
        #maxhomopol_triple_label = tk.Label(self.frame, text="Máximo homopolímeros triples")
        maxhomopol_triple_label = ctk.CTkLabel(self.frame, text="Máximo homopolímeros triples:", text_color="#000000")
        maxhomopol_triple_label.grid(row=16, column=0, padx=5, sticky="e")

        #maxhomopol_triple_spinbox = tk.Spinbox(self.frame, from_=2, to=100, increment=1, width=5)
        maxhomopol_triple_spinbox = IntegerSelector(self.frame, default_value=5, min_value=2, max_value=1000, increment=1)
        maxhomopol_triple_spinbox.grid(row=16, column=1, padx= 5)
        #maxhomopol_triple_spinbox.delete(0, "end")
        #maxhomopol_triple_spinbox.insert(0, "4")

        texto_maxhomopol_simple = "Los homopolímeros simples son secuencias donde se repite una base\nde manera consecutiva: si el máximo de homopolímeros simples\nes 3, se aceptaría una subsecuencia 'GGG' pero no una\nsubsecuencia 'GGGG'"
        tooltip_maxhomopol_simple = ToolTip(maxhomopol_simple_label, texto_maxhomopol_simple)
        texto_maxhomopol_double = "Los homopolímeros dobles se refieren a las subsecuencias de dos\nbases que se repiten más de una vez, como por\nejemplo AGAGAG. Si el máximo de homopolímeros dobles\nes 3, se aceptaría una subsecuencia 'AGAGAG' pero no una\nsubsecuencia 'AGAGAGAG'"
        tooltip_maxhomopol_double = ToolTip(maxhomopol_double_label, texto_maxhomopol_double)
        texto_maxhomopol_triple = "Los homopolímeros triples se refieren a las subsecuencias de tres\nbases que se repiten más de una vez, como por\nejemplo AGCAGCAGC. Si el máximo de homopolímeros triples\nes 3, se aceptaría una subsecuencia 'AGCAGCAGC' pero no una\nsubsecuencia 'AGCAGCAGCAGC'"
        tooltip_maxhomopol_triple = ToolTip(maxhomopol_triple_label, texto_maxhomopol_triple)

        def param_multiplex():
            if multiplex_var.get() == "on":
                mindg_label.grid(row=15, column=2, padx=5, sticky="e")
                mindg_spinbox.grid(row=15, column=3, padx= 5)
                maxdt_label.grid(row=16, column=2, padx=5, sticky="e")
                maxdt_spinbox.grid(row=16, column=3, padx= 5)
            else:
                mindg_label.grid_forget()
                mindg_spinbox.grid_forget()
                maxdt_label.grid_forget()
                maxdt_spinbox.grid_forget()


        multiplex_var = ctk.StringVar(value="off")
        checkbox_multiplex = ctk.CTkCheckBox(self.frame, text="Multiplexar sondas", command=param_multiplex, variable=multiplex_var, onvalue="on", offvalue="off", fg_color="#000000", text_color="#000000", hover_color="#7a7a7a")
        checkbox_multiplex.deselect()
        checkbox_multiplex.grid(row=14, column=2, padx= 5,sticky='e')

        mindg_label = ctk.CTkLabel(self.frame, text="Mínimo delta G heterodimerización:", text_color="#000000")
        mindg_spinbox = IntegerSelector(self.frame, default_value=-20000, min_value=-999999, max_value=1000000, increment=1000)
    
        maxdt_label = ctk.CTkLabel(self.frame, text="Máxima diferencia de Tm (°C):", text_color="#000000")
        maxdt_spinbox = IntegerSelector(self.frame, default_value=5, min_value=0, max_value=100, increment=1)
        
        texto_multiplex = "La multiplexación consiste en agrupar las sondas\npara ser procesadas en conjunto. Los criterios para\nagrupar son: diferencia de temperatura de melting\ny heterodimerización.\nNOTA: El tiempo de carga aumentará considerablemente."
        tooltip_multiplex = ToolTip(checkbox_multiplex, texto_multiplex)

        texto_mindg = "Para que las sondas se procesen en conjunto, es muy\nnecesario que no se produzca dimerización entre las sondas.\nEsto quiere decir que las secuencias generadas\ndeben tener muy baja probabilidad de que se complementen\nentre sí. Es por esto que se establece un límite para\nla diferencia de energía producto de la dimerización."
        tooltip_mindg = ToolTip(mindg_label, texto_mindg)
        texto_maxdt = "Diferencia de temperatura de melting: Este factor\nes determinante debido a que es muy necesario que las\nsondas funcionen de una manera similar. Si tienen\nTM muy diferentes, algunas sondas pueden no funcionar\neficientemente a la temperatura de reacción, mientras\nque otras lo harán de manera óptima."
        tooltip_maxdt = ToolTip(maxdt_label, texto_maxdt)
        

        def ejecutar():
            dict_params = {}
            dict_params['minlen'] = int(minlen_spinbox.valor.get())
            dict_params['maxlen'] = int(maxlen_spinbox.valor.get())
            dict_params['tmmin'] = int(tmmin_spinbox.valor.get())
            dict_params['tmmax'] = int(tmmax_spinbox.valor.get())
            dict_params['gcmin'] = int(gcmin_spinbox.valor.get())
            dict_params['gcmax'] = int(gcmax_spinbox.valor.get())
            dict_params['mindist'] = int(mindist_spinbox.valor.get())
            dict_params['maxdist'] = int(maxdist_spinbox.valor.get())
            dict_params['minoverlap'] = int(minoverlap_spinbox.valor.get())
            dict_params['maxoverlap'] = int(maxoverlap_spinbox.valor.get())
            dict_params['dgmin_homodim'] = int(dgmin_homodim_spinbox.valor.get())
            dict_params['dgmin_hairpin'] = int(dgmin_hairpin_spinbox.valor.get())
            dict_params['maxhomopol_simple'] = int(maxhomopol_simple_spinbox.valor.get())
            dict_params['maxhomopol_double'] = int(maxhomopol_double_spinbox.valor.get())
            dict_params['maxhomopol_triple'] = int(maxhomopol_triple_spinbox.valor.get())
            dict_params['multiplex'] = (multiplex_var.get() == "on")
            dict_params['mindg'] = int(mindg_spinbox.valor.get())
            dict_params['maxdt'] = int(maxdt_spinbox.valor.get())
            app.mostrar_pantalla("carga", seqrecord=self.seqrecord, transcripciones=self.transcripciones, dict_params=dict_params, filepath=self.filepath)

        #volver_button = tk.Button(self.frame, text="Volver", command=self.volver_a_transcripciones)
        #continuar_button = tk.Button(self.frame, text="Continuar a la ejecución", command=ejecutar)
        volver_button = ctk.CTkButton(master=self.frame, text="Volver", fg_color="#7a7a7a", corner_radius=30, width=70 ,command=self.volver_a_transcripciones)
        continuar_button = ctk.CTkButton(master=self.frame, text="Continuar a la ejecución", fg_color="#404040", corner_radius=30,command=ejecutar)
        volver_button.grid(row=20, column=1, pady=20)
        continuar_button.grid(row=20, column=2, pady=20)

        #self.root.geometry("800x1000")


    def volver_a_transcripciones(self):
        app.mostrar_pantalla("transcripciones", seqrecord=self.seqrecord, filepath=self.filepath)


class PantallaCarga(PantallaInicial):
    def __init__(self, root, seqrecord, transcripciones, dict_params, filepath):
        self.seqrecord = seqrecord
        self.transcripciones = transcripciones
        self.dict_params = dict_params
        self.filepath = filepath
        super().__init__(root)


    def crear_interfaz(self):

        self.load_label = ctk.CTkLabel(self.frame, text="Buscando y verificando sondas, por favor espere...", text_color="#000000")
        self.load_label.pack(padx=10, pady=10)
        self.progress_bar = ttk.Progressbar(self.frame, mode="determinate", maximum=100, length=300)
        self.progress_bar.pack(padx=10, pady=10)
        
        if self.dict_params['multiplex']:
            #multiplex_label = ctk.CTkLabel(self.frame, text="Multiplexando sondas, por favor espere...", text_color="#000000")
            #multiplex_label.pack(padx=10, pady=10)
            
            self.progress_multiplex = ttk.Progressbar(self.frame, mode="determinate", maximum=100, length=300)
            self.progress_multiplex.pack(padx=10, pady=10)

        progreso_queue = queue.Queue()

        tarea_thread = threading.Thread(target=self.ejecutar_diseno, args=(progreso_queue,))
        tarea_thread.start()

        actualizar_progreso_thread = threading.Thread(target=self.actualizar_progreso, args=(progreso_queue,tarea_thread,self.dict_params['multiplex'],))
        actualizar_progreso_thread.start()

        #reporte_label = ctk.CTkLabel(self.frame, text="Generando reporte e imagen, por favor espere...", text_color="#000000")
        #reporte_label.pack(padx=10, pady=10)

        self.root.geometry("600x300")
    
    def ejecutar_diseno(self, progreso_queue):
        inicio = time.time()
        df = diseno.probe_designer(self.seqrecord,
                                    self.transcripciones,
                                    progreso_queue,
                                    minlen=self.dict_params['minlen'],
                                    maxlen=self.dict_params['maxlen'],
                                    tmmin=self.dict_params['tmmin'],
                                    tmmax=self.dict_params['tmmax'],
                                    gcmin=self.dict_params['gcmin'],
                                    gcmax=self.dict_params['gcmax'],
                                    mindist=self.dict_params['mindist'],
                                    maxdist=self.dict_params['maxdist'],
                                    minoverlap=self.dict_params['minoverlap'],
                                    maxoverlap=self.dict_params['maxoverlap'],
                                    dgmin_homodim=self.dict_params['dgmin_homodim'],
                                    dgmin_hairpin=self.dict_params['dgmin_hairpin'],
                                    maxhomopol_simple=self.dict_params['maxhomopol_simple'],
                                    maxhomopol_double=self.dict_params['maxhomopol_double'],
                                    maxhomopol_triple=self.dict_params['maxhomopol_triple'],
                                    ismultiplex=self.dict_params['multiplex'],
                                    mindg=self.dict_params['mindg'],
                                    maxdt=self.dict_params['maxdt'])
        fin = time.time()
        filename = diseno.generate_xlsx(df=df,
                                        name=self.seqrecord.id,
                                        genes=str(diseno.get_all_genes(self.seqrecord))[1:-1],
                                        minlen=self.dict_params['minlen'],
                                        maxlen=self.dict_params['maxlen'],
                                        tmmin=self.dict_params['tmmin'],
                                        tmmax=self.dict_params['tmmax'],
                                        gcmin=self.dict_params['gcmin'],
                                        gcmax=self.dict_params['gcmax'],
                                        mindist=self.dict_params['mindist'],
                                        maxdist=self.dict_params['maxdist'],
                                        minoverlap=self.dict_params['minoverlap'],
                                        maxoverlap=self.dict_params['maxoverlap'],
                                        dgmin_homodim=self.dict_params['dgmin_homodim'],
                                        dgmin_hairpin=self.dict_params['dgmin_hairpin'],
                                        maxhomopol_simple=self.dict_params['maxhomopol_simple'],
                                        maxhomopol_double=self.dict_params['maxhomopol_double'],
                                        maxhomopol_triple=self.dict_params['maxhomopol_triple'],
                                        multiplex=self.dict_params['multiplex'],
                                        maxdt=self.dict_params['maxdt'],
                                        mindg=self.dict_params['mindg'],
                                        tiempo_diseno=int(fin-inicio))
        
        imagen.plot_isoforms(self.seqrecord,
                             self.transcripciones,
                             imagen.get_splicings(self.seqrecord, self.transcripciones),
                             filename)
        
        #folderpath = os.path.join(os.getcwd(),'sondas',filename)
        folderpath = resource_path(os.path.join('sondas',filename))
        
        shutil.copy(self.filepath, folderpath)

        #df1 = pd.read_excel(os.path.join(folderpath, self.seqrecord.id+'.xlsx'), sheet_name='Resumen', header=None)
        #df2 = pd.read_excel(os.path.join(folderpath, self.seqrecord.id+'.xlsx'), sheet_name='Parámetros iniciales', header=None)
        #df3 = pd.read_excel(os.path.join(folderpath, self.seqrecord.id+'.xlsx'), sheet_name='Grupos', header=None)
        
        app.mostrar_pantalla(nombre="final", seqrecord=self.seqrecord, transcripciones=self.transcripciones, df=df, dict_params=self.dict_params, filepath=self.filepath, folderpath=folderpath)#, df1=df1, df2=df2, df3=df3)
        
    
    def actualizar_progreso(self, progreso_queue, tarea_thread, multiplex):
        while True:
            try:
                progreso_actual = progreso_queue.get_nowait()
                if not multiplex:
                    if progreso_actual < 95:
                        self.progress_bar["value"] = progreso_actual
                    else:
                        self.progress_bar["value"] = 100
                        self.load_label.configure(text="Generando reporte e imagen, por favor espere...")
                else:
                    if progreso_actual < 100:
                        self.progress_bar["value"] = progreso_actual
                    else:
                        self.progress_bar["value"] = 100
                        self.load_label.configure(text="Multiplexando sondas, por favor espere...")
                        self.progress_multiplex["value"] = progreso_actual - 100
                        if progreso_actual > 195:
                                self.load_label.configure(text="Generando reporte e imagen, esto podría tomar unos minutos...")
                root.update_idletasks()
            except queue.Empty:
                if tarea_thread.is_alive():
                    continue
                break


class PantallaFinal(PantallaInicial):
    def __init__(self, root, seqrecord, transcripciones, df, dict_params, filepath, folderpath, df1, df2, df3):
        self.seqrecord = seqrecord
        self.transcripciones = transcripciones
        self.df = df
        self.dict_params = dict_params
        self.filepath = filepath
        self.folderpath = folderpath
        self.df1 = df1
        self.df2 = df2
        self.df3 = df3
        super().__init__(root)

    def crear_interfaz(self):
    
        sondas = []
        for index, row in self.df.iterrows():
            if row['sonda'] != 'AAA' and row['sonda'] not in sondas:
                sondas.append(row['sonda'])

        df_resumen = self.df.drop_duplicates(subset='sonda')
        df_resumen = df_resumen[df_resumen['sonda'] != "AAA"]

        #fuente_negrita = font.Font(weight="bold")
        #resultado_label = tk.Label(self.frame, text=f"El diseño se ejecutó con éxito. Número de sondas: {len(sondas)}", font=fuente_negrita)
        fuente_negrita = ctk.CTkFont(family='Helvetica', weight='bold')
        resultado_label = ctk.CTkLabel(self.frame, text=f"Sondas diseñadas para {str(diseno.get_all_genes(self.seqrecord))[1:-1]}", text_color="#000000", font=fuente_negrita)
        resultado_label.grid(row=0, column=0, columnspan=2, padx=20, pady=20)

        # img = Image.open(os.path.join(folderpath,self.seqrecord.id+".png"))
        # tkimg = ImageTk.PhotoImage(img)
        # label_imagen = tk.Label(self.frame, image=tkimg)
        # label_imagen.image = tkimg
        # label_imagen.pack(padx=10, pady=10)

        #carpeta_label = tk.Label(self.frame, text=f"El reporte y la imagen se almacenaron en\n" + folderpath)
        carpeta_label = ctk.CTkLabel(self.frame, text=f"El reporte, la imagen y el archivo GenBank se almacenaron en\n" + self.folderpath, text_color="#000000")
        carpeta_label.grid(row=1, column=0, columnspan=2, padx=20, pady=20)

        def abrir_carpeta():
            sistema_operativo = platform.system()
            if sistema_operativo == 'Windows':
                os.system(f'explorer "{self.folderpath}"')
            elif sistema_operativo == 'Darwin':  # macOS
                os.system(f'open "{self.folderpath}"')
            elif sistema_operativo == 'Linux':
                os.system(f'xdg-open "{self.folderpath}"')

        boton_carpeta = ctk.CTkButton(self.frame, text="Abrir carpeta", corner_radius=30, fg_color="#404040", command=abrir_carpeta)
        boton_carpeta.grid(row=2, column=0, columnspan=2, padx=20, pady=20)

        img_frame = ttk.Frame(self.frame, width=900, height=700)
        img_frame.grid(row=3, column=0, padx=10, pady=20)
        Zoom_Advanced(mainframe=img_frame, path=os.path.join(self.folderpath, self.seqrecord.id+'.png'))
        #ScrollView(img_frame)

        tabview = ttk.Notebook(self.frame, width=400, height=300)
        tab1 = ttk.Frame(tabview)
        tab2 = ttk.Frame(tabview)
        tabview.add(tab1, text='Resumen')
        tabview.add(tab2, text='Parámetros')
        tabview.grid(row=3, column=1, padx=10, pady=20)
        
        #tabview_frame = ttk.Frame(self.frame, width=500, height=500)
        #tabview = ctk.CTkTabview(tabview_frame)
        #tabview.configure(height=600)
        #tab1 = tabview.add('Resumen')
        #tab2 = tabview.add('Parámetros')
        #tab3 = tabview.add('Grupos')
        #tabview.grid(row=2, column=1, padx=10, pady=20)



        #wb_reporte = openpyxl.load_workbook(os.path.join(self.folderpath, self.seqrecord.id+'.xlsx'), read_only=True)


        def load_data():
            df1 = pd.read_excel(os.path.join(self.folderpath, self.seqrecord.id+'.xlsx'), sheet_name='Resumen', header=None)
            df2 = pd.read_excel(os.path.join(self.folderpath, self.seqrecord.id+'.xlsx'), sheet_name='Parámetros iniciales', header=None)
            df3 = pd.read_excel(os.path.join(self.folderpath, self.seqrecord.id+'.xlsx'), sheet_name='Grupos', header=None)
            return df1, df2, df3


        def on_load_data():
            #df1, df2, df3 = self.df1, self.df2, self.df3

            ###### Inicio Tab 1: Resumen #########

            canvas1 = tk.Canvas(tab1)
            y_scrollbar = tk.Scrollbar(tab1, orient='vertical', command=canvas1.yview)
            y_scrollbar.pack(side='right', fill='y')
            canvas1.configure(yscrollcommand=y_scrollbar.set)
            canvas1.pack(side='left', fill='both', expand=True)

            #hoja = wb_reporte['Resumen']

            #ctk.CTkLabel(tab1, text=hoja['A1'].value, text_color="#000000").grid(row=0, column=0)
            #ctk.CTkLabel(tab1, text=df1.iat[0,0], text_color="#000000").grid(row=0, column=0, padx=5, sticky='e')
            tk.Label(canvas1,text="Genes: ").grid(row=0, column=0, sticky='e')
            tk.Label(canvas1,text=str(diseno.get_all_genes(self.seqrecord))[1:-1]).grid(row=0, column=1, sticky='w', padx=10)


            canvas1.grid_rowconfigure(1, minsize=10)

            #ctk.CTkLabel(tab1, text=df1.iat[2,0], text_color="#000000").grid(row=2, column=0, sticky='e')
            #ctk.CTkLabel(tab1, text=df1.iat[2,1], text_color="#000000").grid(row=2, column=1, sticky='w')
            tk.Label(canvas1,text='Número deseable de sondas').grid(row=2, column=0, sticky='e')
            tk.Label(canvas1,text=self.df.shape[0]).grid(row=2, column=1, sticky='w', padx=10)

            #ctk.CTkLabel(tab1, text=df1.iat[3,0], text_color="#000000").grid(row=3, column=0, sticky='e')
            #ctk.CTkLabel(tab1, text=df1.iat[3,1], text_color="#000000").grid(row=3, column=1, sticky='w')
            tk.Label(canvas1,text='Número de sondas válidas').grid(row=3, column=0, sticky='e')
            tk.Label(canvas1,text=self.df[self.df['sonda'] != "AAA"].shape[0]).grid(row=3, column=1, sticky='w', padx=10)

            #ctk.CTkLabel(tab1, text=df1.iat[4,0], text_color="#000000").grid(row=4, column=0, sticky='e')
            #ctk.CTkLabel(tab1, text=df1.iat[4,1], text_color="#000000").grid(row=4, column=1, sticky='w')
            tk.Label(canvas1,text='Número de sondas diferentes').grid(row=4, column=0, sticky='e')
            tk.Label(canvas1,text=df_resumen.shape[0]).grid(row=4, column=1, sticky='w', padx=10)

            canvas1.grid_rowconfigure(5, minsize=10)

            #ctk.CTkLabel(tab1, text=df1.iat[6,0], text_color="#000000").grid(row=6, column=0, sticky='e')
            #ctk.CTkLabel(tab1, text=df1.iat[6,1], text_color="#000000").grid(row=6, column=1, sticky='w')
            tk.Label(canvas1,text='Tm promedio').grid(row=6, column=0, sticky='e')
            tk.Label(canvas1,text=str("{:.1f}").format(df_resumen['tm'].mean())+'°').grid(row=6, column=1, sticky='w', padx=10)


            #ctk.CTkLabel(tab1, text=df1.iat[7,0], text_color="#000000").grid(row=7, column=0, sticky='e')
            #ctk.CTkLabel(tab1, text=df1.iat[7,1], text_color="#000000").grid(row=7, column=1, sticky='w')
            tk.Label(canvas1,text='%GC promedio').grid(row=7, column=0, sticky='e')
            tk.Label(canvas1,text=str("{:.1f}").format(df_resumen['gc'].mean())+'%').grid(row=7, column=1, sticky='w', padx=10)


            #ctk.CTkLabel(tab1, text=df1.iat[8,0], text_color="#000000").grid(row=8, column=0, sticky='e')
            #ctk.CTkLabel(tab1, text=df1.iat[8,1], text_color="#000000").grid(row=8, column=1, sticky='w')
            tk.Label(canvas1,text='Largo promedio').grid(row=8, column=0, sticky='e')
            tk.Label(canvas1,text=str("{:.1f}").format(df_resumen['largo'].mean())+' nt').grid(row=8, column=1, sticky='w', padx=10)


            canvas1.grid_rowconfigure(9, minsize=10)

            #ctk.CTkLabel(tab1, text=df1.iat[10,0], text_color="#000000").grid(row=10, column=0, sticky='e')
            #ctk.CTkLabel(tab1, text=df1.iat[10,1], text_color="#000000").grid(row=10, column=1, sticky='w')
            tk.Label(canvas1,text='Largo promedio').grid(row=8, column=0, sticky='e')
            tk.Label(canvas1,text=str("{:.1f}").format(df_resumen['largo'].mean())+' nt').grid(row=8, column=1, sticky='w', padx=10)

            if self.dict_params['multiplex']:
                canvas1.grid_rowconfigure(11, minsize=10)

                #ctk.CTkLabel(tab1, text=df1.iat[12,0], text_color="#000000").grid(row=12, column=0, sticky='e')
                #ctk.CTkLabel(tab1, text=df1.iat[12,1], text_color="#000000").grid(row=12, column=1, sticky='w')
                tk.Label(canvas1,text='Número de grupos (multiplex)').grid(row=12, column=0, sticky='e')
                tk.Label(canvas1,text=int(df_resumen['grupo'].max())).grid(row=12, column=1, sticky='w', padx=10)


            ###### Fin Tab 1: Resumen #########
            
            ###### Inicio Tab 2: Parámetros #########

            #hoja = wb_reporte['Parámetros iniciales']
            
            def on_configure(event):
                canvas2.configure(scrollregion=canvas2.bbox('all'))
                canvas_height = frame.winfo_reqheight()
                if canvas_height != canvas2.winfo_height():
                    canvas2.config(height=canvas_height)
                    scrollbar.config(command=canvas2.yview)
            
            scrollbar = tk.Scrollbar(tab2)
            scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
            canvas2 = tk.Canvas(tab2)
            canvas2.pack(side='left', fill='both', expand=True)
            #canvas2.grid(row=0,column=0)
            scrollbar.config(command=canvas2.yview)

            def scroll(event):
                canvas2.yview_scroll(-1*(event.delta//120), "units")
            canvas2.bind_all("<MouseWheel>", scroll)

            frame = tk.Frame(canvas2)
            canvas2.create_window((0, 0), window=frame, anchor='nw')


            frame.bind("<Configure>", on_configure)


            #ctk.CTkLabel(tab2, text=df2.iat[2,0], text_color="#000000").grid(row=0, column=0, sticky='e')
            #ctk.CTkLabel(tab2, text=df2.iat[2,1], text_color="#000000").grid(row=0, column=1, sticky='e')
            tk.Label(frame,text='Largo de sonda mínimo').grid(row=0, column=0, sticky='e')
            tk.Label(frame,text=self.dict_params['minlen']).grid(row=0, column=1, sticky='w', padx=10)

            #ctk.CTkLabel(tab2, text=df2.iat[3,0], text_color="#000000").grid(row=1, column=0, sticky='e')
            #ctk.CTkLabel(tab2, text=df2.iat[3,1], text_color="#000000").grid(row=1, column=1, sticky='e')
            tk.Label(frame,text='Largo de sonda máximo').grid(row=1, column=0, sticky='e')
            tk.Label(frame,text=self.dict_params['maxlen']).grid(row=1, column=1, sticky='w', padx=10)

            frame.grid_rowconfigure(2, minsize=10)

            #ctk.CTkLabel(tab2, text=df2.iat[5,0], text_color="#000000").grid(row=3, column=0, sticky='e')
            #ctk.CTkLabel(tab2, text=df2.iat[5,1], text_color="#000000").grid(row=3, column=1, sticky='e')
            tk.Label(frame,text='Temp melting mínima').grid(row=3, column=0, sticky='e')
            tk.Label(frame,text=self.dict_params['tmmin']).grid(row=3, column=1, sticky='w', padx=10)

            #ctk.CTkLabel(tab2, text=df2.iat[6,0], text_color="#000000").grid(row=4, column=0, sticky='e')
            #ctk.CTkLabel(tab2, text=df2.iat[6,1], text_color="#000000").grid(row=4, column=1, sticky='e')
            tk.Label(frame,text='Temp melting máxima').grid(row=4, column=0, sticky='e')
            tk.Label(frame,text=self.dict_params['tmmax']).grid(row=4, column=1, sticky='w', padx=10)

            frame.grid_rowconfigure(5, minsize=10)

            #ctk.CTkLabel(tab2, text=df2.iat[8,0], text_color="#000000").grid(row=6, column=0, sticky='e')
            #ctk.CTkLabel(tab2, text=df2.iat[8,1], text_color="#000000").grid(row=6, column=1, sticky='e')
            tk.Label(frame,text='%GC mínimo').grid(row=6, column=0, sticky='e')
            tk.Label(frame,text=self.dict_params['gcmin']).grid(row=6, column=1, sticky='w', padx=10)

            #ctk.CTkLabel(tab2, text=df2.iat[9,0], text_color="#000000").grid(row=7, column=0, sticky='e')
            #ctk.CTkLabel(tab2, text=df2.iat[9,1], text_color="#000000").grid(row=7, column=1, sticky='e')
            tk.Label(frame,text='%GC máximo').grid(row=7, column=0, sticky='e')
            tk.Label(frame,text=self.dict_params['gcmax']).grid(row=7, column=1, sticky='w', padx=10)

            frame.grid_rowconfigure(8, minsize=10)

            #ctk.CTkLabel(tab2, text=df2.iat[11,0], text_color="#000000").grid(row=9, column=0, sticky='e')
            #ctk.CTkLabel(tab2, text=df2.iat[11,1], text_color="#000000").grid(row=9, column=1, sticky='e')
            tk.Label(frame,text='Dist mínima al borde del exón').grid(row=9, column=0, sticky='e')
            tk.Label(frame,text=self.dict_params['mindist']).grid(row=9, column=1, sticky='w', padx=10)

            #ctk.CTkLabel(tab2, text=df2.iat[12,0], text_color="#000000").grid(row=10, column=0, sticky='e')
            #ctk.CTkLabel(tab2, text=df2.iat[12,1], text_color="#000000").grid(row=10, column=1, sticky='e')
            tk.Label(frame,text='Dist máxima al borde del exón').grid(row=10, column=0, sticky='e')
            tk.Label(frame,text=self.dict_params['maxdist']).grid(row=10, column=1, sticky='w', padx=10)

            frame.grid_rowconfigure(11, minsize=10)

            #ctk.CTkLabel(tab2, text=df2.iat[14,0], text_color="#000000").grid(row=12, column=0, sticky='e')
            #ctk.CTkLabel(tab2, text=df2.iat[14,1], text_color="#000000").grid(row=12, column=1, sticky='e')
            tk.Label(frame,text='Sobrelape mínimo').grid(row=12, column=0, sticky='e')
            tk.Label(frame,text=self.dict_params['minoverlap']).grid(row=12, column=1, sticky='w', padx=10)

            #ctk.CTkLabel(tab2, text=df2.iat[15,0], text_color="#000000").grid(row=13, column=0, sticky='e')
            #ctk.CTkLabel(tab2, text=df2.iat[15,1], text_color="#000000").grid(row=13, column=1, sticky='e')
            tk.Label(frame,text='Sobrelape máximo').grid(row=13, column=0, sticky='e')
            tk.Label(frame,text=self.dict_params['maxoverlap']).grid(row=13, column=1, sticky='w', padx=10)

            frame.grid_rowconfigure(14, minsize=10)

            #ctk.CTkLabel(tab2, text=df2.iat[17,0], text_color="#000000").grid(row=15, column=0, sticky='e')
            #ctk.CTkLabel(tab2, text=df2.iat[17,1], text_color="#000000").grid(row=15, column=1, sticky='e')
            tk.Label(frame,text='Delta G hairpin mínimo permitido').grid(row=15, column=0, sticky='e')
            tk.Label(frame,text=self.dict_params['dgmin_hairpin']).grid(row=15, column=1, sticky='w', padx=10)

            #ctk.CTkLabel(tab2, text=df2.iat[18,0], text_color="#000000").grid(row=16, column=0, sticky='e')
            #ctk.CTkLabel(tab2, text=df2.iat[18,1], text_color="#000000").grid(row=16, column=1, sticky='e')
            tk.Label(frame,text='Delta G homodimero mínimo permitido').grid(row=16, column=0, sticky='e')
            tk.Label(frame,text=self.dict_params['dgmin_homodim']).grid(row=16, column=1, sticky='w', padx=10)


            frame.grid_rowconfigure(17, minsize=10)

            #ctk.CTkLabel(tab2, text=df2.iat[20,0], text_color="#000000").grid(row=18, column=0, sticky='e')
            #ctk.CTkLabel(tab2, text=df2.iat[20,1], text_color="#000000").grid(row=18, column=1, sticky='e')
            tk.Label(frame,text='Máximo homopolímeros simples permitidos').grid(row=18, column=0, sticky='e')
            tk.Label(frame,text=self.dict_params['maxhomopol_simple']).grid(row=18, column=1, sticky='w', padx=10)

            #ctk.CTkLabel(tab2, text=df2.iat[21,0], text_color="#000000").grid(row=19, column=0, sticky='e')
            #ctk.CTkLabel(tab2, text=df2.iat[21,1], text_color="#000000").grid(row=19, column=1, sticky='e')
            tk.Label(frame,text='Máximo homopolímeros dobles permitidos').grid(row=19, column=0, sticky='e')
            tk.Label(frame,text=self.dict_params['maxhomopol_double']).grid(row=19, column=1, sticky='w', padx=10)

            #ctk.CTkLabel(tab2, text=df2.iat[22,0], text_color="#000000").grid(row=20, column=0, sticky='e')
            #ctk.CTkLabel(tab2, text=df2.iat[22,1], text_color="#000000").grid(row=20, column=1, sticky='e')
            tk.Label(frame,text='Máximo homopolímeros triples permitidos').grid(row=20, column=0, sticky='e')
            tk.Label(frame,text=self.dict_params['maxhomopol_triple']).grid(row=20, column=1, sticky='w', padx=10)

            #if len(hoja['A25'].value) > 2:
            #if not pd.isna(df2.iat[24,0]):
            if self.dict_params['multiplex']:

                frame.grid_rowconfigure(21, minsize=10)

                #ctk.CTkLabel(tab2, text=df2.iat[24,0], text_color="#000000").grid(row=22, column=0, sticky='e')
                tk.Label(frame,text='Criterios de multiplexación').grid(row=22, column=0, sticky='e')

                frame.grid_rowconfigure(23, minsize=10)

                #ctk.CTkLabel(tab2, text=df2.iat[26,0], text_color="#000000").grid(row=24, column=0, sticky='e')
                #ctk.CTkLabel(tab2, text=df2.iat[26,1], text_color="#000000").grid(row=24, column=1, sticky='e')
                tk.Label(frame,text='Diferencia máxima de Temp Melting').grid(row=24, column=0, sticky='e')
                tk.Label(frame,text=self.dict_params['maxdt']).grid(row=24, column=1, sticky='w', padx=10)

                #ctk.CTkLabel(tab2, text=df2.iat[27,0], text_color="#000000").grid(row=25, column=0, sticky='e')
                #ctk.CTkLabel(tab2, text=df2.iat[27,1], text_color="#000000").grid(row=25, column=1, sticky='e')
                tk.Label(frame,text='Mínimo delta G heterodimerización').grid(row=25, column=0, sticky='e')
                tk.Label(frame,text=self.dict_params['mindg']).grid(row=25, column=1, sticky='w', padx=10)
                

            ###### Fin Tab 2: Parámetros #########

            ###### Inicio Tab 3: Grupos #########

            #hoja = wb_reporte['Resumen']

            #if len(hoja['D1'].value) > 2:
            #if not pd.isna(df3.iat[0,0]):
            if self.dict_params['multiplex']:
                tab3 = ttk.Frame(tabview)
                tabview.add(tab3, text='Grupos')
                canvas3 = tk.Canvas(tab3)
                y_scrollbar = tk.Scrollbar(tab3, orient='vertical', command=canvas3.yview)
                y_scrollbar.pack(side='right', fill='y')
                canvas3.configure(yscrollcommand=y_scrollbar.set)
                canvas3.pack(side='left', fill='both', expand=True)

                # ctk.CTkLabel(tab3, text=df3.iat[0,0], text_color="#000000").grid(row=0, column=0)
                # ctk.CTkLabel(tab3, text=df3.iat[0,1], text_color="#000000").grid(row=0, column=1)
                # ctk.CTkLabel(tab3, text=df3.iat[0,2], text_color="#000000").grid(row=0, column=2)
                # ctk.CTkLabel(tab3, text=df3.iat[0,3], text_color="#000000").grid(row=0, column=3)
                # ctk.CTkLabel(tab3, text=df3.iat[0,4], text_color="#000000").grid(row=0, column=4)

                #ctk.CTkLabel(tab3, text="Grupo", text_color="#000000").grid(row=0, column=0)
                #ctk.CTkLabel(tab3, text="N° sondas", text_color="#000000").grid(row=0, column=1, padx=5)
                #ctk.CTkLabel(tab3, text="Tm", text_color="#000000").grid(row=0, column=2)
                #ctk.CTkLabel(tab3, text="%GC", text_color="#000000").grid(row=0, column=3, padx=5)
                #ctk.CTkLabel(tab3, text="Largo", text_color="#000000").grid(row=0, column=4)
                tk.Label(canvas3,text='Grupo').grid(row=0, column=0, padx=20, pady=20)
                tk.Label(canvas3,text='N° Sondas').grid(row=0, column=1, padx=20, pady=20)
                tk.Label(canvas3,text='Tm').grid(row=0, column=2, padx=20, pady=20)
                tk.Label(canvas3,text='%GC').grid(row=0, column=3, padx=20, pady=20)
                tk.Label(canvas3,text='Largo').grid(row=0, column=4, padx=20, pady=20)

                # continua = True
                # i = 1
                # while continua:
                #     ctk.CTkLabel(tab3, text=df3.iat[i,0], text_color="#000000").grid(row=i, column=0)
                #     ctk.CTkLabel(tab3, text=df3.iat[i,1], text_color="#000000").grid(row=i, column=1)
                #     ctk.CTkLabel(tab3, text=df3.iat[i,2], text_color="#000000").grid(row=i, column=2)
                #     ctk.CTkLabel(tab3, text=df3.iat[i,3], text_color="#000000").grid(row=i, column=3)
                #     ctk.CTkLabel(tab3, text=df3.iat[i,4], text_color="#000000").grid(row=i, column=4)
                #     #if hoja['D'+str(i+1)].value is None:
                #     #if pd.isna(df3.iat[i+1,0]):
                #     if len(df3) < i:
                #         continua = False
                #     i += 1

                for i in range(int(df_resumen['grupo'].max())):
                    grupo = i+1
                    tk.Label(canvas3,text=grupo).grid(row=1+grupo, column=0)
                    tk.Label(canvas3,text=df_resumen[df_resumen['grupo'] == grupo]['grupo'].count()).grid(row=1+grupo, column=1)
                    tk.Label(canvas3,text=str("{:.1f}").format(df_resumen[df_resumen['grupo'] == grupo]['tm'].mean())).grid(row=1+grupo, column=2)
                    tk.Label(canvas3,text=str("{:.1f}").format(df_resumen[df_resumen['grupo'] == grupo]['gc'].mean())).grid(row=1+grupo, column=3)
                    tk.Label(canvas3,text=str("{:.1f}").format(df_resumen[df_resumen['grupo'] == grupo]['largo'].mean())).grid(row=1+grupo, column=4)

                    
            ###### Fin Tab 3: Grupos #########

        
        threading.Thread(target=on_load_data).start()


        def abrir_imagen():
            ventana = tk.Toplevel(root)
            ventana.geometry("700x400")
            Zoom_Advanced(ventana, path=os.path.join(self.folderpath, self.seqrecord.id+'.png'))

        #boton_imagen = ctk.CTkButton(self.frame, text="Visualizar transcritos", corner_radius=30, fg_color="#404040", command=abrir_imagen)
        #boton_imagen.grid(pady=20)

        def guardar_xlsx():
            opciones = {
                'defaultextension': '.xlsx',  # Extensión predeterminada del archivo
                'filetypes': [('Archivos Excel', '.xlsx')],  # Tipos de archivos permitidos
                'initialfile': 'reporte_'+self.seqrecord.id+'.xlsx',  # Nombre predeterminado del archivo
                'title': 'Guardar reporte de diseño',  # Título del diálogo
            }
            ruta_archivo = filedialog.asksaveasfilename(**opciones)

            if ruta_archivo:
                archivo_a_copiar = os.path.join(self.folderpath,self.seqrecord.id+".xlsx")
                shutil.copy(archivo_a_copiar, ruta_archivo)

        def abrir_registro():
            self.ventana = tk.Toplevel(root)
            self.ventana.title("Registrar diseño")
            self.ventana.geometry("")

            descripcion_label = tk.Label(self.ventana, text="Descripción del diseño:")
            descripcion_label.grid(row=1, column=0, padx=10, pady=10)

            def limitar_caracteres(text):
                if len(text) > max_caracteres:
                    return False
                return True

            max_caracteres = 70

            descripcion_entry = ttk.Entry(self.ventana, width=40)
            validation = root.register(limitar_caracteres)
            descripcion_entry.config(validate="key", validatecommand=(validation, "%P"))
            descripcion_entry.grid(row=1, column=1, padx=10, pady=10)
            
            boton_registrar = ctk.CTkButton(self.ventana, text="Registrar diseño", corner_radius=30, fg_color="#404040", command=lambda: registrar_diseno(descripcion_entry.get()))
            boton_registrar.grid(row=2, column=0, padx=20, pady=20)

            boton_cancelar = ctk.CTkButton(self.ventana, text="Cancelar", corner_radius=30, fg_color="#7a7a7a", command=cancelar)
            boton_cancelar.grid(row=2, column=1, padx=20, pady=20)


        def registrar_diseno(descripcion):
            try:
                #conn = sqlite3.connect(os.path.join(os.getcwd(), 'databases', 'probesdb.db'))
                conn = sqlite3.connect(resource_path(os.path.join('databases', 'probesdb.db')))
                cursor = conn.cursor()
                cursor.execute('''CREATE TABLE IF NOT EXISTS sondas (
                                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                                    accnum TEXT,
                                    descripcion TEXT,
                                    secuencia TEXT,
                                    genes TEXT,
                                    fecha_hora DATETIME,
                                    carpeta TEXT)''')
                
                cursor.execute("INSERT INTO sondas (accnum, secuencia, descripcion, genes, fecha_hora, carpeta) VALUES (?,?,?,?,?,?)", (str(self.seqrecord.id), str(self.seqrecord.description), descripcion, str(diseno.get_all_genes(self.seqrecord))[1:-1], datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), self.folderpath))
                conn.commit()
                conn.close()
                self.ventana.destroy()
                messagebox.showinfo("Registro completado", "Se registró correctamente el diseño de las sondas para " + str(diseno.get_all_genes(self.seqrecord))[1:-1])
            except Exception as e:
                messagebox.showinfo("Registro no completado", "No se pudo guardar el diseño de las sondas para " + str(diseno.get_all_genes(self.seqrecord))[1:-1] + ": "+e)

        def cancelar():
            self.ventana.destroy()

        frame_botones = ctk.CTkFrame(master=self.frame, bg_color="transparent")
        frame_botones.grid(row=4, column=0, columnspan=2)

        #boton_volver = tk.Button(self.frame, text="Volver al inicio", command=self.volver_a_inicio)
        ctk.CTkButton(frame_botones, text="Volver al inicio", corner_radius=30, fg_color="#7a7a7a", command=self.volver_a_inicio).grid(row=0, column=0, padx=10)

        ctk.CTkButton(frame_botones, text="Exportar a XLSX", corner_radius=30, fg_color="#7a7a7a", command=guardar_xlsx).grid(row=0, column=1, padx=10)

        ctk.CTkButton(frame_botones, text="Registrar diseño", corner_radius=30, fg_color="#404040", command=abrir_registro).grid(row=0, column=2, padx=10)

        self.root.geometry("800x540")

    
    def volver_a_inicio(self):
        app.mostrar_pantalla('inicial')


class PantallaHistorial(PantallaInicial):
    def __init__(self, root):
        super().__init__(root)

    def crear_interfaz(self):
        # if not os.path.exists(os.path.join(os.getcwd(),'databases')):
        #     os.makedirs(os.path.join(os.getcwd(),'databases'))
        if not os.path.exists(resource_path('databases')):
            os.makedirs(resource_path('databases'))
        #conexion = sqlite3.connect(os.path.join(os.getcwd(),"databases", "probesdb.db"))
        conexion = sqlite3.connect(resource_path(os.path.join("databases", "probesdb.db")))  
        cursor = conexion.cursor()
        cursor.execute('''CREATE TABLE IF NOT EXISTS sondas (
                                id INTEGER PRIMARY KEY AUTOINCREMENT,
                                accnum TEXT,
                                descripcion TEXT,
                                secuencia TEXT,
                                genes TEXT,
                                fecha_hora DATETIME,
                                carpeta TEXT)''')
        cursor.execute("SELECT * FROM sondas ORDER BY id DESC")
        registros = cursor.fetchall()
        #conexion.close()

        genes = []
        for registro in registros:
            if registro[4] not in genes:
                genes.append(registro[4])
        
        #fuente_titulo = font.Font(weight="bold", size=16)
        #historial_label = tk.Label(self.frame, text=f"Historial de paneles de sondas", font=fuente_titulo)
        fuente_titulo = ctk.CTkFont(family='Helvetica', size=20, weight='bold')
        historial_label = ctk.CTkLabel(self.frame, text=f"Historial de paneles de sondas", text_color="#000000", font=fuente_titulo)
        historial_label.grid(row=0, column=0, padx=20, pady=20)


        #tabla_frame = ttk.Frame(self.frame)
        #tabla_frame.grid(row=2,column=0, sticky="nsew")

        #data_table = DataTable(tabla_frame)
        #data_table.pack()
        #data_table = DataTable(self.frame)
        #data_table.grid(row=2,column=0)
        scrollable_frame = ScrollableFrame(self.frame)
        scrollable_frame.grid(row=2, column=0, sticky="nsew")
        data_table = DataTable(scrollable_frame.frame)
        data_table.pack()


        #combobox_filtro = ttk.Combobox(self.frame, values=genes, state="readonly", width=maxlen(genes))
        combobox_filtro = ctk.CTkComboBox(self.frame, values=[""]+genes, width=100, command=data_table.mostrar_gen)
        combobox_filtro.grid(row=1, column=0, padx=20, pady=20)

        #boton_filtro = tk.Button(self.frame, text="Filtrar", command=filtrar_descripcion)
        #boton_filtro.grid(row=1, column=1, padx=20, pady=20)

        #tree = ttk.Treeview(self.frame, columns=("ID", "Genes", "Fecha y hora", "Descripción", "Transcritos", "N° Grupos", "Acciones"), show="headings")
        # tree.heading("ID", text="ID")
        # tree.heading("Genes", text="Genes")
        # tree.heading("Descripción", text="Descripción")
        # tree.heading("Fecha y hora", text="Fecha y hora")        
        # tree.heading("Transcritos", text="Transcritos")
        # tree.heading("N° Grupos", text="N° Grupos")
        # tree.grid(row=2,column=0,padx=10,pady=10)

        

            #row_id = tree.insert("", "end", values=(row[0], row[1], row[2], row[3], row[4], row[5]))

            #button_frame = ttk.Frame(tree)
            #button = ttk.Button(button_frame, text="Editar", command=lambda row_id=row_id: mostrar_id(row_id))
            #button.pack(fill="both", expand=True)
            #tree.attach_widget(row_id, "#0", button_frame)

            #button_label = ttk.Button(tree, text="Ver", command=lambda row_id=row[0]: mostrar_id(row_id))
            #button_label.bind("<Button-1>", lambda event: mostrar_id(row[0]))
            #tree.window_create("", window=button_label, anchor="w", tags=("button",))
            #tree.insert(index, "end", values=(None, None, None, None, None, None, button_label))

            #button_label2 = ttk.Button(tree, text="Eliminar", command=lambda row_id=row[0]: eliminar_registro(row_id))
            #button_label2.bind("<Button-1>", lambda event: eliminar_registro(row[0]))
            #tree.window_create("", window=button_label, anchor="w", tags=("button",))
            #tree.insert(index, "end", values=(None, None, None, None, None, None, None, button_label2))

            #tree.insert("", "end", values=(row[0], row[1], row[2], row[3], row[4], row[5], button_label, button_label2))

            #btn_mostrar = tk.Button(self.frame, text="Ver", command=lambda id=row[0]: mostrar_id(id))
            #btn_mostrar.grid(row=tree.index(index), column=3)

            #btn_eliminar = tk.Button(self.frame, text="X", command=lambda id_registro=row[0]: eliminar_registro(id_registro))
            #btn_eliminar.grid(row=tree.index(index), column=4)

        #for col in tree["columns"]:
            #tree.column(col, width=120, stretch=False)


        # #texto_registros = tk.Text(self.frame, wrap=tk.WORD, width=100, height=20)
        # texto_registros = ctk.CTkTextbox(self.frame, width=600, text_color="#000000")
        # #texto_registros.grid(row=2, column=0, padx=20, pady=20)
    
        # for registro in registros:
        #     texto_registros.insert("end", f"ID: {registro[0]}\n")
        #     texto_registros.insert("end", f"Descripción: {registro[1]}\n")
        #     texto_registros.insert("end", f"Secuencia: {registro[2]}\n")
        #     texto_registros.insert("end", f"Genes: {registro[3]}\n")
        #     texto_registros.insert("end", f"Fecha y Hora: {registro[4]}\n")
        #     texto_registros.insert("end", f"Carpeta: {registro[5]}\n\n")
        
        # texto_registros.configure(state="disabled")

        #boton_volver = tk.Button(self.frame, text="Volver", command=self.volver_a_seleccion)
        boton_volver = ctk.CTkButton(self.frame, text="Volver", corner_radius=30, fg_color="#404040", command=self.volver_a_seleccion)
        boton_volver.grid(row=4, column=0,columnspan=2, pady=20)

        self.root.geometry("800x600")

    
    def volver_a_seleccion(self):
        app.mostrar_pantalla("inicial")

class PantallaRegistro(PantallaInicial):
    def __init__(self, root, folderpath, genes, accnum):
        self.folderpath = folderpath
        self.genes = genes
        self.accnum = accnum
        super().__init__(root)

    def crear_interfaz(self):
        fuente_negrita = ctk.CTkFont(family='Helvetica', weight='bold')
        resultado_label = ctk.CTkLabel(self.frame, text=f"Sondas diseñadas para {self.genes}", text_color="#000000", font=fuente_negrita)
        resultado_label.grid(row=0, column=0, columnspan=2, padx=20, pady=20)

        carpeta_label = ctk.CTkLabel(self.frame, text=f"El reporte, la imagen y el archivo GenBank se almacenaron en\n" + self.folderpath, text_color="#000000")
        carpeta_label.grid(row=1, column=0, columnspan=2, padx=20, pady=20)

        img_frame = ttk.Frame(self.frame, width=1300, height=700)
        img_frame.grid(row=2, column=0, padx=10, pady=20)
        Zoom_Advanced(mainframe=img_frame, path=os.path.join(self.folderpath, self.accnum+'.png'))

        tabview = ttk.Notebook(self.frame, width=400, height=400)
        tab1 = ttk.Frame(tabview)
        tab2 = ttk.Frame(tabview)
        tabview.add(tab1, text='Resumen')
        tabview.add(tab2, text='Parámetros')
        tabview.grid(row=2, column=1, padx=10, pady=20)

        xls = pd.ExcelFile(os.path.join(self.folderpath, self.accnum+'.xlsx'))
        hojas_excel = {}
        for nombre_hoja in xls.sheet_names:
            hojas_excel[nombre_hoja] = pd.read_excel(os.path.join(self.folderpath, self.accnum+'.xlsx'), sheet_name=nombre_hoja)

        #wb_reporte = openpyxl.load_workbook(os.path.join(self.folderpath, self.accnum+'.xlsx'), read_only=True)
        #hoja = wb_reporte['Resumen']
            
        def on_load_data():
            canvas1 = tk.Canvas(tab1)
            #y_scrollbar = tk.Scrollbar(tab1, orient='vertical', command=canvas1.yview)
            #y_scrollbar.pack(side='right', fill='y')
            #canvas1.configure(yscrollcommand=y_scrollbar.set)
            canvas1.grid(row=0,column=0)

            ###### Inicio Tab 1: Resumen #########

            a=[1,2,3,5,6,7,9]

            tk.Label(canvas1,text=hojas_excel["Resumen"].columns[0]).grid(row=0, column=0, sticky='e')

            for i in a:
                tk.Label(canvas1,text=hojas_excel["Resumen"].iat[i,0]).grid(row=i+1, column=0, sticky='e')
                tk.Label(canvas1,text=hojas_excel["Resumen"].iat[i,1]).grid(row=i+1, column=1, sticky='w', padx=10)

            canvas1.grid_rowconfigure(1, minsize=10)
            canvas1.grid_rowconfigure(5, minsize=10)
            canvas1.grid_rowconfigure(9, minsize=10)

            multiplex = not hojas_excel["Grupos"].empty
            
            if multiplex:
                canvas1.grid_rowconfigure(11, minsize=10)
                tk.Label(canvas1,text=hojas_excel["Resumen"].iat[11,0]).grid(row=12, column=0, sticky='e')
                tk.Label(canvas1,text=hojas_excel["Resumen"].iat[11,1]).grid(row=12, column=1, sticky='w', padx=10)

            ###### Inicio Tab 2: Parámetros #########
                
            
            #y_scrollbar = tk.Scrollbar(tab2, orient='vertical', command=canvas2.yview)
            #y_scrollbar.grid(row=0,column=1)
            
            def on_configure(event):
                canvas2.configure(scrollregion=canvas2.bbox('all'))
                canvas_height = frame.winfo_reqheight()
                if canvas_height != canvas2.winfo_height():
                    canvas2.config(height=canvas_height)
                    scrollbar.config(command=canvas2.yview)
            
            scrollbar = tk.Scrollbar(tab2)
            scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
            canvas2 = tk.Canvas(tab2)
            canvas2.pack(side='left', fill='both', expand=True)
            #canvas2.grid(row=0,column=0)
            scrollbar.config(command=canvas2.yview)

            def scroll(event):
                canvas2.yview_scroll(-1*(event.delta//120), "units")
            canvas2.bind_all("<MouseWheel>", scroll)

            frame = tk.Frame(canvas2)
            canvas2.create_window((0, 0), window=frame, anchor='nw')


            frame.bind("<Configure>", on_configure)

            a=[1,2,4,5,7,8,10,11,13,14,16,17,19,20,21]

            for i in a:
                tk.Label(frame,text=hojas_excel["Parámetros iniciales"].iat[i,0]).grid(row=i+1, column=0, sticky='e')
                tk.Label(frame,text=hojas_excel["Parámetros iniciales"].iat[i,1]).grid(row=i+1, column=1, sticky='w', padx=10)

            frame.grid_rowconfigure(4, minsize=10)
            frame.grid_rowconfigure(7, minsize=10)
            frame.grid_rowconfigure(10, minsize=10)
            frame.grid_rowconfigure(13, minsize=10)
            frame.grid_rowconfigure(16, minsize=10)
            frame.grid_rowconfigure(19, minsize=10)
            frame.grid_rowconfigure(23, minsize=10)
            
            if multiplex:
                frame.grid_rowconfigure(25, minsize=10)
                tk.Label(frame,text=hojas_excel["Parámetros iniciales"].iat[23,0]).grid(row=24, column=0, sticky='e')
                
                tk.Label(frame,text=hojas_excel["Parámetros iniciales"].iat[25,0]).grid(row=26, column=0, sticky='w', padx=10)
                tk.Label(frame,text=hojas_excel["Parámetros iniciales"].iat[25,1]).grid(row=26, column=1, sticky='w', padx=10)

                tk.Label(frame,text=hojas_excel["Parámetros iniciales"].iat[26,0]).grid(row=27, column=0, sticky='w', padx=10)
                tk.Label(frame,text=hojas_excel["Parámetros iniciales"].iat[26,1]).grid(row=27, column=1, sticky='w', padx=10)

            ###### Inicio Tab 3: Grupos #########
            if multiplex:
                tab3 = ttk.Frame(tabview)
                tabview.add(tab3, text='Grupos')
                canvas3 = tk.Canvas(tab3)
                #y_scrollbar = tk.Scrollbar(tab3, orient='vertical', command=canvas3.yview)
                #y_scrollbar.pack(side='right', fill='y')
                #canvas3.configure(yscrollcommand=y_scrollbar.set)
                #canvas3.grid(row=0,column=0,side='left', fill='both', expand=True)
                canvas3.grid(row=0,column=0)

                tk.Label(canvas3,text='Grupo').grid(row=0, column=0, padx=20, pady=20)
                tk.Label(canvas3,text='N° Sondas').grid(row=0, column=1, padx=20, pady=20)
                tk.Label(canvas3,text='Tm').grid(row=0, column=2, padx=20, pady=20)
                tk.Label(canvas3,text='%GC').grid(row=0, column=3, padx=20, pady=20)
                tk.Label(canvas3,text='Largo').grid(row=0, column=4, padx=20, pady=20)

                for i in range(int(hojas_excel["Grupos"]['GRUPO'].max())):
                    grupo = i+1
                    tk.Label(canvas3,text=hojas_excel["Grupos"].iat[i,0]).grid(row=1+grupo, column=0)
                    tk.Label(canvas3,text=hojas_excel["Grupos"].iat[i,1]).grid(row=1+grupo, column=1)
                    tk.Label(canvas3,text=hojas_excel["Grupos"].iat[i,2]).grid(row=1+grupo, column=2)
                    tk.Label(canvas3,text=hojas_excel["Grupos"].iat[i,3]).grid(row=1+grupo, column=3)
                    tk.Label(canvas3,text=hojas_excel["Grupos"].iat[i,4]).grid(row=1+grupo, column=4)
            
        def guardar_xlsx():
            opciones = {
                'defaultextension': '.xlsx',  # Extensión predeterminada del archivo
                'filetypes': [('Archivos Excel', '.xlsx')],  # Tipos de archivos permitidos
                'initialfile': 'reporte_'+self.accnum+'.xlsx',  # Nombre predeterminado del archivo
                'title': 'Guardar reporte de diseño',  # Título del diálogo
            }
            ruta_archivo = filedialog.asksaveasfilename(**opciones)

            if ruta_archivo:
                archivo_a_copiar = os.path.join(self.folderpath,self.accnum+".xlsx")
                shutil.copy(archivo_a_copiar, ruta_archivo)
        
        self.root.geometry("800x540")
        threading.Thread(target=on_load_data).start()

        def abrir_carpeta():
            sistema_operativo = platform.system()
            if sistema_operativo == 'Windows':
                os.system(f'explorer "{self.folderpath}"')
            elif sistema_operativo == 'Darwin':  # macOS
                os.system(f'open "{self.folderpath}"')
            elif sistema_operativo == 'Linux':
                os.system(f'xdg-open "{self.folderpath}"')

        frame_botones = ctk.CTkFrame(master=self.frame, bg_color="transparent")
        frame_botones.grid(row=3, column=0, columnspan=2)
        
        ctk.CTkButton(frame_botones, text="Volver a historial", corner_radius=30, fg_color="#7a7a7a", command=self.volver_a_inicio).grid(row=0, column=0, padx=10)

        ctk.CTkButton(frame_botones, text="Exportar a XLSX", corner_radius=30, fg_color="#7a7a7a", command=guardar_xlsx).grid(row=0, column=1, padx=10)

        ctk.CTkButton(frame_botones, text="Abrir carpeta", corner_radius=30, fg_color="#404040", command=abrir_carpeta).grid(row=0, column=2, padx=10)


        

    def volver_a_inicio(self):
            app.mostrar_pantalla('historial')

        

class ControladorApp:
    def __init__(self, root):
        self.root = root
        self.pantallas = {
            "inicial": PantallaInicial(root),
            "transcripciones": None,
            "parametros": None,
            "carga": None,
            "historial": None,
            "final": None,
            "registro": None
        }
        self.mostrar_pantalla("inicial")

    def mostrar_pantalla(self, nombre, seqrecord=None, transcripciones=None, dict_params=None, df=None, filepath=None, folderpath=None, df1=None, df2=None, df3=None, genes=None, accnum=None):
        # Ocultar pantalla actual
        if hasattr(self, "pantalla_actual"):
            self.pantalla_actual.frame.pack_forget()
        # Mostrar la pantalla solicitada
        if nombre == "transcripciones":
            self.pantalla_actual = PantallaTranscripciones(self.root, seqrecord, filepath)
            self.pantalla_actual.frame.pack(padx=50, pady=50)
        elif nombre == "parametros":
            self.pantalla_actual = PantallaParametros(self.root, seqrecord, transcripciones, filepath)
            self.pantalla_actual.frame.pack(padx=10, pady=10)
        elif nombre == "carga":
            #self.carga = tk.Toplevel(self.root)
            self.pantalla_actual = PantallaCarga(self.root, seqrecord, transcripciones, dict_params, filepath)
            self.pantalla_actual.frame.pack(padx=50, pady=50)
            #self.root.after(200, lambda: self.disenar_sondas(seqrecord, transcripciones, dict_params, filepath))
        elif nombre == "final":
            self.pantalla_actual = PantallaFinal(self.root, seqrecord, transcripciones, df, dict_params, filepath, folderpath, df1, df2, df3)
            self.pantalla_actual.frame.pack(padx=50, pady=50)
        elif nombre == "historial":
            self.pantalla_actual = PantallaHistorial(self.root)
            self.pantalla_actual.frame.pack(padx=50, pady=50)
            #self.pantalla_actual.frame.grid(row=0,column=0,padx=50, pady=50)
        elif nombre == "registro":
            self.pantalla_actual = PantallaRegistro(self.root, folderpath, genes, accnum)
            self.pantalla_actual.frame.pack(padx=50, pady=50)
        else:
            self.pantalla_actual = self.pantallas[nombre]
            self.pantalla_actual.frame.pack(padx=50, pady=50)
            self.root.geometry("800x450")      


class ToolTip:
    def __init__(self, widget, text, image_path=None):
        self.widget = widget
        self.text = text
        self.image_path = image_path
        self.tooltip = None
        self.widget.bind("<Enter>", self.mostrar_tooltip)
        self.widget.bind("<Leave>", self.ocultar_tooltip)

    def mostrar_tooltip(self, event):
        x, y, _, _ = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 25

        self.tooltip = tk.Toplevel(self.widget)
        self.tooltip.wm_overrideredirect(True)
        self.tooltip.wm_geometry(f"+{x}+{y}")

        if self.image_path:
            imagen = Image.open(self.image_path)
            imagen = ImageTk.PhotoImage(imagen)
            label_imagen = tk.Label(self.tooltip, image=imagen)
            label_imagen.image = imagen
            label_imagen.pack()

        label_texto = tk.Label(self.tooltip, text=self.text, background="lightyellow", relief="solid")
        label_texto.pack()

    def ocultar_tooltip(self, event):
        if self.tooltip:
            self.tooltip.destroy()
            self.tooltip = None


class IntegerSelector(ctk.CTkFrame):
    def __init__(self, master, default_value=0, min_value=None, max_value=None, increment=1, **kwargs):
        super().__init__(master, **kwargs)
        
        self.min_value = min_value
        self.max_value = max_value
        self.increment = increment
        
        self.valor = tk.IntVar()
        self.valor.set(default_value)
        
        #self.btn_disminuir = tk.Button(self, text="-", command=self.disminuir)
        self.btn_disminuir = ctk.CTkButton(master=self, text="-", fg_color="#404040", corner_radius=30, width=10, height=10,command=self.disminuir)
        self.btn_disminuir.pack(side="left")

        self.entry = tk.Entry(self, textvariable=self.valor, validate="key", validatecommand=(self.register(self.validar), "%P"), width=15, justify='center')
        #self.entry = customtkinter.CTkEntry(self, placeholder_text="CTkEntry")
        self.entry.pack(side="left", padx=5)
        
        #self.btn_aumentar = tk.Button(self, text="+", command=self.aumentar)
        self.btn_aumentar = ctk.CTkButton(master=self, text="+", fg_color="#404040", corner_radius=30, width=15, height=10,command=self.aumentar)
        self.btn_aumentar.pack(side="left")

        self.entry.bind("<Up>", lambda event: self.aumentar())
        self.entry.bind("<Down>", lambda event: self.disminuir())
        self.entry.bind("<Right>", lambda event: self.aumentar())
        self.entry.bind("<Left>", lambda event: self.disminuir())
    
    def validar(self, nuevo_valor):
        if nuevo_valor == "" or nuevo_valor.isdigit():
            return True
        else:
            return False

    def disminuir(self):
        nuevo_valor = self.valor.get() - self.increment
        if self.min_value is not None:
            nuevo_valor = max(nuevo_valor, self.min_value)
        self.valor.set(nuevo_valor)

    def aumentar(self):
        nuevo_valor = self.valor.get() + self.increment
        if self.max_value is not None:
            nuevo_valor = min(nuevo_valor, self.max_value)
        self.valor.set(nuevo_valor)


class AutoScrollbar(ttk.Scrollbar):
    ''' A scrollbar that hides itself if it's not needed.
        Works only if you use the grid geometry manager '''
    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            self.grid_remove()
        else:
            self.grid()
            ttk.Scrollbar.set(self, lo, hi)

    def pack(self, **kw):
        raise tk.TclError('Cannot use pack with this widget')

    def place(self, **kw):
        raise tk.TclError('Cannot use place with this widget')


#https://stackoverflow.com/questions/41656176/tkinter-canvas-zoom-move-pan
class Zoom_Advanced(ttk.Frame):
    ''' Advanced zoom of the image '''
    def __init__(self, mainframe, path):
        ''' Initialize the main Frame '''
        ttk.Frame.__init__(self, master=mainframe)
        #self.master.title('Visualizador de transcritos')
        # Vertical and horizontal scrollbars for canvas
        vbar = AutoScrollbar(self.master, orient='vertical')
        hbar = AutoScrollbar(self.master, orient='horizontal')
        vbar.grid(row=0, column=1, sticky='ns')
        hbar.grid(row=1, column=0, sticky='we')
        # Create canvas and put image on it
        self.canvas = tk.Canvas(self.master, highlightthickness=0,
                                xscrollcommand=hbar.set, yscrollcommand=vbar.set)
        self.canvas.grid(row=0, column=0, sticky='nswe')
        self.canvas.update()  # wait till canvas is created
        vbar.configure(command=self.scroll_y)  # bind scrollbars to the canvas
        hbar.configure(command=self.scroll_x)
        # Make the canvas expandable
        self.master.rowconfigure(0, weight=1)
        self.master.columnconfigure(0, weight=1)
        # Bind events to the Canvas
        self.canvas.bind('<Configure>', self.show_image)  # canvas is resized
        self.canvas.bind('<ButtonPress-1>', self.move_from)
        self.canvas.bind('<B1-Motion>',     self.move_to)
        self.canvas.bind('<MouseWheel>', self.wheel)  # with Windows and MacOS, but not Linux
        self.canvas.bind('<Button-5>',   self.wheel)  # only with Linux, wheel scroll down
        self.canvas.bind('<Button-4>',   self.wheel)  # only with Linux, wheel scroll up
        self.image = Image.open(path)  # open image
        self.width, self.height = self.image.size
        self.imscale = 1.0  # scale for the canvaas image
        self.delta = 1.3  # zoom magnitude
        # Put image into container rectangle and use it to set proper coordinates to the image
        self.container = self.canvas.create_rectangle(0, 0, self.width, self.height, width=0)
        self.show_image()

    def scroll_y(self, *args, **kwargs):
        ''' Scroll canvas vertically and redraw the image '''
        self.canvas.yview(*args, **kwargs)  # scroll vertically
        self.show_image()  # redraw the image

    def scroll_x(self, *args, **kwargs):
        ''' Scroll canvas horizontally and redraw the image '''
        self.canvas.xview(*args, **kwargs)  # scroll horizontally
        self.show_image()  # redraw the image

    def move_from(self, event):
        ''' Remember previous coordinates for scrolling with the mouse '''
        self.canvas.scan_mark(event.x, event.y)

    def move_to(self, event):
        ''' Drag (move) canvas to the new position '''
        self.canvas.scan_dragto(event.x, event.y, gain=1)
        self.show_image()  # redraw the image

    def wheel(self, event):
        ''' Zoom with mouse wheel '''
        x = self.canvas.canvasx(event.x)
        y = self.canvas.canvasy(event.y)
        bbox = self.canvas.bbox(self.container)  # get image area
        if bbox[0] < x < bbox[2] and bbox[1] < y < bbox[3]: pass  # Ok! Inside the image
        else: return  # zoom only inside image area
        scale = 1.0
        # Respond to Linux (event.num) or Windows (event.delta) wheel event
        if event.num == 5 or event.delta == -120:  # scroll down
            i = min(self.width, self.height)
            if int(i * self.imscale) < 30: return  # image is less than 30 pixels
            self.imscale /= self.delta
            scale        /= self.delta
        if event.num == 4 or event.delta == 120:  # scroll up
            i = min(self.canvas.winfo_width(), self.canvas.winfo_height())
            if i < self.imscale: return  # 1 pixel is bigger than the visible area
            self.imscale *= self.delta
            scale        *= self.delta
        self.canvas.scale('all', x, y, scale, scale)  # rescale all canvas objects
        self.show_image()

    def show_image(self, event=None):
        ''' Show image on the Canvas '''
        bbox1 = self.canvas.bbox(self.container)  # get image area
        # Remove 1 pixel shift at the sides of the bbox1
        bbox1 = (bbox1[0] + 1, bbox1[1] + 1, bbox1[2] - 1, bbox1[3] - 1)
        bbox2 = (self.canvas.canvasx(0),  # get visible area of the canvas
                 self.canvas.canvasy(0),
                 self.canvas.canvasx(self.canvas.winfo_width()),
                 self.canvas.canvasy(self.canvas.winfo_height()))
        bbox = [min(bbox1[0], bbox2[0]), min(bbox1[1], bbox2[1]),  # get scroll region box
                max(bbox1[2], bbox2[2]), max(bbox1[3], bbox2[3])]
        if bbox[0] == bbox2[0] and bbox[2] == bbox2[2]:  # whole image in the visible area
            bbox[0] = bbox1[0]
            bbox[2] = bbox1[2]
        if bbox[1] == bbox2[1] and bbox[3] == bbox2[3]:  # whole image in the visible area
            bbox[1] = bbox1[1]
            bbox[3] = bbox1[3]
        self.canvas.configure(scrollregion=bbox)  # set scroll region
        x1 = max(bbox2[0] - bbox1[0], 0)  # get coordinates (x1,y1,x2,y2) of the image tile
        y1 = max(bbox2[1] - bbox1[1], 0)
        x2 = min(bbox2[2], bbox1[2]) - bbox1[0]
        y2 = min(bbox2[3], bbox1[3]) - bbox1[1]
        if int(x2 - x1) > 0 and int(y2 - y1) > 0:  # show image if it in the visible area
            x = min(int(x2 / self.imscale), self.width)   # sometimes it is larger on 1 pixel...
            y = min(int(y2 / self.imscale), self.height)  # ...and sometimes not
            image = self.image.crop((int(x1 / self.imscale), int(y1 / self.imscale), x, y))
            imagetk = ImageTk.PhotoImage(image.resize((int(x2 - x1), int(y2 - y1))))
            imageid = self.canvas.create_image(max(bbox2[0], bbox1[0]), max(bbox2[1], bbox1[1]),
                                               anchor='nw', image=imagetk)
            self.canvas.lower(imageid)  # set image into background
            self.canvas.imagetk = imagetk  # keep an extra reference to prevent garbage-collection

class ScrollableFrame(tk.Frame):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        # Configurar el canvas para hacer el frame scrollable
        self.canvas = tk.Canvas(self)
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.scrollbar = ttk.Scrollbar(self, orient=tk.VERTICAL, command=self.canvas.yview)
        self.scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.canvas.configure(yscrollcommand=self.scrollbar.set)

        # Crear un frame dentro del canvas para colocar el contenido
        self.frame = tk.Frame(self.canvas)
        self.canvas.create_window((0, 0), window=self.frame, anchor='nw')

        # Asegurar que el canvas se ajuste al contenido cuando cambie el tamaño
        self.frame.bind("<Configure>", self.resize_canvas)
        self.canvas.bind_all("<MouseWheel>", self._on_mousewheel)
    
    def _on_mousewheel(self, event):
        self.canvas.yview_scroll(int(-1*(event.delta/120)), "units")

        # Evitar el desplazamiento de la ventana principal en Windows
        return "break"

    def resize_canvas(self, event):
        self.canvas.configure(scrollregion=self.canvas.bbox('all'))
        canvas_height = min(self.frame.winfo_reqheight(), 600)  # Establecer el alto máximo
        self.canvas.config(height=canvas_height,width=self.frame.winfo_reqwidth())

class DataTable(tk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self.table_frames = []
        self.create_table()

    def create_table(self):

        #self.conexion = sqlite3.connect(os.path.join(os.getcwd(),"databases", "probesdb.db"))
        self.conexion = sqlite3.connect(resource_path(os.path.join("databases", "probesdb.db")))
        
        self.cursor = self.conexion.cursor()
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS sondas (
                                id INTEGER PRIMARY KEY AUTOINCREMENT,
                                accnum TEXT,
                                descripcion TEXT,
                                secuencia TEXT,
                                genes TEXT,
                                fecha_hora DATETIME,
                                carpeta TEXT)''')
        data=[]
        data.append(['ID', 'Descripción', 'Accession Num', 'Genes', 'Fecha y hora'])
        
        self.cursor.execute("SELECT id, accnum, descripcion, secuencia, genes, fecha_hora FROM sondas ORDER BY id DESC")
        for row in self.cursor.fetchall():
            data.append([row[0], row[2], row[1], row[4], row[5]])

        self.data = data

        # Destruir todos los widgets de la tabla actual
        for frame in self.table_frames:
            frame.destroy()
        
        # Limpiar la lista de frames
        self.table_frames.clear()

        # Crear una nueva tabla con los datos actualizados
        column_widths = []
        for col in zip(*self.data):
            column_widths.append(max(len(str(cell))-1 for cell in col))

        for row_index, row_data in enumerate(self.data):
            row_frame = tk.Frame(self)
            row_frame.grid(row=row_index, column=0, sticky="ew")
            self.table_frames.append(row_frame)

            for col_index, cell_data in enumerate(row_data):
                cell_label = tk.Label(row_frame, text=str(cell_data), padx=5, pady=2, borderwidth=1, relief="solid", width=column_widths[col_index])
                cell_label.grid(row=0, column=col_index, sticky="ew")

            if row_index > 0 :
                ver_button = tk.Button(row_frame, text="Ver", command=lambda row=row_index: self.on_button_click_1(row))
                ver_button.grid(row=0, column=len(row_data), padx=5, pady=2, sticky="e")

                delete_button = tk.Button(row_frame, text="Eliminar", command=lambda row=row_index: self.on_button_click_2(row))
                delete_button.grid(row=0, column=len(row_data)+1, padx=5, pady=2, sticky="e")

            row_frame.columnconfigure(len(row_data), weight=1)

    def on_button_click_1(self, row):
        id_registro = self.data[row][0]
        self.cursor.execute("SELECT accnum, carpeta, genes FROM sondas WHERE id=?", (id_registro,))
        registro = self.cursor.fetchone()
        #print(registro[0])
        #print(registro[1])
        #print(registro[2])
        #print(f"Botón clickeado en la fila {row}")
        app.mostrar_pantalla("registro", folderpath=registro[1], genes=registro[2], accnum=registro[0])


    def on_button_click_2(self, row):
        id_registro = self.data[row][0]
        self.cursor.execute("SELECT carpeta FROM sondas WHERE id=?", (id_registro,))
        carpeta = self.cursor.fetchone()
        respuesta = messagebox.askyesno("Confirmación", "¿Estás seguro(a) que quieres eliminar este diseño? Se eliminará la carpeta "+carpeta[0])
        if respuesta:
            try:
                self.cursor.execute("DELETE FROM sondas WHERE carpeta=?", (carpeta[0],))
                self.conexion.commit()
                shutil.rmtree(carpeta[0])
            except OSError as error:
                print("No se pudo eliminar la carpeta: {error}")
            # Eliminar la fila correspondiente de los datos
            del self.data[row]
            # Actualizar la tabla con los datos actualizados
            self.create_table()

    def mostrar_gen(self, gen):
        if len(gen)>1:
            data=[]
            data.append(['ID', 'Accession Num', 'Descripción', 'Genes', 'Fecha y hora'])
            
            self.cursor.execute("SELECT id, accnum, descripcion, secuencia, genes, fecha_hora FROM sondas WHERE genes = (?) ORDER BY id DESC",(gen,))
            for row in self.cursor.fetchall():
                data.append([row[0], row[1], row[2], row[4], row[5]])

            self.data = data

            # Destruir todos los widgets de la tabla actual
            for frame in self.table_frames:
                frame.destroy()
            
            # Limpiar la lista de frames
            self.table_frames.clear()

            # Crear una nueva tabla con los datos actualizados
            column_widths = []
            for col in zip(*self.data):
                column_widths.append(max(len(str(cell)) for cell in col))

            for row_index, row_data in enumerate(self.data):
                row_frame = tk.Frame(self)
                row_frame.grid(row=row_index, column=0, sticky="ew")
                self.table_frames.append(row_frame)

                for col_index, cell_data in enumerate(row_data):
                    cell_label = tk.Label(row_frame, text=str(cell_data), padx=2, pady=2, borderwidth=1, relief="solid", width=column_widths[col_index])
                    cell_label.grid(row=0, column=col_index, sticky="ew")

                if row_index > 0 :
                    ver_button = tk.Button(row_frame, text="Ver", command=lambda row=row_index: self.on_button_click_1(row))
                    ver_button.grid(row=0, column=len(row_data), padx=5, pady=2, sticky="e")

                    delete_button = tk.Button(row_frame, text="Eliminar", command=lambda row=row_index: self.on_button_click_2(row))
                    delete_button.grid(row=0, column=len(row_data)+1, padx=5, pady=2, sticky="e")

                row_frame.columnconfigure(len(row_data), weight=1)


if __name__ == "__main__":
    ctk.set_appearance_mode("light")
    #root = tk.Tk()
    root = ctk.CTk()
    root.state("zoomed")
    app = ControladorApp(root)
    root.iconbitmap(resource_path(os.path.join("images", "gemini2.ico")))
    root.wm_iconbitmap(resource_path(os.path.join("images", "gemini2.ico")))
    #root.tk.call('wm', 'iconphoto', root._w, tk.PhotoImage(file=os.path.join('images','gemini2.png')))
    root.mainloop()
