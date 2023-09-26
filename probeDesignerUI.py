import descarga
import diseno
import imagen
import tkinter as tk
import customtkinter
import os
import sqlite3
import datetime
import platform
import shutil
import threading
import queue
from diseno import progreso_compartido
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from PIL import Image, ImageTk

class PantallaInicial:
    def __init__(self, root):
        self.root = root
        self.root.title("GEMINi - Diseñador de sondas de hibridación para ARN mensajero")
        #self.root.iconbitmap(os.path.join("images", "gemini2.ico"))
        #self.frame = tk.Frame(root)
        self.frame = customtkinter.CTkFrame(master=root)
        self.frame.pack(padx=50, pady=50)
        self.crear_interfaz()
        #self.root.minsize(800,450)
        #self.root.geometry(f"{self.frame.winfo_reqwidth()+20}x{self.frame.winfo_reqheight()+20}")

    def crear_interfaz(self):
        # Crear widgets de la pantalla de inicio
        #fuente_titulo = font.Font(weight="bold", size=16)
        #titulo_archivo = tk.Label(self.frame, text="Seleccione un archivo GenBank", font=fuente_titulo)
        fuente_titulo = customtkinter.CTkFont(family='Helvetica', size=20, weight='bold')
        titulo_archivo = customtkinter.CTkLabel(self.frame, text="Seleccione un archivo GenBank\n(.gb o .gbk)", text_color="#BE272F", font=fuente_titulo)
        titulo_archivo.pack(padx=100, pady=30)

        # Pantalla de selección de archivo
        #seleccion_button = tk.Button(self.frame, text="Seleccionar Archivo", command=self.seleccionar_archivo)
        seleccion_button = customtkinter.CTkButton(self.frame, text="Seleccionar Archivo", corner_radius=30, fg_color="#BE272F", command=self.seleccionar_archivo)
        seleccion_button.pack(pady=30)

        #historial_button = tk.Button(self.frame, text="Ver historial de sondas", command=self.ver_historial)
        historial_button = customtkinter.CTkButton(self.frame, text="Ver historial de sondas", corner_radius=30, fg_color="#E89434",command=self.ver_historial)
        historial_button.pack(pady=60)
        
        self.root.geometry("800x450")
        
    def seleccionar_archivo(self):
        archivo = filedialog.askopenfilename(filetypes=[("Archivos GenBank", "*.gbk *.gb")])
        if archivo:
            seqrecord = descarga.parse_file_to_seqrecord(archivo)
            app.mostrar_pantalla("transcripciones", seqrecord, filepath=archivo)
    
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
        fuente_titulo = customtkinter.CTkFont(family='Helvetica', size=20, weight='bold')
        transcripciones_label = customtkinter.CTkLabel(self.frame, text="Seleccione las transcripciones para el diseño de sondas", text_color="#BE272F", font=fuente_titulo)
        transcripciones_label.grid(row=0, column=0, columnspan=3, padx=10, pady=20)

        label_izq = customtkinter.CTkLabel(self.frame, text="Transcripciones no seleccionadas", text_color="#BE272F")
        label_izq.grid(row=1, column=0, padx=10, pady=10)

        label_der = customtkinter.CTkLabel(self.frame, text="Transcripciones seleccionadas", text_color="#BE272F")
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
        
        self.lista_disponibles = tk.Listbox(self.frame, selectmode=tk.SINGLE, selectbackground="#BE272F", width=maxlen(self.opciones_disponibles))
        self.lista_disponibles.grid(row=2, column= 0, padx=10, pady=10)

        for opcion in self.opciones_disponibles:
            self.lista_disponibles.insert(tk.END, opcion)

        self.lista_seleccionadas = tk.Listbox(self.frame, selectmode=tk.SINGLE, selectbackground="#BE272F", width=maxlen(self.opciones_disponibles))
        self.lista_seleccionadas.grid(row=2, column= 2, padx=10, pady=10)

        botones_frame = customtkinter.CTkFrame(master=self.frame, bg_color="transparent")
        botones_frame.grid(row=2, column=1)

        #self.boton_seleccionar = tk.Button(botones_frame, text="Seleccionar", command=self.seleccionar)
        self.boton_seleccionar = customtkinter.CTkButton(botones_frame, text=">", corner_radius=30, fg_color="#BE272F", command=self.seleccionar)
        self.boton_seleccionar.pack(pady=10)
        
        #self.boton_eliminar = tk.Button(botones_frame, text="Eliminar", command=self.eliminar)
        self.boton_eliminar = customtkinter.CTkButton(botones_frame, text="<", corner_radius=30, fg_color="#E89434", command=self.eliminar)
        self.boton_eliminar.pack(pady=10)

        #self.boton_seleccionar_todo = tk.Button(botones_frame, text="Seleccionar Todo", command=self.seleccionar_todo)
        self.boton_seleccionar_todo = customtkinter.CTkButton(botones_frame, text=">>>>", corner_radius=30, fg_color="#BE272F", command=self.seleccionar_todo)
        self.boton_seleccionar_todo.pack(pady=10)

        #self.boton_eliminar_todo = tk.Button(botones_frame, text="Eliminar Todo", command=self.eliminar_todo)
        self.boton_eliminar_todo = customtkinter.CTkButton(botones_frame, text="<<<<", corner_radius=30, fg_color="#E89434", command=self.eliminar_todo)
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
        boton_volver = customtkinter.CTkButton(self.frame, text="Volver", corner_radius=30, fg_color="#E89434", command=self.volver_a_seleccion)
        boton_volver.grid(row=3, column=0, pady=20)

        # Botón para obtener las opciones seleccionadas
        #boton_obtener_seleccionadas = tk.Button(self.frame, text="Continuar", command=obtener_seleccion)
        boton_obtener_seleccionadas = customtkinter.CTkButton(self.frame, text="Continuar", corner_radius=30, fg_color="#BE272F", command=obtener_seleccion)
        boton_obtener_seleccionadas.grid(row=3, column=2, pady=20)

        self.root.geometry("800x500")

    def seleccionar(self):
        seleccion = self.lista_disponibles.get(tk.ACTIVE)
        if seleccion not in self.opciones_seleccionadas:
            self.opciones_seleccionadas.append(seleccion)
            self.lista_seleccionadas.insert(tk.END, seleccion)
            self.lista_disponibles.delete(tk.ACTIVE)
            self.opciones_disponibles.remove(seleccion)

    def eliminar(self):
        seleccion = self.lista_seleccionadas.get(tk.ACTIVE)
        if seleccion in self.opciones_seleccionadas:
            self.opciones_disponibles.append(seleccion)
            self.lista_disponibles.insert(tk.END, seleccion)
            self.lista_seleccionadas.delete(tk.ACTIVE)
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
        self.root.geometry('1000x800')
    
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
        nombre_label = customtkinter.CTkLabel(self.frame, text=f"{self.seqrecord.id} | {self.seqrecord.description}", text_color="#BE272F")
        nombre_label.grid(row=0, column=0, columnspan=4)

        genes = diseno.get_all_genes(self.seqrecord)
        #descripcion_label = tk.Label(self.frame, text=f"Descripción: {self.seqrecord.description}")
        descripcion_label = customtkinter.CTkLabel(self.frame, text=f"Genes: {str(genes)[1:-1]}", text_color="#BE272F")
        descripcion_label.grid(row=1, column=0, columnspan=4)

        #largo_label = tk.Label(self.frame, text=f"Largo de la secuencia: {len(self.seqrecord)}")
        largo_label = customtkinter.CTkLabel(self.frame, text=f"Largo de la secuencia: {len(self.seqrecord)}", text_color="#BE272F")
        largo_label.grid(row=2, column=0, columnspan=4)

        locations = []
        for f in self.seqrecord.features:
            if f.type == 'CDS':
                loc = f.location
                if str(loc) not in locations:
                    locations.append(str(loc))

        #transcripciones_label = tk.Label(self.frame, text=f"Se considerará(n) {len(self.transcripciones)} de las {len(locations)} transcripciones.")
        transcripciones_label = customtkinter.CTkLabel(self.frame, text=f"Se considerará(n) {len(self.transcripciones)} de las {len(locations)} transcripciones. ({len(diseno.get_splicings(self.seqrecord, self.transcripciones))} puntos de empalme)", text_color="#BE272F")
        transcripciones_label.grid(row=3, column=0, columnspan=4)

        #fuente_negrita = font.Font(weight="bold")
        #parametros_label = tk.Label(self.frame, text="Ingrese los parámetros para el diseño de las sondas:", font=fuente_negrita)
        fuente_titulo = customtkinter.CTkFont(family='Helvetica', size=20, weight='bold')
        parametros_label = customtkinter.CTkLabel(self.frame, text="Ingrese los parámetros para el diseño de las sondas:", text_color="#BE272F", font=fuente_titulo)
        parametros_label.grid(row=4, column=0, columnspan=3, padx=20, pady=10)

        def mostrar_ayuda():
            ayuda_texto = "Usa el cursor para obtener más información de los parámetros"
            messagebox.showinfo("Ayuda", ayuda_texto)

        #boton_ayuda = tk.Button(self.frame, text="?", command=mostrar_ayuda)
        boton_ayuda = customtkinter.CTkButton(master=self.frame, text="?", fg_color="#E89434", corner_radius=20, width=10, height=10,command=mostrar_ayuda)
        boton_ayuda.grid(row=4, column=3)

        #Largo mínimo
        #minlen_label = tk.Label(self.frame, text="Largo mínimo (nt)")
        minlen_label = customtkinter.CTkLabel(self.frame, text="Largo mínimo (nt):", text_color="#BE272F")
        minlen_label.grid(row=5, column=0, padx=5, pady=10, sticky="e")
        
        #minlen_spinbox = tk.Spinbox(self.frame, from_=10, to=300, increment=1, width=5)
        minlen_spinbox = IntegerSelector(self.frame, default_value=60, min_value=10, max_value=500, increment=1)
        minlen_spinbox.grid(row=5, column=1, pady=10)
        #minlen_spinbox.delete(0, "end")
        #minlen_spinbox.insert(0, "60")


        #Largo máximo
        #maxlen_label = tk.Label(self.frame, text="Largo máximo (nt)")
        maxlen_label = customtkinter.CTkLabel(self.frame, text="Largo máximo (nt):", text_color="#BE272F")
        maxlen_label.grid(row=6, column=0, padx=5, sticky="e")

        #maxlen_spinbox = tk.Spinbox(self.frame, from_=10, to=300, increment=1, width=5)
        maxlen_spinbox = IntegerSelector(self.frame, default_value=120, min_value=10, max_value=500, increment=1)
        maxlen_spinbox.grid(row=6, column=1, padx= 5)
        #maxlen_spinbox.delete(0, "end")
        #maxlen_spinbox.insert(0, "120")

        texto_largo = "Largo de las sondas: Corresponde a la longitud en\nnucleótidos de la sonda a diseñar. Se prioriza\nel diseño de sondas de menor tamaño."
        tooltip_minlen = ToolTip(minlen_label, texto_largo, os.path.join("images","largo.png"))
        tooltip_maxlen = ToolTip(maxlen_label, texto_largo, os.path.join("images","largo.png"))

        #TM mínima
        #tmmin_label = tk.Label(self.frame, text="Temp melting mínima (°C)")
        tmmin_label = customtkinter.CTkLabel(self.frame, text="Temp melting mínima (°C):", text_color="#BE272F")
        tmmin_label.grid(row=5, column=2, padx=5, pady=10, sticky="e")
        
        #tmmin_spinbox = tk.Spinbox(self.frame, from_=10, to=300, increment=1, width=5)
        tmmin_spinbox = IntegerSelector(self.frame, default_value=65, min_value=10, max_value=300, increment=1)
        tmmin_spinbox.grid(row=5, column=3, padx= 5, pady=10)
        #tmmin_spinbox.delete(0, "end")
        #tmmin_spinbox.insert(0, "65")


        #TM máxima
        #tmmax_label = tk.Label(self.frame, text="Temp melting máxima (°C)")
        tmmax_label = customtkinter.CTkLabel(self.frame, text="Temp melting máxima (°C):", text_color="#BE272F")
        tmmax_label.grid(row=6, column=2, padx=5, sticky="e")
        
        #tmmax_spinbox = tk.Spinbox(self.frame, from_=10, to=300, increment=1, width=5)
        tmmax_spinbox = IntegerSelector(self.frame, default_value=80, min_value=10, max_value=300, increment=1)
        tmmax_spinbox.grid(row=6, column=3, padx= 5)
        #tmmax_spinbox.delete(0, "end")
        #tmmax_spinbox.insert(0, "80")

        texto_tm = "La temperatura melting (Tm) o temperatura de fusión de\nuna secuencia se refiere a la temperatura en la cual\nse desnaturalizaría la doble hebra. Más en concreto\nes la temperatura en la cual 50 % de las copias de esa secuencia\npresentes en una reacción se encuentra en forma monocatenaria\ny 50 % en forma bicatenaria, que interactúan con su secuencia\ncomplementaria."
        tooltip_tmmin = ToolTip(tmmin_label, texto_tm)
        tooltip_tmmax = ToolTip(tmmax_label, texto_tm)

        #%GC mínimo
        #gcmin_label = tk.Label(self.frame, text="%GC mínimo")
        gcmin_label = customtkinter.CTkLabel(self.frame, text="%GC mínimo:", text_color="#BE272F")
        gcmin_label.grid(row=8, column=0, padx=5, pady=10, sticky="e")
        
        #gcmin_spinbox = tk.Spinbox(self.frame, from_=0, to=100, increment=1, width=5)
        gcmin_spinbox = IntegerSelector(self.frame, default_value=30, min_value=0, max_value=100, increment=1)
        gcmin_spinbox.grid(row=8, column=1, padx= 5, pady=10)
        #gcmin_spinbox.delete(0, "end")
        #gcmin_spinbox.insert(0, "30")


        #%GC máximo
        #gcmax_label = tk.Label(self.frame, text="%GC máximo")
        gcmax_label = customtkinter.CTkLabel(self.frame, text="%GC máximo:", text_color="#BE272F")
        gcmax_label.grid(row=9, column=0, padx=5, sticky="e")

        #gcmax_spinbox = tk.Spinbox(self.frame, from_=0, to=100, increment=1, width=5)
        gcmax_spinbox = IntegerSelector(self.frame, default_value=70, min_value=0, max_value=100, increment=1)
        gcmax_spinbox.grid(row=9, column=1, padx= 5)
        #gcmax_spinbox.delete(0, "end")
        #gcmax_spinbox.insert(0, "70")

        texto_gc = "Porcentaje de GC: Se refiere a la fracción de bases\nnitrogenadas que son citosinas o guaninas dentro de\nla secuancia de nucleótidos."
        tooltip_gcmin = ToolTip(gcmin_label, texto_gc)
        tooltip_gcmax = ToolTip(gcmax_label, texto_gc)

        #Distancia mínima al borde del exón
        #mindist_label = tk.Label(self.frame, text="Distancia mínima al borde del exón (nt)")
        mindist_label = customtkinter.CTkLabel(self.frame, text="Distancia mínima al borde del exón (nt):", text_color="#BE272F")
        mindist_label.grid(row=8, column=2, padx=5, pady=10, sticky="e")
        
        #mindist_spinbox = tk.Spinbox(self.frame, from_=0, to=500, increment=1, width=5)
        mindist_spinbox = IntegerSelector(self.frame, default_value=0, min_value=0, max_value=500, increment=1)
        mindist_spinbox.grid(row=8, column=3, padx= 5, pady=10)
        #mindist_spinbox.delete(0, "end")
        #mindist_spinbox.insert(0, "0")


        #Distancia máxima al borde del exón
        #maxdist_label = tk.Label(self.frame, text="Distancia máxima al borde del exón (nt)")
        maxdist_label = customtkinter.CTkLabel(self.frame, text="Distancia máxima al borde del exón (nt):", text_color="#BE272F")
        maxdist_label.grid(row=9, column=2, padx=5, sticky="e")
        
        #maxdist_spinbox = tk.Spinbox(self.frame, from_=0, to=500, increment=1, width=5)
        maxdist_spinbox = IntegerSelector(self.frame, default_value=100, min_value=0, max_value=500, increment=1)
        maxdist_spinbox.grid(row=9, column=3, padx= 5)
        #maxdist_spinbox.delete(0, "end")
        #maxdist_spinbox.insert(0, "100")

        texto_dist = "Se refiere a la cantidad de nucleótidos que hay entre\nel punto de empalme y las sondas internas donor y acceptor.\n(No aplica para sondas centrales)"
        tooltip_mindist = ToolTip(mindist_label, texto_dist, os.path.join("images","distancia.png"))
        tooltip_maxdist = ToolTip(maxdist_label, texto_dist, os.path.join("images","distancia.png"))


        #Sobrelape mínimo
        #minoverlap_label = tk.Label(self.frame, text="Sobrelape mínimo (%)")
        minoverlap_label = customtkinter.CTkLabel(self.frame, text="Sobrelape mínimo (%):", text_color="#BE272F")
        minoverlap_label.grid(row=11, column=0, padx=5, pady=10, sticky="e")
        
        #minoverlap_spinbox = tk.Spinbox(self.frame, from_=0, to=100, increment=1, width=5)
        minoverlap_spinbox = IntegerSelector(self.frame, default_value=25, min_value=0, max_value=100, increment=1)
        minoverlap_spinbox.grid(row=11, column=1, padx= 5, pady=10)
        #minoverlap_spinbox.delete(0, "end")
        #minoverlap_spinbox.insert(0, "25")


        #Sobrelape máximo
        #maxoverlap_label = tk.Label(self.frame, text="Sobrelape máximo (%)")
        maxoverlap_label = customtkinter.CTkLabel(self.frame, text="Sobrelape máximo (%):", text_color="#BE272F")
        maxoverlap_label.grid(row=12, column=0, padx=5, sticky="e")

        maxoverlap_spinbox = tk.Spinbox(self.frame, from_=0, to=100, increment=1, width=5)
        maxoverlap_spinbox = IntegerSelector(self.frame, default_value=50, min_value=0, max_value=100, increment=1)
        maxoverlap_spinbox.grid(row=12, column=1, padx= 5)
        #maxoverlap_spinbox.delete(0, "end")
        #maxoverlap_spinbox.insert(0, "50")

        texto_overlap = "El sobrelape u overlap es el porcentaje de la\nsonda interior que es compartido con las sondas\nexteriores, solo aplica para donor y acceptor."
        tooltip_minoverlap = ToolTip(minoverlap_label, texto_overlap, os.path.join("images","sobrelape.png"))
        tooltip_maxoverlap = ToolTip(maxoverlap_label, texto_overlap, os.path.join("images","sobrelape.png"))

        #Delta G mínimo homodimerización
        #dgmin_homodim_label = tk.Label(self.frame, text="Delta G mínimo homodimerización")
        dgmin_homodim_label = customtkinter.CTkLabel(self.frame, text="Delta G mínimo homodimerización:", text_color="#BE272F")
        dgmin_homodim_label.grid(row=11, column=2, padx=5, pady=10, sticky="e")
        
        #dgmin_homodim_spinbox = tk.Spinbox(self.frame, from_=-50000, to=10000, increment=1000, width=5)
        dgmin_homodim_spinbox = IntegerSelector(self.frame, default_value=-10000, min_value=-999999, max_value=1000000, increment=1000)
        dgmin_homodim_spinbox.grid(row=11, column=3, padx= 5, pady=10)
        #dgmin_homodim_spinbox.delete(0, "end")
        #dgmin_homodim_spinbox.insert(0, "-10000")


        #Delta G mínimo hairpin
        #dgmin_hairpin_label = tk.Label(self.frame, text="Delta G mínimo hairpin")
        dgmin_hairpin_label = customtkinter.CTkLabel(self.frame, text="Delta G mínimo hairpin:", text_color="#BE272F")
        dgmin_hairpin_label.grid(row=12, column=2, padx=5, sticky="e")
        
        #dgmin_hairpin_spinbox = tk.Spinbox(self.frame, from_=-50000, to=10000, increment=1000, width=6)
        dgmin_hairpin_spinbox = IntegerSelector(self.frame, default_value=-10000, min_value=-999999, max_value=1000000, increment=1000)
        dgmin_hairpin_spinbox.grid(row=12, column=3, padx= 5)
        #dgmin_hairpin_spinbox.delete(0, "end")
        #dgmin_hairpin_spinbox.insert(0, "-10000")

        texto_dgmin_homodim = "Se refiere a la evaluación de la capacidad\nde la secuencia para formar estructuras de dímeros\nhomólogos. Los dímeros homólogos son la unión de\ndos secuencias de ADN idénticas o muy similares\nentre sí. Básicamente que la sonda no se complemente\ncon sí misma. Esto se verifica con la librería\nprimer3, que calcula la energía máxima entre las\ndos secuencias iguales."
        tooltip_dgmin_homodim = ToolTip(dgmin_homodim_label, texto_dgmin_homodim)
        texto_hairpin_label = "Capacidad de una secuencia para formar estructuras\nde horquilla (hairpin) en sí misma. Una estructura\nde hairpin es una conformación en la que una región\nde la secuencia se pliega hacia atrás y se empareja\ncon su complementaria, formando una estructura en forma de\nhorquilla."
        tooltip_hairpin_label = ToolTip(dgmin_hairpin_label, texto_hairpin_label)

        #Máximo homopolímeros simples
        #maxhomopol_simple_label = tk.Label(self.frame, text="Máximo homopolímeros simples")
        maxhomopol_simple_label = customtkinter.CTkLabel(self.frame, text="Máximo homopolímeros simples:", text_color="#BE272F")
        maxhomopol_simple_label.grid(row=14, column=0, padx=5, sticky="e")
        
        #maxhomopol_simple_spinbox = tk.Spinbox(self.frame, from_=2, to=100, increment=1, width=5)
        maxhomopol_simple_spinbox = IntegerSelector(self.frame, default_value=6, min_value=2, max_value=1000, increment=1)
        maxhomopol_simple_spinbox.grid(row=14, column=1, padx= 5)
        #maxhomopol_simple_spinbox.delete(0, "end")
        #maxhomopol_simple_spinbox.insert(0, "6")


        #Máximo homopolímeros dobles
        #maxhomopol_double_label = tk.Label(self.frame, text="Máximo homopolímeros dobles")
        maxhomopol_double_label = customtkinter.CTkLabel(self.frame, text="Máximo homopolímeros dobles:", text_color="#BE272F")
        maxhomopol_double_label.grid(row=15, column=0, padx=5, pady=10, sticky="e")

        #maxhomopol_double_spinbox = tk.Spinbox(self.frame, from_=2, to=100, increment=1, width=5)
        maxhomopol_double_spinbox = IntegerSelector(self.frame, default_value=5, min_value=2, max_value=1000, increment=1)
        maxhomopol_double_spinbox.grid(row=15, column=1, padx= 5, pady=10)
        #maxhomopol_double_spinbox.delete(0, "end")
        #maxhomopol_double_spinbox.insert(0, "5")

        #Máximo homopolímeros triples
        #maxhomopol_triple_label = tk.Label(self.frame, text="Máximo homopolímeros triples")
        maxhomopol_triple_label = customtkinter.CTkLabel(self.frame, text="Máximo homopolímeros triples:", text_color="#BE272F")
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


        multiplex_var = customtkinter.StringVar(value="off")
        checkbox_multiplex = customtkinter.CTkCheckBox(self.frame, text="Multiplexar sondas", command=param_multiplex, variable=multiplex_var, onvalue="on", offvalue="off", fg_color="#BE272F", text_color="#BE272F", hover_color="#E89434")
        checkbox_multiplex.deselect()
        checkbox_multiplex.grid(row=14, column=2, padx= 5,sticky='e')

        mindg_label = customtkinter.CTkLabel(self.frame, text="Mínimo delta G heterodimerización:", text_color="#BE272F")
        mindg_spinbox = IntegerSelector(self.frame, default_value=-20000, min_value=-999999, max_value=1000000, increment=1000)
    
        maxdt_label = customtkinter.CTkLabel(self.frame, text="Máxima diferencia de Tm (°C):", text_color="#BE272F")
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
        volver_button = customtkinter.CTkButton(master=self.frame, text="Volver", fg_color="#E89434", corner_radius=30, width=70 ,command=self.volver_a_transcripciones)
        continuar_button = customtkinter.CTkButton(master=self.frame, text="Continuar a la ejecución", fg_color="#BE272F", corner_radius=30,command=ejecutar)
        volver_button.grid(row=20, column=1, pady=20)
        continuar_button.grid(row=20, column=2, pady=20)

        self.root.geometry("800x1000")


    def volver_a_transcripciones(self):
        app.mostrar_pantalla("transcripciones", seqrecord=self.seqrecord, filepath=self.filepath)

class PantallaCarga(PantallaInicial):
    def __init__(self, root):
        self.root = root
        self.root.title("Cargando...")
        self.root.geometry("400x200")
        self.frame = ttk.Frame(self.root)
        #self.frame = customtkinter.CTkFrame(master=self.root)
        self.frame.pack(padx=20, pady=20)
        #self.label = ttk.Label(self.frame, text="Cargando, por favor espere...")
        self.label = customtkinter.CTkLabel(self.frame, text="Cargando, por favor espere...", text_color="#BE272F")
        self.label.pack()
        self.progreso = ttk.Progressbar(self.frame, mode="determinate", length=300)
        self.progreso.pack(pady=10)
        #self.actualizar_barra_de_progreso()

    def actualizar_progreso(self, progreso_queue):
        while True:
            try:
                progreso_actual = progreso_queue.get_nowait()
                self.progreso["value"] = progreso_actual
                root.update_idletasks()
            except queue.Empty:
                if tarea_thread.is_alive():
                    continue
                break

    def actualizar_barra_de_progreso(self):
        global progreso_compartido
        self.progreso["value"] = progreso_compartido
        self.frame.update()
        #print(progreso_compartido,'xd')
        if progreso_compartido < 100:            
            self.frame.after(10000, self.actualizar_barra_de_progreso())
    

class PantallaFinal(PantallaInicial):
    def __init__(self, root, seqrecord, transcripciones, df, dict_params, filepath):
        self.seqrecord = seqrecord
        self.transcripciones = transcripciones
        self.df = df
        self.dict_params = dict_params
        self.filepath = filepath
        super().__init__(root)

    def crear_interfaz(self):
        
        filename = diseno.generate_xlsx(df=self.df,
                                        name=self.seqrecord.id,
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
                                        mindg=self.dict_params['mindg'])
        
        imagen.plot_isoforms(self.seqrecord,
                             self.transcripciones,
                             imagen.get_splicings(self.seqrecord, self.transcripciones),
                             filename)
        
        folderpath = os.path.join(os.getcwd(),'sondas',filename)


        shutil.copy(self.filepath, folderpath)

        conn = sqlite3.connect(os.path.join(os.getcwd(), 'databases', 'probesdb.db'))
        cursor = conn.cursor()
        cursor.execute('''CREATE TABLE IF NOT EXISTS sondas (
                            id INTEGER PRIMARY KEY,
                            secuencia TEXT,
                            genes TEXT,
                            fecha_hora DATETIME,
                            carpeta TEXT)''')
        
        cursor.execute("INSERT INTO sondas (secuencia, genes, fecha_hora, carpeta) VALUES (?,?,?,?)", (str(self.seqrecord.id)+" | "+str(self.seqrecord.description), str(diseno.get_all_genes(self.seqrecord))[1:-1], datetime.datetime.now(), folderpath))
        conn.commit()
        conn.close()

        sondas = []
        for index, row in self.df.iterrows():
            if row['sonda'] != 'AAA' and row['sonda'] not in sondas:
                sondas.append(row['sonda'])

        #fuente_negrita = font.Font(weight="bold")
        #resultado_label = tk.Label(self.frame, text=f"El diseño se ejecutó con éxito. Número de sondas: {len(sondas)}", font=fuente_negrita)
        fuente_negrita = customtkinter.CTkFont(family='Helvetica', weight='bold')
        resultado_label = customtkinter.CTkLabel(self.frame, text=f"El diseño se ejecutó con éxito. Número de sondas: {len(sondas)}", text_color="#BE272F", font=fuente_negrita)
        resultado_label.pack(padx=20, pady=40)

        # img = Image.open(os.path.join(folderpath,self.seqrecord.id+".png"))
        # tkimg = ImageTk.PhotoImage(img)
        # label_imagen = tk.Label(self.frame, image=tkimg)
        # label_imagen.image = tkimg
        # label_imagen.pack(padx=10, pady=10)

        #carpeta_label = tk.Label(self.frame, text=f"El reporte y la imagen se almacenaron en\n" + folderpath)
        carpeta_label = customtkinter.CTkLabel(self.frame, text=f"El reporte y la imagen se almacenaron en\n" + folderpath, text_color="#BE272F")
        carpeta_label.pack(padx=20, pady=20)



        def abrir_carpeta():
            sistema_operativo = platform.system()
            if sistema_operativo == 'Windows':
                os.system(f'explorer "{folderpath}"')
            elif sistema_operativo == 'Darwin':  # macOS
                os.system(f'open "{folderpath}"')
            elif sistema_operativo == 'Linux':
                os.system(f'xdg-open "{folderpath}"')

        boton_imagen = customtkinter.CTkButton(self.frame, text="Abrir carpeta", corner_radius=30, fg_color="#BE272F", command=abrir_carpeta)
        boton_imagen.pack(pady=20)

        #boton_volver = tk.Button(self.frame, text="Volver al inicio", command=self.volver_a_inicio)
        boton_volver = customtkinter.CTkButton(self.frame, text="Volver al inicio", corner_radius=30, fg_color="#E89434", command=self.volver_a_inicio)
        boton_volver.pack(pady=20)

        self.root.geometry("1200x800")

    
    def volver_a_inicio(self):
        app.mostrar_pantalla('inicial')

class PantallaHistorial(PantallaInicial):
    def __init__(self, root):
        super().__init__(root)

    def crear_interfaz(self):
        conexion = sqlite3.connect(os.path.join(os.getcwd(),"databases", "probesdb.db"))
        cursor = conexion.cursor()
        cursor.execute('''CREATE TABLE IF NOT EXISTS sondas (
                            id INTEGER PRIMARY KEY,
                            secuencia TEXT,
                            genes TEXT,
                            fecha_hora DATETIME,
                            carpeta TEXT)''')
        cursor.execute("SELECT * FROM sondas ORDER BY id DESC")
        registros = cursor.fetchall()
        conexion.close()


        genes = []
        for registro in registros:
            if registro[2] not in genes:
                genes.append(registro[2])

        def maxlen(a):
            maximo = 0
            for i in a:
                length = len(i)
                if length > maximo:
                    maximo = length
            return maximo
        
        def filtrar_descripcion(value):
            conexion = sqlite3.connect(os.path.join(os.getcwd(),"databases", "probesdb.db"))
            cursor = conexion.cursor()
            cursor.execute("SELECT * FROM sondas WHERE genes = (?) ORDER BY id DESC",(value,))
            registros = cursor.fetchall()
            conexion.close()

            texto_registros.configure(state="normal")
            texto_registros.delete("0.0", "end")
            for registro in registros:
                texto_registros.insert("end", f"ID: {registro[0]}\n")
                texto_registros.insert("end", f"Secuencia: {registro[1]}\n")
                texto_registros.insert("end", f"Genes: {registro[2]}\n")
                texto_registros.insert("end", f"Fecha y Hora: {registro[3]}\n")
                texto_registros.insert("end", f"Carpeta: {registro[4]}\n\n")
            texto_registros.configure(state="disabled")

        
        #fuente_titulo = font.Font(weight="bold", size=16)
        #historial_label = tk.Label(self.frame, text=f"Historial de paneles de sondas", font=fuente_titulo)
        fuente_titulo = customtkinter.CTkFont(family='Helvetica', size=20, weight='bold')
        historial_label = customtkinter.CTkLabel(self.frame, text=f"Historial de paneles de sondas", text_color="#BE272F", font=fuente_titulo)
        historial_label.grid(row=0, column=0, padx=20, pady=20)

        #combobox_filtro = ttk.Combobox(self.frame, values=genes, state="readonly", width=maxlen(genes))
        combobox_filtro = customtkinter.CTkComboBox(self.frame, values=genes, width=100, command=filtrar_descripcion)
        combobox_filtro.grid(row=1, column=0, padx=20, pady=20)

        #boton_filtro = tk.Button(self.frame, text="Filtrar", command=filtrar_descripcion)
        #boton_filtro.grid(row=1, column=1, padx=20, pady=20)

        #texto_registros = tk.Text(self.frame, wrap=tk.WORD, width=100, height=20)
        texto_registros = customtkinter.CTkTextbox(self.frame, width=400, text_color="#BE272F")
        texto_registros.grid(row=2, column=0, padx=20, pady=20)

        for registro in registros:
            texto_registros.insert("end", f"ID: {registro[0]}\n")
            texto_registros.insert("end", f"Secuencia: {registro[1]}\n")
            texto_registros.insert("end", f"Genes: {registro[2]}\n")
            texto_registros.insert("end", f"Fecha y Hora: {registro[3]}\n")
            texto_registros.insert("end", f"Carpeta: {registro[4]}\n\n")
        
        texto_registros.configure(state="disabled")

        #boton_volver = tk.Button(self.frame, text="Volver", command=self.volver_a_seleccion)
        boton_volver = customtkinter.CTkButton(self.frame, text="Volver", corner_radius=30, fg_color="#BE272F", command=self.volver_a_seleccion)
        boton_volver.grid(row=3, column=0,columnspan=2, pady=20)

        self.root.geometry("800x600")

    
    def volver_a_seleccion(self):
        app.mostrar_pantalla("inicial")


class ControladorApp:
    def __init__(self, root):
        self.root = root
        self.pantallas = {
            "inicial": PantallaInicial(root),
            "transcripciones": None,
            "parametros": None,
            "carga": None,
            "historial": None,
            "final": None
        }
        self.mostrar_pantalla("inicial")

    def mostrar_pantalla(self, nombre, seqrecord=None, transcripciones=None, dict_params=None, df=None, filepath=None):
        # Ocultar pantalla actual
        if hasattr(self, "pantalla_actual"):
            self.pantalla_actual.frame.pack_forget()
        # Mostrar la pantalla solicitada
        if nombre == "transcripciones":
            self.pantalla_actual = PantallaTranscripciones(self.root, seqrecord, filepath)
            self.pantalla_actual.frame.pack(padx=50, pady=50)
        elif nombre == "parametros":
            self.pantalla_actual = PantallaParametros(self.root, seqrecord, transcripciones, filepath)
            self.pantalla_actual.frame.pack(padx=50, pady=50)
        elif nombre == "carga":
            self.carga = tk.Toplevel(self.root)
            self.pantalla_actual = PantallaCarga(self.carga)
            self.root.after(200, lambda: self.disenar_sondas(seqrecord, transcripciones, dict_params, filepath))
        elif nombre == "final":
            self.pantalla_actual = PantallaFinal(self.root, seqrecord, transcripciones, df, dict_params, filepath)
            self.pantalla_actual.frame.pack(padx=50, pady=50)
        elif nombre == "historial":
            self.pantalla_actual = PantallaHistorial(self.root)
            self.pantalla_actual.frame.pack(padx=50, pady=50)
        else:
            self.pantalla_actual = self.pantallas[nombre]
            self.pantalla_actual.frame.pack(padx=50, pady=50)
            self.root.geometry("800x450")

    def disenar_sondas(self, seqrecord, transcripciones, dict_params, filepath):
        global progreso_compartido
        progreso_compartido = 0
        progreso_queue = queue.Queue()
        df = diseno.probe_designer(seqrecord,
                                    transcripciones,
                                    progreso_queue,
                                    minlen=dict_params['minlen'],
                                    maxlen=dict_params['maxlen'],
                                    tmmin=dict_params['tmmin'],
                                    tmmax=dict_params['tmmax'],
                                    gcmin=dict_params['gcmin'],
                                    gcmax=dict_params['gcmax'],
                                    mindist=dict_params['mindist'],
                                    maxdist=dict_params['maxdist'],
                                    minoverlap=dict_params['minoverlap'],
                                    maxoverlap=dict_params['maxoverlap'],
                                    dgmin_homodim=dict_params['dgmin_homodim'],
                                    dgmin_hairpin=dict_params['dgmin_hairpin'],
                                    maxhomopol_simple=dict_params['maxhomopol_simple'],
                                    maxhomopol_double=dict_params['maxhomopol_double'],
                                    maxhomopol_triple=dict_params['maxhomopol_triple'],
                                    multiplex=dict_params['multiplex'],
                                    mindg=dict_params['mindg'],
                                    maxdt=dict_params['maxdt'])
        

        self.carga.destroy()  # Cierra la pantalla de carga
        self.mostrar_pantalla(nombre="final", seqrecord=seqrecord, transcripciones=transcripciones, df=df, dict_params=dict_params, filepath=filepath)


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

class IntegerSelector(customtkinter.CTkFrame):
    def __init__(self, master, default_value=0, min_value=None, max_value=None, increment=1, **kwargs):
        super().__init__(master, **kwargs)
        
        self.min_value = min_value
        self.max_value = max_value
        self.increment = increment
        
        self.valor = tk.IntVar()
        self.valor.set(default_value)
        
        #self.btn_disminuir = tk.Button(self, text="-", command=self.disminuir)
        self.btn_disminuir = customtkinter.CTkButton(master=self, text="-", fg_color="#BE272F", corner_radius=30, width=10, height=10,command=self.disminuir)
        self.btn_disminuir.pack(side="left")

        self.entry = tk.Entry(self, textvariable=self.valor, validate="key", validatecommand=(self.register(self.validar), "%P"), width=15, justify='center')
        #self.entry = customtkinter.CTkEntry(self, placeholder_text="CTkEntry")
        self.entry.pack(side="left", padx=5)
        
        #self.btn_aumentar = tk.Button(self, text="+", command=self.aumentar)
        self.btn_aumentar = customtkinter.CTkButton(master=self, text="+", fg_color="#BE272F", corner_radius=30, width=15, height=10,command=self.aumentar)
        self.btn_aumentar.pack(side="left")
    
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

if __name__ == "__main__":
    customtkinter.set_appearance_mode("light")
    #customtkinter.set_default_color_theme("dark-blue")
    #root = tk.Tk()
    root = customtkinter.CTk()
    app = ControladorApp(root)
    root.iconbitmap(os.path.join("images", "gemini2.ico"))
    #root.tk.call('wm', 'iconphoto', root._w, tk.PhotoImage(file=os.path.join('images','gemini2.png')))
    root.mainloop()
