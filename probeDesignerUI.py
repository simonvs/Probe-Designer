import descarga
import diseno
import imagen
import tkinter as tk
import os
import sqlite3
import datetime
import platform
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import font
from PIL import Image, ImageTk

class PantallaInicial:
    def __init__(self, root):
        self.root = root
        self.root.title("Diseñador de sondas de hibridación para ARN mensajero")
        self.frame = tk.Frame(root)
        self.frame.pack()
        self.crear_interfaz()

    def crear_interfaz(self):
        # Crear widgets de la pantalla de inicio
        fuente_titulo = font.Font(weight="bold", size=16)
        titulo_archivo = tk.Label(self.frame, text="Seleccione un archivo GenBank", font=fuente_titulo)
        titulo_archivo.pack(padx=100, pady=30)

        # Pantalla de selección de archivo
        seleccion_button = tk.Button(self.frame, text="Seleccionar Archivo", command=self.seleccionar_archivo)
        seleccion_button.pack(pady=30)

        historial_button = tk.Button(self.frame, text="Ver historial de sondas", command=self.ver_historial)
        historial_button.pack(pady=60)
        
    def seleccionar_archivo(self):
        archivo = filedialog.askopenfilename(filetypes=[("Archivos GenBank", "*.gbk *.gb")])
        if archivo:
            seqrecord = descarga.parse_file_to_seqrecord(archivo)
            app.mostrar_pantalla("transcripciones", seqrecord)
    
    def ver_historial(self):
        app.mostrar_pantalla('historial')


class PantallaTranscripciones(PantallaInicial):
    def __init__(self, root, seqrecord):
        self.seqrecord = seqrecord
        super().__init__(root)

    def crear_interfaz(self):
        # Crear widgets de la pantalla de inicio
        fuente_titulo = font.Font(weight="bold", size=16)
        transcripciones_label = tk.Label(self.frame, text=f"Seleccione las transcripciones para el diseño de sondas", font=fuente_titulo)
        transcripciones_label.grid(row=0, column=0, columnspan=2, padx=10, pady=20)

        frame_combobox = tk.Frame(self.frame)
        frame_combobox.grid(row=1, column=0)

        frame_listbox = tk.Frame(self.frame)
        frame_listbox.grid(row=1, column=1)

        #Obtener lista con transcripciones
        transcripts = []
        for f in self.seqrecord.features:
            if f.type == 'CDS':
                info = f.qualifiers.get('product')[0]
                if str(info) not in transcripts:
                    transcripts.append(str(info))
        
        def maxlen(a):
            maximo = 0
            for i in a:
                length = len(i)
                if length > maximo:
                    maximo = length
            return maximo

        # Crear Combobox para seleccionar opciones
        combobox = ttk.Combobox(frame_combobox, values=transcripts, state="readonly", width=maxlen(transcripts))
        combobox.pack(pady=20)

        # Función para agregar una opción al Listbox si no se ha seleccionado previamente
        def agregar_opcion():
            opcion = combobox.get()
            if opcion and opcion not in lista_seleccion.get(0, tk.END):
                lista_seleccion.insert(tk.END, opcion)

        # Función para seleccionar todas las opciones
        def seleccionar_todas():
            for opcion in transcripts:
                if opcion not in lista_seleccion.get(0, tk.END):
                    lista_seleccion.insert(tk.END, opcion)

        # Función para eliminar las opciones seleccionadas
        def eliminar_seleccionadas():
            seleccionados = lista_seleccion.curselection()
            for seleccionado in seleccionados:
                lista_seleccion.delete(seleccionado)

        
        def obtener_seleccion():
            seleccionados = lista_seleccion.get(0, tk.END)
            #pantalla parametreos
            app.mostrar_pantalla("parametros", self.seqrecord, [seleccionado for seleccionado in seleccionados])

        # Botón para agregar la opción seleccionada al Listbox
        boton_agregar = tk.Button(frame_combobox, text="Agregar Opción", command=agregar_opcion)
        boton_agregar.pack()

        # Listbox para mostrar las opciones seleccionadas
        lista_seleccion = tk.Listbox(frame_listbox, selectmode=tk.MULTIPLE, width=maxlen(transcripts))
        lista_seleccion.pack(pady=10)

        # Botón para seleccionar todas las opciones
        boton_seleccionar_todas = tk.Button(frame_listbox, text="Seleccionar Todas", command=seleccionar_todas)
        boton_seleccionar_todas.pack()

        # Botón para eliminar las opciones seleccionadas
        boton_eliminar_seleccionadas = tk.Button(frame_listbox, text="Eliminar Seleccionadas", command=eliminar_seleccionadas)
        boton_eliminar_seleccionadas.pack()

        # Botón para obtener las opciones seleccionadas
        boton_volver = tk.Button(self.frame, text="Volver", command=self.volver_a_seleccion)
        boton_volver.grid(row=2, column=0, pady=20)

        # Botón para obtener las opciones seleccionadas
        boton_obtener_seleccionadas = tk.Button(self.frame, text="Continuar", command=obtener_seleccion)
        boton_obtener_seleccionadas.grid(row=2, column=1, pady=20)


    def volver_a_seleccion(self):
        app.mostrar_pantalla("inicial")


class PantallaParametros(PantallaInicial):
    def __init__(self, root, seqrecord, transcripciones):
        self.seqrecord = seqrecord
        self.transcripciones = transcripciones
        super().__init__(root)
    
    def crear_interfaz(self):
        #Espaciadores
        self.frame.grid_rowconfigure(7, minsize=20)
        self.frame.grid_rowconfigure(10, minsize=20)
        self.frame.grid_rowconfigure(13, minsize=30)
        self.frame.grid_rowconfigure(17, minsize=20)
        self.frame.grid_columnconfigure(4, minsize=20)

        nombre_label = tk.Label(self.frame, text=f"Nombre de la secuencia: {self.seqrecord.id}")
        nombre_label.grid(row=0, column=0, columnspan=4)

        descripcion_label = tk.Label(self.frame, text=f"Descripción: {self.seqrecord.description}")
        descripcion_label.grid(row=1, column=0, columnspan=4)

        largo_label = tk.Label(self.frame, text=f"Largo de la secuencia: {len(self.seqrecord)}")
        largo_label.grid(row=2, column=0, columnspan=4)

        locations = []
        for f in self.seqrecord.features:
            if f.type == 'CDS':
                loc = f.location
                if str(loc) not in locations:
                    locations.append(str(loc))

        transcripciones_label = tk.Label(self.frame, text=f"Se considerará(n) {len(self.transcripciones)} de las {len(locations)} transcripciones.")
        transcripciones_label.grid(row=3, column=0, columnspan=4)

        fuente_negrita = font.Font(weight="bold")
        parametros_label = tk.Label(self.frame, text="Ingrese los parámetros para el diseño de las sondas:", font=fuente_negrita)
        parametros_label.grid(row=4, column=0, columnspan=2, padx=10, pady=10)

        def mostrar_ayuda():
            ayuda_texto = "Este es un mensaje de ayuda.\nPuedes agregar información adicional aquí."
            messagebox.showinfo("Ayuda", ayuda_texto)

        boton_ayuda = tk.Button(self.frame, text="?", command=mostrar_ayuda)
        boton_ayuda.grid(row=4, column=2)

        #Largo mínimo
        minlen_label = tk.Label(self.frame, text="Largo mínimo (nt)")
        minlen_label.grid(row=5, column=0, padx=5, pady=10, sticky="e")
        
        minlen_spinbox = tk.Spinbox(self.frame, from_=10, to=300, increment=1, width=5)
        minlen_spinbox.grid(row=5, column=1, pady=10)
        minlen_spinbox.delete(0, "end")
        minlen_spinbox.insert(0, "60")


        #Largo máximo
        maxlen_label = tk.Label(self.frame, text="Largo máximo (nt)")
        maxlen_label.grid(row=6, column=0, padx=5, sticky="e")

        maxlen_spinbox = tk.Spinbox(self.frame, from_=10, to=300, increment=1, width=5)
        maxlen_spinbox.grid(row=6, column=1, padx= 5)
        maxlen_spinbox.delete(0, "end")
        maxlen_spinbox.insert(0, "120")

        texto_largo = "Largo de las sondas: Corresponde a la longitud en\nnucleótidos de la sonda a diseñar. Se prioriza\nel diseño de sondas de menor tamaño."
        tooltip_minlen = ToolTip(minlen_label, texto_largo, os.path.join("images","largo.png"))
        tooltip_maxlen = ToolTip(maxlen_label, texto_largo, os.path.join("images","largo.png"))

        #TM mínima
        tmmin_label = tk.Label(self.frame, text="Temp melting mínima (°C)")
        tmmin_label.grid(row=5, column=2, padx=5, pady=10, sticky="e")
        
        tmmin_spinbox = tk.Spinbox(self.frame, from_=10, to=300, increment=1, width=5)
        tmmin_spinbox.grid(row=5, column=3, padx= 5, pady=10)
        tmmin_spinbox.delete(0, "end")
        tmmin_spinbox.insert(0, "65")


        #TM máxima
        tmmax_label = tk.Label(self.frame, text="Temp melting máxima (°C)")
        tmmax_label.grid(row=6, column=2, padx=5, sticky="e")
        
        tmmax_spinbox = tk.Spinbox(self.frame, from_=10, to=300, increment=1, width=5)
        tmmax_spinbox.grid(row=6, column=3, padx= 5)
        tmmax_spinbox.delete(0, "end")
        tmmax_spinbox.insert(0, "80")

        texto_tm = "La temperatura melting (Tm) o temperatura de fusión de\nuna secuencia se refiere a la temperatura en la cual\nse desnaturalizaría la doble hebra. Más en concreto\nes la temperatura en la cual 50 % de las copias de esa secuencia\npresentes en una reacción se encuentra en forma monocatenaria y\n50 % en forma bicatenaria, que interactúan con su secuencia complementaria."
        tooltip_tmmin = ToolTip(tmmin_label, texto_tm)
        tooltip_tmmax = ToolTip(tmmax_label, texto_tm)

        #%GC mínimo
        gcmin_label = tk.Label(self.frame, text="%GC mínimo")
        gcmin_label.grid(row=8, column=0, padx=5, pady=10, sticky="e")
        
        gcmin_spinbox = tk.Spinbox(self.frame, from_=0, to=100, increment=1, width=5)
        gcmin_spinbox.grid(row=8, column=1, padx= 5, pady=10)
        gcmin_spinbox.delete(0, "end")
        gcmin_spinbox.insert(0, "30")


        #%GC máximo
        gcmax_label = tk.Label(self.frame, text="%GC máximo")
        gcmax_label.grid(row=9, column=0, padx=5, sticky="e")

        gcmax_spinbox = tk.Spinbox(self.frame, from_=0, to=100, increment=1, width=5)
        gcmax_spinbox.grid(row=9, column=1, padx= 5)
        gcmax_spinbox.delete(0, "end")
        gcmax_spinbox.insert(0, "70")

        texto_gc = "Porcentaje de GC: Se refiere a la fracción de bases nitrogenadas que son citosinas o guaninas dentro de la secuancia de nucleótidos."
        tooltip_gcmin = ToolTip(gcmin_label, texto_largo, os.path.join("images","distancia.png"))
        tooltip_gcmax = ToolTip(gcmax_label, texto_largo, os.path.join("images","distancia.png"))

        #Distancia mínima al borde del exón
        mindist_label = tk.Label(self.frame, text="Distancia mínima al borde del exón (nt)")
        mindist_label.grid(row=8, column=2, padx=5, pady=10, sticky="e")
        
        mindist_spinbox = tk.Spinbox(self.frame, from_=0, to=500, increment=1, width=5)
        mindist_spinbox.grid(row=8, column=3, padx= 5, pady=10)
        mindist_spinbox.delete(0, "end")
        mindist_spinbox.insert(0, "0")


        #Distancia máxima al borde del exón
        maxdist_label = tk.Label(self.frame, text="Distancia máxima al borde del exón (nt)")
        maxdist_label.grid(row=9, column=2, padx=5, sticky="e")
        
        maxdist_spinbox = tk.Spinbox(self.frame, from_=0, to=500, increment=1, width=5)
        maxdist_spinbox.grid(row=9, column=3, padx= 5)
        maxdist_spinbox.delete(0, "end")
        maxdist_spinbox.insert(0, "100")


        #Sobrelape mínimo
        minoverlap_label = tk.Label(self.frame, text="Sobrelape mínimo (%)")
        minoverlap_label.grid(row=11, column=0, padx=5, pady=10, sticky="e")
        
        minoverlap_spinbox = tk.Spinbox(self.frame, from_=0, to=100, increment=1, width=5)
        minoverlap_spinbox.grid(row=11, column=1, padx= 5, pady=10)
        minoverlap_spinbox.delete(0, "end")
        minoverlap_spinbox.insert(0, "25")


        #Sobrelape máximo
        maxoverlap_label = tk.Label(self.frame, text="Sobrelape máximo (%)")
        maxoverlap_label.grid(row=12, column=0, padx=5, sticky="e")

        maxoverlap_spinbox = tk.Spinbox(self.frame, from_=0, to=100, increment=1, width=5)
        maxoverlap_spinbox.grid(row=12, column=1, padx= 5)
        maxoverlap_spinbox.delete(0, "end")
        maxoverlap_spinbox.insert(0, "50")


        #Delta G mínimo homodimerización
        dgmin_homodim_label = tk.Label(self.frame, text="Delta G mínimo homodimerización")
        dgmin_homodim_label.grid(row=11, column=2, padx=5, pady=10, sticky="e")
        
        dgmin_homodim_spinbox = tk.Spinbox(self.frame, from_=-50000, to=10000, increment=1000, width=5)
        dgmin_homodim_spinbox.grid(row=11, column=3, padx= 5, pady=10)
        dgmin_homodim_spinbox.delete(0, "end")
        dgmin_homodim_spinbox.insert(0, "-10000")


        #Delta G mínimo hairpin
        dgmin_hairpin_label = tk.Label(self.frame, text="Delta G mínimo hairpin")
        dgmin_hairpin_label.grid(row=12, column=2, padx=5, sticky="e")
        
        dgmin_hairpin_spinbox = tk.Spinbox(self.frame, from_=-50000, to=10000, increment=1000, width=6)
        dgmin_hairpin_spinbox.grid(row=12, column=3, padx= 5)
        dgmin_hairpin_spinbox.delete(0, "end")
        dgmin_hairpin_spinbox.insert(0, "-10000")


        #Máximo homopolímeros simples
        maxhomopol_simple_label = tk.Label(self.frame, text="Máximo homopolímeros simples")
        maxhomopol_simple_label.grid(row=14, column=0, padx=5, sticky="e")
        
        maxhomopol_simple_spinbox = tk.Spinbox(self.frame, from_=2, to=100, increment=1, width=5)
        maxhomopol_simple_spinbox.grid(row=14, column=1, padx= 5)
        maxhomopol_simple_spinbox.delete(0, "end")
        maxhomopol_simple_spinbox.insert(0, "6")


        #Máximo homopolímeros dobles
        maxhomopol_double_label = tk.Label(self.frame, text="Máximo homopolímeros dobles")
        maxhomopol_double_label.grid(row=15, column=0, padx=5, pady=10, sticky="e")

        maxhomopol_double_spinbox = tk.Spinbox(self.frame, from_=2, to=100, increment=1, width=5)
        maxhomopol_double_spinbox.grid(row=15, column=1, padx= 5, pady=10)
        maxhomopol_double_spinbox.delete(0, "end")
        maxhomopol_double_spinbox.insert(0, "5")

        #Máximo homopolímeros triples
        maxhomopol_triple_label = tk.Label(self.frame, text="Máximo homopolímeros triples")
        maxhomopol_triple_label.grid(row=16, column=0, padx=5, sticky="e")

        maxhomopol_triple_spinbox = tk.Spinbox(self.frame, from_=2, to=100, increment=1, width=5)
        maxhomopol_triple_spinbox.grid(row=16, column=1, padx= 5)
        maxhomopol_triple_spinbox.delete(0, "end")
        maxhomopol_triple_spinbox.insert(0, "4")

        def ejecutar():
            dict_params = {}
            dict_params['minlen'] = int(minlen_spinbox.get())
            dict_params['maxlen'] = int(maxlen_spinbox.get())
            dict_params['tmmin'] = int(tmmin_spinbox.get())
            dict_params['tmmax'] = int(tmmax_spinbox.get())
            dict_params['gcmin'] = int(gcmin_spinbox.get())
            dict_params['gcmax'] = int(gcmax_spinbox.get())
            dict_params['mindist'] = int(mindist_spinbox.get())
            dict_params['maxdist'] = int(maxdist_spinbox.get())
            dict_params['minoverlap'] = int(minoverlap_spinbox.get())
            dict_params['maxoverlap'] = int(maxoverlap_spinbox.get())
            dict_params['dgmin_homodim'] = int(dgmin_homodim_spinbox.get())
            dict_params['dgmin_hairpin'] = int(dgmin_hairpin_spinbox.get())
            dict_params['maxhomopol_simple'] = int(maxhomopol_simple_spinbox.get())
            dict_params['maxhomopol_double'] = int(maxhomopol_double_spinbox.get())
            dict_params['maxhomopol_triple'] = int(maxhomopol_triple_spinbox.get())
            app.mostrar_pantalla("carga", self.seqrecord, self.transcripciones, dict_params)

        volver_button = tk.Button(self.frame, text="Volver a la selección de transcripciones", command=self.volver_a_transcripciones)
        continuar_button = tk.Button(self.frame, text="Continuar a la ejecución", command=ejecutar)
        volver_button.grid(row=20, column=1, pady=10)
        continuar_button.grid(row=20, column=2, pady=10)


    def volver_a_transcripciones(self):
        app.mostrar_pantalla("transcripciones", self.seqrecord)


class PantallaCarga(PantallaInicial):
    def __init__(self, root):
        self.root = root
        self.root.title("Cargando...")
        self.root.geometry("300x100")
        self.frame = ttk.Frame(self.root)
        self.frame.pack(padx=20, pady=20)
        self.label = ttk.Label(self.frame, text="Cargando, por favor espere...")
        self.label.pack()
        self.progreso = ttk.Progressbar(self.frame, mode="indeterminate", length=200)
        self.progreso.pack()
        self.progreso.start()


class PantallaFinal(PantallaInicial):
    def __init__(self, root, seqrecord, transcripciones, df, dict_params):
        self.seqrecord = seqrecord
        self.transcripciones = transcripciones
        self.df = df
        self.dict_params = dict_params
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
                                        maxhomopol_triple=self.dict_params['maxhomopol_triple'])
        imagen.plot_isoforms(self.seqrecord,
                             self.transcripciones,
                             imagen.get_splicings(self.seqrecord, self.transcripciones),
                             filename)
        
        folderpath = os.path.join(os.getcwd(),'sondas',filename)

        conn = sqlite3.connect(os.path.join(os.getcwd(), 'databases', 'probesdb.db'))
        cursor = conn.cursor()
        cursor.execute('''CREATE TABLE IF NOT EXISTS sondas (
                            id INTEGER PRIMARY KEY,
                            accession TEXT,
                            descripcion TEXT,
                            fecha_hora DATETIME,
                            carpeta TEXT)''')
        
        cursor.execute("INSERT INTO sondas (accession, descripcion, fecha_hora, carpeta) VALUES (?,?,?,?)", (self.seqrecord.id, self.seqrecord.description, datetime.datetime.now(), folderpath))
        conn.commit()
        conn.close()

        sondas = []
        for index, row in self.df.iterrows():
            if row['sonda'] != 'AAA' and row['sonda'] not in sondas:
                sondas.append(row['sonda'])

        fuente_negrita = font.Font(weight="bold")
        resultado_label = tk.Label(self.frame, text=f"El diseño se ejecutó con éxito. Número de sondas: {len(sondas)}", font=fuente_negrita)
        resultado_label.pack(pady=40)

        carpeta_label = tk.Label(self.frame, text=f"El reporte y la imagen se almacenaron en " + folderpath)
        carpeta_label.pack(pady=20)

        boton_volver = tk.Button(self.frame, text="Volver al inicio", command=self.volver_a_inicio)
        boton_volver.pack(pady=20)


        sistema_operativo = platform.system()

        if sistema_operativo == 'Windows':
            os.system(f'explorer "{folderpath}"')
        elif sistema_operativo == 'Darwin':  # macOS
            os.system(f'open "{folderpath}"')
        elif sistema_operativo == 'Linux':
            os.system(f'xdg-open "{folderpath}"')
    
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
                            accession TEXT,
                            descripcion TEXT,
                            fecha_hora DATETIME,
                            carpeta TEXT)''')
        cursor.execute("SELECT * FROM sondas ORDER BY id DESC")
        registros = cursor.fetchall()
        conexion.close()


        descripciones = []
        for registro in registros:
            if registro[2] not in descripciones:
                descripciones.append(registro[2])

        def maxlen(a):
            maximo = 0
            for i in a:
                length = len(i)
                if length > maximo:
                    maximo = length
            return maximo
        
        def filtrar_descripcion():
            conexion = sqlite3.connect(os.path.join(os.getcwd(),"databases", "probesdb.db"))
            cursor = conexion.cursor()
            print(combobox_filtro.get())
            cursor.execute("SELECT * FROM sondas WHERE descripcion = (?) ORDER BY id DESC",(combobox_filtro.get(),))
            registros = cursor.fetchall()
            conexion.close()

            texto_registros.config(state="normal")
            texto_registros.delete("1.0", tk.END)
            for registro in registros:
                texto_registros.insert(tk.END, f"ID: {registro[0]}\n")
                texto_registros.insert(tk.END, f"Accesion Numer: {registro[1]}\n")
                texto_registros.insert(tk.END, f"Descripción: {registro[2]}\n")
                texto_registros.insert(tk.END, f"Fecha y Hora: {registro[3]}\n")
                texto_registros.insert(tk.END, f"Carpeta: {registro[4]}\n\n")
            texto_registros.config(state="disabled")

        
        fuente_titulo = font.Font(weight="bold", size=16)
        historial_label = tk.Label(self.frame, text=f"Historial de paneles de sondas", font=fuente_titulo)
        historial_label.grid(row=0, column=0, columnspan=2, padx=20, pady=20)

        combobox_filtro = ttk.Combobox(self.frame, values=descripciones, state="readonly", width=maxlen(descripciones))
        combobox_filtro.grid(row=1, column=0, padx=20, pady=20)

        boton_filtro = tk.Button(self.frame, text="Filtrar", command=filtrar_descripcion)
        boton_filtro.grid(row=1, column=1, padx=20, pady=20)

        texto_registros = tk.Text(self.frame, wrap=tk.WORD, width=100, height=20)
        texto_registros.grid(row=2, column=0, columnspan=2, padx=20, pady=20)

        for registro in registros:
            texto_registros.insert(tk.END, f"ID: {registro[0]}\n")
            texto_registros.insert(tk.END, f"Accesion Numer: {registro[1]}\n")
            texto_registros.insert(tk.END, f"Descripción: {registro[2]}\n")
            texto_registros.insert(tk.END, f"Fecha y Hora: {registro[3]}\n")
            texto_registros.insert(tk.END, f"Carpeta: {registro[4]}\n\n")
        
        texto_registros.config(state="disabled")

        boton_volver = tk.Button(self.frame, text="Volver", command=self.volver_a_seleccion)
        boton_volver.grid(row=3, column=0,columnspan=2, pady=20)

    
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

    def mostrar_pantalla(self, nombre, seqrecord=None, transcripciones=None, dict_params=None, df=None):
        # Ocultar pantalla actual
        if hasattr(self, "pantalla_actual"):
            self.pantalla_actual.frame.pack_forget()
        # Mostrar la pantalla solicitada
        if nombre == "transcripciones":
            self.pantalla_actual = PantallaTranscripciones(self.root, seqrecord)
            self.pantalla_actual.frame.pack()
        elif nombre == "parametros":
            self.pantalla_actual = PantallaParametros(self.root, seqrecord, transcripciones)
            self.pantalla_actual.frame.pack()
        elif nombre == "carga":
            self.carga = tk.Toplevel(self.root)
            self.pantalla_actual = PantallaCarga(self.carga)
            self.root.after(200, lambda: self.disenar_sondas(seqrecord, transcripciones, dict_params))
        elif nombre == "final":
            self.pantalla_actual = PantallaFinal(self.root, seqrecord, transcripciones, df, dict_params)
            self.pantalla_actual.frame.pack()
        elif nombre == "historial":
            self.pantalla_actual = PantallaHistorial(self.root)
            self.pantalla_actual.frame.pack()
        else:
            self.pantalla_actual = self.pantallas[nombre]
            self.pantalla_actual.frame.pack()

    def disenar_sondas(self, seqrecord, transcripciones, dict_params):
        df = diseno.probe_designer(seqrecord,
                                    transcripciones,
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
                                    maxhomopol_triple=dict_params['maxhomopol_triple'])
        print(len(df))
        self.carga.destroy()  # Cierra la pantalla de carga
        self.mostrar_pantalla(nombre="final", seqrecord=seqrecord, transcripciones=transcripciones, df=df, dict_params=dict_params)  # Muestra la pantalla de configuración

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


if __name__ == "__main__":
    root = tk.Tk()
    app = ControladorApp(root)
    root.mainloop()
