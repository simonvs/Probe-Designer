import descarga
import diseno
import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import font
import time

# Función para simular una tarea que toma tiempo
def tarea_larga():
    time.sleep(3)
    return 0

# Función para mostrar la pantalla de carga
# def pantalla_carga():
#     carga_window = tk.Toplevel(root)
#     carga_window.title("Cargando...")
#     carga_label = tk.Label(carga_window, text='Ejecutando tarea')
#     carga_label.pack()
#     resultado = tarea_larga()
#     carga_window.destroy()
#     mostrar_resultado(resultado)

# Función para mostrar el resultado
# def mostrar_resultado(resultado):
#     resultado_window = tk.Toplevel(root)
#     resultado_window.title("Resultado")
#     resultado_label = tk.Label(resultado_window, text=f"La tarea se ejecutó con éxito. Resultado: {resultado}")
#     resultado_label.pack()
#     resultado_button = tk.Button(resultado_window, text="Cerrar", command=resultado_window.destroy)
#     resultado_button.pack()

# Función para seleccionar el archivo GenBank
def seleccionar_archivo():
    archivo = filedialog.askopenfilename(filetypes=[("Archivos GenBank", "*.gbk *.gb")])
    if archivo:
        seqrecord = descarga.parse_file_to_seqrecord(archivo)
        pantalla_transcripciones(seqrecord)


def pantalla_transcripciones(seqrecord):
    root.withdraw()
    transcripciones_window = tk.Toplevel(root)
    transcripciones_window.title("Diseñador de sondas de hibridación para ARN mensajero")

    transcripciones_label = tk.Label(transcripciones_window, text=f"Seleccione las transcripciones para el diseño de sondas", font=fuente_titulo)
    transcripciones_label.pack()

    #Obtener lista con transcripciones
    transcripts = []
    for f in seqrecord.features:
        if f.type == 'CDS':
            info = f.qualifiers.get('product')[0]
            if str(info) not in transcripts:
                transcripts.append(str(info))
    
    checkboxes = []
    variables = []

    for transcripcion in transcripts:
        var = tk.BooleanVar()
        checkbox = tk.Checkbutton(transcripciones_window, text=transcripcion, variable=var)
        checkbox.pack(pady=20)  # Alineación a la izquierda
        checkboxes.append(checkbox)
        variables.append(var)

    def obtener_seleccion():
        opciones_seleccionadas = []
        for i, var in enumerate(variables):
            if var.get():
                opciones_seleccionadas.append(transcripts[i])
        transcripciones_window.destroy()
        pantalla_parametros(seqrecord, opciones_seleccionadas)

    boton_parametros = tk.Button(transcripciones_window, text="Continuar", command=obtener_seleccion)
    boton_parametros.pack()


def pantalla_parametros(seqrecord, transcripts):
    # Cerrar la pantalla actual
    root.withdraw()

    parametros_window = tk.Toplevel(root)
    parametros_window.title("Diseñador de sondas de hibridación para ARN mensajero")
    #parametros_window.geometry("800x450")

    #Espaciadores
    parametros_window.grid_rowconfigure(7, minsize=20)
    parametros_window.grid_rowconfigure(10, minsize=20)
    parametros_window.grid_rowconfigure(13, minsize=30)
    parametros_window.grid_rowconfigure(17, minsize=20)
    parametros_window.grid_columnconfigure(4, minsize=20)

    nombre_label = tk.Label(parametros_window, text=f"Nombre de la secuencia: {seqrecord.id}")
    nombre_label.grid(row=0, column=0, columnspan=4)

    descripcion_label = tk.Label(parametros_window, text=f"Descripción: {seqrecord.description}")
    descripcion_label.grid(row=1, column=0, columnspan=4)

    largo_label = tk.Label(parametros_window, text=f"Largo de la secuencia: {len(seqrecord)}")
    largo_label.grid(row=2, column=0, columnspan=4)

    locations = []
    for f in seqrecord.features:
        if f.type == 'CDS':
            loc = f.location
            if str(loc) not in locations:
                locations.append(str(loc))

    transcripciones_label = tk.Label(parametros_window, text=f"Se considerará(n) {len(transcripts)} de las {len(locations)} transcripciones.")
    transcripciones_label.grid(row=3, column=0, columnspan=4)

    parametros_label = tk.Label(parametros_window, text="Ingrese los parámetros para el diseño de las sondas:", font=fuente_negrita)
    parametros_label.grid(row=4, column=0, columnspan=2, padx=10, pady=10)

    def mostrar_ayuda():
        ayuda_texto = "Este es un mensaje de ayuda.\nPuedes agregar información adicional aquí."
        messagebox.showinfo("Ayuda", ayuda_texto)

    boton_ayuda = tk.Button(parametros_window, text="?", command=mostrar_ayuda)
    boton_ayuda.grid(row=4, column=2)

    #Largo mínimo
    minlen_label = tk.Label(parametros_window, text="Largo mínimo (nt)")
    minlen_label.grid(row=5, column=0, padx=5, pady=10, sticky="e")
    
    minlen_spinbox = tk.Spinbox(parametros_window, from_=10, to=300, increment=1, width=5)
    minlen_spinbox.grid(row=5, column=1, pady=10)
    minlen_spinbox.delete(0, "end")
    minlen_spinbox.insert(0, "60")


    #Largo máximo
    maxlen_label = tk.Label(parametros_window, text="Largo máximo (nt)")
    maxlen_label.grid(row=6, column=0, padx=5, sticky="e")

    maxlen_spinbox = tk.Spinbox(parametros_window, from_=10, to=300, increment=1, width=5)
    maxlen_spinbox.grid(row=6, column=1, padx= 5)
    maxlen_spinbox.delete(0, "end")
    maxlen_spinbox.insert(0, "120")

    #TM mínima
    tmmin_label = tk.Label(parametros_window, text="Temp melting mínima (°C)")
    tmmin_label.grid(row=5, column=2, padx=5, pady=10, sticky="e")
    
    tmmin_spinbox = tk.Spinbox(parametros_window, from_=10, to=300, increment=1, width=5)
    tmmin_spinbox.grid(row=5, column=3, padx= 5, pady=10)
    tmmin_spinbox.delete(0, "end")
    tmmin_spinbox.insert(0, "65")


    #TM máxima
    tmmax_label = tk.Label(parametros_window, text="Temp melting máxima (°C)")
    tmmax_label.grid(row=6, column=2, padx=5, sticky="e")
    
    tmmax_spinbox = tk.Spinbox(parametros_window, from_=10, to=300, increment=1, width=5)
    tmmax_spinbox.grid(row=6, column=3, padx= 5)
    tmmax_spinbox.delete(0, "end")
    tmmax_spinbox.insert(0, "80")


    #%GC mínimo
    gcmin_label = tk.Label(parametros_window, text="%GC mínimo")
    gcmin_label.grid(row=8, column=0, padx=5, pady=10, sticky="e")
    
    gcmin_spinbox = tk.Spinbox(parametros_window, from_=0, to=100, increment=1, width=5)
    gcmin_spinbox.grid(row=8, column=1, padx= 5, pady=10)
    gcmin_spinbox.delete(0, "end")
    gcmin_spinbox.insert(0, "30")


    #%GC máximo
    gcmax_label = tk.Label(parametros_window, text="%GC máximo")
    gcmax_label.grid(row=9, column=0, padx=5, sticky="e")

    gcmax_spinbox = tk.Spinbox(parametros_window, from_=0, to=100, increment=1, width=5)
    gcmax_spinbox.grid(row=9, column=1, padx= 5)
    gcmax_spinbox.delete(0, "end")
    gcmax_spinbox.insert(0, "70")

    #Distancia mínima al borde del exón
    mindist_label = tk.Label(parametros_window, text="Distancia mínima al borde del exón (nt)")
    mindist_label.grid(row=8, column=2, padx=5, pady=10, sticky="e")
    
    mindist_spinbox = tk.Spinbox(parametros_window, from_=0, to=500, increment=1, width=5)
    mindist_spinbox.grid(row=8, column=3, padx= 5, pady=10)
    mindist_spinbox.delete(0, "end")
    mindist_spinbox.insert(0, "0")


    #Distancia máxima al borde del exón
    maxdist_label = tk.Label(parametros_window, text="Distancia máxima al borde del exón (nt)")
    maxdist_label.grid(row=9, column=2, padx=5, sticky="e")
    
    maxdist_spinbox = tk.Spinbox(parametros_window, from_=0, to=500, increment=1, width=5)
    maxdist_spinbox.grid(row=9, column=3, padx= 5)
    maxdist_spinbox.delete(0, "end")
    maxdist_spinbox.insert(0, "50")


    #Sobrelape mínimo
    minoverlap_label = tk.Label(parametros_window, text="Sobrelape mínimo (%)")
    minoverlap_label.grid(row=11, column=0, padx=5, pady=10, sticky="e")
    
    minoverlap_spinbox = tk.Spinbox(parametros_window, from_=0, to=100, increment=1, width=5)
    minoverlap_spinbox.grid(row=11, column=1, padx= 5, pady=10)
    minoverlap_spinbox.delete(0, "end")
    minoverlap_spinbox.insert(0, "25")


    #Sobrelape máximo
    maxoverlap_label = tk.Label(parametros_window, text="Sobrelape máximo (%)")
    maxoverlap_label.grid(row=12, column=0, padx=5, sticky="e")

    maxoverlap_spinbox = tk.Spinbox(parametros_window, from_=0, to=100, increment=1, width=5)
    maxoverlap_spinbox.grid(row=12, column=1, padx= 5)
    maxoverlap_spinbox.delete(0, "end")
    maxoverlap_spinbox.insert(0, "50")


    #Delta G mínimo homodimerización
    dgmin_homodim_label = tk.Label(parametros_window, text="Delta G mínimo homodimerización")
    dgmin_homodim_label.grid(row=11, column=2, padx=5, pady=10, sticky="e")
    
    dgmin_homodim_spinbox = tk.Spinbox(parametros_window, from_=-20000, to=1000, increment=1, width=5)
    dgmin_homodim_spinbox.grid(row=11, column=3, padx= 5, pady=10)
    dgmin_homodim_spinbox.delete(0, "end")
    dgmin_homodim_spinbox.insert(0, "-8000")


    #Delta G mínimo hairpin
    dgmin_hairpin_label = tk.Label(parametros_window, text="Delta G mínimo hairpin")
    dgmin_hairpin_label.grid(row=12, column=2, padx=5, sticky="e")
    
    dgmin_hairpin_spinbox = tk.Spinbox(parametros_window, from_=-20000, to=1000, increment=1, width=5)
    dgmin_hairpin_spinbox.grid(row=12, column=3, padx= 5)
    dgmin_hairpin_spinbox.delete(0, "end")
    dgmin_hairpin_spinbox.insert(0, "-8000")


    #Máximo homopolímeros simples
    maxhomopol_simple_label = tk.Label(parametros_window, text="Máximo homopolímeros simples")
    maxhomopol_simple_label.grid(row=14, column=0, padx=5, sticky="e")
    
    maxhomopol_simple_spinbox = tk.Spinbox(parametros_window, from_=2, to=100, increment=1, width=5)
    maxhomopol_simple_spinbox.grid(row=14, column=1, padx= 5)
    maxhomopol_simple_spinbox.delete(0, "end")
    maxhomopol_simple_spinbox.insert(0, "6")


    #Máximo homopolímeros dobles
    maxhomopol_double_label = tk.Label(parametros_window, text="Máximo homopolímeros dobles")
    maxhomopol_double_label.grid(row=15, column=0, padx=5, pady=10, sticky="e")

    maxhomopol_double_spinbox = tk.Spinbox(parametros_window, from_=2, to=100, increment=1, width=5)
    maxhomopol_double_spinbox.grid(row=15, column=1, padx= 5, pady=10)
    maxhomopol_double_spinbox.delete(0, "end")
    maxhomopol_double_spinbox.insert(0, "5")

    #Máximo homopolímeros triples
    maxhomopol_triple_label = tk.Label(parametros_window, text="Máximo homopolímeros triples")
    maxhomopol_triple_label.grid(row=16, column=0, padx=5, sticky="e")

    maxhomopol_triple_spinbox = tk.Spinbox(parametros_window, from_=2, to=100, increment=1, width=5)
    maxhomopol_triple_spinbox.grid(row=16, column=1, padx= 5)
    maxhomopol_triple_spinbox.delete(0, "end")
    maxhomopol_triple_spinbox.insert(0, "4")


    def volver_a_seleccion():
        parametros_window.destroy()
        root.deiconify()

    def pantalla_carga():
        carga_window = tk.Toplevel(root)
        carga_window.title("Cargando...")
        carga_label = tk.Label(carga_window, text='Ejecutando tarea, Por favor espere...', font=fuente_titulo)
        carga_label.pack()
        df = diseno.probe_designer(seqrecord,
                                   transcripts,
                                   minlen=int(minlen_spinbox.get()),
                                   maxlen=int(maxlen_spinbox.get()),
                                   tmmin=int(tmmin_spinbox.get()),
                                   tmmax=int(tmmax_spinbox.get()),
                                   gcmin=int(gcmin_spinbox.get()),
                                   gcmax=int(gcmax_spinbox.get()),
                                   mindist=int(mindist_spinbox.get()),
                                   maxdist=int(maxdist_spinbox.get()),
                                   minoverlap=int(minoverlap_spinbox.get()),
                                   maxoverlap=int(maxoverlap_spinbox.get()),
                                   dgmin_homodim=int(dgmin_homodim_spinbox.get()),
                                   dgmin_hairpin=int(dgmin_hairpin_spinbox.get()),
                                   maxhomopol_simple=int(maxhomopol_simple_spinbox.get()),
                                   maxhomopol_double=int(maxhomopol_double_spinbox.get()),
                                   maxhomopol_triple=int(maxhomopol_triple_spinbox.get()))
        #resultado = tarea_larga()
        carga_window.destroy()
        mostrar_resultado(df)

    def mostrar_resultado(df):
        resultado_window = tk.Toplevel(root)
        resultado_window.title("Resultado")
        resultado_label = tk.Label(resultado_window, text=f"La tarea se ejecutó con éxito. Nro sondas: {df.shape[0]}")
        resultado_label.pack()
        resultado_button = tk.Button(resultado_window, text="Cerrar", command=resultado_window.destroy)
        resultado_button.pack()

    # def continuar_a_ejecucion():
    #     # Pasar a la pantalla de carga
    #     pantalla_carga(seqrecord, transcripts)

    volver_button = tk.Button(parametros_window, text="Volver a la selección de secuencia", command=volver_a_seleccion)
    continuar_button = tk.Button(parametros_window, text="Continuar a la ejecución", command=pantalla_carga)
    volver_button.grid(row=20, column=1, pady=10)
    continuar_button.grid(row=20, column=2, pady=10)


# Crear la ventana principal
root = tk.Tk()
root.title("Diseñador de sondas de hibridación para ARN mensajero")
root.geometry("800x450")  # Tamaño fijo 16:9

fuente_titulo = font.Font(weight="bold", size=16)
fuente_negrita = font.Font(weight="bold")

titulo_archivo = tk.Label(root, text="Seleccione un archivo GenBank", font=fuente_titulo)
titulo_archivo.pack(pady=20)

# Pantalla de selección de archivo
seleccion_button = tk.Button(root, text="Seleccionar Archivo", command=seleccionar_archivo)
seleccion_button.pack(pady=20)

historial_button = tk.Button(root, text="Ver historial de sondas")
historial_button.pack(pady=20)

root.mainloop()
