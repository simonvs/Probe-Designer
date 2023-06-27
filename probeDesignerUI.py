import tkinter as tk
from tkinter import filedialog

class Interfaz:
    def __init__(self):
        self.ventana = tk.Tk()
        self.ventana.title("Diseñador de Sondas de Hibridación para ARNm")
        
        self.pantalla_actual = None
        self.secuencias = []
        self.seleccionar_secuencias()
    
    def seleccionar_secuencias(self):
        if self.pantalla_actual is not None:
            self.pantalla_actual.destroy()
        
        self.pantalla_actual = tk.Frame(self.ventana)
        self.pantalla_actual.pack()
        
        etiqueta = tk.Label(self.pantalla_actual, text="Selecciona las Secuencias")
        etiqueta.pack()

        boton_seleccionar = tk.Button(self.pantalla_actual, text="Seleccionar Archivo", command=self.seleccionar_archivo)
        boton_seleccionar.pack()
        
        boton = tk.Button(self.pantalla_actual, text="Cambiar a Pantalla 2", command=self.seleccionar_parametros)
        boton.pack()
    
    def seleccionar_archivo(self):
        tipos_archivo = (("Archivos de texto", "*.txt"), ("Todos los archivos", "*.*"))
        archivo = filedialog.askopenfilename(filetypes=tipos_archivo)
        
        if archivo:
            print("Archivo seleccionado:", archivo)
        else:
            print("No se seleccionó ningún archivo.")

    def seleccionar_parametros(self):
        if self.pantalla_actual is not None:
            self.pantalla_actual.destroy()
        
        self.pantalla_actual = tk.Frame(self.ventana)
        self.pantalla_actual.pack()
        
        etiqueta = tk.Label(self.pantalla_actual, text="Pantalla 2")
        etiqueta.pack()
        
        boton = tk.Button(self.pantalla_actual, text="Cambiar a Pantalla 1", command=self.seleccionar_secuencias)
        boton.pack()

    def ejecutar(self):
        self.ventana.mainloop()

if __name__ == "__main__":
    interfaz = Interfaz()
    interfaz.ejecutar()
