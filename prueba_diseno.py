import tkinter as tk
from tkinter import ttk

class DataTable:
    def __init__(self, parent):
        self.parent = parent
        self.table = ttk.Treeview(self.parent)
        self.table['columns'] = ('ID', 'Nombre', 'Edad')
        self.table.column("#0", width=0, stretch=tk.NO)
        self.table.column("ID", anchor=tk.W, width=100)
        self.table.column("Nombre", anchor=tk.W, width=100)
        self.table.column("Edad", anchor=tk.W, width=100)
        self.table.heading("ID", text="ID")
        self.table.heading("Nombre", text="Nombre")
        self.table.heading("Edad", text="Edad")
        self.table.grid(row=0, column=0, sticky="nsew")

def resize_frame(event):
    canvas.configure(scrollregion=canvas.bbox('all'))
    canvas_height = frame.winfo_reqheight()
    if canvas_height != canvas.winfo_height():
        canvas.config(height=canvas_height)

root = tk.Tk()
root.geometry("400x300")

# Crear un frame contenedor para la tabla
tabla_frame = tk.Frame(root)
tabla_frame.grid(row=0, column=0, sticky="nsew")

# Crear una tabla dentro del frame
data_table = DataTable(tabla_frame)

# Agregar barras de desplazamiento al frame si es necesario
canvas = tk.Canvas(tabla_frame)
canvas.grid(row=0, column=0, sticky="nsew")
scrollbar = ttk.Scrollbar(tabla_frame, orient=tk.VERTICAL, command=canvas.yview)
scrollbar.grid(row=0, column=1, sticky="ns")
canvas.configure(yscrollcommand=scrollbar.set)

frame = tk.Frame(canvas)
canvas.create_window((0, 0), window=frame, anchor='nw')
frame.bind("<Configure>", resize_frame)

root.mainloop()
