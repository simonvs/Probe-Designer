import tkinter as tk
from tkinter import ttk
import sqlite3
import os

from pyparsing import col

# Conectar a la base de datos SQLite3 (creará la base de datos si no existe)
conn = sqlite3.connect(os.path.join(os.getcwd(),"databases", "personas.db"))
cursor = conn.cursor()

# Crear una tabla de ejemplo
cursor.execute('''CREATE TABLE IF NOT EXISTS personas
                  (id INTEGER PRIMARY KEY AUTOINCREMENT, nombre TEXT, edad INTEGER)''')
conn.commit()

# Insertar algunos datos de ejemplo
cursor.executemany("INSERT INTO personas (nombre, edad) VALUES (?, ?)",
                   [("Alice", 25), ("Bob", 30), ("Charlie", 22)])
conn.commit()

# Función para mostrar la edad
def mostrar_edad(edad):
    print(f"La edad es: {edad}")

# Función para eliminar un registro
def eliminar_registro(id_registro):
    cursor.execute("DELETE FROM personas WHERE id=?", (id_registro,))
    conn.commit()
    actualizar_tabla()

def actualizar_tabla():
    # Limpiar la tabla antes de volver a cargar los datos
    for row in tree.get_children():
        tree.delete(row)

    # Obtener los datos de la base de datos y mostrarlos en la tabla
    cursor.execute("SELECT * FROM personas")
    for row in cursor.fetchall():
        # Para cada registro, agregar una fila en la tabla con dos botones
        index = tree.insert("", "end", values=(row[0], row[1], row[2]))
        # Botón para mostrar la edad
        btn_mostrar = tk.Button(root, text="Mostrar Edad", command=lambda edad=row[2]: mostrar_edad(edad))
        btn_mostrar.grid(row=tree.index(index), column=3)
        # Botón para eliminar el registro
        btn_eliminar = tk.Button(root, text="Eliminar", command=lambda id_registro=row[0]: eliminar_registro(id_registro))
        btn_eliminar.grid(row=tree.index(index), column=4)
    for col in tree["columns"]:
        tree.column(col, width=80)
# Configurar la ventana principal
root = tk.Tk()
root.title("Registros de Personas")

# Crear el Treeview para mostrar la tabla
tree = ttk.Treeview(root, columns=("ID", "Nombre", "Edad"), show="headings")
tree.heading("ID", text="ID")
tree.heading("Nombre", text="Nombre")
tree.heading("Edad", text="Edad")
tree.grid(row=0,column=0,pady=10)

# Función para ajustar dinámicamente la altura de las filas
def ajustar_altura(event):
    for row_id in tree.get_children():
        tree.item(row_id, tags="")
        tree.update_idletasks()  # Asegurar que las etiquetas se apliquen
        # Obtener la altura del contenido de la columna "Edad"
        altura_contenido = tree.bbox(row_id, column="Edad")[-1]

        # Obtener la altura de la fila
        altura_fila = tree.bbox(row_id, column="#0")[-1]

        # Configurar la etiqueta de la fila con la altura adecuada
        tree.tag_configure(row_id, height=altura_contenido + altura_fila)


# ...

# Enlazar la función para ajustar la altura al evento <Configure>
tree.bind("<Configure>", ajustar_altura)


# Obtener los datos de la base de datos y mostrarlos en la tabla
cursor.execute("SELECT * FROM personas")
for row in cursor.fetchall():
    # Para cada registro, agregar una fila en la tabla con dos botones
    index = tree.insert("", "end", values=(row[0], row[1], row[2]))
    
    # Botón para mostrar la edad
    btn_mostrar = tk.Button(root, text="Mostrar Edad", command=lambda edad=row[2]: mostrar_edad(edad))
    btn_mostrar.grid(row=tree.index(index), column=3)
    # Botón para eliminar el registro
    btn_eliminar = tk.Button(root, text="Eliminar", command=lambda id_registro=row[0]: eliminar_registro(id_registro))
    btn_eliminar.grid(row=tree.index(index), column=4)

# Configurar el evento de clic en la tabla
tree.bind("<Button-1>", lambda event: btn_mostrar.focus_set())

# Iniciar el bucle principal
root.mainloop()

# Cerrar la conexión a la base de datos al finalizar el programa
conn.close()
