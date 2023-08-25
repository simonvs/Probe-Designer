import tkinter as tk

# Función para obtener las opciones seleccionadas
def obtener_seleccion():
    opciones_seleccionadas = []
    for i, var in enumerate(variables):
        if var.get():
            opciones_seleccionadas.append(opciones[i])
    print("Opciones seleccionadas:", opciones_seleccionadas)

# Crear la ventana principal
root = tk.Tk()
root.title("Selección de Opciones")

# Lista de opciones
opciones = ["Opción 1", "Opción 2", "Opción 3", "Opción 4"]

# Crear una lista para mantener los Checkbuttons y sus variables asociadas
checkboxes = []
variables = []

# Crear los Checkbuttons
for opcion in opciones:
    var = tk.BooleanVar()
    checkbox = tk.Checkbutton(root, text=opcion, variable=var)
    checkbox.pack(anchor="w")  # Alineación a la izquierda
    checkboxes.append(checkbox)
    variables.append(var)

# Botón para obtener las opciones seleccionadas
boton_obtener_seleccion = tk.Button(root, text="Obtener Selección", command=obtener_seleccion)
boton_obtener_seleccion.pack()

root.mainloop()
