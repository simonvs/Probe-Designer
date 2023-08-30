import tkinter as tk

root = tk.Tk()
root.title("Grid dentro de una Celda Grid")

# Crear un Frame principal con grid
frame_principal = tk.Frame(root)
frame_principal.grid(row=0, column=0)

# Crear un Label en la celda (0, 0) del Frame principal
label1 = tk.Label(frame_principal, text="Label 1")
label1.grid(row=0, column=0)

# Crear un Frame secundario dentro de la celda (0, 1) del Frame principal
frame_secundario = tk.Frame(frame_principal)
frame_secundario.grid(row=0, column=1)

# Dentro del Frame secundario, crear elementos con grid
label2 = tk.Label(frame_secundario, text="Label 2")
label2.grid(row=0, column=0)

entry = tk.Entry(frame_secundario)
entry.grid(row=1, column=0)

button = tk.Button(frame_secundario, text="Botón")
button.grid(row=2, column=0)

root.mainloop()
