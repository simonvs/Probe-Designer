import tkinter as tk


class SeleccionarOpciones:
    def __init__(self, root):
        super().__init__(root)
        self.root = root
        self.root.title("Seleccionar Opciones")

        self.opciones_disponibles = ["Opción 1", "Opción 2", "Opción 3", "Opción 4", "Opción 5"]
        self.opciones_seleccionadas = []

        self.lista_disponibles = tk.Listbox(self, selectmode=tk.SINGLE)
        self.lista_disponibles.pack(side=tk.LEFT, padx=10, pady=10)

        for opcion in self.opciones_disponibles:
            self.lista_disponibles.insert(tk.END, opcion)

        self.lista_seleccionadas = tk.Listbox(self, selectmode=tk.SINGLE)
        self.lista_seleccionadas.pack(side=tk.RIGHT, padx=10, pady=10)

        self.boton_seleccionar = tk.Button(self, text="Seleccionar", command=self.seleccionar)
        self.boton_seleccionar.pack()
        
        self.boton_eliminar = tk.Button(self, text="Eliminar", command=self.eliminar)
        self.boton_eliminar.pack()

        self.boton_seleccionar_todo = tk.Button(self, text="Seleccionar Todo", command=self.seleccionar_todo)
        self.boton_seleccionar_todo.pack()

        self.boton_eliminar_todo = tk.Button(self, text="Eliminar Todo", command=self.eliminar_todo)
        self.boton_eliminar_todo.pack()

    def seleccionar(self):
        seleccion = self.lista_disponibles.get(tk.ACTIVE)
        if seleccion not in self.opciones_seleccionadas:
            self.opciones_seleccionadas.append(seleccion)
            self.lista_seleccionadas.insert(tk.END, seleccion)
            self.lista_disponibles.delete(tk.ACTIVE)

    def eliminar(self):
        seleccion = self.lista_seleccionadas.get(tk.ACTIVE)
        if seleccion in self.opciones_seleccionadas:
            self.opciones_seleccionadas.remove(seleccion)
            self.lista_disponibles.insert(tk.END, seleccion)
            self.lista_seleccionadas.delete(tk.ACTIVE)

    def seleccionar_todo(self):
        for opcion in self.opciones_disponibles:
            if opcion not in self.opciones_seleccionadas:
                self.opciones_seleccionadas.append(opcion)
                self.lista_seleccionadas.insert(tk.END, opcion)
        self.lista_disponibles.delete(0, tk.END)

    def eliminar_todo(self):
        for opcion in self.opciones_seleccionadas:
            self.opciones_seleccionadas.remove(opcion)
            self.lista_disponibles.insert(tk.END, opcion)
        self.lista_seleccionadas.delete(0, tk.END)

if __name__ == "__main__":
    root = tk.Tk()
    app = SeleccionarOpciones(root)
    app.pack()
    root.mainloop()
