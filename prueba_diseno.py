import tkinter as tk
from tkinter import ttk

class ScrollableTab(tk.Frame):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)
        
        # Crear un Notebook para las pestañas
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill='both', expand=True)
        
        # Agregar pestañas de ejemplo
        for i in range(5):
            tab_content = tk.Frame(self.notebook)
            self.notebook.add(tab_content, text=f'Tab {i + 1}')
            
            # Crear un Canvas y una barra de desplazamiento para cada pestaña
            canvas = tk.Canvas(tab_content)
            y_scrollbar = tk.Scrollbar(tab_content, orient='vertical', command=canvas.yview)
            y_scrollbar.pack(side='right', fill='y')
            
            # Configurar el Canvas y la barra de desplazamiento
            canvas.configure(yscrollcommand=y_scrollbar.set)
            canvas.pack(side='left', fill='both', expand=True)
            
            # Agregar contenido a cada pestaña (puedes personalizar esto)
            label = tk.Label(canvas, text=f"Contenido de la pestaña {i + 1}")
            canvas.create_window((0, 0), window=label, anchor='nw')
            
            # Configurar el desplazamiento del Canvas
            canvas.update_idletasks()
            canvas.config(scrollregion=canvas.bbox('all'))

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Interfaz con Tabview Scrollable")

    # Crear una instancia de ScrollableTab
    scrollable_tab = ScrollableTab(root)
    scrollable_tab.pack(fill='both', expand=True)

    root.mainloop()
