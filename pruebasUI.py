import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageTk

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

root = tk.Tk()
root.title("Tooltip con Imagen")

label = tk.Label(root, text="Pasa el cursor por encima de mí", cursor="hand2")
label.pack(padx=20, pady=20)

imagen_path = "images/largo.png"  # Reemplaza con la ruta de tu imagen
tooltip = ToolTip(label, "Texto de tooltip aquí", imagen_path)

root.mainloop()
