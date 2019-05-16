#!/usr/bin/python
# Proyecto: BIOIN
# Trabajo de Grado
# Alejandro Valencia
# Juan Jose Varela
# Universidad del Valle

import tkinter as tk
from tkinter import ttk
from tkinter import font
from tkinter import filedialog
from panel import PanelWindow
from buildProyect import BuildProyectWindow
import webbrowser
import pickle


class NewProyectWindow(ttk.Frame):

    def __init__(self, new_proyect_window, main_window):

        # Intancias y configuraciones de la ventana
        self.build_proyect_window = new_proyect_window
        self.main_window = main_window
        super().__init__(self.build_proyect_window)
        self.build_proyect_window.title("BIOIN - NUEVO PROYECTO")
        self.build_proyect_window.geometry("300x250")
        self.place(relwidth=1, relheight=1)
        self.build_proyect_window.resizable(0, 0)
        self.build_proyect_window.protocol("WM_DELETE_WINDOW", self.onClosing)

        # menu de la pantalla de inicio
        self.menubar = tk.Menu(self.build_proyect_window)
        self.fileMenu = tk.Menu(self.menubar, tearoff=0)
        self.fileMenu.add_command(label="Abrir...", command=self.openProyect)

        self.fileMenu.add_separator()

        self.fileMenu.add_command(label="Salir", command=self.onClosing)
        self.menubar.add_cascade(label="Proyecto", menu=self.fileMenu)

        self.helpMenu = tk.Menu(self.menubar, tearoff=0)
        self.helpMenu.add_command(label="Más Ayuda", command=self.openWebHelp)
        self.helpMenu.add_command(label="Acerca de Bioin...", command=self.openWeb)
        self.menubar.add_cascade(label="Ayuda", menu=self.helpMenu)

        self.build_proyect_window.config(menu=self.menubar)

        # definimos fuentes a usar
        self.titleFont = font.Font(family="Times", size=14)

        # Label principal
        self.labelPrincipal = ttk.Label(self, text="Iniciar un nuevo proyecto", font=self.titleFont)
        self.labelPrincipal.place(x=10, y=10)

        # label secundario
        self.labelSecond = ttk.Label(self, text="Seleccione el tipo de proyecto:")
        self.labelSecond.place(x=10, y=40)

        # menu de selección de  tipo de proyecto
        self.selectionMenu = ttk.Combobox(self, state= "readonly", values=["Análisis de ADN", "Análisis de proteinas"])
        self.selectionMenu.current(0)
        self.selectionMenu.place(x=10, y=70)

        # label terciario
        self.labelThird = ttk.Label(self, text="Cargar archivo:")
        self.labelThird.place(x=10, y=120)

        # text para definir la ruta de carga de archivos
        self.routeText = ttk.Entry(self)
        self.routeText.place(x=10, y=150)

        # botón para añadir la ruta por medio de una ventana
        self.routeButton = ttk.Button(self, text="...", width=3, command=self.findFileRoute)
        self.routeButton.place(x=150, y=148)

        # boton de inicio de proyecto
        self.confirmButton = ttk.Button(self, text="Continuar", command=self.onStartProyect)
        self.confirmButton.place(x=10, y=200)

        # boton de regreso
        self.returnButton = ttk.Button(self, text="Regresar", command=self.onClosing)
        self.returnButton.place(x=200, y=200)

    def onStartProyect(self):
        proyectType= self.selectionMenu.get()
        fileRoute= self.routeText.get()
        if proyectType != "" and fileRoute != "":
            self.build_proyect_window.destroy()
            new_window = tk.Tk()
            buildProyect = BuildProyectWindow(new_window, proyectType, fileRoute, self.main_window)
            buildProyect.mainloop()
        else:
            self.top = tk.Toplevel(self.build_proyect_window)
            self.top.title("Alerta")
            tk.Label(self.top, text="No ha seleccionado un archivo para analizar, por favor seleccione un archivo").grid(row=0, column=0, columnspan=2)
            self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
            self.button2.grid(row=1, column=0, padx=5, pady=5)


    def onClosing(self):
        self.top = tk.Toplevel(self.build_proyect_window)
        self.top.title("Salir")

        tk.Label(self.top, text="¿Está seguro?").grid(row=0, column=0, columnspan=2)

        self.button1 = tk.Button(self.top, text="Si, deseo regresar.", command=self.salir)
        self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
        self.button1.grid(row=1, column=0, padx=5, pady=5)
        self.button2.grid(row=1, column=1, padx=5, pady=5)

    def openProyect(self):
        dirRoute = filedialog.askopenfilename()
        if dirRoute != () and dirRoute != '':
            self.build_proyect_window.destroy()
            with open(dirRoute, "br") as archivo:
                proyect = pickle.load(archivo)
            new_window = tk.Tk()
            panelwindow = PanelWindow(new_window, proyect, self.main_window)
            panelwindow.mainloop()

    def findFileRoute(self):
        dirRoute = filedialog.askopenfilename()
        if dirRoute != () and dirRoute != '':
            self.routeText.insert(0, dirRoute)

    def salir(self):
        self.main_window.deiconify()
        self.top.destroy()
        self.build_proyect_window.destroy()

    def cancelar(self):
        self.top.destroy()

    def openWeb(self):
        url = "http://bioinformatica.univalle.edu.co/"
        webbrowser.open(url)

    def openWebHelp(self):
        url = "http://bioinformatica.univalle.edu.co/"
        webbrowser.open(url)
