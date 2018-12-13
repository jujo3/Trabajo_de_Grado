#!/usr/bin/python
# Proyecto: BIOIN
# Trabajo de Grado
# Alejandro Valencia
# Juan Jose Varela
# Universidad del Valle

import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
from tkinter import font
from proyect import Proyect
import webbrowser
import pickle
import os
import sys


class PanelWindow(ttk.Frame):

    def __init__(self, panel_window, proyect, main_window):

        # Intancias y configuraciones de la ventana
        self.panel_window = panel_window
        self.proyect = proyect
        self.main_window = main_window
        super().__init__(self.panel_window)

        # Se guarda el proyecto
        self.saveProyect()

        # menu de la pantalla de inicio
        self.menubar = tk.Menu(self.panel_window)
        self.fileMenu = tk.Menu(self.menubar, tearoff=0)
        self.fileMenu.add_command(label="Crear nuevo proyecto", command=self.newProyect)
        self.fileMenu.add_command(label="Abrir...", command=self.openProyect)
        self.fileMenu.add_command(label="Guardar...", command=self.saveProyect)
        self.fileMenu.add_command(label="Guardar como...", command=self.saveProyectAs)

        self.fileMenu.add_separator()

        self.fileMenu.add_command(label="Salir", command=self.onClosing)
        self.menubar.add_cascade(label="Proyecto", menu=self.fileMenu)

        self.helpMenu = tk.Menu(self.menubar, tearoff=0)
        self.helpMenu.add_command(label="Más Ayuda", command=self.openWebHelp)
        self.helpMenu.add_command(label="Acerca de Bioin...", command=self.openWeb)
        self.menubar.add_cascade(label="Ayuda", menu=self.helpMenu)

        self.panel_window.config(menu=self.menubar)

        # definimos fuentes a usar
        self.titleFont = font.Font(family="Times", size=14)

        # Label principal
        self.labelPrincipal = ttk.Label(self, text="Panel principal:", font=self.titleFont)
        self.labelPrincipal.place(x=10, y=10)

        # label secundario
        self.labelSecond = ttk.Label(self, text="Desde este panel puede revisar y chequear el progreso de su proyecto:")
        self.labelSecond.place(x=10, y=40)

        # variables de entorno
        add = 0

        # Implantación de los pasos:

        if self.proyect.steps_list[0]:
            self.labelStep = ttk.Label(self, text="Paso 1", font=self.titleFont).place(x=10, y=70 + add)

            self.labelDescription = ttk.Label(self, text="Esta es la descripcion del paso 1.").place(x=10, y=90 + add)

            self.button = ttk.Button(self, text="Evaluar").place(x=10, y=110 + add)

            self.progress = ttk.Progressbar(self).place(x=90, y=110 + add, width=200)

            self.buttonReport = ttk.Button(self, text="Leer reporte").place(x=300, y=110 + add)

            add += 100

        if self.proyect.steps_list[1]:
            self.labelStep = ttk.Label(self, text="Paso 2", font=self.titleFont).place(x=10, y=70 + add)

            self.labelDescription = ttk.Label(self, text="Esta es la descripcion del paso 2.").place(x=10, y=90 + add)

            self.button = ttk.Button(self, text="Evaluar").place(x=10, y=110 + add)

            self.progress = ttk.Progressbar(self).place(x=90, y=110 + add, width=200)

            self.buttonReport = ttk.Button(self, text="Leer reporte").place(x=300, y=110 + add)

            add += 100

        if self.proyect.steps_list[2]:
            self.labelStep = ttk.Label(self, text="Paso 3", font=self.titleFont).place(x=10, y=70 + add)

            self.labelDescription = ttk.Label(self, text="Esta es la descripcion del paso 3.").place(x=10, y=90 + add)

            self.button = ttk.Button(self, text="Evaluar").place(x=10, y=110 + add)

            self.progress = ttk.Progressbar(self).place(x=90, y=110 + add, width=200)

            self.buttonReport = ttk.Button(self, text="Leer reporte").place(x=300, y=110 + add)

            add += 100

        if self.proyect.steps_list[3]:
            self.labelStep = ttk.Label(self, text="Paso 4", font=self.titleFont).place(x=10, y=70 + add)

            self.labelDescription = ttk.Label(self, text="Esta es la descripcion del paso 4.").place(x=10, y=90 + add)

            self.button = ttk.Button(self, text="Evaluar").place(x=10, y=110 + add)

            self.progress = ttk.Progressbar(self).place(x=90, y=110 + add, width=200)

            self.buttonReport = ttk.Button(self, text="Leer reporte").place(x=300, y=110 + add)

            add += 100

        if self.proyect.steps_list[4]:
            self.labelStep = ttk.Label(self, text="Paso 5", font=self.titleFont).place(x=10, y=70 + add)

            self.labelDescription = ttk.Label(self, text="Esta es la descripcion del paso 5.").place(x=10, y=90 + add)

            self.button = ttk.Button(self, text="Evaluar").place(x=10, y=110 + add)

            self.progress = ttk.Progressbar(self).place(x=90, y=110 + add, width=200)

            self.buttonReport = ttk.Button(self, text="Leer reporte").place(x=300, y=110 + add)

            add += 100

        self.labelStep = ttk.Label(self, text="Progreso general", font=self.titleFont).place(x=10, y=70 + add)

        self.progress = ttk.Progressbar(self).place(x=10, y=110 + add, width=450)

        self.buttonReport = ttk.Button(self, text="Leer reporte").place(x=500, y=110 + add)

        # Ajustes de pantalla
        add += 150
        screen_size = "800x"+str(add)

        self.panel_window.title("BIOIN - PANEL PRINCIPAL")
        self.panel_window.geometry(screen_size)
        self.place(relwidth=1, relheight=1)
        self.panel_window.resizable(0, 0)
        self.panel_window.protocol("WM_DELETE_WINDOW", self.onClosing)

    def newProyect(self):
        self.top = tk.Toplevel(self.panel_window)
        self.top.title("Iniciar un nuevo proyecto")

        tk.Label(self.top, text="¿Está seguro?").grid(row=0, column=0, columnspan=2)

        self.button1 = tk.Button(self.top, text="Si, deseo iniciar un nuevo proyecto.", command=self.startNewProyect)
        self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
        self.button1.grid(row=1, column=0, padx=5, pady=5)
        self.button2.grid(row=1, column=1, padx=5, pady=5)

    def startNewProyect(self):
        self.main_window.deiconify()
        self.top.destroy()
        self.panel_window.destroy()

    def openProyect(self):
        dirRoute = filedialog.askopenfilename()
        if dirRoute!=():
            self.panel_window.destroy()
            with open(dirRoute, "br") as archivo:
                proyect = pickle.load(archivo)
            new_window = tk.Tk()
            panelwindow = PanelWindow(new_window, proyect, self.main_window)
            panelwindow.mainloop()

    def saveProyect(self):
        with open(self.proyect.route, "bw") as archivo:
            pickle.dump(self.proyect, archivo)

    def saveProyectAs(self):
        dirRoute = filedialog.asksaveasfilename()
        os.mkdir(dirRoute)
        self.proyect.route = dirRoute+"/archivo.bin"
        with open(self.proyect.route, "bw") as archivo:
            pickle.dump(self.proyect, archivo)

    def onClosing(self):

        self.top = tk.Toplevel(self.panel_window)
        self.top.title("Salir")

        tk.Label(self.top, text="¿Está seguro?").grid(row=0, column=0, columnspan=2)

        self.button1 = tk.Button(self.top, text="Si, deseo salir.", command=self.salir)
        self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
        self.button1.grid(row=1, column=0, padx=5, pady=5)
        self.button2.grid(row=1, column=1, padx=5, pady=5)

    def salir(self):
        self.top.destroy()
        self.panel_window.destroy()
        sys.exit()

    def cancelar(self):
        self.top.destroy()

    def openWeb(self):
        url = "http://bioinformatica.univalle.edu.co/"
        webbrowser.open(url)

    def openWebHelp(self):
        url = "http://bioinformatica.univalle.edu.co/"
        webbrowser.open(url)
