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
from step import Step
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

        # Creacion de Tabs
        tabControl = ttk.Notebook(self)

        # Implantación de los pasos:
        if self.proyect.steps_list[0]:

            self.tab1 = ttk.Frame(tabControl)
            tabControl.add(self.tab1, text="Ensamblaje")

            self.labelStep = ttk.Label(self.tab1, text="Ensamblaje", font=self.titleFont).place(x=10, y=50)

            self.labelDescription = ttk.Label(self.tab1, text="Esta es la descripcion del paso 1.").place(x=10, y=90)

            self.button = ttk.Button(self.tab1, text="Evaluar").place(x=10, y=110)

            self.progress = ttk.Progressbar(self.tab1).place(x=90, y=110, width=200)

            self.buttonReport = ttk.Button(self.tab1, text="Leer reporte").place(x=300, y=110)

        if self.proyect.steps_list[1]:

            self.tab2 = ttk.Frame(tabControl)
            tabControl.add(self.tab2, text="Alineamiento")

            self.labelStep = ttk.Label(self.tab2, text="Alineamiento", font=self.titleFont).place(x=10, y=50)

            self.labelDescription = ttk.Label(self.tab2, text="Esta es la descripcion del paso 2.").place(x=10, y=90)

            self.button = ttk.Button(self.tab2, text="Evaluar").place(x=10, y=110)

            self.progress = ttk.Progressbar(self.tab2).place(x=90, y=110, width=200)

            self.buttonReport = ttk.Button(self.tab2, text="Leer reporte").place(x=300, y=110)

        if self.proyect.steps_list[2]:

            self.tab3 = ttk.Frame(tabControl)
            tabControl.add(self.tab3, text="Predictor")

            self.labelStep = ttk.Label(self.tab3, text="Predictor", font=self.titleFont).place(x=10, y=50)

            self.labelDescription = ttk.Label(self.tab3, text="Esta es la descripcion del paso 3.").place(x=10, y=90)

            self.button = ttk.Button(self.tab3, text="Evaluar").place(x=10, y=110)

            self.progress = ttk.Progressbar(self.tab3).place(x=90, y=110, width=200)

            self.buttonReport = ttk.Button(self.tab3, text="Leer reporte").place(x=300, y=110)

        if self.proyect.steps_list[3]:

            self.tab4 = ttk.Frame(tabControl)
            tabControl.add(self.tab4, text="Genome Browser")

            self.labelStep = ttk.Label(self.tab4, text="Genome Browser", font=self.titleFont).place(x=10, y=50)

            self.labelDescription = ttk.Label(self.tab4, text="Esta es la descripcion del paso 4.").place(x=10, y=90)

            self.button = ttk.Button(self.tab4, text="Evaluar").place(x=10, y=110)

            self.progress = ttk.Progressbar(self.tab4).place(x=90, y=110, width=200)

            self.buttonReport = ttk.Button(self.tab4, text="Leer reporte").place(x=300, y=110)

        if self.proyect.steps_list[4]:

            self.tab5 = ttk.Frame(tabControl)
            tabControl.add(self.tab5, text="Filogenia")

            self.labelStep = ttk.Label(self.tab5, text="Filogenia", font=self.titleFont).place(x=10, y=50)

            self.labelDescription = ttk.Label(self.tab5, text="Esta es la descripcion del paso 5.").place(x=10, y=90)

            self.button = ttk.Button(self.tab5, text="Evaluar").place(x=10, y=110)

            self.progress = ttk.Progressbar(self.tab5).place(x=90, y=110, width=200)

            self.buttonReport = ttk.Button(self.tab5, text="Leer reporte").place(x=300, y=110)

        self.tab6 = ttk.Frame(tabControl)
        tabControl.add(self.tab6, text="Reporte Final")

        self.labelStep = ttk.Label(self.tab6, text="Progreso general", font=self.titleFont).place(x=10, y=50)

        self.progress = ttk.Progressbar(self.tab6).place(x=10, y=110, width=450)

        self.buttonReport = ttk.Button(self.tab6, text="Leer reporte").place(x=500, y=110)

        # Ajustes de pantalla
        tabControl.pack(expan=1, fill="both")
        screen_size = "1024x768"

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
        if dirRoute != () and dirRoute != '':
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
        if dirRoute != () and dirRoute != '':
            os.mkdir(dirRoute)
            self.proyect.steps = []
            if self.proyect.steps_list[0]:
                os.mkdir(dirRoute+"/Ensamblaje")
                self.proyect.steps += [Step("Ensamblaje", "script", "config", dirRoute+"/Ensamblaje")]
            if self.proyect.steps_list[1]:
                os.mkdir(dirRoute+"/Alineamiento")
                self.proyect.steps += [Step("Alineamiento", "script", "config", dirRoute+"/Alineamiento")]
            if self.proyect.steps_list[2]:
                os.mkdir(dirRoute+"/Predictor")
                self.proyect.steps += [Step("Predictor", "script", "config", dirRoute+"/Predictor")]
            if self.proyect.steps_list[3]:
                os.mkdir(dirRoute+"/GenomeBrowser")
                self.proyect.steps += [Step("GenomeBrowser", "script", "config", dirRoute+"/GenomeBrowser")]
            if self.proyect.steps_list[4]:
                os.mkdir(dirRoute+"/Filogenia")
                self.proyect.steps += [Step("Filogenia", "script", "config", dirRoute+"/Filogenia")]
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
