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
from nucleotidesModule.aligners.aligner import Aligner
from nucleotidesModule.assemblers.assembler import Assembler
from nucleotidesModule.genePredictor.predictor import Predictor
from nucleotidesModule.genomeBrowsers.browser import Browser
from nucleotidesModule.orfFinder.orfFinder import OrfFinder
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
        # self.labelPrincipal = ttk.Label(self, text="Panel principal:", font=self.titleFont)
        # self.labelPrincipal.place(x=10, y=10)

        # label secundario
        # self.labelSecond = ttk.Label(self, text="Desde este panel puede revisar
        # y chequear el progreso de su proyecto:")
        # self.labelSecond.place(x=10, y=40)

        # Creacion de Tabs
        tabControl = ttk.Notebook(self)

        if self.proyect.type == "Análisis de ADN":

            if proyect.steps[0] == "Ensamblaje y alineamiento":
                # Creacion de objetos:
                self.alignTool = Aligner(self.proyect)
                self.assembleTool = Assembler(self.proyect.sequenceRoute, self.proyect.dirRoute + "/Ensamblaje")
                self.visualTool = Browser(self.proyect)

                # Implantación de los pasos:
                self.tab1 = ttk.Frame(tabControl)
                tabControl.add(self.tab1, text="Ensamblaje")

                self.labelStep = ttk.Label(self.tab1, text="Ensamblaje", font=self.titleFont).place(x=10, y=50)

                self.labelDescription = ttk.Label(self.tab1,
                                                  text="Desde este panel puedes revisar y ejecutar el proceso de Ensamblaje").place(
                    x=10, y=90)

                self.buttonAssemble = ttk.Button(self.tab1, text="Evaluar", command=self.evaluateAssemble)
                self.buttonAssemble.place(x=10, y=110)

                # self.configButtonAlign = ttk.Button(self.tab1, text="Configurar", command=self.configAlign)
                # .place(x=10, y=150)

                self.buttonReportAssemble = ttk.Button(self.tab1, text="Leer reporte", command=self.readAssembleReport)
                self.buttonReportAssemble.place(x=10, y=180)
                if not self.assembleTool.fileExist():
                    self.buttonReportAssemble.state(["disabled"])

                # -----------------------------------------------------------------------------------------------------

                self.tab2 = ttk.Frame(tabControl)
                tabControl.add(self.tab2, text="Alineamiento")

                self.labelStep = ttk.Label(self.tab2, text="Alineamiento con Genoma de Referencia", font=self.titleFont).place(x=10, y=50)

                self.labelDescription = ttk.Label(self.tab2,
                                                  text="Desde este panel puedes revisar y ejecutar el proceso de Alineamiento").place(
                    x=10, y=90)

                self.buttonAlign = ttk.Button(self.tab2, text="Evaluar", command=self.evaluateAlign)
                self.buttonAlign.place(x=10, y=110)
                if not self.assembleTool.fileExist():
                    self.buttonAlign.state(["disabled"])

                # self.configAssemble = ttk.Button(self.tab2, text="Configurar").place(x=10, y=150)

                self.buttonReportAlign = ttk.Button(self.tab2, text="Leer reporte",
                                                         command=self.readAlignReport)
                if not self.alignTool.fileExist():
                    self.buttonReportAlign.state(["disabled"])
                self.buttonReportAlign.place(x=10, y=180)

                # -----------------------------------------------------------------------------------------------------

                self.tab3 = ttk.Frame(tabControl)
                tabControl.add(self.tab3, text="Visualizacion")

                self.labelStep = ttk.Label(self.tab3, text="Visualizacion", font=self.titleFont).place(x=10, y=50)

                self.labelDescription = ttk.Label(self.tab3,
                                                  text="Desde este panel puedes revisar y ejecutar el proceso de visualizacion").place(
                    x=10, y=90)

                self.buttonVisual = ttk.Button(self.tab3, text="Abrir Visualizador (Kablammo)", command=self.evaluateVisual)
                self.buttonVisual.place(x=10, y=110)
                if not self.alignTool.fileExist():
                    self.buttonVisual.state(["disabled"])

                # self.progress = ttk.Button(self.tab3, text="Configurar").place(x=10, y=150)

                # self.buttonReportAlign = ttk.Button(self.tab3, text="Leer reporte", command=self.readAlignReport,
                #                                    state="disabled")
                # self.buttonReportAlign.place(x=10, y=180)

            if proyect.steps[0] == "Ensamblaje, homologia":
                # Creacion de objetos:
                self.alignTool = Aligner(self.proyect)
                self.assembleTool = Assembler(self.proyect.sequenceRoute, self.proyect.dirRoute + "/Ensamblaje")
                self.visualTool = Browser(self.proyect)
                self.orfTool = OrfFinder(self.proyect)

                # Implantación de los pasos:
                self.tab1 = ttk.Frame(tabControl)
                tabControl.add(self.tab1, text="Ensamblaje")

                self.labelStep = ttk.Label(self.tab1, text="Ensamblaje", font=self.titleFont).place(x=10, y=50)

                self.labelDescription = ttk.Label(self.tab1,
                                                  text="Desde este panel puedes revisar y ejecutar el proceso de Ensamblaje").place(
                    x=10, y=90)

                self.buttonAssemble = ttk.Button(self.tab1, text="Evaluar", command=self.evaluateAssemble)
                self.buttonAssemble.place(x=10, y=110)

                # self.configButtonAlign = ttk.Button(self.tab1, text="Configurar", command=self.configAlign).place(x=10, y=150)

                self.buttonReportAssemble = ttk.Button(self.tab1, text="Leer reporte", command=self.readAssembleReport)
                self.buttonReportAssemble.place(x=10, y=180)
                if not self.assembleTool.fileExist():
                    self.buttonReportAssemble.state(["disabled"])

                # -----------------------------------------------------------------------------------------------------

                self.tab3 = ttk.Frame(tabControl)
                tabControl.add(self.tab3, text="ORFoma")

                self.labelStep = ttk.Label(self.tab3, text="ORFoma", font=self.titleFont).place(x=10, y=50)

                self.labelDescription = ttk.Label(self.tab3,
                                                  text="Desde este panel puedes revisar y ejecutar el proceso de encontrar ORFs").place(
                    x=10, y=90)

                self.buttonORF = ttk.Button(self.tab3, text="Evaluar", command=self.evaluateORFoma)
                self.buttonORF.place(x=10, y=110)
                if not self.assembleTool.fileExist():
                    self.buttonORF.state(["disabled"])

                # self.progress = ttk.Button(self.tab3, text="Configurar").place(x=10, y=150)

                self.buttonReportORF = ttk.Button(self.tab3, text="Leer reporte", command=self.readORFReport)
                self.buttonReportORF.place(x=10, y=180)
                if not self.orfTool.fileExist():
                    self.buttonReportORF.state(["disabled"])

                # -----------------------------------------------------------------------------------------------------

                self.tab3 = ttk.Frame(tabControl)
                tabControl.add(self.tab3, text="Homologia")

                self.labelStep = ttk.Label(self.tab3, text="Homologia", font=self.titleFont).place(x=10, y=50)

                self.labelDescription = ttk.Label(self.tab3,
                                                  text="Desde este panel puedes revisar y ejecutar el proceso de Homologia").place(
                    x=10, y=90)

                self.buttonAlign = ttk.Button(self.tab3, text="Evaluar", command=self.evaluateAlign2)
                self.buttonAlign.place(x=10, y=110)
                if not self.orfTool.fileExist():
                    self.buttonAlign.state(["disabled"])

                # self.progress = ttk.Button(self.tab3, text="Configurar").place(x=10, y=150)

                self.buttonReportAlign = ttk.Button(self.tab3, text="Leer reporte", command=self.readAlignReport)
                self.buttonReportAlign.place(x=10, y=180)
                if not self.alignTool.fileExist():
                    self.buttonReportAlign.state(["disabled"])

                # -----------------------------------------------------------------------------------------------------

                self.tab3 = ttk.Frame(tabControl)
                tabControl.add(self.tab3, text="Visualizacion")

                self.labelStep = ttk.Label(self.tab3, text="Visualizacion", font=self.titleFont).place(x=10, y=50)

                self.labelDescription = ttk.Label(self.tab3,
                                                  text="Desde este panel puedes revisar y ejecutar el proceso de visualizacion").place(
                    x=10, y=90)

                self.buttonVisual = ttk.Button(self.tab3, text="Abrir Visualizador (Kablammo)", command=self.evaluateVisual)
                self.buttonVisual.place(x=10, y=110)
                if not self.alignTool.fileExist():
                    self.buttonVisual.state(["disabled"])

            if proyect.steps[0] == "otro":
                # Creacion de objetos:
                self.alignTool = Aligner(self.proyect)
                self.assembleTool = Assembler(self.proyect.sequenceRoute, self.proyect.dirRoute + "/Ensamblaje")
                self.visualTool = Browser(self.proyect)
                self.predictionTool = Predictor(self.proyect)

                # Implantación de los pasos:
                self.tab1 = ttk.Frame(tabControl)
                tabControl.add(self.tab1, text="Ensamblaje")

                self.labelStep = ttk.Label(self.tab1, text="Ensamblaje", font=self.titleFont).place(x=10, y=50)

                self.labelDescription = ttk.Label(self.tab1,
                                                  text="Desde este panel puedes revisar y "
                                                       "ejecutar el proceso de Ensamblaje").place(
                    x=10, y=90)

                self.buttonAssemble = ttk.Button(self.tab1, text="Evaluar", command=self.evaluateAssemble)
                self.buttonAssemble.place(x=10, y=110)

                # self.configButtonAlign = ttk.Button(self.tab1, text="Configurar", command=self.configAlign)
                # .place(x=10, y=150)

                self.buttonReportAssemble = ttk.Button(self.tab1, text="Leer reporte", command=self.readAssembleReport)
                self.buttonReportAssemble.place(x=10, y=180)
                if not self.assembleTool.fileExist():
                    self.buttonReportAssemble.state(["disabled"])

                # -----------------------------------------------------------------------------------------------------

                self.tab3 = ttk.Frame(tabControl)
                tabControl.add(self.tab3, text="Prediccion")

                self.labelStep = ttk.Label(self.tab3, text="Prediccion", font=self.titleFont).place(x=10, y=50)

                self.labelDescription = ttk.Label(self.tab3,
                                                  text="Desde este panel puedes revisar y ejecutar"
                                                       " el proceso de Prediccion de genes").place(
                    x=10, y=90)

                self.buttonPrediction = ttk.Button(self.tab3, text="Evaluar", command=self.evaluatePrediction)
                self.buttonPrediction.place(x=10, y=110)
                if not self.assembleTool.fileExist():
                    self.buttonPrediction.state(["disabled"])

                # self.progress = ttk.Button(self.tab3, text="Configurar").place(x=10, y=150)

                self.buttonReportPrediction = ttk.Button(self.tab3, text="Leer reporte",
                                                         command=self.readPredictionReport)
                self.buttonReportPrediction.place(x=10, y=180)
                if not self.predictionTool.fileExist():
                    self.buttonReportPrediction.state(["disabled"])

                # -----------------------------------------------------------------------------------------------------

                self.tab3 = ttk.Frame(tabControl)
                tabControl.add(self.tab3, text="Homologia")

                self.labelStep = ttk.Label(self.tab3, text="Homologia", font=self.titleFont).place(x=10, y=50)

                self.labelDescription = ttk.Label(self.tab3,
                                                  text="Desde este panel puedes revisar y "
                                                       "ejecutar el proceso de Homologia").place(
                    x=10, y=90)

                self.buttonAlign = ttk.Button(self.tab3, text="Evaluar", command=self.evaluateAlign3)
                self.buttonAlign.place(x=10, y=110)
                if not self.predictionTool.fileExist():
                    self.buttonAlign.state(["disabled"])

                # self.progress = ttk.Button(self.tab3, text="Configurar").place(x=10, y=150)

                self.buttonReportAlign = ttk.Button(self.tab3, text="Leer reporte", command=self.readAlignReport)
                self.buttonReportAlign.place(x=10, y=180)
                if not self.alignTool.fileExist():
                    self.buttonReportAlign.state(["disabled"])

                # -----------------------------------------------------------------------------------------------------

                self.tab3 = ttk.Frame(tabControl)
                tabControl.add(self.tab3, text="Visualizacion")

                self.labelStep = ttk.Label(self.tab3, text="Visualizacion", font=self.titleFont).place(x=10, y=50)

                self.labelDescription = ttk.Label(self.tab3,
                                                  text="Desde este panel puedes revisar y "
                                                       "ejecutar el proceso de visualizacion").place(
                    x=10, y=90)

                self.buttonVisual = ttk.Button(self.tab3, text="Abrir Visualizador (Kablammo)",
                                               command=self.evaluateVisual)
                self.buttonVisual.place(x=10, y=110)
                if not self.alignTool.fileExist():
                    self.buttonVisual.state(["disabled"])

        elif self.proyect.type == "Análisis de proteinas":

            if proyect.steps[0] == "Proteina1":
                # Creacion de objetos:

                # Implantación de los pasos:
                self.tab1 = ttk.Frame(tabControl)
                tabControl.add(self.tab1, text="Ensamblaje")

                self.labelStep = ttk.Label(self.tab1, text="Ensamblaje", font=self.titleFont).place(x=10, y=50)

                self.labelDescription = ttk.Label(self.tab1,
                                                  text="Desde este panel puedes revisar y ejecutar el proceso de Ensamblaje").place(
                    x=10, y=90)

                self.buttonAssemble = ttk.Button(self.tab1, text="Evaluar", command=self.evaluateAssemble)
                self.buttonAssemble.place(x=10, y=110)

                # self.configButtonAlign = ttk.Button(self.tab1, text="Configurar", command=self.configAlign)
                # .place(x=10, y=150)

                self.buttonReportAssemble = ttk.Button(self.tab1, text="Leer reporte", command=self.readAssembleReport)
                self.buttonReportAssemble.place(x=10, y=180)
                if not self.assembleTool.fileExist():
                    self.buttonReportAssemble.state(["disabled"])

                self.tab2 = ttk.Frame(tabControl)
                tabControl.add(self.tab2, text="Alineamiento")

                self.labelStep = ttk.Label(self.tab2, text="Alineamiento con Genoma de Referencia", font=self.titleFont).place(x=10, y=50)

                self.labelDescription = ttk.Label(self.tab2,
                                                  text="Desde este panel puedes revisar y ejecutar el proceso de Alineamiento").place(
                    x=10, y=90)

                self.buttonAlign = ttk.Button(self.tab2, text="Evaluar", command=self.evaluateAlign)
                self.buttonAlign.place(x=10, y=110)
                if not self.assembleTool.fileExist():
                    self.buttonAlign.state(["disabled"])

                # self.configAssemble = ttk.Button(self.tab2, text="Configurar").place(x=10, y=150)

                self.buttonReportAlign = ttk.Button(self.tab2, text="Leer reporte",
                                                         command=self.readAlignReport)
                if not self.alignTool.fileExist():
                    self.buttonReportAlign.state(["disabled"])
                self.buttonReportAlign.place(x=10, y=180)

                self.tab3 = ttk.Frame(tabControl)
                tabControl.add(self.tab3, text="Visualizacion")

                self.labelStep = ttk.Label(self.tab3, text="Visualizacion", font=self.titleFont).place(x=10, y=50)

                self.labelDescription = ttk.Label(self.tab3,
                                                  text="Desde este panel puedes revisar y ejecutar el proceso de visualizacion").place(
                    x=10, y=90)

                self.buttonVisual = ttk.Button(self.tab3, text="Abrir Visualizador (Kablammo)", command=self.evaluateVisual)
                self.buttonVisual.place(x=10, y=110)
                if not self.alignTool.fileExist():
                    self.buttonVisual.state(["disabled"])

                # self.progress = ttk.Button(self.tab3, text="Configurar").place(x=10, y=150)

                # self.buttonReportAlign = ttk.Button(self.tab3, text="Leer reporte", command=self.readAlignReport,
                #                                    state="disabled")
                # self.buttonReportAlign.place(x=10, y=180)
            self.tab4 = ttk.Frame(tabControl)
            tabControl.add(self.tab4, text="Predicción")

            self.labelStep = ttk.Label(self.tab4, text="Predicción", font=self.titleFont).place(x=10, y=50)

            self.labelDescription = ttk.Label(self.tab4,
                                              text="Desde este panel puedes revisar y ejecutar el proceso de predicción").place(
                x=10, y=90)

            self.buttonPrediction = ttk.Button(self.tab4, text="Evaluar", command=self.evaluatePrediction)
            self.buttonPrediction.place(x=10, y=110)

            # self.configPrediction = ttk.Button(self.tab4, text="Configurar").place(x=10, y=150)

            self.buttonReportPrediction = ttk.Button(self.tab4, text="Leer reporte", command=self.readPredictionReport,
                                                     state="disabled")
            self.buttonReportPrediction.place(x=10, y=180)

            self.tab5 = ttk.Frame(tabControl)
            tabControl.add(self.tab5, text="Filogenia")

            self.labelStep = ttk.Label(self.tab5, text="Filogenia", font=self.titleFont).place(x=10, y=50)

            self.labelDescription = ttk.Label(self.tab5,
                                              text="Desde este panel puedes revisar y ejecutar el proceso de filogenia").place(
                x=10, y=90)

            self.button = ttk.Button(self.tab5, text="Evaluar").place(x=10, y=110)

            self.progress = ttk.Button(self.tab5, text="Configurar").place(x=10, y=150)

            self.buttonReport = ttk.Button(self.tab5, text="Leer reporte").place(x=10, y=180)

            self.tab6 = ttk.Frame(tabControl)
            tabControl.add(self.tab6, text="Reporte Final")

            self.labelStep = ttk.Label(self.tab6, text="Progreso general", font=self.titleFont).place(x=10, y=50)

            self.progress = ttk.Progressbar(self.tab6).place(x=10, y=110, width=450)

            self.buttonReport = ttk.Button(self.tab6, text="Leer reporte").place(x=10, y=150)

        # Ajustes de pantalla
        tabControl.pack(expan=1, fill="both")
        screen_size = "600x300"

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

        if self.proyect.type == "Análisis de ADN":

            if self.proyect.steps[0] == "Ensamblaje y alineamiento":

                dirRoute = filedialog.asksaveasfilename()
                os.mkdir(dirRoute)

                steps = ["Ensamblaje y alineamiento"]
                os.mkdir(dirRoute+"/Ensamblaje")
                os.mkdir(dirRoute+"/Ensamblaje/Reads")
                steps += [dirRoute+"/Ensamblaje"]

                os.mkdir(dirRoute+"/Homologia")
                os.mkdir(dirRoute+"/Homologia/GR")
                steps += [dirRoute+"/Homologia"]

                if self.proyect.sequenceRoute != "" and self.proyect.genomeRefRoute != "":
                    fileName = self.proyect.sequenceRoute.split("/")[-1]
                    os.rename(self.proyect.sequenceRoute, dirRoute+"/Ensamblaje/Reads/" + fileName)
                    refGenomaFileName = self.proyect.genomeRefRoute.split("/")[-1]
                    os.rename(self.proyect.genomeRefRoute, dirRoute+"/Homologia/GR/"+refGenomaFileName)
                    proyect = Proyect(self.proyect.type, steps, dirRoute + "/archivo.bin", dirRoute+"/Ensamblaje/Reads/"
                                      + fileName, dirRoute+"/Homologia/GR/"+refGenomaFileName, dirRoute)
                    self.panel_window.destroy()
                    new_window = tk.Tk()
                    panelwindow = PanelWindow(new_window, proyect, self.main_window)
                    panelwindow.mainloop()

                else:
                    self.top = tk.Toplevel(self.panel_window)
                    self.top.title("Alerta")
                    tk.Label(self.top,
                             text="No ha seleccionado un archivo para analizar, por favor seleccione un archivo").grid(
                        row=0,
                        column=0,
                        columnspan=2)
                    self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
                    self.button2.grid(row=1, column=0, padx=5, pady=5)

            elif self.proyect.steps[0] == "Ensamblaje, homologia":

                dirRoute = filedialog.asksaveasfilename()
                os.mkdir(dirRoute)

                steps = ["Ensamblaje, homologia"]
                os.mkdir(dirRoute + "/Ensamblaje")
                os.mkdir(dirRoute+"/Ensamblaje/Reads")
                steps += [dirRoute + "/Ensamblaje"]

                os.mkdir(dirRoute+"/ORFoma")
                steps += [dirRoute + "/ORFoma"]

                os.mkdir(dirRoute + "/Homologia")
                os.mkdir(dirRoute+"/Homologia/GR")
                steps += [dirRoute + "/Homologia"]

                if self.proyect.sequenceRoute != "" and self.proyect.genomeRefRoute != "":
                    fileName = self.proyect.sequenceRoute.split("/")[-1]
                    os.rename(self.proyect.sequenceRoute, dirRoute + "/Ensamblaje/Reads/" + fileName)
                    refGenomaFileName = self.proyect.genomeRefRoute.split("/")[-1]
                    os.rename(self.proyect.genomeRefRoute, dirRoute + "/Homologia/GR/" + refGenomaFileName)
                    proyect = Proyect(self.proyect.type, steps, dirRoute + "/archivo.bin",
                                      dirRoute + "/Ensamblaje/Reads/"
                                      + fileName, dirRoute + "/Homologia/GR/" + refGenomaFileName, dirRoute)
                    self.panel_window.destroy()
                    new_window = tk.Tk()
                    panelwindow = PanelWindow(new_window, proyect, self.main_window)
                    panelwindow.mainloop()

                else:
                    self.top = tk.Toplevel(self.panel_window)
                    self.top.title("Alerta")
                    tk.Label(self.top,
                             text="No ha seleccionado un archivo para analizar, por favor seleccione un archivo").grid(
                        row=0,
                        column=0,
                        columnspan=2)
                    self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
                    self.button2.grid(row=1, column=0, padx=5, pady=5)

            elif self.proyect.steps[0] == "otro":

                dirRoute = filedialog.asksaveasfilename()
                os.mkdir(dirRoute)

                steps = ["otro"]
                os.mkdir(dirRoute + "/Ensamblaje")
                os.mkdir(dirRoute+"/Ensamblaje/Reads")
                steps += [dirRoute + "/Ensamblaje"]

                os.mkdir(dirRoute + "/Prediccion")
                steps += [dirRoute + "/Prediccion"]

                os.mkdir(dirRoute + "/Homologia")
                os.mkdir(dirRoute+"/Homologia/GR")
                steps += [dirRoute + "/Homologia"]

                if self.proyect.sequenceRoute != "" and self.proyect.genomeRefRoute != "":
                    fileName = self.proyect.sequenceRoute.split("/")[-1]
                    os.rename(self.proyect.sequenceRoute, dirRoute + "/Ensamblaje/Reads/" + fileName)
                    refGenomaFileName = self.proyect.genomeRefRoute.split("/")[-1]
                    os.rename(self.proyect.genomeRefRoute, dirRoute + "/Homologia/GR/" + refGenomaFileName)
                    proyect = Proyect(self.proyect.type, steps, dirRoute + "/archivo.bin",
                                      dirRoute + "/Ensamblaje/Reads/"
                                      + fileName, dirRoute + "/Homologia/GR/" + refGenomaFileName, dirRoute)
                    self.panel_window.destroy()
                    new_window = tk.Tk()
                    panelwindow = PanelWindow(new_window, proyect, self.main_window)
                    panelwindow.mainloop()

                else:
                    self.top = tk.Toplevel(self.panel_window)
                    self.top.title("Alerta")
                    tk.Label(self.top,
                             text="No ha seleccionado un archivo para analizar, por favor seleccione un archivo").grid(
                        row=0,
                        column=0,
                        columnspan=2)
                    self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
                    self.button2.grid(row=1, column=0, padx=5, pady=5)

        if self.proyect.type == "Análisis de proteinas":

            if self.proyect.steps[0] == "Proteina1":

                dirRoute = filedialog.asksaveasfilename()
                os.mkdir(dirRoute)

                steps = ["Proteina1"]
                os.mkdir(dirRoute+"/Ensamblaje")
                steps += [dirRoute+"/Ensamblaje"]

                os.mkdir(dirRoute+"/Prediccion")
                steps += [dirRoute+"/Prediccion"]

                os.mkdir(dirRoute+"/Homologia")
                steps += [dirRoute+"/Homologia"]

                if self.proyect.sequenceRoute != "" and self.proyect.genomeRefRoute != "":
                    fileName = self.proyect.sequenceRoute.split("/")[-1]
                    os.rename(self.proyect.sequenceRoute, dirRoute + "/Ensamblaje/Reads/" + fileName)
                    refGenomaFileName = self.proyect.genomeRefRoute.split("/")[-1]
                    os.rename(self.proyect.genomeRefRoute, dirRoute + "/Homologia/GR/" + refGenomaFileName)
                    proyect = Proyect(self.proyect.type, steps, dirRoute + "/archivo.bin",
                                      dirRoute + "/Ensamblaje/Reads/"
                                      + fileName, dirRoute + "/Homologia/GR/" + refGenomaFileName, dirRoute)
                    self.panel_window.destroy()
                    new_window = tk.Tk()
                    panelwindow = PanelWindow(new_window, proyect, self.main_window)
                    panelwindow.mainloop()

                else:
                    self.top = tk.Toplevel(self.panel_window)
                    self.top.title("Alerta")
                    tk.Label(self.top,
                             text="No ha seleccionado un archivo para analizar, por favor seleccione un archivo").grid(
                        row=0,
                        column=0,
                        columnspan=2)
                    self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
                    self.button2.grid(row=1, column=0, padx=5, pady=5)

            elif self.proyect.steps[0] == "Proteina2":

                dirRoute = filedialog.asksaveasfilename()
                os.mkdir(dirRoute)

                steps = ["Proteina2"]
                os.mkdir(dirRoute + "/Ensamblaje")
                steps += [dirRoute + "/Ensamblaje"]

                os.mkdir(dirRoute + "/Homologia")
                steps += [dirRoute + "/Homologia"]

                if self.proyect.sequenceRoute != "" and self.proyect.genomeRefRoute != "":
                    fileName = self.proyect.sequenceRoute.split("/")[-1]
                    os.rename(self.proyect.sequenceRoute, dirRoute + "/Ensamblaje/Reads/" + fileName)
                    refGenomaFileName = self.proyect.genomeRefRoute.split("/")[-1]
                    os.rename(self.proyect.genomeRefRoute, dirRoute + "/Homologia/GR/" + refGenomaFileName)
                    proyect = Proyect(self.proyect.type, steps, dirRoute + "/archivo.bin",
                                      dirRoute + "/Ensamblaje/Reads/"
                                      + fileName, dirRoute + "/Homologia/GR/" + refGenomaFileName, dirRoute)
                    self.panel_window.destroy()
                    new_window = tk.Tk()
                    panelwindow = PanelWindow(new_window, proyect, self.main_window)
                    panelwindow.mainloop()

                else:
                    self.top = tk.Toplevel(self.panel_window)
                    self.top.title("Alerta")
                    tk.Label(self.top,
                             text="No ha seleccionado un archivo para analizar, por favor seleccione un archivo").grid(
                        row=0,
                        column=0,
                        columnspan=2)
                    self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
                    self.button2.grid(row=1, column=0, padx=5, pady=5)

            elif self.proyect.steps[0] == "Proteina3":

                dirRoute = filedialog.asksaveasfilename()
                os.mkdir(dirRoute)

                steps = ["Proteina3"]
                os.mkdir(dirRoute + "/Ensamblaje")
                steps += [dirRoute + "/Ensamblaje"]

                os.mkdir(dirRoute + "/Prediccion")
                steps += [dirRoute + "/Prediccion"]

                os.mkdir(dirRoute + "/Homologia")
                steps += [dirRoute + "/Homologia"]

                if self.proyect.sequenceRoute != "" and self.proyect.genomeRefRoute != "":
                    fileName = self.proyect.sequenceRoute.split("/")[-1]
                    os.rename(self.proyect.sequenceRoute, dirRoute + "/Ensamblaje/Reads/" + fileName)
                    refGenomaFileName = self.proyect.genomeRefRoute.split("/")[-1]
                    os.rename(self.proyect.genomeRefRoute, dirRoute + "/Homologia/GR/" + refGenomaFileName)
                    proyect = Proyect(self.proyect.type, steps, dirRoute + "/archivo.bin",
                                      dirRoute + "/Ensamblaje/Reads/"
                                      + fileName, dirRoute + "/Homologia/GR/" + refGenomaFileName, dirRoute)
                    self.panel_window.destroy()
                    new_window = tk.Tk()
                    panelwindow = PanelWindow(new_window, proyect, self.main_window)
                    panelwindow.mainloop()

                else:
                    self.top = tk.Toplevel(self.panel_window)
                    self.top.title("Alerta")
                    tk.Label(self.top,
                             text="No ha seleccionado un archivo para analizar, por favor seleccione un archivo").grid(
                        row=0,
                        column=0,
                        columnspan=2)
                    self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
                    self.button2.grid(row=1, column=0, padx=5, pady=5)

        else:
            self.top = tk.Toplevel(self.panel_window)
            self.top.title("Alerta")
            tk.Label(self.top,
                     text="Porfavor seleccione un solo grupo de pasos").grid(row=0, column=0, columnspan=2)
            self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
            self.button2.grid(row=1, column=0, padx=5, pady=5)

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
        url = "https://bioin.innovarweb.com/"
        webbrowser.open(url)

    def openWebHelp(self):
        url = "https://bioin.innovarweb.com/"
        webbrowser.open(url)

    def evaluateAlign(self):
        self.alignTool.ejecutCommand()
        if self.alignTool.fileExist():
            self.buttonReportAlign.state(["!disabled"])
            self.buttonVisual.state(["!disabled"])

    def evaluateAlign2(self):
        self.alignTool.ejecutCommand2()
        if self.alignTool.fileExist():
            self.buttonReportAlign.state(["!disabled"])
            self.buttonVisual.state(["!disabled"])

    def evaluateAlign3(self):
        self.alignTool.ejecutCommand3()
        if self.alignTool.fileExist():
            self.buttonReportAlign.state(["!disabled"])
            self.buttonVisual.state(["!disabled"])

    def configAlign(self):
        return

    def readAlignReport(self):
        self.topAlign = tk.Toplevel(self.panel_window)
        self.topAlign.title("Reporte Alineamiento")

        self.alignerTextBox = tk.Text(self.topAlign, height=100, width=100)
        self.alignerTextBox.pack()
        f = open(self.proyect.dirRoute + "/Homologia/output.xml", "r")
        self.alignerTextBox.insert(tk.INSERT, f.buffer.read())

    def evaluateAssemble(self):
        self.assembleTool.ejecutCommand()
        if self.assembleTool.fileExist():
            self.buttonReportAssemble.state(["!disabled"])
            self.buttonAlign.state(["!disabled"])

    def configAssemble(self):
        return

    def readAssembleReport(self):
        self.topAssemble = tk.Toplevel(self.panel_window)
        self.topAssemble.title("Reporte Ensamblaje")

        self.assembleTextBox = tk.Text(self.topAssemble, height=100, width=100)
        self.assembleTextBox.pack()
        f = open(self.proyect.dirRoute + "/Ensamblaje/output_assembly/output_d_results/output_out.padded.fasta", "r")
        self.assembleTextBox.insert(tk.INSERT, f.buffer.read())

    def evaluateVisual(self):
        self.visualTool.ejecutCommand()

    def readGBrowserReport(self):
        return

    def evaluatePrediction(self):
        self.predictionTool.ejecutCommand()
        if self.predictionTool.fileExist():
            self.buttonReportPrediction.state(["!disabled"])

    def configPrediction(self):
        return

    def readPredictionReport(self):
        self.topPredictor = tk.Toplevel(self.panel_window)
        self.topPredictor.title("Reporte Prediccion")

        self.predictorTextBox = tk.Text(self.topPredictor, height=100, width=100)
        self.predictorTextBox.pack()
        f = open(self.proyect.dirRoute + "/Prediccion/output.fasta", "r")
        self.predictorTextBox.insert(tk.INSERT, f.buffer.read())

    def evaluateORFoma(self):
        self.orfTool.ejecutCommand()
        if self.orfTool.fileExist():
            self.buttonReportORF.state(["!disabled"])
            self.buttonAlign.state(["!disabled"])

    def readORFReport(self):
        self.topORF = tk.Toplevel(self.panel_window)
        self.topORF.title("Reporte ORFoma")

        self.orfTextBox = tk.Text(self.topORF, height=100, width=100)
        self.orfTextBox.pack()
        f = open(self.proyect.dirRoute + "/ORFoma/output.fasta", "r")
        self.orfTextBox.insert(tk.INSERT, f.buffer.read())
