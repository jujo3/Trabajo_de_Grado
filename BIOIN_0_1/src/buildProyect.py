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
from proyect import Proyect
from step import Step
import webbrowser
import os
import pickle


class BuildProyectWindow(ttk.Frame):

    def __init__(self, build_proyect_window, proyectType, fileRoute, file2route, main_window):

        # Intancias y configuraciones de la ventana
        self.build_proyect_window = build_proyect_window
        self.proyectType = proyectType
        self.fileRoute = fileRoute
        self.file2route = file2route
        self.main_window = main_window
        super().__init__(self.build_proyect_window)
        self.build_proyect_window.title("BIOIN - NUEVO PROYECTO")
        self.build_proyect_window.geometry("400x250")
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
        self.labelSecond = ttk.Label(self, text="Seleccione los pasos:")
        self.labelSecond.place(x=10, y=40)

        # booleans:

        if self.proyectType == "Análisis de ADN":
            # check buttons
            self.step1Value = tk.BooleanVar(self)
            self.step1 = ttk.Checkbutton(self, text="Ensamblaje y alineamiento con genoma de referencia", variable=self.step1Value)
            self.step1.place(x=15, y=70)

            # botón para añadir la ruta por medio de una ventana
            # self.configButton1 = ttk.Button(self, text="...", width=3)
            # self.configButton1.place(x=150, y=68)

            self.step2Value = tk.BooleanVar(self)
            self.step2 = ttk.Checkbutton(self, text="Ensamblaje y homologia con respecto a una base de datos", variable=self.step2Value)
            self.step2.place(x=15, y=100)

            # botón para añadir la ruta por medio de una ventana
            # self.configButton2 = ttk.Button(self, text="...", width=3)
            # self.configButton2.place(x=150, y=98)

            self.step3Value = tk.BooleanVar(self)
            self.step3 = ttk.Checkbutton(self, text="Ensamblaje y prediccion de genes", variable=self.step3Value)
            self.step3.place(x=15, y=130)

            # botón para añadir la ruta por medio de una ventana
            # self.configButton3 = ttk.Button(self, text="...", width=3)
            # self.configButton3.place(x=150, y=138)

        elif self.proyectType == "Análisis de proteinas":
            # check buttons
            self.step4Value = tk.BooleanVar(self)
            self.step4 = ttk.Checkbutton(self, text="Proteina 1", variable=self.step4Value)
            self.step4.place(x=15, y=70)

            # botón para añadir la ruta por medio de una ventana
            # self.configButton1 = ttk.Button(self, text="...", width=3)
            # self.configButton1.place(x=150, y=68)

            self.step5Value = tk.BooleanVar(self)
            self.step5 = ttk.Checkbutton(self, text="Proteina 2", variable=self.step5Value)
            self.step5.place(x=15, y=100)

            # botón para añadir la ruta por medio de una ventana
            # self.configButton2 = ttk.Button(self, text="...", width=3)
            # self.configButton2.place(x=150, y=98)

            self.step6Value = tk.BooleanVar(self)
            self.step6 = ttk.Checkbutton(self, text="Proteina 3", variable=self.step6Value)
            self.step6.place(x=15, y=130)

            # botón para añadir la ruta por medio de una ventana
            # self.configButton3 = ttk.Button(self, text="...", width=3)
            # self.configButton3.place(x=150, y=138)

        # boton de inicio de proyecto
        self.confirmButton = ttk.Button(self, text="Iniciar Proyecto", command=self.onStartProyect)
        self.confirmButton.place(x=10, y=160)

        # boton de regreso
        self.confirmButton = ttk.Button(self, text="Regresar", command=self.onClosing)
        self.confirmButton.place(x=200, y=160)

    def onStartProyect(self):

        steps = []
        dirRoute = ""

        if self.proyectType == "Análisis de ADN":

            if self.step1Value.get() and not self.step2Value.get() and not self.step3Value.get():

                dirRoute = filedialog.asksaveasfilename()
                os.mkdir(dirRoute)

                steps = ["Ensamblaje y alineamiento"]
                os.mkdir(dirRoute+"/Ensamblaje")
                os.mkdir(dirRoute+"/Ensamblaje/Reads")
                steps += [dirRoute+"/Ensamblaje"]

                os.mkdir(dirRoute+"/Homologia")
                os.mkdir(dirRoute+"/Homologia/GR")
                steps += [dirRoute+"/Homologia"]

                if self.fileRoute != "" and self.file2route != "":
                    fileName = self.fileRoute.split("/")[-1]
                    os.rename(self.fileRoute, dirRoute+"/Ensamblaje/Reads/" + fileName)
                    refGenomaFileName = self.file2route.split("/")[-1]
                    os.rename(self.file2route, dirRoute+"/Homologia/GR/"+refGenomaFileName)
                    proyect = Proyect(self.proyectType, steps, dirRoute + "/archivo.bin", dirRoute+"/Ensamblaje/Reads/"
                                      + fileName, dirRoute+"/Homologia/GR/"+refGenomaFileName, dirRoute)
                    self.build_proyect_window.destroy()
                    new_window = tk.Tk()
                    panelwindow = PanelWindow(new_window, proyect, self.main_window)
                    panelwindow.mainloop()

                else:
                    self.top = tk.Toplevel(self.build_proyect_window)
                    self.top.title("Alerta")
                    tk.Label(self.top,
                             text="No ha seleccionado un archivo para analizar, por favor seleccione un archivo").grid(
                        row=0,
                        column=0,
                        columnspan=2)
                    self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
                    self.button2.grid(row=1, column=0, padx=5, pady=5)

            elif self.step2Value.get() and not self.step1Value.get() and not self.step3Value.get():

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

                if self.fileRoute != "":
                    sequenceRoute = self.fileRoute
                    fileName = self.fileRoute.split("/")[-1]
                    os.rename(self.fileRoute, dirRoute + "/Ensamblaje/Reads/" + fileName)
                    refGenomaFileName = self.file2route.split("/")[-1]
                    os.rename(self.file2route, dirRoute + "/Homologia/GR/" + refGenomaFileName)
                    proyect = Proyect(self.proyectType, steps, dirRoute + "/archivo.bin", sequenceRoute,
                                      dirRoute+"/Homologia/GR/" + refGenomaFileName, dirRoute)

                    self.build_proyect_window.destroy()
                    new_window = tk.Tk()
                    panelwindow = PanelWindow(new_window, proyect, self.main_window)
                    panelwindow.mainloop()

                else:
                    self.top = tk.Toplevel(self.build_proyect_window)
                    self.top.title("Alerta")
                    tk.Label(self.top,
                             text="No ha seleccionado un archivo para analizar, por favor seleccione un archivo").grid(
                        row=0,
                        column=0,
                        columnspan=2)
                    self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
                    self.button2.grid(row=1, column=0, padx=5, pady=5)

            elif self.step3Value.get() and not self.step1Value.get() and not self.step2Value.get():

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

                if self.fileRoute != "":
                    sequenceRoute = self.fileRoute
                    fileName = self.fileRoute.split("/")[-1]
                    os.rename(self.fileRoute, dirRoute + "/Ensamblaje/Reads/" + fileName)
                    refGenomaFileName = self.file2route.split("/")[-1]
                    os.rename(self.file2route, dirRoute + "/Homologia/GR/" + refGenomaFileName)
                    proyect = Proyect(self.proyectType, steps, dirRoute + "/archivo.bin", sequenceRoute,
                                      dirRoute+"/Homologia/GR/" + refGenomaFileName, dirRoute)

                    self.build_proyect_window.destroy()
                    new_window = tk.Tk()
                    panelwindow = PanelWindow(new_window, proyect, self.main_window)
                    panelwindow.mainloop()

                else:
                    self.top = tk.Toplevel(self.build_proyect_window)
                    self.top.title("Alerta")
                    tk.Label(self.top,
                             text="No ha seleccionado un archivo para analizar, por favor seleccione un archivo").grid(
                        row=0,
                        column=0,
                        columnspan=2)
                    self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
                    self.button2.grid(row=1, column=0, padx=5, pady=5)

        if self.proyectType == "Análisis de proteinas":

            if self.step4Value.get() and not self.step5Value.get() and not self.step6Value.get():

                dirRoute = filedialog.asksaveasfilename()
                os.mkdir(dirRoute)

                steps = ["Proteina1"]
                os.mkdir(dirRoute+"/Ensamblaje")
                steps += [dirRoute+"/Ensamblaje"]

                os.mkdir(dirRoute+"/Prediccion")
                steps += [dirRoute+"/Prediccion"]

                os.mkdir(dirRoute+"/Homologia")
                steps += [dirRoute+"/Homologia"]

                if self.fileRoute != "":
                    sequenceRoute = self.fileRoute
                    proyect = Proyect(self.proyectType, steps, dirRoute + "/archivo.bin", sequenceRoute, dirRoute)

                    self.build_proyect_window.destroy()
                    new_window = tk.Tk()
                    panelwindow = PanelWindow(new_window, proyect, self.main_window)
                    panelwindow.mainloop()

                else:
                    self.top = tk.Toplevel(self.build_proyect_window)
                    self.top.title("Alerta")
                    tk.Label(self.top,
                             text="No ha seleccionado un archivo para analizar, por favor seleccione un archivo").grid(
                        row=0,
                        column=0,
                        columnspan=2)
                    self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
                    self.button2.grid(row=1, column=0, padx=5, pady=5)

            elif self.step5Value.get() and not self.step4Value.get() and not self.step6Value.get():

                dirRoute = filedialog.asksaveasfilename()
                os.mkdir(dirRoute)

                steps = ["Proteina2"]
                os.mkdir(dirRoute + "/Ensamblaje")
                steps += [dirRoute + "/Ensamblaje"]

                os.mkdir(dirRoute + "/Homologia")
                steps += [dirRoute + "/Homologia"]

                if self.fileRoute != "":
                    sequenceRoute = self.fileRoute
                    proyect = Proyect(self.proyectType, steps, dirRoute + "/archivo.bin", sequenceRoute, dirRoute)

                    self.build_proyect_window.destroy()
                    new_window = tk.Tk()
                    panelwindow = PanelWindow(new_window, proyect, self.main_window)
                    panelwindow.mainloop()

                else:
                    self.top = tk.Toplevel(self.build_proyect_window)
                    self.top.title("Alerta")
                    tk.Label(self.top,
                             text="No ha seleccionado un archivo para analizar, por favor seleccione un archivo").grid(
                        row=0,
                        column=0,
                        columnspan=2)
                    self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
                    self.button2.grid(row=1, column=0, padx=5, pady=5)

            elif self.step6Value.get() and not self.step4Value.get() and not self.step5Value.get():

                dirRoute = filedialog.asksaveasfilename()
                os.mkdir(dirRoute)

                steps = ["Proteina3"]
                os.mkdir(dirRoute + "/Ensamblaje")
                steps += [dirRoute + "/Ensamblaje"]

                os.mkdir(dirRoute + "/Prediccion")
                steps += [dirRoute + "/Prediccion"]

                os.mkdir(dirRoute + "/Homologia")
                steps += [dirRoute + "/Homologia"]

                if self.fileRoute != "":
                    sequenceRoute = self.fileRoute
                    proyect = Proyect(self.proyectType, steps, dirRoute + "/archivo.bin", sequenceRoute, dirRoute)

                    self.build_proyect_window.destroy()
                    new_window = tk.Tk()
                    panelwindow = PanelWindow(new_window, proyect, self.main_window)
                    panelwindow.mainloop()

                else:
                    self.top = tk.Toplevel(self.build_proyect_window)
                    self.top.title("Alerta")
                    tk.Label(self.top,
                             text="No ha seleccionado un archivo para analizar, por favor seleccione un archivo").grid(
                        row=0,
                        column=0,
                        columnspan=2)
                    self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
                    self.button2.grid(row=1, column=0, padx=5, pady=5)

        else:
            self.top = tk.Toplevel(self.build_proyect_window)
            self.top.title("Alerta")
            tk.Label(self.top,
                     text="Porfavor seleccione un solo grupo de pasos").grid(row=0, column=0, columnspan=2)
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
