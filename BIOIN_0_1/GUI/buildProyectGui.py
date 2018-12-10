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
from panelGui import PanelWindow
from proyect import Proyect
import webbrowser
import os
import pickle


class BuildProyectWindow(ttk.Frame):

    def __init__(self, build_proyect_window, main_window):

        # Intancias y configuraciones de la ventana
        self.build_proyect_window = build_proyect_window
        self.main_window = main_window
        super().__init__(self.build_proyect_window)
        self.build_proyect_window.title("BIOIN - NUEVO PROYECTO")
        self.build_proyect_window.geometry("500x400")
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

        # check buttons
        self.step1Value = tk.BooleanVar(self)
        self.step1 = ttk.Checkbutton(self, text="Paso 1", variable=self.step1Value)
        self.step1.place(x=15, y=70)

        self.step2Value = tk.BooleanVar(self)
        self.step2 = ttk.Checkbutton(self, text="Paso 2", variable=self.step2Value)
        self.step2.place(x=15, y=90)

        self.step3Value = tk.BooleanVar(self)
        self.step3 = ttk.Checkbutton(self, text="Paso 3", variable=self.step3Value)
        self.step3.place(x=15, y=110)

        self.step4Value = tk.BooleanVar(self)
        self.step4 = ttk.Checkbutton(self, text="Paso 4", variable=self.step4Value)
        self.step4.place(x=15, y=130)

        self.step5Value = tk.BooleanVar(self)
        self.step5 = ttk.Checkbutton(self, text="Paso 5", variable=self.step5Value)
        self.step5.place(x=15, y=150)

        # boton de inicio de proyecto
        self.confirmButton = ttk.Button(self, text="Iniciar Proyecto", command=self.onStartProyect)
        self.confirmButton.place(x=10, y=190)

        # boton de regreso
        self.confirmButton = ttk.Button(self, text="Regresar", command=self.onClosing)
        self.confirmButton.place(x=200, y=190)

    def onStartProyect(self):
        stepsList = []
        stepsList += [self.step1Value.get()]
        stepsList += [self.step2Value.get()]
        stepsList += [self.step3Value.get()]
        stepsList += [self.step4Value.get()]
        stepsList += [self.step5Value.get()]

        dirRoute = filedialog.asksaveasfilename()
        os.mkdir(dirRoute)
        proyect = Proyect(stepsList, dirRoute+"/archivo.bin")

        self.build_proyect_window.destroy()
        new_window = tk.Tk()
        panelwindow = PanelWindow(new_window, proyect, self.main_window)
        panelwindow.mainloop()

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
        if dirRoute != ():
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