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
from buildProyect import BuildProyectWindow
from panel import PanelWindow
import webbrowser
import sys
import pickle


class MainWindow(ttk.Frame):

    def __init__(self, main_window):

        # Intancias y configuraciones de la ventana
        self.main_window = main_window
        super().__init__(self.main_window)
        self.main_window.title("BIOIN")
        self.main_window.geometry("500x400")
        self.main_window.resizable(0, 0)
        self.main_window.protocol("WM_DELETE_WINDOW", self.onClosing)
        self.place(relwidth=1, relheight=1)

        # menu de la pantalla de inicio
        self.menubar = tk.Menu(self.main_window)
        self.fileMenu = tk.Menu(self.menubar, tearoff=0)
        self.fileMenu.add_command(label="Nuevo Proyecto", command=self.newProyectWindow)
        self.fileMenu.add_command(label="Abrir...", command=self.openProyect)

        self.fileMenu.add_separator()

        self.fileMenu.add_command(label="Salir", command=self.onClosing)
        self.menubar.add_cascade(label="Proyecto", menu=self.fileMenu)

        self.helpMenu = tk.Menu(self.menubar, tearoff=0)
        self.helpMenu.add_command(label="Más Ayuda", command=self.openWebHelp)
        self.helpMenu.add_command(label="Acerca de Bioin...", command=self.openWeb)
        self.menubar.add_cascade(label="Ayuda", menu=self.helpMenu)

        main_window.config(menu=self.menubar)

        # Cargar la imagen del logo
        self.logo = tk.PhotoImage(file="res/images/logo.png")
        self.logo = self.logo.subsample(2)

        # Label del Logo
        self.logoLabel = ttk.Label(self, image=self.logo)
        self.logoLabel.place(x=100, y=-100)

        # Botones de la pantalla de inicio
        self.startNewProyectButton = ttk.Button(self, text="Nuevo Proyecto", width=15, command=self.newProyectWindow)
        self.startNewProyectButton.place(x=200, y=210)

        self.loadProyectButton = ttk.Button(self, text="Cargar Proyecto", width=15, command=self.openProyect)
        self.loadProyectButton.place(x=200, y=240)

        self.settingsButton = ttk.Button(self, text="Ajustes", width=15)
        self.settingsButton.place(x=200, y=270)

        self.infoButton = ttk.Button(self, text="Más Información", width=15, command=self.openWeb)
        self.infoButton.place(x=200, y=300)

    # Funciones de la interfaz
    def newProyectWindow(self):
        self.main_window.withdraw()
        new_window = tk.Tk()
        buildProyectWindow = BuildProyectWindow(new_window, self.main_window)
        buildProyectWindow.mainloop()

    def onClosing(self):
        self.top = tk.Toplevel(self.main_window)
        self.top.title("Salir")

        tk.Label(self.top, text="¿Está seguro?").grid(row=0, column=0, columnspan=2)

        self.button1 = tk.Button(self.top, text="Si, deseo salir.", command=self.salir)
        self.button2 = tk.Button(self.top, text="Cancelar", command=self.cancelar)
        self.button1.grid(row=1, column=0, padx=5, pady=5)
        self.button2.grid(row=1, column=1, padx=5, pady=5)

    def openProyect(self):
        dirRoute = filedialog.askopenfilename()
        if dirRoute != () and dirRoute != '':
            self.main_window.withdraw()
            with open(dirRoute, "br") as archivo:
                proyect = pickle.load(archivo)
            new_window = tk.Tk()
            panelwindow = PanelWindow(new_window, proyect, self.main_window)
            panelwindow.mainloop()

    def salir(self):
        self.top.destroy()
        self.main_window.destroy()
        sys.exit()

    def cancelar(self):
        self.top.destroy()

    def openWeb(self):
        url = "http://bioinformatica.univalle.edu.co/"
        webbrowser.open(url)

    def openWebHelp(self):
        url = "http://bioinformatica.univalle.edu.co/"
        webbrowser.open(url)


# Main donde todo inicia

if __name__ == "__main__":
    root = tk.Tk()
    app = MainWindow(root)
    root.mainloop()
