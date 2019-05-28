#!/usr/bin/python
# Proyecto: BIOIN
# Trabajo de Grado
# Alejandro Valencia R.
# Juan Jose Varela V.
# Universidad del Valle

# librerias a importar
import subprocess as sp
import shlex as sl
import webbrowser
import threading


class Browser:

    def __init__(self, project):
        self.project = project

    def worker(self):
        comando2 = "python3 -m http.server 8081"
        args2 = sl.split(comando2)
        sp.call(args2)

    def ejecutCommand(self):
        # comando = "cd /nucleotidesModule/genomeBrowsers/kablammo"
        # args = sl.split(comando)
        # sp.call(args)
        # montaje de servidor en segundo plano para no detener la interfaz
        thread = threading.Thread(target=self.worker)
        thread.start()
        url = "http://localhost:8081/nucleotidesModule/genomeBrowsers/kablammo/"
        webbrowser.open(url)
        # guardamos en una variable el comando a ejecutar
        # comando = "./nucleotidesModule/genomeBrowsers/JBrowse-1.16.4-desktop-linux-x64/JBrowse-1.16.4-desktop"

        # convertimos el string en una lista para poder pasar de manera adecuada los comandos desde python
        # args = sl.split(comando)

        # ejecutamos la función call de subprocess que permite ejecutar comandos desde la temrinal
        # sp.call(args)
