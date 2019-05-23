#!/usr/bin/python
# Proyecto: BIOIN
# Trabajo de Grado
# Alejandro Valencia R.
# Juan Jose Varela V.
# Universidad del Valle

# librerias a importar
import subprocess as sp
import shlex as sl
import os


class Aligner:

    def __init__(self, proyect):
        self.proyect = proyect

    def ejecutCommand(self):
        # guardamos en una variable el comando a ejecutar
        comando = "./nucleotidesModule/aligners/clustalw2 -ALIGN -INFILE="
        comando += self.inputfile + " -OUTFILE=" + self.outputfile + "/output"

        # convertimos el string en una lista para poder pasar de manera adecuada los comandos desde python
        args = sl.split(comando)

        # ejecutamos la función call de subprocess que permite ejecutar comandos desde la temrinal
        sp.call(args)

    def fileExist(self):
        path = self.outputfile + "/output"
        return os.path.exists(path)

