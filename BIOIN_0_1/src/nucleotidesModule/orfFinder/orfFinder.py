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


class OrfFinder:

    def __init__(self, proyect):
        self.proyect = proyect

    def ejecutCommand(self):
        # guardamos en una variable el comando a ejecutar
        comando = "./nucleotidesModule/orfFinder/ORFfinder -in " + self.proyect.dirRoute +\
                  "/Ensamblaje/output_assembly/output_d_results/output_out.padded.fasta -out " + \
                  self.proyect.dirRoute + "/ORFoma/output.fasta"

        # convertimos el string en una lista para poder pasar de manera adecuada los comandos desde python
        args = sl.split(comando)

        # ejecutamos la función call de subprocess que permite ejecutar comandos desde la temrinal
        sp.call(args)

    def fileExist(self):
        path = self.proyect.dirRoute + "/ORFoma/output.fasta"
        return os.path.exists(path)
