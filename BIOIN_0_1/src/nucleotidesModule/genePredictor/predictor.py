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


class Predictor:

    def __init__(self, project):
        self.project = project

    def ejecutCommand(self):
        # guardamos en una variable el comando a ejecutar
        comando = "./nucleotidesModule/genePredictor/glimmer3.02/scripts/g3-iterated.csh " + self.project.dirRoute +\
                  "/Ensamblaje/output_assembly/output_d_results/output_out.padded.fasta " + self.project.dirRoute + "/Prediccion/output"

        # convertimos el string en una lista para poder pasar de manera adecuada los comandos desde python
        args = sl.split(comando)

        # ejecutamos la función call de subprocess que permite ejecutar comandos desde la temrinal
        sp.call(args)

        # traducimos los resultados a formato fasta
        comando2 = "./nucleotidesModule/genePredictor/glimerTranslate.py " \
                   + self.project.dirRoute + "/Prediccion/output.predict " \
                   + self.project.dirRoute + "/Ensamblaje/output_assembly/output_d_results/output_out.padded.fasta " \
                   + self.project.dirRoute + "/Prediccion/output.fasta"
        args2 = sl.split(comando2)
        sp.call(args2)

    def fileExist(self):
        path = self.project.dirRoute + "/Prediccion/output.fasta"
        return os.path.exists(path)
