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
        comando2 = "./nucleotidesModule/genePredictor/glimmer3.02/bin/build-icm " + self.project.genomeRefRoute
        print(comando2)
        args2 = sl.split(comando2)
        sp.call(args2)

        # guardamos en una variable el comando a ejecutar
        comando = "./nucleotidesModule/genePredictor/glimmer3.02/bin/glimmer3 " + self.project.dirRoute +\
                  "/Ensamblaje/output_assembly/output_d_results/output_out.padded.fasta " + \
                  self.project.dirRoute + "/Prediccion/output.icm " + self.project.dirRoute + "/Prediccion/output"

        # convertimos el string en una lista para poder pasar de manera adecuada los comandos desde python
        args = sl.split(comando)

        # ejecutamos la función call de subprocess que permite ejecutar comandos desde la temrinal
        sp.call(args)

    def fileExist(self):
        path = self.project.dirRoute + "/Prediccion/output.detail"
        return os.path.exists(path)
