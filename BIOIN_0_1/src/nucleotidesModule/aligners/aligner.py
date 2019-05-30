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
        comando = "./nucleotidesModule/aligners/blastn -query " + self.proyect.dirRoute + \
                  "/Ensamblaje/output_assembly/output_d_results/output_out.padded.fasta -subject " + \
                  self.proyect.genomeRefRoute + " -out " + self.proyect.dirRoute + "/Homologia/output.xml -outfmt 5"

        # convertimos el string en una lista para poder pasar de manera adecuada los comandos desde python
        args = sl.split(comando)

        # ejecutamos la función call de subprocess que permite ejecutar comandos desde la temrinal
        sp.call(args)

    def fileExist(self):
        path = self.proyect.dirRoute + "/Homologia/output.xml"
        return os.path.exists(path)

    def ejecutCommand2(self):
        # guardamos en una variable el comando a ejecutar
        comando = "./nucleotidesModule/aligners/blastp -query " + self.proyect.dirRoute + \
                  "/ORFoma/output.fasta -db nr -remote -out " + \
                  self.proyect.dirRoute + "/Homologia/output.xml -outfmt 5"

        # convertimos el string en una lista para poder pasar de manera adecuada los comandos desde python
        args = sl.split(comando)

        # ejecutamos la función call de subprocess que permite ejecutar comandos desde la temrinal
        sp.call(args)

    def ejecutCommand3(self):
        # guardamos en una variable el comando a ejecutar
        comando = "./nucleotidesModule/aligners/blastx -query " + self.proyect.dirRoute + \
                  "/Prediccion/output.fasta -db nr -remote -out " + \
                  self.proyect.dirRoute + "/Homologia/output.xml -outfmt 5"

        # convertimos el string en una lista para poder pasar de manera adecuada los comandos desde python
        args = sl.split(comando)

        # ejecutamos la función call de subprocess que permite ejecutar comandos desde la temrinal
        sp.call(args)
