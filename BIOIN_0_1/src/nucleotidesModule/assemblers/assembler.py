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


class Assembler:

    def __init__(self, inputfile, outputfile):
        self.inputfile = inputfile
        self.outputfile = outputfile

    def ejecutCommand(self):

        # creamos el archivo manifest para ejecutar mira:
        f = open(self.outputfile + "/manifest.txt", "w+")
        manifest = '''project = output
job = genome,denovo,draft
parameters = -NW:cmrnl=no
readgroup = project 
data = ruta
technology = solexa'''

        manifest = manifest.replace("ruta", self.inputfile)
        f.write(manifest)
        f.close()

        # guardamos en una variable el comando a ejecutar
        comando = "./nucleotidesModule/assemblers/mira " + self.outputfile + "/manifest.txt"

        # convertimos el string en una lista para poder pasar de manera adecuada los comandos desde python
        args = sl.split(comando)

        # ejecutamos la función call de subprocess que permite ejecutar comandos desde la temrinal
        sp.call(args)

        # movemos la carpeta resultado al directorio correspondiente:
        os.rename("output_assembly", self.outputfile + "/output_assembly")

    def fileExist(self):
        path = self.outputfile + "/output_assembly/output_d_results/output_out.padded.fasta"
        return os.path.exists(path)
