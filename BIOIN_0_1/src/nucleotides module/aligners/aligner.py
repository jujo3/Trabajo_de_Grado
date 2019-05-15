#!/usr/bin/python
# Proyecto: BIOIN
# Trabajo de Grado
# Alejandro Valencia R.
# Juan Jose Varela V.
# Universidad del Valle

# librerias a importar
import subprocess as sp
import shlex as sl


class Aligner:

    def __init__(self, name, script, config, route):
        self.name = name
        self.script = script
        self.config = config
        self.route = route

    def ejecutCommand(self):
        # guardamos en una variable el comando a ejecutar
        comando = "sudo"

        # convertimos el string en una lista para poder pasar de manera adecuada los comandos desde python
        args = sl.split(comando)

        # ejecutamos la función call de subprocess que permite ejecutar comandos desde la temrinal
        sp.call(args)
