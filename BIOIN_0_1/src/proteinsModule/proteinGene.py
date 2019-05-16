#!/usr/bin/python
# Proyecto: BIOIN
# Trabajo de Grado
# Alejandro Valencia R.
# Juan Jose Varela V.
# Universidad del Valle

# librerias a importar
import subprocess as sp
import shlex as sl


class ProteinGene:

    def __init__(self, name, script, config, route):
        self.name = name
        self.script = script
        self.config = config
        self.route = route