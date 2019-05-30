#!/usr/bin/python
# Proyecto: BIOIN
# Trabajo de Grado
# Alejandro Valencia
# Juan Jose Varela
# Universidad del Valle


class Proyect:

    def __init__(self, type, steps, route, sequenceRoute, genomeRefRoute, dirRoute):
        self.type = type
        self.route = route
        self.steps = steps
        self.sequenceRoute = sequenceRoute
        self.genomeRefRoute = genomeRefRoute
        self.dirRoute = dirRoute
