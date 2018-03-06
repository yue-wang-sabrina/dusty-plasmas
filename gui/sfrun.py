# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 15:41:01 2017

@author: sabrina.yue.wang
"""
import sys
import os
from PyQt5 import QtCore, QtGui, QtWidgets
from sf import Ui_Dialog
from PyQt5.QtCore import QThread, pyqtSignal, pyqtSlot, QObject
from PyQt5.QtWidgets import QApplication, QMainWindow, QInputDialog, QFileDialog, QStatusBar
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy
import pickle
import math
import importlib
import os.path

from analysis_dust import BEffectsAnalysis
from utils.utils import generate_particle_equilibrium_positions, prepare_modified_b_field
from plots import dustplots

from IPython import get_ipython


ipython = get_ipython()
ipython.magic('load_ext autoreload')
ipython.magic('autoreload 2')


class SpectrumCanvas(FigureCanvas):
    '''
    Provides ability for displaying figures and animations on gui
    '''
    def __init__(self, parent=None, width=7.31, height=4.21, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
    
    


class SFGui(Ui_Dialog): #Setting up/Connecting the gui buttons and connecting them to their associated functions
    def __init__(self, dialog):
        Ui_Dialog.__init__(self)
        self.setupUi(dialog)
        self.graph = SpectrumCanvas(self.graphicsView)
        self.graph.setObjectName("graph")
        self.error=QtWidgets.QGraphicsScene()
        self.pushButton.clicked.connect(self.simtime)
        self.pushButton_2.clicked.connect(self.particlenumber)
        self.pushButton_2.clicked.connect(self.run)
        self.time = 0;
        self.particlenumber = 0;
 
    def addInputTextToListbox(self): #Add user input
        txt = self.myTextInput.text()

    def simtime(self):
        text2,ok2 = QInputDialog.getText(dialog, 'User Input','Enter simulation time in seconds')
        self.time=float(text2.split()[0])

    def particlenumber(self):
        text, ok = QInputDialog.getText(dialog, 'User Input','Enter particle number')
        self.particlenumber = float(text.split()[0])
    def run(self):
        beffect1 = BEffectsAnalysis()
        beffect1.create_particles(
            numparticles=10,
            initpositions=generate_particle_equilibrium_positions()
        )
        beffect1.create_pairs()
        beffect1.interact_and_iterate(
            iterationsB=100,
            init_iterations=100,
            method='Gibs',
            modified_b_field=prepare_modified_b_field()
        )
beffect1.sort_positions_of_particles()

    


        
 
if __name__ == '__main__':
	app = QtWidgets.QApplication(sys.argv)
	dialog = QtWidgets.QDialog()
 
	prog = SFGui(dialog)
 
	dialog.show()
	sys.exit(app.exec_())
    

    