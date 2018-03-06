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
import importlib
import os.path
import matplotlib.pyplot as plt
import matplotlib.animation
import mpl_toolkits.mplot3d.axes3d as p3

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

from pylab import rcParams

from msci.analysis.analysis_dust import BEffectsAnalysis
from msci.utils.utils import generate_particle_equilibrium_positions, prepare_modified_b_field
import msci.analysis.constants as const

from IPython import get_ipython

ipython = get_ipython()
ipython.magic('load_ext autoreload')
ipython.magic('autoreload 2')


class SpectrumCanvas(FigureCanvas):
    '''
    Provides ability for displaying figures and animations on gui
    '''

    def __init__(self, parent=None, width=7.21, height=5.91, dpi=100):
        self.fig = Figure(figsize=(width*3, height*3), dpi=dpi)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def display(self, dustanalysis):
        self.fig.clf()
        test = self.fig.add_subplot(111, projection='3d')

        def _update_graph(n_iter):
            data = dustanalysis.positions_df[dustanalysis.positions_df['time'] == n_iter]
            point.set_data(data.x, data.y)
            point.set_3d_properties(data.z)
            return point,

        # test.view_init(elev=90., azim=90)
        if min(dustanalysis.position_array[:, 0]) == max(dustanalysis.position_array[:, 0]):
            test.set_xlim([-10 * dustanalysis.const.lambdaD, 10 * dustanalysis.const.lambdaD])
        else:
            test.set_xlim([min(dustanalysis.position_array[:, 0]), max(dustanalysis.position_array[:, 0])])
        if min(dustanalysis.position_array[:, 1]) == max(dustanalysis.position_array[:, 1]):
            test.set_ylim([-10 * dustanalysis.const.lambdaD, 10 * dustanalysis.const.lambdaD])
        else:
            test.set_ylim([min(dustanalysis.position_array[:, 1]), max(dustanalysis.position_array[:, 1])])
        if min(dustanalysis.position_array[:, 2]) == max(dustanalysis.position_array[:, 2]):
            test.set_zlim([-10 * dustanalysis.const.lambdaD, 10 * dustanalysis.const.lambdaD])
        else:
            test.set_zlim(
                [min(dustanalysis.position_array[:, 2]) - 0.5 * min(dustanalysis.position_array[:, 2]),
                 max(dustanalysis.position_array[:, 2])])

        test.set_xlim([-dustanalysis.modified_b_field['rmax'] * 1.5, dustanalysis.modified_b_field['rmax'] * 1.5])
        test.set_ylim([-dustanalysis.modified_b_field['rmax'] * 1.5, dustanalysis.modified_b_field['rmax'] * 1.5])
        test.set_xlabel("x (m)")
        test.set_ylabel("y (m)")

        data = dustanalysis.positions_df[dustanalysis.positions_df['time'] == 0]
        point, = test.plot(data.x, data.y, data.z, linestyle="", marker=".")

        self.ani = matplotlib.animation.FuncAnimation(
            self.fig,
            _update_graph,
            frames=dustanalysis.iterationsB + dustanalysis.init_iterations,
            interval=1,
            blit=False
        )
        self.draw()


class SFGui(Ui_Dialog):  # Setting up/Connecting the gui buttons and connecting them to their associated functions
    def __init__(self, dialog):
        Ui_Dialog.__init__(self)
        self.setupUi(dialog)
        self.graph = SpectrumCanvas(self.graphicsView)
        self.graph.setObjectName("graph")
        self.pushButton.clicked.connect(self.simtime)
        self.pushButton_2.clicked.connect(self.particlenumber)
        self.pushButton_3.clicked.connect(self.run)
        self.time = 0;
        self.particlenumber = 0;

    # def addInputTextToListbox(self):  # Add user input
    #     txt = self.myTextInput.text()

    def simtime(self):
        text2, ok2 = QInputDialog.getText(dialog, 'User Input', 'Enter simulation time in seconds')
        self.time = float(text2.split()[0])

    def particlenumber(self):
        text, ok = QInputDialog.getText(dialog, 'User Input', 'Enter particle number')
        self.particlenumber = int(text.split()[0])

    def run(self):
        self.beffect1 = BEffectsAnalysis()
        self.beffect1.create_particles(
            numparticles=self.particlenumber,
            initpositions=generate_particle_equilibrium_positions()
        )
        self.beffect1.create_pairs()
        self.beffect1.interact_and_iterate(
            iterationsB=100,
            init_iterations=int(self.time / const.dt),
            method='Gibs',
            modified_b_field=prepare_modified_b_field()
        )
        self.beffect1.sort_positions_of_particles()
        self.graph.display(self.beffect1)

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    dialog = QtWidgets.QDialog()

    prog = SFGui(dialog)

    dialog.show()
    sys.exit(app.exec_())
