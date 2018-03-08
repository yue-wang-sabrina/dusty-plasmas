# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 15:41:01 2017

@author: sabrina.yue.wang
"""
import sys
import os
from PyQt5 import QtCore, QtGui, QtWidgets
from sftest import Ui_Dialog
from PyQt5.QtCore import QThread, pyqtSignal, pyqtSlot, QObject
from PyQt5.QtWidgets import QApplication, QMainWindow, QInputDialog, QFileDialog, QStatusBar
import matplotlib

matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy
import pickle
import math
import simulationsensor
from simulationsensor import dt
import importlib
import os.path


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

    def display(self):  # Creates interactive animation (including EMTS and sensor fusion reconstruction)
        self.fig.clf()
        test = self.fig.add_subplot(111, projection='3d')
        if os.path.exists("simulatedsensor.obj"):
            import fusion as fuse
            importlib.reload(fuse)
            import EMTS as em
            importlib.reload(em)
            self.anilen = len(fuse.orientationvectorz)
            test.set_title("Reconstruction from simulated sensor data")
            self.maximum = max(
                [max((numpy.array(fuse.globalposition)[:, 0])), max((numpy.array(fuse.globalposition)[:, 1])),
                 max((numpy.array(fuse.globalposition)[:, 2]))])
            self.minimum = min(
                [min(numpy.array(fuse.globalposition)[:, 0]), min(numpy.array(fuse.globalposition)[:, 1]),
                 min(numpy.array(fuse.globalposition)[:, 2])])
            # self.maximum2=max([max((numpy.array(em.EMTSx))),max((numpy.array(em.EMTSy))),max((numpy.array(em.EMTSz)))])
            # self.minimum2=min([min((numpy.array(em.EMTSx))),min((numpy.array(em.EMTSy))),min((numpy.array(em.EMTSz)))])
            # self.maximum=max(self.maximum,self.maximum2)
            # self.minimum=min(self.minimum,self.minimum2)
            if abs(self.maximum) < 10 ** (-3) and abs(self.minimum) < 10 ** (-3):
                self.maximum = 1
                self.minimum = -1
            test.set_xlim([self.minimum, self.maximum])
            test.set_ylim([self.minimum, self.maximum])
            test.set_zlim([self.minimum, self.maximum])
            test.quiver3D(fuse.globalposition[0][0], fuse.globalposition[0][1], fuse.globalposition[0][2],
                          fuse.orientationvectorx[0][0], fuse.orientationvectorx[0][1], fuse.orientationvectorx[0][2],
                          length=abs(self.maximum * 0.5), arrow_length_ratio=0.5, color='r')
            test.quiver3D(fuse.globalposition[0][0], fuse.globalposition[0][1], fuse.globalposition[0][2],
                          fuse.orientationvectory[0][0], fuse.orientationvectory[0][1], fuse.orientationvectory[0][2],
                          length=abs(0.5 * self.maximum), arrow_length_ratio=0.5, color='k')
            test.quiver3D(fuse.globalposition[0][0], fuse.globalposition[0][1], fuse.globalposition[0][2],
                          fuse.orientationvectorz[0][0], fuse.orientationvectorz[0][1], fuse.orientationvectorz[0][2],
                          length=abs(0.5 * self.maximum), arrow_length_ratio=0.5, color='b')
            # test.quiver3D(em.EMTSx[0],em.EMTSy[0],em.EMTSz[0],-em.EMTSorientationvectorx[0][0],-em.EMTSorientationvectorx[0][1],-em.EMTSorientationvectorx[0][2],length=abs(self.maximum*0.6),arrow_length_ratio=0.5,color='b')
            # test.quiver3D(em.EMTSx[0],em.EMTSy[0],em.EMTSz[0],-em.EMTSorientationvectory[0][0],-em.EMTSorientationvectory[0][1],-em.EMTSorientationvectory[0][2],length=abs(self.maximum*0.6),arrow_length_ratio=0.5,color='k')
            # test.quiver3D(em.EMTSx[0],em.EMTSy[0],em.EMTSz[0],-em.EMTSorientationvectorz[0][0],-em.EMTSorientationvectorz[0][1],-em.EMTSorientationvectorz[0][2],length=abs(self.maximum*0.6),arrow_length_ratio=0.5,color='r')
            test.set_xlabel('x')
            test.set_ylabel('y')
            test.set_zlabel('z')

            def update_graph(p):
                test.clear()
                test.plot([item[0] for item in fuse.globalposition], [item[1] for item in fuse.globalposition],
                          [item[2] for item in fuse.globalposition], color='r', label='Using orientation updates')
                # test.plot(em.EMTSx,em.EMTSy,em.EMTSz,color='b',label='EMTS measurements')
                test.quiver3D(fuse.globalposition[p][0], fuse.globalposition[p][1], fuse.globalposition[p][2],
                              fuse.orientationvectorx[p][0], fuse.orientationvectorx[p][1],
                              fuse.orientationvectorx[p][2], length=abs((self.maximum - self.minimum) * 0.5),
                              arrow_length_ratio=0.5, color='r')
                test.quiver3D(fuse.globalposition[p][0], fuse.globalposition[p][1], fuse.globalposition[p][2],
                              fuse.orientationvectory[p][0], fuse.orientationvectory[p][1],
                              fuse.orientationvectory[p][2], length=abs((self.maximum - self.minimum) * 0.5),
                              arrow_length_ratio=0.5, color='k')
                test.quiver3D(fuse.globalposition[p][0], fuse.globalposition[p][1], fuse.globalposition[p][2],
                              fuse.orientationvectorz[p][0], fuse.orientationvectorz[p][1],
                              fuse.orientationvectorz[p][2], length=abs((self.maximum - self.minimum) * 0.5),
                              arrow_length_ratio=0.5, color='b')
                # test.quiver3D(em.EMTSx[p],em.EMTSy[p],em.EMTSz[p],-em.EMTSorientationvectorx[p][0],-em.EMTSorientationvectorx[p][1],-em.EMTSorientationvectorx[p][2],length=abs((self.maximum-self.minimum)*0.6),arrow_length_ratio=0.5,color='b')
                # test.quiver3D(em.EMTSx[p],em.EMTSy[p],em.EMTSz[p],-em.EMTSorientationvectory[p][0],-em.EMTSorientationvectory[p][1],-em.EMTSorientationvectory[p][2],length=abs((self.maximum-self.minimum)*0.6),arrow_length_ratio=0.5,color='k')
                # test.quiver3D(em.EMTSx[p],em.EMTSy[p],em.EMTSz[p],-em.EMTSorientationvectorz[p][0],-em.EMTSorientationvectorz[p][1],-em.EMTSorientationvectorz[p][2],length=abs((self.maximum-self.minimum)*0.6),arrow_length_ratio=0.5,color='b')

                test.set_xlim([self.minimum, self.maximum])
                test.set_ylim([self.minimum, self.maximum])
                test.set_zlim([self.minimum, self.maximum])
                test.set_xlabel('x')
                test.set_ylabel('y')
                test.set_zlabel('z')
                test.set_title("Reconstruction from simulated sensor data")
                test.legend()

            self.ani = matplotlib.animation.FuncAnimation(self.fig, update_graph, frames=self.anilen, interval=10,
                                                          blit=False)
            self.draw()

        else:
            self.error.addText("Error: No data exists to plot!")

    def clear(self):  # clear figure from gui
        self.fig.clf()


class SFGui(Ui_Dialog):  # Setting up/Connecting the gui buttons and connecting them to their associated functions
    def __init__(self, dialog):
        Ui_Dialog.__init__(self)
        self.setupUi(dialog)
        self.graph = SpectrumCanvas(self.graphicsView)
        self.graph.setObjectName("graph")
        self.error = QtWidgets.QGraphicsScene()
        self.graphicsView_2.setScene(self.error)
        self.pushButton.clicked.connect(self.run)
        self.pushButton_2.clicked.connect(self.linearacc)
        self.pushButton_3.clicked.connect(self.rotacc)
        self.pushButton_6.clicked.connect(self.initorient)
        self.pushButton_4.clicked.connect(self.deleteitem)
        self.pushButton_5.clicked.connect(self.clearall)
        self.comboBox.activated[str].connect(self.combobox)
        self.allacc = []  # All inputs are added to this list to be evaluated by simulationsensor.py
        self.dt = dt  # Time step
        self.count = 0  # TCount of how many individual motions are inputed to be evaluated
        os.system('AllClasses.py')
        os.system('simulationsensor.py')

    def addInputTextToListbox(self):  # Add user input
        txt = self.myTextInput.text()
        self.listWidget.addItem(txt)

    def run(self):  # Run simulation of all inputs shown in the listWidget
        if len(self.allacc) == 0:
            self.error.addText("Error: No data for analysis")
        else:
            importlib.reload(simulationsensor)
            for i in self.allacc:  # Add all the inputs into the simulationsensor.py individually depending on whether it is a rotation or acceleration
                if i[1] == 0:  # Linear acceleration
                    simulationsensor.data.updatelinaccfromGUI(i[0])
                    simulationsensor.appenddatamotion(simulationsensor.data, flatness=0)
                elif i[1] == 1:  # Angular acc
                    simulationsensor.data.updateangfromGUI(i[0])
                    simulationsensor.appenddatarot(simulationsensor.data)
            filehandler = open(b"simulatedsensor.obj", "wb")
            pickle.dump(simulationsensor.ax, filehandler)
            pickle.dump(simulationsensor.ay, filehandler)
            pickle.dump(simulationsensor.az, filehandler)
            pickle.dump(simulationsensor.phidot, filehandler)
            pickle.dump(simulationsensor.thetadot, filehandler)
            pickle.dump(simulationsensor.psidot, filehandler)
            pickle.dump(simulationsensor.compx, filehandler)
            pickle.dump(simulationsensor.compy, filehandler)
            pickle.dump(simulationsensor.compz, filehandler)
            pickle.dump(simulationsensor.xcompx, filehandler)
            pickle.dump(simulationsensor.xcompy, filehandler)
            pickle.dump(simulationsensor.xcompz, filehandler)
            pickle.dump(simulationsensor.xax, filehandler)
            pickle.dump(simulationsensor.xay, filehandler)
            pickle.dump(simulationsensor.xaz, filehandler)
            filehandler.close()

            simulationsensor.EMTSgetdata()
            filehandler2 = open(b"EMTSsimulation.obj", "wb")
            pickle.dump(simulationsensor.EMTSx, filehandler2)
            pickle.dump(simulationsensor.EMTSy, filehandler2)
            pickle.dump(simulationsensor.EMTSz, filehandler2)
            pickle.dump(simulationsensor.EMTSroll, filehandler2)
            pickle.dump(simulationsensor.EMTSpitch, filehandler2)
            pickle.dump(simulationsensor.EMTSyaw, filehandler2)
            pickle.dump(simulationsensor.EMTSorientationvectorx, filehandler2)
            pickle.dump(simulationsensor.EMTSorientationvectory, filehandler2)
            pickle.dump(simulationsensor.EMTSorientationvectorz, filehandler2)
            filehandler2.close()

            # Reset
            self.allacc = []
            self.count = 0
            self.graph.display()

    def combobox(
            self):  # This refers to the dropdown list with all the simple motions. Here they are pre-made so you don't have to retype them everytime.
        rotatex = [2 * math.pi, 0, 0]
        rotatexm = [-2 * math.pi, 0, 0]
        rotatey = [0, 2 * math.pi, 0]
        rotateym = [0, -2 * math.pi, 0]
        rotatez = [0, 0, 2 * math.pi]
        rotatezm = [0, 0, -2 * math.pi]
        accx = [1, 0, 0]
        accxm = [-1, 0, 0]
        accy = [0, 1, 0]
        accym = [0, -1, 0]
        accz = [0, 0, 1]
        acczm = [0, 0, -1]
        nochange = [0, 0, 0]
        if self.comboBox.currentText() == "Rotate x 90deg":
            timenew = 0.5
            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(rotatex) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([rotatex, 1, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(rotatexm) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([rotatexm, 1, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(nochange) + "x %i" % int(1) + ", %i" % self.count)
            self.allacc.append([nochange, 1, self.count])
        elif self.comboBox.currentText() == "Rotate x -90deg":
            timenew = 0.5
            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(rotatexm) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([rotatexm, 1, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(rotatex) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([rotatex, 1, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(nochange) + "x %i" % int(1) + ", %i" % self.count)
            self.allacc.append([nochange, 1, self.count])
        elif self.comboBox.currentText() == "Rotate y 90deg":
            timenew = 0.5
            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(rotatey) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([rotatey, 1, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(rotateym) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([rotateym, 1, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(nochange) + "x %i" % int(1) + ", %i" % self.count)
            self.allacc.append([nochange, 1, self.count])
        elif self.comboBox.currentText() == "Rotate y -90deg":
            timenew = 0.5
            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(rotateym) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([rotateym, 1, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(rotatey) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([rotatey, 1, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(nochange) + "x %i" % int(1) + ", %i" % self.count)
            self.allacc.append([nochange, 1, self.count])
        elif self.comboBox.currentText() == "Rotate z 90deg":
            timenew = 0.5
            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(rotatez) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([rotatez, 1, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(rotatezm) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([rotatezm, 1, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(nochange) + "x %i" % int(1) + ", %i" % self.count)
            self.allacc.append([nochange, 1, self.count])
        elif self.comboBox.currentText() == "Rotate z -90deg":
            timenew = 0.5
            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(rotatezm) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([rotatezm, 1, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(rotatez) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([rotatez, 1, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Rotational acceleration added" + str(nochange) + "x %i" % int(1) + ", %i" % self.count)
            self.allacc.append([nochange, 1, self.count])
        elif self.comboBox.currentText() == "Move x 2m":
            timenew = 1.
            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accx) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accx, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accxm) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accxm, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accxm) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accxm, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accx) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accx, 0, self.count])
        elif self.comboBox.currentText() == "Move -x 2m":
            timenew = 1.
            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accxm) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accxm, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accx) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accx, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accx) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accx, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accxm) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accxm, 0, self.count])
        elif self.comboBox.currentText() == "Move y 2m":
            timenew = 1.
            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accy) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accy, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accym) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accym, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accym) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accym, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accy) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accy, 0, self.count])
        elif self.comboBox.currentText() == "Move -y 2m":
            timenew = 1.
            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accym) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accym, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accy) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accy, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accy) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accy, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accym) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accym, 0, self.count])
        elif self.comboBox.currentText() == "Move z 2m":
            timenew = 1.
            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accz) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accz, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(acczm) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([acczm, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(acczm) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([acczm, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accz) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accz, 0, self.count])
        elif self.comboBox.currentText() == "Move -z 2m":
            timenew = 1.
            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(acczm) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([acczm, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accz) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accz, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(accz) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([accz, 0, self.count])

            self.count += 1
            self.listWidget.addItem(
                "Linear acceleration added:" + str(acczm) + "x %i" % int(timenew / self.dt) + ", %i" % self.count)
            for k in numpy.arange(int(timenew / self.dt)):
                self.allacc.append([acczm, 0, self.count])

    def linearacc(self):  # Button for User input linear jerk and time period in seconds
        text2, ok2 = QInputDialog.getText(dialog, 'User Input', 'Enter linear jerk (e.g. 0 0 1):')
        text2check = text2.split()
        if len(text2check) == 3:
            if isinstance(eval(text2check[0]), float) or isinstance(eval(text2check[0]), int) and isinstance(
                    eval(text2check[1]), float) or isinstance(eval(text2check[1]), int) and isinstance(
                    eval(text2check[2]), float) or isinstance(eval(text2check[2]),
                                                              int):  # If all 3 inputs are either floats or integers
                text2 = [eval(b) for b in text2.split()]  # input jerk separated into 3 components
                self.error.clear()
                time2, check2 = QInputDialog.getText(dialog, 'User input', 'Enter time (s) in steps of dt=%s' % self.dt)
                if len(time2.split()) == 1:
                    if isinstance(eval(time2.split()[0]), float) or isinstance(eval(time2.split()[0]),
                                                                               int):  # If the input is either a float or integer
                        timenew = eval(time2.split()[0])  # Time in seconds the jerk lasts for
                        self.count += 1  # Add one to counter for number of independent motions
                        self.listWidget.addItem("Linear acceleration added:" + str(text2) + "x %i" % int(
                            timenew / self.dt) + ", %i" % self.count)  # Display this motion to listwidget
                        for k in numpy.arange(
                                int(timenew / self.dt)):  # timenew/self.dt = no. iterations the jerk lasts for
                            self.allacc.append([text2, 0,
                                                self.count])  # Add jerk to list containing all user inputs. The 0 assigns the acceleration as linear.

            else:
                self.error.addText("Error: Please enter 3 linear \n acceleration values separated by spaces")
                self.graphicsView_2.show()
        else:
            self.error.addText("Error: Please enter 3 rotational acceleration values separated by spaces")
            self.graphicsView_2.show()

    def rotacc(
            self):  # Button for user input rotational acceleration and time period in seconds. Same method as that for linear jerk
        text, ok = QInputDialog.getText(dialog, 'User Input', 'Enter rotational acceleration (e.g. 1 0 0)')
        textcheck = text.split()
        if len(textcheck) == 3:
            if isinstance(eval(textcheck[0]), float) or isinstance(eval(textcheck[0]), int) and isinstance(
                    eval(textcheck[1]), float) or isinstance(eval(textcheck[1]), int) and isinstance(eval(textcheck[2]),
                                                                                                     float) or isinstance(
                    eval(textcheck[2]), int):
                text = [eval(b) for b in text.split()]
                self.error.clear()
                time, check = QInputDialog.getText(dialog, 'User input', 'Enter time (s) in steps of dt=%s' % self.dt)
                if len(time.split()) == 1:
                    if isinstance(eval(time.split()[0]), float) or isinstance(eval(time.split()[0]), int):
                        timenew = eval(time.split()[0])
                        self.count += 1
                        self.listWidget.addItem("Rotational acceleration added" + str(text) + "x %i" % int(
                            timenew / self.dt) + ", %i" % self.count)
                        for k in numpy.arange(int(timenew / self.dt)):
                            self.allacc.append([text, 1, self.count])  # the 1 assigns the acceleration as rotational

            else:
                self.error.addText("Error: Please enter 3 rotational acceleration values separated by spaces")
                self.graphicsView_2.show()
        else:
            self.error.addText("Error: Please enter 3 rotational acceleration values separated by spaces")
            self.graphicsView_2.show()

    def initorient(self):  # Set up initial orientation, required before doing any other motion
        self.count += 1
        self.allacc.append([[0, 0, 0], 1, self.count])
        self.listWidget.addItem("Initial orientation and position added: origin" + ", %i" % self.count)

    def deleteitem(self):  # Delete's an item off from List Widget
        item = self.listWidget.selectedItems()
        itemtype = self.listWidget.currentItem().text().split()
        indexval = int(itemtype[-1])  # The index is the self.count value
        deletelist = []
        for j in self.allacc:
            if j[2] == indexval:  # Find the matching self.count in Listwidget with the selected item
                deletelist.append(j)
            else:
                pass
        if len(deletelist) == 0:
            self.error.addText("Did not delete anything!")
        else:
            self.allacc = [x for x in self.allacc if x not in deletelist]  # Delete the selected item from ListWidget
        if not item:
            return
        for i in item:
            self.listWidget.takeItem(self.listWidget.row(i))

    def clearall(
            self):  # Attempt at resetting, but still problematic... Can tell when it does not work : it just crashes
        self.count = 0
        self.allacc = []
        self.graph.clear()
        self.listWidget.clear()
        self.error.clear()
        os.remove("simulatedsensor.obj")
        os.remove("EMTSsimulation.obj")


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    dialog = QtWidgets.QDialog()

    prog = SFGui(dialog)

    dialog.show()
    sys.exit(app.exec_())
