# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'sftest.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(817, 688)
        self.pushButton = QtWidgets.QPushButton(Dialog)
        self.pushButton.setGeometry(QtCore.QRect(40, 180, 161, 23))
        self.pushButton.setObjectName("pushButton")
        self.pushButton_2 = QtWidgets.QPushButton(Dialog)
        self.pushButton_2.setGeometry(QtCore.QRect(40, 50, 161, 23))
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_3 = QtWidgets.QPushButton(Dialog)
        self.pushButton_3.setGeometry(QtCore.QRect(40, 80, 161, 23))
        self.pushButton_3.setObjectName("pushButton_3")
        self.listWidget = QtWidgets.QListWidget(Dialog)
        self.listWidget.setGeometry(QtCore.QRect(210, 20, 321, 181))
        self.listWidget.setObjectName("listWidget")
        self.pushButton_6 = QtWidgets.QPushButton(Dialog)
        self.pushButton_6.setGeometry(QtCore.QRect(40, 20, 161, 23))
        self.pushButton_6.setObjectName("pushButton_6")
        self.progressBar = QtWidgets.QProgressBar(Dialog)
        self.progressBar.setGeometry(QtCore.QRect(10, 660, 118, 23))
        self.progressBar.setProperty("value", 24)
        self.progressBar.setObjectName("progressBar")
        self.graphicsView = QtWidgets.QGraphicsView(Dialog)
        self.graphicsView.setGeometry(QtCore.QRect(40, 220, 731, 421))
        self.graphicsView.setObjectName("graphicsView")
        self.pushButton_4 = QtWidgets.QPushButton(Dialog)
        self.pushButton_4.setGeometry(QtCore.QRect(40, 110, 161, 23))
        self.pushButton_4.setObjectName("pushButton_4")
        self.graphicsView_2 = QtWidgets.QGraphicsView(Dialog)
        self.graphicsView_2.setGeometry(QtCore.QRect(540, 20, 261, 51))
        self.graphicsView_2.setObjectName("graphicsView_2")
        self.pushButton_5 = QtWidgets.QPushButton(Dialog)
        self.pushButton_5.setGeometry(QtCore.QRect(40, 140, 161, 23))
        self.pushButton_5.setObjectName("pushButton_5")
        self.comboBox = QtWidgets.QComboBox(Dialog)
        self.comboBox.setGeometry(QtCore.QRect(540, 80, 161, 21))
        self.comboBox.setObjectName("comboBox")
        self.comboBox.addItem("Rotate x 90deg")
        self.comboBox.addItem("Rotate x -90deg")
        self.comboBox.addItem("Rotate y 90deg")
        self.comboBox.addItem("Rotate y -90deg")
        self.comboBox.addItem("Rotate z 90deg")
        self.comboBox.addItem("Rotate z -90deg")
        self.comboBox.addItem("Move x 2m")
        self.comboBox.addItem("Move -x 2m")
        self.comboBox.addItem("Move y 2m")
        self.comboBox.addItem("Move -y 2m")
        self.comboBox.addItem("Move z 2m")
        self.comboBox.addItem("Move z -2m")
        self.comboBox.addItem("Select simple motion")

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.pushButton.setText(_translate("Dialog", "Run Simulation"))
        self.pushButton_2.setText(_translate("Dialog", "Add linear jerk"))
        self.pushButton_3.setText(_translate("Dialog", "Add rotational acceleration"))
        self.pushButton_6.setText(_translate("Dialog", "Add initial orientation"))
        self.pushButton_4.setText(_translate("Dialog", "Delete item"))
        self.pushButton_5.setText(_translate("Dialog", "Clear all"))
        self.comboBox.setItemText(0, _translate("Dialog", "Select simple motion"))
        self.comboBox.setItemText(1, _translate("Dialog", "Rotate x 90deg"))
        self.comboBox.setItemText(2, _translate("Dialog", "Rotate x -90deg"))
        self.comboBox.setItemText(3, _translate("Dialog", "Rotate y 90deg"))
        self.comboBox.setItemText(4, _translate("Dialog", "Rotate y -90deg"))
        self.comboBox.setItemText(5, _translate("Dialog", "Rotate z 90deg"))
        self.comboBox.setItemText(6, _translate("Dialog", "Rotate z -90deg"))
        self.comboBox.setItemText(7, _translate("Dialog", "Move x 2m"))
        self.comboBox.setItemText(8, _translate("Dialog", "Move -x 2m"))
        self.comboBox.setItemText(9, _translate("Dialog", "Move y 2m"))
        self.comboBox.setItemText(10, _translate("Dialog", "Move -y 2m"))
        self.comboBox.setItemText(11, _translate("Dialog", "Move z 2m"))
        self.comboBox.setItemText(12, _translate("Dialog", "Move -z 2m"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())

