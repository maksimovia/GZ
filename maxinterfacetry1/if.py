# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'if.ui'
#
# Created by: PyQt5 UI code generator 5.15.7
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(822, 266)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.pict = QtWidgets.QLabel(self.centralwidget)
        self.pict.setGeometry(QtCore.QRect(0, 0, 701, 261))
        self.pict.setText("")
        self.pict.setPixmap(QtGui.QPixmap("схема.png"))
        self.pict.setObjectName("pict")
        self.Tmark1 = QtWidgets.QLabel(self.centralwidget)
        self.Tmark1.setGeometry(QtCore.QRect(0, 190, 47, 13))
        self.Tmark1.setObjectName("Tmark1")
        self.T_yx_gas_value7 = QtWidgets.QLabel(self.centralwidget)
        self.T_yx_gas_value7.setGeometry(QtCore.QRect(0, 210, 47, 13))
        self.T_yx_gas_value7.setObjectName("T_yx_gas_value7")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(110, 190, 47, 13))
        self.label.setObjectName("label")
        self.T_yx_gas_value6 = QtWidgets.QLabel(self.centralwidget)
        self.T_yx_gas_value6.setGeometry(QtCore.QRect(110, 200, 47, 13))
        self.T_yx_gas_value6.setObjectName("T_yx_gas_value6")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(180, 190, 47, 13))
        self.label_2.setObjectName("label_2")
        self.T_yx_gas_value5 = QtWidgets.QLabel(self.centralwidget)
        self.T_yx_gas_value5.setGeometry(QtCore.QRect(180, 200, 47, 13))
        self.T_yx_gas_value5.setObjectName("T_yx_gas_value5")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(250, 190, 47, 13))
        self.label_3.setObjectName("label_3")
        self.T_yx_gas_value4 = QtWidgets.QLabel(self.centralwidget)
        self.T_yx_gas_value4.setGeometry(QtCore.QRect(250, 200, 47, 13))
        self.T_yx_gas_value4.setObjectName("T_yx_gas_value4")
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(310, 190, 47, 13))
        self.label_4.setObjectName("label_4")
        self.T_yx_gas_value3 = QtWidgets.QLabel(self.centralwidget)
        self.T_yx_gas_value3.setGeometry(QtCore.QRect(310, 200, 47, 13))
        self.T_yx_gas_value3.setObjectName("T_yx_gas_value3")
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(390, 190, 47, 13))
        self.label_5.setObjectName("label_5")
        self.T_yx_gas_value2 = QtWidgets.QLabel(self.centralwidget)
        self.T_yx_gas_value2.setGeometry(QtCore.QRect(390, 200, 47, 13))
        self.T_yx_gas_value2.setObjectName("T_yx_gas_value2")
        self.label_6 = QtWidgets.QLabel(self.centralwidget)
        self.label_6.setGeometry(QtCore.QRect(480, 190, 47, 13))
        self.label_6.setObjectName("label_6")
        self.T_yx_gas_value1 = QtWidgets.QLabel(self.centralwidget)
        self.T_yx_gas_value1.setGeometry(QtCore.QRect(480, 200, 47, 13))
        self.T_yx_gas_value1.setObjectName("T_yx_gas_value1")
        self.label_7 = QtWidgets.QLabel(self.centralwidget)
        self.label_7.setGeometry(QtCore.QRect(710, 120, 71, 16))
        self.label_7.setObjectName("label_7")
        self.T_yx_gas_value1_2 = QtWidgets.QLabel(self.centralwidget)
        self.T_yx_gas_value1_2.setGeometry(QtCore.QRect(710, 140, 47, 13))
        self.T_yx_gas_value1_2.setObjectName("T_yx_gas_value1_2")
        self.pushButton_2 = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_2.setGeometry(QtCore.QRect(720, 70, 91, 41))
        self.pushButton_2.setStyleSheet("background-color: rgb(255, 0, 0);")
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton.setGeometry(QtCore.QRect(720, 10, 91, 41))
        self.pushButton.setStyleSheet("background-color: rgb(92, 255, 67);")
        self.pushButton.setObjectName("pushButton")
        self.label_8 = QtWidgets.QLabel(self.centralwidget)
        self.label_8.setGeometry(QtCore.QRect(70, 60, 47, 13))
        self.label_8.setObjectName("label_8")
        self.T_GPK_IN = QtWidgets.QLabel(self.centralwidget)
        self.T_GPK_IN.setGeometry(QtCore.QRect(70, 70, 47, 13))
        self.T_GPK_IN.setObjectName("T_GPK_IN")
        self.label_9 = QtWidgets.QLabel(self.centralwidget)
        self.label_9.setGeometry(QtCore.QRect(70, 20, 61, 16))
        self.label_9.setObjectName("label_9")
        self.T_GPK_OUT_2 = QtWidgets.QLabel(self.centralwidget)
        self.T_GPK_OUT_2.setGeometry(QtCore.QRect(70, 30, 47, 13))
        self.T_GPK_OUT_2.setObjectName("T_GPK_OUT_2")
        self.label_10 = QtWidgets.QLabel(self.centralwidget)
        self.label_10.setGeometry(QtCore.QRect(0, 90, 81, 16))
        self.label_10.setObjectName("label_10")
        self.G_GPK = QtWidgets.QLabel(self.centralwidget)
        self.G_GPK.setGeometry(QtCore.QRect(0, 110, 47, 13))
        self.G_GPK.setObjectName("G_GPK")
        self.label_11 = QtWidgets.QLabel(self.centralwidget)
        self.label_11.setGeometry(QtCore.QRect(140, 60, 61, 16))
        self.label_11.setObjectName("label_11")
        self.T_IND_IN = QtWidgets.QLabel(self.centralwidget)
        self.T_IND_IN.setGeometry(QtCore.QRect(140, 70, 47, 13))
        self.T_IND_IN.setObjectName("T_IND_IN")
        self.T_IND_OUT = QtWidgets.QLabel(self.centralwidget)
        self.T_IND_OUT.setGeometry(QtCore.QRect(140, 30, 47, 13))
        self.T_IND_OUT.setObjectName("T_IND_OUT")
        self.label_12 = QtWidgets.QLabel(self.centralwidget)
        self.label_12.setGeometry(QtCore.QRect(140, 20, 61, 16))
        self.label_12.setObjectName("label_12")
        self.T_PPND_IN = QtWidgets.QLabel(self.centralwidget)
        self.T_PPND_IN.setGeometry(QtCore.QRect(210, 70, 47, 13))
        self.T_PPND_IN.setObjectName("T_PPND_IN")
        self.label_13 = QtWidgets.QLabel(self.centralwidget)
        self.label_13.setGeometry(QtCore.QRect(210, 60, 61, 16))
        self.label_13.setObjectName("label_13")
        self.label_14 = QtWidgets.QLabel(self.centralwidget)
        self.label_14.setGeometry(QtCore.QRect(210, 20, 61, 16))
        self.label_14.setObjectName("label_14")
        self.T_PPND_OUT = QtWidgets.QLabel(self.centralwidget)
        self.T_PPND_OUT.setGeometry(QtCore.QRect(210, 30, 47, 13))
        self.T_PPND_OUT.setObjectName("T_PPND_OUT")
        self.T_VE_IN = QtWidgets.QLabel(self.centralwidget)
        self.T_VE_IN.setGeometry(QtCore.QRect(280, 70, 47, 13))
        self.T_VE_IN.setObjectName("T_VE_IN")
        self.label_15 = QtWidgets.QLabel(self.centralwidget)
        self.label_15.setGeometry(QtCore.QRect(280, 60, 61, 16))
        self.label_15.setObjectName("label_15")
        self.label_16 = QtWidgets.QLabel(self.centralwidget)
        self.label_16.setGeometry(QtCore.QRect(280, 20, 61, 16))
        self.label_16.setObjectName("label_16")
        self.T_VE_OUT = QtWidgets.QLabel(self.centralwidget)
        self.T_VE_OUT.setGeometry(QtCore.QRect(280, 30, 47, 13))
        self.T_VE_OUT.setObjectName("T_VE_OUT")
        self.T_IVD_IN = QtWidgets.QLabel(self.centralwidget)
        self.T_IVD_IN.setGeometry(QtCore.QRect(350, 70, 47, 13))
        self.T_IVD_IN.setObjectName("T_IVD_IN")
        self.label_17 = QtWidgets.QLabel(self.centralwidget)
        self.label_17.setGeometry(QtCore.QRect(350, 60, 61, 16))
        self.label_17.setObjectName("label_17")
        self.label_18 = QtWidgets.QLabel(self.centralwidget)
        self.label_18.setGeometry(QtCore.QRect(350, 20, 61, 16))
        self.label_18.setObjectName("label_18")
        self.T_IVD_OUT = QtWidgets.QLabel(self.centralwidget)
        self.T_IVD_OUT.setGeometry(QtCore.QRect(350, 30, 47, 13))
        self.T_IVD_OUT.setObjectName("T_IVD_OUT")
        self.T_PPVD_IN = QtWidgets.QLabel(self.centralwidget)
        self.T_PPVD_IN.setGeometry(QtCore.QRect(420, 70, 47, 13))
        self.T_PPVD_IN.setObjectName("T_PPVD_IN")
        self.label_19 = QtWidgets.QLabel(self.centralwidget)
        self.label_19.setGeometry(QtCore.QRect(420, 60, 61, 16))
        self.label_19.setObjectName("label_19")
        self.label_20 = QtWidgets.QLabel(self.centralwidget)
        self.label_20.setGeometry(QtCore.QRect(420, 20, 61, 16))
        self.label_20.setObjectName("label_20")
        self.T_PPVD_OUT = QtWidgets.QLabel(self.centralwidget)
        self.T_PPVD_OUT.setGeometry(QtCore.QRect(420, 30, 47, 13))
        self.T_PPVD_OUT.setObjectName("T_PPVD_OUT")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 822, 21))
        self.menubar.setDefaultUp(False)
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.Tmark1.setText(_translate("MainWindow", "T ух.газ. "))
        self.T_yx_gas_value7.setText(_translate("MainWindow", "Значение, С"))
        self.label.setText(_translate("MainWindow", "T ух. газ"))
        self.T_yx_gas_value6.setText(_translate("MainWindow", "Значение, С"))
        self.label_2.setText(_translate("MainWindow", "T ух. газ"))
        self.T_yx_gas_value5.setText(_translate("MainWindow", "Значение, С"))
        self.label_3.setText(_translate("MainWindow", "T ух. газ"))
        self.T_yx_gas_value4.setText(_translate("MainWindow", "Значение, С"))
        self.label_4.setText(_translate("MainWindow", "T ух. газ"))
        self.T_yx_gas_value3.setText(_translate("MainWindow", "Значение, С"))
        self.label_5.setText(_translate("MainWindow", "T ух. газ"))
        self.T_yx_gas_value2.setText(_translate("MainWindow", "Значение, С"))
        self.label_6.setText(_translate("MainWindow", "T ух. газ"))
        self.T_yx_gas_value1.setText(_translate("MainWindow", "Значение, С"))
        self.label_7.setText(_translate("MainWindow", "T окр среды"))
        self.T_yx_gas_value1_2.setText(_translate("MainWindow", "Нагрузка, %"))
        self.pushButton_2.setText(_translate("MainWindow", "Стоп"))
        self.pushButton.setText(_translate("MainWindow", "Старт"))
        self.label_8.setText(_translate("MainWindow", "T ГПК вх"))
        self.T_GPK_IN.setText(_translate("MainWindow", "Значение"))
        self.label_9.setText(_translate("MainWindow", "T ГПК вых"))
        self.T_GPK_OUT_2.setText(_translate("MainWindow", "Значение"))
        self.label_10.setText(_translate("MainWindow", "Расход"))
        self.G_GPK.setText(_translate("MainWindow", "Значение"))
        self.label_11.setText(_translate("MainWindow", "T ИНД вх"))
        self.T_IND_IN.setText(_translate("MainWindow", "Значение"))
        self.T_IND_OUT.setText(_translate("MainWindow", "Значение"))
        self.label_12.setText(_translate("MainWindow", "T ИНД вых"))
        self.T_PPND_IN.setText(_translate("MainWindow", "Значение"))
        self.label_13.setText(_translate("MainWindow", "T ППНД вх"))
        self.label_14.setText(_translate("MainWindow", "T ППНД вых"))
        self.T_PPND_OUT.setText(_translate("MainWindow", "Значение"))
        self.T_VE_IN.setText(_translate("MainWindow", "Значение"))
        self.label_15.setText(_translate("MainWindow", "T ВЭ вх"))
        self.label_16.setText(_translate("MainWindow", "T ВЭ вых"))
        self.T_VE_OUT.setText(_translate("MainWindow", "Значение"))
        self.T_IVD_IN.setText(_translate("MainWindow", "Значение"))
        self.label_17.setText(_translate("MainWindow", "T ИВД вх"))
        self.label_18.setText(_translate("MainWindow", "T ИВД вых"))
        self.T_IVD_OUT.setText(_translate("MainWindow", "Значение"))
        self.T_PPVD_IN.setText(_translate("MainWindow", "Значение"))
        self.label_19.setText(_translate("MainWindow", "T ППВД вх"))
        self.label_20.setText(_translate("MainWindow", "T ППВД вых"))
        self.T_PPVD_OUT.setText(_translate("MainWindow", "Значение"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
