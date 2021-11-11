import os
import sys
from PyQt5 import uic,QtWidgets
from newProject import NewProject
from compute import Compute
from warning import Warning
from PyQt5.QtCore import *
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QTreeView, QSplashScreen
from PyQt5.QtGui import QIcon, QPixmap, QMovie
from data import Data
from viewgraphs import Graphs
from CONFIG import logo, wel, logo_picture, warning_window
from time import sleep


PATH, _ = uic.loadUiType('gui/main.ui')

class Main(QtWidgets.QMainWindow, PATH):
    """
    Main page of Agdas
    """
    def __init__(self):
        super().__init__()
        self.splashScreen()
        print(logo)
        print(wel)
        self.setupUi(self)
        self.setWindowIcon(QIcon(logo_picture))

        # show opened/closed picture
        pixmap = QPixmap('picture/unchecked.png')
        pixmap=pixmap.scaled(16,16)
        self.check.setPixmap(pixmap)

        # connect buttons with method
        self.newProject.triggered.connect(self.new_Project)
        self.computing.triggered.connect(self.Computing)
        self.viewdata.triggered.connect(self.viewData)
        self.viewGraphs.triggered.connect(self.viewgraphs)
        self.closeProject.triggered.connect(self.closeproject)
        self.openProject.triggered.connect(self.open_project)

        # self.Import_data.triggered.connect(self.importdata)

    def splashScreen(self):
        splash_pix = QPixmap(logo_picture)
        splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
        splash.setWindowFlags(Qt.WindowStaysOnTopHint | Qt.FramelessWindowHint)
        splash.show()
        # splash.showMessage("<h1><font color='white'>Welcome to the Agdas!</font></h1>", Qt.AlignTop | Qt.AlignCenter, Qt.black)
        sleep(1)
        splash.hide()

    def open_project(self):

        path=QtWidgets.QFileDialog.getExistingDirectory()
        self.newProjectWin=NewProject().accept()

    def new_Project(self):
        # show opened/closed picture
        self.newProjectWin=NewProject()
        try:
            self.newProjectWin.pathDir
            pixmap = QPixmap('picture/checked.png')
            self.check.setPixmap(pixmap)
        except AttributeError:
            pass

    def Computing(self):
        try:
            self.result=Compute(self.newProjectWin.path, stationData=self.newProjectWin.stationData, instrumentData=self.newProjectWin.instrumentData,
        processingResults=self.newProjectWin.processingResults, gravityCorrections=self.newProjectWin.gravityCorrections,
        header2=self.newProjectWin.header2,rawlines=self.newProjectWin.rawlines, header1=self.newProjectWin.header1, projDirPath=self.newProjectWin.pathDir, setFile=self.newProjectWin.setFile)
        except AttributeError:
            Warning(error=warning_window['import_data'],icon='critical', title='Warning')

    def viewData(self):

        try:
            Data(self.newProjectWin.pathDir)
        except AttributeError:
            Warning(error=warning_window['import_data'],icon='critical', title='Warning')

    def viewgraphs(self):
        try:
            Graphs(self.newProjectWin.pathDir+'/Graphs')
        except AttributeError:
            Warning(error=warning_window['import_data'],icon='critical', title='Warning')

    def closeproject(self):
        try:
            del self.newProjectWin
            pixmap = QPixmap('picture/unchecked.png')
            pixmap=pixmap.scaled(16,16)
            self.check.setPixmap(pixmap)
        except AttributeError:
            Warning(error=warning_window['project'],icon='critical', title='Warning')
        # self.__del__()


    # def __del__(self):
    #     pass


if __name__ == "__main__":
    app=QtWidgets.QApplication([])
    main=Main()
    main.show()
    app.exec_()
