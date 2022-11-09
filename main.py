import threading
import winapps
from PyQt5 import uic, QtWidgets
from newProject import NewProject
from compute import Compute
from warning import Warning
from PyQt5.QtCore import *
from PyQt5.QtWidgets import QSplashScreen
from PyQt5.QtGui import QIcon, QPixmap
from data import Data
from viewgraphs import Graphs
from CONFIG import logo, wel, logo_picture, warning_window, picture_unchecked
from time import sleep
from functions import printDict, read_last_version, get_date_last_version
import os, subprocess,platform

script_path = os.path.dirname(os.path.realpath(__file__))

PATH, _ = uic.loadUiType(script_path + r'\gui\main.ui')


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
        pixmap = QPixmap(picture_unchecked)
        pixmap = pixmap.scaled(16, 16)
        self.check.setPixmap(pixmap)

        # connect buttons with method
        self.newProject.triggered.connect(self.new_Project)
        # self.Process_data_batch.triggered.connect(self.new_Project_batch)
        self.computing.triggered.connect(self.Computing)
        self.viewdata.triggered.connect(self.viewData)
        self.viewGraphs.triggered.connect(self.viewgraphs)
        self.closeProject.triggered.connect(self.closeproject)
        self.openProject.triggered.connect(self.open_project)

        self.newproject.clicked.connect(self.new_Project)

        self.finals_folder()

        self.last_version()

        # self.Import_data.triggered.connect(self.importdata)

    def last_version(self):
        last = get_date_last_version()
        file = read_last_version()

        if last != file:
            Warning(error=warning_window['commits'], icon='', title='Warning')

    def finals_folder(self):
        """
        This method just creates 'finals' folder
        @return:
        """
        try:
            os.mkdir('finals')
        except FileExistsError:
            pass

    def set_info(self):
        """
        Set some information on the main page
        :return:
        """

        self.campaign.setText(self.newProjectWin.stationData['ProjName'])
        self.sitecode.setText(self.newProjectWin.stationData['SiteCode'])
        self.trasnferheight.setText('{} {}'.format(self.newProjectWin.stationData['transferHeight'], self.newProjectWin.units['transferHeight']))
        self.gradient.setText('{} {}'.format(self.newProjectWin.stationData['gradient'], self.newProjectWin.units['gradient']))
        self.totaldrops.setText(str(int(self.newProjectWin.processingResults['setsCollected']) * int(self.newProjectWin.processingResults['dropsInSet'])))
        self.gravity.setText('{} {}'.format(self.newProjectWin.processingResults['gravity'], self.newProjectWin.units['gravity']))
        self.setscatter.setText('{} {}'.format(self.newProjectWin.processingResults['setScatter'], self.newProjectWin.units['setScatter']))

    @staticmethod
    def splashScreen():
        """
        Create splash screen with apple during starting of Agdas
        :return:
        """
        splash_pix = QPixmap(logo_picture)
        splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
        splash.setWindowFlags(Qt.WindowStaysOnTopHint | Qt.FramelessWindowHint)
        splash.show()
        splash.showMessage("<h1><font color='white'>Welcome to the Agdas!</font></h1>", Qt.AlignTop | Qt.AlignCenter,
                           Qt.black)
        sleep(1)
        splash.hide()

    def open_project(self):
        """

        :return:
        """
        path = QtWidgets.QFileDialog.getExistingDirectory()
        self.newProjectWin = NewProject().accept()

    def new_Project(self):
        """
        This method serves for creating of new project
        :return:
        """
        # Open newproject dialog
        self.newProjectWin = NewProject()

        # If newproject was created successfully
        if self.newProjectWin.succ:
            # Load check picture to main page
            pixmap = QPixmap(os.path.dirname(os.path.realpath(__file__)) + r'\picture\checked.png')
            self.check.setPixmap(pixmap)

            # Set name of campaign
            self.set_info()

            # Info about successfully upload
            self.succ_upload.setText('Following input files have been uploaded successfully!')

            self.label_2.clear()

    def Computing(self):
        """
        OThis method serves for creating of Compute
        :return:
        """
        try:
            # Open Compute dialog
            self.result = Compute(self.newProjectWin.path, stationData=self.newProjectWin.stationData,
                                  instrumentData=self.newProjectWin.instrumentData,
                                  processingResults=self.newProjectWin.processingResults,
                                  gravityCorrections=self.newProjectWin.gravityCorrections,
                                  header2=self.newProjectWin.header2, rawlines=self.newProjectWin.rawlines,
                                  header1=self.newProjectWin.header1, projDirPath=self.newProjectWin.pathDir,
                                  setFile=self.newProjectWin.setFile)
            # self.result.Run()

        except AttributeError:
            Warning(error=warning_window['import_data'], icon='critical', title='Warning')

    def viewData(self):
        """
        This method is for viewing of processed data by DB Browser for SQLite.
        :return:
        """

        try:
            # For windows system
            if 'Win' in platform.platform():
                # Find install location of DB browser
                for item in winapps.list_installed():
                    if 'DB Browser' in item.name:
                        browser_path = item.install_location
                        browser_name = item.name

                # Open DB browser in another thread
                command = ["{}".format(os.path.join(browser_path, browser_name)), "{}".format(os.path.join(self.newProjectWin.pathDir, 'data.db'))]
                subprocess.Popen(command)
            else:
                Warning(error='Support for Linux is in process!', icon='critical', title='Warning')

        # If DB browser is not installed
        except UnboundLocalError:
            url_link = "<a href=\"https://sqlitebrowser.org/dl/\">Click this link to go to DB Browser for SQLite " \
                       "Download</a> "
            Warning(error=url_link, icon='critical', title='Warning')

        # If project is not opened
        except AttributeError:
            Warning(error='Opened project is necessary!', icon='critical', title='Warning')

    def viewgraphs(self):
        try:
            Graphs(self.newProjectWin.pathDir + '/Graphs')
        except AttributeError:
            Warning(error=warning_window['import_data'], icon='critical', title='Warning')

    def closeproject(self):
        """
        This method delete self.newProjectWin variable and change main page
        :return:
        """
        try:
            del self.newProjectWin
            pixmap = QPixmap(os.path.dirname(os.path.realpath(__file__)) + r'\picture\unchecked.png')
            pixmap = pixmap.scaled(16, 16)
            self.check.setPixmap(pixmap)

            self.campaign.clear()

            self.succ_upload.clear()

            self.label_2.setText('Upload data for measurement campaign')

        except AttributeError:
            Warning(error=warning_window['project'], icon='critical', title='Warning')
        # self.__del__()

    # def __del__(self):
    #     pass


if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    main = Main()
    main.show()
    app.exec_()
