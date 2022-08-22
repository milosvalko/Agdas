import sys, os
from PyQt5 import uic, QtWidgets
from glob import glob
from classes import projectFile, rawFile
from warning import Warning
from sumarize import Sumarize
from PyQt5.QtGui import QIcon
from CONFIG import logo_picture

# PATH, _ = uic.loadUiType('gui/newProject.ui')
PATH, _ = uic.loadUiType(os.path.join(os.path.dirname(__file__), 'gui/newProject.ui'))


class NewProject(QtWidgets.QDialog, PATH):
    """
    Creating new project
    """

    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.setWindowIcon(QIcon(logo_picture))

        # connect buttons with methods
        self.loadFiles.clicked.connect(self.load_files)
        self.projDir.clicked.connect(self.saveDir)

        # info if upload of input files was successfully
        self.succ = False

        self.show()
        self.exec()

    def load_files(self):
        """
        Get files direction
        """
        self.path = QtWidgets.QFileDialog.getExistingDirectory()
        self.name.setText(self.path)

    def saveDir(self):
        """
        Get direction where will results
        """
        self.pathDir = QtWidgets.QFileDialog.getExistingDirectory()
        self.projDirPath.setText(self.pathDir)

    def accept(self):

        if len(self.name.toPlainText()) > 0 and len(self.projDirPath.toPlainText()) > 0:
            try:

                self.path = self.name.toPlainText()
                self.pathDir = self.projDirPath.toPlainText()

                for self.rawfilepath in glob('{}/*.raw.txt'.format(self.path)):
                    self.rawfile = rawFile(self.rawfilepath)
                    header1 = self.rawfile.rawHeader1()
                    self.header2 = self.rawfile.rawHeader2()
                    self.header1 = self.rawfile.rawHeader1()
                    self.rawlines = self.rawfile.rawLines()

                for self.projectfile in glob('{}/*.project.txt'.format(self.path)):
                    a = projectFile(self.projectfile)
                    a.read()
                    self.stationData, self.instrumentData, self.processingResults, self.gravityCorrections, self.names, self.units = a.createDictionary()

                self.setFile = glob('{}/*.set.txt'.format(self.path))[0]

                self.close()

                Sumarize(header1, self.stationData, self.instrumentData, self.processingResults, self.gravityCorrections, self.names, self.units)

                self.succ = True

                try:
                    os.mkdir(self.pathDir + '/Graphs')
                    os.mkdir(self.pathDir + '/Files')
                except FileExistsError:
                    pass

            except IndexError:
                Warning(error='Bad direction with files!', icon='critical', title='Warning')
        else:
            Warning(error='Something is missing', icon='critical', title='Warning')


if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    main = NewProject()
    main.show()
    app.exec_()
