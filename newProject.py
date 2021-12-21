import sys, os
from PyQt5 import  uic, QtWidgets
from glob import glob
from classes import projectFile, rawFile
from warning import Warning
from sumarize import Sumarize
from PyQt5.QtGui import QIcon
from CONFIG import logo_picture


# PATH, _ = uic.loadUiType('gui/newProject.ui')
PATH, _ = uic.loadUiType(os.path.join(os.path.dirname(__file__), 'gui/newProject.ui'))

class NewProject(QtWidgets.QDialog,PATH):
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

        self.show()
        self.exec()

    def load_files(self):
        """
        Get files direction
        """
        self.path=QtWidgets.QFileDialog.getExistingDirectory()
        # self.path = r'c:\Users\Jakub\Desktop\pecny\data\files_wetzel'

    def saveDir(self):
        """
        Get direction where will results
        """
        self.pathDir=QtWidgets.QFileDialog.getExistingDirectory()
        # self.pathDir = r'c:\Users\Jakub\Desktop\pecny\data\res'
        self.projDirPath.setText(self.pathDir)

    def accept(self):
        # for self.rawfilepath in glob('{}/*.raw.txt'.format(self.path)):
        #     rawfile=rawFile(self.rawfilepath)
        #     print(rawfile.rawHeader1())
        # self.pathDir = r'c:\Users\Jakub\Desktop\pecny\data\res'
        # self.path = r'c:\Users\Jakub\Desktop\pecny\data\files_wetzel'
        try:
            for self.rawfilepath in glob('{}/*.raw.txt'.format(self.path)):
                self.rawfile=rawFile(self.rawfilepath)
                header1=self.rawfile.rawHeader1()
                self.header2=self.rawfile.rawHeader2()
                self.header1=self.rawfile.rawHeader1()
                self.rawlines=self.rawfile.rawLines()


            for self.projectfile in glob('{}/*.project.txt'.format(self.path)):
                a=projectFile(self.projectfile)
                a.read()
                self.stationData, self.instrumentData, self.processingResults, self.gravityCorrections = a.createDictionary()

            self.setFile = glob('{}/*.set.txt'.format(self.path))[0]





            self.close()
            Sumarize(header1, self.stationData, self.instrumentData, self.processingResults, self.gravityCorrections)


        except AttributeError:
            Warning(error='Something is missing',icon='critical', title='Warning')

        except IndexError:
            Warning(error='Bad direction with files!',icon='critical', title='Warning')


        try:
            os.mkdir(self.pathDir+'/Graphs')
            os.mkdir(self.pathDir+'/Files')
        except FileExistsError:
            pass





if __name__ == "__main__":
    app=QtWidgets.QApplication([])
    main=NewProject()
    main.show()
    app.exec_()
