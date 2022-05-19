import sys, os
from PyQt5 import uic, QtWidgets
from glob import glob
from classes import projectFile, rawFile
from warning import Warning
from PyQt5.QtGui import QIcon
from CONFIG import logo_picture, warning_window
from compare import Compare

# PATH, _ = uic.loadUiType('gui/newProject.ui')
PATH, _ = uic.loadUiType(os.path.join(os.path.dirname(__file__), r'gui\comparison.ui'))

class Comparison(QtWidgets.QDialog, PATH):
    """
    Compare campagnes
    """

    def __init__(self, delimiter =','):
        super().__init__()
        self.setupUi(self)
        self.setWindowIcon(QIcon(logo_picture))

        # # connect buttons with methods
        self.source.clicked.connect(self.source_path1)
        self.source2.clicked.connect(self.source_path2)
        self.target.clicked.connect(self.target_path)

        self.delimiter = delimiter

        self.show()
        self.exec()

    def source_path1(self):
        self.path1 = QtWidgets.QFileDialog.getExistingDirectory()

    def source_path2(self):
        self.path2 = QtWidgets.QFileDialog.getExistingDirectory()

    def target_path(self):
        self.comparison_path = QtWidgets.QFileDialog.getExistingDirectory()

    def accept(self):

        try:
            Compare(self.path1, self.path2, self.comparison_path, delimiter=self.delimiter)
            Warning(error=warning_window['comparison'], icon='information', title='Done')
        except AttributeError:
            Warning(error='Something is missing!', icon='warning', title='Warning')
        self.close()


if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    main = Comparison()
    main.show()
    app.exec_()