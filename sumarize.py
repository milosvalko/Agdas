import sys
from PyQt5 import uic, QtWidgets
from PyQt5.QtWidgets import QLabel
from PyQt5.QtGui import QIcon
from CONFIG import logo_picture
import os

script_path = os.path.dirname(os.path.realpath(__file__))

PATH, _ = uic.loadUiType(script_path + '\gui\sumarize.ui')

class Sumarize(QtWidgets.QDialog, PATH):
    """
    Show information from project file after creating new project
    """

    def __init__(self, rawheader, stationData, instrumentData, processingResults, gravityCorrections, names, units):
        super().__init__()
        self.setupUi(self)
        self.setWindowIcon(QIcon(logo_picture))

        self.rawheader = rawheader
        self.stationData = stationData
        self.instrumentData = instrumentData
        self.processingResults = processingResults
        self.gravityCorrections = gravityCorrections
        self.names = names
        self.units = units

        self.accept_button.clicked.connect(self.accept)

        self.load()

        self.show()
        self.exec()

    def accept(self):
        self.close()

    def load(self):
        """
        Load dictionaries into text browsers
        """
        self.station.setText(self.string(self.stationData))
        self.instrument.setText(self.string(self.instrumentData))
        self.processing.setText(self.string(self.processingResults))
        self.gravity.setText(self.string(self.gravityCorrections))


    def string(self, dictionary):
        String = ''
        for keys in dictionary.keys():
            try:
                unit = self.units[keys]
            except KeyError:
                unit = ' '

            line = self.names[keys] + ' ' + dictionary[keys] + ' ' + unit + '\n'
            String = String + line

        return String


if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    sumarize = Sumarize()
    sumarize.show()
    app.exec_()
