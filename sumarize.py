import sys
from PyQt5 import uic,QtWidgets
from PyQt5.QtWidgets import QLabel
from PyQt5.QtGui import QIcon


PATH, _ = uic.loadUiType('gui/sumarize.ui')


class Sumarize(QtWidgets.QDialog,PATH):
    """
    Show information from project file after creating new project
    """
    def __init__(self,rawheader, stationData, instrumentData, processingResults, gravityCorrections):
        super().__init__()
        self.setupUi(self)
        self.setWindowIcon(QIcon('picture/logo.svg'))

        self.rawheader = rawheader
        self.stationData = stationData
        self.instrumentData = instrumentData
        self.processingResults = processingResults
        self.gravityCorrections = gravityCorrections

        self.load()

        self.show()
        self.exec()

    def load(self):
        """
        Load dictionaries into text browsers
        """
        self.raw_tab.setText(self.raw(self.rawheader))
        self.station.setText(self.string(self.stationData))
        self.instrument.setText(self.string(self.instrumentData))
        self.processing.setText(self.string(self.processingResults))
        self.gravity.setText(self.string(self.gravityCorrections))

    def raw(self,rawheader):
        String=''
        for i in range(len(rawheader)):
            if i%2==1:
                n='\n'
                d=''
            else:
                n=''
                d='         '
            line=rawheader[i]+d+n
            String=String+line
        return String

    def string(self,dictionary):
        String=''
        for keys in dictionary.keys():
            line=keys+':        '+dictionary[keys]+'\n'
            String=String+line

        return String


if __name__ == "__main__":
    app=QtWidgets.QApplication([])
    sumarize=Sumarize()
    sumarize.show()
    app.exec_()
