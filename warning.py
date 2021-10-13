import sys
from PyQt5 import uic,QtWidgets
from PyQt5.QtWidgets import QLabel

class Warning():
    """Generate warning window with 3 parameters:
    error - definition of problem
    icon - picture in window
    title - name of window
    """
    def __init__(self, error, icon, title):

        msgBox=QtWidgets.QMessageBox()
        if icon=='critical':
            msgBox.setIcon(msgBox.Critical)

        msgBox.setText(error)
        msgBox.setWindowTitle(title)
        msgBox.exec()

if __name__ == "__main__":
    app=QtWidgets.QApplication([])
    warning=Warning()
    warning.show()
    app.exec_()
