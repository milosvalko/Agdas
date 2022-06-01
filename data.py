from PyQt5 import uic, QtWidgets
from PyQt5.QtSql import QSqlDatabase, QSqlQuery, QSqlQueryModel
from PyQt5.QtGui import QIcon
from CONFIG import logo_picture
import os

script_path = os.path.dirname(os.path.realpath(__file__))

PATH, _ = uic.loadUiType(script_path + '\gui\data.ui')


class Data(QtWidgets.QDialog, PATH):

    def __init__(self, matr):
        super().__init__()
        self.setupUi(self)
        self.setWindowIcon(QIcon(logo_picture))

        self.matr_db = matr

        # open database with results
        db = QSqlDatabase.addDatabase("QSQLITE", "db")
        db.setDatabaseName(self.matr_db + '/data.db')
        db.open()

        # create model and load data
        projectModel = QSqlQueryModel()
        projectModel.setQuery('select * from results', db)
        self.view.setModel(projectModel)

        db.close()

        self.show()
        self.exec()


if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    data = Data()
    data.show()
    app.exec_()
