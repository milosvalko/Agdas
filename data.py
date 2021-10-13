from PyQt5 import uic,QtWidgets
from PyQt5.QtSql import QSqlDatabase,QSqlQuery,QSqlQueryModel
from PyQt5.QtGui import QIcon

PATH, _ = uic.loadUiType('gui/data.ui')

class Data(QtWidgets.QDialog,PATH):

    def __init__(self, matr=r'c:\Users\Jakub\Desktop\pecny\pyAgdasGui\data\res'):
        super().__init__()
        self.setupUi(self)
        self.setWindowIcon(QIcon('picture/logo.svg'))

        self.matr_db=matr

        # open database with results
        db = QSqlDatabase.addDatabase("QSQLITE","db")
        db.setDatabaseName(self.matr_db+'/data.db')
        db.open()

        # create model and load data
        projectModel = QSqlQueryModel()
        projectModel.setQuery('select * from results',db)
        self.view.setModel(projectModel)

        db.close()

        self.show()
        self.exec()


if __name__ == "__main__":
    app=QtWidgets.QApplication([])
    data=Data()
    data.show()
    app.exec_()
