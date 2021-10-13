from PyQt5 import uic,QtWidgets
from PyQt5.QtSql import QSqlDatabase,QSqlQuery,QSqlQueryModel
import os
from PyQt5.QtGui import QIcon, QPixmap

PATH, _ = uic.loadUiType('gui/graph.ui')

class Graphs(QtWidgets.QDialog,PATH):

    def __init__(self, res_path):
        super().__init__()
        self.setupUi(self)
        self.setWindowIcon(QIcon('picture/logo.svg'))
        self.res_path=res_path

        self.graphs_list.addItems(os.listdir(res_path))

        self.graphs_list.currentItemChanged.connect(self.showG)
        self.reload.clicked.connect(self.Reload)



        self.show()
        self.exec()

    def showG(self):
        gr=self.graphs_list.currentItem().text()
        # gr='/effective_height.png'


        pixmap = QPixmap(self.res_path+'/'+gr)
        pixmap=pixmap.scaled(731,591)
        self.graph_label.setPixmap(pixmap)

    def Reload(self):
        self.__init__(self.res_path)
        self.close()


if __name__ == "__main__":
    app=QtWidgets.QApplication([])
    data=Graphs()
    data.show()
    app.exec_()
