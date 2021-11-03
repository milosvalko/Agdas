from PyQt5 import uic,QtWidgets
from PyQt5.QtSql import QSqlDatabase,QSqlQuery,QSqlQueryModel
import os
from PyQt5.QtGui import QIcon, QPixmap
from CONFIG import logo_picture

PATH, _ = uic.loadUiType('gui/graph.ui')

class Graphs(QtWidgets.QDialog,PATH):

    def __init__(self, res_path):
        super().__init__()
        self.setupUi(self)
        self.setWindowIcon(QIcon(logo_picture))
        self.res_path=res_path

        self.setItem()

        self.graphs_list.currentItemChanged.connect(self.showG)
        self.reload.clicked.connect(self.Reload)



        self.show()
        self.exec()

    def setItem(self):

        it=[]
        for i in os.listdir(self.res_path):
            if i[-1] == 'g':
                it.append(i)

        self.graphs_list.addItems(it)

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
