from install_requiremets_windows import Install


def run_app():
    """
    Run pyAGDAS
    :return:
    """
    from PyQt5 import QtWidgets
    from main import Main

    app = QtWidgets.QApplication([])
    main = Main()
    main.show()
    app.exec_()


try:
    run_app()

except ModuleNotFoundError:
    Install()
    run_app()
