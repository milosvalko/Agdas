from install_requiremets_windows import Install

def run_app():
    """
    Run pyAGDAS.
    """
    from PyQt5 import QtWidgets
    from main import Main

    app = QtWidgets.QApplication([])
    main = Main()
    main.show()
    app.exec_()

def try_run():
    """
    Try to run Agdas. If some library is not installed, run install_requiremets_windows and open Agdas again.
    """
    try:
        run_app()

    except ModuleNotFoundError:
        Install()
        run_app()

try_run()