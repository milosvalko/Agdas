import sys, glob, os, urllib.request, csv
from PyQt5 import uic, QtWidgets, QtCore
from PyQt5.QtWidgets import QLabel
from time import sleep, time
from warning import Warning
from classes import Fall, projectFile, rawFile, dropFile, estim, matr_db, res_final, Graph
from CONFIG import getFG5X, getFG5, matrDatabase, statistic, separator, headers, logo_picture, round_line_ind, \
    warning_window
import sqlite3 as sql
from datetime import datetime, timedelta
from time import time
import matplotlib.pyplot as plt
from PyQt5.QtGui import QIcon, QPixmap
from math import sin, cos, pi, sqrt, floor
from functions import allan, roundList, date_to_mjd, rssq, movingAverage
import numpy as np
from scipy.stats import t
import scipy.signal as sig

PATH, _ = uic.loadUiType('gui/compute.ui')


class Compute(QtWidgets.QDialog, PATH):

    def __init__(self, path, stationData, instrumentData, processingResults, gravityCorrections, header2, rawlines,
                 header1, projDirPath, setFile):
        super().__init__()
        self.setupUi(self)
        self.setWindowIcon(QIcon(logo_picture))

        # set values to widgets
        self.gravimeter_box.addItems(['FG5X', 'FG5'])
        # self.prescale.addItems(['100','800','1000'])
        self.multiplex.setText(processingResults['multiplex'])
        self.scaleFactor.setText(processingResults['scaleFactor'])
        self.preScale.setText(str(int(processingResults['scaleFactor']) * int(processingResults['multiplex'])))
        self.frminT.setText(processingResults['fringeStart'])
        self.frmaxT.setText(str(int(processingResults['processedFringes']) + int(processingResults['fringeStart'])))
        self.grad.setText(stationData['gradient'])
        self.setPrescale(10)
        self.set_gravimeter()
        self.fmodf.setText(str(self.gravimeter['fmodf']))
        self.lpar.setText(str(self.gravimeter['Lpar']))
        self.poleCorr_file.setText(gravityCorrections['polarMotion'])
        self.gravimeter_box.currentTextChanged.connect(self.set_gravimeter)

        # connect buttons with method
        self.run.clicked.connect(self.Run)
        self.allDrop.stateChanged.connect(self.numDrops)
        self.downloadPoleCorr.clicked.connect(self.downloadPole)
        self.statistics.setDisabled(True)
        self.numDrop.valueChanged.connect(self.DisStat)
        self.split_set.stateChanged.connect(self.disabledSplit)
        self.sets_choose.activated.connect(self.currentSet)
        # self.num_sets.valueChanged.connect(self.splitSet)

        # make class values
        self.path = path
        self.stationData = stationData
        self.instrumentData = instrumentData
        self.processingResults = processingResults
        self.gravityCorrections = gravityCorrections
        self.columns_rawfile = header2
        self.raw_lines = rawlines
        self.header1 = header1
        self.projDirPath = projDirPath
        self.setFile = setFile

        self.setDelimiter(',')

        self.ndrops = len(self.raw_lines)
        # self.run=True

        self.kalpha.setText(str(50))

        # self.Prescale()
        self.numDrops()
        self.drop()  # open and load drop file

        self.show()
        self.exec()

    def set_gravimeter(self):

        grav = self.gravimeter_box.currentText()

        if grav == 'FG5X':
            self.gravimeter = getFG5X(self.ps)

        if grav == 'FG5':
            self.gravimeter = getFG5(self.ps)

    def setPrescale(self, ps):
        """
        Setter of prescale
        :param ps:
        """
        self.ps = ps

    def setDelimiter(self, delimiter):
        """
        Setter of delimiter for printing text files
        :param delimiter: character, for example ;
        """
        self.delimiter = delimiter

    def currentSet(self):
        """
        Show count of drops in sets, when split on sets is giving by the user
        """
        s = self.numDrop.value() / int(self.sets_choose.currentText())
        self.sets_view.setText(str(int(s)))

    def disabledSplit(self):
        """
        Change QComboBox Enabled/Disabled by click on the checkbox
        """
        a = self.sets_choose.isEnabled()
        self.sets_choose.setDisabled(a)

    def DisStat(self):
        """
        This method set possible choice of sets to QComboBox
        """
        if self.numDrop.value() > 0:
            self.statistics.setDisabled(False)
        else:
            self.statistics.setDisabled(True)

        self.sets_choose.clear()
        x = self.numDrop.value()
        d = []
        for i in range(x):
            if x % (i + 1) == 0 and i + 1 != 1 and i + 1 != x:
                d.append(str(i + 1))

        self.sets_choose.addItems(d)

    def numDrops(self):
        """
        Set all drops for computing
        """
        check = self.allDrop.isChecked()
        self.numDrop.setRange(1, self.ndrops)

        if check == True:
            self.numDrop.setValue(self.ndrops)
            self.numDrop.setDisabled(True)

        if check == False:
            self.numDrop.setDisabled(False)
            self.numDrop.setValue(1)

    def downloadPole(self):
        '''
        Download pole coordinates from IERS and compute corrections for each drop
        '''

        url = 'https://datacenter.iers.org/data/csv/finals2000A.all.csv'
        try:
            urllib.request.urlretrieve(url, os.getcwd() + '/finals/finals2000A.all.csv')
        except urllib.error.URLError:
            Warning(error=warning_window['internet'], icon='critical', title='Warning')

        # open and load file from IERS
        file = open(os.getcwd() + '/finals/finals2000A.all.csv', 'r')
        # reader = csv.DictReader(file, delimiter=';')
        reader = csv.reader(file, delimiter=';')
        rows = list(reader)
        file.close()

        # date of first day
        a = datetime(int(self.lines[2].split()[4]), 1, 1) + timedelta(int(self.lines[2].split()[3]))

        month = (str(a.month))
        if len(month) == 1:
            month = '0' + month

        day = (str(a.day))
        if len(day) == 1:
            day = '0' + day

        date = [str(a.year), month, day]

        # get unique DOY
        doys = []
        years = []
        for l in self.lines:
            doy = int(l.split()[3])
            if doy not in doys:
                doys.append(doy)
                years.append(int(l.split()[4]))

        # get index of date in finals file
        i = 0
        for row in reversed(rows):
            if row[1:4] == date:
                today = row
                break
            i += 1

        # coordinates of pole from finals to interpolation
        x = []
        y = []
        for j in range(-i - 1, -i + len(doys)):
            x.append(rows[j][5])
            y.append(rows[j][7])

        # coordinates of station
        fi = float(self.stationData['lat']) * pi / 180
        lam = float(self.stationData['long']) * pi / 180

        # compute pole corrections
        self.dg = []
        for l in self.lines:
            line = l.split()

            doy = line[3]

            xInt = [float(x[doys.index(int(doy))]), float(x[doys.index(int(doy)) + 1])]
            yInt = [float(y[doys.index(int(doy))]), float(y[doys.index(int(doy)) + 1])]

            time = (line[2].split(':'))
            time = int(time[2]) / 3600 + int(time[1]) / 60 + int(time[0])

            ypole = time * (yInt[1] - yInt[0]) / 24 + yInt[0]
            xpole = time * (xInt[1] - xInt[0]) / 24 + xInt[0]

            self.dg.append(-19.139 * sin(2 * fi) * (xpole * cos(lam) - ypole * sin(lam)))
        # print(self.dg)
        self.poleCorrIERS.setText('<{}; {}>'.format(str(round(min(self.dg), 2)), str(round(max(self.dg), 2))))
        # print(self.dg)

    def defineSets(self):
        """
        Generate user define split to sets
        """

        drops = self.numDrop.value()
        set = int(self.sets_choose.currentText())
        # self.nset = set
        p = floor(drops / set)

        drop1 = [i + 1 for i in range(p)]

        self.drop_in_set = []
        self.sets = []
        for i in range(0, set):
            self.sets.extend([i + 1] * p)
            self.drop_in_set.extend(drop1)

        self.sets.extend([set + 1] * (drops % p))
        self.drop_in_set.extend([i + 1 for i in range(drops - set * p)])

    def drop(self):
        """
        Open and read drop file
        """
        for dropfile in glob.glob(self.path + '\*.drop.txt'):
            d = dropFile(dropfile)
            self.columns_dropfile = d.dropHeader4()
            self.lines = d.dropLines()

    def Run(self):
        """
        In this method is manage the whole calculating
        """
        # set color of RUN button on green
        self.run.setStyleSheet("background-color : green")
        # delete information about time of calculating from last computing
        self.calc_time.clear()

        # clear the logging window?
        if self.clwin.isChecked():
            self.logWindow.clear()

        # define user sets
        if self.split_set.isChecked():
            try:
                self.defineSets()
            except ValueError:
                Warning(error=warning_window['split_set'], icon='critical', title='Warning')
                self.run.setStyleSheet("background-color :#f0f0f0;")
                return

        # count of the fringes
        self.nfringe = int(self.header1[-1])

        # load kalpha
        kalpha = float(self.kalpha.toPlainText())

        frminss = self.gravimeter['frminss']
        self.frmaxss = self.gravimeter['frmaxss']

        # create estim file with head
        if self.files.isChecked():
            estim = res_final(path=self.projDirPath, header=headers['estim'].format(self.delimiter),
                              name=self.stationData['ProjName'] + '_' + 'estim', delimiter=self.delimiter)
            estim_grad = res_final(path=self.projDirPath, header=headers['estim_grad'].format(self.delimiter),
                                   name=self.stationData['ProjName'] + '_' + 'estimgrad', delimiter=self.delimiter)

        # create database for save every measuring
        self.matr_connection = matr_db(self.projDirPath + '/data.db')

        # measuring time of run computing
        self.t = time()
        # count of the drops
        self.ndrop = self.numDrop.value()

        if self.split_set.isChecked():
            self.nset = int(self.sets_choose.currentText())
        else:
            self.nset = int(self.processingResults['setsCollected'])

        files = False
        atm = []  # atmospheric correction data
        time_gr = []  # data for x axis for atm/baro/tides data
        baro = []  # barometric correction data
        tides = []  # tide correction data
        # self.all_xef=[]
        self.allRes = np.zeros((self.ndrop, self.frmaxss))  # numpy array for residuals of all drops
        self.v0 = []  # list of v0 values
        self.g0 = []  # list of g0 values
        self.resgradsum4 = np.zeros(
            (self.ndrop, int(self.processingResults['totalFringes'])))  # All residuals from gradient estimation fit
        self.ssresAr = []  # Standard deviations of fits with gradient
        self.m0grad4Sig = []  # Standard deviations of gradient estimation fit

        # loop for all drops
        for i in range(0, self.ndrop):

            # ===========================================================================#
            # create drop dictionary with data from dropfile
            drop_line = self.lines[i].split()
            drop = dict(zip(self.columns_dropfile, drop_line))
            # ===========================================================================#
            # load users define split to set
            if self.split_set.isChecked():
                drop['Set'] = self.sets[i]
                drop['Drp'] = self.drop_in_set[i]

            # ===========================================================================#
            # create raw dictionary with data from rawfile
            raw_line = self.raw_lines[i].split()
            d1 = dict(zip(self.columns_rawfile[0:5], raw_line[0:5]))
            d1['ftime'] = raw_line[5:5 + self.nfringe]
            d2 = dict(zip(self.columns_rawfile[6:], raw_line[5 + self.nfringe:10 + self.nfringe]))
            raw = d1 | d2
            # ===========================================================================#

            # ===========================================================================#
            # get wave length of laser
            laser = 'I{}'.format(drop['LaserLock'])
            Lambda = self.instrumentData[laser]
            # ===========================================================================#

            # ===========================================================================#
            # compute of LST
            fall = Fall()
            # Setters
            fall.setFringe(raw['ftime'])
            fall.setLambda(Lambda)
            fall.setScaleFactor(self.processingResults['scaleFactor'])
            fall.setMultiplex(self.processingResults['multiplex'])
            # fall.setGradient(self.stationData['gradient'])
            fall.setGradient(float(self.grad.toPlainText()))
            # fall.setModulFreq(self.instrumentData['modulFreq'])
            fall.setModulFreq(float(self.fmodf.toPlainText()))
            # fall.setLpar(self.FG5X['Lpar'])
            fall.setLpar(float(self.lpar.toPlainText()))
            fall.setRubiFreq(self.instrumentData['rubiFreq'])
            fall.setFrRange(int(self.frminT.toPlainText()), int(self.frmaxT.toPlainText()))
            # fall.setFrRange(frmin,frmax)
            fall.setFRssRange(self.frmaxss, frminss)
            fall.setKpar(self.kpar.isChecked())
            fall.setPcable(self.gravimeter['Pcable'])
            fall.setAcable(self.gravimeter['Acable'])
            fall.setLcable(self.gravimeter['Lcable'])
            if self.kdis.isChecked():
                fall.checkKDIS()
            if self.kimp.isChecked():
                fall.checkKIMP()
            if self.ksae.isChecked():
                fall.checkKSAE()
            # Compute fit by least squere method
            fall.LST()
            fall.effectiveHeight()
            fall.effectiveHeightTop()
            fall.effectivePosition()
            fall.gTop()
            fall.gTopCor(drop['Tide'], drop['Load'], drop['Baro'], drop['Polar'])
            self.tt = fall.tt
            self.Lambda = fall.Lambda
            # ===========================================================================#
            self.v0.append(fall.x_grad[0][1] * 1e-6)
            self.g0.append(fall.x_grad[0][2])

            # ===========================================================================#
            self.resgradsum4[i, :] = fall.resgrad4
            self.m0grad4Sig.append(fall.m0grad4)

            # ===========================================================================#
            if self.files.isChecked():
                # create line for estim file
                estim_line = self.estimLine(fall.x_grad[0], fall.std_grad, drop['Set'], drop['Drp'], fall.m02_grad)
                # print line to the estim file
                estim.printResult(line=roundList(estim_line, round_line_ind['estim']))

                # Create line of estimgrad file
                estim_grad_line = [drop['Set'], drop['Drp']]
                estim_grad_line.extend(fall.x_grad[0])
                if self.kpar.isChecked() == False:
                    estim_grad_line.extend(['-', '-'])
                estim_grad_line.extend(fall.x[0])
                if self.kpar.isChecked() == False:
                    estim_grad_line.extend(['-', '-'])
                estim_grad_line.extend(fall.xgrad4[0][:3])
                estim_grad_line.append(fall.stdGradX[2])

                estim_grad.printResult(line=roundList(estim_grad_line, round_line_ind['estim_grad']))

            date = datetime(int(drop['Year']), 1, 1) + timedelta(int(drop['DOY']) - 1)

            # decision if measuring is accepted by kalpha
            accepted = True
            self.ssresAr.append(fall.ssres)
            if fall.ssres > kalpha:
                accepted = False

            # polar correction from file
            if self.useFilePoleCorr.isChecked():
                Polar = float(self.poleCorr_file.toPlainText())

            # polar correction calculated for each drop from IERS file
            if self.useIERSPoleCorr.isChecked():
                Polar = self.dg[i]

            # residuals added into matrix with all residuals
            self.allRes[i, :] = fall.res_grad1[0:self.frmaxss]

            # transfer residuals to string for adding into database
            res = ', '.join(str(r) for r in fall.res_grad1[0:self.frmaxss])
            # date for database
            date_database = drop['Year'] + ':' + str(date.month).zfill(2) + ':' + str(date.day).zfill(2) + ' ' + drop[
                'Time']
            day_ = drop['Time'].split(':')
            day_ = (int(day_[-1]) / 3600 + int(day_[1]) / 60 + int(day_[0])) / 24 + date.day
            date_mjd = date_to_mjd(int(drop['Year']), date.month, day_)

            # line for database
            try:
                matr_drop = [i + 1, drop['Set'], drop['Drp'], date_database, date_mjd, fall.x_grad[0][0],
                             fall.x_grad[0][1], fall.x_grad[0][3], fall.x_grad[0][4], fall.x_grad[0][5],
                             fall.x_grad[0][6], fall.x_grad[0][7], fall.x_grad[0][8], fall.g0_Gr,
                             - fall.gradient * fall.Grad,
                             float(drop['Tide']) * 10, float(drop['Load']) * 10, float(drop['Baro']) * 10, Polar * 10,
                             fall.gTopCor, fall.g0,
                             fall.h * 1e-6, fall.Grad * 1e-6, fall.xgrad4[0][2], fall.m0gradient, fall.std,
                             fall.xef[0][3], fall.ssres, accepted, res]


            except UnboundLocalError:
                Warning(error=warning_window['pole_corr'], icon='critical', title='Warning')
                break

            except IndexError:
                matr_drop = [i + 1, drop['Set'], drop['Drp'], date_database, date_mjd, fall.x_grad[0][0],
                             fall.x_grad[0][1], fall.x_grad[0][3], fall.x_grad[0][4], fall.x_grad[0][5],
                             fall.x_grad[0][6], 0.0, 0.0, fall.g0_Gr, - fall.gradient * fall.Grad,
                             float(drop['Tide']) * 10, float(drop['Load']) * 10, float(drop['Baro']) * 10, Polar * 10,
                             fall.gTopCor, fall.g0,
                             fall.h * 1e-6, fall.Grad * 1e-6, fall.xgrad4[0][2], fall.m0gradient, fall.std,
                             fall.xef[0][3], fall.ssres, accepted, res]

            # send line to database
            self.matr_connection.insert(matrDatabase['insert'].format(*matr_drop))

            # send message to logging window
            mess = (
                'Drop: {} >> g: {:.2f} m0: {:.2f}'.format(str(i + 1).rjust(len(str(self.ndrop))), round(fall.g0_Gr, 2),
                                                          round(fall.m02_grad[0], 2)))
            self.logWindow.append(mess)
            self.progressBar.setValue(float((i + 1) * 100 / self.ndrop))

            # data for atmsorefic graphs
            atm.append(float(drop['Baro']))
            baro.append(float(drop['Pres']))
            tides.append(float(drop['Tide']))
            time1 = drop['Time'].split(':')
            time1 = int(time1[2]) / 3600 + int(time1[1]) / 60 + int(time1[0])
            if i > 0:
                if time1 < time_gr[-1]:
                    time1 += 24
            time_gr.append(time1)

            QtCore.QCoreApplication.processEvents()
            # QtCore.QCoreApplication.sendPostedEvents()

        # commit data to database
        self.matr_connection.commit()

        # close connection with database
        # self.matr_connection.close()

        # Compute statistics
        if self.statistics.isChecked():
            self.rejectBySigma()
            self.meanResidualsBySets()
            self.sensitivity()
            self.fourier()
            self.compute_normres()
            self.print_allanFile()
            self.ressets_res()
            # self.harmonic()

        # print results with gradient to estim file
        if self.files.isChecked():
            try:
                self.writeDropsFile()
            except AttributeError:
                Warning(error='Cannot write file due statistic is not computed', icon='critical', title='Warning')

            try:
                self.write_res_final()

            except AttributeError:
                Warning(error=warning_window['cannot_wrtite_file'], icon='critical', title='Warning')

        # Create graph
        if self.graph_save.isChecked():
            g = Graph(path=self.projDirPath + '/Graphs', name='atm_corr', project=self.stationData['ProjName'],
                      show=self.open_graphs.isChecked(), x_label='Time /h', y_label='Correction /μGal',
                      title='Atmosferic correction')
            g.plotXY(x=[time_gr], y=[atm], mark=['b+'], columns_name=['atm_corr'])
            g.saveSourceData()
            g.save()

            g = Graph(path=self.projDirPath + '/Graphs', name='atm_press', project=self.stationData['ProjName'],
                      show=self.open_graphs.isChecked(), x_label='Time /h', y_label='Recorder pressure /hPa',
                      title='Atmosferic pressure')
            g.plotXY(x=[time_gr], y=[baro], mark=['b+'], columns_name=['atm_press'])
            g.saveSourceData()
            g.save()

            g = Graph(path=self.projDirPath + '/Graphs', name='tides', project=self.stationData['ProjName'],
                      show=self.open_graphs.isChecked(), x_label='Time /h', y_label='Tides /μGal',
                      title='Tidal acceleration')
            g.plotXY(x=[time_gr], y=[tides], mark=['b+'], columns_name=['tides'])
            g.saveSourceData()
            g.save()

            self.Graph_EffHeight_CorToEffHeight(project=self.stationData['ProjName'])
            self.graphGravityChange()
            self.graphRes()
            self.graphSetG()
            self.graphParasitic()
            self.graphHistogramAccDrops()
            self.graphHistogramAccDropsNorm()
            self.graphEffectiveHeights2()
            self.graphSensitivityStd()
            self.graphVGG()
            self.graphResidualsBySets()
            self.allResGraph()
            self.graphResiduals()
            self.graphResidualsGradient()
            self.graphSpectrumParts()
            self.graphSpectrumRatio()

            r = self.matr_connection.get('select gTopCor from results where Accepted = 1')
            data = [i[0] for i in r]
            title = 'Drop data'
            ylabel = r'$-/nm.s^2$'
            name = 'allan1'
            self.graphAllan1(data, title, ylabel, name)

            title = 'Drop data - normalized'
            ylabel = r'$g-g_0  /nm.s^2$'
            name = 'allan1_normalized'
            self.graphAllan1(data=self.normres, title=title, ylabel=ylabel, name=name)

            r = self.matr_connection.get('select Gradient from results where Accepted = 1')
            data = [i[0] for i in r]
            title = 'Gradient'
            ylabel = r'Data - Median $-/nm.s^2/mm$'
            name = 'allan3'
            self.graphAllan1(data, title, ylabel, name)

        # close files
        if self.files.isChecked():
            estim.close()
            estim_grad.close()

        # Change color of Run button
        self.run.setStyleSheet("background-color:#f0f0f0;")
        # Set time of run calculation
        self.calc_time.setText('Calculation time: {:.2f} s'.format(time() - self.t))
        self.calc_time.setStyleSheet('color: red; font-size: 10pt')

    def ressets_res(self):

        prescale = int(self.processingResults['multiplex']) * int(self.processingResults['scaleFactor'])

        nset = self.matr_connection.get('select max(Set1) from results')
        nset = nset[0][0]

        tfit = self.meanResSets[1, :]

        self.v0m_bysets = []
        self.g0m_bysets = []
        self.tkor = []
        self.zzh = np.zeros((self.gravimeter['frmaxplot'], nset))
        for i in range(nset):

            a = self.matr_connection.get('select v0_withGR from results where Set1 = {}'.format(i + 1))
            self.v0m_bysets.append(np.median(a) * 1e-9)

            a = self.matr_connection.get('select g0_Gr from results where Set1 = {} '.format(i + 1))
            self.g0m_bysets.append(np.median(a) * 1e-9)

            self.tkor.append(-self.v0m_bysets[-1] / self.g0m_bysets[-1])

            self.zzh[0, i] = 0.5 * 9.809 * (tfit[0] - self.tkor[-1]) ** 2

            for j in range(1, self.gravimeter['frmaxplot']):
                self.zzh[j, i] = self.zzh[0, i] + self.Lambda * 1e-9 / 2 * prescale * j

    def allanGraph(self, a, tau, path):
        """
        Create loglog graph of allan standard deviations

        @param a: result of allan function
        @param tau: list with counts of intervals
        @param path: path for saving graphs
        """
        p = plt
        p.loglog(tau[:len(a)], [i[0] for i in a], '.r', ms=20)
        p.loglog(tau[:len(a)], [a[0][0] / np.sqrt(tau[i]) for i in range(len(a))], '-r')
        p.plot([tau[:len(a)], tau[:len(a)]], [[a[i][0] - a[i][0] / np.sqrt(a[i][2]) for i in range(len(a))],
                                              [a[i][0] + a[i][0] / np.sqrt(a[i][2]) for i in range(len(a))]], '-k',
               lw=2)
        p.plot([[i * 0.95 for i in tau[:len(a)]], [i * 1.05 for i in tau[:len(a)]]],
               [[a[i][0] + a[i][0] / np.sqrt(a[i][2]) for i in range(len(a))],
                [a[i][0] + a[i][0] / np.sqrt(a[i][2]) for i in range(len(a))]], '-k', lw=2)
        p.plot([[i * 0.95 for i in tau[:len(a)]], [i * 1.05 for i in tau[:len(a)]]],
               [[a[i][0] - a[i][0] / np.sqrt(a[i][2]) for i in range(len(a))],
                [a[i][0] - a[i][0] / np.sqrt(a[i][2]) for i in range(len(a))]], '-k', lw=2)
        p.legend(['Mean values', 'White noise'])
        p.title('Allan deviation - normalized data')
        p.xlabel('Drop number (n)')
        p.ylabel(r'σ (n)/nm.s$^-2$')
        p.savefig(path)
        plt.close()

    def graphAllan1(self, data, title, ylabel, name):
        """
        Create graph with drop data

        @param data: list of plotting data
        @param title: title of the graph
        @param ylabel:
        @param name: name of the graph
        """
        median = np.median(data)
        data_centered = [i - median for i in data]  # data centered on median
        median1 = np.median(np.abs(data_centered)) / 0.6745  # median from absolute value of residuals
        min = np.min(data_centered)
        max = np.max(data_centered)

        p = plt
        p.plot(data_centered, '.', ms=6)
        p.plot([0, len(data_centered)], [5 * median1, 5 * median1], 'r', lw=0.5)
        p.plot([0, len(data_centered)], [-5 * median1, -5 * median1], 'r', lw=0.5)
        p.title(title)
        p.xlabel('Drop #')
        p.ylabel(ylabel)

        path = self.projDirPath + '/Graphs/'
        project = self.stationData['ProjName']
        p.savefig(path + project + '_' + name + '.png', format='png', transparent=False)
        p.close()

    def graphSpectrumRatio(self):
        """
        Compute 2 subplot in one graph
        1) FFT of average and average of FFT
        2) their ratio
        """
        frmin = self.gravimeter['frmin']
        frmax = self.gravimeter['frmax']
        nforfft = self.gravimeter['nforfft']
        indexpad = self.matr_connection.get('select count(*) from results where Accepted = 1')[0]

        fs = nforfft / (2 * (self.tin[nforfft - 1] - self.tin[0]))
        frk = 2 * fs / (nforfft - 3)
        fr = np.arange(0, fs + 1, frk)

        n = 1
        # ratio of arrays
        ratio = [self.yfdMean[0, i] / self.yffa[i] for i in range(self.yffa.shape[0])]

        # ploting of graph
        p, (ax1, ax2) = plt.subplots(2, 1)

        ax1.loglog(fr[n:], self.yffa[n:], '-r', fr[n:], self.yfdMean[0, n:], '-b', lw=0.5)
        ax1.set(title='Spectrum comparison', xlabel='Frequency /Hz', ylabel='Amplitude /nm')
        ax1.legend(['FFT of average', 'Average of FFT'])

        ax2.loglog(fr, ratio, color=(0.64, 0, 1), lw=0.5)
        ax2.loglog([np.min(fr), fs], [np.sqrt(indexpad), np.sqrt(indexpad)], '-k', lw=0.5)
        ax2.set(xlabel='Frequency /Hz', ylabel='Ratio')
        ax2.legend(['Average of FFT/FFT of average', 'sqrt of #drops'])

        p.tight_layout()

        # save graph
        path = self.projDirPath + '/Graphs/'
        name = 'spectrum_ratio'
        project = self.stationData['ProjName']
        p.savefig(path + project + '_' + name + '.png', format='png', transparent=False)
        plt.close()

    def graphSpectrumParts(self):
        """

        """
        frmin = self.gravimeter['frmin']
        frmax = self.gravimeter['frmax']
        nforfft = self.gravimeter['nforfft']

        self.tin = np.linspace(self.tt[frmin - 1], self.tt[frmax - 1], nforfft)
        tin1 = np.linspace(self.tt[frmin - 1], self.tt[int(frmin + (frmax - frmin + self.ps) / 2 - 1)], nforfft)
        tin2 = np.linspace(self.tt[int(frmin + (frmax - frmin + self.ps) / 2)], self.tt[frmax - 1], nforfft)

        res1 = np.interp(tin1, self.tt[:self.gravimeter['frmaxplot']], self.meanRes[0, :])
        res2 = np.interp(tin2, self.tt[:self.gravimeter['frmaxplot']], self.meanRes[0, :])

        yres1 = 2 / nforfft * np.fft.fft(res1)
        yres1a = np.abs(yres1[:int((nforfft - 1) / 2)])

        yres2 = 2 / nforfft * np.fft.fft(res2)
        yres2a = np.abs(yres2[:int((nforfft - 1) / 2)])

        fs = nforfft / (2 * (self.tin[nforfft - 1] - self.tin[0]))

        fs1 = nforfft / (2 * (tin1[nforfft - 1] - tin1[0]))
        frk1 = 2 * fs1 / (nforfft - 3)
        fr1 = [i * frk1 for i in range(0, int((nforfft - 1) / 2))]

        fs2 = nforfft / (2 * (tin2[nforfft - 1] - tin2[0]))
        frk2 = 2 * fs2 / (nforfft - 3)
        fr2 = [i * frk2 for i in range(0, int((nforfft - 1) / 2))]

        frx = range(5, int(1e4) + 1)

        yres1x = np.interp(frx, fr1, yres1a)
        yres2x = np.interp(frx, fr2, yres2a)

        ratio = [yres1x[i] / yres2x[i] for i in range(yres1x.shape[0])]

        p, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 10))
        ax1.loglog(frx, yres1x, '-r', frx, yres2x, '-b', lw=0.5)
        ax1.set(title='Spectrum comparison for 2 parts of drops', xlabel='Frequency /Hz', ylabel='Amplitude /nm')
        ax1.legend(['1. half of drop', '2. half of drop'])

        ax2.loglog(frx, ratio, color=(0.64, 0, 1), lw=0.5)
        ax2.loglog([np.min(frx), fs], [1, 1], '-k', lw=0.5)
        ax2.set(xlabel='Frequency /Hz', ylabel='Ratio')
        ax2.legend(['1. half of drop / 2. half of drop'])

        p.tight_layout()

        path = self.projDirPath + '/Graphs/'
        name = 'spectrum_parts'
        project = self.stationData['ProjName']
        p.savefig(path + project + '_' + name + '.png', format='png', transparent=False)
        plt.close()

    def graphResidualsGradient(self):

        # self.resgradsum4Mean = np.loadtxt('resgradsum4x.csv', delimiter = ';')
        yl = 0.75
        n = 4

        aa, bb = sig.butter(n, self.kcutoff, btype='low')
        self.resgradsum4filt = sig.filtfilt(aa, bb, self.resgradsum4Mean)

        resm0 = [0, 0]
        tt0 = [self.tt[0], self.tt[self.nfringe - 1]]

        ylim = [-yl, yl]
        xlim = [self.tt[self.gravimeter['frmin']], self.tt[self.gravimeter['frmin']]]
        xxlim = [self.tt[self.gravimeter['frmax']], self.tt[self.gravimeter['frmax']]]

        import matplotlib as mpb
        mpb.use('ps')

        p = mpb.pyplot
        p.rcParams['figure.figsize'] = (25, 10)

        p.plot(tt0, resm0, '-k')
        p.plot(xlim, ylim, '-b')
        p.plot(xxlim, ylim, '-b')

        # p.plot(tt[:frmaxplot], resmm, '-k', lw = 2)
        p.plot(self.tt[:self.gravimeter['frmaxplot']], self.resgradsum4Mean[0, :self.gravimeter['frmaxplot']], '-k',
               lw=1)
        p.plot(self.tt[:self.gravimeter['frmaxplot']], self.resgradsum4filt[0, :self.gravimeter['frmaxplot']], '-',
               lw=3,
               color=(1, 0, 1))
        p.text(xlim[0] + 0.001, -yl, 'Start fringe', color='b')
        p.text(xxlim[0] + 0.001, -yl, 'Final fringe', color='b')

        p.title('Mean residuals - gravity gradient estimation')
        p.ylabel('Residuals /nm')
        p.xlabel('Time /s')

        path = self.projDirPath + '/Graphs/'
        name = 'residuals_gradient'
        project = self.stationData['ProjName']
        p.savefig(path + project + '_' + name + '.png', format='png', transparent=False)
        p.close()

        mpb.use('Qt5Agg')

    def graphResiduals(self):

        yl = 0.75  # ylimit
        resm0 = [0, 0]
        tt0 = [self.tt[0], self.tt[self.nfringe - 1]]

        ylim = [-yl, yl]
        xlim = [self.tt[self.gravimeter['frmin']], self.tt[self.gravimeter['frmin']]]
        xxlim = [self.tt[self.gravimeter['frmax']], self.tt[self.gravimeter['frmax']]]

        import matplotlib as mpb
        mpb.use('ps')

        siz = 15  # range for y axes

        # r = self.matr_connection.get('select Accepted from results')

        p = mpb.pyplot
        p.rcParams['figure.figsize'] = (25, 10)

        p.plot(tt0, resm0, '-k')
        p.plot(xlim, ylim, '-b')
        p.plot(xxlim, ylim, '-b')
        for i in self.meanResSets:
            p.plot(self.tt[:self.gravimeter['frmaxplot']], i, '-', lw=1, color='0.75')

        p.plot(self.tt[:self.gravimeter['frmaxplot']], self.meanRes[0, :], '-k', lw=2)
        p.plot(self.tinc, self.yn, '-', lw=3, color='tab:pink')
        p.text(xlim[0] + 0.001, -yl, 'Start fringe', color='b')
        p.text(xxlim[0] + 0.001, -yl, 'Final fringe', color='b')

        p.title('Stacked residuals')
        p.ylabel('Residuals /nm')
        p.xlabel('Time /s')

        path = self.projDirPath + '/Graphs/'
        name = 'residuals'
        project = self.stationData['ProjName']
        p.savefig(path + project + '_' + name + '.png', format='png', transparent=False)
        p.close()

        mpb.use('Qt5Agg')

    def allResGraph(self):
        import matplotlib as mpb
        mpb.use('ps')

        siz = 15  # range for y axes

        r = self.matr_connection.get('select Accepted from results')

        p = mpb.pyplot
        p.rcParams['figure.figsize'] = (20, 13)
        p.ylim([-siz, siz])
        p.plot([self.gravimeter['frmin'], self.gravimeter['frmin']], [-siz, siz], '-b', lw=1)
        p.plot([self.gravimeter['frmax'], self.gravimeter['frmax']], [-siz, siz], '-b', lw=1)
        p.title('Residuals for all drops')
        p.ylabel('Residuals /nm')
        p.xlabel('Fringe #')
        x = range(1, self.gravimeter['frmaxplot'] + 1)
        for i in range(len(r)):
            acc = r[i]

            if acc:
                # black color
                p.plot(x, self.allRes[i, :], '-k', lw=0.2)
            else:
                # grey color
                p.plot(x, self.allRes[i, :], '-', color="0.5", lw=0.2)

        path = self.projDirPath + '/Graphs/'
        name = 'resid_all'
        project = self.stationData['ProjName']
        p.savefig(path + project + '_' + name + '.png', format='png', transparent=False)
        p.close()

        mpb.use('Qt5Agg')

    def graphResidualsBySets(self):
        """
        Printing shifted residuals by sets to graph
        """

        frmaxplot = self.gravimeter['frmaxplot']  # range of data
        frmin = self.gravimeter['frmin']  # start fringe
        frmax = self.gravimeter['frmax']  # final fringe
        x = self.tt[:frmaxplot]  # data for x axis

        j = 1
        X = []
        Y = []
        mark = []
        col_name = []
        lw = []
        XX = []  # shifted lines
        YY = []
        markk = []
        lww = []
        text_x = []  # description of lines, Set 1-nset
        text_y = []
        text_color = []
        for l in self.meanResSets:
            X.append(x)
            y = [k + j for k in l]
            Y.append(y)
            mark.append('-k')
            col_name.append('Set{}'.format(j))
            lw.append(0.3)
            text_x.append(0.265)
            text_y.append(j)
            text_color.append('k')

            XX.append([0, x[-1]])
            YY.append([j, j])
            markk.append('-b')
            lww.append(0.1)

            j += 1

        # start fringe
        XX.append([x[frmin], x[frmin]])
        YY.append([0.5, self.nset + 0.5])
        markk.append('-b')
        lww.append(0.3)

        # final fringe
        XX.append([x[frmax], x[frmax]])
        YY.append([0.5, self.nset + 0.5])
        markk.append('-b')
        lww.append(0.3)

        g = Graph(path=self.projDirPath + '/Graphs', name='residuals_shifted', project=self.stationData['ProjName'],
                  show=self.open_graphs.isChecked(), x_label='Time /s', y_label='Shifted Residuals /nm',
                  title='Set residuals', winsize=(15, 10))
        g.plotXY(x=X, y=Y, mark=mark, columns_name=col_name, lw=lw)
        g.plotXY(x=XX, y=YY, mark=markk, columns_name=col_name, lw=lww)
        g.text(x=[x[frmin], x[frmax]], y=[0.3, 0.3], t=['Start fringe', 'Final fringe'], c=['b', 'b'])
        g.text(x=text_x, y=text_y, t=col_name, c=text_color)
        g.saveSourceData()
        g.save()

        del X, Y, XX, YY, x

    def graphVGG(self):
        """
        Create graph of gradients
        """

        res = self.matr_connection.get('select Gradient, GradientLSTm0 from results where Accepted = 1')
        grad = [i[0] for i in res]

        moving_average, moving_avg_x = movingAverage(grad, n=50)

        vggp3 = np.mean(grad)  # average value of gradients
        mggp3 = np.std(grad)  # standart deviation of gradients
        xlim = [0, self.ndrop]  # xrange
        ylim = [vggp3, vggp3]  # yrange for mean
        yylim = [vggp3 - 3 * mggp3, vggp3 - 3 * mggp3]  # yrange for -3sigma
        yyylim = [vggp3 + 3 * mggp3, vggp3 + 3 * mggp3]  # yrange for +3sigma

        m0 = []  # vector of m0 values
        x = []  # x range for plotting data
        cumulative_average = []
        for i in range(1, len(res) + 1):
            m0.append(res[i - 1][1])
            x.append(i)
            cumulative_average.append(sum(grad[:i]) / len(grad[:i]))

        g = Graph(path=self.projDirPath + '/Graphs', name='vgg', project=self.stationData['ProjName'],
                  show=self.open_graphs.isChecked(), x_label='Drop #', y_label=r'$VGG /nm.s^2/mm$',
                  title='Estimated VGGs', winsize=(13, 8))
        g.error_bar(x, grad, m0, 'r', ms=5, capsize=5)
        g.plotXY(x=[x, moving_avg_x, xlim, xlim, xlim], y=[cumulative_average, moving_average, ylim, yylim, yyylim],
                 mark=['k-', '-b', '-p', '-y', '-y'],
                 columns_name=['cumulative_mean', 'moving_average', 'mean', '-3sigma_range', '+3sigma_range'],
                 legend=['Cumulative average', 'Moving average', 'Average vgg-value', '+3σ range', '-3σ range'],
                 lw=[3, 3, 1, 1, 1])
        g.saveSourceData()
        g.save()

    def graphSensitivityStd(self):
        """
        Variability of set g-values on the choice of first and final fringe
        """
        sensa_tx = self.gravimeter['sensa_tx']
        sensa_tn = self.gravimeter['sensa_tn']
        sensa_bx = self.gravimeter['sensa_bx']
        sensa_bn = self.gravimeter['sensa_bn']

        # y values
        celk = sensa_tx - sensa_tn
        dglrms = rssq(self.dglc) / np.sqrt(celk + 1)

        celk = sensa_bx - sensa_bn
        dgrrms = rssq(self.dgrc) / np.sqrt(celk + 1)

        # xrange
        ts = np.linspace(1, self.nset, self.nset)

        # legend
        l1 = r'$RMS_S({}-{}) = {:.1f}  nm.s^2$'.format(sensa_tn, sensa_tx, np.mean(dglrms))
        l2 = r'$RMS_F({}-{}) = {:.1f}  nm.s^2$'.format(sensa_bn, sensa_bx, np.mean(dgrrms))

        g = Graph(path=self.projDirPath + '/Graphs', name='sensitivity_std', project=self.stationData['ProjName'],
                  show=self.open_graphs.isChecked(), x_label='Set #', y_label=r'Standart deviation $[nm.s^2]$',
                  title='Variability of set g-values on the choice of first and final fringe')
        g.plotXY(x=[ts, ts], y=[dglrms, dgrrms], mark=['k+-', 'r+-'], columns_name=['left', 'right'], legend=[l1, l2],
                 lw=[1, 1])
        g.saveSourceData()
        g.save()

    def graphGravityChange(self):
        Y = []
        X = []
        l = []
        cn = []
        m = []
        lw = []

        # x range
        tttt = np.linspace(self.gravimeter['sens_bn'], self.gravimeter['sens_bx'] - 1,
                           self.gravimeter['sens_bx'] - self.gravimeter['sens_bn'])

        # sensitivity data
        for i in range(len(self.dgr)):
            # g.plotXY(x=[tttt], y=[dgr[i,:]], mark=['C'+str((i)%10)+ '-'], columns_name=['Set ' + str(i+1)], legend =['Set ' + str(i+1)])
            X.append(tttt)
            Y.append(self.dgr[i, :])
            l.append('Set ' + str(i + 1))
            cn.append('Set ' + str(i + 1))
            m.append('C' + str((i) % 10) + '-')
            lw.append(0.3)

        X.append(tttt)
        Y.append(self.dgrm.T)
        l.append('Mean')
        cn.append('Mean')
        m.append('k-')
        lw.append(1)

        g = Graph(path=self.projDirPath + '/Graphs', name='sensitivity_bottom', project=self.stationData['ProjName'],
                  show=self.open_graphs.isChecked(), x_label='Final Fringe #', y_label='',
                  title='Gravity change due to choice of the last fringe')
        g.plotXY(x=[tttt], y=[[0 for i in range(len(tttt))]], mark=['b-'], columns_name='xx', legend='', lw=[0.3])
        g.plotXY(x=[[self.gravimeter['frmax'], self.gravimeter['frmax']]], y=[[-10, 10]], mark=['b-'],
                 columns_name='xx', legend='',
                 lw=[0.3])
        g.plotXY(x=X, y=Y, mark=m, columns_name=cn, legend=l, lw=lw)
        g.saveSourceData()
        g.save()

        del X
        del Y

    def graphSetG(self):

        g0 = 1000 * (floor(self.gfinal / 1000))
        x = [0, self.nset + 1]

        res = self.matr_connection.get(statistic['mean:vxv'])
        mean_by_set = [i[2] - g0 for i in res]

        # Set gravity at top of the drop
        g = Graph(path=self.projDirPath + '/Graphs', name='set_g', project=self.stationData['ProjName'],
                  show=self.open_graphs.isChecked(), x_label='Set #',
                  y_label='Set gravity -{:.2f} /nm.s^(-2)'.format(g0), title='Set gravity at top of the drop')
        g.plotXY(x=[x, x, x],
                 y=[[self.gfinal - g0, self.gfinal - g0], [self.gfinal - g0 - self.gstd, self.gfinal - g0 - self.gstd],
                    [self.gfinal - g0 + self.gstd, self.gfinal - g0 + self.gstd]], mark=['b-', 'g-', 'g-'],
                 columns_name=['mean', 'mean-1σ', 'mean+1σ'], legend=['Set g-values', 'Avegare g-value', '1σ range'])
        g.error_bar(range(1, x[1]), mean_by_set, self.stodchmod, 'r')
        g.saveSourceData()
        g.save()

        # Standart deviation for set g-values
        g = Graph(path=self.projDirPath + '/Graphs', name='set_std', project=self.stationData['ProjName'],
                  show=self.open_graphs.isChecked(), x_label='Set #', y_label='Set standart deviation /nm.s^(-2)',
                  title='Standart deviation for set g-values')
        # g.plotXY(x=[x, x, x], y=[[gfinal-g0, gfinal-g0], [gfinal-g0-gstd, gfinal-g0-gstd], [gfinal-g0+gstd, gfinal-g0+gstd]], mark=['b-', 'g-', 'g-'], columns_name=['Sine component', 'Cosine component'], legend =['Set g-values', 'Avegare g-value', '1 range'])
        g.error_bar(range(1, x[1]), self.stodch, self.stodchs, 'r')
        g.saveSourceData()
        g.save()

    def graphRes(self):

        rms = self.matr_connection.get('select n, ssres from results')
        std = [r[1] for r in rms]
        n = range(1, len(std) + 1)

        g = Graph(path=self.projDirPath + '/Graphs', name='resid_RMS', project=self.stationData['ProjName'],
                  show=self.open_graphs.isChecked(), x_label='Drop #', y_label='Standart deviation /nm',
                  title='Standart deviation of the residuals for each drop')
        g.plotXY(x=[n], y=[std], mark=['-g'], columns_name=['rms'], legend=[])
        g.saveSourceData()
        g.save()

        # acc = self.matr_connection.get('select Accepted from results')
        # g=Graph(path=self.projDirPath+'/Graphs', name='resid_all', project = self.stationData['ProjName'], show=self.open_graphs.isChecked(), x_label='Fringe #', y_label = 'Residual /nm', title='Residuals for each drop')
        # x=range(1, self.frmaxss + 1)
        # for i in range(len(acc)):
        #     if acc[i][0] == 1:
        #         # g.plotXY(x=[x], y=[self.allRes[i, :]], mark=['-k'], columns_name=[], legend =[])
        #         g.gr.plot(x, self.allRes[i, :], '-k')
        #
        # # g.saveSourceData()
        # g.save()

    def graphParasitic(self):

        res = self.matr_connection.get('select e_withGR, f_withGR from results where Accepted = 1')

        e = [i[0] for i in res]
        f = [i[1] for i in res]
        x = range(1, len(e) + 1)

        g = Graph(path=self.projDirPath + '/Graphs', name='parasitic', project=self.stationData['ProjName'],
                  show=self.open_graphs.isChecked(), x_label='Drop #', y_label='sin/cos amplitude /nm',
                  title='Sine/Cosine amplitudes of the parasitic wave with L = {:.3f} m'.format(
                      float(self.lpar.toPlainText()) / 1e9))
        g.plotXY(x=[x, x], y=[e, f], mark=['r-', 'g-'], columns_name=['Sine component', 'Cosine component'],
                 legend=['Sine component', 'Cosine component'])
        g.saveSourceData()
        g.save()

    def graphHistogramAccDrops(self):

        r = self.matr_connection.get('select gTopCor from results where Accepted = 1')
        r = [i[0] for i in r]

        g = Graph(path=self.projDirPath + '/Graphs', name='histogram', project=self.stationData['ProjName'],
                  x_label='Drop gravity - final g/nm.s^{-2}', y_label='Frequency', title='Histogram of accepted drops',
                  show=self.open_graphs.isChecked())
        g.histogram(r, fit=True)
        g.saveSourceData()
        g.save()

    def graphHistogramAccDropsNorm(self):
        g = Graph(path=self.projDirPath + '/Graphs', name='histogram_norm', project=self.stationData['ProjName'],
                  x_label='Drop gravity - final g/normalized nm.s^{-2}', y_label='Frequency',
                  title='Histogram of accepted drops (normalized)', show=self.open_graphs.isChecked())
        g.histogram(self.normres, fit=True)
        g.saveSourceData()
        g.save()

    def graphEffectiveHeights2(self):

        y = self.matr_connection.get('select EffHeight + CorToEffHeight from results')
        y = [i[0] for i in y]
        x = [i for i in range(1, self.ndrop + 1)]
        g = Graph(path=self.projDirPath + '/Graphs', name='effective_height2', project=self.stationData['ProjName'],
                  x_label='Drop #', y_label='Effective measurement height /mm',
                  title='Effective measurement height from top of the drop', show=self.open_graphs.isChecked())
        g.plotXY(x=[x], y=[y], mark=['-b'], columns_name=['effective_height'])
        g.save()

    def estimLine(self, X, std, set, drop, m0):
        """
        Print result of drop to estim file
        """
        dropResult = [set, drop, m0[0]]
        for i in range(len(X)):
            dropResult.append(X[i])
            dropResult.append(std[i])

        dropResult[3] /= 1e6
        dropResult[4] /= 1e6
        dropResult[5] /= 1e6
        dropResult[6] /= 1e6

        if self.kpar.isChecked() == False:
            dropResult.extend(['-', '-', '-', '-'])

        return dropResult

    def writeDropsFile(self):

        r = self.matr_connection.get(
            'SELECT Set1, Drop1, Date, g0_Gr, CorrToTop, Tide, Load, Baro, Polar, gTopCor, g0, EffHeight, CorToEffHeight, Accepted from results')
        a = res_final(path=self.projDirPath, header=headers['drops'].format(self.delimiter),
                      name=self.stationData['ProjName'] + '_' + 'drops', delimiter=self.delimiter)

        for i in r:
            k = list(i[:4])

            k.append(self.stdodchpadu[i[0] - 1])

            for j in range(4, len(i)):
                k.append(i[j])

            k = roundList(k, round_line_ind['drops'])
            a.printResult(line=k)

        a.close()

    def printMatlog(self):
        r = self.matr_connection.get(matrDatabase['matlog'])

        a = open(self.setFile, 'r')

        a = a.read().splitlines()
        press = [a[i].split()[18] for i in range(4, len(a))]

        self.vgg_median_bysets = []

        a = res_final(path=self.projDirPath, header=headers['matlog'].format(self.delimiter),
                      name=self.stationData['ProjName'] + '_' + 'matlogsets', delimiter=self.delimiter)
        it = 0
        for i in r:
            tst = abs((t.cdf(self.vv[it] / self.mm, len(r) - 1) - 0.5) * 200)

            vgg = np.median(self.matr_connection.get('select vgg from results where Set1 = {}'.format(i[0])))

            self.vgg_median_bysets.append(vgg)

            line = [self.stationData['ProjName'], i[0], int(i[1]), int(i[2]), int(i[3]), int(i[4]), int(i[5]),
                    int(i[6]), i[-1], self.stationData['gradient'], i[7], self.stodch[it], self.dglrms[it],
                    self.dgrrms[it], i[8], float(self.stationData['actualHeight']) / 100, press[it], vgg, tst]
            a.printResult(line=roundList(line, round_line_ind['matlogsets']))
            it += 1

        a.close()

    def print_allanFile(self):

        file = res_final(path=self.projDirPath, header=headers['allan'].format(self.delimiter),
                         name=self.stationData['ProjName'] + '_' + 'allan', delimiter=self.delimiter)
        # file.write('n;ALLAN1;STD1;ALLAN2;STD2;ALLAN3;STD3\n')
        tau = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800,
               900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500]

        r = self.matr_connection.get('select gTopCor, Gradient from results where Accepted = 1')
        gTopCor = [i[0] for i in r]
        grad = [i[1] for i in r]

        # print(self.normres)
        a1 = allan(gTopCor, tau)
        # np.savetxt('gTopCor.csv', gTopCor, delimiter = ';')
        a2 = allan(self.normres, tau)
        a3 = allan(grad, tau)

        for i in range(len(a1)):
            line = [tau[i], a1[i][0], a1[i][1], a2[i][0], a2[i][1], a3[i][0], a3[i][1]]
            file.printResult(line=roundList(line, round_line_ind['allan']))

        file.close()

        if self.graph_save.isChecked():
            self.allanGraph(a2, tau, self.projDirPath + '/Graphs/' + self.stationData['ProjName'] + '_allan_deviation')

            self.allanGraph(a3, tau, self.projDirPath + '/Graphs/' + self.stationData['ProjName'] + '_allan_gradient')

    def compute_normres(self):
        self.logWindow.append(separator)
        self.logWindow.append('Compute allan standart deviation')
        QtCore.QCoreApplication.processEvents()

        r = self.matr_connection.get('select max(Set1) from results')
        g = self.matr_connection.get('select avg(gTopCor) from results where Accepted = 1 group by Set1')

        # r=[i[0] for i in r]
        r = r[0][0]
        nset = (r)
        # self.nset = nset

        ksmooth = self.gravimeter['ksmooth']
        self.stodch = []  # mean error of sets
        self.stdodchpadu = []  # mean error of drops
        count = 0
        weight = []  # weight of sets
        gfinal = 0  # weight mean of g by sets
        sumweight = 0
        self.vv = []  # mean correction
        tst = []  # what is this
        self.normres = []
        self.stodchmod = []
        self.stodchs = []
        for i in range(r):
            d = self.matr_connection.get('select gTopCor from results where Accepted = 1 and Set1 = {}'.format(i + 1))

            count += len(d)
            self.stodch.append(np.std(d) / np.sqrt(len(d)))
            self.stodchs.append(self.stodch[-1] / np.sqrt(2 * len(d)))
            stodchmod = self.stodch[-1] * self.stodch[-1] + ksmooth * ksmooth
            self.stodchmod.append(sqrt(stodchmod))
            self.stdodchpadu.append(self.stodch[-1] * np.sqrt(len(d)))
            # print('STD {}'.format(stdodchpadu[-1]))
            weight.append(100 / stodchmod)
            sumweight += weight[-1]
            gfinal += g[i][0] * weight[-1]

        self.gfinal = gfinal / sumweight

        # if nset>1:
        self.mm = 0
        for i in range(r):
            self.vv.append((g[i][0] - self.gfinal) * np.sqrt(weight[i]))
            self.mm += self.vv[-1] * self.vv[-1]
            # print(weight[i])

        self.mm = np.sqrt(self.mm / (nset - 1))

        # if count>=1:
        gstd = np.std(self.vv) / np.sqrt(sumweight)
        self.gstd = gstd

        gtop = self.matr_connection.get('select gTopCor, Set1 from results where Accepted = 1')

        for i in gtop:
            # print((i[0]-self.gfinal)/stdodchpadu[i[1]-1]*gstd*np.sqrt(count))
            self.normres.append((i[0] - self.gfinal) / self.stdodchpadu[i[1] - 1] * gstd * np.sqrt(count))

        if self.files.isChecked():
            self.print_allanFile()
            self.printMatlog()

    def write_res_final(self):
        prescale = int(self.processingResults['multiplex']) * int(self.processingResults['scaleFactor'])
        # it=0

        # self.FG5X['frmaxplot']
        self.v0m = np.median(self.v0) / 1e3
        self.g0m = np.median(self.g0) / 10e9
        v0mg0mkor = (-self.v0m / self.g0m) / 10

        self.tinc = np.linspace(self.tt[0], self.tt[self.gravimeter['frmaxplot']], self.gravimeter['nforfft'])

        a = res_final(path=self.projDirPath, header=headers['residuals_final'].format(self.delimiter),
                      name=self.stationData['ProjName'] + '_' + 'residuals_final', delimiter=self.delimiter)
        by_sets = res_final(path=self.projDirPath, header=headers['residuals_sets'].format(self.delimiter),
                            name=self.stationData['ProjName'] + '_' + 'residuals_sets', delimiter=self.delimiter)
        resgradsum = res_final(path=self.projDirPath, header=headers['resgradsum'].format(self.delimiter),
                               name=self.stationData['ProjName'] + '_' + 'resgradsum', delimiter=self.delimiter)

        a1000 = res_final(path=self.projDirPath, header=headers['residuals_final1000'].format(self.delimiter),
                          name=self.stationData['ProjName'] + '_' + 'residuals_final1000', delimiter=self.delimiter)

        # =======================================================================
        resyyy = np.interp(self.tinc, self.tt[:self.gravimeter['frmaxplot']], self.meanRes[0, :])

        coff = int(self.coff.toPlainText())
        n = 4
        self.kcutoff = 2 * coff * (self.tt[self.gravimeter['frmaxplot'] - 1] - self.tt[0]) / self.gravimeter['nforfft']

        aa, bb = sig.butter(n, self.kcutoff, btype='low')
        self.yn = sig.filtfilt(aa, bb, resyyy)

        resmmi = np.interp(self.tt[:self.gravimeter['frmaxplot']], self.tinc, self.yn)

        # =======================================================================

        self.resgradsm4filt = sig.filtfilt(aa, bb, self.resgradsum4Mean)

        # =======================================================================

        # data for round value to print into file
        round_ind_bysets = [[1, 5], [2, 5], [3, 5]]
        round_ind_bysets.extend([[i, 6] for i in range(4, len(self.meanResSets[:, 0]) + 1)])
        for it in range(self.gravimeter['frmaxplot']):

            z = (it) * self.Lambda / 2 * 1e-9 * prescale
            line = [it + 1, z, self.tt[it], self.tt[it] - v0mg0mkor, self.meanRes[0, it], resmmi[it]]
            a.printResult(line=roundList(line, round_line_ind['residuals_final']))

            line = [it + 1, z, self.tt[it], self.tt[it] - v0mg0mkor]
            line.extend(self.meanResSets[:, it])
            by_sets.printResult(line=roundList(line, round_ind_bysets))

            line = [it + 1, z, self.tt[it], self.tt[it] - v0mg0mkor, self.resgradsum4Mean[0, it],
                    self.resgradsm4filt[0, it]]
            resgradsum.printResult(roundList(line, round_line_ind['resgradsum']))

            if (it + 1) % 10 == 5:
                line = [it + 1, z, self.tt[it], self.tt[it] - v0mg0mkor,
                        np.mean(self.resgradsum4Mean[0, it - 4:it + 6]), resmmi[it]]
                a1000.printResult(line=roundList(line, round_line_ind['residuals_final1000']))

        a.close()
        by_sets.close()
        resgradsum.close()
        a1000.close()

    def rejectBySigma(self):
        """
        Reject drops by sigma
        Reject drops by median
        """
        # Print to logWindow
        self.logWindow.append(separator)
        self.logWindow.append('Reject drops with rejsigma>3*std')
        QtCore.QCoreApplication.processEvents()

        # Get mean and v*v by sets from database
        mean = self.matr_connection.get(statistic['mean:vxv'])
        # Get gTopCor from databse
        res = self.matr_connection.get('select Set1, Drop1, gTopCor, Accepted from results where Accepted = 1')

        n = int(self.processingResults['dropsInSet'])
        # Get rejsigma from GUI
        rejsigma = float(self.rejsigma.toPlainText())

        # Median of
        resMed = np.median(self.ssresAr)

        m0grad4Med = np.median(self.m0grad4Sig)
        self.resgradsum4Mean = np.zeros((1, int(self.processingResults['totalFringes'])))

        kalpha = float(self.kalpha.toPlainText())
        kalpha_2 = float(self.kalpha_2.toPlainText())

        std = []
        mean1 = []
        for i in mean:
            std.append(sqrt(i[3] / (n - 1)))
            mean1.append(i[2])

        it = 0
        grad4Acc = 0
        for j in res:
            set_std = std[j[0] - 1]
            set_mean = mean1[j[0] - 1]

            # accepted if m0 < 1.5*median m0
            if self.m0grad4Sig[it] < (1 + kalpha_2 / 100) * m0grad4Med:
                self.resgradsum4Mean[0, :] += self.resgradsum4[it, :]
                grad4Acc += 1

            # accepted if m0 < 1.5*median m0
            if self.ssresAr[it] > (1 + kalpha / 100) * resMed:
                update = matrDatabase['updateAcc'].format(j[0], j[1])
                self.matr_connection.insert(update)

            # accepted if v < sigma*std
            if abs(j[2] - set_mean) > rejsigma * set_std:
                update = matrDatabase['updateAcc'].format(j[0], j[1])
                self.matr_connection.insert(update)

                # print(j)
            it += 1

        self.resgradsum4Mean = self.resgradsum4Mean / grad4Acc

        self.matr_connection.commit()

    def meanResidualsBySets(self):
        """
        Compute mean residuals by sets.
        """
        self.logWindow.append(separator)
        self.logWindow.append('Compute mean residuals by sets')
        QtCore.QCoreApplication.processEvents()

        # initialization of arrays for averages
        self.meanResSets = np.zeros((self.nset, self.frmaxss))
        self.meanRes = np.zeros((1, self.frmaxss))
        # get is drop accepted set values from database
        d = self.matr_connection.get('select Accepted, Set1 from results')
        self.d = d
        c = self.matr_connection.get('''select count(*) from results
        where Accepted = 1
        group by Set1''')
        self.count = c

        self.allAcc = 0
        for i in c:
            self.allAcc += i[0]

        it = 0
        for i in d:
            if i[0] == 1:
                # mean the residuals
                self.meanResSets[i[1] - 1, :] += self.allRes[it, :] / c[i[1] - 1][0]
                self.meanRes[0, :] += self.allRes[it, :]

            it += 1

        self.meanRes = self.meanRes / self.allAcc

    def fourier(self):
        """
        Compute fourier transformation for residuals: all res., mean res. by set and mean res.
        """
        self.logWindow.append(separator)
        self.logWindow.append('Fourier transformation')
        QtCore.QCoreApplication.processEvents()

        # split tt on nforfft parts
        tin = np.linspace(self.tt[self.gravimeter['frmin'] - 1], self.tt[self.gravimeter['frmax'] - 1],
                          self.gravimeter['nforfft'])
        ttx = self.tt[self.gravimeter['frmin'] - 1:self.gravimeter['frmax']]

        x = int((self.gravimeter['nforfft'] - 1) / 2)
        # arrays with results
        yfd = np.zeros((len(self.allRes), self.gravimeter['nforfft']), dtype=complex)
        yfdMeanBySet = np.zeros((self.nset, x))
        self.yfdMean = np.zeros((1, x))
        yfs = np.zeros((self.nset, self.gravimeter['nforfft']), dtype=complex)
        yfsa = np.zeros((self.nset, x))  # by set

        it = 0
        # fourier transformation for all drops
        for ress in self.allRes:
            # for accepted drops
            if self.d[it][0] == 1:
                ress = ress[self.gravimeter['frmin'] - 1: self.gravimeter['frmax']]

                # transformation of ress = f(ttx) => tin
                resd = np.interp(tin, ttx, ress)
                resd = resd - np.mean(resd)

                # fft
                fft = 2 / self.gravimeter['nforfft'] * np.fft.fft(resd)

                # fft for all residuals
                yfd[it, :] = fft

                set = self.d[it][1]

                l = np.absolute(fft[0:x]) / self.count[set - 1][0]
                yfdMeanBySet[set - 1, :] += np.real(l)

                l = np.absolute(fft[0:x] / self.allAcc)
                self.yfdMean[0, :] += np.real(l)

            it += 1

        # fourier transformation for mean residuals by set
        for i in range(self.nset):
            ress = np.interp(tin, ttx, self.meanResSets[i, self.gravimeter['frmin'] - 1:self.gravimeter['frmax']])

            ressm = ress - np.mean(ress)

            fft = 2 / self.gravimeter['nforfft'] * np.fft.fft(ressm)

            yfs[i, :] = fft

            yfsa[i, :] = np.real(np.absolute(fft[0:x]))

        # spectrum for mean all residuals
        resf = np.interp(tin, ttx, self.meanRes[0, self.gravimeter['frmin'] - 1:self.gravimeter['frmax']])
        resfm = resf - np.mean(resf)
        yff = 2 / self.gravimeter['nforfft'] * np.fft.fft(resfm)
        self.yffa = np.real(np.absolute(yff[0:x]))

        # writing data to text files
        if self.files.isChecked():
            tins = np.linspace(0, 10000, 2151)
            fs = self.gravimeter['nforfft'] / (2 * (tin[-1] - tin[0]))
            frk = 2 * fs / (self.gravimeter['nforfft'] - 3)
            # k=0
            # fr=[0]
            # for i in yffa:
            #     fr.append(fr[-1]+frk)
            fr = [i * frk for i in range(len(self.yffa))]
            yffas = np.interp(tins, fr, self.yffa)
            self.yfdamm = np.interp(tins, fr, self.yfdMean[0, :])

            # np.savetxt()
            a = res_final(path=self.projDirPath, header=headers['spectrum'].format(self.delimiter),
                          name=self.stationData['ProjName'] + '_' + 'spectrum', delimiter=self.delimiter)
            for i in range(len(tins)):
                line = [tins[i], yffas[i], self.yfdamm[i]]
                a.printResult(line=roundList(line, round_line_ind['spectrum']))

            a.close()

    def sensitivity(self):
        """
        Computing of sensitivity
        """
        self.logWindow.append(separator)
        self.logWindow.append('Compute sensitivity')
        QtCore.QCoreApplication.processEvents()

        sens_tn = self.gravimeter['sens_tn']
        sens_tx = self.gravimeter['sens_tx']

        # initialization of arrays for sensitivity
        dgl = np.zeros((self.nset, sens_tx - sens_tn + 1))
        self.dgr = np.zeros((self.nset, sens_tx - sens_tn + 1))

        dglm = np.zeros((1, sens_tx - sens_tn + 1))
        self.dgrm = np.zeros((1, sens_tx - sens_tn + 1))

        # sensitivity is calculated for averages by sets
        for i in range(self.nset):
            for j in range(sens_tn, sens_tx + 1):
                # data for fitting by parabola on left side of drop
                x = self.tt[j - 1: self.gravimeter['frmax']]
                y = self.meanResSets[i, j - 1:self.gravimeter['frmax']]
                # fitting by parabola
                koef = np.polyfit(x, y, deg=2)
                # storing quadratic coefficient of equation of fitted parabola
                dgl[i, j - 1] = koef[0] * 2

                # data for fitting by parabola on right side of drop
                x = self.tt[self.gravimeter['frmin'] - 1: j + self.gravimeter['sens_bn'] - 1]
                y = self.meanResSets[i, self.gravimeter['frmin'] - 1: j + self.gravimeter['sens_bn'] - 1]
                # fitting by parabola
                koef = np.polyfit(x, y, deg=2)
                # storing quadratic coefficient of equation of fitted parabola
                self.dgr[i, j - 1] = koef[0] * 2

            dglm += dgl[i, :]
            self.dgrm += self.dgr[i, :]

        # average
        dglm /= self.nset
        self.dgrm /= self.nset
        # ==============================================================================

        rozd = self.gravimeter['sensa_tn'] - self.gravimeter['sens_tn']
        celk = self.gravimeter['sensa_tx'] - self.gravimeter['sensa_tn']

        self.dglc = dgl[:, rozd:rozd + 1 + celk]

        self.dglrms = np.sqrt(np.sum(np.square(self.dglc.transpose()), axis=0)) / np.sqrt(celk + 1)

        # ==============================================================================
        rozd = self.gravimeter['sensa_bn'] - self.gravimeter['sens_bn']
        celk = self.gravimeter['sensa_bx'] - self.gravimeter['sensa_bn']

        self.dgrc = self.dgr[:, rozd:rozd + 1 + celk]

        self.dgrrms = np.sqrt(np.sum(np.square(self.dgrc.transpose()), axis=0)) / np.sqrt(celk + 1)
        # ==============================================================================

    def Graph_EffHeight_CorToEffHeight(self, project):

        # Get data from database
        res1 = self.matr_connection.get('select EffHeight, CorToEffHeight from results')
        # Open result file
        # file=open(self.projDirPath+'/Graphs/'+project+'_'+'effective_height_corr.csv', 'w')
        file = res_final(path=self.projDirPath, header=headers['effective_height_corr'].format(self.delimiter),
                         name=self.stationData['ProjName'] + '_' + 'effective_height_corr', files='/Graphs/',
                         delimiter=self.delimiter)

        # Create lists print graph and print lines of result file
        res = []
        res2 = []
        x = range(0, len(res1))
        for i in x:
            res.append(res1[i][0])
            res2.append(res1[i][1])

            line = [i + 1, res[i], res2[i]]
            file.printResult(line=roundList(line, round_line_ind['effHeightCorr_Graph']))

        file.close()

        # Print graph
        t = x
        data1 = res
        data2 = res2
        fig, ax1 = plt.subplots()

        color = 'tab:blue'
        ax1.set_xlabel('Drop')
        ax1.set_ylabel('Effective measurement height /mm', color=color)
        ax1.plot(t, data1, color=color, linewidth=0.5)
        ax1.tick_params(axis='y', labelcolor=color)

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        color = 'tab:red'
        ax2.set_ylabel('Top of the drop /mm', color=color)  # we already handled the x-label with ax1
        ax2.plot(t, data2, color=color, linewidth=0.5)
        ax2.tick_params(axis='y', labelcolor=color)

        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        fig.savefig(self.projDirPath + '/Graphs/' + project + '_' + 'effective_height.png', dpi=250)
        plt.close(fig)

    def end(self):
        self.logWindow.append(separator)
        self.logWindow.append('Done')
        # self.logWindow.append('Time of run computing: {} s'.format(round(time()-t)))
        QtCore.QCoreApplication.processEvents()


if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    compute = Compute()
    compute.show()
    app.exec_()

# def harmonic(self):
#     """
#     Effect of harmonica on VGG a change g
#     """
#     Lmin = self.FG5X['Lmin']
#     Lmax = self.FG5X['Lmax']
#     frmin=self.FG5X['frmin'] #start fringe
#     frmax=self.FG5X['frmax'] #final fringe
#     frmaxplot = self.FG5X['frmaxplot'] #range of data
#     hefm = self.matr_connection.get('select avg(EffHeight + CorToEffHeight) from results where Accepted = 1')
#     hefm = hefm[0][0]
#     gamma = float(self.stationData['gradient']) #Gradient from file
#     nset = self.nset
#     # self.v0m=np.median(self.v0)/1e3
#     # self.g0m=np.median(self.g0)/10e9
#
#     # tfit = np.loadtxt('tfit.csv', delimiter = ';')
#     tfit = self.meanResSets[1,:]
#
#     ressets2 = np.zeros((frmaxplot, nset))
#     A0 = np.ndarray((frmax-frmin+1, 3))
#     AH = np.ndarray((frmax-frmin+1, 3))
#     AL1 = np.ndarray((frmax-frmin+1, 3))
#     LA = np.ndarray((frmax-frmin+1, 3))
#     harmonic_res = np.ndarray((10*Lmax-10*Lmin+1, nset))
#     T0_m = np.ndarray((10*Lmax-10*Lmin+1, nset))
#     A0_m = np.ndarray((10*Lmax-10*Lmin+1, nset))
#     B0_m = np.ndarray((10*Lmax-10*Lmin+1, nset))
#
#     xg=[]
#     xvgg=[]
#     mx0=[]
#     mxL1=[]
#     dg=[]
#     mdg=[]
#     vgg=[]
#     mvgg=[]
#     # minimum_index_array = []
#     # minimum_array=[]
#     LLmin=[] #wave lengths m
#     T0min=[] #zeros member
#     A0min=[] #cosinus member
#     B0min=[] #sinus member
#     harm_amp=[] #amplitude of harmonic
#     harm_phase=[] #phase of harmonic
#     t=time()
#     for n in range(nset):
#
#         #==========================================================#
#         A0[:, 0] = 1
#         A0[:, 1] = tfit[frmin-1: frmax]
#         A0[:, 2] = 0.5*tfit[frmin-1: frmax]*tfit[frmin-1: frmax]
#         b0=self.meanResSets[n, frmin-1:frmax]*1e-9
#
#         x, covar, m02, std, stdX, res, m0=Fall.computeLST(A0, b0, frmin-1, frmax)
#         mx0.append(m0*1e8)
#         xg.append(x[0][2])
#         #==========================================================#
#
#
#         ressets2[:, n] = (self.meanResSets[n, :frmaxplot] - x[0][0] - x[0][1]*tfit[:frmaxplot] - 0.5*xg[n]*tfit[:frmaxplot]*tfit[:frmaxplot])*1e-9
#
#         #==========================================================#
#
#         AL1[:, 0] = 1
#         AL1[:, 1] = tfit[frmin-1: frmax]
#         AL1[:, 2] = self.v0m_bysets[n]/6*tfit[frmin-1: frmax]**3 + self.g0m_bysets[n]/24*tfit[frmin-1: frmax]**4 - hefm/2*tfit[frmin-1: frmax]*tfit[frmin-1: frmax] - hefm*gamma/24*tfit[frmin-1: frmax]**4
#         bL1 = ressets2[frmin-1:frmax, n]
#
#         x1, covar, m02, std, stdX, res, m0=Fall.computeLST(AL1, bL1, frmin-1, frmax)
#         xvgg.append((gamma + x1[0][2])*1e6)
#         mxL1.append(m0*1e6)
#         #==========================================================#
#
#         #calculating residuals and their std in range 3-16 cm
#         LLL=[]
#         for i in range(10*Lmin, 10*Lmax+1):
#             LL=i/1000
#             LLL.append(LL)
#
#             LA[:, 0] = 1
#             LA[:, 1] = np.sin(2*np.pi*self.zzh[frmin-1: frmax, n]/LL)
#             LA[:, 2] = np.cos(2*np.pi*self.zzh[frmin-1: frmax, n]/LL)
#
#             b = ressets2[frmin-1: frmax, n]
#
#             x2, covar, m02, std, stdX, res, m0=Fall.computeLST(LA, b, frmin-1, frmax)
#
#             T0_m[i-10*Lmin, n] = x2[0][0]*1e9
#             A0_m[i-10*Lmin, n] = x2[0][1]*1e9
#             B0_m[i-10*Lmin, n] = x2[0][2]*1e9
#             harmonic_res[i-10*Lmin, n] = np.std(res)*1e9
#
#         #searching of minimal residuals std
#         minimum = min(harmonic_res[:, n])
#         minimum_index = np.argmin(harmonic_res[:, n])
#
#         T0min.append(T0_m[minimum_index, n])
#         A0min.append(A0_m[minimum_index, n])
#         B0min.append(B0_m[minimum_index, n])
#         LLmin.append(Lmin/100 - 0.001 + minimum_index/1000)
#         # print(B0_m[minimum_index, n])
#
#         amplitude = np.sqrt(A0min[-1]*A0min[-1] + B0min[-1]*B0min[-1])
#         harm_amp.append(amplitude)
#
#         phase = np.arctan2(B0min[-1], A0min[-1])*(180/np.pi)
#         if phase < 0:
#             harm_phase.append(phase + 360)
#         else:
#             harm_phase.append(phase)
#
#         #==========================================================#
#
#         bH = ressets2[frmin-1:frmax, n] - (T0min[n] - A0min[n]*np.sin(2*np.pi*self.zzh[frmin-1:frmax, n]/LLmin[n]) - B0min[n]*np.cos(2*np.pi*self.zzh[frmin-1:frmax, n]/LLmin[n]))/1e9
#         x1, covar, m02, std, stdX, res, m0=Fall.computeLST(A0, bH, frmin-1, frmax)
#         dg.append(x1[0][2]*1e8)
#         mdg.append(m0*1e8)
#
#         #==========================================================#
#         x1, covar, m02, std, stdX, res, m0=Fall.computeLST(AL1, res, frmin-1, frmax)
#         vgg.append((gamma + x1[0][2])*1e6)
#         mvgg.append(m0*1e6)
#
#
#     if self.files.isChecked():
#
#         a=res_final(path=self.projDirPath, header=headers['vgg_per_sets0'].format(self.delimiter), name=self.stationData['ProjName']+'_'+'vgg_per_sets0', delimiter = self.delimiter)
#         for i in range(len(xvgg)):
#             line=[i+1, xvgg[i], mxL1[i], xg[i]*1e8, mx0[i]]
#
#             a.printResult(line = roundList(line, round_line_ind['vgg_per_sets0']))
#
#         a.close()
#
#
#         a=res_final(path=self.projDirPath, header=False, name=self.stationData['ProjName']+'_'+'ressets_per', delimiter = self.delimiter)
#         for j in range(self.FG5X['frmaxplot']):
#             line=[j+1, tfit[j], tfit[j]-self.tkor[0], self.zzh[j,0]]
#
#             line1=[]
#             for n in range(nset):
#                 ll = ressets2[j,n]*1e9-T0min[n]-A0min[n]*np.sin(2*n.pi*self.zzh[i,n]/LLmin[n])-B0min[n]*np.cos(2*n.pi*self.zzh[j,n]/LLmin[n])
#                 line1.append(ll)
#
#             line.extend(line1)
#
#             a.printResult(line = line)
#
#         a.close()
