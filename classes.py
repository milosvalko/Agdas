import numpy as np
from numpy import random
import os, glob
from functions import printDict
import csv
import sqlite3 as sql
from CONFIG import matrDatabase, SAE
from warning import Warning
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.stats import norm
import subprocess
import scipy.interpolate as interp


class Fall():

    def __init__(self):
        self.ksol = 1
        self.c = 2.99792458e+17
        self.keys = ['z0', 'v0', 'g0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6']
        self.kdis = False
        self.kimp = False
        self.ksae = False

    def setFringe(self, times):
        """Short summary.

        Parameters
        ----------
        times : type
            Description of parameter `times`.
            -Convert fringe in list from string to float
        Returns
        -------
        type
            -Set to class variable


        """

        self.fringe = np.float_((times))
        # self.fringe=[float(i) for i in times]

    def setLambda(self, Lambda):
        """Short summary.

        Parameters
        ----------
        Lambda : type
            Description of parameter `Lambda`.
            -
        Returns
        -------
        type
            Description of returned object.
            -Set to class variable
        """
        self.Lambda = float(Lambda)

    def setScaleFactor(self, scaleFactor):
        self.scaleFactor = float(scaleFactor)

    def setMultiplex(self, multiplex):
        self.multiplex = float(multiplex)

    def setGradient(self, grad):
        self.gradient = -100 * float(grad) * 1e-8

    def setModulFreq(self, fmod):
        self.fmod = float(fmod)

    def setLpar(self, Lpar):
        self.Lpar = Lpar

    def setRubiFreq(self, freq):
        self.rubiFreq = float(freq)

    def setFrRange(self, frmin, frmax):
        self.frmin = frmin - 1
        self.frmax = frmax

    def setFRssRange(self, frmaxss, frminss):
        self.frmaxss = frmaxss
        self.frminss = frminss - 1

    def setKpar(self, kpar):
        self.kpar = kpar

    def checkKDIS(self):
        self.kdis = True

    def checkKIMP(self):
        self.kimp = True

    def checkKSAE(self):
        self.ksae = True

    def setLcable(self, Lcable):
        self.Lcable = Lcable

    def setAcable(self, Acable):
        self.Acable = Acable

    def setPcable(self, Pcable):
        self.Pcable = Pcable

    @staticmethod
    def computeLST(A, z, frmin, frmax):
        """Short summary.

        Parameters
        ----------
        A : type
            Description of parameter `A`.
            -matrix A
        z : type
            Description of parameter `z`.
            -vector of meassurement

        Returns
        -------
        type
            Description of returned object.

        """

        x = np.linalg.lstsq(A, z, rcond=None)  # solution of LST
        covar = np.matrix(np.dot(A.T, A)).I  # covariance matrix
        m02 = x[1] / (frmax - frmin - A.shape[1] - 1)  # m02=res*res'/(frmax-frmin-k)
        stdX = np.sqrt(np.multiply(np.diag(covar), m02))  # standart devivations of x
        std = np.dot(stdX, stdX)
        res = np.subtract(z, np.matmul(A, x[0]))  # res=z-A*X
        m0 = np.sqrt(np.dot(res[frmin:frmax], res[frmin:frmax]) / (frmax - frmin) * (covar[2, 2]))

        return x, covar, m02, std, stdX, res, m0

    def LST(self):
        """
        Compute 'z0','v0','g0','a1','a2','a3','a4','a5','a6' by least squares method.
        For computing LST with non-zero gradient call method - LST(grad=True)
        For computing LST with zero gradient call method - LST(grad=False)
        """

        # count of fringe use for computing
        nfringe = len(self.fringe)

        # interpolated SAE for measuring
        if self.ksae:
            asae = interp.pchip([SAE[i] * 1e9 for i in range(0, len(SAE), 2)],
                                [SAE[i] * 10 for i in range(1, len(SAE), 2)])

        # initialization of z, tt, A
        z1 = np.zeros(nfringe)  # vector of z coordinates
        self.tt = np.zeros(nfringe)  # vector of fringe time

        A1 = np.zeros((nfringe, 7 + 2 * self.kpar))  # matrix to LST without gradient
        A_grad1 = np.zeros((nfringe, 7 + 2 * self.kpar))  # matrix to LST with gradient
        A4 = np.zeros((nfringe, 3 + 2 * self.kpar))  # matrix to LST with modulation

        dv_list = []
        asae_list = []
        vsae_list = []
        zsae_list = []
        dz_list = []

        # loop for all fringes
        for i in range(nfringe):

            # vector of z
            z1[i] = ((self.Lambda / 2 * (i) * self.multiplex * self.scaleFactor))
            # correct fringe
            self.tt[i] = (((self.fringe[i]) * (1e7 / self.rubiFreq) + self.ksol * z1[i] / self.c))

            # interpolated SAE
            if self.ksae:
                asae_interp = asae(z1[i])
                asae_list.append(asae_interp)
                if i == 0:
                    dv = asae_interp / 2 * self.tt[0]
                    z1[i] -= dv

                    dv_list.append(dv)
                    vsae_list.append(dv)
                    dz = vsae_list[0] / 2 * self.tt[0]
                    dz_list.append(dz)
                    zsae_list.append(dz * 1e-8)

                if i > 0:
                    dv = (asae_list[i - 1] + asae_list[i]) / 2 * (self.tt[i] - self.tt[i - 1])
                    vsae_list.append(vsae_list[-1] + dv)
                    dz_list.append((vsae_list[i - 1] + vsae_list[i]) / 2 * (self.tt[i] - self.tt[i - 1]))
                    zsae_list.append(zsae_list[i - 1] + dz_list[-1])
                    z1[i] -= zsae_list[-1]
                    print(zsae_list[-1])

            # corrections from dispersion and distortion
            fcefr = self.tt[i] * 9.8093 / (self.Lambda / 2 / 10e8)
            if self.kdis:
                z1[i] += (178.2 - 4.47 * np.log10(fcefr / 10e5) + 1.3 * (np.log10(fcefr / 10e5)) ** 2) / 36 * 9.8093 * \
                         self.tt[i] * self.Lcable

            if self.kimp:
                z1[i] += 9.8093 * self.tt[i] * 10e8 / (2 * np.pi * fcefr) * np.arcsin(
                    self.Acable * np.sin(4 * np.pi * fcefr * self.Lcable / 0.66 / 299792458 + self.Pcable))

            # arguments for sin and cos functions
            arg1 = 2 * np.pi * self.fmod * self.tt[i]
            arg2 = 2 * np.pi * self.fmod * 2 * z1[i] / self.c

            # derivations
            z_z0 = 1

            # with gradient
            z_v0_grad = self.tt[i] + (self.gradient / 6) * (self.tt[i]) ** 3
            z_g0_grad = (self.tt[i] ** 2) / 2 + (self.gradient / 24) * self.tt[i] ** 4

            # without gradient
            z_v0 = self.tt[i]
            z_g0 = (self.tt[i] ** 2) / 2

            # A matrix elements
            z_a1 = np.sin(arg1) * np.cos(arg2)
            z_a2 = np.cos(arg1) * np.cos(arg2)
            z_a3 = np.sin(arg1) * np.sin(arg2)
            z_a4 = np.cos(arg1) * np.sin(arg2)

            if self.kpar:
                z_a5 = np.sin(2 * np.pi * z1[i] / self.Lpar)
                z_a6 = np.cos(2 * np.pi * z1[i] / self.Lpar)

                line = [z_z0, z_v0, z_g0, z_a3, z_a1, z_a4, z_a2, z_a5, z_a6]

                line_grad = [z_z0, z_v0_grad, z_g0_grad, z_a3, z_a1, z_a4, z_a2, z_a5, z_a6]

                A4[i, 3] = z_a5
                A4[i, 4] = z_a6

            else:
                line = [z_z0, z_v0, z_g0, z_a3, z_a1, z_a4, z_a2]
                line_grad = [z_z0, z_v0_grad, z_g0_grad, z_a3, z_a1, z_a4, z_a2]

            # line of A matrix
            # line=[z_z0, z_v0, z_g0, z_a1, z_a2, z_a3, z_a4, z_a5, z_a6]

            A_grad1[i, :] = line_grad

            A4[i, 0] = 1
            A4[i, 1] = self.tt[i]

            A1[i, :] = line

        # separate
        A_grad = A_grad1[self.frmin:self.frmax, :]
        A = A1[self.frmin:self.frmax, :]
        z = z1[self.frmin:self.frmax]
        A2 = A[:, :7 + 2 * self.kpar]

        # A4=np.zeros([self.frmax-self.frmin+1,5])
        # A4[:,0:2]=A1[:,0:2]

        # Fit without gradient
        self.x, covar, self.m02, self.std, self.stdX, res, m0withouGR = Fall.computeLST(A=A, z=z, frmin=self.frmin,
                                                                                        frmax=self.frmax)

        # Fit with gradient
        self.x_grad, covar_grad, self.m02_grad, stdstd, self.std_grad, res_grad, m0withgradient = Fall.computeLST(
            A=A_grad, z=z, frmin=self.frmin, frmax=self.frmax)

        self.xef, xefCovar, xefM02, xefStd, stdXX, xefRes, m00 = Fall.computeLST(A2, z, frmin=self.frmin,
                                                                                 frmax=self.frmax)

        self.res_grad1 = np.subtract(z1, np.matmul(A_grad1, self.x_grad[0]))  # residuals for all fringes

        ress = self.res_grad1[self.frminss:self.frmaxss]  # residuals for correct interval
        self.ssres = np.sqrt(np.dot(ress, ress) / (self.frmaxss - self.frminss + 1))  # is drop accepted value

        # LST with modulation
        A4[:, 2] = [(self.x[0][1] / 6 * self.tt[i] ** 3 + self.x[0][2] / 24 * self.tt[i] ** 4 - (
                self.x[0][2] - self.x_grad[0][2]) * self.tt[i] ** 2 / (self.gradient * 2)) / 1e6 for i in
                    range(nfringe)]
        zgrad = np.zeros((nfringe))
        zgrad[:] = [z1[i] - self.x[0][2] / 2 * self.tt[i] ** 2 for i in range(nfringe)]
        modkor = np.matmul(A_grad1[:, 3:], self.x[0][3:])
        zgrad = np.subtract(zgrad, modkor.transpose())
        AA4 = A4[self.frmin:self.frmax, :]
        zgrad4 = zgrad[self.frmin:self.frmax]

        # Gradient estimation
        self.xgrad4, covarXgrad4, self.m0grad4, stdGrad4, self.stdGradX, self.Resgrad4, self.m0gradient = Fall.computeLST(
            A=AA4, z=zgrad4, frmin=self.frmin, frmax=self.frmax)

        self.resgrad4 = np.subtract(zgrad, np.matmul(A4, self.xgrad4[0]))

        # printDict(dict(zip(self.keys,x_grad[0])))
        # printDict(dict(zip(self.keys,x[0])))

        self.g0_Gr = self.x_grad[0][2]
        self.z0 = self.x[0][0]
        self.v0 = self.x[0][1]

        self.g0 = self.x[0][2]

    def effectiveHeight(self):
        """
        Effective height of measuring
        """
        self.h = -(self.g0_Gr - self.g0) / self.gradient

    def effectiveHeightTop(self):
        """
        Effective height of measuring in top
        """
        self.Grad = self.v0 ** 2 / (2 * self.g0_Gr)
        self.htop = self.h + self.Grad

    def effectivePosition(self):
        """
        Effective position of free fall
        """
        self.effectiveZ = self.h + self.z0

    def gTop(self):
        """
        G in effective height
        """
        self.gTop = self.g0_Gr - self.gradient * self.Grad

    def gTopCor(self, tide, load, baro, polar):
        """
        Load correction from tide, load, baro and polar to gTop
        """
        self.gTopCor = self.gTop + 10 * float(tide) + 10 * float(load) + 10 * float(baro) + 10 * float(polar)


class projectFile():

    def __init__(self, projectfile):
        self.projectfile = projectfile

    def read(self):
        """
        Read project file word by word and create 'self.file' list
        """
        self.file = []

        with open(self.projectfile, "r") as f:
            for line in f:
                self.file.extend(line.split())

    def createDictionary(self):
        """
        Read 'self.file' list and create dictionaries:
        --stationData
        --instrumentData
        --processingResults
        --gravityCorrections
        """
        stationData = {}
        instrumentData = {}
        processingResults = {}
        gravityCorrections = {}
        for i in range(0, len(self.file)):

            if self.file[i] == 'Project' and self.file[i + 1] == 'Name:':
                stationData['ProjName'] = self.file[i + 2]

            if self.file[i] == 'Name:' and self.file[i - 1] == 'Data':
                stationData['name'] = self.file[i + 1]

            if self.file[i] == 'Code:':
                stationData['SiteCode'] = self.file[i + 1]

            if self.file[i] == 'Lat:':
                stationData['lat'] = self.file[i + 1]

            if self.file[i] == 'Long:':
                stationData['long'] = self.file[i + 1]

            if self.file[i] == 'Elev:':
                stationData['elev'] = self.file[i + 1]

            if self.file[i] == 'Height:' and self.file[i - 1] == 'Setup':
                stationData['setupHeight'] = self.file[i + 1]

            if self.file[i] == 'Height:' and self.file[i - 1] == 'Transfer' and self.file[i - 2] == 'cm':
                stationData['transferHeight'] = self.file[i + 1]

            if self.file[i] == 'Height:' and self.file[i - 1] == 'Actual':
                stationData['actualHeight'] = self.file[i + 1]

            if self.file[i] == 'Gradient:' and self.file[i - 1] == 'cm':
                stationData['gradient'] = self.file[i + 1]

            if self.file[i] == 'Pressure:' and self.file[i - 1] == 'Air':
                stationData['airPressure'] = self.file[i + 1]

            if self.file[i] == 'Factor:' and self.file[i - 1] == 'Admittance':
                stationData['barometricFactor'] = self.file[i + 1]

            if self.file[i] == 'Coord:':
                stationData['polarX'] = self.file[i + 1]
                stationData['polarY'] = self.file[i + 3]

            if self.file[i] == 'Filename:' and self.file[i - 1] == 'Potential':
                stationData['potentialFile'] = self.file[i + 1]

            if self.file[i] == 'Filename:' and self.file[i - 1] == 'Factor':
                stationData['deltaFactorFile'] = self.file[i + 1]

            if self.file[i] == 'Type:' and self.file[i - 1] == 'Meter':
                instrumentData['meterType'] = self.file[i + 1]

            if self.file[i] == 'S/N:':
                instrumentData['meterS/N'] = self.file[i + 1]

            if self.file[i] == 'Height:' and self.file[i - 1] == 'Factory':
                instrumentData['factoryHeight'] = self.file[i + 1]

            if self.file[i] == 'Frequency:' and self.file[i - 1] == 'Rubidium':
                instrumentData['rubiFreq'] = self.file[i + 1]

            if self.file[i] == 'Laser:' and self.file[i + 3] == 'ID:':
                instrumentData['laser'] = self.file[i + 1]

            if self.file[i] == 'ID:':
                instrumentData['ID'] = self.file[i + 1]
                instrumentData['ID_V'] = self.file[i + 4]

            if self.file[i] == 'IE:':
                instrumentData['IE'] = self.file[i + 1]
                instrumentData['IE_V'] = self.file[i + 4]

            if self.file[i] == 'IF:':
                instrumentData['IF'] = self.file[i + 1]
                instrumentData['IF_V'] = self.file[i + 4]

            if self.file[i] == 'IG:':
                instrumentData['IG'] = self.file[i + 1]
                instrumentData['IG_V'] = self.file[i + 4]

            if self.file[i] == 'IH:':
                instrumentData['IH'] = self.file[i + 1]
                instrumentData['IH_V'] = self.file[i + 4]

            if self.file[i] == 'II:':
                instrumentData['II'] = self.file[i + 1]
                instrumentData['II_V'] = self.file[i + 4]

            if self.file[i] == 'IJ:':
                instrumentData['IJ'] = self.file[i + 1]
                instrumentData['IJ_V'] = self.file[i + 4]

            if self.file[i] == 'Frequency:' and self.file[i - 1] == 'Modulation':
                instrumentData['modulFreq'] = self.file[i + 1]

            if self.file[i] == 'Date:':
                processingResults['date'] = self.file[i + 1]

            if self.file[i] == 'Time:':
                processingResults['time'] = self.file[i + 1]

            if self.file[i] == 'DOY:':
                processingResults['doy'] = self.file[i + 1]

            if self.file[i] == 'Year:':
                processingResults['year'] = self.file[i + 1]

            if self.file[i] == 'Offset':
                processingResults['timeOffset'] = self.file[i + 4]

            if self.file[i] == 'Gravity:':
                processingResults['gravity'] = self.file[i + 1]

            if self.file[i] == 'Scatter:':
                processingResults['setScatter'] = self.file[i + 1]

            if self.file[i] == 'Precision:':
                processingResults['precision'] = self.file[i + 1]

            if self.file[i] == 'Uncertainty:' and self.file[i - 1] == 'Total':
                processingResults['totalUncertainty'] = self.file[i + 1]

            if self.file[i] == 'Collected:':
                processingResults['setsCollected'] = self.file[i + 1]

            if self.file[i] == 'Processed:' and self.file[i - 1] == 'Sets':
                processingResults['setsProcessed'] = self.file[i + 1]

            if self.file[i] == 'Processed:' and self.file[i - 1] == '#s':
                processingResults['processedSets'] = self.file[i + 1]

            if self.file[i] == 'Processed:' and self.file[i - 1] == 'NOT' and self.file[i - 2] == 'Sets':
                processingResults['numNotProcessed'] = self.file[i + 1]

            if self.file[i] == 'Drops/Set:':
                processingResults['dropsInSet'] = self.file[i + 1]

            if self.file[i] == 'Accepted:':
                processingResults['acceptedDrops'] = self.file[i + 1]

            if self.file[i] == 'Rejected:':
                processingResults['rejectedDrops'] = self.file[i + 1]

            if self.file[i] == 'Acquired:':
                processingResults['totalFringes'] = self.file[i + 1]

            if self.file[i] == 'Start:':
                processingResults['fringeStart'] = self.file[i + 1]

            if self.file[i] == 'Fringes:':
                processingResults['processedFringes'] = self.file[i + 1]

            if self.file[i] == 'Multiplex:':
                processingResults['multiplex'] = self.file[i + 1]

            if self.file[i] == 'Factor:' and self.file[i - 1] == 'Scale':
                processingResults['scaleFactor'] = self.file[i + 1]

            if self.file[i] == '(ETGTAB):':
                gravityCorrections['earthTide'] = self.file[i + 1]

            if self.file[i] == 'Motion:' and self.file[i - 1] == 'Polar' and i < 460:
                gravityCorrections['polarMotion'] = self.file[i + 1]

            if self.file[i] == 'Pressure:' and self.file[i - 1] == 'Barometric':
                gravityCorrections['baroPress'] = self.file[i + 1]

            if self.file[i] == 'Height:' and self.file[i + 2] == 'µGal':
                gravityCorrections['transferHeight'] = self.file[i + 1]

            if self.file[i] == 'Xo:':
                gravityCorrections['referenceXo'] = self.file[i + 1]

        return stationData, instrumentData, processingResults, gravityCorrections


class rawFile():

    def __init__(self, rawfile):
        self.rawfile = rawfile
        self.read()

    def read(self):
        """
        Read raw file from direction
        """
        raw = open(self.rawfile, 'r')
        self.raw_lines = raw.read().splitlines()
        raw.close()

    def rawHeader2(self):
        """
        Return second line of header from raw file
        """
        return self.raw_lines[1].split()

    def rawHeader1(self):
        """
        Return first line of header from raw file
        """
        return self.raw_lines[0].split()

    def rawLines(self):
        """
        Return lines of raw file
        """
        return self.raw_lines[2:]


class dropFile():

    def __init__(self, dropfile):
        self.dropfile = dropfile
        self.read()

    def read(self):
        """
        Read drop file
        """
        drop = open(self.dropfile, 'r')
        self.drop_lines = drop.read().splitlines()
        drop.close()

    def dropHeader4(self):
        """
        Return fourth header from drop file
        """
        return self.drop_lines[3].split()

    def dropLines(self):
        """
        Return lines of drop file
        """
        return self.drop_lines[4:]


class estim():

    def __init__(self, path, name):
        header = 'Set {0} Drop {0} m0 {0} z0 {0} z0-std {0} v0 {0} v0-std {0} g0 {0} g0-std {0} a {0} a-std {0} b {0} b-std {0} c {0} c-std {0} d {0} d-std {0} e {0} e-std {0} f {0} f-std '.format(
            ';')
        units = '   {0}    {0} {0} mm {0} mm {0} mm.s-1 {0} mm.s-1 {0} nm.s-2 {0} nm.s-2 {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm   '.format(
            ';')

        try:
            os.remove(path + '/Files/estim.csv')
        except FileNotFoundError:
            pass

        try:
            self.f = open(path + '/Files/' + name + '_' + 'estim.csv', 'w', newline='')
        except PermissionError:
            Warning(error='Close ' + name + '_estim.csv' + '!', icon='critical', title='Warning')
            self.f = open(path + '/Files/' + name + '_' + 'estim.csv', 'w', newline='')

        self.f.write(header + '\n')
        self.f.write(units + '\n')
        # self.f.close()

    def printResult(self, X, std, set, drop, m0):
        """
        Print result of drop to estim file
        """
        dropResult = [set, drop, m0[0]]
        for i in range(len(X)):
            dropResult.append(X[i])
            dropResult.append(std[i])

        # dropResult[3] = round(dropResult[3], 4)
        # dropResult[4] = round(dropResult[4], 14)
        # dropResult[5] = round(dropResult[5], 4)
        # dropResult[6] = round(dropResult[6], 15)
        # dropResult[8] = round(dropResult[8], 3)

        writer = csv.writer(self.f, delimiter=';')
        writer.writerow(dropResult)


class res_final():

    def __init__(self, path, header, name, files='/Files/', delimiter=','):

        self.delimiter = delimiter
        # try:
        #     os.remove(path+files+name+'.csv')
        # except FileNotFoundError:
        #     pass
        # except PermissionError:
        #     Warning(error='Close ' + name + '.csv' + '!',icon='critical', title='Warning')
        #     os.remove(path+files+name+'.csv')

        try:
            self.f = open(path + files + name + '.csv', 'w', newline='')
        except PermissionError:
            Warning(error='Close ' + name + '.csv' + '!', icon='critical', title='Warning')
            self.f = open(path + files + name + '.csv', 'w', newline='')

        if isinstance(header, str):
            self.f.write(header + '\n')
        # self.f.write(units+'\n')

    def printResult(self, line):
        # line='{};{};{};{};{};{}'.format(str(i), str(z), str(t), str(tt), str(value), str(fil_value))
        # line=[i, z, t, tt, value, fil_value]

        # print(line)
        writer = csv.writer(self.f, delimiter=self.delimiter)
        writer.writerow(line)

    def close(self):
        self.f.close()


class matr_db():

    def __init__(self, path):
        """
        Create connection with matr database
        """
        self.matr_db = sql.connect(path)
        self.cursor = self.matr_db.cursor()
        try:
            self.cursor.execute(matrDatabase['schema'])
        except:
            pass
        self.insert('DELETE FROM results')

    def insert(self, query):
        """
        Execute insert query
        """
        self.cursor.execute(query)

    def get(self, query):
        """
        Return data from database
        """
        self.cursor.execute(query)
        res = self.cursor.fetchall()
        return res

    def commit(self):
        """
        Commit data to database
        """
        self.matr_db.commit()

    def close(self):
        """
        Close connection with database
        """
        self.matr_db.close()


class Graph():

    def __init__(self, path, name, project, x_label, y_label, title, show, winsize=(11, 4)):
        # Create graph
        self.gr = plt
        self.gr.rcParams['figure.dpi'] = 500
        self.gr.rcParams['figure.figsize'] = winsize

        self.path = path
        self.name = name
        self.project = project
        self.show = show
        self.x_label = x_label
        self.y_label = y_label
        self.title = title

    def plotXY(self, x, y, mark, columns_name=[], legend=[], lw=0):
        # Plot XY data
        self.x = x
        self.y = y
        self.columns_name = columns_name

        if lw == 0:
            for i in range(len(x)):
                self.gr.plot(x[i], y[i], mark[i], lw=0.3)
        else:
            for i in range(len(x)):
                self.gr.plot(x[i], y[i], mark[i], lw=lw[i])

        if len(legend) > 0:
            self.gr.legend(legend, loc='upper right')

    def saveSourceData(self):
        # Save source data to csv file
        err = False
        try:
            self.x.append(self.x_err)
            self.y.append(self.y_err)
            self.columns_name.append('err')
            err = True
        except:
            pass

        hist = False
        try:
            self.x.append(self.hist_data)
            # self.columns_name.append('hist')
            hist = True
        except:
            pass

        # Write graph source data to file
        try:
            d = open(self.path + '/' + self.project + '_' + self.name + '.csv', 'w', encoding='utf-8')
        except PermissionError:
            Warning(error='Close ' + self.project + '_' + self.name + '.csv' + '!', icon='critical', title='Warning')
            d = open(self.path + '/' + self.project + '_' + self.name + '.csv', 'w', encoding='utf-8')

        h = ''
        try:
            for i in self.columns_name:
                h += i + '_x'
                h += ';'
                h += i + '_y'
                h += ';'
            if err:
                h += 'err' + ';'
            if hist:
                h += 'hist'
            d.write(h + '\n')

            n = max([len(i) for i in self.x])

            for i in range(n):
                l = []
                for j in range(len(self.y)):

                    try:
                        l.append(self.x[j][i])
                    except IndexError:
                        self.x[j].append('')
                        l.append(self.x[j][i])

                    try:
                        l.append(self.y[j][i])
                    except IndexError:
                        self.y[j].append('')
                        l.append(self.y[j][i])

                ll = ''
                for k in l:
                    ll += str(k)
                    ll += ';'

                if err:
                    ll += str(self.yerr[i]) + ';'
                if hist:
                    ll += str(self.hist_data[i])
                ll += '\n'

                d.write(ll)
        except AttributeError:
            pass

        d.close()

    def histogram(self, hist_data, fit):
        # Create histogram to plot
        self.hist_data = hist_data

        # Count of histogram columns by Sturges rule
        # bins=np.floor(1+3.32*np.log(len(hist_data)))
        bins = np.floor(1 + 5 * np.log(len(hist_data)))

        a = self.gr.hist(hist_data, edgecolor='black', bins=int(bins))

        if fit:
            bins = np.linspace(min(hist_data), max(hist_data), 100)

            mu, sigma = norm.fit(hist_data)

            fit = norm.pdf(bins, mu, sigma)

            self.gr.plot(bins, fit * max(a[0]) / max(fit), 'r')

    def error_bar(self, x_err, y_err, yerr, color_err, ms=10, capsize=5):
        # Create error bar
        self.x_err = x_err
        self.y_err = y_err
        self.yerr = yerr

        self.gr.errorbar(x_err, y_err, yerr, marker='o', color=color_err, ms=ms, linestyle='', capsize=capsize,
                         zorder=0)

    def text(self, x, y, t, c):
        # Text to plot
        for i in range(len(x)):
            self.gr.text(x[i], y[i], t[i], color=c[i])

    def save(self):
        # Save plot to direction as png
        self.gr.xlabel(self.x_label)
        self.gr.ylabel(self.y_label)
        self.gr.title(self.title)
        self.gr.savefig(self.path + '/' + self.project + '_' + self.name)
        if self.show:
            # self.gr.show()
            ret_code = subprocess.call(['start', self.path + '/' + self.project + '_' + self.name + '.png'], shell=True)

        self.gr.close()
