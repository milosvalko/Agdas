import numpy as np
from numpy import random
import os, glob
import csv
import sqlite3 as sql
from CONFIG import matrDatabase, SAE
from warning import Warning
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
from scipy.stats import norm
import subprocess
import scipy.interpolate as interp
import scipy.stats

# import matplotlib
# font = {'size' : 18}
# matplotlib.rc('font', **font)

script_path = os.path.dirname(os.path.realpath(__file__))

class Fall():
    """
    Class for fitting of the drop
    """

    def __init__(self):
        self.c = 2.99792458e+17
        self.keys = ['z0', 'v0', 'g0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6']
        self.kdis = False
        self.kimp = False
        self.ksae = False

    def set_ksol_k(self, ksol_k: bool):
        self.ksol_k = ksol_k

    def set_ksol(self, ksol: float):
        self.ksol = ksol

    def setFringe(self, times: list):
        """
        Set measured fringes

        Parameters
        ----------
        times : list
            list of fringes, fringes are representation by string
        """

        self.fringe = np.float_((times))

    def setLambda(self, Lambda: str):
        """
        Set wave length

        Parameters
        ----------
        Lambda : str
            Wave length of laser
        """
        self.Lambda = float(Lambda)

    def setScaleFactor(self, scaleFactor: str):
        """
        Set GuideCard Scale Factor from Project file

        Parameters
        ----------
        scaleFactor : str
        """
        self.scaleFactor = float(scaleFactor)

    def setMultiplex(self, multiplex: str):
        """
        Set GuideCard Multiplex from Project file

        Parameters
        ----------
        multiplex : str
        """
        self.multiplex = float(multiplex)

    def setGradient(self, grad: str):
        """
        Set Gradient from Project file

        Parameters
        ----------
        grad : str
        """
        self.gradient = -100 * float(grad) * 1e-8

    def setModulFreq(self, fmod: float):
        """
        Set Modulation frequency from GUI - fmodf

        Parameters
        ----------
        fmod : float
        """
        self.fmod = fmod

    def setLpar(self, Lpar: float):
        """
        Set Parasitic wave length from GUI - lpar

        Parameters
        ----------
        Lpar : float
        """
        self.Lpar = Lpar

    def setRubiFreq(self, freq: str):
        """
        Set Rubidium Frequency from Project file - Instrument data

        Parameters
        ----------
        freq : str
        """
        self.rubiFreq = float(freq)

    def setFrRange(self, frmin: int, frmax: int):
        """
        Set first and last fringe

        Parameters
        ----------
        frmin : int
            index of first fringe
        frmax : int
            index of last fringe
        """
        self.frmin = frmin - 1
        self.frmax = frmax

    def setKpar(self, kpar: bool):
        """
        The method set if parasitic wavelength correction should be included

        Parameters
        ----------
        kpar : bool
        """
        self.kpar = kpar

    def checkKDIS(self):
        """
        Enable cable dispersion correction
        """
        self.kdis = True

    def checkKIMP(self):
        """
        Enable impedance mismatch correction
        """
        self.kimp = True

    def checkKSAE(self):
        """
        Enable self-attraction correction
        """
        self.ksae = True

    def setLcable(self, Lcable: float):
        """
        Set length of L cabel - GUI - lcable_ar

        Parameters
        ----------
        Lcable : float

        """
        self.Lcable = Lcable

    def setAcable(self, Acable: float):
        """
        Set A cable from GONFIG

        Parameters
        ----------
        Acable : float

        """
        self.Acable = Acable

    def setPcable(self, Pcable):
        """
        Set P cable from GONFIG

        Parameters
        ----------
        Pcable : float

        """
        self.Pcable = Pcable

    @staticmethod
    def computeLST(A: np.ndarray, z: np.ndarray, frmin: int, frmax: int):
        """

        Parameters
        ----------
        A : np.ndarray
            matrix of derivations
        z : np.ndarray
            vector of measuring
        frmin : int
            first fringe
        frmax : int
            last fringe

        Returns
        -------
        x : tuple
            return value from numpy.linalg.lstsq()
        covar : np.ndarray
            covariance matrix
        m02 : float
            m02=res*res'/(frmax-frmin-k)
        std : float
            standard deviation of fit
        stdX : standard deviation of fitted unknow values
        res : np.ndarray
            residuals for fringes, res=z-A*X
        m0 : float
            m0 computed from fringes in [frmin, frmax] range
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
        # """
        # Compute 'z0','v0','g0','a1','a2','a3','a4','a5','a6' by least squares method.
        # For computing LST with non-zero gradient call method - LST(grad=True)
        # For computing LST with zero gradient call method - LST(grad=False)
        # """
        """
        This method compute all of fits by least square method and set all of important variable as variable of class Fall
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
            z1[i] = (self.Lambda / 2 * (i) * self.multiplex * self.scaleFactor)
            # correct fringe
            self.tt[i] = ((self.fringe[i]) * (1e7 / self.rubiFreq) + self.ksol * z1[i] / self.c)

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
                    # print(zsae_list[-1])

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
        self.x_grad, covar_grad, self.m02_grad, self.stdstd, self.std_grad, res_grad, m0withgradient = Fall.computeLST(
            A=A_grad, z=z, frmin=self.frmin, frmax=self.frmax)

        # from sandbox1 import max_ind
        # self.ind_covar = max_ind(covar_grad)

        self.xef, xefCovar, xefM02, xefStd, stdXX, xefRes, m00 = Fall.computeLST(A2, z, frmin=self.frmin,
                                                                                 frmax=self.frmax)

        self.res_grad1 = np.subtract(z1, np.matmul(A_grad1, self.x_grad[0]))  # residuals for all fringes

        ress = res[self.frmin:self.frmax]  # residuals for correct interval
        self.ssres = np.sqrt(np.dot(ress, ress) / (self.frmax - self.frmin + 1))  # is drop accepted value

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
        Effective height of measuring at top of the drop
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
        self.gTopCor = self.gTop + 10 * float(tide) + 10 * float(load) + 10 * float(baro) + 10 * polar

class projectFile():
    """
    Read Project file and create dictionary --stationData
        --instrumentData
        --processingResults
        --gravityCorrections
        --self.names_summary
        --self.units
    """

    def __init__(self, projectfile):
        """

        Parameters
        ----------
        projectfile : str
            path to Project file
        """
        self.projectfile = projectfile

    def read(self):
        """
        Read project file word by word and create 'self.file' list.
        """
        self.file = []
        self.names = []

        file = open(self.projectfile, 'r')
        self.file_lines = file.read().splitlines()
        file.close()

        i = 0
        for line in self.file_lines:
            line_split = line.split()
            self.file.extend(line_split)
            ind = [i for j in range(len(line_split))]
            self.names.extend(ind)
            i += 1

    def insert_to_names(self, dict: dict, i: str):
        """
        This method collects name of lines from project file for summary

        Parameters
        ----------
        dict : dict
            dictionary of values
        i : index of line in raw file

        Returns
        -------

        """
        ind = self.file_lines[self.names[i]].index(':')
        name = self.file_lines[self.names[i]][:ind + 1]
        self.names_summary[list(dict.keys())[-1]] = name

    def insert_to_units(self, dict, i):
        """
        This method collects units of values for summary

        Parameters
        ----------
        dict : dict
            dictionary of values
        i : index of line in raw file

        Returns
        -------

        """
        self.units[list(dict.keys())[-1]] = self.file[i+2]

    def createDictionary(self):
        """
        This method create dictionaries with data from Raw file

        Returns
        -------
        stationData : dict
        instrumentData : dict
        processingResults : dict
        gravityCorrections : dict
        self.names_summary : dict
        self.units : dict
        """

        stationData = {}
        instrumentData = {}
        processingResults = {}
        gravityCorrections = {}
        self.names_summary = {}
        self.units = {}
        for i in range(0, len(self.file)):

            if self.file[i] == 'Project' and self.file[i + 1] == 'Name:':
                stationData['ProjName'] = self.file[i + 2]
                self.insert_to_names(stationData, i)

            if self.file[i] == 'Name:' and self.file[i - 1] == 'Data':
                stationData['name'] = self.file[i + 1]
                self.insert_to_names(stationData, i)

            if self.file[i] == 'Code:':
                stationData['SiteCode'] = self.file[i + 1]
                self.insert_to_names(stationData, i)

            if self.file[i] == 'Lat:':
                stationData['lat'] = self.file[i + 1]
                self.names_summary['lat'] = 'Latitude:'
                self.units['lat'] = 'deg'

            if self.file[i] == 'Long:':
                stationData['long'] = self.file[i + 1]
                self.names_summary['long'] = 'Longitude:'
                self.units['long'] = 'deg'

            if self.file[i] == 'Elev:':
                stationData['elev'] = self.file[i + 1]
                self.names_summary['elev'] = 'Elevation:'
                self.units['elev'] = 'm'

            if self.file[i] == 'Height:' and self.file[i - 1] == 'Setup':
                stationData['setupHeight'] = self.file[i + 1]
                self.insert_to_names(stationData, i)
                self.insert_to_units(stationData, i)

            if self.file[i] == 'Height:' and self.file[i - 1] == 'Transfer' and self.file[i - 2] == 'cm':
                stationData['transferHeight'] = self.file[i + 1]
                self.insert_to_names(stationData, i)
                # self.insert_to_units(stationData, i)
                self.units['transferHeight'] = 'cm'

            if self.file[i] == 'Height:' and self.file[i - 1] == 'Actual':
                stationData['actualHeight'] = self.file[i + 1]
                self.insert_to_names(stationData, i)
                self.insert_to_units(stationData, i)

            if self.file[i] == 'Gradient:' and self.file[i - 1] == 'cm':
                stationData['gradient'] = self.file[i + 1]
                self.insert_to_names(stationData, i)
                self.insert_to_units(stationData, i)

            if self.file[i] == 'Pressure:' and self.file[i - 1] == 'Air':
                stationData['airPressure'] = self.file[i + 1]
                self.insert_to_names(stationData, i)
                self.insert_to_units(stationData, i)

            if self.file[i] == 'Factor:' and self.file[i - 1] == 'Admittance':
                stationData['barometricFactor'] = self.file[i + 1]
                self.insert_to_names(stationData, i)

            if self.file[i] == 'Coord:':
                stationData['polarX'] = self.file[i + 1]
                self.names_summary['polarX'] = 'X polar coordinate:'
                self.units['polarX'] = '"'

                stationData['polarY'] = self.file[i + 3]
                self.names_summary['polarY'] = 'Y polar coordinate:'
                self.units['polarY'] = '"'

            if self.file[i] == 'Filename:' and self.file[i - 1] == 'Potential':
                stationData['potentialFile'] = self.file[i + 1]
                self.insert_to_names(stationData, i)

            if self.file[i] == 'Filename:' and self.file[i - 1] == 'Factor':
                stationData['deltaFactorFile'] = self.file[i + 1]
                self.insert_to_names(stationData, i)

            if self.file[i] == 'Type:' and self.file[i - 1] == 'Meter':
                instrumentData['meterType'] = self.file[i + 1]
                self.insert_to_names(instrumentData, i)

            if self.file[i] == 'S/N:':
                instrumentData['meterS/N'] = self.file[i + 1]
                self.insert_to_names(instrumentData, i)

            if self.file[i] == 'Height:' and self.file[i - 1] == 'Factory':
                instrumentData['factoryHeight'] = self.file[i + 1]
                self.insert_to_names(instrumentData, i)
                self.insert_to_units(instrumentData, i)

            if self.file[i] == 'Frequency:' and self.file[i - 1] == 'Rubidium':
                instrumentData['rubiFreq'] = self.file[i + 1]
                self.insert_to_names(instrumentData, i)
                self.insert_to_units(instrumentData, i)

            if self.file[i] == 'Laser:' and self.file[i + 3] == 'ID:':
                instrumentData['laser'] = self.file[i + 1]
                self.insert_to_names(instrumentData, i)

            if self.file[i] == 'ID:':
                instrumentData['ID'] = self.file[i + 1]
                self.insert_to_names(instrumentData, i)
                self.insert_to_units(instrumentData, i)

                # instrumentData['ID_V'] = self.file[i + 4]
                # self.insert_to_names(instrumentData, i)
                # self.insert_to_units(instrumentData, i)

            if self.file[i] == 'IE:':
                instrumentData['IE'] = self.file[i + 1]
                self.insert_to_names(instrumentData, i)
                self.insert_to_units(instrumentData, i)

                # instrumentData['IE_V'] = self.file[i + 4]
                # self.insert_to_names(instrumentData, i)
                # self.insert_to_units(instrumentData, i)

            if self.file[i] == 'IF:':
                instrumentData['IF'] = self.file[i + 1]
                self.insert_to_names(instrumentData, i)
                self.insert_to_units(instrumentData, i)

                # instrumentData['IF_V'] = self.file[i + 4]
                # self.insert_to_names(instrumentData, i)
                # self.insert_to_units(instrumentData, i)

            if self.file[i] == 'IG:':
                instrumentData['IG'] = self.file[i + 1]
                self.insert_to_names(instrumentData, i)
                self.insert_to_units(instrumentData, i)

                # instrumentData['IG_V'] = self.file[i + 4]
                # self.insert_to_names(instrumentData, i)
                # self.insert_to_units(instrumentData, i)

            if self.file[i] == 'IH:':
                instrumentData['IH'] = self.file[i + 1]
                self.insert_to_names(instrumentData, i)
                self.insert_to_units(instrumentData, i)

                # instrumentData['IH_V'] = self.file[i + 4]
                # self.insert_to_names(instrumentData, i)
                # self.insert_to_units(instrumentData, i)

            if self.file[i] == 'II:':
                instrumentData['II'] = self.file[i + 1]
                self.insert_to_names(instrumentData, i)
                self.insert_to_units(instrumentData, i)

                # instrumentData['II_V'] = self.file[i + 4]
                # self.insert_to_names(instrumentData, i)
                # self.insert_to_units(instrumentData, i)

            if self.file[i] == 'IJ:':
                instrumentData['IJ'] = self.file[i + 1]
                self.insert_to_names(instrumentData, i)
                self.insert_to_units(instrumentData, i)

                # instrumentData['IJ_V'] = self.file[i + 4]
                # self.insert_to_names(instrumentData, i)
                # self.insert_to_units(instrumentData, i)

            if self.file[i] == 'Frequency:' and self.file[i - 1] == 'Modulation':
                instrumentData['modulFreq'] = self.file[i + 1]
                self.insert_to_names(instrumentData, i)
                self.insert_to_units(instrumentData, i)

            if self.file[i] == 'Date:':
                processingResults['date'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)

            if self.file[i] == 'Time:':
                processingResults['time'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)

            if self.file[i] == 'DOY:':
                processingResults['doy'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)

            if self.file[i] == 'Year:':
                processingResults['year'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)

            if self.file[i] == 'Offset':
                processingResults['timeOffset'] = self.file[i + 4]
                self.names_summary['timeOffset'] = 'Time Offset (D h:m:s):'

            if self.file[i] == 'Gravity:':
                processingResults['gravity'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)
                self.insert_to_units(processingResults, i)

            if self.file[i] == 'Scatter:':
                processingResults['setScatter'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)
                self.insert_to_units(processingResults, i)

            if self.file[i] == 'Precision:':
                processingResults['precision'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)
                self.insert_to_units(processingResults, i)

            if self.file[i] == 'Uncertainty:' and self.file[i - 1] == 'Total':
                processingResults['totalUncertainty'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)
                self.insert_to_units(processingResults, i)

            if self.file[i] == 'Collected:':
                processingResults['setsCollected'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)

            if self.file[i] == 'Processed:' and self.file[i - 1] == 'Sets':
                processingResults['setsProcessed'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)

            if self.file[i] == 'Processed:' and self.file[i - 1] == '#s':
                processingResults['processedSets'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)

            if self.file[i] == 'Processed:' and self.file[i - 1] == 'NOT' and self.file[i - 2] == 'Sets':
                processingResults['numNotProcessed'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)

            if self.file[i] == 'Drops/Set:':
                processingResults['dropsInSet'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)

            if self.file[i] == 'Accepted:':
                processingResults['acceptedDrops'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)

            if self.file[i] == 'Rejected:':
                processingResults['rejectedDrops'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)

            if self.file[i] == 'Acquired:':
                processingResults['totalFringes'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)

            if self.file[i] == 'Start:':
                processingResults['fringeStart'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)

            if self.file[i] == 'Fringes:':
                processingResults['processedFringes'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)

            if self.file[i] == 'Multiplex:':
                processingResults['multiplex'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)

            if self.file[i] == 'Factor:' and self.file[i - 1] == 'Scale':
                processingResults['scaleFactor'] = self.file[i + 1]
                self.insert_to_names(processingResults, i)

            if self.file[i] == '(ETGTAB):':
                gravityCorrections['earthTide'] = self.file[i + 1]
                self.insert_to_names(gravityCorrections, i)
                self.insert_to_units(gravityCorrections, i)

            if self.file[i] == 'Motion:' and self.file[i - 1] == 'Polar' and self.file[i+3] == 'Barometric':
                gravityCorrections['polarMotion'] = self.file[i + 1]
                self.insert_to_names(gravityCorrections, i)
                self.insert_to_units(gravityCorrections, i)

            if self.file[i] == 'Pressure:' and self.file[i - 1] == 'Barometric':
                gravityCorrections['baroPress'] = self.file[i + 1]
                self.insert_to_names(gravityCorrections, i)
                self.insert_to_units(gravityCorrections, i)

            if self.file[i] == 'Height:' and self.file[i + 2] == 'µGal':
                gravityCorrections['transferHeightGal'] = self.file[i + 1]
                self.insert_to_names(gravityCorrections, i)
                self.insert_to_units(gravityCorrections, i)

            if self.file[i] == 'Xo:':
                gravityCorrections['referenceXo'] = self.file[i + 1]
                self.insert_to_names(gravityCorrections, i)
                self.insert_to_units(gravityCorrections, i)

        return stationData, instrumentData, processingResults, gravityCorrections, self.names_summary, self.units


class rawFile():
    """
    Read Raw file
    """

    def __init__(self, rawfile: str):
        """

        Parameters
        ----------
        rawfile : str
            path of rawfile
        """
        self.rawfile = rawfile
        self.read()

    def read(self):
        """
        Read raw file from direction.
        """
        raw = open(self.rawfile, 'r')
        self.raw_lines = raw.read().splitlines()
        raw.close()

    def rawHeader2(self):
        """
        Return second line of header from raw file

        Returns
        -------
        header2 : list
        """
        return self.raw_lines[1].split()

    def rawHeader1(self):
        """
        Return first line of header from raw file

        Returns
        -------
        header1 : list
        """
        return self.raw_lines[0].split()

    def rawLines(self):
        """
        Return lines of raw file
        Returns
        -------
        lines : list
            All lines except headers
        """

        return self.raw_lines[2:]


class dropFile():
    """
    Read Drop file
    """

    def __init__(self, dropfile):
        """

        Parameters
        ----------
        dropfile : str
            path of Drop file
        """

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
        Returns
        -------
        header4 : list
            4. header
        """
        return self.drop_lines[3].split()

    def dropLines(self):
        """
        Return lines of drop file
        Returns
        -------
        lines : list
            All lines except headers
        """

        return self.drop_lines[4:]


class estim():
    """
    Write estim file
    """

    def __init__(self, path, name):
        """

        Parameters
        ----------
        path : str
            file path
        name : str

        """
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
    """
    Class for writing outputs
    """

    def __init__(self, path: str, header: list, name: str, files='/Files/', delimiter=','):
        """

        Parameters
        ----------
        path : str
            Path of output
        header : list
            Headers splitted by delimiter
        name : str
            Project name
        files
        delimiter
        """

        self.delimiter = delimiter

        try:
            self.f = open(path + files + name + '.csv', 'w', newline='')
        except PermissionError:
            Warning(error='Close ' + name + '.csv' + '!', icon='critical', title='Warning')
            self.f = open(path + files + name + '.csv', 'w', newline='')

        if isinstance(header, str):
            self.f.write(header + '\n')
        # self.f.write(units+'\n')

    def printResult(self, line: list):
        """
        Write row (as list) of output file.

        Parameters
        ----------
        line : list
            Row of output file

        """

        writer = csv.writer(self.f, delimiter=self.delimiter)
        writer.writerow(line)

    def write_line(self, line: str):
        """
        Write row (as string) of output file.

        Parameters
        ----------
        line : str
            line of output file formated by delimiter
        """
        self.f.write(line + '\n')

    def close(self):
        """
        Close connection with the file
        """
        self.f.close()


class matr_db():
    """
    Class for work with output database
    """

    def __init__(self, path: str):
        """
        Create connection with matr database.

        Parameters
        ----------
        path : str
            path of output database
        """

        self.matr_db = sql.connect(path)
        self.matr_db.enable_load_extension(True)
        self.matr_db.load_extension(os.path.join(script_path, 'math.dll'))
        self.cursor = self.matr_db.cursor()
        try:
            self.cursor.execute(matrDatabase['schema'])
        except:
            pass
        self.insert('DELETE FROM results')

    def insert(self, query: str):
        """
        Execute query.

        Parameters
        ----------
        query : str
            sql query

        """
        self.cursor.execute(query)

    def get(self, query: str):
        """
        Return database request.

        Parameters
        ----------
        query : str
            query

        Returns
        -------
        res : list
            list of tuples with request
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
    """
    Class for graph plotting
    """

    def __init__(self, path: str, name: str, project: str, x_label: str, y_label: str, title: str, show: bool, winsize=(11, 4)):
        """

        Parameters
        ----------

        path : str
            path of graph
        name : str
            name of graph
        project : str
            project name
        x_label : str
            label of x ax
        y_label : str
            label of y ax
        title : str
            title of graph
        show : bool
            open graph after plot
        winsize : tuple
            size of graph
        """
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

    def decimal(self, decimal_number: int):
        """
        Set decimal number for y ax

        Parameters
        ----------
        decimal_number : int
            decimal number

        """
        ax = self.gr.gca()
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.{}f'.format(decimal_number)))

    def x_ax_int(self):
        """
        Set integers at x ax
        """

        ax = self.gr.gca()
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    def plotXY(self, x: list, y: list, mark: list, columns_name=[], legend=[], lw=0):
        """

        Parameters
        ----------
        x : list
            list of lists, x = [x1, ..., x2]
        y : list
            list of lists, y = [y1, ..., y2]
        mark : list
            list of markers
        columns_name
        legend : list
            list of legends
        lw : float
            line width

        """
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

    def histogram(self, hist_data: list, fit: bool):
        """
        Plot histogram from data in list

        Parameters
        ----------
        hist_data : list
            input data
        fit : bool
            True - fit by Gaussian curve

        """
        # Create histogram to plot
        self.hist_data = hist_data

        # Count of histogram columns by Sturges rule
        # bins=np.floor(1+3.32*np.log(len(hist_data)))
        bins = np.floor(1 + 5 * np.log(len(hist_data)))

        _, bins, _ = self.gr.hist(hist_data, int(bins), density=False, alpha=1, edgecolor='black')

        # fitting graph by Gaussian curve
        if fit:
            mean = np.mean(hist_data)
            variance = np.var(hist_data)
            sigma = np.sqrt(variance)
            dx = bins[1] - bins[0]
            scale = len(hist_data) * dx
            best_fit_line = scipy.stats.norm.pdf(bins, mean, sigma)
            self.gr.plot(bins, best_fit_line*scale, '-r')

    def error_bar(self, x_err: list, y_err: list, yerr: list, color_err: str, ms=10, capsize=5):
        """
        Plot error bar graph

        Parameters
        ----------
        x_err : list
            x coordinates
        y_err : list
            y coordinates
        yerr : list
            length of error bar
        color_err : str
            color of points
        ms : float
            marker size
        capsize : float
            size of

        """
        # Create error bar
        self.x_err = x_err
        self.y_err = y_err
        self.yerr = yerr

        self.gr.errorbar(x_err, y_err, yerr, marker='o', color=color_err, ms=ms, linestyle='', capsize=capsize,
                         zorder=0)

    def text(self, x: list, y: list, t: list, c: list):
        """
        Plot text to graph

        Parameters
        ----------
        x : list
            list of lists x coordinates, x = [x1, ..., x2]
        y : list
            list of lists y coordinates, y = [y1, ..., y2]
        t : list
            list of lists texts, t = [t1, ..., t2]
        c : list
            list of lists colors, c = [c1, ..., c2]

        """
        # Text to plot
        for i in range(len(x)):
            self.gr.text(x[i], y[i], t[i], color=c[i])

    def ylim(self, lim: list):
        """
        Set y limits

        Parameters
        ----------
        lim : list
            y limits = [lim1,  lim2]

        """
        self.gr.ylim(lim)

    def save(self):
        """
        Save graph to path
        """
        # Save plot to direction as png
        self.gr.xlabel(self.x_label)
        self.gr.ylabel(self.y_label)
        self.gr.title(self.title)
        self.gr.savefig(self.path + '/' + self.project + '_' + self.name)
        if self.show:
            # self.gr.show()
            ret_code = subprocess.call(['start', self.path + '/' + self.project + '_' + self.name + '.png'], shell=True)

        self.gr.close()


class Compare_gsoft_agdas():
    """
    This class compares values of gravity from G software and Agdas.
    The class generates compare_gsoft_agdas.csv file in Files folder with mean, standard deviation,
    max difference and differences between Agdas an g software
    """

    def __init__(self, path: str, vgg: str, project: str):
        """

        Parameters
        ----------
        path : str
            path of output folder
        vgg : str
            Gradient from Project file
        project :
            Project name
        """
        self.gsoft = [] #
        self.agdas = []
        self.set = []
        self.drp = []
        self.vgg = float(vgg) # gradient from project file
        self.path = os.path.join(path, 'Files', project + 'compare_gsoft_agdas.csv')
        self.path_hist = os.path.join(path, 'Graphs')
        self.project = project

    def add_gsoft(self, drop: dict, hef: float):
        """
        Reduce corrections of gsoft results and add to gsoft list

        Parameters
        ----------
        drop : dict
            line of drop file in dictionary structure
        hef : float
            g in effective height

        """
        self.set.append(drop['Set'])
        self.drp.append(drop['Drp'])

        # reduce corrections
        ffa_gsoft = float(drop['Gravity']) - float(drop['Tide']) - float(drop['Load']) - float(
            drop['Baro']) - float(
            drop['Polar']) - float(drop['Transfer'])
        # compute g in effective height
        ffa_gsoft_hef = ffa_gsoft - self.vgg * hef / 10

        self.gsoft.append(ffa_gsoft_hef)

    def add_agdas(self, EfH: float):
        """
        Add g to agdas
        Parameters
        ----------
        EfH : float
            g in effective height

        """
        self.agdas.append(EfH / 10)

    def print_file(self, acc: list, delimiter: str):
        """
        Write file with differences between Agdas and g-soft

        Parameters
        ----------
        acc : list
            list of tuples, bool, true - accepted
        delimiter : str
            delimiter for writing output file

        """

        # compute absolute value of differences for accepted drops
        self.diff_abs = [abs(self.gsoft[i] - self.agdas[i]) for i in range(len(self.gsoft)) if acc[i][0] == 1]
        self.diff_ = [self.gsoft[i] - self.agdas[i] for i in range(len(self.gsoft))]
        self.diff__ = [(self.gsoft[i] - self.agdas[i]) for i in range(len(self.gsoft)) if acc[i][0] == 1]

        # create summary
        summ = """max_diff_acc [nm.s-2]{}{:.2f}\nmin_diff_acc [nm.s-2]{}{:.2f}\navg_diff_acc [nm.s-2]{}{:.2f}\nstd_diff_acc [nm.s-2]{}{:.2f}\n
        """.format(delimiter, np.max(self.diff__)*10, delimiter, np.min(self.diff__)*10, delimiter, np.mean(self.diff__)*10, delimiter, np.std(self.diff__, ddof=1)*10)

        # create header
        # set ; drop ; acc ; gsoft ; agdas ; diff
        header = 'set{0}drop{0}acc{0}gsoft [nm.s-2]{0}agdas [nm.s-2]{0}diff\n'.format(delimiter)
        line = '{}{}{}{}{}{}{:.3f}{}{:.3f}{}{:.3f}\n'

        # open a create file
        file = open(self.path, 'w')
        file.write(summ)
        file.write(header)

        for i in range(len(self.gsoft)):
            line_ = line.format(self.set[i], delimiter, self.drp[i], delimiter, acc[i][0], delimiter,
                                self.gsoft[i]*10, delimiter, self.agdas[i]*10, delimiter, self.diff_[i]*10)
            file.write(line_)

        file.close()

    def print_histogram(self, graph_lang: list):
        """
        Plot histogram of differences

        Parameters
        ----------
        graph_lang : list
            differences

        """
        g = Graph(path=self.path_hist, name=self.project + '_histogram_compare_gsoft_pyagdas', project='',
                  x_label=graph_lang['histogram_diff']['xlabel'], y_label=graph_lang['histogram_diff']['ylabel'],
                  title=graph_lang['histogram_diff']['title'],
                  show=False)
        g.histogram(self.diff__*10, fit=True)
        # g.saveSourceData()
        g.save()
