import os
import numpy as np
from glob import glob
import sqlite3 as sql

# -----------------------------------------------------------------
path1 = r'c:\Users\Jakub\Desktop\pecny\data\italie\agdas_808_parOFF'
path2 = r'c:\Users\Jakub\Desktop\pecny\data\italie\agdas_808_parOFF\res\Files'
compare_path = r'c:\Users\Jakub\Desktop\pecny\data\italie\agdas_808_parOFF\comp'

# path1 = r'c:\Users\Jakub\Desktop\pecny\data\x153_glb\agdas_807zaloha'
# path2 = r'c:\Users\Jakub\Desktop\pecny\data\x153_glb\res1\Files'
# compare_path = r'c:\Users\Jakub\Desktop\pecny\data\x153_glb\res1'
# delimiter = ','

# path1 = r'c:\Users\Jakub\Desktop\pecny\data\566_glb\agdas_808'
# path2 = r'c:\Users\Jakub\Desktop\pecny\data\566_glb\res\Files'
# compare_path = r'c:\Users\Jakub\Desktop\pecny\data\566_glb\res'

# path1 = r'c:\Users\Jakub\Desktop\pecny\data\566_glb\agdas_808'
# path2 = r'c:\Users\Jakub\Desktop\pecny\data\566_glb\res\Files'
# compare_path = r'c:\Users\Jakub\Desktop\pecny\data\566_glb\res'
# delimiter = ','

# path1 = r'c:\Users\Jakub\Desktop\pecny\data\x172_wetzel\files_wetzel'
# path2 = r'c:\Users\Jakub\Desktop\pecny\data\x172_wetzel\res1\Files'
# compare_path = r'c:\Users\Jakub\Desktop\pecny\data\x172_wetzel\res1'
# delimiter = ','

# path1 = r'c:\Users\Jakub\Desktop\pecny\data\575_glb\agdas_902'
# path2 = r'c:\Users\Jakub\Desktop\pecny\data\575_glb\res\Files'
# compare_path = r'c:\Users\Jakub\Desktop\pecny\data\575_glb\res'
delimiter = ','

# path1 = r'c:\Users\Jakub\Desktop\pecny\data\572_glb\agdas_808'
# path2 = r'c:\Users\Jakub\Desktop\pecny\data\572_glb\res_m0\Files'
# compare_path = r'c:\Users\Jakub\Desktop\pecny\data\572_glb\res_m0'

# path1 = r'c:\Users\Jakub\Desktop\pecny\data\573_glb\agdas_901'
# path2 = r'c:\Users\Jakub\Desktop\pecny\data\573_glb\res\Files'
# compare_path = r'c:\Users\Jakub\Desktop\pecny\data\573_glb\res'

# path1 = r'c:\Users\Jakub\Desktop\pecny\data\x179_\agdas_901'
# path2 = r'c:\Users\Jakub\Desktop\pecny\data\x179_\res\Files'
# compare_path = r'c:\Users\Jakub\Desktop\pecny\data\x179_\res'

# path1 = r'c:\Users\Jakub\Desktop\pecny\data\x163_glb\agdas_808'
# path2 = r'c:\Users\Jakub\Desktop\pecny\data\x163_glb\res\Files'
# compare_path = r'c:\Users\Jakub\Desktop\pecny\data\x163_glb\res'

# path1 = r'c:\Users\Jakub\Desktop\pecny\data\x165_glb\agdas_808_A_parOFF'
# path2 = r'c:\Users\Jakub\Desktop\pecny\data\x165_glb\res_off\Files'
# compare_path = r'c:\Users\Jakub\Desktop\pecny\data\x165_glb\res_off'
#
# path1 = r'c:\Users\Jakub\Desktop\pecny\data\x165_glb\agdas_808_A_parON'
# path2 = r'c:\Users\Jakub\Desktop\pecny\data\x165_glb\res_on\Files'
# compare_path = r'c:\Users\Jakub\Desktop\pecny\data\x165_glb\res_on'

# path1 = r'c:\Users\Jakub\Desktop\pecny\data\x168_gyula\agdas_808'
# path2 = r'c:\Users\Jakub\Desktop\pecny\data\x168_gyula\res\Files'
# compare_path = r'c:\Users\Jakub\Desktop\pecny\data\x168_gyula\res'


# Compare(path1, path2, compare_path, delimiter)
out = True
matlab_nonacc = True
allan = False
summary = True


# -----------------------------------------------------------------
class Compare():

    def __init__(self, path1: str, path2: str, compare_path: str, delimiter: str):

        self.path1 = path1
        self.path2 = path2
        self.compare_path = compare_path
        self.delimiter = delimiter

        self.matlab_acc = 0

        self.Compare_dir()

    def summary(self, path_database, accepted):

        # ==================================================================#
        database = sql.connect(path_database)
        cursor = database.cursor()
        cursor.execute('select count(*) from results where Accepted = 1')
        acc = cursor.fetchall()
        cursor.execute('select count(*) from results')
        all = cursor.fetchall()

        cursor.execute('select n from results where Accepted = 0')
        rejected_python = cursor.fetchall()
        # ==================================================================#
        acc_matlab = 0
        for i in accepted:
            if i:
                acc_matlab += 1
        # ==================================================================#
        matlog = glob(self.path1 + '/*matlog.csv')
        matlog = open(matlog[0], 'r')
        matlog_r = matlog.read().splitlines()
        matlog.close()
        gravimeter = matlog_r[2].split(delimiter)[1]

        path = os.path.join(self.compare_path, 'summary.txt')
        summary_file = open(path, 'w')
        output = """Gravimeter: {}     
        Rejected/Accepted/All\nMatlab: {}   /   {}/   {}\nPython: {}   /   {}/   {}
        Rejected Matlab / Python \n
        """
        output = output.format(gravimeter, self.matlab_acc, all[0][0] - self.matlab_acc, all[0][0],
                               all[0][0] - acc[0][0], acc[0][0],
                               all[0][0])
        summary_file.write(output)

        if len(rejected_python) > len(self.rejected_matlab):
            d = len(rejected_python) - len(self.rejected_matlab)
            [self.rejected_matlab.append('-') for i in range(d)]
        else:
            d = len(self.rejected_matlab) - len(rejected_python)
            [rejected_python.append('-') for i in range(d)]

        for i in range(len(rejected_python)):
            line = '      {}   /   {}'.format(self.rejected_matlab[i], rejected_python[i][0]) + '\n'
            summary_file.write(line)

        summary_file.close()

    def Compare_dir(self):

        # create direction "comparisons"
        self.compare_path = os.path.join(self.compare_path, 'Comparisions')
        try:
            os.mkdir(self.compare_path)
        except:
            pass

        # read all files with ".csv" extension
        files1 = glob(self.path1 + '\*.csv')
        files2 = glob(self.path2 + '\*.csv')

        if allan:
            del files1[0]

        # get accepted
        path_database = os.path.realpath(os.path.join(os.path.dirname(path2), '.', 'data.db'))
        path_drops_file_matlab = glob(self.path1 + '\*drops.csv')[0]
        accepted = self.get_accepted(path_database, path_drops_file_matlab, delimiter)

        # compare each file
        for i in range(len(files1)):

            file1 = files1[i].split('\\')[-1]

            for j in range(len(files2)):

                file2 = files2[j].split('\\')[-1]

                if file1 == file2:
                    if 'matlog.' in file1:
                        headers1 = 2
                        headers2 = 2
                        delimiter1 = ','
                        delimiter2 = ','
                    else:
                        headers1, delimiter1 = Compare.get_headers_delimiter(files1[i])

                        headers2, delimiter2 = Compare.get_headers_delimiter(files2[j])

                    if 'estim' in file1 or 'drops' in file1:
                        comparision = Compare.compare2files(files1[i], files2[j], self.delimiter, delimiter1,
                                                            delimiter2,
                                                            headers1, headers2, accepted)

                    else:
                        comparision = Compare.compare2files(files1[i], files2[j], self.delimiter, delimiter1,
                                                            delimiter2,
                                                            headers1, headers2)

                    output_file = os.path.join(self.compare_path, 'Comp_' + file1)

                    if out:
                        f = open(output_file, 'w')
                        f.write(comparision)
                        f.close()

        if summary:
            self.summary(path_database, accepted)

    @staticmethod
    def match_columns(l1: list, l2: list):

        # matches - position of column from first file with corresponding column from second file
        matches = []

        # prohibited column to comparing
        prohibited = ['', 'Time', 'Campaign', 'Date', 'Accepted', 'Set', 'Drop', 'Year', 'Month', 'Day', 'Hour',
                      'Minute', 'Second', 'Campaign', 'Gravimeter-type', 'Sitename', 'Sitecode', 'Gravimeter-SN']

        # used_columns is used for 'estimgrad' because there are columns with the same name
        # key is name of column and value is frequency
        used_columns = {}
        # find corresponding columns
        for i in range(len(l1)):

            # if column was already used, increase frequency
            if l1[i].strip() in used_columns:
                used_columns[l1[i].strip()] += 1
            else:
                used_columns[l1[i].strip()] = 0

            # decide where shall start with columns
            if l1[i].strip() in used_columns and used_columns[l1[i].strip()] == 0:
                s = 0

            if l1[i].strip() in used_columns and used_columns[l1[i].strip()] == 1:
                s = 10

            if l1[i].strip() in used_columns and used_columns[l1[i].strip()] == 2:
                s = 20

            for j in range(s, len(l2)):

                if l1[i].strip().lower() == l2[j].strip().lower() and l1[i].strip() not in prohibited:
                    matches.append([i, j])
                    break

        return matches

    # @staticmethod
    def get_accepted(self, path_database: str, path_drops_file_matlab: str, delimiter: str):
        """

        @param path_database: path to database with results from pyAgdas
        @param path_drops_file_matlab: path to 'drops' file from matlab
        @param delimiter: delimiter of 'path_drops_file_matlab' file
        @return: list with boolean values if drop was accepted or not
        """

        database = sql.connect(path_database)
        cursor = database.cursor()
        cursor.execute('select Accepted, Set1, Drop1 from results')
        acc = cursor.fetchall()

        file = open(path_drops_file_matlab, 'r')
        file_r = file.read().splitlines()
        file.close()

        acc_index = file_r[0].split(delimiter).index('Acc')

        j = 0
        accepted = []
        self.rejected_matlab = []
        for i in range(len(file_r) - 1):

            a = file_r[i].split(delimiter)[acc_index]

            if a.strip() == '0':
                self.matlab_acc += 1
                self.rejected_matlab.append(i + 1)

            if a != 'Acc':

                acc_python = acc[j][0]

                if acc_python == 0 or a.strip() == '0':
                    accepted.append(False)
                else:
                    accepted.append(True)

                j += 1

        return accepted

    @staticmethod
    def compare2files(path1: str, path2: str, delimiter: str, delimiter1: str, delimiter2: str, headers1: int,
                      headers2: int, acc=True):

        # open files
        file1 = open(path1, 'r')
        file2 = open(path2, 'r')

        # read files
        file1_r = file1.read().splitlines()
        file2_r = file2.read().splitlines()

        # get headers
        if 'estimgrad' in path1:
            a = file1_r[1].split(delimiter1)
            b = file2_r[1].split(delimiter2)

        else:
            a = file1_r[0].split(delimiter1)
            b = file2_r[0].split(delimiter2)

        # units
        c = file1_r[1].split(delimiter1)

        # this is because 'residuals_sets' file has columns without names
        if 'residuals_sets' in path1:
            # create matches for residuals sets
            count_columns = len(file1_r[1].split(delimiter1))

            matches = []

            [matches.append([i, i]) for i in range(count_columns)]

            a = file1_r[0].split(delimiter1)
            b = file2_r[0].split(delimiter2)

            [a.append('-') for i in range(len(a), count_columns)]
            [b.append('-') for i in range(len(b), count_columns)]

        else:
            # find matches
            matches = Compare.match_columns(a, b)

        if 'matlogsets' in path1:
            matches.append([12, 12])
            matches.append([13, 13])

        # l1_p = len(file1_r[header1 + 2].split(delimiter1))
        n = len(file1_r) - headers1

        if acc == True:
            acc = []
            [acc.append(True) for i in range(n)]

        # array with differences between matlab and python
        diff = np.zeros((n, len(matches)))

        # find differences
        for i in range(n):

            l1 = file1_r[headers1 + i].split(delimiter1)
            l2 = file2_r[headers2 + i].split(delimiter2)

            if acc[i]:
                j = 0
                for ind in matches:
                    try:
                        d = abs(float(l1[ind[0]]) - float(l2[ind[1]]))

                        diff[i, j] = d

                        j += 1
                        # else:
                        #     diff[i, j] = 0
                    except ValueError:
                        pass
            else:
                j = 0
                for ind in matches:
                    try:
                        d = abs(float(l1[ind[0]]) - float(l2[ind[1]]))

                        diff[i, j] = 0

                        j += 1

                    except ValueError:
                        pass

        maxima = []
        # minima = []
        means = []
        maxi_ind = []
        # find maximal differences for each column
        for i in range(len(matches)):
            ind_max = np.unravel_index(np.argmax(diff[:, i], axis=None), diff[:, i].shape)
            # ind_min = np.unravel_index(np.argmin(diff[:, i], axis=None), diff[:, i].shape)

            maxi = diff[ind_max[0], i]
            # mini = diff[ind_min[0], i]
            mean = np.mean(diff[:, i])

            maxima.append(maxi)
            # minima.append(mini)
            means.append(mean)

            maxi_ind.append(ind_max)

        header_out1 = '{}'.format(delimiter)
        header_out2 = '{}'.format(delimiter)
        differences_max = 'max_diff{}'.format(delimiter)
        # differences_min = 'min_diff{}'.format(delimiter)
        differences_mean = 'mean_diff{}'.format(delimiter)
        differences_max_ind = 'max_number_of_line{}'.format(delimiter)

        # create output
        for i in range(len(matches)):
            header_out1 += a[matches[i][0]].strip() + delimiter  # add name of columns
            header_out2 += b[matches[i][1]].strip() + delimiter  # add units of columns

            differences_max += str('{:.5f}'.format(maxima[i])) + delimiter
            # differences_min += str('{:.5f}'.format(minima[i])) + delimiter
            differences_mean += str('{:.5f}'.format(means[i])) + delimiter
            differences_max_ind += str('{}'.format(maxi_ind[i][0] + 1)) + delimiter

        output = header_out1[:-1] + '\n'
        output += header_out2[:-1] + '\n'
        output += differences_max[:-1] + '\n'
        # output += differences_min[:-1] + '\n'
        output += differences_mean[:-1] + '\n'
        output += differences_max_ind[:-1]

        file1.close()
        file2.close()

        return output

    @staticmethod
    def get_headers_delimiter(path: str):
        # return count of headers and used delimiter

        # open and read file
        file = open(path, 'r')
        file_r = file.read().splitlines()

        # find delimiter
        if ';' in file_r[0]:
            delimiter = ';'

        if ',' in file_r[0]:
            delimiter = ','

        headers = 0
        # find count of headers
        for i in range(5):
            try:
                first = file_r[i].split(delimiter)[0]
            except:
                break
            try:
                float(first)
            except ValueError:
                headers += 1

        file.close()

        return headers, delimiter


# -----------------------------------------------------------------#
Compare(path1, path2, compare_path, delimiter)
# files_p = open(r'c:\Users\Jakub\Desktop\pecny\data\compare_path.txt', 'r')
# files = files_p.read().splitlines()
# files_p.close()
# for i in range(0, len(files), 3):
#     l_split = files[i].split()
#
#     path1 = l_split[0]
#     path2 = files[i + 1]
#     compare_path = files[i + 2]
#
#     if len(l_split) == 2:
#         exec(l_split[1])
#     print(path2)
#     Compare(path1, path2, compare_path, delimiter=',')
#     allan = False
# -----------------------------------------------------------------#
