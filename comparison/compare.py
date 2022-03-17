import os
import numpy as np
from glob import glob
import sqlite3 as sql

# -----------------------------------------------------------------
# path1 = r'c:\Users\Jakub\Desktop\pecny\agdas_807zaloha'
# path2 = r'c:\Users\Jakub\Desktop\pecny\data\x153_glb\res1\Files'
# compare_path = r'c:\Users\Jakub\Desktop\pecny\data\x153_glb\res1'
# delimiter = ','

# path1 = r'c:\Users\Jakub\Desktop\pecny\data\566_glb\agdas_808'
# path2 = r'c:\Users\Jakub\Desktop\pecny\data\566_glb\res\Files'
# compare_path = r'c:\Users\Jakub\Desktop\pecny\data\566_glb\res'
# delimiter = ','

# path1 = r'c:\Users\Jakub\Desktop\pecny\data\x172_wetzel\files_wetzel'
# path2 = r'c:\Users\Jakub\Desktop\pecny\data\x172_wetzel\res1\Files'
# compare_path = r'c:\Users\Jakub\Desktop\pecny\data\x172_wetzel\res1'
# delimiter = ','

# Compare(path1, path2, compare_path, delimiter)
out = True
# -----------------------------------------------------------------

class Compare():

    def __init__(self, path1: str, path2: str, compare_path: str, delimiter: str):

        self.path1 = path1
        self.path2 = path2
        self.compare_path = compare_path
        self.delimiter = delimiter

        self.Compare_dir()

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

        # get accepted
        path_database = os.path.realpath(os.path.join(os.path.dirname(path2), '.', 'data.db'))
        path_drops_file_matlab = glob(self.path2 + '\*drops.csv')[0]
        accepted = Compare.get_accepted(path_database, path_drops_file_matlab, delimiter)

        # compare each file
        for i in range(len(files1)):

            file1 = files1[i].split('\\')[-1]

            for j in range(len(files2)):

                file2 = files2[j].split('\\')[-1]

                if file1 == file2:
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

    @staticmethod
    def match_columns(l1: list, l2: list):

        # matches - position of column from first file with corresponding column from second file
        matches = []

        # prohibited column to comparing
        prohibited = ['', 'Time', 'Campaign', 'Date', 'Accepted', 'Set', 'Drop', 'Year', 'Month', 'Day', 'Hour',
                      'Minute', 'Second']

        # used_columns is used for 'estimgrad' because there are columns with the same name
        # key is name of column and value is frequency
        used_columns = {}
        # find corresponding columns
        for i in range(len(l1)):

            # if columns was already used, increase frequency
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

                if l1[i].strip() == l2[j].strip() and l1[i].strip() not in prohibited:
                    matches.append([i, j])
                    break

        return matches

    @staticmethod
    def get_accepted(path_database: str, path_drops_file_matlab: str, delimiter: str):
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

        j = 0
        accepted = []
        for i in range(len(file_r)-1):

            a = file_r[i].split(delimiter)[-1]

            if a != 'Acc':

                acc_python = acc[j][0]

                if acc_python == 0 or a.strip() == '0':
                    accepted.append(False)
                else:
                    accepted.append(True)

                j += 1

        return accepted

    @staticmethod
    def compare2files(path1: str, path2: str, delimiter: str, delimiter1: str, delimiter2: str, headers1: int, headers2: int, acc=True):

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
            ind_min = np.unravel_index(np.argmin(diff[:, i], axis=None), diff[:, i].shape)

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

            first = file_r[i].split(delimiter)[0]

            try:
                float(first)
            except ValueError:
                headers += 1

        file.close()

        return headers, delimiter


# -----------------------------------------------------------------#
Compare(path1, path2, compare_path, delimiter)
# -----------------------------------------------------------------#
