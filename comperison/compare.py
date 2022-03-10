import os
import numpy as np
from glob import glob

#-----------------------------------------------------------------
path1 = r'c:\Users\Jakub\Desktop\pecny\agdas_807zaloha'
path2 = r'c:\Users\Jakub\Desktop\pecny\data\x153_glb\res1\Files'
compare_path = r'c:\Users\Jakub\Desktop\pecny\data\x153_glb\res1'
delimiter = ','
# Compare(path1, path2, compare_path, delimiter)
#-----------------------------------------------------------------

class Compare():

    def __init__(self, path1, path2, compare_path, delimiter):

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

        # compare each file
        for i in range(len(files1)):

            file1 = files1[i].split('\\')[-1]

            for j in range(len(files2)):

                file2 = files2[j].split('\\')[-1]

                if file1 == file2:

                    headers1, delimiter1 = Compare.get_headers_delimiter(files1[i])

                    headers2, delimiter2 = Compare.get_headers_delimiter(files2[j])

                    comparision = Compare.compare2files(files1[i], files2[j], self.delimiter, delimiter1, delimiter2, headers1, headers2)

                    output_file = os.path.join(self.compare_path, 'Comp_' + file1)

                    f = open(output_file, 'w')
                    f.write(comparision)
                    f.close()

    @staticmethod
    def match_columns(l1, l2):

        # matches - position of column from first file with corresponding column from second file
        matches = []

        # prohibited column to comparing
        prohibited = ['', 'Time', 'Campaign', 'Date', 'Accepted', 'Set', 'Drop', 'Year', 'Month', 'Day', 'Hour', 'Minute', 'Second']

        # find corresponding columns
        for i in range(len(l1)):

            for j in range(len(l2)):

                if l1[i].strip() == l2[j].strip() and l1[i].strip() not in prohibited:

                    matches.append([i, j])
                    break

        return matches

    @staticmethod
    def compare2files(path1, path2, delimiter, delimiter1, delimiter2, headers1, headers2):

        # open files
        file1 = open(path1, 'r')
        file2 = open(path2, 'r')

        # read files
        file1_r = file1.read().splitlines()
        file2_r = file2.read().splitlines()

        # get headers
        a = file1_r[0].split(delimiter1)
        b = file2_r[0].split(delimiter2)

        # units
        c = file1_r[1].split(delimiter1)

        # find matches
        matches = Compare.match_columns(a, b)

        # l1_p = len(file1_r[header1 + 2].split(delimiter1))
        n = len(file1_r) - headers1

        diff = np.zeros((n, len(matches)))

        # find differences
        for i in range(n):

            l1 = file1_r[headers1 + i].split(delimiter1)
            l2 = file2_r[headers2 + i].split(delimiter2)

            j = 0
            for ind in matches:

                d = abs(float(l1[ind[0]]) - float(l2[ind[1]]))

                diff[i, j] = d

                j += 1


        maxima = []
        minima = []
        means = []
        # find maximal differences for each column
        for i in range(len(matches)):

            maxi = max(diff[:, i])
            mini = min(diff[:, i])
            mean = np.mean(diff[:, i])

            maxima.append(maxi)
            minima.append(mini)
            means.append(mean)


        header_out1 = ','
        header_out2 = ','
        differences_max = 'max_diff,'
        differences_min = 'min_diff,'
        differences_mean = 'mean_diff,'
        # create output
        for i in range(len(matches)):

            header_out1 += a[matches[i][0]].strip() + delimiter # add name of columns
            header_out2 += b[matches[i][1]].strip() + delimiter # add units of columns

            differences_max += str('{:.5f}'.format(maxima[i])) + delimiter
            differences_min += str('{:.5f}'.format(minima[i])) + delimiter
            differences_mean += str('{:.5f}'.format(means[i])) + delimiter

        output = header_out1[:-1] + '\n' + header_out2[:-1] + '\n' + differences_max[:-1] + '\n' + differences_min[:-1] + '\n' + differences_mean[:-1]

        file1.close()
        file2.close()

        return output

    @staticmethod
    def get_headers_delimiter(path):
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


#-----------------------------------------------------------------#
Compare(path1, path2, compare_path, delimiter)
#-----------------------------------------------------------------#
