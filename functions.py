import json
import numpy as np
from math import floor, trunc
import math, requests, os


def allan(data: list, tau: list):
    """
    Allan standard deviation.

    Parameters
    ----------
    data : list
        list with data [x1,...xn]
    tau : list
        list with numbers of intervals

    Returns
    -------
    res : list
        [[allan, std, count], ...]
    """

    res = []
    for f in tau:
        k = 0
        int = []
        s = 0
        it = 0
        while k < len(data):
            se = data[k:k + f]

            if len(se) == f:
                int.append(np.mean(se))

                if it > 0:
                    s += (int[-1] - int[-2]) ** 2

            k += f
            it += 1

        if len(int) - 1 > 0:
            v = np.sqrt(1 / (2 * (len(int) - 1)) * s)
            err = v / np.sqrt(floor(len(data) / f))
            res.append([v, err, len(int) - 1])

    return res


def printDict(dict: dict):
    """
    Just print dictionary in json format in form:
    {
    key_1 : value_1
    .
    .
    key_n : value_n
    }

    Parameters
    ----------
    dict : dict

    """
    print(json.dumps(dict, indent=5))


def roundList(l: list, index: list):
    """
    Round list for printing into text files.

    Parameters
    ----------
    l : list
        line for rounding
    index : list
        information for rounding, [[position1, number_of_decimal_places1],[position2, number_of_decimal_places2], ...]

    Returns
    -------
    l : list
        rounded list
    """

    for j in index:
        a = '{:.' + str(j[1]) + 'f}'
        if l[j[0]] == '-':
            l[j[0]] = '-'
        else:
            l[j[0]] = a.format(l[j[0]])

    return l


def date_to_mjd(year, month, day):
    """
    :Author: Matt Davis
    :Website: http://github.com/jiffyclub

    Change to returning modified Julian date

    Convert a date to Julian Day.

    Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet',
        4th ed., Duffet-Smith and Zwart, 2011.

    Parameters
    ----------
    year : int
        Year as integer. Years preceding 1 A.D. should be 0 or negative.
        The year before 1 A.D. is 0, 10 B.C. is year -9.

    month : int
        Month as integer, Jan = 1, Feb. = 2, etc.

    day : float
        Day, may contain fractional part.

    Returns
    -------
    mjd : float
        modified Julian Day

    Examples
    --------
    Convert 6 a.m., February 17, 1985 to Julian Day

    >>> date_to_mjd(1985,2,17.25)
    2446113.75

    """
    if month == 1 or month == 2:
        yearp = year - 1
        monthp = month + 12
    else:
        yearp = year
        monthp = month

    # this checks where we are in relation to October 15, 1582, the beginning
    # of the Gregorian calendar.
    if ((year < 1582) or
            (year == 1582 and month < 10) or
            (year == 1582 and month == 10 and day < 15)):
        # before start of Gregorian calendar
        B = 0
    else:
        # after start of Gregorian calendar
        A = trunc(yearp / 100.)
        B = 2 - A + trunc(A / 4.)

    if yearp < 0:
        C = trunc((365.25 * yearp) - 0.75)
    else:
        C = trunc(365.25 * yearp)

    D = trunc(30.6001 * (monthp + 1))

    jd = B + C + D + day + 1720994.5

    return jd - 2400000.5


def mjd_to_jd(mjd):
    """
    :Author: Matt Davis
    :Website: http://github.com/jiffyclub
    Convert Modified Julian Day to Julian Day.

    Parameters
    ----------
    mjd : float
        Modified Julian Day

    Returns
    -------
    jd : float
        Julian Day


    """
    return mjd + 2400000.5


def jd_to_date(jd):
    """
    :Author: Matt Davis
    :Website: http://github.com/jiffyclub
    Convert Julian Day to date.

    Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet',
        4th ed., Duffet-Smith and Zwart, 2011.

    Parameters
    ----------
    jd : float
        Julian Day

    Returns
    -------
    year : int
        Year as integer. Years preceding 1 A.D. should be 0 or negative.
        The year before 1 A.D. is 0, 10 B.C. is year -9.

    month : int
        Month as integer, Jan = 1, Feb. = 2, etc.

    day : float
        Day, may contain fractional part.

    Examples
    --------
    Convert Julian Day 2446113.75 to year, month, and day.

    >>> jd_to_date(2446113.75)
    (1985, 2, 17.25)

    """
    jd = jd + 0.5

    F, I = math.modf(jd)
    I = int(I)

    A = math.trunc((I - 1867216.25) / 36524.25)

    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I

    C = B + 1524

    D = math.trunc((C - 122.1) / 365.25)

    E = math.trunc(365.25 * D)

    G = math.trunc((C - E) / 30.6001)

    day = C - E + F - math.trunc(30.6001 * G)

    if G < 13.5:
        month = G - 1
    else:
        month = G - 13

    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715

    hour = (day - math.floor(day)) * 24
    minut = (hour - math.floor(hour)) * 60
    sec = (minut - math.floor(minut)) * 60

    return year, month, math.floor(day), math.floor(hour), math.floor(minut), math.floor(sec)


def rssq(x):
    """
    Root sum of squares.

    Parameters
    ----------
    x : np.array

    Returns
    -------
    res : list

    """
    res = []

    for i in range(x.shape[0]):

        square_sum = 0
        for j in range(x.shape[1]):
            square_sum += x[i, j] * x[i, j]

        # print(np.sqrt(square_sum))
        res.append(np.sqrt(square_sum))

    return res


def movingAverage(x: list, n=50):
    """
    Return moving average with floating window.

    Parameters
    ----------
    x : list
        list with data
    n : int
        range for moving average (kernel), with n = 50 will be average from 101 values

    Returns
    -------
    res : list
        moving averages
    plot_range : list
        range for plotting

    """

    res = []
    plot_range = []
    for i in range(n, len(x) - n):
        m = np.mean(x[i - n:i + n + 1])
        res.append(m)
        plot_range.append(i)

    return res, plot_range


def get_date_last_version():
    """
    Check date of last commit.

    Returns
    -------
    date of last commit, 0 in case that connection failed

    """

    repo = 'jakubsimek97/Agdas'
    try:
        r = requests.get('https://api.github.com/repos/{0}/commits?per_page=1'.format(repo))
        return r.json()[0]['commit']['author']['date']
    except requests.exceptions.ConnectionError:
        return 0


def write_last_version(version_date: str):
    """
    Write last version of repository to last_version.txt

    Parameters
    ----------
    version_date : str
        last version date

    """

    script_path = os.path.dirname(os.path.realpath(__file__))
    path = os.path.join(script_path, 'last_version.txt')
    f = open(path, 'w')
    f.write(version_date)
    f.close()

def read_last_version():
    """
    Return date of last version from last_version.txt.

    Returns
    -------
    r : str
        date of last version

    """
    script_path = os.path.dirname(os.path.realpath(__file__))
    path = os.path.join(script_path, 'last_version.txt')
    f = open(path, 'r')
    r = f.read()
    f.close()
    return r
