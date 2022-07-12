import numpy as np
import os

def getFG5X(ps):
    FG5X = {
        'Lpar': 0.122,  # [m] parasitic wavelength
        'frmin': 15 * ps,  # first fringe
        'frmax': 940 * ps,  # end fringe
        'fmodf': 8333.355,  # [Hz] modulation frequency
        'frmaxss': 1040 * ps,
        'frminss': 1,
        'frmaxplot': 1040 * ps,
        # 'sens_tn': 1,
        # 'sens_tx': 150 * ps,
        # 'sensa_tn': 3 * ps,
        # 'sensa_tx': 80 * ps,
        # 'sensa_bn': 900 * ps,
        # 'sensa_bx': 1000 * ps,
        'nforfft': 4501,
        'ksmooth': 3,
        'sens_bx': 1040 * ps,
        'Lmin': 3,
        'Lmax': 16,
        'Lcable': 2,
        'Acable': 0.004,
        'Pcable': np.pi / 2,
        'valenv': 0.005
    }
    # FG5X['sens_bn'] = FG5X['frmaxss'] - FG5X['sens_tx']
    return FG5X


def getFG5(ps):
    FG5 = {
        'Lpar': 0.044,  # [m] parasitic wavelength
        'frmin': 30 * ps,  # first fringe
        'frmax': 629 * ps,  # end fringe
        'fmodf': 8333.251,  # [Hz] modulation frequency
        'frmaxplot': 650 * ps,
        # 'sens_tn': 1,
        # 'sens_tx': 100 * ps,
        # 'sensa_tn': 15 * ps,
        # 'sensa_tx': 60 * ps,
        # 'sensa_bn': 550 * ps,
        # 'sensa_bx': 640 * ps,
        'nforfft': 3293,
        'ksmooth': 3,
        'sens_bx': 650 * ps,
        'Lmin': 3,
        'Lmax': 16,
        'Lcable': 3.7,
        'Acable': 0.004,
        'Pcable': np.pi / 2,
        'frminss': 1,
        'frmaxss': 650 * ps,
        'valenv': 0.002
    }
    # FG5['sens_bn'] = FG5['frmaxplot'] - 150 * ps
    # FG5['frmaxss'] = FG5['frmax'] + 10
    # FG5['frminss'] = FG5['frmin'] - 10
    return FG5


matrDatabase = {
    'schema': '''CREATE TABLE results (
             n INTEGER,
             m0 REAL,
             Set1 INTEGER,
             Drop1 INTEGER,
             Date TEXT,
             mjd REAL,
             z0_withGR REAL,
             v0_withGR REAL,
             a_withGR REAL,
             b_withGR REAL,
             c_withGR REAL,
             d_withGR REAL,
             e_withGR REAL,
             f_withGR REAL,
             g0_Gr REAL,
             CorrToTop REAL,
             Tide REAL,
             Load REAL,
             Baro REAL,
             Polar REAL,
             gTopCor REAL,
             g0 REAL,
             EffHeight REAL,
             CorToEffHeight REAL,
             Gradient REAL,
             GradientLSTm0 REAL,
             std REAL,
             vgg REAL,
             ssres REAL,
             Accepted INTEGER,
             Res TEXT)''',
    'insert': '''INSERT INTO results ( n, m0, Set1, Drop1, Date, mjd, z0_withGR, v0_withGR, a_withGR, b_withGR,
    c_withGR, d_withGR, e_withGR, f_withGR, g0_Gr, CorrToTop, Tide, Load, Baro, Polar, gTopCor, g0, EffHeight,
    CorToEffHeight, Gradient, GradientLSTm0, std, vgg, ssres, Accepted, Res) values({},{},{},{},"{}",{},{},{},{},{},
    {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},"{}")''',

    'updateAcc': '''UPDATE results
                SET Accepted = 0
                WHERE Set1 = {} and Drop1 = {}''',

    'matlog': '''select Set1, substr(Date,0,5), substr(Date, 6,2), substr(Date, 9,2), round(avg(substr(Date, 12,2))), round(avg(substr(Date, 15,2))), round(avg(substr(Date, 18,2))), avg(gTopCor), count(*), avg(mjd)
            from results
            where Accepted = 1
            group by Set1'''

}

statistic = {
    'mean:vxv': '''with tab as (select avg(gTopCor) as mean, Set1, Drop1 from results where Accepted = 1 group by Set1),
tab1 as (select r.Set1, r.Drop1, tab.mean as mean, r.gTopCor as g from results as r join tab on (tab.Set1 = r.Set1 ) WHERE r.Accepted =1)
select Set1, Drop1, tab1.mean, sum((mean-g)*(mean-g)) as vxv, count(*) from tab1  group by tab1.Set1'''
}

logo = r'''
     /| |‾ ‾ ‾ ‾  |‾ ‾ ‾ \       /| |‾ ‾ ‾ ‾|
    / | |         |       \     / | |
   /  | |    _ _  |        |   /  | |
  /---| |       | |        |  /---|  ‾ ‾ ‾ ‾|
 /    | |       | |       /  /    |         |
/     | |_ _ _ _| |_ _ _ /  /     | |_ _ _ _|
'''

wel = '''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Welcome to the graphical user interface of Agdas
Agdas homepage:
RIGTC homepage:             www.vugtk.cz
G. O. Pecny homepage:       www.pecny.cz
'''

separator = '==============================='
logo_picture = os.path.dirname(os.path.realpath(__file__)) + r'\picture\logo.png'
picture_unchecked = os.path.dirname(os.path.realpath(__file__)) + r'\picture\unchecked.png'

headers = {
    'estim': 'Set {0} Drop {0} Date/Time {0} m0 {0} z0 {0} z0-std {0} v0 {0} v0-std {0} g0 {0} g0-std {0} a {0} a-std {0} b {0} b-std {0} c {0} c-std {0} d {0} d-std {0} e {0} e-std {0} f {0} f-std \n   {0}  {0}  {0} {0} mm {0} mm {0} mm.s-1 {0} mm.s-1 {0} nm.s-2 {0} nm.s-2 {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm   ',
    'drops': ' Set{0} Drop{0} Date{0} g"(t=0s){0} STD{0} TOD{0} Tide{0} Load{0} Baro{0} Polar{0} g(TOD){0} g(Ef.H){0} Ef.H1{0} c.EfH{0} Acc',
    'matlogsets': 'Campaign{0} Set{0} Year{0} Month{0} Day{0} Hour{0} Minute{0} Second{0} MJD{0} VGG_inp{0} g{0} g_std{0} STD-Start{0}STD-Final{0}Accepted{0} Top height{0} Pressure{0} VGG{0} T-stat',
    'allan': 'n{0}ALLAN1{0}STD1{0}ALLAN2{0}STD2{0}ALLAN3{0}STD3',
    'residuals_final': 'Fringe{0} z [m]{0}Time [s]{0}Time Top [s]{0}Value [nm]{0}Filtered value [nm]',
    'residuals_final1000': 'Fringe{0} z [m]{0} Time [s]{0} Time Top [s]{0} resid [nm]{0} Filtered resid [nm]',
    'residuals_sets': 'Fringe{0} Time [s]{0}Time Top [s]{0}Value [nm]',
    'spectrum': 'Frequency [Hz]{0}Avr res [nm]{0}Avr spec [nm]',
    'estim_grad': """   {0}    {0}Fit with gradient{0}{0}{0}{0}{0}{0}{0}{0}{0}Fit without gradient{0}{0}{0}{0}{0}{0}{0}{0}{0}Gradient estimation{0}{0}
                        Set{0}Drop{0}z0{0}v0{0}g0{0}a{0}b{0}c{0}d{0}a_par{0}b_par{0}z0{0}v0{0}g0{0}a{0}b{0}c{0}d{0}a_par{0}b_par{0}z0{0}v0{0}vgg{0}vgg-std
                        {0}    {0}nm{0}nm.s-1{0}nm.s-2{0}nm{0}nm{0}nm{0}nm{0}nm{0}nm{0}nm{0}nm.s-1{0}nm.s-2{0}nm{0}nm{0}nm{0}nm{0}nm{0}nm{0}nm{0}nm.s-1{0}nm.s-2/mm{0}nm.s-2/mm""",
    'resgradsum': 'Fringe{0}z [m]{0}Time [s]{0}Time TOP [s]{0}Resgradsum4_mean [nm]{0}Filtered mean [nm]',
    'vgg_per_sets0': 'Set {0}	vgg0 [uGal/cm] {0}	mvgg0 [uGal/cm] {0}	dg0 [uGal] {0}	mdg0 [uGal]',
    'effective_height_corr': 'Drop{0} EffHeight{0} CorToEffHeight',
    'matlog': """Campaign{0}Gravimeter-type{0}Gravimeter-SN{0}Sitename{0}Sitecode{0}Latitude{0}Longitude{0}Elevation{0}P_norm{0}baro admit.{0}VGG{0}WEO-IE{0}WEO-fmod{0}Clock-10MHz{0}L_parasit{0}Ksol{0}Ksae{0}Kdis{0}Kimp{0}kpar{0}L_TTLcable{0}Start Fringe{0}Final Fringe{0}x_pole{0}y_pole{0}Sets{0}Drops/Set{0}Year{0}Month{0}Day{0}Hour{0}Minute{0}MJD{0}Duration{0}Ave_P{0}dP(MAX-MIN){0}Ave_tide{0}Tide(MAX-MIN){0}Tstud{0}Drops_accept{0}Ef.h t0{0}Ef.h.TOD{0}STD{0}H.ef.ins{0}g@H.ef.ins{0}STD{0}STD-Start{0}STD-final{0}VGG_AG{0}STD
    {0}{0}{0}{0}{0}deg{0}deg{0}m{0}hPa{0}uGal/hPa{0}uGal/m{0}nm{0}Hz{0}Hz{0}m{0}{0}0/1{0}0/1{0}0/1{0}0/1{0}m{0}{0}{0}arcsec{0}arcsec{0}{0}{0}{0}{0}{0}{0}{0}{0}hour{0}hPa{0}hPa{0}uGal{0}uGal{0}%{0}{0}mm{0}mm{0}mm{0}m{0}uGal{0}uGal{0}uGal{0}uGal{0}uGal/cm{0}uGal/cm"""

}

round_line_ind = {
    'drops': [[3, 2], [4, 2], [5, 2], [6, 2], [7, 2], [8, 2], [9, 2], [10, 2], [11, 2], [12, 3], [13, 3]],
    'allan': [[i, 5] for i in range(1, 7)],
    'residuals_final': [[1, 8], [2, 5], [3, 5], [4, 6]],
    'residuals_sets': [[1, 5], [2, 5], [3, 5]],
    'residuals_final1000': [[1, 8], [2, 5], [3, 5], [4, 6], [5, 6]],
    'spectrum': [[0, 4], [1, 4], [2, 4]],
    'estim': [[2, 3], [3, 5], [4, 14], [5, 10], [6, 14], [7, 3], [8, 3], [9, 4], [10, 4], [11, 4], [12, 4], [13, 4],
              [14, 4], [15, 4], [16, 4]],
    'effHeightCorr_Graph': [[1, 5], [2, 5]],
    'estim_grad': [[2, 5], [3, 5], [4, 5], [5, 5], [6, 5], [7, 5], [8, 5], [9, 5], [10, 5], [11, 5], [12, 5], [13, 5],
                   [14, 5], [15, 5], [16, 5], [17, 5], [18, 5], [19, 5], [20, 5], [21, 5], [22, 5], [23, 5]],
    'resgradsum': [[1, 5], [2, 5], [3, 5], [4, 5], [5, 5]],
    'vgg_per_sets0': [[1, 5], [2, 5], [3, 5], [4, 5]],
    'matlogsets': [[8, 5], [10, 5], [11, 5], [12, 5], [13, 5], [15, 4], [17, 5], [18, 5]]
}
# increase index of floated number due to adding date/time to estim file
for i in range(len(round_line_ind['estim'])):
    round_line_ind['estim'][i] = [round_line_ind['estim'][i][0] + 1, round_line_ind['estim'][i][1]]

warning_window = {
    'import_data': 'Import data behind',
    'project': 'The project does not exist',
    'internet': 'Internet connection fail',
    'split_set': 'Choose count of sets',
    'pole_corr': 'Choose polar correction',
    'cannot_wrtite_file': 'Cannot write file due statistic is not computed',
    'comparison': 'Files successfully compared'
}

# colors = ['b', 'g', 'r', ]

SAE = [0, 0, 0.006, 0.10, 0.077, 0.92, 0.087, 1, 0.177, 1.65, 0.200, 2, 0.215, 2.20, 0.228, 2, 0.243, 1.70, 0.262, 2,
       0.277, 2.250, 0.317, 3, 0.322, 3.11, 0.335, 3.38, 0.362, 4, 0.377, 4.34, 0.397, 4.74]

tau = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800,
               900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500]

sh1 = """"BG" - Bulgarian
"CS" - Czech
"DA" - Danish
"DE" - German
"EL" - Greek
"EN-GB" - English
"ES" - Spanish
"ET" - Estonian
"FI" - Finnish
"FR" - French
"HU" - Hungarian
"ID" - Indonesian
"IT" - Italian
"JA" - Japanese
"LT" - Lithuanian
"LV" - Latvian
"NL" - Dutch
"PL" - Polish
"PT-BR" - Brazilian Portuguese
"RO" - Romanian
"RU" - Russian
"SK" - Slovak
"SL" - Slovenian
"SV" - Swedish
"TR" - Turkish
"ZH" - Chinese"""

languages = {}
for i in sh1.splitlines():
    s, l = i.split('- ')
    languages[s.split('"')[1]] = l
