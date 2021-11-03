
def getFG5X(ps):

    FG5X={
    'Lpar' : 0.122*1e9, # [nm] parasitic wavelength
    'frmin' : 15*ps, # first fringe
    'frmax' : 940*ps, # end fringe
    'fmodf' : 8333.355, # [Hz] modulation frequency
    'frmaxss' : 1040*ps,
    'frminss' : 1,
    'frmaxplot' : 1040*ps,
    'sens_tn' : 1,
    'sens_tx' : 150*ps,
    'sensa_tn' : 3*ps,
    'sensa_tx' : 80*ps,
    'sensa_bn' : 900*ps,
    'sensa_bx' : 1000*ps,
    'nforfft' : 4501,
    'ksmooth' : 3
    }
    FG5X['sens_bn'] = FG5X['frmaxss'] - FG5X['sens_tx']
    return FG5X


matrDatabase={
    'schema' : '''CREATE TABLE results (
             Set1 INTEGER,
             Drop1 INTEGER,
             Date TEXT,
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
             std REAL,
             vgg REAL,
             Accepted INTEGER,
             Res TEXT)''',
    'insert' : '''INSERT INTO results (
            Set1,
            Drop1,
            Date,
            g0_Gr,
            CorrToTop,
            Tide,
            Load,
            Baro,
            Polar,
            gTopCor,
            g0,
            EffHeight,
            CorToEffHeight,
            Gradient,
            std,
            vgg,
            Accepted,
            Res) values({},{},"{}",{},{},{},{},{},{},{},{},{},{},{},{},{},{},"{}")''',

    'updateAcc' : '''UPDATE results
                SET Accepted = 0
                WHERE Set1 = {} and Drop1 = {}''',

    'matlog' : '''select Set1, substr(Date,0,5), substr(Date, 6,2), substr(Date, 9,2), round(avg(substr(Date, 12,2))), round(avg(substr(Date, 15,2))), round(avg(substr(Date, 18,2))), avg(gTopCor), count(*)
            from results
            where Accepted = 1
            group by Set1'''

}

statistic={
    'mean:vxv':'''with tab as (select avg(gTopCor) as mean, Set1, Drop1 from results group by Set1),
tab1 as (select r.Set1, r.Drop1, tab.mean as mean, r.gTopCor as g from results as r join tab on (tab.Set1 = r.Set1 ) WHERE r.Accepted =1)
select Set1, Drop1, tab1.mean, sum((mean-g)*(mean-g)) as vxv from tab1 group by tab1.Set1'''
}

logo=r'''
     /| |‾ ‾ ‾ ‾  |‾ ‾ ‾ \       /| |‾ ‾ ‾ ‾|
    / | |         |       \     / | |
   /  | |    _ _  |        |   /  | |
  /---| |       | |        |  /---|  ‾ ‾ ‾ ‾|
 /    | |       | |       /  /    |         |
/     | |_ _ _ _| |_ _ _ /  /     | |_ _ _ _|
'''

wel='''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Welcome to grafical user interface of Agdas
Agdas homepage:
RIGTC homepage:             www.vugtk.cz
G. O. Pecny homepage:       www.pecny.cz
'''

separator='================================='
logo_picture = 'picture/logo.png'


headers={
    'estim' : 'Set {0} Drop {0} m0 {0} z0 {0} z0-std {0} v0 {0} v0-std {0} g0 {0} g0-std {0} a {0} a-std {0} b {0} b-std {0} c {0} c-std {0} d {0} d-std {0} e {0} e-std {0} f {0} f-std \n   {0}    {0} {0} mm {0} mm {0} mm.s-1 {0} mm.s-1 {0} nm.s-2 {0} nm.s-2 {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm {0} nm   ',
    'drops' : ' Set{0} Drop{0} Date{0} g"(t=0s){0} STD{0} TOD{0} Tide{0} Load{0} Baro{0} Polar{0} g(TOD){0} g(Ef.H){0} Ef.H1{0} c.EfH{0} Acc',
    'matlog' : 'Campaign{0} Set{0} Year{0} Month{0} Day{0} Hour{0} Minute{0} Second{0} MJD{0} VGG_inp{0} g{0} g_std{0} STD-Start{0}STD-Final{0}Accepted{0} Top height{0} Pressure{0} VGG{0} T-stat',
    'allan' : 'n{0}ALLAN1{0}STD1{0}ALLAN2{0}STD2{0}ALLAN3{0}STD3',
    'residuals_final' : 'Fringe{0} z [m]{0}Time [s]{0}Time Top [s]{0}Value [nm]{0}Filtered value [nm]',
    'residuals_sets' : 'Fringe{0} Time [s]{0}Time Top [s]{0}Value [nm]',
    'spectrum' : 'Frequency [Hz]{0}Avr res [nm]{0}Avr spec [nm]',
    'estim_grad' : """   {0}    {0}Fit with gradient{0}{0}{0}{0}{0}{0}{0}{0}{0}Fit without gradient{0}{0}{0}{0}{0}{0}{0}{0}{0}Gradient estimation{0}{0}
                        Set{0}Drop{0}z0{0}v0{0}g0{0}a{0}b{0}c{0}d{0}a_par{0}b_par{0}z0{0}v0{0}g0{0}a{0}b{0}c{0}d{0}a_par{0}b_par{0}z0{0}v0{0}vgg{0}vgg-std
                        {0}    {0}nm{0}nm.s-1{0}nm.s-2{0}nm{0}nm{0}nm{0}nm{0}nm{0}nm{0}nm{0}nm.s-1{0}nm.s-2{0}nm{0}nm{0}nm{0}nm{0}nm{0}nm{0}nm{0}nm.s-1{0}nm.s-2/mm{0}nm.s-2/mm"""
}

round_line_ind={
    'drops' : [[3,2],[4,2],[5,2],[6,2],[7,2],[8,2],[9,2],[10,2],[11,2],[12,3],[13,3]],
    'allan' : [[i, 5] for i in range(1,7)],
    'residuals_final' : [[1,8], [2, 5], [3, 5], [4, 6]],
    'residuals_sets' : [[1,5],[2,5],[3,5]],
    'spectrum' : [[1,4],[2,4],[3,4]],
    'estim' : [[2,3], [3,4], [4,6], [5,4], [6,6], [7,3], [8,3], [9,4], [10,4], [11,4], [12,4], [13,4], [14,4], [15,4], [16,4]],
    'effHeightCorr_Graph' : [[1,5],[2,5]],
    'estim_grad' : [[2, 5], [3, 5], [4, 5], [5, 5], [6, 5], [7, 5], [8, 5], [9, 5], [10, 5], [11, 5], [12, 5], [13, 5], [14, 5], [15, 5], [16, 5], [17, 5], [18, 5], [19, 5], [20, 5], [21, 5], [22, 5], [23, 5]]
}
