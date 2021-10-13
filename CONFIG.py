
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
