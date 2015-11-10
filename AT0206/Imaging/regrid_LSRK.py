
'''
Regrid to LSRK
'''


mstransform(vis='M33_b_c.ms/', outputvis='M33_b_c_LSRK.ms', regridms=True,
            phasecenter="J2000 01h33m50.904 +30d39m35.79",
            restfreq="1420.40575177MHz", outframe='LSRK', combinespws=True)
