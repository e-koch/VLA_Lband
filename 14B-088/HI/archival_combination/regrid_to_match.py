
'''
Regrid the data to match AT0206
'''

from tasks import mstransform

mstransform(vis='14B-088_HI_LSRK.ms.contsub',
            outputvis='14B-088_HI_LSRK_AT0206_regrid.ms.contsub',
            regridms=True,
            phasecenter="J2000 01h33m50.904 +30d39m35.79",
            restfreq="1420.40575177MHz", mode='frequency',
            nchan=255, start='1420.582kHz', width='6.103kHz')
