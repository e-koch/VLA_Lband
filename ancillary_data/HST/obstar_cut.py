
'''
Cut out OB stars from the photometry files
'''

from paths import hst_phat_path
import pandas as pd

df_phot = pd.read_hdf(hst_phat_path("halfbrick_6filt_phot/14610_M33-B01_1.phot.hdf5"),
                      key='data')
n_gst = df_phot.filter(regex='_gst').sum(axis=1)
n_gst.head()

df_good = df[n_gst >= 4]
print('Number of stars with 4+ GST measurements:', df_good.shape[0])
# Make cuts
obstars = df_phot.query('(f336w_f475w < -0.5) & (f336w_vega < 18.12) & (f336w_f475w > -2) & ()')

obstars.to_hdf(hst_phat_path('halfbrick_6filt_phot/14610_M33-B01_1.phot.OBcut.hdf5',
                             no_check=False),
               key='data', mode='w')