
'''
Cut out OB stars from the photometry files
'''

from paths import hst_phat_path
import pandas as pd

df_phot = pd.read_hdf(hst_phat_path("halfbrick_6filt_phot/14610_M33-B01_1.phot.hdf5"),
                      key='data')
n_gst = df_phot.filter(regex='_gst').sum(axis=1)
n_gst.head()

df_good = df_phot[n_gst >= 4]
print('Number of stars with 4+ GST measurements:', df_good.shape[0])

# Define some colour columns
df_good = df_good.assign(f336w_f475w_vega = df_good.f336w_vega - df_good.f475w_vega)

# Make cuts
ostars = df_good.query('(f336w_f475w_vega < -0.5) & (f336w_vega < 18.12) & (f336w_f475w_vega > -2)')

ostars.to_hdf(hst_phat_path('halfbrick_6filt_phot/14610_M33-B01_1.phot.Ocut.hdf5',
                            no_check=True),
              key='data', mode='w')