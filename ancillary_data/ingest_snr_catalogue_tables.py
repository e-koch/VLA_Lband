
'''
Combine tables of recent publications with SNR catalogues.

Paths only working for travelling drive
'''

from astropy.table import Table, Column
from astropy.table import join as join_tab
import astropy.units as u
from astropy.coordinates import SkyCoord


mmt_overview = Table.read("/Volumes/Travel_Data/M33_2/MMT_SNR_catalogue_long18.txt", format='ascii')
mmt_fluxes = Table.read("/Volumes/Travel_Data/M33_2/MMT_SNR_catalogue_long18_fluxes.txt", format='ascii')

xmm_tab2 = Table.read("/Volumes/Travel_Data/M33_2/table2_m33_xmmsnrs.txt", format='ascii')
xmm_tab3 = Table.read("/Volumes/Travel_Data/M33_2/table3_m33_xmmsnrs.txt", format='ascii')

mmt_all = join_tab(mmt_overview, mmt_fluxes, keys='ID')

xmm_all = join_tab(xmm_tab2, xmm_tab3,
                   keys=list(set(xmm_tab2.colnames) & set(xmm_tab3.colnames)))

# Combine into a single coordinate axis defined as a skycoord

coord_strings = ["{0} {1} {2} +{3} {4} {5}"
                 .format(rah, ram, ras, decd, decm, decs) for rah, ram, ras, decd, decm, decs
                 in mmt_all['RAh', 'RAm', 'RAs', 'DEd', 'DEm', 'DEs']]

coords = SkyCoord(coord_strings, frame='icrs', unit=(u.hourangle, u.deg))

mmt_all["RA"] = Column([out_str.split(" ")[0] for out_str in coords.to_string('hmsdms')])
mmt_all["Dec"] = Column([out_str.split(" ")[1] for out_str in coords.to_string('hmsdms')])

del mmt_all['RAh'], mmt_all['RAm'], mmt_all['RAs'], mmt_all['DEd'], mmt_all['DEm'], mmt_all['DEs']

xmm_all.rename_column("XMMID", "XMM")

comb_table = join_tab(mmt_all, xmm_all, keys=['XMM'],
                      table_names=['MMT', 'XMM'])

mmt_all.write("/Volumes/Travel_Data/M33_2/MMT_SNR_catalogue_long18_combined.txt",
              format='ascii', overwrite=True)

xmm_all.write("/Volumes/Travel_Data/M33_2/XMM_SNR_catalogue_garofali17_combined.txt",
              format='ascii', overwrite=True)

comb_table.write("/Volumes/Travel_Data/M33_2/MMT_XMM_SNR_catalogue_combined.txt",
                 format='ascii', overwrite=True)
