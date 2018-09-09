
'''
Re-weight the 14B data prior to imaging with the 17B data
'''

from tasks import statwt, initweights

# Calculate weights per observations
# for i in range(12):

#     statwt(vis='14B-088_HI_LSRK.ms.contsub', dorms=False,
#            field='M33*', fitspw="0:0~350;500~600;1750~2000",
#            minsamp=2, datacolumn='data', observation=str(i))


# # And on the 17B data

# for i in range(17):

#     print("On 17B Obs {}".format(i))

#     statwt(vis='17B-162_HI_spw_0_LSRK.mms.contsub', dorms=False,
#            field='M33*', fitspw="0:1043~1393;1543~1643;2793~3043",
#            minsamp=2, datacolumn='data', observation=str(i))

# There remains a much larger weighting in the C config data when using
# statwt which just doesn't seem right...
# Put everything on the same scale by just initializing the weights to match

initweights(vis='14B-088_HI_LSRK.ms.contsub', wtmode='nyq')
initweights(vis='17B-162_HI_spw_0_LSRK.ms.contsub', wtmode='nyq')

initweights(vis='14B-088_HI_LSRK.ms', wtmode='nyq')
initweights(vis='17B-162_HI_spw_0_LSRK.ms', wtmode='nyq')
