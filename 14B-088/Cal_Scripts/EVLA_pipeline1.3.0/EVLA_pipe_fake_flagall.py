
'''
On mixed setups, flagall is already run. This defines the variables set
during that script so it doesn't need to be run multiple times.
'''

logprint("Starting EVLA_pipe_fake_flagall.py",
         logfileout='logs/flagall.log')
time_list = runtiming('flagall', 'start')
QA2_flagall = 'Pass'

logprint("These value are fake! You should have run the actual flagall script"
         " already in the initial pipeline run! Check those logs for actual"
         " flagging fractions.")

start_total = 1.0
start_flagged = 0.0

init_on_source_vis = 1.0

afterzero_total = 1.0
afterzero_flagged = 0.0

zero_flagged = 0.0

aftershadow_total = 1.0
aftershadow_flagged = 0.0

shadow_flagged = 0.0

flagdata_list=[]
cmdreason_list=[]


frac_flagged_on_source1 = 0.0

logprint("Finished EVLA_pipe_fake_flagall.py", logfileout='logs/flagall.log')
logprint("QA2 score: "+QA2_flagall, logfileout='logs/flagall.log')
time_list = runtiming('flagall', 'end')

pipeline_save()
