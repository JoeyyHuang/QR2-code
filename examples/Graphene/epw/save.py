#
# Post-processing script QE --> EPW 
# 14/07/2015 - Samuel Ponce
#

import numpy as np
import os

# Enter the number of irr. q-points 
user_input = input('Enter the prefix used for PH calculations (e.g. diam)\n')
prefix = str(user_input)

# Enter the number of irr. q-points 
user_input = input('Enter the number of irreducible q-points\n')
nqpt = user_input
try:
  nqpt = int(user_input)
except ValueError:
  raise Exception('The value you enter is not an integer!')

os.system('mkdir save')

for iqpt in np.arange(1,nqpt+1):
  label = str(iqpt)

  os.system('cp ../phonon/dyn'+str(iqpt)+' save/'+prefix+'.dyn_q'+label)
  if (iqpt == 1):
    os.system('cp ../phonon/tmp/_ph0/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
    os.system('cp -r ../phonon/tmp/_ph0/'+prefix+'.phsave save/')
    # if PAW
    # os.system('cp ../phonon/tmp/_ph0/'+prefix+'.dvscf_paw1 save/'+prefix+'.dvscf_paw_q'+label)
  else:
    os.system('cp ../phonon/tmp/_ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
    # if PAW
    # os.system('cp ../phonon/tmp/_ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf_paw1 save/'+prefix+'.dvscf_paw_q'+label)
