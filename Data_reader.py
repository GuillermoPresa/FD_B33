import scipy.io
import numpy as np
mat = scipy.io.loadmat('/Users/frederik/GitHub/FD_B33/matlab.mat')
#print(mat)
#print(mat.get('flightdata')[0][0][1][0][0][0])
#print(np.concatenate((mat.get('flightdata')[0][0][1][0][0][0]), axis=None ))
flightdata = {'vane_AOA': np.concatenate((mat.get('flightdata')[0][0][0][0][0][0]), axis=None),
              'elevator_dte': np.concatenate((mat.get('flightdata')[0][0][1][0][0][0]), axis=None),
              'column_fe': np.concatenate((mat.get('flightdata')[0][0][2][0][0][0]), axis=None),
              'lh_engine_FMF':np.concatenate((mat.get('flightdata')[0][0][3][0][0][0]), axis=None),
              'rh_engine_FMF':np.concatenate((mat.get('flightdata')[0][0][4][0][0][0]), axis=None),
              'lh_engine_itt':np.concatenate((mat.get('flightdata')[0][0][5][0][0][0]), axis=None),
              'rh_engine_itt':np.concatenate((mat.get('flightdata')[0][0][6][0][0][0]), axis=None),
              'lh_engine_OP':np.concatenate((mat.get('flightdata')[0][0][7][0][0][0]), axis=None),
              'rh_engine_OP':np.concatenate((mat.get('flightdata')[0][0][8][0][0][0]), axis=None),
              'lh_engine_fan_N1':np.concatenate((mat.get('flightdata')[0][0][9][0][0][0]), axis=None),
              'lh_engine_turbine':np.concatenate((mat.get('flightdata')[0][0][10][0][0][0]), axis=None),
              'rh_engine_fan_N1':np.concatenate((mat.get('flightdata')[0][0][11][0][0][0]), axis=None),
              'rh_engine_turbine':np.concatenate((mat.get('flightdata')[0][0][12][0][0][0]), axis=None),
              'lh_engine_FU':np.concatenate((mat.get('flightdata')[0][0][13][0][0][0]), axis=None),
              'rh_engine_FU':np.concatenate((mat.get('flightdata')[0][0][14][0][0][0]), axis=None),
              'delta_a':np.concatenate((mat.get('flightdata')[0][0][15][0][0][0]), axis=None),
              'delta_e':np.concatenate((mat.get('flightdata')[0][0][16][0][0][0]), axis=None),
              'delta_r':np.concatenate((mat.get('flightdata')[0][0][17][0][0][0]), axis=None),
              'Gps_date':np.concatenate((mat.get('flightdata')[0][0][18][0][0][0]), axis=None),
              'Gps_utcSec':np.concatenate((mat.get('flightdata')[0][0][19][0][0][0]), axis=None),
              'Ahrs1_Roll':np.concatenate((mat.get('flightdata')[0][0][20][0][0][0]), axis=None),
              'Ahrs1_Pitch':np.concatenate((mat.get('flightdata')[0][0][21][0][0][0]), axis=None),
              'Fms1_trueHeading':np.concatenate((mat.get('flightdata')[0][0][21][0][0][0]), axis=None),
              'Gps_lat':np.concatenate((mat.get('flightdata')[0][0][22][0][0][0]), axis=None),
              'Gps_long':np.concatenate((mat.get('flightdata')[0][0][23][0][0][0]), axis=None),
              'Ahrs1_bRollRate':np.concatenate((mat.get('flightdata')[0][0][24][0][0][0]), axis=None)}

print(flightdata['vane_AOA'])



