import scipy.io
import numpy as np
# 'matlab.mat'
mat = scipy.io.loadmat('FTISxprt-20200306_flight1.mat')
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
              'Column_Se':np.concatenate((mat.get('flightdata')[0][0][9][0][0][0]), axis=None),
              'lh_engine_fan_N1':np.concatenate((mat.get('flightdata')[0][0][10][0][0][0]), axis=None),
              'lh_engine_turbine_N2':np.concatenate((mat.get('flightdata')[0][0][11][0][0][0]), axis=None),
              'rh_engine_fan_N1':np.concatenate((mat.get('flightdata')[0][0][12][0][0][0]), axis=None),
              'rh_engine_turbine_N2':np.concatenate((mat.get('flightdata')[0][0][13][0][0][0]), axis=None),
              'lh_engine_FU':np.concatenate((mat.get('flightdata')[0][0][14][0][0][0]), axis=None),
              'rh_engine_FU':np.concatenate((mat.get('flightdata')[0][0][15][0][0][0]), axis=None),
              'delta_a':np.concatenate((mat.get('flightdata')[0][0][16][0][0][0]), axis=None),
              'delta_e':np.concatenate((mat.get('flightdata')[0][0][17][0][0][0]), axis=None),
              'delta_r':np.concatenate((mat.get('flightdata')[0][0][18][0][0][0]), axis=None),
              'Gps_date':np.concatenate((mat.get('flightdata')[0][0][19][0][0][0]), axis=None),
              'Gps_utcSec':np.concatenate((mat.get('flightdata')[0][0][20][0][0][0]), axis=None),
              'Ahrs1_Roll':np.concatenate((mat.get('flightdata')[0][0][21][0][0][0]), axis=None),
              'Ahrs1_Pitch':np.concatenate((mat.get('flightdata')[0][0][22][0][0][0]), axis=None),
              'Fms1_trueHeading':np.concatenate((mat.get('flightdata')[0][0][23][0][0][0]), axis=None),
              'Gps_lat':np.concatenate((mat.get('flightdata')[0][0][24][0][0][0]), axis=None),
              'Gps_long':np.concatenate((mat.get('flightdata')[0][0][25][0][0][0]), axis=None),
              'Ahrs1_bRollRate':np.concatenate((mat.get('flightdata')[0][0][26][0][0][0]), axis=None),
              'Ahrs1_bPitchRate':np.concatenate((mat.get('flightdata')[0][0][27][0][0][0]), axis=None),
              'Ahrs1_bYawRate':np.concatenate((mat.get('flightdata')[0][0][28][0][0][0]), axis=None),
              'Ahrs1_bLongAcc':np.concatenate((mat.get('flightdata')[0][0][29][0][0][0]), axis=None),
              'Ahrs1_bLatAcc':np.concatenate((mat.get('flightdata')[0][0][30][0][0][0]), axis=None),
              'Ahrs1_bNormAcc':np.concatenate((mat.get('flightdata')[0][0][31][0][0][0]), axis=None),
              'Ahrs1_aHdgAcc':np.concatenate((mat.get('flightdata')[0][0][32][0][0][0]), axis=None),
              'Ahrs1_xHdgAcc':np.concatenate((mat.get('flightdata')[0][0][33][0][0][0]), axis=None),
              'Ahrs1_VertAcc':np.concatenate((mat.get('flightdata')[0][0][34][0][0][0]), axis=None),
              'Dadc1_sat':np.concatenate((mat.get('flightdata')[0][0][35][0][0][0]), axis=None),
              'Dadc1_tat':np.concatenate((mat.get('flightdata')[0][0][36][0][0][0]), axis=None),
              'Dadc1_alt':np.concatenate((mat.get('flightdata')[0][0][37][0][0][0]), axis=None),
              'Dadc1_bcAlt':np.concatenate((mat.get('flightdata')[0][0][38][0][0][0]), axis=None),
              'Dadc1_bcAltMb':np.concatenate((mat.get('flightdata')[0][0][39][0][0][0]), axis=None),
              'Dadc1_mach':np.concatenate((mat.get('flightdata')[0][0][40][0][0][0]), axis=None),
              'Dadc1_cas':np.concatenate((mat.get('flightdata')[0][0][41][0][0][0]), axis=None),
              'Dadc1_tas':np.concatenate((mat.get('flightdata')[0][0][42][0][0][0]), axis=None),
              'Dadc1_altRate':np.concatenate((mat.get('flightdata')[0][0][43][0][0][0]), axis=None),
              'measurement_running':np.concatenate((mat.get('flightdata')[0][0][44][0][0][0]), axis=None),
              'measurement_n_rdy':np.concatenate((mat.get('flightdata')[0][0][45][0][0][0]), axis=None),
              'display_graph_state':np.concatenate((mat.get('flightdata')[0][0][46][0][0][0]), axis=None),
              'display_active_screen':np.concatenate((mat.get('flightdata')[0][0][47][0][0][0]), axis=None),
              'time':np.concatenate((mat.get('flightdata')[0][0][48][0][0][0]), axis=None)}


#print(flightdata['Dadc1_alt'])





