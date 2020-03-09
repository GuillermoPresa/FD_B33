import scipy.io
mat = scipy.io.loadmat('matlab.mat')
print(mat)
print(mat.get('flightdata')[0][0][4][0][0][0])
