import numpy as np
import scipy as sp
import math
	
atoms = []
coordinates=[]
resn = []
resid = []
	
#create a dummy plane		
dummy_CA = [0, 0, 0]
dummy_O = [ 0, 2.4, 0]
atoms.append('CA5')
coordinates.append([2.9, 2.4,0.0])
resn.append('RE1')
resid.append(1)
atoms.append('CA')
coordinates.append([2.9, 0.0,0.0])
resn.append('RE1')
resid.append(1)

length = 2
i =0
	
while i < len(length):	
	#create a dummy plane for first atom	
	if i = 0:
		atoms.append('O')
		resid.append(1)
		resn.append('RE1')
		coord_O = [dummy_O[1]+dummy_CA[0],dummy_O[0], dummy_O[2]]
		coordinates.append(coord_O)
	else:
		pass
		
print(coordinates)