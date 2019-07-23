writePdbAA = {'A' : {
					'Res_name': 'ALA', 
					'atoms':['N','CA','C','O','CB'],
					'atom_type':['N','C','C','O','C']},
			  'Ala' : {
					'Res_name': 'ALA', 
					'atoms':['N','CA','C','O','CB'],
					'atom_type':['N','C','C','O','C']},
					
			  'R' : {
					'Res_name': 'ARG', 
					'atoms':['N','CA','C','O','CB','CG','CD','NE','CZ','NH1','NH2'],
					'atom_type':['N','C','C','O','C','C','C','N','C','N','N']},		
			
		      'N' : {
					'Res_name': 'ASN', 
					'atoms':['N','CA','C','O','CB','CG','OD1','ND2'],
					'atom_type':['N','C','C','O','C','C','O','N']},
					
			  'D' : {
					'Res_name': 'ASP', 
					'atoms':['N','CA','C','O','CB','CG','OD1','OD2'],
					'atom_type':['N','C','C','O','C','C','O','O']},
			 
			  'C' : {
					'Res_name': 'CYS', 
					'atoms':['N','CA','C','O','CB','SG'],
					'atom_type':['N','C','C','O','C','S']},
					
			  'E' : {
					'Res_name': 'GLU', 
					'atoms':['N','CA','C','O','CB','CG','CD','OE1','OE2'],
					'atom_type':['N','C','C','O','C','C','C','O','O']},
			 
			  'Q' : {
					'Res_name': 'GLN', 
					'atoms':['N','CA','C','O','CB','CG','CD','OE1','NE2'],
					'atom_type':['N','C','C','O','C','C','C','O','N']},
					
			  'G' : {
					'Res_name': 'GLY', 
					'atoms':['N','CA','C','O'],
					'atom_type':['N','C','C','O']},
					
			  'H' : {
					'Res_name': 'HIS', 
					'atoms':['N','CA','C','O','CB','CG','CD2','ND1','CE1','NE2'],
					'atom_type':['N','C','C','O','C','C','C','N','C','N']},
					
			  'I' : {
					'Res_name': 'ILE', 
					'atoms':['N','CA','C','O','CB','CG1','CG2','CD1'],
					'atom_type':['N','C','C','O','C','C','C','C']},
					
			  'L' : {
					'Res_name': 'LEU', 
					'atoms':['N','CA','C','O','CB','CG','CD1','CD2'],
					'atom_type':['N','C','C','O','C','C','C','C']},
				
			  'K' : {
					'Res_name': 'LYS', 
					'atoms':['N','CA','C','O','CB','CG','CD','CE','NZ'],
					'atom_type':['N','C','C','O','C','C','C','C','N']},
					
			  'M' : {
					'Res_name': 'MET', 
					'atoms':['N','CA','C','O','CB','CG','SD','CE'],
					'atom_type':['N','C','C','O','C','C','S','C']},
			
			  'F' : {
					'Res_name': 'PHE', 
					'atoms':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C']},
				
			  'P' : {
					'Res_name': 'PRO', 
					'atoms':['N','CA','C','O','CB','CG','CD'],
					'atom_type':['N','C','C','O','C','C','C']},
			 
			  'S' : {
					'Res_name': 'SER', 
					'atoms':['N','CA','C','O','CB','OG'],
					'atom_type':['N','C','C','O','C','O']},
			 
			  'T' : {
					'Res_name': 'THR', 
					'atoms':['N','CA','C','O','CB', 'OG1','CG2'],
					'atom_type':['N','C','C','O','C','O','C']},
					
			  'W' : {
					'Res_name': 'TRP', 
					'atoms':['N','CA','C','O','CB','CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2'],
					'atom_type':['N','C','C','O','C','C','C','C','N','C','C','C','C','C']},
			
			  'Y' : {
					'Res_name': 'TYR', 
					'atoms':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ', 'OH'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','O']},
			  
			  'V' : {
					'Res_name': 'VAL', 
					'atoms':['N','CA','C','O','CB','CG1','CG2'],
					'atom_type':['N','C','C','O','C','C','C']},
					
##unnatural amino acids
	
			  'Pheiodo' : {
					'Res_name': 'PHI', 
					'atoms':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ', 'I'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','I'],
					'link_rec':[(0,1),(1,2), (2,3), (3,4,4), (2,5), (5,6),(6,7), (6,8,8), (7,9,9),(8,10),(9,11), (10,11,11),(11,12)]},
			  'Phebromo' : {
					'Res_name': 'PBR', 
					'atoms':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ', 'BR'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','Br'],
					'link_rec':[(1,2), (2,3), (3,4,4), (2,5), (5,6),(6,7), (6,8,8), (7,9,9),(8,10),(9,11), (10,11,11),(11,12)]},
			 
			  'Phefluoro' : {
					'Res_name': 'P5F', 
					'atoms':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ', 'F1','F2','F3','F4','F5'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','F','F','F','F','F','F'],
					'link_rec':[(1,2), (2,3), (3,4,4), (2,5), (5,6),(6,7), (6,8,8), (7,9,9),(8,10),(9,11), (10,11,11),(7,13),(8,12),(9,15),(10,14),(11,16)]},
		      'Phe(I)' : {
					'Res_name': 'PHI', 
					'atoms':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ', 'I'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','I'],
					'link_rec':[(1,2), (2,3), (3,4,4), (2,5), (5,6),(6,7), (6,8,8), (7,9,9),(8,10),(9,11), (10,11,11),(11,12)]},
			  'Phe(Br)' : {
					'Res_name': 'PBR', 
					'atoms':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ', 'BR'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','Br'],
					'link_rec':[(1,2), (2,3), (3,4,4), (2,5), (5,6),(6,7), (6,8,8), (7,9,9),(8,10),(9,11), (10,11,11),(11,12)]},
			 
			  'Phe(5F)' : {
					'Res_name': 'P5F', 
					'atoms':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ', 'F1','F2','F3','F4','F5'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','F','F','F','F','F','F'],
					'link_rec':[(1,2), (2,3), (3,4,4), (2,5), (5,6),(6,7), (6,8,8), (7,9,9),(8,10),(9,11), (10,11,11),(7,13),(8,12),(9,15),(10,14),(11,16)]},
				
			  'O' : {
					'Res_name': 'ORN', 
					'atoms':['N','CA','C','O','CB','CG','CD','NE'],
					'atom_type':['N','C','C','O','C','C','C','N+1'],
					'link_rec': [(1,2), (2,3), (3,4,4),(2,5),(5,6),(6,7),(7,8)]},
					
			  'Nle' : {
					'Res_name': 'NLE', 
					'atoms':['N','CA','C','O','CB','CG','CD','CE'],
					'atom_type':['N','C','C','O','C','C','C','C'],
					'link_rec': [(1,2), (2,3), (3,4,4),(2,5),(5,6),(6,7),(7,8)]},
			
   			  'Nva' : {
					'Res_name': 'NVA', 
					'atoms':['N','CA','C','O','CB','CG','CD'],
					'atom_type':['N','C','C','O','C','C','C'],
					'link_rec': [(1,2), (2,3), (3,4,4),(2,5),(5,6),(6,7)]},
			  'Pra' : {
					'Res_name': 'PRA', 
					'atoms':['N','CA','C','O','CB','CG','CD'],
					'atom_type':['N','C','C','O','C','C','C'],
					'link_rec': [(1,2), (2,3), (3,4,4),(2,5),(5,6),(6,7,7),(6,7)]},
			  'Cha' : {
					'Res_name': 'CHA', 
					'atoms':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (2,5), (5,6),(6,7), (6,8), (7,9),(8,10),(9,11), (10,11)]},
			
			  'Hao' : {
					'Res_name': 'HAO', 
					'atoms':['N','N02','C','O','C08','O06','C06','C05','OM','CM','C04','C03','C02','C07','N01','C01','O01'],
					'atom_type':['N','N','C','O','C','O','C','C','O','C','C','C','C','C','N','C','O'],
					'link_rec': [(1,2),(2,5),(5,6,6), (5,7),(7,8,8), (7,14),(8,9), (9,10), (8, 11), (11,12,12), (12,13), (13,14,14), (13, 15), (15,16), (16,17,17),(16,3),(3,4,4) ]},
					
			  
##Misc
			  'Ac':{
					'Res_name': 'ACE', 
					'atoms':['O','C','CH3'],
					'atom_type':['O','C','C']},
##NMe Amino Acids					
			   
			  'Ala(NMe)' : {
					'Res_name': 'MAA', 
					'atoms':['N','CA','C','O','CB','CN'],
					'atom_type':['N','C','C','O','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (2,5), (1,6)]},
			
			  'Gly(NMe)' : {
					'Res_name': 'SAR', 
					'atoms':['N','CA','C','O','CN'],
					'atom_type':['N','C','C','O','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5)]},
			  'Ile(NMe)' : {
					'Res_name': 'IML', 
					'atoms':['N','CA','C','O','CB','CG1','CG2','CD1','CN'],
					'atom_type':['N','C','C','O','C','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (2,5), (5,6), (5,7), (6,8), (1,9)]},
					
			  'Leu(NMe)' : {
					'Res_name': 'MLE', 
					'atoms':['N','CA','C','O','CB','CG','CD1','CD2','CN'],
					'atom_type':['N','C','C','O','C','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (2,5), (5,6), (6,7), (6,8), (1,9)]},
			
			  'Nle(NMe)' : {
					'Res_name': 'MNL', 
					'atoms':['N','CA','C','O','CB','CG','CD','CE','CN'],
					'atom_type':['N','C','C','O','C','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (2,5), (5,6), (6,7), (7,8), (1,9)]},
					
			  'Nva(NMe)' : {
					'Res_name': 'MNV', 
					'atoms':['N','CA','C','O','CB','CG','CD','CN'],
					'atom_type':['N','C','C','O','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (2,5), (5,6), (6,7), (1,8)]},
					
			  'Phe(NMe)' : {
					'Res_name': 'MEA', 
					'atoms':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ','CN'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (2,5), (5,6),(6,7), (6,8,8), (7,9,9),(8,10),(9,11), (10,11,11),(1,12)]},
					
			  'Tyr(NMe)' : {
					'Res_name': 'YNM', 
					'atoms':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ', 'OH', 'CN'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','O','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (2,5), (5,6),(6,7), (6,8,8), (7,9,9),(8,10),(9,11), (10,11,11),(11,12),(1,13)]},
			  
			  'Val(NMe)' : {
					'Res_name': 'MNV', 
					'atoms':['N','CA','C','O','CB','CG1','CG2','CN'],
					'atom_type':['N','C','C','O','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (2,5), (5,6),(5,7),(1,8)]},
				
##Peptoids	
			  
			  
			  }

					
atom_positions_beta_up = {'N':[0.000, 0.000, 0.000],
					'CA':[1.221, 0.870, 0.000],
					'CH3':[1.221, 0.870, 0.000],
					'C':[2.442, 0.000, 0.000],
					'O':[2.442,-1.200, 0.000], 
					'CB':[1.221, 1.74, 1.221],
					'CG':[1.221,3.250,1.221],
					'SG':[1.221,3.450,1.221],
					'CG1':[2.442,2.610,1.221],
					'CG2':[0.200, 2.610, 1.421],
					'OG':[1.221,3.250,1.221],
					'OG1':[1.221, 0.300, 2.121],
					'CD':[0.000,4.120,1.221],
					'SD':[0.221,4.220,1.221],
					'CD1':[2.442,4.120,1.221],
					'CD2':[0.000,4.120,1.221],
					'ND1':[2.442,4.120,1.221],
					'ND2':[2.442,4.120,1.221],
					'OD1':[0.000,4.120,1.221],
					'OD2':[2.442,4.120,1.221],
					'CE':[0.000,5.620,1.221],
					'CE1':[2.442,5.341,1.221],
					'CE2':[0.200,5.341,1.221],
					'CE3':[-1.221,3.550,1.600],
					'OE1':[-1.170,3.750,1.421],
					'OE2':[0.200,5.250,1.221],
					'NE1':[1.800,5.541,1.221],
					'NE2':[0.551,5.341,1.221],
					'NE':[0.000,5.620,1.221],
					'NZ':[-1.221,6.490,1.221],
					'CZ2':[-1.000, 6.200, 1.500],
					'CZ3':[-2.200, 4.120, 1.500],
					'CZ':[1.221,6.320,1.221],
					'CH2':[-2.200,5.520,1.500],
					'NH1':[2.521,6.490,1.221],
					'NH2':[1.221,7.990,1.221],
					'OH':[1.221,7.820,1.221],
					'I':[1.221,8.020,1.221],
					'BR':[1.221,8.020,1.221],
					'F1':[-1.000,3.250,1.421],
					'F2':[3.44,3.250,1.421],
					'F3':[-1.000,6.320,1.221],
					'F4':[3.44,6.320,1.221],
					'F5':[1.221,7.820,1.221],
					'CN':[0.000, -1.450, 0.600]}
					
##Peptoid starting positions

peptoid_cis = {     'N':[0.000, 0.000, 0.000],
					'CA':[0.000, 1.450, 0.000],
					'CH3':[0.000, 1.450, 0.000],
					'C':[1.221, 2.250, 0.000],
					'O':[1.221, 3.470, 0.000],
					'CA5':[1.221, -0.870, 0.000],
					
					'CB':[1.221,-2.370,0.000],
					'CB1':[-1.250,-2.250,0.000],
					'CB2':[1.250,-2.250,0.000],
					
					'CG':[2.441,-3.170,0.000],
					'OG':[2.441,-3.170,0.000],
					
					'CG1':[2.441,-3.170,0.000],
					'CG2':[1.441,-3.170,0.000],
					
					#'CG1':[-2.400, -1.450, 0.000],
					'CG3':[-1.250,-3.750,0.000],
					'CG4':[1.250,-3.750,0.000],
					
					'NG':[2.400,-3.570,0.000],
					'OD1':[3.663,-2.370,0.500],
					'OD2':[2.441,-4.390,0.500],
					
					'CD':[2.441,-4.390, 0.000],
					'CD1':[3.663,-2.370,0.500],
					'CD2':[2.441,-4.390,-0.500],
					'CD3':[-2.450,-4.570,0.000],
					'CD4':[-3.600,-2.250,0.000],
					'CD5':[0.000, -4.550, 0.000],
					
					'CE':[3.663,-5.190, 0.000],
					'CE1':[4.863,-3.170,0.500],
					'CE2':[3.663,-5.190,-0.500],
					'CE3':[1.250,-3.750,0.000],
					'CE5':[-3.600,-3.700,0.000],
					'CE6':[0.000,-6.050, 0.000],
					'NE':[-2.500,-6.050,0.000],
					
					'CZ1':[4.884,-4.390,0.000],
					'CZ':[3.663, -6.700,0.000],
					'OZ':[3.663, -6.700,0.000],
					'SZ':[-3.750, -6.850,0.000],
					
					'CZ2':[1.250, -6.850,0.000],
					'CZ3':[1.250, -6.700,0.000],
					'CLH':[-1.250, -8.200,0.000],
					'OD1':[-2.450,-4.800,0.000],
					'OD2':[-0.050,-4.800,0.000],
					
					'CH':[4.884, -7.650,0.000],
					'CH1':[-1.250, -8.150,0.000],
					'CT':[4.884, -9.150,0.000],
					'CI':[6.105, -9.950,0.000],
					'OI':[6.105, -9.950,0.000],
					'CK':[6.105, -11.450,0.000],
					'CL':[-6.250, -12.950,0.000]
					
			  }
peptoid_cis_up = {  'N':[2.985, 0.112, -0.055],
					'CA':[3.944, -0.695, 0.707],
					'CH3':[0.000, 1.450, 0.000],
					'C':[4.683,   0.050,   1.769],
					'O':[3.963,   0.961,   2.243],
					'CA5':[3.560  , 0.907 , -1.257],
					
					'CB':[3.550,   0.042, -2.554],
					'CB1':[4.141,-0.136,-2.101],
					'CB2':[4.335,2.066,-0.648],
					
					'CG':[3.743 ,  0.841 , -3.766],
					'OG':[4.339,   2.601,  -2.540],
					
					'CG1':[3.743 ,  0.841 , -3.766],
					'CG2':[2.743 ,  0.841 , -3.766],
					
					'CG3':[3.529,-1.254,-2.579],
					'CG4':[5.438,0.193,-2.757],
					
					'NG':[3.743 ,  0.841 , -3.766],
					'OD1':[2.611, 1.331,-4.282],
					'OD2':[4.844, 1.077,-4.246],
					
					'CD':[4.799,4.064, -2.459],
					'CD1':[3.173, 2.868,-3.498],
					'CD2':[5.446, 2.184,-2.928],
					'CD3':[4.015,-2.035,-3.631],
					'CD4':[5.928,-0.587,-3.81],
					'CD5':[0.000, -4.550, 0.000],
					
					'CE':[5.259,4.564, -3.735],
					'CE1':[3.634,2.959,-4.817],
					'CE2':[5.887,2.558,-4.293],
					'CE3':[5.219,-1.702,-4.244],
					'CE5':[-3.600,-3.700,0.000],
					'CE6':[0.000,-6.050, 0.000],
					'NE':[-2.500,-6.050,0.000],
					
					'CZ1':[4.966,2.960,-5.243],
					'CZ':[5.719,6.027,-3.654],
					'OZ':[5.719,6.027,-3.654],
					'SZ':[-3.750, -6.850,0.000],
					
					'CZ2':[1.250, -6.850,0.000],
					'CZ3':[1.250, -6.700,0.000],
					'CLH':[5.541,3.188,-7.009],
					
					'CH':[6.179, 6.527,-4.930],
					'CH1':[-1.250, -8.150,0.000],
					'CT':[6.639, 7.990, -4.849],
					'CI':[7.099, 8.490, -6.125],
					'OI':[7.099, 8.490, -6.125],
					'CK':[7.559, 9.953, -6.044],
					'CL':[8.019, 10.453, -7.320]
					
			  }			

peptoid_cis_up_2 = {  'N':[2.985, 0.112, -0.055],
					'CA':[3.944,  -0.695, 0.707],
					'CH3':[0.000, 1.450, 0.000],
					'C':[4.683,   0.050,   1.769],
					'O':[3.879,   -0.737,   2.445],
					'CA5':[3.560  , 0.907 , -1.257],
					
					'CB':[3.550,   0.042, -2.554],
					'CB1':[4.141,-0.136,-2.101],
					'CB2':[4.335,2.066,-0.648],
					
					'CG':[3.743 ,  0.841 , -3.766],
					'OG':[4.339,   2.601,  -2.540],
					
					'CG1':[-2.400, -1.450, 0.000],
					
					'CG3':[3.529,-1.254,-2.579],
					'CG4':[5.438,0.193,-2.757],
					
					'NG':[3.743 ,  0.841 , -3.766],
					'OD1':[2.611, 1.331,-4.282],
					'OD2':[4.844, 1.077,-4.246],
					
					'CD':[4.799,4.064, -2.459],
					'CD1':[3.173, 2.868,-3.498],
					'CD2':[5.446, 2.184,-2.928],
					'CD3':[4.015,-2.035,-3.631],
					'CD4':[5.928,-0.587,-3.81],
					'CD5':[0.000, -4.550, 0.000],
					
					'CE':[5.259,4.564, -3.735],
					'CE1':[3.634,2.959,-4.817],
					'CE2':[5.887,2.558,-4.293],
					'CE3':[5.219,-1.702,-4.244],
					'CE5':[-3.600,-3.700,0.000],
					'CE6':[0.000,-6.050, 0.000],
					'NE':[-2.500,-6.050,0.000],
					
					'CZ1':[4.966,2.960,-5.243],
					'CZ':[5.719,6.027,-3.654],
					'OZ':[5.719,6.027,-3.654],
					'SZ':[-3.750, -6.850,0.000],
					
					'CZ2':[1.250, -6.850,0.000],
					'CZ3':[1.250, -6.700,0.000],
					'CLH':[5.541,3.188,-7.009],
					
					'CH':[6.179, 6.527,-4.930],
					'CH1':[-1.250, -8.150,0.000],
					'CT':[6.639, 7.990, -4.849],
					'CI':[7.099, 8.490, -6.125],
					'OI':[7.099, 8.490, -6.125],
					'CK':[7.559, 9.953, -6.044],
					'CL':[8.019, 10.453, -7.320]
					
			  }			  			  
peptoid_cis_down = { 'N':[5.989,  0.036,   2.236],
					'CA':[6.983 , -0.695 ,  1.381],
					'CH3':[0.000, 1.450, 0.000],
					'C':[7.460 ,  0.000 ,  0.119],
					'O':[6.678 ,  0.686 , -0.548],
					'CA5':[6.499,   0.722 ,  3.221],
					'CB':[6.123 ,  0.294  , 4.701],
					#'CB':[6.835 ,  -0.096  , 4.531],
					#'CB':[5.841 ,  0.152  , 4.593],
					'CB1':[7.48,-0.411,4.198],
					'CB2':[7.331,1.652,2.599],
					'CG':[7.103, 0.773, 5.836],
					'CG1':[7.103, 0.773, 5.836],
					'CG2':[6.103, 0.773, 5.836],
					#'CG':[7.548 ,  0.770 ,  5.417],
					#'CG':[6.604 ,  0.710 ,  5.769],
					'OG':[7.548 ,  0.770 ,  5.417],
					
					#'CG1':[-2.400, -1.450, 0.000],
					'CG3':[8.852,-0.657,4.051],
					'CG4':[6.851,-0.875,5.37],
					
					'NG':[-1.809,3.170,5.097],
					'OD1':[6.639  , 2.097 ,  5.989],
					'OD2':[7.313,  -0.143,   6.631],
					
					'CD':[7.848,0.136, 6.659],
					'CD1':[6.595, 1.767, 6.660],
					'CD2':[8.384, 0.183, 5.911],
					'CD3':[9.570,-1.345,5.035],
					'CD4':[7.564,-1.570,6.348],
					'CD5':[0.000, -4.550, 0.000],
					
					'CE':[8.410,1.169, 7.630],
					'CE1':[7.402, 2.180, 7.742],
					'CE2':[9.151, 0.545, 7.017],
					'CE3':[8.901,  -0.796, 4.840],
					'CE5':[-3.600,-3.700,0.000],
					'CE6':[0.000,-6.050, 0.000],
					'NE':[6.639,0.478, 8.260],
					
					'CZ1':[8.692, 1.576, 7.903],
					'CZ':[8.866  , 0.514  , 8.815],
					'OZ':[8.866  , 0.514  , 8.815],
					'SZ':[-3.750, -6.850,0.000],
					
					'CZ2':[10.227, 0.072, 7.213],
					'CZ3':[1.250, -6.700,0.000],
					'CLH':[8.686  , 1.769  , 9.211],
					
					#'CH':[9.327, 1.458,9.784],
					'CH':[9.606, 2.045,9.050],
					'CH1':[9.606, 2.045,9.050],
					'CT':[9.917, 0.713, 10.976],
					'CI':[10.383, 1.654, 11.943],
					'OI':[10.383, 1.654, 11.943],
					'CK':[10.961, 1.004, 13.066],
					'CL':[6.744, -0.218, 15.733]
					
			  }
			  
peptoid_cis_down_2 = { 'N':[5.989,  0.036,   2.236],
					'CA':[6.983 , -0.695 ,  1.381],
					'CH3':[0.000, 1.450, 0.000],
					'C':[7.460 ,  0.000 ,  0.119],
					'O':[6.750 ,  -0.771 , -0.566],
					'CA5':[6.499,   0.722 ,  3.221],
					'CB':[6.123 ,  0.294  , 4.701],
					#'CB':[6.835 ,  -0.096  , 4.531],
					#'CB':[5.841 ,  0.152  , 4.593],
					'CB1':[7.48,-0.411,4.198],
					'CB2':[7.331,1.652,2.599],
					'CG':[7.103, 0.773, 5.836],
					#'CG':[7.548 ,  0.770 ,  5.417],
					#'CG':[6.604 ,  0.710 ,  5.769],
					'OG':[7.548 ,  0.770 ,  5.417],
					
					'CG1':[-2.400, -1.450, 0.000],
					'CG3':[8.852,-0.657,4.051],
					'CG4':[6.851,-0.875,5.37],
					
					'NG':[-1.809,3.170,5.097],
					'OD1':[6.639  , 2.097 ,  5.989],
					'OD2':[7.313,  -0.143,   6.631],
					
					'CD':[7.848,0.136, 6.659],
					'CD1':[6.595, 1.767, 6.660],
					'CD2':[8.384, 0.183, 5.911],
					'CD3':[9.570,-1.345,5.035],
					'CD4':[7.564,-1.570,6.348],
					'CD5':[0.000, -4.550, 0.000],
					
					'CE':[8.410,1.169, 7.630],
					'CE1':[7.402, 2.180, 7.742],
					'CE2':[9.151, 0.545, 7.017],
					'CE3':[8.901,  -0.796, 4.840],
					'CE5':[-3.600,-3.700,0.000],
					'CE6':[0.000,-6.050, 0.000],
					'NE':[6.639,0.478, 8.260],
					
					'CZ1':[8.692, 1.576, 7.903],
					'CZ':[8.866  , 0.514  , 8.815],
					'OZ':[8.866  , 0.514  , 8.815],
					'SZ':[-3.750, -6.850,0.000],
					
					'CZ2':[10.227, 0.072, 7.213],
					'CZ3':[1.250, -6.700,0.000],
					'CLH':[8.686  , 1.769  , 9.211],
					
					#'CH':[9.327, 1.458,9.784],
					'CH':[9.606, 2.045,9.050],
					'CH1':[9.606, 2.045,9.050],
					'CT':[9.917, 0.713, 10.976],
					'CI':[10.383, 1.654, 11.943],
					'OI':[10.383, 1.654, 11.943],
					'CK':[10.961, 1.004, 13.066],
					'CL':[6.744, -0.218, 15.733]
					
			  }			  
peptoid_trans_up = {  'N':[3.743, 0.100, -0.053],
					'CA':[5.093, 0.266, 0.514],
					'CH3':[0.000, 1.450, 0.000],
					'C':[6.305,   0.029,   -0.340],
					'O':[6.419,   -1.145,  -0.733],
					'CA5':[2.964  , -1.058 , 0.418],
					
					'CB':[3.875,   -2.186, 1.093],
					'CB1':[-1.250,-2.250,0.000],
					'CB2':[1.250,-2.250,0.000],
					
					'CG':[3.101,   -3.409,  1.234],
					'OG':[4.339,   2.601,  -2.540],
					
					'CG1':[-2.400, -1.450, 0.000],
					'CG3':[-1.250,-3.750,0.000],
					'CG4':[1.250,-3.750,0.000],
					
					'NG':[3.101,   -3.409,  1.234],
					'OD1':[1.882,-4.526,1.341],
					'OD2':[3.412, -3.451,1.309],
					
					'CD':[4.799,4.064, -2.459],
					'CD1':[-0.328,1.027,-4.421],
					'CD2':[1.857,-0.019,-4.286],
					'CD3':[-2.450,-4.570,0.000],
					'CD4':[-3.600,-2.250,0.000],
					'CD5':[0.000, -4.550, 0.000],
					
					'CE':[5.259,4.564, -3.735],
					'CE1':[0.046,1.649,-5.615],
					'CE2':[2.227,0.609,-5.48],
					'CE3':[1.250,-3.750,0.000],
					'CE5':[-3.600,-3.700,0.000],
					'CE6':[0.000,-6.050, 0.000],
					'NE':[-2.500,-6.050,0.000],
					
					'CZ1':[1.322,1.440,-6.143],
					'CZ':[5.719,6.027,-3.654],
					'OZ':[5.719,6.027,-3.654],
					'SZ':[-3.750, -6.850,0.000],
					
					'CZ2':[1.250, -6.850,0.000],
					'CZ3':[1.250, -6.700,0.000],
					'CLH':[-1.250, -8.200,0.000],
					
					'CH':[6.179, 6.527,-4.930],
					'CH1':[-1.250, -8.150,0.000],
					'CT':[6.639, 7.990, -4.849],
					'CI':[7.099, 8.490, -6.125],
					'OI':[7.099, 8.490, -6.125],
					'CK':[7.559, 9.953, -6.044],
					'CL':[8.019, 10.453, -7.320]
					
			  }			  
peptoid_trans_down = { 'N':[7.245,  1.049,   -0.666],
					'CA':[8.373 , 0.525 ,  -1.511],
					'CH3':[0.000, 1.450, 0.000],
					'C':[9.755 ,  1.084,  -0.896],
					'O':[10.262, 2.133, -1.356],
					'CA5':[7.091,   2.399,  -0.146],
					
					'CB':[7.504,  2.529, 1.297],
					'CB1':[-1.250,-2.250,0.000],
					'CB2':[1.250,-2.250,0.000],
					
					'CG':[7.348,  4.008, 1.744],
					'OG':[6.604 ,  0.710 ,  5.769],
					
					'CG1':[-2.400, -1.450, 0.000],
					'CG3':[-1.250,-3.750,0.000],
					'CG4':[1.250,-3.750,0.000],
					
					'NG':[-1.809,3.170,5.097],
					'OD1':[5.382,2.015,-3.215],
					'OD2':[6.912,0.693,-4.752],
					
					'CD':[5.516,1.036, 7.084],
					'CD1':[6.003, 4.450,  1.926],
					'CD2':[8.393,  4.864,  2.240],
					'CD3':[-2.450,-4.570,0.000],
					'CD4':[-3.600,-2.250,0.000],
					'CD5':[0.000, -4.550, 0.000],
					
					'CE':[6.639,0.478, 8.260],
					'CE1':[5.683, 5.694, 2.409],
					'CE2':[8.033, 6.179,  2.617],
					'CE3':[7.313,  -1.211,   6.805],
					'CE5':[-3.600,-3.700,0.000],
					'CE6':[0.000,-6.050, 0.000],
					'NE':[6.639,0.478, 8.260],
					
					'CZ1':[6.716, 6.566,2.630],
					'CZ':[5.551  , 0.804  , 9.575],
					'OZ':[5.551  , 0.804  , 9.575],
					'SZ':[-3.750, -6.850,0.000],
					
					'CZ2':[8.606  , -0.271  , 7.820],
					'CZ3':[1.250, -6.700,0.000],
					'CLH':[-1.250, -8.200,0.000],
					'OD1':[-0.991,3.948,5.581],
					'OD2':[-1.667,1.808,5.257],
					
					'CH':[6.647, 0.246,10.751],
					'CH1':[8.650, 2.170,9.271],
					'CT':[5.586, 0.572, 12.066],
					'CI':[6.709, 0.014, 13.242],
					'OI':[6.709, 0.014, 13.242],
					'CK':[5.621, 0.340, 14.557],
					'CL':[6.744, -0.218, 15.733]
					
			  }



peptoid_sigma = {	'N':[-3.100,0.000, 1.214],
					'CA':[0.55, -0.094, 1.857],
					'CH3':[3.974, -0.094, 1.857],
					'C':[5.151, -0.698, 1.134],
					'O':[5.157,-1.915, 0.951],
					'CA5':[2.052, -1.648, 1.869],
					
					'CB':[0.993,-1.063,2.822],
					
					'CG':[0.263,-2.122,3.477],
					'OG':[0.263,-2.122,3.477],
					
					'CG1':[-2.400, -1.450, 0.000],
					'NG':[0.263,-2.122,3.477],
					'OG1':[-2.400,-2.250,0.000],
					'OG2':[-1.202,3.487,-0.416],
					
					'CD':[-0.796,-1.537,4.43],
					'CD1':[-2.450,-4.500,0.000],
					'CD2':[-0.050,-4.500,0.000],
					'CD3':[-2.450,-4.570,0.000],
					'CD4':[-3.600,-2.250,0.000],
					
					'CE':[-1.526,-2.596,5.085],
					'CE1':[-2.450,-5.900,0.000],
					'CE2':[-0.050,-5.900,0.000],
					'CE3':[-3.600,-3.700,0.000],
					'NE':[-1.526,-2.596,5.085],
					
					'CZ':[-2.585, -2.011,6.038],
					'OZ':[-2.585, -2.011,6.038],
					'CZ1':[-1.250, -6.700,0.000],
					'CLH':[-1.250, -8.200,0.000],
					'OD1':[-2.450,-4.800,0.000],
					'OD2':[-0.050,-4.800,0.000],
					
					'CH':[-3.315, -3.070,6.693],
					'CT':[-4.374, -2.485,7.646],
					'CI':[-5.104, -3.544,8.301],
					'OI':[-5.104, -3.544,8.301],
					'CK':[-6.163, -2.959,9.254],
					'CL':[-6.893, -4.018,9.909],
					}
					
hao = {				'N':[0.000, 0.000, 0.000],
					'N02':[1.221, 0.870, 0.000],
					'C08':[2.442, 0.000, 0.000],
					'O06':[2.442,-1.200, 0.000],
					'C06':[3.663, 0.870, 0.000],
					'C05':[3.663, 2.370, 0.000],
					'OM':[2.442, 3.200, 0.000],
					'CM':[2.442, 4.600, 0.000],
					'C04':[4.884, 3.170, 0.000],
					'C07':[4.884, 0.000, 0.000],
					'C03':[6.106, 2.300, 0.000],
					'C02':[6.106, 0.870, 0.000],
					'N01':[7.307, 0.00, 0.000],
					'C01':[8.508, 0.870, 0.000],
					'O01':[8.508, 2.270, 0.000],
					'C':[9.708, 0.000, 0.000],
					'O':[9.708,-1.200, 0.000]}
					
peptoid_tube = {	'N':[0.000, 0.000, 0.000],
					'CA':[1.221, 0.870, 0.000],
					'CH3':[1.221, 0.870, 0.000],
					'C':[2.442, 0.000, 0.000],
					'O':[2.442,-1.200, 0.000],
					'CA5':[0.000, -1.450, 0.000],
					
					'CB':[-1.250,-2.250,0.000],
					
					'CG':[-1.250,-3.750,0.000],
					'OG':[-1.250,-3.750,0.000],
					
					'CG1':[-2.400, -1.450, 0.000],
					'NG':[-1.250,-3.750,0.000],
					'OG1':[-2.400,-2.250,0.000],
					'OG2':[-1.250,-4.100,0.000],
					
					'CD':[-2.450,-4.500,0.000],
					'CD1':[-2.450,-4.500,0.000],
					'CD2':[-0.050,-4.500,0.000],
					'CD3':[-2.450,-4.570,0.000],
					'CD4':[-3.600,-2.250,0.000],
					
					'CE':[-2.450,-5.900,0.000],
					'CE1':[-2.450,-5.900,0.000],
					'CE2':[-0.050,-5.900,0.000],
					'CE3':[-3.600,-3.700,0.000],
					'NE':[-2.450,-6.000,0.000],
					
					'CZ':[-3.650, -6.700,0.000],
					'OZ':[-3.650, -6.700,0.000],
					'CZ1':[-1.250, -6.700,0.000],
					'OD1':[-2.450,-4.800,0.000],
					'OD2':[-0.050,-4.800,0.000],
					
					'CH':[-3.650, -8.200,0.000],
					'CT':[-4.850, -8.800,0.000],
					'CI':[-4.850, -10.300,0.000],
					'OI':[-4.850, -10.300,0.000],
					'CK':[-6.050, -11.100,0.000],
					'CL':[-6.050, -12.600,0.000],
					}