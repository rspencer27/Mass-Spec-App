from urllib.request import urlopen


class Import_Peptoid_Names(object):

	def __init__(self):
		pass
	def get_peptoid_names(self):
		peptoid_names=[]
		while True:
			try:
				f = urlopen('http://www.peptoids.org/peptoids_website_additional_files/peptoid_Names.txt')
				peptoid_names = eval(f.read())
				f.close()
				break
			except:
				print ("Could not open custom_Names file. Please upload file.")
				peptoid_names = ['Nab', 'Nae', 'Nbsa', 'Nbn', 'Namd', 'Nce', 'Ncm', 'Ncp', 'Ncpr', 'Ndpe','Nfu','Nhe','Nia','Nme','Nmp','Nbu','Npe','Npip','Npp','Ntrp','Ntyr','Nipr','Nph','Nspe','Nime','Npyr', 'Nte','Ndc','Ncpe','Ntm','Npeo','Npem','Npep'] 
				break
		return peptoid_names
		
	def get_peptoid_mass(self):
		peptoid_masses={}
		while True:	
			try:
				f1 = urlopen('http://www.peptoids.org/peptoids_website_additional_files/peptoid_Masses2.txt')
				peptoid_masses = eval(f1.read())
				f1.close()
				break
			except:
				print ("Could not open custom_Masses file. Please upload file.")
				peptoid_masses={
								'Nab' : [6,12,2,1,0,0,0,0,0], 
								'Nae' : [4,8,2,1,0,0,0,0,0], 
								'Nbsa' : [10,12,2,3,1,0,0,0,0], 
								'Nbn' : [9,9,1,1,0,0,0,0,0], 
								'Namd' : [4,6,2,2,0,0,0,0,0], 
								'Nce' : [5,7,1,3,0,0,0,0,0], 
								'Ncm' : [4,5,1,3,0,0,0,0,0], 
								'Ncp' : [7,11,1,1,0,0,0,0,0], 
								'Ncpr' : [6,9,1,1,0,0,0,0,0] , 
								'Ndpe' : [16,15,1,1,0,0,0,0,0],
								'Nfu' : [7,7,1,2,0,0,0,0,0] ,
								'Nhe' : [4,7,1,2,0,0,0,0,0],
								'Nia' : [7,13,1,1,0,0,0,0,0],
								'Nme' : [5,9,1,2,0,0,0,0,0],
								'Nmp' : [8,8,2,2,0,0,0,0,0],
								'Nbu' : [6,11,1,1,0,0,0,0,0],
								'Npe' : [10,11,1,1,0,0,0,0,0],
								'Npeo': [11,13,1,1,0,0,0,0,0],
								'Npem': [11,13,1,1,0,0,0,0,0],
								'Npep': [11,13,1,1,0,0,0,0,0],
								'Npip' : [10,9,1,3,0,0,0,0,0],
								'Npp' :[9,14,2,2,0,0,0,0,0],
								'Ntrp' : [12,12,2,1,0,0,0,0,0],
								'Ntyr' : [10,11,1,2,0,0,0,0,0],
								'Nipr' :[5,9,1,1,0,0,0,0,0],
								'Nph' : [8,7,1,1,0,0,0,0,0] ,
								'Nspe' : [10,11,1,1,0,0,0,0,0] ,
								'Nime' : [7,9,3,1,0,0,0,0,0],
								'Npyr': [8,8,2,1,0,0,0,0,0],
								'Ndc': [12,23,1,1,0,0,0,0,0],
								'Nte': [9,17,1,4,0,0,0,0,0],
								'Ntm' : [6,11,1,2,1,0,0,0,0],
								'Ncpe': [9,17,1,4,0,0,1,0,0]
								}
				break
		return peptoid_masses
	
	def get_peptoid_comments(self):
		peptoid_comments=[]
		while True:	
			try:
				f2 = urlopen('http://www.peptoids.org/peptoids_website_additional_files/peptoid_Comments.txt')
				peptoid_comments = eval(f2.read())
				f2.close()
				break
			except:
				print("Could not open custom_Comments file. Please upload file")
				peptoid_comments=['aminobutyl', 'aminoethyl', 'aminoethyl-benzenesulfonamide', 'benzylamine', 'glycinamide', 'beta-alanine', 'glycine', 'cyclopentylamine', 'cyclopropanemethylamine','2,2diphenylethanamine','furfurylamine','ethanolamine','isoamyl','2-methoxyethylamine','methoxypyridyl','n-butyl','phenethylamine','piperonylamine','propylpyrrolidinone','tryptamine','tyramine','isopropylamine','aniline','(S)-(-)-alpha-methylbenzylamine','histamine','3-aminomethyl-pyridine', '2-[2-(2-methoxyethoxy)ethoxy]-ethanamine', 'decylamine','para-chloro phenethylamine', 'thiol-amine', 'Npe-o-methyl','Npe-m-methyl','Npe-p-methyl' ]
				break
		return peptoid_comments
	
	def get_peptoid_coordinates(self):
		peptoid_coordinates={}
		while True:	
			try:
				f2 = urlopen('http://www.peptoids.org/peptoids_website_additional_files/peptoid_Coordinates.txt')
				peptoid_coordinates = eval(f2.read())
				f2.close()
				break
			except:
				print("Could not open custom_Comments file. Please upload file")
				peptoid_coordinates= {
			'Nab' : {
					'Res_name': 'NAB', 
					'atoms':['N','CA','C','O','CA5','CB', 'CG', 'CD', 'NE'],
					'atom_type':['N','C','C','O','C','C','C','C','N+1'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (7,8),(8,9)]},
					
			  'Nbn' : {
					'Res_name': 'NBN', 
					'atoms':['N','CA','C','O','CA5','CB', 'CG','CG1','CD3','CD4','CE3'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6), (6,7,7) ,(6,8),(7,9), (8,10,10), (9,11,11),(10,11)]},
			 
			  'Namd' : {
					'Res_name': 'NMD', 
					'atoms':['N','CA','C','O','CA5','CB', 'NG', 'OG1'],
					'atom_type':['N','C','C','O','C','C','N','O'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (6,8,8)]},
			 
			  'Ncm' : {
					'Res_name': 'NCM', 
					'atoms':['N','CA','C','O','CA5','CB', 'OG1', 'OG2'],
					'atom_type':['N','C','C','O','C','C','O','O'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (6,8,8)]},
					
			  'Npe' : {
					'Res_name': 'NPE', 
					'atoms':['N','CA','C','O','CA5','CB', 'CG','CD1','CD2','CE1','CE2','CZ1'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (7,8), (7,9,9),(8,10,10),(10,12),(9,11),(11,12,12)]},
			  'Ncpe' : {
					'Res_name': 'NCPE', 
					'atoms':['N','CA','C','O','CA5','CB', 'CG','CD1','CD2','CE1','CE2','CZ1','CLH'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','C','Cl'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (7,8), (7,9,9),(8,10,10),(10,12),(9,11),(11,12,12), (12,13)]},
			  'Nae' : {
					'Res_name': 'NAE', 
					'atoms':['N','CA','C','O','CA5','CB', 'NG'],
					'atom_type':['N','C','C','O','C','C','N1+'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7)]},
			
			  'Nce' : {
					'Res_name': 'NCE', 
					'atoms':['N','CA','C','O','CA5','CB', 'CG','OD1','OD2'],
					'atom_type':['N','C','C','O','C','C','C','O','O1-'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (7,8,8), (7,9)]},
			 
			  'Ndc' : {
					'Res_name': 'NDC', 
					'atoms':['N','CA','C','O','CA5','CB', 'CG', 'CD', 'CE', 'CZ', 'CH','CT','CI','CK'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (7,8), (8,9),(9,10),(10,11),(11,12), (12,13),(13,14)]},
			  
			  'Nte' : {
					'Res_name': 'NTE', 
					'atoms':['N','CA','C','O','CA5','CB', 'OG', 'CD', 'CE', 'OZ', 'CH','CT','OI','CK'],
					'atom_type':['N','C','C','O','C','C','O','C','C','O','C','C','O','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (7,8), (8,9),(9,10),(10,11),(11,12), (12,13),(13,14)]},
			  
			  'Nia' : {
					'Res_name': 'NIA', 
					'atoms':['N','CA','C','O','CA5','CB', 'CG','CD1','CD2'],
					'atom_type':['N','C','C','O','C','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (7,8), (7,9)]},
			
			  'Ncp' : {
					'Res_name': 'NCP', 
					'atoms':['N','CA','C','O','CA5','CB1', 'CB2','CG3','CG4'],
					'atom_type':['N','C','C','O','C','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6), (5,7), (6,8), (7,9), (8,9)]},
					
			  'Ntm' : {
					'Res_name': 'NTM', 
					'atoms':['N','CA','C','O','CA5','CB', 'OG', 'CD', 'CE', 'SZ'],
					'atom_type':['N','C','C','O','C','C','O','C','C','S'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (7,8), (8,9),(9,10)]},
			  'Nme' : {
					'Res_name': 'NME', 
					'atoms':['N','CA','C','O','CA5','CB', 'OG', 'CD'],
					'atom_type':['N','C','C','O','C','C','O','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (7,8)]},
			  'Npeo' : {
					'Res_name': 'NPEO', 
					'atoms':['N','CA','C','O','CA5','CB', 'CG','CD1','CD2','CE1','CE2','CZ1', 'CE3'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (7,8), (7,9,9),(8,10,10),(10,12),(9,11),(11,12,12), (9,13)]},
			  'Npem' : {
					'Res_name': 'NPEM', 
					'atoms':['N','CA','C','O','CA5','CB', 'CG','CD1','CD2','CE1','CE2','CZ1', 'CZ2'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (7,8), (7,9,9),(8,10,10),(10,12),(9,11),(11,12,12), (11,13)]},
			  'Npep' : {
					'Res_name': 'NPEP', 
					'atoms':['N','CA','C','O','CA5','CB', 'CG','CD1','CD2','CE1','CE2','CZ1', 'CH1'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (7,8), (7,9,9),(8,10,10),(10,12),(9,11),(11,12,12), (12,13)]}		
					}
				break
		return peptoid_coordinates
	
	def get_peptoid_atoms(self):
		peptoid_atoms={}
		while True:	
			try:
				f2 = urlopen('http://www.peptoids.org/peptoids_website_additional_files/peptoid_atoms.txt')
				peptoid_atoms = eval(f2.read())
				f2.close()
				break
			except:
				print("Could not open custom_Comments file. Please upload file")
				peptoid_atoms = {	        
				    'N':[0.000, 0.000, 0.000],
					'CA':[1.221, 0.870, 0.000],
					'CH3':[1.221, 0.870, 0.000],
					'C':[2.442, 0.000, 0.000],
					'O':[2.442,-1.220, 0.000],
					'CA5':[0.000, -1.450, 0.000],
					
					'CB':[-1.250,-2.250,0.000],
					'CB1':[-1.250,-2.250,0.000],
					'CB2':[1.250,-2.250,0.000],
					
					'CG':[-1.250,-3.750,0.000],
					'OG':[-1.250,-3.750,0.000],
					
					'CG1':[-2.400, -1.450, 0.000],
					'CG3':[-1.250,-3.750,0.000],
					'CG4':[1.250,-3.750,0.000],
					
					'NG':[-1.250,-3.750,0.000],
					'OG1':[-2.400,-2.250,0.000],
					'OG2':[-1.250,-3.750,0.000],
					
					'CD':[-2.500,-4.550,0.000],
					'CD1':[-2.450,-4.500,0.000],
					'CD2':[-0.050,-4.500,0.000],
					'CD3':[-2.450,-4.570,0.000],
					'CD4':[-3.600,-2.250,0.000],
					
					'CE':[-2.500,-6.050,0.000],
					'CE1':[-2.450,-5.900,0.000],
					'CE2':[-0.050,-5.900,0.000],
					'CE3':[1.250,-3.750,0.000],
					'CE5':[-3.600,-3.700,0.000],
					'NE':[-2.500,-6.050,0.000],
					
					'CZ':[-3.750, -6.850,0.000],
					'CZ1':[-1.250, -6.700,0.000],
					'OZ':[-3.750, -6.850,0.000],
					'SZ':[-3.750, -6.850,0.000],
					
					'CZ2':[1.250, -6.700,0.000],
					'CLH':[-1.250, -8.200,0.000],
					'OD1':[-2.450,-4.800,0.000],
					'OD2':[-0.050,-4.800,0.000],
					
					'CH':[-3.750, -8.350,0.000],
					'CH1':[-1.250, -8.150,0.000],
					'CT':[-5.000, -9.150,0.000],
					'CI':[-5.000, -10.650,0.000],
					'OI':[-5.000, -10.650,0.000],
					'CK':[-6.250, -11.450,0.000],
					'CL':[-6.250, -12.950,0.000],
					}
				break
		return peptoid_atoms
		
