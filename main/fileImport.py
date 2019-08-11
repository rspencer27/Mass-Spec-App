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
				peptoid_names = ['Nab', 'Nae', 'Nbsa', 'Nbn', 'Namd', 'Nce', 'Ncm', 'Ncp', 'Ncpr', 'Ndpe','Nfu','Nhe','Nia','Nme','Nmp','Nbu','Npe','Npip','Npp','Ntrp','Ntyr','Nipr','Nph','Nspe','Nrpe','Nime','Npyr', 'Nte','Ndc','Ncpe','Ntm','Nreh','Nseh','Npm','Netm','Nib']
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
								'Npip' : [10,9,1,3,0,0,0,0,0],
								'Npp' :[9,14,2,2,0,0,0,0,0],
								'Ntrp' : [12,12,2,1,0,0,0,0,0],
								'Ntyr' : [10,11,1,2,0,0,0,0,0],
								'Nipr' :[5,9,1,1,0,0,0,0,0],
								'Nph' : [8,7,1,1,0,0,0,0,0] ,
								'Nspe' : [10,11,1,1,0,0,0,0,0] ,
								'Nrpe' : [10,11,1,1,0,0,0,0,0] ,
								'Nime' : [7,9,3,1,0,0,0,0,0],
								'Npyr': [8,8,2,1,0,0,0,0,0],
								'Ndc': [12,23,1,1,0,0,0,0,0],
								'Nte': [9,17,1,4,0,0,0,0,0],
								'Ncpe': [10,10,1,1,0,0,1,0,0],
								'Ntm' : [6,11,1,2,1,0,0,0,0],
								'Nseh' : [6,11,1,2,1,0,0,0,0],
								'Nreh' : [6,11,1,2,1,0,0,0,0],
								'Npm' : [3,6,1,2,1,0,1,0,0],
								'Netm' : [3,6,1,1,1,0,0,0,0],
								'Nib' : [6,11,1,1,0,0,0,0,0]
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
				peptoid_comments=['aminobutyl', 'aminoethyl', 'aminoethyl-benzenesulfonamide', 'benzylamine', 'glycinamide', 'beta-alanine', 'glycine', 'cyclopentylamine', 'cyclopropanemethylamine','2,2diphenylethanamine','furfurylamine','ethanolamine','isoamyl','2-methoxyethylamine','methoxypyridyl','n-butyl','phenethylamine','piperonylamine','propylpyrrolidinone','tryptamine','tyramine','isopropylamine','aniline','(S)-(-)-alpha-methylbenzylamine', '(R)-(-)-alpha-methylbenzylamine','histamine','3-aminomethyl-pyridine', '2-[2-(2-methoxyethoxy)ethoxy]-ethanamine', 'decylamine','para-chloro phenethylamine','Nreh','Nseh','Npm','Netm','Nib']
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
					'Res_name': 'NAMD', 
					'atoms':['N','CA','C','O','CA5','CB', 'NG', 'OG'],
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
			'Nspe' : {
					'Res_name': 'NSPE', 
					'atoms':['N','CA','C','O','CA5','CB2','CB1', 'CG3','CG4','CD3','CD4','CE3'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(5,7), (7,8,8),(8,10),(10,12,12),(12,11),(11,9,9), (9,7)]},
			'Nrpe' : {
					'Res_name': 'NRPE', 
					'atoms':['N','CA','C','O','CA5','CB2','CB1', 'CG3','CG4','CD3','CD4','CE3'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(5,7), (7,8,8),(8,10),(10,12,12),(12,11),(11,9,9), (9,7)]},
			'Ncpe' : {
					'Res_name': 'NPC', 
					'atoms':['N','CA','C','O','CA5','CB', 'CG','CD1','CD2','CE1','CE2','CZ1','CLH'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','C','Cl'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (7,8), (7,9,9),(8,10,10),(10,12),(9,11),(11,12,12), (12,13)]},
			  'Nae' : {
					'Res_name': 'NAE', 
					'atoms':['N','CA','C','O','CA5','CB', 'NG'],
					'atom_type':['N','C','C','O','C','C','N1+'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7)]},
			 'Nhe' : {
					'Res_name': 'NHE', 
					'atoms':['N','CA','C','O','CA5','CB', 'OG'],
					'atom_type':['N','C','C','O','C','C','O'],
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
			'Nib' : {
					'Res_name': 'NIB', 
					'atoms':['N','CA','C','O','CA5','CB', 'CG1','CG2'],
					'atom_type':['N','C','C','O','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7),(6,8)]},
			
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
			  'Nreh' : {
					'Res_name': 'NREH', 
					'atoms':['N','CA','C','O','CA5','CB', 'CG1','CD1','CE','CZ', 'CG2', 'CD2'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (7,8), (8,9),(9,10),(6,11),(11,12)]},
			  'Nseh' : {
					'Res_name': 'NSEH', 
					'atoms':['N','CA','C','O','CA5','CB', 'CG1','CD1','CE','CZ', 'CG2', 'CD2'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (7,8), (8,9),(9,10),(6,11),(11,12)]},
			  'Npm': {
					'Res_name': 'NPM', 
					'atoms':['N','CA','C','O','CA5','PB', 'OB1','OB2','OB3'],
					'atom_type':['N','C','C','O','C','P','O','O','O'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (6,8), (6,9,9)]},
			'Netm': {
					'Res_name': 'NETM', 
					'atoms':['N','CA','C','O','CA5','CB', 'NG','CD1','CD2', 'CD3'],
					'atom_type':['N','C','C','O','C','C','N','C','C','C'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (7,8), (7,9), (7,10)]},
			'Ntyr' : {
					'Res_name': 'NTYR', 
					'atoms':['N','CA','C','O','CA5','CB', 'CG','CD1','CD2','CE1','CE2','CZ1','OH'],
					'atom_type':['N','C','C','O','C','C','C','C','C','C','C','C','O'],
					'link_rec':[(1,2), (2,3), (3,4,4), (1,5), (5,6),(6,7), (7,8), (7,9,9),(8,10,10),(10,12),(9,11),(11,12,12), (12,13)]}
				}
				break
		return peptoid_coordinates
		
