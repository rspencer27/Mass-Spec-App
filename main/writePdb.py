from main.PDB_coord import *
from tkinter.filedialog import asksaveasfilename
import math
from main.fileImport import Import_Peptoid_Names

ipn = Import_Peptoid_Names()
chain = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f','g','h','i','j','k','l',\
		'm','n','o','p','q','r','s','t','u','v','w','x','y','z','AA','BB','CC','DD','EE','FF','GG','HH','II','JJ','KK','LL','MM','NN','OO','PP','QQ','RR','SS','TT',\
		'UU','VV','WW','XX','YY','ZZ','aa','bb','cc','dd','ee','ff','gg','hh','ii','jj','kk','ll','mm','nn','oo','pp','qq','rr','ss','tt','uu','vv','ww','xx','yy','zz',\
		'A1','A2','A3','A4','A5','A6','A7','A8','A9','B1','B2','B3','B4','B5','B6','B7','B8','B9','C1','C2','C3','C4','C5','C6','C7','C8','C9','D1','D2','D3','D4','D5',\
		'D6','D7','D8','D9','E1','E2','E3','E4','E5','E6','E7','E8','E9']
		

class write_Pdb(object):
	
	def __init__(self,sequence):
		pass
	
	def write_beta_sheet(self,sequence, start_residue, chain_ID):
		fname = asksaveasfilename(filetypes=(("pdb", "*.pdb"),
                                           ("All files", "*.*") ),
								  defaultextension='.pdb')
		try:							   
			fobj = open(fname, 'w')
		except:
				pass
		print("Sequence file to be converted is")
		print(sequence)
		amino_acids = ['Ac','A','R','N','D','C','E','N','G','H','I','L','K','M','F','P','S','T','W', 'Y', 'Q','V',
						'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Ser', 'Thr', 'Trp', 'Tyr', 'Gln','Val']
		locations = []
		locations_temp=[]
		coordinate=[]
		coordinate2 = []
		temp_coordinate=[]
		total_atoms = 1
		peptoids = ipn.get_peptoid_names()
		writePdbAA.update(ipn.get_peptoid_coordinates())
		if start_residue == '':
			print("Empty residue")
			resid_num = 1
			start_residue = 1
		else:
			start_residue = int(start_residue)
			resid_num = int(start_residue)
		if chain_ID == '':
			chain_ID = 'A'
		else:
			pass
		for i,k in enumerate(sequence):
			res_name = writePdbAA.get(k)['Res_name']
			
			atoms = writePdbAA.get(k)['atoms']
			
			atom_type = writePdbAA.get(k)['atom_type']
			sr = start_residue + i
			print('Residue Number')
			print(sr)
			if i%2 == 0:
				for l,m in enumerate(atoms):
					if k not in peptoids:
						if k == 'Hao':
							temp_coordinate = hao.get(m)
						else:
							temp_coordinate = atom_positions_beta_up.get(m)
					else:
						temp_coordinate = peptoid_sigma.get(m)
					for n,p in enumerate(temp_coordinate):
						coordinate.append(p)
					if m == 'CG' and k == 'P':
						coordinate[0] = 1.221
						coordinate[1] = 0.000
						coordinate[2] = 1.921
					elif m == 'CD' and k == 'P':
						coordinate[0] = 0.000
						coordinate[1] = -0.870
						coordinate[2] = 1.521
					else:
						pass
					if k in peptoids:
						coordinate[0] = coordinate[0] + (resid_num*6.178/2)
						coordinate[1] = coordinate[1] + (resid_num*0.097/2)
						coordinate[2] = coordinate[2] + (resid_num*0.705/2)
					else:	
						coordinate[0] = coordinate[0] + resid_num*3.675
					if atom_type[l] == 'Br':
						fobj.write('HETATM' + str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + ' '+chain_ID + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*8 + str(atom_type[l]).rjust(4,' ')+ '\n')
					elif k not in amino_acids:
						fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + ' '+chain_ID + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*11 + str(atom_type[l])+ '\n')
					else:
						fobj.write('ATOM' +  str(total_atoms).rjust(7,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + ' '+chain_ID + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*11 + str(atom_type[l])+ '\n')
					del coordinate[:]
					total_atoms+=1
			else:
				for l,m in enumerate(atoms):
					if k not in peptoids:
						if k == 'Hao':
							temp_coordinate = hao.get(m)
						else:
							temp_coordinate = atom_positions_beta_up.get(m)
					else:
						temp_coordinate = peptoid.get(m)
					for n,p in enumerate(temp_coordinate):
						coordinate.append(p)
					print(k)
					if k in peptoids:
						if m == 'N':
							coordinate[0] = 0.000
							coordinate[1] = 0.000
							coordinate[2] = 0.000
						elif m == 'CA':
							coordinate[0] = 1.104
							coordinate[1] = -0.537
							coordinate[2] = -0.738
						elif m == 'C':
							coordinate[0] = 2.371
							coordinate[1] = -0.013
							coordinate[2] = 0.049
						elif m == 'O':
							coordinate[0] = 2.978
							coordinate[1] = 0.932
							coordinate[2] = -0.454
						elif m == 'CA5':
							coordinate[0] = 0.034
							coordinate[1] = 1.471
							coordinate[2] = -0.017
						elif m == 'CB':
							coordinate[0] = -0.589
							coordinate[1] = 2.180
							coordinate[2] = -0.972
						elif m == 'CG' or m == 'NG' or m =='OG':
							coordinate[0] = -0.555
							coordinate[1] = 3.487
							coordinate[2] = -0.989
						elif m == 'CD':
							coordinate[0] = -1.178
							coordinate[1] = 4.196
							coordinate[2] = -1.944
						elif m == 'CE':
							coordinate[0] = -1.144
							coordinate[1] = 5.503
							coordinate[2] = -1.961
						elif m == 'CZ' or m =='OZ':
							coordinate[0] = -1.767
							coordinate[1] = 6.212
							coordinate[2] = -2.916
						elif m == 'CH':
							coordinate[0] = -1.733
							coordinate[1] = 7.519
							coordinate[2] = -2.933
						elif m == 'CT':
							coordinate[0] = -2.356
							coordinate[1] = 8.228
							coordinate[2] = -3.888
						elif m == 'CI' or m == 'OI':
							coordinate[0] = -2.322
							coordinate[1] = 9.535
							coordinate[2] = -3.905
						elif m == 'CK':
							coordinate[0] = -2.945
							coordinate[1] = 10.244
							coordinate[2] = -4.860
						elif m == 'CL':
							coordinate[0] = -2.911
							coordinate[1] = 11.551
							coordinate[2] = -4.877
						else:
							pass
					else:
						if m == 'N':
							coordinate[0] = 0.000
							coordinate[1] = 1.121
						elif m == 'CA':
							coordinate[1] = 0.000
						elif m == 'C':
							coordinate[1] = 1.121
						elif m == 'O':
							coordinate[1] = 2.121
						elif m == 'CG' and k == 'P':
							coordinate[0] = 0.000
							coordinate[1] = -0.870
							coordinate[2] = -1.521
						elif m == 'CD' and k == 'P':
							coordinate[0] = 0.000
							coordinate[1] = 0.870
							coordinate[2] = -1.521
						else:
							coordinate[1] = coordinate[1]*(-1) + 1.170
							coordinate[2] = coordinate[2]*(-1)
					if k in peptoids:
						coordinate[0] = coordinate[0] + (resid_num*6.178/2) + 6.178/2
						coordinate[1] = coordinate[1] + (resid_num*0.097/2) + 0.097/2
						coordinate[2] = coordinate[2] + (resid_num*0.705/2) + 0.705/2
					else:
						coordinate[0] = coordinate[0] + resid_num*3.675
					if atom_type[l] == 'Br':
						fobj.write('HETATM' + str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + ' '+chain_ID + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*8 + str(atom_type[l]).rjust(4,' ')+ '\n')
					elif k not in amino_acids:
						fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + ' '+chain_ID + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*11 + str(atom_type[l])+ '\n')
					else:
						fobj.write('ATOM' +  str(total_atoms).rjust(7,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + ' '+chain_ID + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*11 + str(atom_type[l])+ '\n')
					del coordinate[:]
					total_atoms+=1
			if k == 'Hao':
				resid_num +=3
			else:
				resid_num +=1
		total_atoms2 = 0
		
		for start,res in enumerate(sequence):
			atoms2 = writePdbAA.get(res)['atoms']
			if res not in amino_acids:
				print('Looking for residue')
				print(res)
				link_rec = writePdbAA.get(res)['link_rec']
				for x,y in enumerate(link_rec):
					if len(y) > 2:
						link1 = y[0] + total_atoms2
						link2 = y[1] + total_atoms2
						link3 = y[2] + total_atoms2
						fobj.write('CONECT' + str(link1).rjust(5,' ') + str(link2).rjust(5, ' ') + str(link3).rjust(5, ' ') +'\n')
					else:
						link1 = y[0] + total_atoms2
						link2 = y[1] + total_atoms2
						fobj.write('CONECT' + str(link1).rjust(5,' ') + str(link2).rjust(5, ' ') + '\n')
				fobj.write('CONECT' + str(total_atoms2 + 3).rjust(5,' ') + str(total_atoms2 + len(atoms2)+1).rjust(5, ' ') + '\n')
				total_atoms2 = total_atoms2 + len(atoms2)
				print(total_atoms2)
			else:
				total_atoms2 = total_atoms2 + len(atoms2)
		fobj.write('END')		
		fobj.close()
		
	def write_alpha_helix(self,sequence, start_residue, chain_ID):
		fname = asksaveasfilename(filetypes=(("pdb", "*.pdb"),
                                           ("All files", "*.*") ),
								  defaultextension='.pdb')
		try:							   
			fobj = open(fname, 'w')
		except:
				pass
		print("Sequence file to be converted is")
		print(sequence)
		writePdbAA.update(ipn.get_peptoid_coordinates())
		amino_acids = ['Ac','A','R','N','D','C','E','N','G','H','I','L','K','M','F','P','S','T','W', 'Y', 'Q','V',
						'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Ser', 'Thr', 'Trp', 'Tyr', 'Gln','Val']
		locations = []
		locations_temp=[]
		coordinate=[]
		coordinate2 = []
		temp_coordinate=[]
		total_atoms = 1
		residue_count = 0
		turn_number = 7.5
		if start_residue == '':
			print("Empty residue")
			start_residue = 1
		else:
			pass
		if chain_ID == '':
			chain_ID = 'A'
		else:
			pass
		for i,k in enumerate(sequence):
			sr = int(start_residue) + i
			res_name = writePdbAA.get(k)['Res_name']
			
			print(k)
			
			atoms = writePdbAA.get(k)['atoms']
			
			atom_type = writePdbAA.get(k)['atom_type']
			
			for l,m in enumerate(atoms):
				coordinate=[0.00,0.00,0.00]
				if m == 'N':
					coordinate[0] = math.cos((720/turn_number)*residue_count)*1.9
					coordinate[1] = math.sin((720/turn_number)*residue_count)*1.9
					coordinate[2] = 0.000 + 1.400*residue_count
				elif m == 'CA' or m =='CH3':
					coordinate[0] = math.cos((720/turn_number)*residue_count+384)*2.2
					coordinate[1] = math.sin((720/turn_number)*residue_count+384)*2.2
					coordinate[2] = 1.300 + 1.400*residue_count
				elif m == 'C':
					coordinate[0] = math.cos((720/turn_number)*residue_count+818)*1.9
					coordinate[1] = math.sin((720/turn_number)*residue_count+818)*1.9
					coordinate[2] = 2.400 + 1.400*residue_count
				elif m == 'O':
					coordinate[0] = math.cos((720/turn_number)*residue_count+818)*1.9
					coordinate[1] = math.sin((720/turn_number)*residue_count+818)*1.9
					coordinate[2] = 3.800 + 1.400*residue_count
				elif m == 'CB':
					coordinate[0] = math.cos((720/turn_number)*residue_count+384)*3.4
					coordinate[1] = math.sin((720/turn_number)*residue_count+384)*3.4
					coordinate[2] = 0.900 + 1.400*residue_count
				elif m == 'CG' or m =='SG':
					if k == 'P':
						coordinate[0] = math.cos((720/turn_number)*residue_count+390)*4.0
						coordinate[1] = math.sin((720/turn_number)*residue_count+390)*4.0
						coordinate[2] = 0.200 + 1.400*residue_count
					else:
						coordinate[0] = math.cos((720/turn_number)*residue_count+384)*4.6
						coordinate[1] = math.sin((720/turn_number)*residue_count+384)*4.6
						coordinate[2] = 0.900 + 1.400*residue_count
				elif m == 'CG1' or m == 'OG' or m == 'OG1':
					coordinate[0] = math.cos((720/turn_number)*residue_count+378)*4.2
					coordinate[1] = math.sin((720/turn_number)*residue_count+378)*4.2
					coordinate[2] = 2.000 + 1.400*residue_count
				elif m == 'CG2':
					coordinate[0] = math.cos((720/turn_number)*residue_count+378)*4.2
					coordinate[1] = math.sin((720/turn_number)*residue_count+378)*4.2
					coordinate[2] = -0.200 + 1.400*residue_count
				elif m == 'CN':
					coordinate[0] = math.cos((720/turn_number)*residue_count)*3.4
					coordinate[1] = math.sin((720/turn_number)*residue_count)*3.4
					coordinate[2] = -0.200 + 1.400*residue_count
				elif m == 'CD1' or m == 'CD' or m == 'OD1':
					if k == 'I':
						coordinate[0] = math.cos((720/turn_number)*residue_count+378)*5.6
						coordinate[1] = math.sin((720/turn_number)*residue_count+378)*5.6
						coordinate[2] = 2.000 + 1.400*residue_count
					elif k == 'P':
						coordinate[0] = math.cos((720/turn_number)*residue_count)*3.4
						coordinate[1] = math.sin((720/turn_number)*residue_count)*3.4
						coordinate[2] = -0.200 + 1.400*residue_count
					elif k =='W':
						coordinate[0] = math.cos((720/turn_number)*residue_count+384)*5.8
						coordinate[1] = math.sin((720/turn_number)*residue_count+384)*5.8
						coordinate[2] = 2.100 + 1.400*residue_count
					else:
						coordinate[0] = math.cos((720/turn_number)*residue_count+384)*5.8
						coordinate[1] = math.sin((720/turn_number)*residue_count+384)*5.8
						coordinate[2] = -0.200 + 1.400*residue_count
				elif m == 'F2':
					coordinate[0] = math.cos((720/turn_number)*residue_count+384)*5.8
					coordinate[1] = math.sin((720/turn_number)*residue_count+384)*5.8
					coordinate[2] = -1.600 + 1.400*residue_count
				elif m == 'ND1' or m == 'ND2' or m =='CD2' or m =='OD2' or m == 'SD':
					if k == 'W':
						coordinate[0] = math.cos((720/turn_number)*residue_count+384)*5.8
						coordinate[1] = math.sin((720/turn_number)*residue_count+384)*5.8
						coordinate[2] = -0.200 + 1.400*residue_count
					else:
						coordinate[0] = math.cos((720/turn_number)*residue_count+384)*5.8
						coordinate[1] = math.sin((720/turn_number)*residue_count+384)*5.8
						coordinate[2] = 2.100 + 1.400*residue_count
				elif m == 'F1':
					coordinate[0] = math.cos((720/turn_number)*residue_count+384)*5.8
					coordinate[1] = math.sin((720/turn_number)*residue_count+384)*5.8
					coordinate[2] = 3.500 + 1.400*residue_count
				elif m == 'NE2' or m =='NE' or m == 'OE2' or m == 'CE2':
					if k == 'F' or k == 'Y' or k == 'Phe(I)' or k == 'Phe(Br)' or k == 'Phe(5F)':
						coordinate[0] = math.cos((720/turn_number)*residue_count+384)*7.2
						coordinate[1] = math.sin((720/turn_number)*residue_count+384)*7.2
						coordinate[2] = 2.100 + 1.400*residue_count
					elif k == 'W':
						coordinate[0] = math.cos((720/turn_number)*residue_count+384)*7.2
						coordinate[1] = math.sin((720/turn_number)*residue_count+384)*7.2
						coordinate[2] = 0.100 + 1.400*residue_count
					else:
						coordinate[0] = math.cos((720/turn_number)*residue_count+384)*7.2
						coordinate[1] = math.sin((720/turn_number)*residue_count+384)*7.2
						coordinate[2] = 0.000 + 1.400*residue_count
				elif m == 'F3':
					coordinate[0] = math.cos((720/turn_number)*residue_count+384)*7.2
					coordinate[1] = math.sin((720/turn_number)*residue_count+384)*7.2
					coordinate[2] = 3.500 + 1.400*residue_count
				elif m == 'OE1':
					coordinate[0] = math.cos((720/turn_number)*residue_count+384)*5.8
					coordinate[1] = math.sin((720/turn_number)*residue_count+384)*5.8
					coordinate[2] = -1.600 + 1.400*residue_count
				elif m == 'CE1' or m == 'CE' or m == 'NE1':
					if k == 'K' or k == 'M':
						coordinate[0] = math.cos((720/turn_number)*residue_count+384)*7.2
						coordinate[1] = math.sin((720/turn_number)*residue_count+384)*7.2
						coordinate[2] = 1.000 + 1.400*residue_count
					elif k == 'W':
						coordinate[0] = math.cos((720/turn_number)*residue_count+384)*7.2
						coordinate[1] = math.sin((720/turn_number)*residue_count+384)*7.2
						coordinate[2] = 1.700 + 1.400*residue_count
					elif k == 'F' or k == 'Y' or k == 'Phe(I)' or k == 'Phe(Br)' or k == 'Phe(5F)':
						coordinate[0] = math.cos((720/turn_number)*residue_count+384)*7.2
						coordinate[1] = math.sin((720/turn_number)*residue_count+384)*7.2
						coordinate[2] = -0.200 + 1.400*residue_count
					else:
						coordinate[0] = math.cos((720/turn_number)*residue_count+384)*7.2
						coordinate[1] = math.sin((720/turn_number)*residue_count+384)*7.2
						coordinate[2] = 1.900 + 1.400*residue_count
				elif m == 'F4':
					coordinate[0] = math.cos((720/turn_number)*residue_count+384)*7.2
					coordinate[1] = math.sin((720/turn_number)*residue_count+384)*7.2
					coordinate[2] = -1.600 + 1.400*residue_count
				elif m == 'CE3':
					coordinate[0] = math.cos((720/turn_number)*residue_count+384)*5.0
					coordinate[1] = math.sin((720/turn_number)*residue_count+384)*5.0
					coordinate[2] = -1.500 + 1.400*residue_count
				elif m == 'CZ2':
					coordinate[0] = math.cos((720/turn_number)*residue_count+384)*8.2
					coordinate[1] = math.sin((720/turn_number)*residue_count+384)*8.2
					coordinate[2] = -1.500 + 1.400*residue_count
				elif m == 'F5':
					coordinate[0] = math.cos((720/turn_number)*residue_count+384)*9.6
					coordinate[1] = math.sin((720/turn_number)*residue_count+384)*9.6
					coordinate[2] = 0.900 + 1.400*residue_count
				elif m == 'CZ3':
					coordinate[0] = math.cos((720/turn_number)*residue_count+384)*5.5
					coordinate[1] = math.sin((720/turn_number)*residue_count+384)*5.5
					coordinate[2] = -2.700 + 1.400*residue_count
				elif m == 'CH2':
					coordinate[0] = math.cos((720/turn_number)*residue_count+384)*7.4
					coordinate[1] = math.sin((720/turn_number)*residue_count+384)*7.4
					coordinate[2] = -2.700 + 1.400*residue_count
				elif m == 'CZ' or m == 'NZ':
					coordinate[0] = math.cos((720/turn_number)*residue_count+384)*8.4
					coordinate[1] = math.sin((720/turn_number)*residue_count+384)*8.4
					coordinate[2] = 0.900 + 1.400*residue_count
				elif m == 'NH1':
					coordinate[0] = math.cos((720/turn_number)*residue_count+384)*10.0
					coordinate[1] = math.sin((720/turn_number)*residue_count+384)*10.0
					coordinate[2] = -0.200 + 1.400*residue_count
				elif m == 'NH2':
					coordinate[0] = math.cos((720/turn_number)*residue_count+384)*10.0
					coordinate[1] = math.sin((720/turn_number)*residue_count+384)*10.0
					coordinate[2] = 2.000 + 1.400*residue_count
				elif m == 'OH' or m == 'BR' or m == 'I':
					coordinate[0] = math.cos((720/turn_number)*residue_count+384)*9.6
					coordinate[1] = math.sin((720/turn_number)*residue_count+384)*9.6
					coordinate[2] = 1.000 + 1.400*residue_count
				else:
					coordinate[1] = coordinate[1]*(-1) + 0.870
					coordinate[2] = coordinate[2]*(-1)
					coordinate[0] = coordinate[0] + i*3.675
				print(coordinate)
				if atom_type[l] == 'Br':
					fobj.write('HETATM' + str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + ' '+chain_ID + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*8 + str(atom_type[l]).rjust(4,' ')+ '\n')
				elif k not in amino_acids:
					fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + ' '+chain_ID + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*11 + str(atom_type[l])+ '\n')
				else:
					fobj.write('ATOM' +  str(total_atoms).rjust(7,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + ' '+chain_ID + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*11 + str(atom_type[l])+ '\n')
				del coordinate[:]
				total_atoms+=1
			residue_count +=1
		total_atoms2 = 0
		for start,res in enumerate(sequence):
			atoms2 = writePdbAA.get(res)['atoms']
			if res not in amino_acids:
				print('Looking for residue')
				print(res)
				link_rec = writePdbAA.get(res)['link_rec']
				for x,y in enumerate(link_rec):
					if len(y) > 2:
						link1 = y[0] + total_atoms2
						link2 = y[1] + total_atoms2
						link3 = y[2] + total_atoms2
						fobj.write('CONECT' + str(link1).rjust(5,' ') + str(link2).rjust(5, ' ') + str(link3).rjust(5, ' ') +'\n')
					else:
						link1 = y[0] + total_atoms2
						link2 = y[1] + total_atoms2
						fobj.write('CONECT' + str(link1).rjust(5,' ') + str(link2).rjust(5, ' ') + '\n')
				fobj.write('CONECT' + str(total_atoms2 + 3).rjust(5,' ') + str(total_atoms2 + len(atoms2)+1).rjust(5, ' ') + '\n')
				total_atoms2 = total_atoms2 + len(atoms2)
				print(total_atoms2)
			else:
				total_atoms2 = total_atoms2 + len(atoms2)
		fobj.write('END')		
		fobj.close()
	
	def write_totally_tubular_par(self,sequence):
		fname = asksaveasfilename(filetypes=(("pdb", "*.pdb"),
                                           ("All files", "*.*") ),
								  defaultextension='.pdb')
		try:							   
			fobj = open(fname, 'w')
		except:
				pass
		print(len(sequence))
		sequence2 = sequence[1:]
		writePdbAA.update(ipn.get_peptoid_coordinates())
		start_pos=0
		sr = 1
		Z_spacing = 26
		number_of_rings = 8
		repeats = 3
		l=0
		for q in range(0,number_of_rings):
			if q%2 == 0:
				start_pos =0
				degrees = (360/((len(sequence)-1)*3*repeats+repeats*2))
			else:
				start_pos = ((len(sequence)-1)*3)/2
				degrees = (-360/((len(sequence)-1)*3*repeats+repeats*2))
			for p in range (0, repeats):
				radius = (3.1 * (len(sequence)-1)*repeats+repeats*2)/(2*math.pi)
				
				#radians = (degrees*math.pi*i)/180
				print(radius)
				print(degrees)
				chain_ID = chain[q+p+q*3]
				
				for i,k in enumerate(sequence2):
					#sr = int(start_residue) + i
					res_name = writePdbAA.get(k)['Res_name']
			
					atoms = writePdbAA.get(k)['atoms']
			
					atom_type = writePdbAA.get(k)['atom_type']
					total_atoms =0
					 
					
					if i%2 == 0:
						for l,m in enumerate(atoms):
							coordinate=[0.00,0.00,0.00]
							if m == 'N':
								coordinate[0] = math.cos((degrees*math.pi*start_pos)/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*start_pos)/180)*radius
								coordinate[2] = 0.000+Z_spacing*q
							elif m == 'CA' or m =='CH3':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos+1))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos+1))/180)*radius
								coordinate[2] = 1.200+Z_spacing*q
							elif m == 'C':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos+2))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos+2))/180)*radius
								coordinate[2] = 0.000+Z_spacing*q
							elif m == 'O':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos+2))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos+2))/180)*radius
								coordinate[2] = -1.400+Z_spacing*q
							elif m == 'CA5':
								coordinate[0] = math.cos((degrees*math.pi*start_pos)/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*start_pos)/180)*radius
								coordinate[2] = -1.400+Z_spacing*q
							elif m == 'CB':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-1))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-1))/180)*radius
								coordinate[2] = -2.200+Z_spacing*q
							elif m == 'CG' or m == 'NG' or m == 'OG':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-1))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-1))/180)*radius
								coordinate[2] = -3.600+Z_spacing*q
							elif m == 'CD':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-2))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-2))/180)*radius
								coordinate[2] = -4.400+Z_spacing*q
							elif m == 'CE':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-2))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-2))/180)*radius
								coordinate[2] = -5.800+Z_spacing*q
							elif m == 'CZ' or m == 'OZ':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-3))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-3))/180)*radius
								coordinate[2] = -6.600+Z_spacing*q
							elif m == 'CH':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-3))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-3))/180)*radius
								coordinate[2] = -8.100+Z_spacing*q
							elif m == 'CT':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-4))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-4))/180)*radius
								coordinate[2] = -8.900+Z_spacing*q
							elif m == 'CI' or m =='OI':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-4))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-4))/180)*radius
								coordinate[2] = -10.400+Z_spacing*q
							elif m == 'CK':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-5))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-5))/180)*radius
								coordinate[2] = -11.200+Z_spacing*q
							else:
								pass
							fobj.write('ATOM' +  str(total_atoms).rjust(7,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() +chain_ID.ljust(3,' ') + str(sr).rjust(3,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*11 + str(atom_type[l])+ '\n')
							del coordinate[:]
							total_atoms+=1
					else:
						for l,m in enumerate(atoms):
							coordinate=[0.00,0.00,0.00]
							if m == 'N':
								coordinate[0] = math.cos((degrees*math.pi*start_pos)/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*start_pos)/180)*radius
								coordinate[2] = 1.200+Z_spacing*q
							elif m == 'CA' or m =='CH3':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos+1))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos+1))/180)*radius
								coordinate[2] = 0.000+Z_spacing*q
							elif m == 'C':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos+2))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos+2))/180)*radius
								coordinate[2] = 1.200+Z_spacing*q
							elif m == 'O':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos+2))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos+2))/180)*radius
								coordinate[2] = 2.600+Z_spacing*q
							elif m == 'CA5':
								coordinate[0] = math.cos((degrees*math.pi*start_pos)/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*start_pos)/180)*radius
								coordinate[2] = 2.200+Z_spacing*q
							elif m == 'CB':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-1))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-1))/180)*radius
								coordinate[2] = 3.000+Z_spacing*q
							elif m == 'CG' or m == 'NG' or m == 'OG':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-1))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-1))/180)*radius
								coordinate[2] =  4.500+Z_spacing*q
							elif m == 'CD':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-2))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-2))/180)*radius
								coordinate[2] = 5.300+Z_spacing*q
							elif m == 'CE':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-2))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-2))/180)*radius
								coordinate[2] = 6.800+Z_spacing*q
							elif m == 'CZ' or m == 'OZ':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-3))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-3))/180)*radius
								coordinate[2] = 7.600+Z_spacing*q
							elif m == 'CH':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-3))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-3))/180)*radius
								coordinate[2] = 9.100+Z_spacing*q
							elif m == 'CT':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-4))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-4))/180)*radius
								coordinate[2] = 9.900+Z_spacing*q
							elif m == 'CI' or m =='OI':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-4))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-4))/180)*radius
								coordinate[2] = 11.400+Z_spacing*q
							elif m == 'CK':
								coordinate[0] = math.cos((degrees*math.pi*(start_pos-5))/180)*radius
								coordinate[1] = math.sin((degrees*math.pi*(start_pos-5))/180)*radius
								coordinate[2] = 12.200+Z_spacing*q
							else:
								pass
							fobj.write('ATOM' +  str(total_atoms).rjust(7,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() +chain_ID.ljust(3,' ') + str(sr).rjust(3,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*11 + str(atom_type[l])+ '\n')
							del coordinate[:]
							total_atoms+=1
					start_pos +=3
					sr +=1
				start_pos +=2
			l+=1
		fobj.write('END')
		fobj.close()
			#residue_count +=1
			
	def write_totally_tubular_slip(self,sequence):
		fname = asksaveasfilename(filetypes=(("pdb", "*.pdb"),
                                           ("All files", "*.*") ),
								  defaultextension='.pdb')
		try:							   
			fobj = open(fname, 'w')
		except:
				pass
		print(len(sequence))
		sequence2 = sequence[1:]
		writePdbAA.update(ipn.get_peptoid_coordinates())
		start_pos=0
		sr = 1
		Z_spacing = 26
		number_of_rings = 10
		repeats = 6
		compression = -80
		total_atoms=1
		tilt=0
		for q in range(0,number_of_rings):
			if q%2 == 0:
				ring_orientation = 1
			else:
				ring_orientation = 1
			for p in range (0, repeats):
				sr = 1
				radius = (3.1 * (len(sequence)-1)*(repeats/2)+(repeats/2))/(2*math.pi)
				chain_ID = chain[q+p+q*5]
				rotation_x = -math.sin(((60+p*60)*math.pi)/180)*7 + math.cos(((0+p*60)*math.pi)/180)*7
				rotation_y = -math.cos(((60+p*60)*math.pi)/180)*7 - math.sin(((0+p*60)*math.pi)/180)*7
				if p%2 == 0:
					start_pos = p*(len(sequence)*3/2)
					degrees = (360/((len(sequence)-1)*3*(repeats/2)))
					orientation = 1
					for i,k in enumerate(sequence2):
						res_name = writePdbAA.get(k)['Res_name']
						atoms = writePdbAA.get(k)['atoms']
						atom_type = writePdbAA.get(k)['atom_type']
						if i%2 == 0:
							for l,m in enumerate(atoms):
								coordinate=[0.00,0.00,0.00]
								if m == 'N':
									coordinate[0] = math.sin((degrees*math.pi*start_pos)/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*start_pos)/180)*radius+rotation_y
									coordinate[2] = 0.000+Z_spacing*q
								elif m == 'CA' or m =='CH3':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+1))/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+1))/180)*radius+rotation_y
									coordinate[2] = 1.200+Z_spacing*q
								elif m == 'C':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+2))/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+2))/180)*radius+rotation_y
									coordinate[2] = 0.000+Z_spacing*q
								elif m == 'O':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+2))/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+2))/180)*radius+rotation_y
									coordinate[2] = -1.400+Z_spacing*q
								elif m == 'CA5':
									coordinate[0] = math.sin((degrees*math.pi*start_pos)/180)*(radius+tilt*1)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*start_pos)/180)*(radius+tilt*1)+rotation_y
									coordinate[2] = -1.400+Z_spacing*q
								elif m == 'CB':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+tilt*2)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+tilt*2)+rotation_y
									coordinate[2] = -2.200+Z_spacing*q
								elif m == 'CG' or m == 'NG' or m == 'OG':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+tilt*3)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+tilt*3)+rotation_y
									coordinate[2] = -3.600+Z_spacing*q
								elif m == 'CD':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+tilt*4)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+tilt*4)+rotation_y
									coordinate[2] = -4.400+Z_spacing*q
								elif m == 'CE':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+tilt*5)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+tilt*5)+rotation_y
									coordinate[2] = -5.800+Z_spacing*q
								elif m == 'CZ' or m == 'OZ':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+tilt*6)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+tilt*6)+rotation_y
									coordinate[2] = -6.600+Z_spacing*q
								elif m == 'CH':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+tilt*7)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+tilt*7)+rotation_y
									coordinate[2] = -8.100+Z_spacing*q
								elif m == 'CT':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+tilt*8)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+tilt*8)+rotation_y
									coordinate[2] = -8.900+Z_spacing*q
								elif m == 'CI' or m =='OI':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+tilt*9)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+tilt*9)+rotation_y
									coordinate[2] = -10.400+Z_spacing*q
								elif m == 'CK':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-5*orientation*ring_orientation))/180)*(radius+tilt*10)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-5*orientation*ring_orientation))/180)*(radius+tilt*10)+rotation_y
									coordinate[2] = -11.200+Z_spacing*q
								else:
									pass
								fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() +' '+chain_ID.ljust(2,' ') + str(sr).rjust(3,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*11 + str(atom_type[l])+ '\n')
								del coordinate[:]
								total_atoms+=1
							sr +=1
						else:
							for l,m in enumerate(atoms):
								
								coordinate=[0.00,0.00,0.00]
								if m == 'N':
									coordinate[0] = math.sin((degrees*math.pi*start_pos)/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*start_pos)/180)*radius+rotation_y
									coordinate[2] = 1.200+Z_spacing*q 
								elif m == 'CA' or m =='CH3':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+1))/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+1))/180)*radius+rotation_y
									coordinate[2] = 0.000+Z_spacing*q
								elif m == 'C':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+2))/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+2))/180)*radius+rotation_y
									coordinate[2] = 1.200+Z_spacing*q
								elif m == 'O':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+2))/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+2))/180)*radius+rotation_y
									coordinate[2] = 2.600+Z_spacing*q
								elif m == 'CA5':
									coordinate[0] = math.sin((degrees*math.pi*start_pos)/180)*(radius+tilt*-1)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*start_pos)/180)*(radius+tilt*-1)+rotation_y
									coordinate[2] = 2.200+Z_spacing*q
								elif m == 'CB':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+tilt*-2)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+tilt*-2)+rotation_y
									coordinate[2] = 3.000+Z_spacing*q
								elif m == 'CG' or m == 'NG' or m == 'OG':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+tilt*-3)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+tilt*-3)+rotation_y
									coordinate[2] =  4.500+Z_spacing*q
								elif m == 'CD':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+tilt*-4)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+tilt*-4)+rotation_y
									coordinate[2] = 5.300+Z_spacing*q
								elif m == 'CE':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+tilt*-5)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+tilt*-5)+rotation_y
									coordinate[2] = 6.800+Z_spacing*q
								elif m == 'CZ' or m == 'OZ':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+tilt*-6)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+tilt*-6)+rotation_y
									coordinate[2] = 7.600+Z_spacing*q
								elif m == 'CH':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+tilt*-7)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+tilt*-7)+rotation_y
									coordinate[2] = 9.100+Z_spacing*q
								elif m == 'CT':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+tilt*-8)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+tilt*-8)+rotation_y
									coordinate[2] = 9.900+Z_spacing*q
								elif m == 'CI' or m =='OI':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+tilt*-9)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+tilt*-9)+rotation_y
									coordinate[2] = 11.400+Z_spacing*q
								elif m == 'CK':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-5*orientation*ring_orientation))/180)*(radius+tilt*-10)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-5*orientation*ring_orientation))/180)*(radius+tilt*-10)+rotation_y
									coordinate[2] = 12.200+Z_spacing*q
								else:
									pass
								fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() +' '+chain_ID.ljust(2,' ') + str(sr).rjust(3,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*11 + str(atom_type[l])+ '\n')
								del coordinate[:]
								total_atoms+=1
							sr +=1
						start_pos +=3
					start_pos = start_pos/2
					fobj.write('TER\n')
				else:
					start_pos = p*((len(sequence)-1)*3)+(len(sequence)*3/2)
					degrees = (-360/((len(sequence)-1)*3*(repeats/2)))
					#start_pos = p*(len(sequence)*3/2)
					#degrees = (360/((len(sequence)-1)*3*(repeats/2)))
					orientation = 1
					slip =0.5*Z_spacing
					for i,k in enumerate(sequence2):
						res_name = writePdbAA.get(k)['Res_name']
						atoms = writePdbAA.get(k)['atoms']
						atom_type = writePdbAA.get(k)['atom_type']
						if i%2 == 0:
							for l,m in enumerate(atoms):
								coordinate=[0.00,0.00,0.00]
								if m == 'N':
									coordinate[0] = math.sin((degrees*math.pi*start_pos)/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*start_pos)/180)*radius+rotation_y
									coordinate[2] = 1.200+Z_spacing*q+slip
								elif m == 'CA' or m =='CH3':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+1))/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+1))/180)*radius+rotation_y
									coordinate[2] = 0.000+Z_spacing*q+slip
								elif m == 'C':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+2))/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+2))/180)*radius+rotation_y
									coordinate[2] = 1.200+Z_spacing*q+slip
								elif m == 'O':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+2))/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+2))/180)*radius+rotation_y
									coordinate[2] = 2.600+Z_spacing*q+slip
								elif m == 'CA5':
									coordinate[0] = math.sin((degrees*math.pi*start_pos)/180)*(radius+tilt*-1)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*start_pos)/180)*(radius+tilt*-1)+rotation_y
									coordinate[2] = 2.200+Z_spacing*q+slip
								elif m == 'CB':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+tilt*-2)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+tilt*-2)+rotation_y
									coordinate[2] = 3.000+Z_spacing*q+slip
								elif m == 'CG' or m == 'NG' or m == 'OG':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+tilt*-3)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+tilt*-3)+rotation_y
									coordinate[2] =  4.500+Z_spacing*q+slip
								elif m == 'CD':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+tilt*-4)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+tilt*-4)+rotation_y
									coordinate[2] = 5.300+Z_spacing*q+slip
								elif m == 'CE':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+tilt*-5)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+tilt*-5)+rotation_y
									coordinate[2] = 6.800+Z_spacing*q+slip
								elif m == 'CZ' or m == 'OZ':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+tilt*-6)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+tilt*-6)+rotation_y
									coordinate[2] = 7.600+Z_spacing*q+slip
								elif m == 'CH':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+tilt*-7)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+tilt*-7)+rotation_y
									coordinate[2] = 9.100+Z_spacing*q+slip
								elif m == 'CT':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+tilt*-8)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+tilt*-8)+rotation_y
									coordinate[2] = 9.900+Z_spacing*q+slip
								elif m == 'CI' or m =='OI':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+tilt*-9)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+tilt*-9)+rotation_y
									coordinate[2] = 11.400+Z_spacing*q+slip
								elif m == 'CK':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-5*orientation*ring_orientation))/180)*(radius+tilt*-10)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-5*orientation*ring_orientation))/180)*(radius+tilt*-10)+rotation_y
									coordinate[2] = 12.200+Z_spacing*q+slip
								else:
									pass
								
								fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() +' '+chain_ID.ljust(2,' ') + str(sr).rjust(3,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*11 + str(atom_type[l])+ '\n')
								del coordinate[:]
								total_atoms+=1
							sr +=1
						else:
							for l,m in enumerate(atoms):
								coordinate=[0.00,0.00,0.00]
								if m == 'N':
									coordinate[0] = math.sin((degrees*math.pi*start_pos)/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*start_pos)/180)*radius+rotation_y
									coordinate[2] = 0.000+Z_spacing*q+slip
								elif m == 'CA' or m =='CH3':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+1))/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+1))/180)*radius+rotation_y
									coordinate[2] = 1.200+Z_spacing*q+slip
								elif m == 'C':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+2))/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+2))/180)*radius+rotation_y
									coordinate[2] = 0.000+Z_spacing*q+slip
								elif m == 'O':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+2))/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+2))/180)*radius+rotation_y
									coordinate[2] = -1.400+Z_spacing*q+slip
								elif m == 'CA5':
									coordinate[0] = math.sin((degrees*math.pi*start_pos)/180)*(radius+tilt*1)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*start_pos)/180)*(radius+tilt*1)+rotation_y
									coordinate[2] = -1.400+Z_spacing*q+slip
								elif m == 'CB':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+tilt*2)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+tilt*2)+rotation_y
									coordinate[2] = -2.200+Z_spacing*q+slip
								elif m == 'CG' or m == 'NG' or m == 'OG':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+tilt*3)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+tilt*3)+rotation_y
									coordinate[2] = -3.600+Z_spacing*q+slip
								elif m == 'CD':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+tilt*4)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+tilt*4)+rotation_y
									coordinate[2] = -4.400+Z_spacing*q+slip
								elif m == 'CE':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+tilt*5)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+tilt*5)+rotation_y
									coordinate[2] = -5.800+Z_spacing*q+slip
								elif m == 'CZ' or m == 'OZ':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+tilt*6)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+tilt*6)+rotation_y
									coordinate[2] = -6.600+Z_spacing*q+slip
								elif m == 'CH':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+tilt*7)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+tilt*7)+rotation_y
									coordinate[2] = -8.100+Z_spacing*q+slip
								elif m == 'CT':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+tilt*8)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+tilt*8)+rotation_y
									coordinate[2] = -8.900+Z_spacing*q+slip
								elif m == 'CI' or m =='OI':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+tilt*9)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+tilt*9)+rotation_y
									coordinate[2] = -10.400+Z_spacing*q+slip
								elif m == 'CK':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-5*orientation*ring_orientation))/180)*(radius+tilt*10)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-5*orientation*ring_orientation))/180)*(radius+tilt*10)+rotation_y
									coordinate[2] = -11.200+Z_spacing*q+slip
								else:
									pass
								fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() +' '+chain_ID.ljust(2,' ') + str(sr).rjust(3,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*11 + str(atom_type[l])+ '\n')
								del coordinate[:]
								total_atoms+=1
							sr +=1
						start_pos +=3
					start_pos = start_pos/2
					fobj.write('TER\n')
				
		total_atoms2 = 0
		for q in range(0,number_of_rings):
			
			for p in range (0, repeats):
				link = 0
				for start,res in enumerate(sequence2):
					atoms2 = writePdbAA.get(res)['atoms']
					print('Looking for residue')
					print(res)
					link_rec = writePdbAA.get(res)['link_rec']
					for x,y in enumerate(link_rec):
						if len(y) > 2:
							link1 = y[0] + total_atoms2
							link2 = y[1] + total_atoms2
							link3 = y[2] + total_atoms2
							fobj.write('CONECT' + str(link1).rjust(5,' ') + str(link2).rjust(5, ' ') + str(link3).rjust(5, ' ') +'\n')
						else:
							link1 = y[0] + total_atoms2
							link2 = y[1] + total_atoms2
							fobj.write('CONECT' + str(link1).rjust(5,' ') + str(link2).rjust(5, ' ') + '\n')
					if link < (len(sequence)-2):
						fobj.write('CONECT' + str(total_atoms2 + 3).rjust(5,' ') + str(total_atoms2 + len(atoms2)+1).rjust(5, ' ') + '\n')
						link +=1
					else:
						pass
					total_atoms2 = total_atoms2 + len(atoms2)
					print(total_atoms2)
					
		print(sr)
		fobj.write('END')
		fobj.close()
		
		
	def write_totally_tubular_side(self,sequence):
		fname = asksaveasfilename(filetypes=(("pdb", "*.pdb"),
                                           ("All files", "*.*") ),
								  defaultextension='.pdb')
		try:							   
			fobj = open(fname, 'w')
		except:
				pass
		print(len(sequence))
		sequence2 = sequence[1:]
		writePdbAA.update(ipn.get_peptoid_coordinates())
		start_pos=0
		sr = 1
		Z_spacing = 10
		number_of_rings = 20
		repeats = 6
		compression = -80
		total_atoms=1
		for q in range(0,number_of_rings):
			if q%2 == 0:
				ring_orientation = 1
			else:
				ring_orientation = 1
			for p in range (0, repeats):
				print(p)
				sr = 1
				radius = (3.1 * (len(sequence)-1)*(repeats/2)+(repeats/2))/(2*math.pi)
				chain_ID = chain[q+p+q*5]
				#rotation_x = -math.sin(((60+p*60)*math.pi)/180)*7 + math.cos(((0+p*60)*math.pi)/180)*7
				#rotation_y = -math.cos(((60+p*60)*math.pi)/180)*7 - math.sin(((0+p*60)*math.pi)/180)*7
				rotation_x=0
				rotation_y=0
				if p%2 == 0:
					start_pos = p*(len(sequence)*3/2)
					degrees = (360/((len(sequence)-1)*3*(repeats/2)))
					orientation = 1
					for i,k in enumerate(sequence2):
						res_name = writePdbAA.get(k)['Res_name']
						atoms = writePdbAA.get(k)['atoms']
						atom_type = writePdbAA.get(k)['atom_type']
						if i%2 == 0:
							for l,m in enumerate(atoms):
								coordinate=[0.00,0.00,0.00]
								if m == 'N':
									coordinate[0] = math.sin((degrees*math.pi*start_pos)/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*start_pos)/180)*radius+rotation_y
									coordinate[2] = 0.000+Z_spacing*q
								elif m == 'CA' or m =='CH3':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+1))/180)*(radius+1.200)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+1))/180)*(radius+1.200)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'C':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+2))/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+2))/180)*radius+rotation_y
									coordinate[2] = 0.000+Z_spacing*q
								elif m == 'O':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+2))/180)*(radius-1.400)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+2))/180)*(radius-1.400)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CA5':
									coordinate[0] = math.sin((degrees*math.pi*start_pos)/180)*(radius-1.400)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*start_pos)/180)*(radius-1.400)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CB':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius-2.200)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius-2.200)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CG' or m == 'NG' or m == 'OG':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius-3.600)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius-3.600)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CD':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius-4.400)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius-4.400)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CE':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius-5.800)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius-5.800)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CZ' or m == 'OZ':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius-6.600)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius-6.600)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CH':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius-8.100)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius-8.100)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CT':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius-8.900)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius-8.900)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CI' or m =='OI':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius-10.400)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius-10.400)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CK':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-5*orientation*ring_orientation))/180)*(radius-11.200)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-5*orientation*ring_orientation))/180)*(radius-11.200)+rotation_y
									coordinate[2] = Z_spacing*q
								else:
									pass
								fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() +' '+chain_ID.ljust(2,' ') + str(sr).rjust(3,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*11 + str(atom_type[l])+ '\n')
								del coordinate[:]
								total_atoms+=1
							sr +=1
						else:
							for l,m in enumerate(atoms):
								
								coordinate=[0.00,0.00,0.00]
								if m == 'N':
									coordinate[0] = math.sin((degrees*math.pi*start_pos)/180)*(radius+1.200)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*start_pos)/180)*(radius+1.200)+rotation_y
									coordinate[2] = Z_spacing*q 
								elif m == 'CA' or m =='CH3':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+1))/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+1))/180)*radius+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'C':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+2))/180)*(radius+1.200)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+2))/180)*(radius+1.200)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'O':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+2))/180)*(radius+2.600)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+2))/180)*(radius+2.600)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CA5':
									coordinate[0] = math.sin((degrees*math.pi*start_pos)/180)*(radius+2.200)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*start_pos)/180)*(radius+2.200)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CB':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+3.000)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+3.000)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CG' or m == 'NG' or m == 'OG':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+4.500)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+4.500)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CD':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+5.300)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+5.300)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CE':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+6.800)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+6.800)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CZ' or m == 'OZ':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+7.600)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+7.600)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CH':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+9.100)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+9.100)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CT':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+9.900)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+9.900)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CI' or m =='OI':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+11.400)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+11.400)+rotation_y
									coordinate[2] = Z_spacing*q
								elif m == 'CK':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-5*orientation*ring_orientation))/180)*(radius+12.200)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-5*orientation*ring_orientation))/180)*(radius+12.200)+rotation_y
									coordinate[2] = Z_spacing*q
								else:
									pass
								fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() +' '+chain_ID.ljust(2,' ') + str(sr).rjust(3,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*11 + str(atom_type[l])+ '\n')
								del coordinate[:]
								total_atoms+=1
							sr +=1
						start_pos +=3
					start_pos = start_pos/2
					fobj.write('TER\n')
				else:
					start_pos = p*((len(sequence)-1)*3)+(len(sequence)*3/2)
					degrees = (-360/((len(sequence)-1)*3*(repeats/2)))
					orientation = 1
					slip =0.5*Z_spacing
					for i,k in enumerate(sequence2):
						res_name = writePdbAA.get(k)['Res_name']
						atoms = writePdbAA.get(k)['atoms']
						atom_type = writePdbAA.get(k)['atom_type']
						if i%2 == 0:
							for l,m in enumerate(atoms):
								coordinate=[0.00,0.00,0.00]
								if m == 'N':
									coordinate[0] = math.sin((degrees*math.pi*start_pos)/180)*(radius+1.200)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*start_pos)/180)*(radius+1.200)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CA' or m =='CH3':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+1))/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+1))/180)*radius+rotation_y
									coordinate[2] = 0.000+Z_spacing*q+slip
								elif m == 'C':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+2))/180)*(radius+1.200)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+2))/180)*(radius+1.200)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'O':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+2))/180)*(radius+2.600)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+2))/180)*(radius+2.600)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CA5':
									coordinate[0] = math.sin((degrees*math.pi*start_pos)/180)*(radius+2.200)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*start_pos)/180)*(radius+2.200)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CB':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+3.000)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+3.000)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CG' or m == 'NG' or m == 'OG':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+4.500)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius+4.500)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CD':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+5.300)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+5.300)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CE':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+6.800)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius+6.800)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CZ' or m == 'OZ':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+7.600)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+7.600)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CH':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+9.100)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius+9.100)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CT':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+9.900)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+9.900)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CI' or m =='OI':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+11.400)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius+11.400)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CK':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-5*orientation*ring_orientation))/180)*(radius+12.200)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-5*orientation*ring_orientation))/180)*(radius+12.200)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								else:
									pass
								
								fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() +' '+chain_ID.ljust(2,' ') + str(sr).rjust(3,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*11 + str(atom_type[l])+ '\n')
								del coordinate[:]
								total_atoms+=1
							sr +=1
						else:
							for l,m in enumerate(atoms):
								coordinate=[0.00,0.00,0.00]
								if m == 'N':
									coordinate[0] = math.sin((degrees*math.pi*start_pos)/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*start_pos)/180)*radius+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CA' or m =='CH3':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+1))/180)*(radius+1.200)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+1))/180)*(radius+1.200)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'C':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+2))/180)*radius+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+2))/180)*radius+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'O':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos+2))/180)*(radius-1.400)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos+2))/180)*(radius-1.400)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CA5':
									coordinate[0] = math.sin((degrees*math.pi*start_pos)/180)*(radius-1.400)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*start_pos)/180)*(radius-1.400)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CB':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius-2.200)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius-2.200)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CG' or m == 'NG' or m == 'OG':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius-3.600)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-1*orientation*ring_orientation))/180)*(radius-3.600)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CD':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius-4.400)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius-4.400)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CE':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius-5.800)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-2*orientation*ring_orientation))/180)*(radius-5.800)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CZ' or m == 'OZ':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius-6.600)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius-6.600)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CH':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius-8.100)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-3*orientation*ring_orientation))/180)*(radius-8.100)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CT':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius-8.900)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius-8.900)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CI' or m =='OI':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius-10.400)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-4*orientation*ring_orientation))/180)*(radius-10.400)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								elif m == 'CK':
									coordinate[0] = math.sin((degrees*math.pi*(start_pos-5*orientation*ring_orientation))/180)*(radius-11.200)+rotation_x
									coordinate[1] = math.cos((degrees*math.pi*(start_pos-5*orientation*ring_orientation))/180)*(radius-11.200)+rotation_y
									coordinate[2] = Z_spacing*q+slip
								else:
									pass
								fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() +' '+chain_ID.ljust(2,' ') + str(sr).rjust(3,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*11 + str(atom_type[l])+ '\n')
								del coordinate[:]
								total_atoms+=1
							sr +=1
						start_pos +=3
					start_pos = start_pos/2
					fobj.write('TER\n')
				
		total_atoms2 = 0
		for q in range(0,number_of_rings):
			
			for p in range (0, repeats):
				link = 0
				for start,res in enumerate(sequence2):
					atoms2 = writePdbAA.get(res)['atoms']
					print('Looking for residue')
					print(res)
					link_rec = writePdbAA.get(res)['link_rec']
					for x,y in enumerate(link_rec):
						if len(y) > 2:
							link1 = y[0] + total_atoms2
							link2 = y[1] + total_atoms2
							link3 = y[2] + total_atoms2
							fobj.write('CONECT' + str(link1).rjust(5,' ') + str(link2).rjust(5, ' ') + str(link3).rjust(5, ' ') +'\n')
						else:
							link1 = y[0] + total_atoms2
							link2 = y[1] + total_atoms2
							fobj.write('CONECT' + str(link1).rjust(5,' ') + str(link2).rjust(5, ' ') + '\n')
					if link < (len(sequence)-2):
						fobj.write('CONECT' + str(total_atoms2 + 3).rjust(5,' ') + str(total_atoms2 + len(atoms2)+1).rjust(5, ' ') + '\n')
						link +=1
					else:
						pass
					total_atoms2 = total_atoms2 + len(atoms2)
					print(total_atoms2)
					
		print(sr)
		fobj.write('END')
		fobj.close()
		
		
	def what_the_sheet_par(self,sequence,length,height,width, l_adjust):
		fname = asksaveasfilename(filetypes=(("pdb", "*.pdb"),
                                           ("All files", "*.*") ),
								  defaultextension='.pdb')
		try:							   
			fobj = open(fname, 'w')
		except:
				pass
		print(len(sequence))
		sequence2 = sequence[1:]
			#residue_count +=1
		locations = []
		writePdbAA.update(ipn.get_peptoid_coordinates())
		locations_temp=[]
		coordinate=[]
		coordinate2 = []
		temp_coordinate=[]
		total_atoms = 1
		peptoids = ipn.get_peptoid_names()
		start_residue = '1'
		amino_acids = ['Ac','A','R','N','D','C','E','N','G','H','I','L','K','M','F','P','S','T','W', 'Y', 'Q','V',
						'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Ser', 'Thr', 'Trp', 'Tyr', 'Gln','Val']
		if start_residue == '':
			print("Empty residue")
			resid_num = 1
			start_residue = 1
		else:
			start_residue = int(start_residue)
			resid_num = int(start_residue)
		#chain_ID = 'A'
		#if chain_ID == '':
		#	chain_ID = 'A'
		#else:
		#	pass
		print(sequence)
		length_strands = int(length)
		width_strands = int(width)
		height_strands = int(height)
		res_num_temp=1
		res_num_temp2=1
		total_chains=1
		total_chains2=1
		for r in range(0, width_strands):
			if r%2 == 0:
				start_pos = 0
			else:
				start_pos = (len(sequence)/2)*3.52+3.80
			for q in range(0, length_strands):
				chain_ID = ' '
				if q == 0:
					adjust = 0
				else:
					adjust = float(l_adjust)*q
				print(adjust)
				for i,k in enumerate(sequence2):
					res_name = writePdbAA.get(k)['Res_name']
			
					atoms = writePdbAA.get(k)['atoms']
			
					atom_type = writePdbAA.get(k)['atom_type']
					sr = start_residue + i
					if i%2 == 0:
						for l,m in enumerate(atoms):
							if k not in peptoids:
								if k == 'Hao':
									temp_coordinate = hao.get(m)
								else:
									temp_coordinate = atom_positions_beta_up.get(m)
							else:
								temp_coordinate = peptoid.get(m)
							for n,p in enumerate(temp_coordinate):
								coordinate.append(p)
							if m == 'CG' and k == 'P':
								coordinate[0] = 1.221
								coordinate[1] = 0.000
								coordinate[2] = 1.921
							elif m == 'CD' and k == 'P':
								coordinate[0] = 0.000
								coordinate[1] = -0.870
								coordinate[2] = 1.521
							else:
								pass
							if (resid_num-1)%((len(sequence2)*length_strands)) == 0:
								res_num_temp = 1
								coordinate[0] = coordinate[0] + res_num_temp*3.675+ adjust+start_pos
							else:
								coordinate[0] = coordinate[0] + res_num_temp*3.675+ adjust+start_pos
							coordinate[2] = coordinate[2] + r*4
							print(res_num_temp)
							if atom_type[l] == 'Br' or atom_type[l] == 'Cl':
								fobj.write('HETATM' + str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'U'+str(total_chains)+ str(atom_type[l]).rjust(4,' ')+ '\n')
							elif k not in amino_acids:
								fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'U'+str(total_chains)+ str(atom_type[l]).rjust(4,' ')+'\n')
							else:
								fobj.write('ATOM' +  str(total_atoms).rjust(7,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'U'+str(total_chains)+ str(atom_type[l]).rjust(4,' ')+'\n')
							del coordinate[:]
							total_atoms+=1
					else:
						for l,m in enumerate(atoms):
							if k not in peptoids:
								if k == 'Hao':
									temp_coordinate = hao.get(m)
								else:
									temp_coordinate = atom_positions_beta_up.get(m)
							else:
								temp_coordinate = peptoid.get(m)
							for n,p in enumerate(temp_coordinate):
								coordinate.append(p)
							print(k)
							if m == 'N':
								coordinate[0] = 0.000
								coordinate[1] = 1.121
							elif m == 'CA':
								coordinate[1] = 0.000
							elif m == 'C':
								coordinate[1] = 1.121
							elif m == 'O':
								coordinate[1] = 2.141
							elif m == 'CG' and k == 'P':
								coordinate[0] = 0.000
								coordinate[1] = -0.870
								coordinate[2] = -1.521
							elif m == 'CD' and k == 'P':
								coordinate[0] = 0.000
								coordinate[1] = 0.870
								coordinate[2] = -1.521
							else:
								coordinate[1] = coordinate[1]*(-1) + 1.170
								coordinate[2] = coordinate[2]*(-1)
							if (resid_num-1)%((len(sequence2)*length_strands)) == 0:
								res_num_temp = 1
								coordinate[0] = coordinate[0] + res_num_temp*3.675+ adjust+start_pos
							else:
								coordinate[0] = coordinate[0] + res_num_temp*3.675+ adjust+start_pos
							coordinate[2] = coordinate[2] + r*4
							if atom_type[l] == 'Br' or atom_type[l] == 'Cl':
								fobj.write('HETATM' + str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'U'+str(total_chains)+ str(atom_type[l]).rjust(4,' ')+ '\n')
							elif k not in amino_acids:
								fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'U'+str(total_chains)+ str(atom_type[l]).rjust(4,' ')+'\n')
							else:
								fobj.write('ATOM' +  str(total_atoms).rjust(7,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'U'+str(total_chains)+ str(atom_type[l]).rjust(4,' ')+'\n')
							del coordinate[:]
							total_atoms+=1
					if k == 'Hao':
						resid_num +=3
					else:
						resid_num +=1
						res_num_temp +=1
				total_atoms2 = 0
				total_chains+=1
				for i,k in enumerate(sequence2):
					res_name = writePdbAA.get(k)['Res_name']
			
					atoms = writePdbAA.get(k)['atoms']
			
					atom_type = writePdbAA.get(k)['atom_type']
					sr = start_residue + i
					if i%2 != 0:
						for l,m in enumerate(atoms):
							if k not in peptoids:
								if k == 'Hao':
									temp_coordinate = hao.get(m)
								else:
									temp_coordinate = atom_positions_beta_up.get(m)
							else:
								temp_coordinate = peptoid.get(m)
							for n,p in enumerate(temp_coordinate):
								coordinate.append(p)
							if m == 'CG' and k == 'P':
								coordinate[0] = 1.221
								coordinate[1] = 0.000
								coordinate[2] = 1.921
							elif m == 'CD' and k == 'P':
								coordinate[0] = 0.000
								coordinate[1] = -0.870
								coordinate[2] = 1.521
							else:
								pass
							if res_num_temp2 ==((len(sequence2)*length_strands)+1):
								res_num_temp2 = 1
								coordinate[0] = coordinate[0] + res_num_temp2*3.675+ adjust+start_pos
							else:
								coordinate[0] = coordinate[0] + res_num_temp2*3.675+ adjust+start_pos
							coordinate[1] = coordinate[1] + height_strands
							coordinate[2] = coordinate[2] + r*4
							if atom_type[l] == 'Br' or atom_type[l] == 'Cl':
								fobj.write('HETATM' + str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'D'+str(total_chains2)+ str(atom_type[l]).rjust(4,' ')+ '\n')
							elif k not in amino_acids:
								fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'D'+str(total_chains2)+ str(atom_type[l]).rjust(4,' ')+ '\n')
							else:
								fobj.write('ATOM' +  str(total_atoms).rjust(7,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'D'+str(total_chains2)+ str(atom_type[l]).rjust(4,' ')+ '\n')
							del coordinate[:]
							total_atoms+=1
					else:
						for l,m in enumerate(atoms):
							if k not in peptoids:
								if k == 'Hao':
									temp_coordinate = hao.get(m)
								else:
									temp_coordinate = atom_positions_beta_up.get(m)
							else:
								temp_coordinate = peptoid.get(m)
							for n,p in enumerate(temp_coordinate):
								coordinate.append(p)
							print(k)
							if m == 'N':
								coordinate[0] = 0.000
								coordinate[1] = 1.121
							elif m == 'CA':
								coordinate[1] = 0.000
							elif m == 'C':
								coordinate[1] = 1.121
							elif m == 'O':
								coordinate[1] = 2.141
							elif m == 'CG' and k == 'P':
								coordinate[0] = 0.000
								coordinate[1] = -0.870
								coordinate[2] = -1.521
							elif m == 'CD' and k == 'P':
								coordinate[0] = 0.000
								coordinate[1] = 0.870
								coordinate[2] = -1.521
							else:
								coordinate[1] = coordinate[1]*(-1) + 1.170
								coordinate[2] = coordinate[2]*(-1)
							if res_num_temp2 ==((len(sequence2)*length_strands)+1):
								res_num_temp2 = 1
								coordinate[0] = coordinate[0] + res_num_temp2*3.675+ adjust+start_pos
							else:
								coordinate[0] = coordinate[0] + res_num_temp2*3.675+ adjust+start_pos
							coordinate[1] = coordinate[1] + height_strands
							coordinate[2] = coordinate[2] + r*4
							if atom_type[l] == 'Br' or atom_type[l] == 'Cl':
								fobj.write('HETATM' + str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'D'+str(total_chains2)+ str(atom_type[l]).rjust(4,' ')+ '\n')
							elif k not in amino_acids:
								fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'D'+str(total_chains2)+ str(atom_type[l]).rjust(4,' ')+ '\n')
							else:
								fobj.write('ATOM' +  str(total_atoms).rjust(7,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'D'+str(total_chains2)+ str(atom_type[l]).rjust(4,' ')+ '\n')
							del coordinate[:]
							total_atoms+=1
						
					if k == 'Hao':
						resid_num +=3
					else:
						resid_num +=1
						res_num_temp2 +=1
				total_atoms2 = 0
				total_chains2+=1
		for q in range(0, length_strands*2*width_strands):
			link = 0
			for start,res in enumerate(sequence2):
				atoms2 = writePdbAA.get(res)['atoms']
				if res not in amino_acids:
					link_rec = writePdbAA.get(res)['link_rec']
					for x,y in enumerate(link_rec):
						if len(y) > 2:
							link1 = y[0] + total_atoms2
							link2 = y[1] + total_atoms2
							link3 = y[2] + total_atoms2
							fobj.write('CONECT' + str(link1).rjust(5,' ') + str(link2).rjust(5, ' ') + str(link3).rjust(5, ' ') +'\n')
						else:
							link1 = y[0] + total_atoms2
							link2 = y[1] + total_atoms2
							fobj.write('CONECT' + str(link1).rjust(5,' ') + str(link2).rjust(5, ' ') + '\n')
					if link < (len(sequence2)-1):
						fobj.write('CONECT' + str(total_atoms2 + 3).rjust(5,' ') + str(total_atoms2 + len(atoms2)+1).rjust(5, ' ') + '\n')
						link +=1
					else:
						pass
					total_atoms2 = total_atoms2 + len(atoms2)
					print(total_atoms2)
				else:
					total_atoms2 = total_atoms2 + len(atoms2)
					link +=1
		fobj.write('END')		
		fobj.close()
		
	def what_the_sheet_anti(self,sequence,length,height,width, l_adjust):
		fname = asksaveasfilename(filetypes=(("pdb", "*.pdb"),
                                           ("All files", "*.*") ),
								  defaultextension='.pdb')
		try:							   
			fobj = open(fname, 'w')
		except:
				pass
		sequence2 = sequence[1:]
			#residue_count +=1
		locations = []
		writePdbAA.update(ipn.get_peptoid_coordinates())
		locations_temp=[]
		coordinate=[]
		coordinate2 = []
		temp_coordinate=[]
		total_atoms = 1
		peptoids = ipn.get_peptoid_names()
		start_residue = '1'
		amino_acids = ['Ac','A','R','N','D','C','E','N','G','H','I','L','K','M','F','P','S','T','W', 'Y', 'Q','V',
						'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Ser', 'Thr', 'Trp', 'Tyr', 'Gln','Val']
		if start_residue == '':
			print("Empty residue")
			resid_num = 1
			start_residue = 1
		else:
			start_residue = int(start_residue)
			resid_num = int(start_residue)
		#chain_ID = 'A'
		#if chain_ID == '':
		#	chain_ID = 'A'
		#else:
		#	pass
		length_strands = int(length)
		width_strands = int(width)
		height_strands = int(height)
		res_num_temp=1
		res_num_temp2=1
		total_chains=1
		total_chains2=1
		for r in range(0, width_strands):
			if r%2 == 0:
				start_pos = 0
			else:
				start_pos = (len(sequence))*3.52*length_strands + (len(sequence)/2)*3.52+3.68
			for q in range(0, length_strands):
				chain_ID = ' '
				if q == 0:
					adjust = 0
				else:
					adjust = float(l_adjust)*q
				for i,k in enumerate(sequence2):
					res_name = writePdbAA.get(k)['Res_name']
			
					atoms = writePdbAA.get(k)['atoms']
			
					atom_type = writePdbAA.get(k)['atom_type']
					sr = start_residue + i
					if i%2 == 0:
						for l,m in enumerate(atoms):
							if k not in peptoids:
								if k == 'Hao':
									temp_coordinate = hao.get(m)
								else:
									temp_coordinate = atom_positions_beta_up.get(m)
							else:
								temp_coordinate = peptoid.get(m)
							for n,p in enumerate(temp_coordinate):
								coordinate.append(p)
							#if r%2 == 0:
							#	pass
							#else:
							#	coordinate[0] = coordinate[0]*-1
							if m == 'CG' and k == 'P':
								coordinate[0] = 1.221
								coordinate[1] = 0.000
								coordinate[2] = 1.921
							elif m == 'CD' and k == 'P':
								coordinate[0] = 0.000
								coordinate[1] = -0.870
								coordinate[2] = 1.521
							else:
								pass
							if (resid_num-1)%((len(sequence2)*length_strands)) == 0:
								res_num_temp = 1
								if r%2 == 0:
									coordinate[0] = coordinate[0] + res_num_temp*3.675+ adjust+start_pos
								else:
									if m =='N' or m =='CA' or m == 'O' or m == 'C' :
										coordinate[0] = coordinate[0]*-1 - res_num_temp*3.675 - adjust+start_pos
									else:
										coordinate[0] = coordinate[0] - res_num_temp*3.675 - adjust+start_pos
							else:
								if r%2 == 0:
									coordinate[0] = coordinate[0] + res_num_temp*3.675+ adjust+start_pos
								else:
									if m =='N' or m =='CA' or m == 'O' or m == 'C' :
										coordinate[0] = coordinate[0]*-1 - res_num_temp*3.675 - adjust+start_pos
									else:
										coordinate[0] = coordinate[0] - res_num_temp*3.675 - adjust+start_pos
							coordinate[2] = coordinate[2] + r*4
							if atom_type[l] == 'Br' or atom_type[l] == 'Cl':
								fobj.write('HETATM' + str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'U'+str(total_chains)+ str(atom_type[l]).rjust(4,' ')+ '\n')
							elif k not in amino_acids:
								fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'U'+str(total_chains)+ str(atom_type[l]).rjust(4,' ')+'\n')
							else:
								fobj.write('ATOM' +  str(total_atoms).rjust(7,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'U'+str(total_chains)+ str(atom_type[l]).rjust(4,' ')+'\n')
							del coordinate[:]
							total_atoms+=1
					else:
						for l,m in enumerate(atoms):
							if k not in peptoids:
								if k == 'Hao':
									temp_coordinate = hao.get(m)
								else:
									temp_coordinate = atom_positions_beta_up.get(m)
							else:
								temp_coordinate = peptoid.get(m)
							for n,p in enumerate(temp_coordinate):
								coordinate.append(p)
							if r%2 == 0:
								if m == 'N':
									coordinate[0] = 0.000
									coordinate[1] = 1.121
								elif m == 'CA':
									coordinate[1] = 0.000
								elif m == 'C':
									coordinate[1] = 1.121
								elif m == 'O':
									coordinate[1] = 2.341
								elif m == 'CG' and k == 'P':
									coordinate[0] = 0.000
									coordinate[1] = -0.870
									coordinate[2] = -1.521
								elif m == 'CD' and k == 'P':
									coordinate[0] = 0.000
									coordinate[1] = 0.870
									coordinate[2] = -1.521
								else:
									coordinate[1] = coordinate[1]*(-1) + 1.170
									coordinate[2] = coordinate[2]*(-1)
							else:
								if m == 'N':
									coordinate[0] = 0.000
									coordinate[1] = 1.121
								elif m == 'CA':
									coordinate[0] = coordinate[0]*-1
									coordinate[1] = 0.000
								elif m == 'C':
									coordinate[0] = coordinate[0]*-1
									coordinate[1] = 1.121
								elif m == 'O':
									coordinate[0] = coordinate[0]*-1
									coordinate[1] = 2.341
								elif m == 'CG' and k == 'P':
									coordinate[0] = 0.000
									coordinate[1] = -0.870
									coordinate[2] = -1.521
								elif m == 'CD' and k == 'P':
									coordinate[0] = 0.000
									coordinate[1] = 0.870
									coordinate[2] = -1.521
								else:
									coordinate[1] = coordinate[1]*(-1) + 1.170
									coordinate[2] = coordinate[2]*(-1)
							if (resid_num-1)%((len(sequence2)*length_strands)) == 0:
								res_num_temp = 1
								if r%2 == 0:
									coordinate[0] = coordinate[0] + res_num_temp*3.675+ adjust+start_pos
								else:
									coordinate[0] = coordinate[0] - res_num_temp*3.675 - adjust+start_pos
							else:
								if r%2 == 0:
									coordinate[0] = coordinate[0] + res_num_temp*3.675+ adjust+start_pos
								else:
									coordinate[0] = coordinate[0] - res_num_temp*3.675 - adjust+start_pos
							coordinate[2] = coordinate[2] + r*4
							if atom_type[l] == 'Br' or atom_type[l] == 'Cl':
								fobj.write('HETATM' + str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'U'+str(total_chains)+ str(atom_type[l]).rjust(4,' ')+ '\n')
							elif k not in amino_acids:
								fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'U'+str(total_chains)+ str(atom_type[l]).rjust(4,' ')+'\n')
							else:
								fobj.write('ATOM' +  str(total_atoms).rjust(7,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'U'+str(total_chains)+ str(atom_type[l]).rjust(4,' ')+'\n')
							del coordinate[:]
							total_atoms+=1
					if k == 'Hao':
						resid_num +=3
					else:
						resid_num +=1
						res_num_temp +=1
				total_atoms2 = 0
				total_chains+=1
				for i,k in enumerate(sequence2):
					res_name = writePdbAA.get(k)['Res_name']
			
					atoms = writePdbAA.get(k)['atoms']
			
					atom_type = writePdbAA.get(k)['atom_type']
					sr = start_residue + i
					if i%2 != 0:
						for l,m in enumerate(atoms):
							if k not in peptoids:
								if k == 'Hao':
									temp_coordinate = hao.get(m)
								else:
									temp_coordinate = atom_positions_beta_up.get(m)
							else:
								temp_coordinate = peptoid.get(m)
							for n,p in enumerate(temp_coordinate):
								coordinate.append(p)
							#if r%2 == 0:
							#	pass
							#else:
							#	coordinate[0]==coordinate[0]*-1
							if m == 'CG' and k == 'P':
								coordinate[0] = 1.221
								coordinate[1] = 0.000
								coordinate[2] = 1.921
							elif m == 'CD' and k == 'P':
								coordinate[0] = 0.000
								coordinate[1] = -0.870
								coordinate[2] = 1.521
							else:
								pass
							if res_num_temp2 ==((len(sequence2)*length_strands)+1):
								res_num_temp2 = 1
								if r%2 == 0:
									coordinate[0] = coordinate[0] + res_num_temp2*3.675+ adjust+start_pos
								else:
									if m =='N' or m =='CA' or m == 'O' or m == 'C' :
										coordinate[0] = coordinate[0]*-1 - res_num_temp2*3.675 - adjust+start_pos
									else:
										coordinate[0] = coordinate[0] - res_num_temp2*3.675 - adjust+start_pos
							else:
								if r%2 == 0:
									coordinate[0] = coordinate[0] + res_num_temp2*3.675+ adjust+start_pos
								else:
									if m =='N' or m =='CA' or m == 'O' or m == 'C' :
										coordinate[0] = coordinate[0]*-1 - res_num_temp2*3.675 - adjust+start_pos
									else:
										coordinate[0] = coordinate[0] - res_num_temp2*3.675 - adjust+start_pos
							coordinate[1] = coordinate[1] + height_strands
							coordinate[2] = coordinate[2] + r*4
							if atom_type[l] == 'Br' or atom_type[l] == 'Cl':
								fobj.write('HETATM' + str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'D'+str(total_chains2)+ str(atom_type[l]).rjust(4,' ')+ '\n')
							elif k not in amino_acids:
								fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'D'+str(total_chains2)+ str(atom_type[l]).rjust(4,' ')+ '\n')
							else:
								fobj.write('ATOM' +  str(total_atoms).rjust(7,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'D'+str(total_chains2)+ str(atom_type[l]).rjust(4,' ')+ '\n')
							del coordinate[:]
							total_atoms+=1
					else:
						for l,m in enumerate(atoms):
							if k not in peptoids:
								if k == 'Hao':
									temp_coordinate = hao.get(m)
								else:
									temp_coordinate = atom_positions_beta_up.get(m)
							else:
								temp_coordinate = peptoid.get(m)
							for n,p in enumerate(temp_coordinate):
								coordinate.append(p)
							if r%2 == 0:
								if m == 'N':
									coordinate[0] = 0.000
									coordinate[1] = 1.121
								elif m == 'CA':
									coordinate[1] = 0.000
								elif m == 'C':
									coordinate[1] = 1.121
								elif m == 'O':
									coordinate[1] = 2.341
								elif m == 'CG' and k == 'P':
									coordinate[0] = 0.000
									coordinate[1] = -0.870
									coordinate[2] = -1.521
								elif m == 'CD' and k == 'P':
									coordinate[0] = 0.000
									coordinate[1] = 0.870
									coordinate[2] = -1.521
								else:
									coordinate[1] = coordinate[1]*(-1) + 1.170
									coordinate[2] = coordinate[2]*(-1)
							else:
								
								if m == 'N':
									coordinate[0] = 0.000
									coordinate[1] = 1.121
								elif m == 'CA':
									coordinate[0]=coordinate[0]*-1
									coordinate[1] = 0.000
								elif m == 'C':
									coordinate[0]=coordinate[0]*-1
									coordinate[1] = 1.121
								elif m == 'O':
									coordinate[0]=coordinate[0]*-1
									coordinate[1] = 2.341
								elif m == 'CG' and k == 'P':
									coordinate[0] = 0.000
									coordinate[1] = -0.870
									coordinate[2] = -1.521
								elif m == 'CD' and k == 'P':
									coordinate[0] = 0.000
									coordinate[1] = 0.870
									coordinate[2] = -1.521
								else:
									coordinate[1] = coordinate[1]*(-1) + 1.170
									coordinate[2] = coordinate[2]*(-1)
							if res_num_temp2 ==((len(sequence2)*length_strands)+1):
								res_num_temp2 = 1
								if r%2 == 0:
									coordinate[0] = coordinate[0] + res_num_temp2*3.675+ adjust+start_pos
								else:
									coordinate[0] = coordinate[0] - res_num_temp2*3.675 - adjust+start_pos
							else:
								if r%2 == 0:
									coordinate[0] = coordinate[0] + res_num_temp2*3.675+ adjust+start_pos
								else:
									coordinate[0] = coordinate[0] - res_num_temp2*3.675 - adjust+start_pos
							coordinate[1] = coordinate[1] + height_strands
							coordinate[2] = coordinate[2] + r*4
							if atom_type[l] == 'Br' or atom_type[l] == 'Cl':
								fobj.write('HETATM' + str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'D'+str(total_chains2)+ str(atom_type[l]).rjust(4,' ')+ '\n')
							elif k not in amino_acids:
								fobj.write('HETATM' +  str(total_atoms).rjust(5,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'D'+str(total_chains2)+ str(atom_type[l]).rjust(4,' ')+ '\n')
							else:
								fobj.write('ATOM' +  str(total_atoms).rjust(7,' ')+ '  ' + str(atoms[l]).ljust(4,' ') + str(res_name).strip() + str(chain_ID).rjust(2,' ') + str(sr).rjust(4,' ')+ str(format(coordinate[0],'.3f')).rjust(12,' ') + str(format(coordinate[1],'.3f')).rjust(8,' ') + str(format(coordinate[2],'.3f')).rjust(8,' ')+ '  1.00 30.00' + ' '*6 + 'D'+str(total_chains2)+ str(atom_type[l]).rjust(4,' ')+ '\n')
							del coordinate[:]
							total_atoms+=1
						
					if k == 'Hao':
						resid_num +=3
					else:
						resid_num +=1
						res_num_temp2 +=1
				total_atoms2 = 0
				total_chains2+=1
		for q in range(0, length_strands*2*width_strands):
			link = 0
			for start,res in enumerate(sequence2):
				atoms2 = writePdbAA.get(res)['atoms']
				if res not in amino_acids:
					link_rec = writePdbAA.get(res)['link_rec']
					for x,y in enumerate(link_rec):
						if len(y) > 2:
							link1 = y[0] + total_atoms2
							link2 = y[1] + total_atoms2
							link3 = y[2] + total_atoms2
							fobj.write('CONECT' + str(link1).rjust(5,' ') + str(link2).rjust(5, ' ') + str(link3).rjust(5, ' ') +'\n')
						else:
							link1 = y[0] + total_atoms2
							link2 = y[1] + total_atoms2
							fobj.write('CONECT' + str(link1).rjust(5,' ') + str(link2).rjust(5, ' ') + '\n')
					if link < (len(sequence2)-1):
						fobj.write('CONECT' + str(total_atoms2 + 3).rjust(5,' ') + str(total_atoms2 + len(atoms2)+1).rjust(5, ' ') + '\n')
						link +=1
					else:
						pass
					total_atoms2 = total_atoms2 + len(atoms2)
				else:
					total_atoms2 = total_atoms2 + len(atoms2)
					link +=1
		fobj.write('END')		
		fobj.close()
		print("Finished writing file")
			
		