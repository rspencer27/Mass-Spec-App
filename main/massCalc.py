import re
from main.fileImport import *
from copy import deepcopy

ipn = Import_Peptoid_Names()

class Sequence_Convert(object):
	'''Main class for converting input string to list'''
    
    
	def __init__(self, sequence):
		self.sequence = sequence
        

	def sequence_to_list(self,sequence):
		stripped_sequence = str(sequence).replace(' ','').replace('(NMe)','nme').replace('(Br)','bromo').replace('(I)','iodo').replace('(5F)','5fluoro').replace('TFAsalt','Tfasalt').replace('TFAester','Tfaester').strip()
		sequence = re.findall('[A-Z][^A-Z]*', stripped_sequence)
		return sequence
    
	def triple_to_single(self,sequence):
		sequence_converted=[]
		temp = self.sequence_to_list(sequence)
		for i, n in enumerate(temp):
			if n in single_letter_AA:
				sequence_converted.append(single_letter_AA.get(n))
			else:
				sequence_converted.append(n)
		return sequence_converted

class Mass_Spec_Calc(object):
	'''This calculates the mass spec based on the converted sequence'''
    
	def __init__(self,sequence):
		self.sequence = sequence
    
	def mol_formula(self,sequence):
		carbons=0;hydrogens=0;nitrogens=0;oxygens=0;sulfurs=0;fluorines=0;chlorines=0;bromines=0;iodines=0
		for i, n in enumerate(sequence):
			carbons = carbons + exactMassAA.get(n)[0]
			hydrogens = hydrogens + exactMassAA.get(n)[1]
			nitrogens = nitrogens + exactMassAA.get(n)[2]
			oxygens = oxygens + exactMassAA.get(n)[3]
			sulfurs = sulfurs + exactMassAA.get(n)[4]
			fluorines = fluorines + exactMassAA.get(n)[5]
			chlorines = chlorines + exactMassAA.get(n)[6]
			bromines = bromines + exactMassAA.get(n)[7]
			iodines = iodines + exactMassAA.get(n)[8]
		return carbons, hydrogens, nitrogens, oxygens, sulfurs, fluorines, chlorines,bromines, iodines
    
	def exact_mass_calc(self, mol_form):
		exact_mass = mol_form[0]*12.00000 + mol_form[1]*1.007825 + mol_form[2]*14.003074 + mol_form[3]*15.994915 + mol_form[4]*31.972072 + mol_form[5]*18.998403 +mol_form[6]*34.968852+ mol_form[7]*78.918336 + mol_form[8]*126.904477
		return exact_mass
	
	def molecular_weight(self, mol_form):
		mol_weight = exact_mass = mol_form[0]*12.011 + mol_form[1]*1.0079 + mol_form[2]*14.0067 + mol_form[3]*15.9994 + mol_form[4]*32.065 + mol_form[5]*18.998403 + mol_form[6]*35.453 + mol_form[7]*79.904 + mol_form[8]*126.904477
		return mol_weight
		
	def mass_spec_peaks(self, linear_mass, sequence):
		'''This generates a dictionary based on the possible ion combination and returns m/z masses'''
		H = 1.007825; Na = 22.989770; K = 38.963708
		i = 0; j = 0; k = 0; l = 1
		masses = {}
		mass_calc = 0
		TFA_ester = sequence.count('S') + sequence.count('T') + sequence.count('Ser') + sequence.count('Thr')
		for l in range(1,6):
			for k in range (0,4):
				for j in range (0,4):
					for i in range (0,6):
						if i==0 and j==0 and k==0:
							mass_label = "Blank"
							masses[0] = mass_label								
						else:
							for p in range(0, TFA_ester+1):
								mass_temp = ((linear_mass)*l+ 95.99012*p + i*H +j*Na +k*K)/(i+j+k)
								mass_calc = float(format(mass_temp,'.5f'))
								if mass_calc in masses:
									pass
								else:
									if p == 0:
										if i == 0 and j == 0 and k!=0:
											mass_label = "[%sM + %sK]+%s" % (l,k,k)
											masses[mass_calc] = mass_label
										elif i == 0 and k == 0 and j !=0:
											mass_label = "[%sM + %sNa]+%s"%(l,j,j)
											masses[mass_calc] = mass_label
										elif i != 0 and j == 0 and k == 0:
											mass_label = "[%sM + %sH]+%s" %(l,i,i)
											masses[mass_calc] = mass_label
										elif i != 0 and j !=0 and k == 0:
											mass_label = "[%sM + %sH + %sNa]+%s" %(l,i,j, i+j)
											masses[mass_calc] = mass_label
										elif i !=0 and j ==0 and k != 0:
											mass_label = "[%sM + %sH + %sK]+%s" %(l,i,k, i+k)
											masses[mass_calc] = mass_label
										elif i ==0 and j !=0 and k != 0:
											mass_label = "[%sM + %sNa + %sK]+%s" %(l,j,k, j+k)
											masses[mass_calc] = mass_label
										else:
											mass_label = "[%sM + %sH +%sNa + %sK]+%s" %(l,i,j,k, i+j+k)
											masses[mass_calc] = mass_label
									else:
										if i == 0 and j == 0 and k!=0:
											mass_label = "[%sM + %sK + %sTFA ester]+%s" % (l,k,p,k)
											masses[mass_calc] = mass_label
										elif i == 0 and k == 0 and j !=0:
											mass_label = "[%sM + %sNa + %sTFA ester]+%s"%(l,j,p,j)
											masses[mass_calc] = mass_label
										elif i != 0 and j == 0 and k == 0:
											mass_label = "[%sM + %sH + %sTFA ester]+%s" %(l,i,p,i)
											masses[mass_calc] = mass_label
										elif i != 0 and j !=0 and k == 0:
											mass_label = "[%sM + %sH + %sNa + %sTFA ester]+%s" %(l,i,j,p,i+j)
											masses[mass_calc] = mass_label
										elif i !=0 and j ==0 and k != 0:
											mass_label = "[%sM + %sH + %sK + %sTFA ester]+%s" %(l,i,k,p,i+k)
											masses[mass_calc] = mass_label
										elif i ==0 and j !=0 and k != 0:
											mass_label = "[%sM + %sNa + %sK + %sTFA ester]+%s" %(l,j,k,p,j+k)
											masses[mass_calc] = mass_label
										else:
											mass_label = "[%sM + %sH +%sNa + %sK + %sTFA ester]+%s" %(l,i,j,k,p,i+j+k)
											masses[mass_calc] = mass_label
       
		return masses
    
	def mass_spec_compare(self,peak,masses):
		answer = None
		for mass, ion in masses.items():
			if abs((peak) - (mass)) < 0.30:
				answer = ion + '   :   ' + str(format(mass, '.5f'))
				break
			else:
				pass
		if answer == None:
			for mass, ion in masses.items():
				if abs((peak) - (mass)) < 0.50:
					answer = ion + '   :   ' + str(format(mass, '.5f'))
					break
				else:
					pass
		if answer == None:
			answer = 'No Match'
			return answer
		else:
			return answer
	
	def mass_spec_compare_del_add(self,peak,masses):
		answer = None
		answer_del_add=[]
		for mass, ion in masses.items():
			if abs((peak) - (mass)) < 0.50:
				answer = ion + '   :   ' + str(format(mass, '.5f'))
				answer_del_add.append(answer)		
			else:
				pass
		if answer == None:
			answer = []
			answer.append('No Match')
			return answer
		else:
			return answer_del_add
                        
	def mass_spec_del_add(self, linear_mass, sequence):
		H = 1.007825; Na = 22.989770; K = 38.963708
		i = 0; j = 0; k = 0; l = 1
		masses = {}
		masses2={}
		mass_calc = 0
		print(sequence)
		for l in range(1,2):
			for k in range(0,3):
				for j in range(0,3):
					for i in range(0,4):
						if i==0 and j==0 and k==0:
							mass_label = "Blank"
							masses[0] = mass_label
						else:
							for m in range(0,len(sequence)-1):
								z = Mass_Spec_Calc.exact_mass_calc('self',exactMassAA.get(sequence[m]))
								dd_truncate = sequence[:]
								dd_truncate.remove(sequence[m])
								temp_amino_acid = sequence[m]
								mass_del = ((float(linear_mass)-float(z))*l+i*H +j*Na +k*K)/(i+j+k)
								mass_add = ((float(linear_mass)+float(z))*l+i*H +j*Na +k*K)/(i+j+k)
								if i == 0 and j == 0 and k!=0:
									mass_label_del = "[%sM + %sK -%s]+%s" %(l,k,sequence[m],k)
									mass_label_add = "[%sM + %sK +%s]+%s" %(l,k,sequence[m],k)
									masses[mass_del] = mass_label_del
									masses[mass_add] = mass_label_add
								elif i == 0 and k == 0 and j !=0:
									mass_label_del = "[%sM + %sNa -%s]+%s" %(l,j,sequence[m],j)
									mass_label_add = "[%sM + %sNa +%s]+%s" %(l,j,sequence[m],j)
									masses[mass_del] = mass_label_del
									masses[mass_add] = mass_label_add
								elif i != 0 and j == 0 and k == 0:
									mass_label_del = "[%sM + %sH -%s]+%s" %(l,i,sequence[m],i)
									mass_label_add = "[%sM + %sH +%s]+%s" %(l,i,sequence[m],i)
									masses[mass_del] = mass_label_del
									masses[mass_add] = mass_label_add
								elif i != 0 and j !=0 and k == 0:
									mass_label_del = "[%sM + %sH + %sNa -%s]+%s" %(l,i,j,sequence[m],i+j)
									mass_label_add = "[%sM + %sH + %sNa +%s]+%s" %(l,i,j,sequence[m],i+j)
									masses[mass_del] = mass_label_del
									masses[mass_add] = mass_label_add
								elif i !=0 and j ==0 and k != 0:
									mass_label_del = "[%sM + %sH + %sK -%s]+%s" %(l,i,k,sequence[m],i+k)
									mass_label_add = "[%sM + %sH + %sK +%s]+%s" %(l,i,k,sequence[m],i+k)
									masses[mass_del] = mass_label_del
									masses[mass_add] = mass_label_add
								elif i ==0 and j !=0 and k != 0:
									mass_label_del = "[%sM + %sNa + %sK -%s]+%s" %(l,j,k,sequence[m],j+k)
									mass_label_add = "[%sM + %sNa + %sK +%s]+%s" %(l,j,k,sequence[m],j+k)
									masses[mass_del] = mass_label_del
									masses[mass_add] = mass_label_add
								else:
									mass_label_del = "[%sM + %sH + %sNa + %sK -%s]+%s" %(l,i,j,k,sequence[m],i+j+k)
									mass_label_add = "[%sM + %sH + %sNa + %sK +%s]+%s" %(l,i,j,k,sequence[m],i+j+k)
									masses[mass_del] = mass_label_del
									masses[mass_add] = mass_label_add	
								for y in range (0, len(dd_truncate)-1):
									x = Mass_Spec_Calc.exact_mass_calc('self',exactMassAA.get(dd_truncate[y]))
									mass_ddel = ((float(linear_mass)-float(z)-float(x))*l+i*H +j*Na +k*K)/(i+j+k)
									mass_dadd = ((float(linear_mass)+float(z)+float(x))*l+i*H +j*Na +k*K)/(i+j+k)
									mass_del_add = ((float(linear_mass)-float(z)+float(x))*l+i*H +j*Na +k*K)/(i+j+k)
									temp_amino_acid2 = dd_truncate[y]
									mass_del_add_temp = float(format(mass_del_add,'.5f'))
									if i == 0 and j == 0 and k!=0:
										mass_label_ddel = "[%sM + %sK -%s -%s]+%s" %(l,k,temp_amino_acid,temp_amino_acid2,k)
										mass_label_dadd = "[%sM + %sK +%s +%s]+%s" %(l,k,temp_amino_acid,temp_amino_acid2,k)
										masses2[mass_ddel] = mass_label_ddel
										masses2[mass_dadd] = mass_label_dadd
									elif i == 0 and k == 0 and j !=0:
										mass_label_ddel = "[%sM + %sNa -%s -%s]+%s" %(l,j,temp_amino_acid,temp_amino_acid2,j)
										mass_label_dadd = "[%sM + %sNa +%s +%s]+%s" %(l,j,temp_amino_acid,temp_amino_acid2,j)
										masses2[mass_ddel] = mass_label_ddel
										masses2[mass_dadd] = mass_label_dadd
									elif i != 0 and j == 0 and k == 0:
										mass_label_ddel = "[%sM + %sH -%s -%s]+%s" %(l,i,temp_amino_acid,temp_amino_acid2,i)
										mass_label_dadd = "[%sM + %sH +%s +%s]+%s" %(l,i,temp_amino_acid,temp_amino_acid2,i)
										masses2[mass_ddel] = mass_label_ddel
										masses2[mass_dadd] = mass_label_dadd
									elif i != 0 and j !=0 and k == 0:
										mass_label_ddel = "[%sM + %sH + %sNa -%s -%s]+%s" %(l,i,j,temp_amino_acid,temp_amino_acid2,i+j)
										mass_label_dadd = "[%sM + %sH + %sNa +%s +%s]+%s" %(l,i,j,temp_amino_acid,temp_amino_acid2,i+j)
										masses2[mass_ddel] = mass_label_ddel
										masses2[mass_dadd] = mass_label_dadd
									elif i !=0 and j ==0 and k != 0:
										mass_label_ddel = "[%sM + %sH + %sK -%s -%s]+%s" %(l,i,k,temp_amino_acid,temp_amino_acid2,i+k)
										mass_label_dadd = "[%sM + %sH + %sK +%s +%s]+%s" %(l,i,k,temp_amino_acid,temp_amino_acid2,i+k)
										masses2[mass_ddel] = mass_label_ddel
										masses2[mass_dadd] = mass_label_dadd
									elif i ==0 and j !=0 and k != 0:
										mass_label_ddel = "[%sM + %sNa + %sK -%s -%s]+%s" %(l,j,k,temp_amino_acid,temp_amino_acid2,j+k)
										mass_label_dadd = "[%sM + %sNa + %sK +%s +%s]+%s" %(l,j,k,temp_amino_acid,temp_amino_acid2,j+k)
										masses2[mass_ddel] = mass_label_ddel
										masses2[mass_dadd] = mass_label_dadd
									else:
										mass_label_ddel = "[%sM + %sH + %sNa + %sK -%s -%s]+%s" %(l,i,j,k,temp_amino_acid,temp_amino_acid2,i+j+k)
										mass_label_dadd = "[%sM + %sH + %sNa + %sK +%s +%s]+%s" %(l,i,j,k,temp_amino_acid,temp_amino_acid2,i+j+k)
										masses2[mass_ddel] = mass_label_ddel
										masses2[mass_dadd] = mass_label_dadd
									if mass_del_add_temp in masses:
										pass
									elif z == x :
										pass
									else:
										if i == 0 and j == 0 and k!=0:
											mass_label_del_add = "[%sM + %sK -%s +%s]+%s" %(l,k,temp_amino_acid,temp_amino_acid2,k)
											masses2[mass_del_add] = mass_label_del_add
										elif i == 0 and k == 0 and j !=0:
											mass_label_del_add = "[%sM + %sNa -%s +%s]+%s" %(l,j,temp_amino_acid,temp_amino_acid2,j)
											masses2[mass_del_add] = mass_label_del_add
										elif i != 0 and j == 0 and k == 0:
											mass_label_del_add = "[%sM + %sH -%s +%s]+%s" %(l,i,temp_amino_acid,temp_amino_acid2,i)
											masses2[mass_del_add] = mass_label_del_add
										elif i != 0 and j !=0 and k == 0:
											mass_label_del_add = "[%sM + %sH + %sNa -%s +%s]+%s" %(l,i,j,temp_amino_acid,temp_amino_acid2,i+j)
											masses2[mass_del_add] = mass_label_del_add
										elif i !=0 and j ==0 and k != 0:
											mass_label_del_add = "[%sM + %sH + %sK -%s +%s]+%s" %(l,i,k,temp_amino_acid,temp_amino_acid2,i+k)
											masses2[mass_del_add] = mass_label_del_add
										elif i ==0 and j !=0 and k != 0:
											mass_label_del_add = "[%sM + %sNa + %sK -%s +%s]+%s" %(l,j,k,temp_amino_acid,temp_amino_acid2,j+k)
											masses2[mass_del_add] = mass_label_del_add
										else:
											mass_label_del_add = "[%sM + %sH + %sNa + %sK -%s +%s]+%s" %(l,i,j,k,temp_amino_acid,temp_amino_acid2,i+j+k)
											masses2[mass_del_add] = mass_label_del_add
		return masses, masses2 
       

        

single_letter_AA = {'Ala' : 'A', 'Arg' : 'R', 'Asn' : 'N','Asp' : 'D','Cys' : 'C','Glu' : 'E','Gln' : 'Q','Gly' : 'G',
                    'His': 'H','Ile': 'I','Leu' : 'L','Lys' : 'K','Met' : 'M','Phe' : 'F','Pro' : 'P' ,'Ser': 'S','Thr' : 'T',
                    'Trp' : 'W', 'Tyr':'Y', 'Val' : 'V','Acid' : 'OH', 'Amide' : 'NH2', 'C-Term' : 'C-Term','Cyclic' : 'Cyclic',
                    'Pheiodo' : 'Phe(I)', 'Phebromo' : 'Phe(Br)', 'Phe5fluoro' : 'Phe(5F)', 'Orn' : 'O','Hao' : 'Hao',
                    'Alanme' : 'Ala(NMe)', 'Glynme' : 'Gly(NMe)','Ilenme' : 'Ile(NMe)', 'Leunme' : 'Leu(NMe)',
                    'Phenme' :  'Phe(NMe)', 'Valnme' : 'Val(NMe)','Ac' : 'Ac', 'Tfasalt' : 'TFA salt', 'Tfaester' : 'TFA ester',
                    'Nle' : 'Nle', 'Nva' : 'Nva', 'Chg' : 'Chg','Nlenme' : 'Nle(NMe)', 'Tyrnme' : 'Tyr(NMe)',
					'Nvanme' : 'Nva(NMe)', 'Fmoc' : 'Fmoc', 'Disulfide' : 'Disulfide', 'Cha':'Cha', 'Tol':'Tol'}

# List for residues {Name, [carbons, hydrogens, nitrogens, oxygens, sulfurs, fluorines, chlorines, bromines, iodines]}
exactMassAA = {'A' : [3,5,1,1,0,0,0,0,0], 'R' : [6,12,4,1,0,0,0,0,0], 'N' : [4,6,2,2,0,0,0,0,0], 'D' : [4,5,1,3,0,0,0,0,0],'C' : [3,5,1,1,1,0,0,0,0],
               'E' : [5,7,1,3,0,0,0,0,0], 'Q' : [5,8,2,2,0,0,0,0,0],   'G' : [2,3,1,1,0,0,0,0,0], 'H' : [6,7,3,1,0,0,0,0,0],'I' : [6,11,1,1,0,0,0,0,0],
               'L' : [6,11,1,1,0,0,0,0,0],'K' : [6,12,2,1,0,0,0,0,0],  'M' : [5,9,1,1,1,0,0,0,0], 'F' : [9,9,1,1,0,0,0,0,0],'P' : [5,7,1,1,0,0,0,0,0],
               'S' : [3,5,1,2,0,0,0,0,0], 'T' : [4,7,1,2,0,0,0,0,0],   'W' : [11,10,2,1,0,0,0,0,0], 'Y': [9,9,1,2,0,0,0,0,0], 'V' : [5,9,1,1,0,0,0,0,0],
               'OH' : [0,2,0,1,0,0,0,0,0], 'NH2' : [0,3,1,0,0,0,0,0,0], 'C-Term' : [0,3,1,1,0,0,0,0,0],'Cyclic' : [0,0,0,0,0,0,0,0,0], 'Phe(I)' : [9,8,1,1,0,0,0,0,1],
			   'Phe(Br)' : [9,8,1,1,0,0,0,1,0], 'Phe(5F)' : [9,4,1,1,0,5,0,0,0],'O' :[5,10,2,1,0,0,0,0,0] ,'Hao' : [10,9,3,4,0,0,0,0,0], 'Ala(NMe)' : [4,7,1,1,0,0,0,0,0],
			   'Gly(NMe)' : [3,5,1,1,0,0,0,0,0], 'Ile(NMe)' : [7,13,1,1,0,0,0,0,0], 'Leu(NMe)' : [7,13,1,1,0,0,0,0,0], 'Phe(NMe)' : [10,11,1,1,0,0,0,0,0],
			   'Val(NMe)' : [6,11,1,1,0,0,0,0,0],'Ac' : [2,2,0,1,0,0,0,0,0] , 'TFA salt' : [2,1,0,2,0,3,0,0,0], 'TFA ester' : [2,0,0,1,0,3,0,0,0],
			   'Nle' : [5,9,1,1,0,0,0,0,0], 'Nva' : [5,9,1,1,0,0,0,0,0], 'Chg' : [8,13,1,1,0,0,0,0,0],'Nle(NMe)' : [7,13,1,1,0,0,0,0,0], 'Tyr(NMe)' : [10,11,1,1,0,0,0,0,0],
			   'Nva(NMe)' : [6,11,1,1,0,0,0,0,0], 'Fmoc' : [15,10,0,2,0,0,0,0,0], 'Disulfide' : [0,-2,0,0,0,0,0,0,0], 'Pra': [5,5,1,1,0,0,0,0,0], 'Cha':[9,15,1,1,0,0,0,0,0], 'Tol':[10,11,1,1,0,0,0,0,0]}



					
exactMassAA.update(ipn.get_peptoid_mass())
		   



