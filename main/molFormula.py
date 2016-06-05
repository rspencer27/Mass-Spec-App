
class Pretty_Mol_Formula(object):
	
	def __init__(self, mol_formula):
		self.mol_formula = mol_formula
	
	def make_pretty(self,mol_formula):
		mol_form_temp = []
		for i,k in enumerate(mol_formula):
			if k != 0:
				if i == 0:
					mol_form_temp.append('C')
					mol_form_temp.append(k)
				elif i == 1:
					mol_form_temp.append('H')
					mol_form_temp.append(k)
				elif i == 2:
					mol_form_temp.append('N')
					mol_form_temp.append(k)
				elif i == 3:
					mol_form_temp.append('O')
					mol_form_temp.append(k)
				elif i == 4:
					mol_form_temp.append('S')
					mol_form_temp.append(k)
				elif i == 5:
					mol_form_temp.append('F')
					mol_form_temp.append(k)
				elif i == 6:
					mol_form_temp.append('Cl')
					mol_form_temp.append(k)
				elif i == 7:
					mol_form_temp.append('Br')
					mol_form_temp.append(k)
				elif i == 8:
					mol_form_temp.append('I')
					mol_form_temp.append(k)
			else:
				pass
		return mol_form_temp
			