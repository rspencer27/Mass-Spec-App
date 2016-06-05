from tkinter.filedialog import asksaveasfilename
from main.massCalc import *
from main.molFormula import *
import textwrap

class printReport(object):
	''' Main class for printing a report '''
	def __init__(self):
		pass
		
	def write_Report (self, sequence_print, peaks):
		convert = Sequence_Convert('sequence')
		mass_spec_calc = Mass_Spec_Calc('mass_calc')
		pretty_formula = Pretty_Mol_Formula('mol_formula')
		fname = asksaveasfilename(filetypes=(("txt", "*.txt"),
                                           ("All files", "*.*") ),
								  defaultextension='.txt')
		try:							   
			fobj = open(fname, 'w')
		except:
				pass
		
		converted_sequence = convert.triple_to_single(''.join(str(t) for t in sequence_print))
		seq = converted_sequence[1:]
		seq.append(converted_sequence[0])
		print(converted_sequence)
		mol_formula = list(mass_spec_calc.mol_formula(converted_sequence))

		pretty_mol_formula = pretty_formula.make_pretty(mol_formula)
		
		mol_formula_pretty = ' '.join(str(e) for e in pretty_mol_formula)
		exact_mass2 = mass_spec_calc.exact_mass_calc(mol_formula)
		exact_mass_formatted = str(format(exact_mass2, '.5f'))
		mol_weight = mass_spec_calc.molecular_weight(mol_formula)
		mol_weight_formatted = str(format(mol_weight,'.5f'))
		count = -1
		for i in sequence_print:
			count +=1
		
		sequence_count = str(count)
		
		pretty = textwrap.wrap('-'.join(str(t) for t in seq),width=60)
		polymer = seq[:-1]
		p = dict((i,polymer.count(i)) for i in polymer)
		polymer_nice = str(p)[1:-1]
		
		fobj.write('-'*80+'\n')
		fobj.write(' '*29+ 'Sequence and MS Report' + ' '*29)
		fobj.write('-'*80+'\n'*3)
		fobj.write(' '*5 + 'Sequence: ')
		for i, text in enumerate(pretty):
			if i == 0:
				fobj.write(text + '\n')
			else:
				fobj.write(' '*15 + text + '\n')
		fobj.write(' '*5 + 'Residues: ' + sequence_count + '\n')
		fobj.write('  Mol Formula: ' + mol_formula_pretty + '\n')
		fobj.write('   Exact Mass: ' + exact_mass_formatted + '\n')
		fobj.write('   Mol Weight: ' + mol_weight_formatted + '\n')
		fobj.write('      Repeats: ' + polymer_nice + '\n')
		print('Peaks')
		print(peaks)
		if not peaks:
			print('No peaks checked')
		else:
			fobj.write('-'*80+'\n')
			fobj.write(' '*35+ 'Peak Report' + ' '*34)
			fobj.write('-'*80+'\n'*2)
			
			for number, peak in enumerate(peaks):
				fobj.write(' '*5 + 'Peak Input: ' + str(peak) + '\n')
				peak = float(peak)
				check_value = mass_spec_calc.mass_spec_compare(peak,mass_spec_calc.mass_spec_peaks(exact_mass2,sequence_print))
				check_value_temp = mass_spec_calc.mass_spec_del_add(exact_mass2,converted_sequence)
				check_val_add_del=[]
				check_val_add_del.append(mass_spec_calc.mass_spec_compare_del_add(peak,check_value_temp[0]))
				check_val_dd=[]
				check_val_dd.append(mass_spec_calc.mass_spec_compare_del_add(peak,check_value_temp[1]))
				
				add_del = ', '.join(str(e) for e in check_val_add_del[0])
				dd = ', '.join(str(e) for e in check_val_dd[0])
				fobj.write(' '*5 +'Peak Match: ' + check_value + '\n')
				fobj.write(' '*5 +'Add or Del: ' + add_del + '\n')
				fobj.write(' '*2 +'Other Matches: ' + dd + '\n'*2)
				fobj.write('-'*80+'\n'*2)
		fobj.close()
		