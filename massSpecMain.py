import sys
import os.path
from main.massCalc import *
from main.molFormula import *
from main.fileImport import Import_Peptoid_Names
from main.writePdb import *
from main.printReport import *
from main.remove_neg_log import *
import Pmw as Pmw
import tkinter as tk
from tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename


#written by Ryan K. Spencer

resin_choices = ['Chlorotrityl', 'Rink', 'Wang', 'Cyclic','Expressed']

sequence = []

final_sequence=['Acid']
peak = 0
exact_mass = 0
mol_formula=[]
convert = Sequence_Convert('sequence')
mass_spec_calc = Mass_Spec_Calc('mass_calc')
pretty_formula = Pretty_Mol_Formula('mol_formula')
ipn = Import_Peptoid_Names()
wPdb = write_Pdb('self')
report = printReport()
global user_peaks
user_peaks = []
global converted_sequence


class MassSpecApp(object):
	
	def __init__(self, master, **kwargs):
		
		label = tk.Label(root,
								anchor="w", text="Choose an Option")
		label.grid(row=0, column=0)
		label = tk.Label(root,
								anchor="e", text="Sequence: ")
		label.grid(row=2, column=1, padx=5,sticky=E)
		
		label = tk.Label(root,
								anchor="e", text="Molecular Formula: ")
		label.grid(row=5, column=1, padx=5,sticky=E)
		
		label = tk.Label(root,
								anchor="e", text="Exact Mass: ")
		label.grid(row=6, column=1, padx=5,sticky=E)
		
		label = tk.Label(root,
								anchor="e", text="Molecular Weight: ")
		label.grid(row=7, column=1, padx=5,sticky=E)
		
		label = tk.Label(root,
								anchor="e", text="Mass Spec Match: ")
		label.grid(row=8, column=1, padx=5,sticky=E)
		
		label = tk.Label(root,
								anchor="e", text="Residues: ")
		label.grid(row=9, column=1, padx=5,sticky=E)
		
		label = tk.Label(root,
								anchor="e", text="Repeats: ")
		label.grid(row=10, column=1, padx=5,sticky=E)
		
		label = tk.Label(root,
								anchor="e", text="Single Del or Add: ")
		label.grid(row=11, column=1, padx=5,sticky=E)
		label = tk.Label(root,
								anchor="e", text="Other Matches: ")
		label.grid(row=12, column=1, padx=5,sticky=E)
		
		label = tk.Label(root,
								anchor="e", text=" "*200)
		label.grid(row=1, column=2,rowspan =3)
		
		self.check_mass = tk.Entry(root)
		self.check_mass.grid(row=len(resin_choices)+3, column=0, pady=5, padx=5)
		
		btn = tk.Button(root, width = 15,text="Check Mass", command=self.check_mass_peak)
		btn.grid(row=len(resin_choices)+4, column=0, pady=5)
		
		btn = tk.Button(root, width = 15,text="Clear Sequence", command=onclick('Clear Sequence'))
		btn.grid(row=len(resin_choices)+6, column=0, pady =5)
		
		btn = tk.Button(root, width = 20,text="Open Sequence Builder", command=onclick('Open Builder'))
		btn.grid(row=len(resin_choices)+7, column=0, pady =5)
		
		label = tk.Label(root,
								anchor="e", text="Quick References: ")
		label.grid(row=15, column=1, padx=5,sticky=E)
	
		#Makes a menu for the App
		menubar = Menu(root)
		filemenu = Menu(menubar, tearoff=0)
		filemenu.add_command(label="Load Sequence", command=onclick('Load Sequence'))
		filemenu.add_command(label="Save Sequence", command=onclick('Save Sequence'))
		filemenu.add_command(label="Print Report", command=onclick('Print Report'))
		filemenu.add_separator()
		filemenu.add_command(label="Exit", command=root.destroy)
		menubar.add_cascade(label="File", menu = filemenu)
		
		buildermenu = Menu(menubar, tearoff=0)
		buildermenu.add_command(label="Open Sequence Builder", command=onclick('Open Builder'))
		buildermenu.add_command(label="Add Custom Buttons", command = onclick('Add Custom Buttons'))
		buildermenu.add_command(label="Delete Custom Buttons", command = onclick('Delete Custom Buttons'))
		menubar.add_cascade(label="Open Builder", menu = buildermenu)
		
		pdbmenu = Menu(menubar, tearoff=0)
		pdbmenu.add_command(label="Write Pdb file", command=onclick('Open PDB'))
		#pdbmenu.add_command(label="Totally Tubular", command=onclick('Totally Tubular'))
		pdbmenu.add_command(label="What the Sheet?", command=onclick('What the Sheet'))
		pdbmenu.add_command(label="Convert PDB to String", command=onclick('Convert PDB'))
		menubar.add_cascade(label="Write to PDB", menu = pdbmenu)
		
		convertmenu = Menu(menubar, tearoff=0)
		convertmenu.add_command(label='Remove Neg Log', command=onclick('Remove Neg Log'))
		convertmenu.add_command(label='Buffer Subtract', command=onclick('Subtract Buffer'))
		#convertmenu.add_command(label='Convert PDB to String', command=self.about_program)
		menubar.add_cascade(label="Conversions", menu = convertmenu)
		
		aboutmenu = Menu(menubar, tearoff=0)
		aboutmenu.add_command(label='Help', command=self.help_menu)
		aboutmenu.add_command(label='About', command=self.about_program)
		menubar.add_cascade(label="About", menu = aboutmenu)
		
		root.config(menu=menubar)
		v = StringVar()
		v.set('Chlorotrityl')
		keyboardButtons=[]
		for i, k in enumerate(resin_choices):
			radio = tk.Radiobutton(root, text=k, value=k, variable=v,
								command=onclick(k))
			radio.grid(row=i+1, column=0, padx=5, sticky=W)
			keyboardButtons.append(radio)
	def help_menu(self):
		top = Toplevel()
		top.wm_title('Help Menu')
		help1 = 'To load a peptide/peptoid sequence select File -> Load Sequence. You may load any txt file format' + \
				' containing your sequence of interest. Sequences can be labeled with single letters or triple letters such as "A" or "Ala".' + \
				' Unnatural building blocks can be entered with their corresponding three or more letter codes such as "Hao" or "Phe(I)".'

		help2 = 'To build a new sequence, select "Sequence Builder" or press the "Sequence Builder" button on the main window.' + \
				' This will generate a second window showing all the current building blocks the program can handle. To build a sequence either press ' + \
				'the button corresponding to the desired residue or type a sequence at the bottom of the Sequence Builder window ' + \
				'and hit the "Enter Sequence" button. The sequence is built N -> C. Information about the exact mass, molecular weight, ' +\
				'number of residues, and residue repeats are shown in the main window. These values will automatically adjust while ' +\
				'you build your sequence.'
		
		help6 = 'To add a custom button to the "Sequence Builder" select Open Builder -> Add Custom Buttons. This will open a small window where' +\
				' you may enter a three or four letter code and the molecular formula of your residue. You may use any combination of letters or numbers'+\
				' for the entry code ex. Aw3. For the molecular formula enter the number of carbons (C), hydrogens (H), oxygens (O), nitrogens (N), sulfurs (S), fluorines (F),' +\
				' chlorines (Cl), bromines (Br), and iodines(I). You may enter the molecular formula in any order, ex. C 12 H 13 N 1 O 2 F 3. Enter the molecular formula of' +\
				' the acid and free amine species. For peptoids enter the molecular formula of the fragment including the backbone atoms. Custom buttons will' +\
				' appear on the far right side of the "Sequence Builder" menu.'
			
		help7=	'To delete a custom button select Open Builder -> Delete Custom Buttons. Click on the button you wish to delete then reopen the "Sequence Builder".'
				
		help3 = 'To save a sequence select File -> Save Sequence. This will save your sequence in a txt string format so that you may' +\
				' reload it at another time'
				
		help4 = 'To check a Mass Spec peak you must choose an "Option" on the left side of the main window. If an option is not chosen' +\
				' the program will assume the C-terminal end is an acid. Enter an observed MS peak into the entry window above the ' +\
				'"Check Mass" button. Your entry peak will be compared to the expected values of multiple MS adducts. Matches will be displayed' +\
				' in the main window. The program also tracks the peaks you have entered so that you can print out those peaks and the matches' +\
				' by selecting File -> Print Report.'
				
		help5 = 'To build a "crude" PDB structure based on your sequence you may select "Write to PDB". You will have two simple building options' +\
				' to choose from: alpha-helix or beta-sheet. All 20 natural amino acids can be built in either conformation. Other residues such as ' +\
				'N-Me amino acids, Phe(I), and some of the unnatural amino acids can also be built in either conformation. Peptoids and "Hao" are limited' +\
				'to only a beta-sheet conformation. You may also choose a starting residue number (default is 1) and residue chain ID (default A).' +\
				' The coordinates generated should be further minimized in other programs such as Pymol or Macromodel. **Updated** The Mass Spec Calc ' +\
				'is able to generate nanosheet structures using sequences containing natural amino acids and a handful of peptoids.'
				
		label = tk.Label(top,
							 text=help1, wraplength=500, justify=LEFT)
		label.grid(row=0, column=0, padx=10, pady=10)
		
		label = tk.Label(top,
							 text=help2, wraplength=500, justify=LEFT)
		label.grid(row=1, column=0, padx=10, pady=10)
		label = tk.Label(top,
							 text=help6, wraplength=500, justify=LEFT)
		label.grid(row=2, column=0, padx=10, pady=10)
		label = tk.Label(top,
							 text=help7, wraplength=500, justify=LEFT)
		label.grid(row=3, column=0, padx=10, pady=10)
		label = tk.Label(top,
							 text=help3, wraplength=500, justify=LEFT)
		label.grid(row=4, column=0, padx=10, pady=10)
		label = tk.Label(top,
							 text=help4, wraplength=500, justify=LEFT)
		label.grid(row=5, column=0, padx=10, pady=10)
		label = tk.Label(top,
							 text=help5, wraplength=500, justify=LEFT)
		label.grid(row=5, column=0, padx=10, pady=10)
	
	def about_program(self):
		top = Toplevel()
		top.wm_title('About')
		about1 = 'This program was developed as a tool to aid the everyday peptide/peptoid chemist analyze and visualize their products.'

		about2 = 'The program was written by Ryan K. Spencer.'
		
		about3 = 'Please send any bugs/questions/comments to ryan.spencer.79@gmail.com or rspencer@uci.edu.'
				
		label = tk.Label(top,
							 text=about1, wraplength=500, justify=LEFT)
		label.grid(row=0, column=0, padx=10, pady=10)
		label = tk.Label(top,
							 text=about2, wraplength=500, justify=LEFT)
		label.grid(row=1, column=0, padx=10, pady=10)
		label = tk.Label(top,
							 text=about3, wraplength=500, justify=LEFT)
		label.grid(row=2, column=0, padx=10, pady=10)
	
	def check_mass_peak(self):
			peak = self.check_mass.get()
			user_peaks.append(peak)
			print(user_peaks)
			print (peak)
			peak = float(peak)
			if sequence == None:
				print("Please enter a sequence first")
			else:
				check_value = mass_spec_calc.mass_spec_compare(peak, mass_spec_calc.mass_spec_peaks(exact_mass,final_sequence2))
				check_value_temp = mass_spec_calc.mass_spec_del_add(exact_mass,final_sequence2)
				check_val_add_del=[]
				check_val_add_del.append(mass_spec_calc.mass_spec_compare_del_add(peak,check_value_temp[0]))
				check_val_dd=[]
				check_val_dd.append(mass_spec_calc.mass_spec_compare_del_add(peak,check_value_temp[1]))
				
				add_del = ', '.join(str(e) for e in check_val_add_del[0])
				dd = ', '.join(str(e) for e in check_val_dd[0])
			label = tk.Label(root,
								anchor="e", text=check_value + ' '*50)
			label.grid(row=8, column=2, padx=5,sticky=W)
			label = tk.Label(root,
								anchor="e", text=add_del + ' '*50)
			label.grid(row=11, column=2, padx=5,sticky=W)
			label4 = tk.Label(root,
								anchor="e", text=' '*1500, wraplength=500, justify=LEFT)
			label4.grid(row=13, column=2, padx=5,sticky=W)
			label4 = tk.Label(root,
								anchor="e", text=' '*1500, wraplength=500, justify=LEFT)
			label4.grid(row=13, column=2, padx=5,sticky=W)
			label4 = tk.Label(root,
								anchor="e", text=' '*500 + dd + ' '*500, wraplength=500, justify=LEFT)
			label4.grid(row=12, column=2, padx=5,sticky=W)
	
			
			

standard_AA = ['Ala','Arg','Asn','Asp','Cys','Glu','Gln','Gly','His','Ile','Leu','Lys','Met','Phe','Pro','Ser','Thr','Trp', 'Tyr', 'Val']
single_letters_AA =['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
unnatural_AA = ['Phe(I)', 'Phe(Br)', 'Phe(5F)', 'Orn','Hao', 'Nle', 'Nva', 'Chg','Cha','Pra']
misc = ['Ac', 'Fmoc', 'TFA salt', 'Disulfide', 'Tol']
N_me_AA=['Ala(NMe)', 'Gly(NMe)', 'Ile(NMe)', 'Leu(NMe)', 'Nle(NMe)','Phe(NMe)','Tyr(NMe)', 'Val(NMe)', 'Nva(NMe)']
peptoid_names = ipn.get_peptoid_names()
peptoid_comments = ipn.get_peptoid_comments()

class Sequence_Builder(object):
	'''This is the main class for building a peptide using buttons'''
	def __init__(self, master, **kwargs):
		top = Toplevel()
		top.wm_title('Sequence Builder')
		keyboardButtons =[]
		if len(standard_AA)>10:
			label = tk.Label(top,
								anchor="w", padx=5, pady=5, text="Standard AA")
			label.grid(row=0, column=0, columnspan=2)
			AA1 = standard_AA[0:10]
			AA2 = standard_AA[10:]
			for i, k in enumerate(AA1):
				button = tk.Button(top, text=k, width=7, relief='raised',
							command=onclick(k), state=tk.NORMAL)
				button.grid(row=i+1, padx=5,pady =5, column=0)
				keyboardButtons.append(button)	
			for i, k in enumerate(AA2):
				button = tk.Button(top, text=k, width=7, relief='raised',
							command=onclick(k), state=tk.NORMAL)
				button.grid(row=i+1, padx=5,pady =5,column=1)
				keyboardButtons.append(button)
		else:
			pass
		btn = tk.Button(top, width = 15,text="Delete Last Residue", command=onclick("Delete Last"))
		btn.grid(row=12, column=0, columnspan = 2, pady =5, padx =5)
		
		btn = tk.Button(top, width = 15,text="Enter Sequence", command=onclick('Enter Sequence'))
		btn.grid(row=13, column=0, columnspan = 2, pady =5, padx =5)
		
		
		global enter_sequence
		enter_sequence = tk.Entry(top, width = 80)
		enter_sequence.grid(row=13, column=2, columnspan = 6,pady=5, padx=5)
		
		label = tk.Label(top,
								anchor="w", padx=5, pady=5, text="Unnatural AA")
		label.grid(row=0, column=2)
		
		for i, k in enumerate(unnatural_AA):
				button = tk.Button(top, text=k, width=7, relief='raised',
							command=onclick(k), state=tk.NORMAL)
				button.grid(row=i+1, padx=5,pady =5,column=2)
				keyboardButtons.append(button)
		
		label = tk.Label(top,
								anchor="w", padx=5, pady=5, text="Misc.")
		label.grid(row=0, column=6)
		for i, k in enumerate(misc):
				button = tk.Button(top, text=k, width=7, relief='raised',
							command=onclick(k), state=tk.NORMAL)
				button.grid(row=i+1, padx=5,pady =5,column=6)
				keyboardButtons.append(button)
		
		label = tk.Label(top,
								anchor="w", padx=5, pady=5, text="N-Me AA")
		label.grid(row=0, column=7)
		for i, k in enumerate(N_me_AA):
				button = tk.Button(top, text=k, width=7, relief='raised',
							command=onclick(k), state=tk.NORMAL)
				button.grid(row=i+1, padx=5,pady =5,column=7)
				keyboardButtons.append(button)
				
		
		
		if len(peptoid_names) > 20:
			balloon = Pmw.Balloon()
			peptoid1 = peptoid_names[0:10]
			peptoid_comment1 = peptoid_comments[0:10]
			peptoid2 = peptoid_names[10:20]
			peptoid_comment2 = peptoid_comments[10:20]
			peptoid3 = peptoid_names[20:]
			peptoid_comment3 = peptoid_comments[20:]
			label = tk.Label(top,
								anchor="w", padx=5, pady=5, text="Peptoids")
			label.grid(row=0, column=3, columnspan=3)
			for i, k in enumerate(peptoid1):
				button = tk.Button(top, text=k, width=7, relief='raised',
							command=onclick(k), state=tk.NORMAL)
				button.grid(row=i+1, padx=5,pady =5,column=3)
				keyboardButtons.append(button)
				
				balloon.bind(button, peptoid_comment1[i])
			for i, k in enumerate(peptoid2):
				button = tk.Button(top, text=k, width=7, relief='raised',
							command=onclick(k), state=tk.NORMAL)
				button.grid(row=i+1, padx=5,pady =5,column=4)
				keyboardButtons.append(button)
				balloon.bind(button, peptoid_comment2[i])
			for i, k in enumerate(peptoid3):
				button = tk.Button(top, text=k, width=7, relief='raised',
							command=onclick(k), state=tk.NORMAL)
				button.grid(row=i+1, padx=5,pady =5,column=5)
				keyboardButtons.append(button)
				balloon.bind(button, peptoid_comment3[i])
		else:
			pass
		if os.path.isfile('custom_buttons.txt'):
			print('Opening custom file')
			fobj2 = open('custom_buttons.txt','r')
			custom_names_user = eval(fobj2.readline())
			custom_atoms_user = eval(fobj2.readline())
			fobj2.close()
			if custom_names_user == []:
				pass
			else:
				label = tk.Label(top,
								anchor="w", padx=5, pady=5, text="Custom")
				label.grid(row=0, column=8) 
				for i,k in enumerate(custom_names_user):
					button = tk.Button(top, text=k, width=7, relief='raised',
								command=onclick(k), state=tk.NORMAL)
					button.grid(row=i+1, padx=5,pady =5,column=8)
					keyboardButtons.append(button)
				exactMassAA.update(custom_atoms_user)
		else:
			print('Did not open file')
			pass

class buffer_subtract(object):
	def __init__(self,master, **kwargs):
		top = Toplevel()
		top.wm_title('Subtract Buffer')
		btn = tk.Button(top, text="Select Sample SAXS", relief='raised',
						command=self.Sample(fname='Sample'), state=tk.NORMAL)
		btn.grid(row=2, column=0, columnspan=2,pady =5, padx =5)
		
		btn = tk.Button(top, text="Select Buffer SAXS", relief='raised',
						command=self.Sample(fname='Buffer'), state=tk.NORMAL)
		btn.grid(row=3, column=0, columnspan=2,pady =5, padx =5)
		global adjustment
		adjustment = tk.Entry(top, width = 10)
		adjustment.grid(row=4, column=0, pady=5, padx=5)
		btn = tk.Button(top, text="Write Buffer Subtraction SAXS", relief='raised',
						command=self.Sample(fname='Subtract'), state=tk.NORMAL)
		btn.grid(row=5, column=0, columnspan=2,pady =5, padx =5)
		
	def Sample(self,fname):
		def sample_sub():
			print('Clicked Something')
			if fname == 'Sample':
				fnameSample = askopenfilename(filetypes=(("dat", "*.dat"),
														("ndat","*.ndat"),
														("All files", "*.*") ),
														defaultextension='.dat')
				print('Opening .dat file')
				try:
					global fobjSample
					fobjSample = open(fnameSample, 'r')
					fnameSample.close()
				except:
					pass
			elif fname == 'Buffer':
				fnameBuffer = askopenfilename(filetypes=(("dat", "*.dat"),
														("ndat","*.ndat"),
														("All files", "*.*") ),
														defaultextension='.dat')
				print('Opening .dat file')
				try:
					global fobjBuffer
					fobjBuffer = open(fnameBuffer, 'r')
					fnameBuffer.close()
				except:
					pass
			elif fname == 'Subtract':
				Sample_dict={}
				Buffer_dict={}
				fnameSubtract = asksaveasfilename(filetypes=(("dat", "*.dat"),
														("ndat","*.ndat"),
														("All files", "*.*") ),
														defaultextension='.dat')
				try:
					fobjSubtract = open(fnameSubtract, 'w')
				except:
					pass
				i = 0
				for line in fobjSample:
					line_strip = line.strip()
					line_split=line_strip.split()
					Sample_dict.update({i:line_split})
					i+=1
				i = 0
				fobjSample.seek(0)
				for line in fobjBuffer:
					line_strip = line.strip()
					line_split=line_strip.split()
					Buffer_dict.update({i:line_split})
					i+=1
				fobjBuffer.seek(0)
				k = 0
				
				while k <i:
					Sample_temp = Sample_dict.get(k)
					Buffer_temp = Buffer_dict.get(k)
					difference = float(Sample_temp[1]) - float(Buffer_temp[1])*float(adjustment.get())
					fobjSubtract.write(Sample_temp[0] + '\t' + str(difference) + '\t' + Sample_temp[2] + '\n')
					k +=1
					
				fobjSubtract.close()	
		return sample_sub

class Write_nanosheet(object):
	def __init__(self, master, **kwargs):
		top = Toplevel()
		top.wm_title('Generate Nanosheet')
		label = tk.Label(top,
					anchor = "w", padx=5, pady=5, text="Please Enter Appropriate Values")
		label.grid(row = 0, column = 0, columnspan = 2)
		options=[]
		Choices = [("Parallel", "para"),
				   ("Antiparallel", "anti")]
		sel = StringVar()
		sel.set("para")
		i = 0
		global sheet_choice
		sheet_choice = ['para']
		for text, choice in Choices:
			radio = tk.Radiobutton(top, text=text, variable=sel, value=choice, command = self.write_sheets(choice))
			radio.grid(row=1, column=i+1, padx=5, sticky=W)
			options.append(radio)
			i += 1
		label = tk.Label(top,
					anchor = "w", padx=5, pady=5, text="Sheet length # units (ex. 3): ")
		label.grid(row = 2, column = 0, sticky = E)
		global length_strands
		length_strands = tk.Entry(top, width = 7)
		length_strands.grid(row=2, column=1, pady=5, padx=5)
		label = tk.Label(top,
					anchor = "w", padx=5, pady=5, text="Sheet width # units (ex. 10): ")
		label.grid(row = 3, column = 0, sticky=E)
		global width_strands
		width_strands = tk.Entry(top, width = 7)
		width_strands.grid(row=3, column=1, pady=5, padx=5)
		label = tk.Label(top,
					anchor = "w", padx=5, pady=5, text="Distance between layers (ex. 20): ")
		label.grid(row = 4, column = 0, sticky=E)
		global height_strands
		height_strands = tk.Entry(top, width = 7)
		height_strands.grid(row=4, column=1, pady=5, padx=5)
		global adjust_l
		adjust_l = tk.Entry(top, width = 7)
		adjust_l.grid(row=5, column=1, pady=5, padx=5)
		label = tk.Label(top,
					anchor = "w", padx=5, pady=5, text="Distance between monomer ends (ex. 5): ")
		label.grid(row = 5, column = 0, sticky=E)
		btn = tk.Button(top, text="Write Nanosheet PDB", relief='raised',
						command=self.write_sheets(choice='Write Sheet'), state=tk.NORMAL)
		btn.grid(row=6, column=0, columnspan=2,pady =5, padx =5)
	
	def write_sheets(self, choice):
		def sheets():
			if choice == 'para':
				sheet_choice[0]='para'
				print('You selected para')
			elif choice == 'anti':
				sheet_choice[0]='anti'
				print('You selected anti')
			elif choice == 'Write Sheet':
				length=length_strands.get()
				width=width_strands.get()
				height=height_strands.get()
				l_adjust = adjust_l.get()
				if l_adjust == '':
					l_adjust = 0
				else:
					pass
				if sheet_choice[0] == 'para':
					wPdb.what_the_sheet_par(sequence_temp,length,height,width, l_adjust)
				elif sheet_choice[0] == 'anti':
					wPdb.what_the_sheet_anti(sequence_temp,length,height,width, l_adjust)
			else:
				pass
				
		return sheets

	
class Write_to_Pdb(object):
	'''Opens an additional window for writing sequence to a PDB'''
	def __init__(self, master, **kwargs):
		top = Toplevel()
		top.wm_title('Write to PDB')
		label = tk.Label(top,
					anchor = "w", padx=5, pady=5, text="Please select an Option")
		label.grid(row = 0, column = 0, columnspan = 2)
		keyboardButtons2=[]
		Choices = [("Alpha-Helix", "alpha"),
				   ("Beta-Sheet", "beta")]
		sel = StringVar()
		sel.set("alpha")
		i = 1
		global PDB_choice
		PDB_choice = ['alpha']
		for text, choice in Choices:
			radio = tk.Radiobutton(top, text=text, variable=sel, value=choice, command = self.set_sele(choice))
			radio.grid(row=i+1, column=0, padx=5, sticky=W)
			keyboardButtons2.append(radio)
			i += 1
		label = tk.Label(top,
					anchor = "w", padx=5, pady=5, text="Starting Residue #: ")
		label.grid(row = 4, column = 0, sticky = E)
		label = tk.Label(top,
					anchor = "w", padx=5, pady=5, text="Chain ID (A,B,C..): ")
		label.grid(row = 5, column = 0, sticky=E)
		
		global start_residue
		start_residue = tk.Entry(top, width = 7)
		start_residue.grid(row=4, column=1, pady=5, padx=5)
		
		global chain_ID
		chain_ID = tk.Entry(top, width = 7)
		chain_ID.grid(row=5, column=1, pady=5, padx=5)
		
		btn = tk.Button(top, text="Write Sequence to PDB", relief='raised',
						command=self.set_sele(choice='Write PDB'), state=tk.NORMAL)
		btn.grid(row=6, column=0, columnspan=2,pady =5, padx =5)
		
	
	def set_sele(self,choice):
		def sele():
			if choice == 'alpha':
				PDB_choice[0] = 'alpha'
				print("You Selected Alpha")
			elif choice == 'beta':
				PDB_choice[0] = 'beta'
				print("You selected Beta")
			elif choice == 'Write PDB':
				pdb_sequence = ''.join(str(e) for e in sequence_temp[1:])
				print(start_residue.get())
				print(chain_ID.get())
				sr = start_residue.get()
				if PDB_choice[0] == 'alpha':
					print("Writing alpha-helix")
					wPdb.write_alpha_helix(convert.triple_to_single(pdb_sequence), sr, chain_ID.get())
				elif PDB_choice[0] == 'beta':
					print("Writing beta-sheet")
					wPdb.write_beta_sheet(convert.triple_to_single(pdb_sequence), sr, chain_ID.get())
				else:
					pass
			else:
				pass
		return sele
		
class delete_custom_buttons(object):
	def __init__(self, master, **kwargs):
		top = Toplevel()
		if os.path.isfile('custom_buttons.txt'):
			fobj2 = open('custom_buttons.txt','r')
			custom_names2 = eval(fobj2.readline())
			custom_atoms2 = eval(fobj2.readline())
			fobj2.close()
			label = tk.Label(top,
				anchor = "w", padx=5, pady=5, text="Delete Custom Buttons")
			label.grid(row = 0, column = 0, columnspan =3)
			label = tk.Label(top,
				anchor = "w", padx=5, pady=5, text="Please click once on the fragment you wish to delete.")
			label.grid(row = 1, column = 0, columnspan =3)
			for i, k in enumerate(custom_names2):
				button = tk.Button(top, text=k, width=7, relief='raised',
							command= self.delete_button(k,custom_names2,custom_atoms2),state=tk.NORMAL)
				button.grid(row=i+2, padx=5,pady =5,column=0)
		else:
			label = tk.Label(top,
				anchor = "w", padx=5, pady=5, text="You have no custom buttons to delete.")
			label.grid(row = 0, column = 0)
	
	def delete_button(self,k, custom_names2, custom_atoms2):
		def del_button():
			if k in custom_names2:
				custom_names2.remove(k)
				print('Removing Button')
				fname = 'custom_buttons.txt' 
				fobj = open(fname,'w')
				fobj.write('[')
				count = 0
				if len(custom_names2) == 1:
					fobj.write('\''+custom_names2[0]+'\'')
				else:
					for m,n in enumerate(custom_names2):
						if count < len(custom_names2)-1:
							fobj.write('\''+str(n)+'\'')
							fobj.write(',')
						else:
							fobj.write('\''+str(n)+'\'')
						count+=1
				fobj.write(']\n')
				count = 0
				fobj.write('{')	
				for m,n in enumerate(custom_names2):
					key = n
					value = custom_atoms2.get(key)
					fobj.write('\''+key +'\':[')
					for i,j in enumerate(value):
						if i < len(value)-1:
							fobj.write(str(j)+',')
						else:
							fobj.write(str(j))
					fobj.write(']')
					if count < len(custom_names2)-1:
						fobj.write(',')
					else:
						pass
					count+=1
				fobj.write('}')
				fobj.close()
		return del_button
		
class convert_PDB_to_string(object):
	'''This converts a formated PDB file to a single string containing \n'''
	def __init__(self, master, **kwargs):
		try:
			fname = askopenfilename(filetypes=(("pdb", "*.pdb"),
                                           ("All files", "*.*") ))
			fobj = open(fname+'converted_PDB.txt','w')
		except:
			pass
		initial_string=[]
		fobj.write('\'')
		with open(fname) as openfileobject:
			for line in openfileobject:
				line2 = line.replace("\n","")
				initial_string.append(line2)
		full_string = ('\\'+'n').join(str(e) for e in initial_string)
		fobj.write(full_string)
		fobj.write('\'')
		fobj.close()
		openfileobject.close()

class test(object):
	def __init__(self):
		pass

class add_custom_buttons(object):
	'''This is the subprogram that allows you to add additional residue fragments'''
	def __init__(self, master, **kwargs):
		top = Toplevel()
		try:
			if os.path.isfile('custom_buttons.txt'):
				print('Found the custom buttons file, appending.')
				fobj2 = open('custom_buttons.txt','r')
					
			else:
				print('Did not find the file, making one')
				fname = 'custom_buttons.txt' 
				fobj = open(fname,'w')
				fobj.write('[]\n')
				fobj.write('{}\n')
				fobj.close()
				fobj2 = open(fname, 'r')
		except:
			print('Could not Open/Create a file.')
		global custom_names
		custom_names = eval(fobj2.readline())
		global custom_atoms
		custom_atoms = eval(fobj2.readline())
		fobj2.close
		#build GUI for custom entries
		label = tk.Label(top,
					anchor = "w", padx=5, pady=5, text="Add Custom Residues")
		label.grid(row = 0, column = 0, columnspan=2) 
		label = tk.Label(top,
					anchor = "w", padx=5, pady=5, text="Desired Residue Code: ")
		label.grid(row = 1, column = 0, sticky = E) 
		global letter_ID
		letter_ID = tk.Entry(top, width = 7)
		letter_ID.grid(row=1, column=1, pady=5, padx=5)
		submit_button = tk.Button(top, padx =5, pady =5, text = "Submit Custom Residue", command = self.add_custom)
		submit_button.grid(row=3, column=0, pady=5, padx=5, columnspan=2)
		global mol_form_new
		mol_form_new = tk.Entry(top, width = 15)
		mol_form_new.grid(row=2, column=1, pady=5, padx=5)
		label = tk.Label(top,
					anchor = "w", padx=5, pady=5, text="Molecular Formula: ")
		label.grid(row = 2, column = 0, sticky = E)
		
	
	def add_custom(self):
		custom_entry = letter_ID.get()
		if custom_entry in custom_names or custom_entry in standard_AA or custom_entry in unnatural_AA or custom_entry in misc or custom_entry in N_me_AA or custom_entry in peptoid_names or custom_entry in single_letters_AA:
			self.popupmsg('Entry already exists, try another name.')
		else:
			print(type(custom_names))
			new_entry = str(custom_entry)
			custom_names.append(new_entry)
			custom_formula = mol_form_new.get()
			chopped_formula = re.findall('[A-Z][^A-Z]*', custom_formula.replace(' ',''))
			print(chopped_formula)
			carbons =0; hydrogens =0; nitrogens =0; oxygens = 0; sulfurs = 0; fluorines = 0; chlorines=0; bromines = 0; iodines = 0
			for i, number in enumerate(chopped_formula):
				#further_chopped = re.findall('([A-Z])([0-9]+)', number)
				further_chopped = re.findall('([a-zA-Z])([0-9]+)', number)
				print(further_chopped)
				if further_chopped[0][0] == 'C':
					carbons = int(further_chopped[0][1])
				elif further_chopped[0][0] == 'H':
					hydrogens = int(further_chopped[0][1])-2
				elif further_chopped[0][0] == 'N':
					nitrogens = int(further_chopped[0][1])
				elif further_chopped[0][0] == 'O':
					oxygens = int(further_chopped[0][1])-1
				elif further_chopped[0][0] == 'S':
					sulfurs = int(further_chopped[0][1])
				elif further_chopped[0][0] == 'F':
					fluorines = int(further_chopped[0][1])
				elif further_chopped[0][0] == 'l':
					chlorines = int(further_chopped[0][1]) 
				elif further_chopped[0][0] == 'r':
					bromines = int(further_chopped[0][1]) 
				elif further_chopped[0][0] == 'I':
					iodines = int(further_chopped[0][1])
				else:
					pass
			exact_mass_custom = {str(custom_entry):[carbons, hydrogens, nitrogens, oxygens, sulfurs, fluorines, chlorines, bromines, iodines]}
			print(chopped_formula)
			custom_atoms.update(exact_mass_custom)
			print(custom_atoms)
			custom_atoms2 = {key:value for key, value in custom_atoms.items()
								if value != 'temp'}
			if 'temp' in custom_names:
				custom_names.remove('temp')
			else:
				pass
			print (custom_names)
			fname = 'custom_buttons.txt' 
			fobj = open(fname,'w')
			fobj.write('[')
			count = 0
			if len(custom_names) == 1:
				fobj.write('\''+str(custom_names[0])+'\'')
			else:
				for m,n in enumerate(custom_names):
					if count < len(custom_names)-1:
						fobj.write('\''+str(n)+'\'')
						fobj.write(',')
					else:
						fobj.write('\''+str(n)+'\'')
					count+=1
			fobj.write(']\n')
			count = 0
			fobj.write('{')	
			for key, value in custom_atoms2.items():
			
				fobj.write('\''+key +'\':[')
				for i,j in enumerate(value):
					if i < len(value)-1:
						fobj.write(str(j)+',')
					else:
						fobj.write(str(j))
				fobj.write(']')
				if count < len(custom_atoms2)-1:
					fobj.write(',')
				else:
					pass
				count+=1
			fobj.write('}')
			fobj.close()
			print(custom_entry)
	def popupmsg(self,msg):
		popup = tk.Tk()
		popup.wm_title("Residue Not Added")
		label = tk.Label(popup, text=msg)
		label.pack(side="top", fill="x", pady=10)
		B1 = tk.Button(popup, text="Okay", command = popup.destroy)
		B1.pack()
		popup.mainloop()
		
class Totally_tubular(object):
		def __init__(self):
			top = Toplevel()
			label = tk.Label(top,
					anchor = "w", padx=5, pady=5, text="Totally Tubular")
			
	
		
	
		
	
sequence_temp=['Acid']			


def onclick(k):
	def click():
		global sequence_temp
		if k == 'Chlorotrityl' or k == 'Wang' or k == 'Expressed':
			sequence_temp[0] = 'Acid'
		elif k == 'Rink':
			sequence_temp[0] = 'Amide'
		elif k == 'Cyclic':
			sequence_temp[0] = 'Cyclic'
		elif k == 'Open Builder':
			Sequence_Builder('top')
		elif k == 'Add Custom Buttons':
			add_custom_buttons('top')
		elif k == 'Delete Custom Buttons':
			delete_custom_buttons('top')
		elif k == 'Remove Neg Log':
			remove_neg_log()
		elif k == 'Subtract Buffer':
			buffer_subtract('top')
		elif k == 'Open PDB':
			Write_to_Pdb('top')
		elif k == 'Convert PDB':
			convert_PDB_to_string('top')
		elif k == 'Totally Tubular':
			wPdb.write_totally_tubular_par(sequence_temp)
		elif k == 'What the Sheet':
			Write_nanosheet('top')
		elif k == 'Enter Sequence':
			sequence_temp = sequence_temp + convert.sequence_to_list(enter_sequence.get())
		elif k == 'Print Report':
			report.write_Report(sequence_temp,user_peaks)
		elif k == 'Delete Last':
			i = len(sequence_temp)-1
			if i>0:
				del sequence_temp[i]
				label = tk.Label(root,
								anchor="e", text=' '*200)
				label.grid(row=1, column=2)
				label = tk.Label(root,
								anchor="e", text=' '*200)
				label.grid(row=2, column=2)
				label = tk.Label(root,
								anchor="e", text=' '*200)
				label.grid(row=3, column=2)
				label = tk.Label(root,
								anchor="e", text=' '*200)
				label.grid(row=4, column=2)
			else:
				pass
		elif k == 'Save Sequence':
			fname = asksaveasfilename(filetypes=(("txt", "*.txt"),
                                           ("All files", "*.*") ),
									  defaultextension='.txt')
			try:							   
				fobj = open(fname, 'w')
				fobj.write(str(sequence_temp[1:]).replace('[','').replace('\'','').replace(']','').replace(',',''))
				fobj.close()
			except:
				pass
		elif k == 'Clear Sequence':
			i = len(sequence_temp)-1
			while i>0:
				del sequence_temp[i]
				i -=1
			label = tk.Label(root,
								anchor="e", text=' '*200)
			label.grid(row=1, column=2)
			label = tk.Label(root,
								anchor="e", text=' '*200)
			label.grid(row=2, column=2)
			label = tk.Label(root,
								anchor="e", text=' '*200)
			label.grid(row=3, column=2)
			label = tk.Label(root,
								anchor="e", text=' '*200)
			label.grid(row=4, column=2)
			label = tk.Label(root,
								anchor="e", text=' '*200)
			label.grid(row=9, column=2)
			label = tk.Label(root,
								anchor="e", text=' '*200)
			label.grid(row=10, column=2)
			label = tk.Label(root,
								anchor="e", text=' '*200)
			label.grid(row=11, column=2)
			label = tk.Label(root,
								anchor="e", text=' '*200)
			label.grid(row=12, column=2)
			label = tk.Label(root,
								anchor="e", text=' '*200)
			label.grid(row=13, column=2)
			label = tk.Label(root,
								anchor="e", text=' '*200)
			label.grid(row=14, column=2)
			
		elif k == 'Load Sequence':
			fname = askopenfilename(filetypes=(("txt", "*.txt"),
                                           ("All files", "*.*") ))
			if fname:
				try:
					fobj = open(fname)
					sequence_temp = sequence_temp + convert.sequence_to_list((str(fobj.read().strip())))
					fobj.close()
				except:                 
					showerror("Open Source File")
		else:
			sequence_temp.append(k)
		
		sequence = ''.join(str(e) for e in sequence_temp)
		converted_sequence = convert.triple_to_single(sequence)
		global final_sequence2
		final_sequence2 = []
		final_sequence2 = converted_sequence
		
		mol_formula = list(mass_spec_calc.mol_formula(converted_sequence))

		pretty_mol_formula = pretty_formula.make_pretty(mol_formula)
		
		pretty_sequence = '-'.join(str(e) for e in converted_sequence[1:])
		mol_formula_pretty = ' '.join(str(e) for e in pretty_mol_formula)
		global exact_mass
		exact_mass = mass_spec_calc.exact_mass_calc(mol_formula)
		exact_mass_formatted = str(format(exact_mass, '.5f'))
		H = 1.007825; Na = 22.989770; K = 38.963708
		quick_reference = '[1M + 1H]+1 : ' + str(format((exact_mass + H), '.5f')) + ' '*5 + ' [1M + 1Na]+1 : ' + str(format((exact_mass + Na),'.5f')) +\
		' '*5 + ' [1M + 2H]+2 : ' + str(format((exact_mass+2*H)/2, '.5f')) + ' '*5 + '[1M + 1H + 1Na]+2 : ' + str(format((exact_mass+H+Na)/2, '.5f')) +\
		' '*5 +'[1M + 3H]+3 : ' + str(format((exact_mass+3*H)/3, '.5f')) + ' '*5 + '[1M + 2H + 1Na]+3 : ' + str(format((exact_mass+2*H+Na)/3, '.5f')) +\
		' '*5 + '[1M + 1H + 2Na]+3 : ' + str(format((exact_mass+H+2*Na)/3, '.5f'))+' '*5 + '[1M + 4H]+4 : ' + str(format((exact_mass+4*H)/4, '.5f')) +\
		' '*5 + '[1M + 3H + 1Na]+4 : ' + str(format((exact_mass+3*H+Na)/4, '.5f'))
		
		mol_weight = mass_spec_calc.molecular_weight(mol_formula)
		mol_weight_formatted = str(format(mol_weight,'.5f'))
		count = -1
		for i in sequence_temp:
			count +=1
		
		sequence_count = str(count)
		
		polymer = converted_sequence[1:]
		p = dict((i,polymer.count(i)) for i in polymer)
		polymer_nice = str(p)[1:-1]
		mol_form_label = tk.Label(root,
								anchor="w", pady=5, text=str(polymer_nice)+" "*60)
		mol_form_label.grid(row=10, column=2, sticky=W)
		
		label = tk.Label(root,
								anchor="e", text=pretty_sequence, wraplength=500, justify=LEFT)
		label.grid(row=1, column=2,sticky=W, rowspan=3)
		
		label = tk.Label(root,
								anchor="e", text=mol_formula_pretty + ' '*100)
		label.grid(row=5, column=2, padx=5,sticky=W)
		
		label = tk.Label(root,
								anchor="e", text=exact_mass_formatted + ' '*100)
		label.grid(row=6, column=2, padx=5,sticky=W)
		
		label = tk.Label(root,
								anchor="e", text=mol_weight_formatted + ' '*100)
		label.grid(row=7, column=2, padx=5,sticky=W)
		
		label = tk.Label(root,
								anchor="e", text=sequence_count + ' '*100, wraplength=500)
		label.grid(row=9, column=2, padx=5,sticky=W)
		
		label = tk.Label(root,
								anchor="e", text=quick_reference, wraplength = 475)
		label.grid(row=15, column=2, padx=5,  columnspan=2, sticky=W)
		
	
	return click




print ("Exact Mass: " + format(exact_mass,'.5f'))




try:
    type(peak) == float
except:
    print("Please enter a real number")


root = Tk()

app = MassSpecApp(root, title='Mass Spec Calc v1.3')
root.wm_title("Ultimate Awesome of Awesomeness Mass Spec. Calc. + Builder v1.3")
root.mainloop()