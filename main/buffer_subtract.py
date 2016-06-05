from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename
import re


class remove_neg_log(object):
	def __init__(self):
		
		fname = askopenfilename(filetypes=(("dat", "*.dat"),
										   ("ndat","*.ndat"),
                                           ("All files", "*.*") ),
								  defaultextension='.dat')
		print('Opening .dat file')
		try:
			fobj = open(fname, 'r')
			fobj2 = open(fname+'neg_remove.dat','w')
			initial_string=[]
			with open(fname) as openfileobject:
				for line in openfileobject:
					print(line)
					line2 = line.strip()
					print(type(line2))
					initial_string=line2.split()
					print(initial_string)
					print (type(initial_string))
					if float(initial_string[1]) > 0:
						fobj2.write(initial_string[0] + '     ' + initial_string[1] + '     ' + initial_string[2] + '\n')
					else:
						pass
			fobj.close()
			fobj2.close()
		except:
			pass
			