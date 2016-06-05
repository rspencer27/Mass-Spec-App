from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename


class create_CONECT_table(object):
	def __init__(self):
		fname = askopenfilename(filetypes=(("pdb", "*.pdb"),
                                           ("All files", "*.*") ),
								  defaultextension='.pdb')
		try:
			fobj = open(fname, 'r')
			with open(fname) as openfileobject:
				for line in openfileobject:
				line2 = line.replace("\n","")
				initial_string.append(line2)
			line = fobj.readline()
			totat_input.append(line)
			