import tkinter as tk
from tkinter import *

AA_names1 = ['Dof','Arg','Asn','Asp','Cys','Glu','Gln','Gly','His','Ile']
AA_names2 = ['Leu','Lys','Met','Phe','Pro','Ser','Thr','Trp', 'Tyr', 'Val']

class Sequence_Builder2(object):
	'''This is the main class for building a peptide using buttons'''
	def __init__(self, master, **kwargs):
		top = Toplevel()
		keyboardButtons =[]
		label =  tk.Label(top,
								anchor="e", text="Standard Amino Acids")
		label.grid(row=1, column=0,columnspan=2, sticky=NSEW)
		for i, k in enumerate(AA_names1):
			button = tk.Button(top, text=k, width=7, relief='raised',
								command=onclick(k), state=tk.NORMAL)
			button.grid(row=i+1, column=0)
			keyboardButtons.append(button)