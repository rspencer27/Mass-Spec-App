from main.massCalc import *

class write_Pdb(object):
	
	def __init__(self,sequence):
		pass
	
	def write_pdb_to_file(self,sequence):
		print(sequence)
		atoms = 0
		links = 0
		positions_temp = []
		positions = []
		locations = []
		locations_temp=[]
		for i,k in enumerate(sequence):
			atoms += writeMolAA.get(k)['atoms']
			links += writeMolAA.get(k)['links']
			positions_temp = writeMolAA.get(k)['positions']['atom']
			locations_temp = writeMolAA.get(k)['positions']['locations']
			positions = positions + positions_temp
			locations = locations_temp*i + locations
		
		print(atoms)
		print(links)
		print(positions)
		print(locations)