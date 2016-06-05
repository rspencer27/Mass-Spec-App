# Mass-Spec-App
Python code for the Mass Spec App GUI

This is the Python source code for the Mass Spec App currently in developemet. 

Current App calculates the exact mass and molecular weight of a peptide/peptoid sequence. Multiply charged adducts are calculated and compared to user input. 

Crude models can be written in PDB format directly from the input sequence. Users may output the sequence as an alpha-helix or beta-sheet. **NOTE** not all coordinates are available for PDB generation. 

Two setup.py files are included to create executable files for either Windows or Mac. To use, run the appropriate setup.py using py2exe (Windows) or py2app (Mac).

The GUI is based on the Tkinter libraries. Python code is written in Python 3.

Future developement plans:

1) Create an interface for users to create custom residues from PDB coordinates.

2) Add additional comments to the code for clarity.
