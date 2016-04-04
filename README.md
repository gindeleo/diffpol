# diffpol
A small tool to calculate the polarization profile in ABO3 perovskites from a CASTEP .cell files (or .cif files).

The polarization is calculated for each individual B-cation centered conventional unit cell based on atomic displacements and the Born effectrive charges. Born effective charges and cation names may need to be updated or added. There is no library included at the moment.
The script works with defects (vacancies) but no polarization is calculated for defective unit cells. 

For .cif files use cif2cell to convert to the CASTEP .cell format.
https://sourceforge.net/projects/cif2cell/

created by Oliver Gindele  (UCL)

