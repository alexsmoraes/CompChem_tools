#!/usr/bin/env python3

from pathlib import Path
import subprocess
import os
import numpy as np
from shutil import copy2, move, rmtree
import argparse

parser = argparse.ArgumentParser(description='Gets the N best crest conformers and optimize them using xTB at GFN2 level. Also saves energies.csv with all files and energies.',
                                  epilog='Written by Dr. Alex S. Moraes')
parser.add_argument('--fol', metavar='FOLDER', type=str, default='crest_best_xtbopt', help='specify output folder name (default=crest_best_xtbopt)')
parser.add_argument('-n', '--nbest', metavar='N', type=int, default=5, help='specify the number of best conformers to get (default=5)')
parser.add_argument('-c', '--chrg', metavar='CHRG', type=int, default=0, help='specify the charges of the molecules (default=0)')
parser.add_argument('--uhf', metavar='UHF', type=int, default=0, help='specify the number of unpaired electrons for Unrestricted Hartree-Fock calculation (default=0)')
parser.add_argument('-s', '--solvent', metavar='SOL', help='specify the solvent (check xTB documentation for solvent options) (default=None)')
parser.add_argument('--keepout', action='store_true', help='keeps xtb outputs')
parser.add_argument('--crest_fol', metavar='C_FOL', help='specify the folder with the crest files to analyze (default=current folder)')

args = parser.parse_args()

n_best = args.nbest
outfol = Path(args.fol)
chrg = args.chrg
uhf = args.uhf 
sol = args.solvent
keepout = args.keepout
crestfol = args.crest_fol

Eh2eV = 27.2114 # Conversion factor from hartree to eletronvolt

if crestfol:
	allcrest = sorted(list(Path(crestfol).rglob('crest_conformers.xyz')))
else:
	allcrest = sorted(list(Path().rglob('crest_conformers.xyz')))

energies = ['xyzfile,energy(Eh),energy(eV)\n']

allxyz = outfol / 'allxyz'
allbest = outfol / 'allbest'
allout = outfol / 'allout'

if outfol.is_dir():
	rmtree(outfol)

os.mkdir(outfol)
os.mkdir(allxyz)
os.mkdir(allbest)
os.mkdir(allout)

for mol_file in allcrest:
	curr_path = Path()
	mol_name = os.path.basename(mol_file.parent)
	mol = open(mol_file).readlines()
	n_atoms = int(mol[0].strip())
	idx_i = 0
	idx_f = idx_i + n_atoms + 2
	en = 1e20
	for i in range(n_best):
		if len(mol[idx_i:idx_f]) != 0:
			nmol = str(i+1).zfill(2)
			name = f'{mol_name}_{nmol}.xyz'
			with open(allxyz / name,'w+') as f:
				for line in mol[idx_i:idx_f]:
					f.write(line)
			idx_i = idx_f
			idx_f = idx_i + n_atoms + 2
			outname = allout / name.replace('xyz','out')
			if sol:
				cmd = f"xtb {allxyz / name} --opt --chrg {chrg} --uhf {uhf} --gfn 2 -alpb {sol} > {outname}"
			else:
				cmd = f"xtb {allxyz / name} --opt --chrg {chrg} --uhf {uhf} --gfn 2 > {outname}"
			print(cmd)
			opt_xtb_gfn2 = subprocess.run([cmd], shell=True, capture_output=True)
			if not keepout:
				os.remove(outname)
			aux = float(open('xtbopt.xyz').readlines()[1].split()[1].strip())
			if aux < en:
				en = aux
				best = name.replace(f'{nmol}.', 'xtb_best.')
			energies.append(f'{name},{en},{en*Eh2eV}\n')
			os.rename('xtbopt.xyz',name)
			move(name, allxyz / name)
	copy2(allxyz / name, allbest / best)

if not keepout:
	os.rmdir(allout)

for file in ['charges', 'wbo', 'xtbopt.log', 'xtbrestart', 'xtbtopo.mol']:
	os.remove(file)

with open(outfol / 'energies.csv','w+') as f:
	for line in energies:
		f.write(line)
