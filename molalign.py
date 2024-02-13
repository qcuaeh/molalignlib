# MolAlignLib
# Copyright (C) 2022 José M. Vásquez

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from os import path
from math import inf
from argparse import ArgumentParser, SUPPRESS
from ase.io import read, write
from molalignlib import remap_atoms

def molalign():

    parser = ArgumentParser()
    parser.add_argument('file', nargs='*')
    parser.add_argument('-remap', action='store_true')
    parser.add_argument('-mirror', action='store_true')
    parser.add_argument('-stats', action='store_true')
    parser.add_argument('-fast', action='store_true')
    parser.add_argument('-mass', action='store_true')
    parser.add_argument('-test', action='store_true')
    parser.add_argument('-rec', type=int, default=1)
    parser.add_argument('-trials', type=int, default=SUPPRESS)
    parser.add_argument('-count', type=int, default=SUPPRESS)
    parser.add_argument('-tol', type=float, default=SUPPRESS)
    parser.add_argument('-out', type=str, default='aligned.xyz')
    args = parser.parse_args()

    if len(args.file) == 0:
        raise SystemExit('Error: Missing arguments')
    if len(args.file) == 1:
        raise SystemExit('Error: Too few arguments')
    elif len(args.file) == 2:
        mol0 = read(args.file[0], index=0)
        mol1 = read(args.file[1], index=0)
    else:
        raise SystemExit('Error: Too many arguments')

    if not path.splitext(args.out)[1]:
        raise SystemExit('Error: Output file must have an extension')

    if args.mass:
        mol0.get_weights = mol0.get_masses
        mol1.get_weights = mol1.get_masses

    if args.mirror:
        mol1.mirror_x()

    if args.remap:
        assignments = remap_atoms(
            mol0,
            mol1,
            **vars(args)
        )
        write(args.out, mol0, comment='Reference')
        minrmsd = inf
        for i, assignment in enumerate(assignments, start=1):
            auxmol = mol1[assignment.order]
            rmsd = auxmol.align_to(mol0)
            minrmsd = min(minrmsd, rmsd)
            comment = 'RMSD {:.4f}'.format(rmsd)
            write(args.out, auxmol, append=True, comment=comment)
        if not args.stats:
            print('Optimized RMSD = {:.4f}'.format(minrmsd))
    else:
        rmsd = mol1.align_to(mol0)
        print('RMSD = {:.4f}'.format(rmsd))
        write(args.out, mol0, comment='Reference')
        comment='RMSD {:.4f}'.format(rmsd)
        write(args.out, mol1, append=True, comment=comment)

if __name__ == '__main__':
    molalign()
