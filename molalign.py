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
from molalignlib import assign_atoms

def molalign():

    parser = ArgumentParser()
    parser.add_argument('filelist', nargs='+')
    parser.add_argument('-sort', dest='reorder', action='store_true', default=SUPPRESS)
    parser.add_argument('-mass', action='store_true')
    parser.add_argument('-mirror', action='store_true')
    parser.add_argument('-stats', action='store_true')
    parser.add_argument('-fast', action='store_true')
    parser.add_argument('-test', action='store_true')
    parser.add_argument('-rec', type=int, default=1)
    parser.add_argument('-trials', type=int, default=SUPPRESS)
    parser.add_argument('-tol', type=float, default=SUPPRESS)
    parser.add_argument('-count', type=int, default=SUPPRESS)
    parser.add_argument('-out', type=str, default='aligned.xyz')
    args = parser.parse_args()

    if len(args.filelist) == 1:
        print('Error: Too few arguments')
    elif len(args.filelist) == 2:
        atoms0 = read(args.filelist[0], index=0)
        atoms1 = read(args.filelist[1], index=0)
    else:
        raise SystemExit('Error: Too many arguments')

    if not path.splitext(args.out)[1]:
        raise SystemExit('Error: Output file must have an extension')

    if args.mass:
        atoms0.get_weights = atoms0.get_masses
        atoms1.get_weights = atoms1.get_masses

    if args.mirror:
        atoms1.mirror_x()

    if args.reorder:
        assignments = assign_atoms(
            atoms0,
            atoms1,
            **vars(args)
        )
        write(args.out, atoms0, comment='Reference')
        minrmsd = inf
        for i, assignment in enumerate(assignments, start=1):
            atoms2 = atoms1[assignment.order]
            rmsd = atoms2.align_to(atoms0)
            minrmsd = min(minrmsd, rmsd)
            comment = 'Map={0} RMSD={1:.4f}'.format(i, rmsd)
            write(args.out, atoms2, append=True, comment=comment)
        if not args.stats:
            print('Optimized RMSD = {:.4f}'.format(minrmsd))
    else:
        rmsd = atoms1.align_to(atoms0)
        print('RMSD {:.4f}'.format(rmsd))
        write(args.out, atoms0, comment='Reference')
        comment='RMSD={:.4f}'.format(rmsd)
        write(args.out, atoms1, append=True, comment=comment)

if __name__ == '__main__':
    molalign()
