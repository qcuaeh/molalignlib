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
from argparse import ArgumentParser
from ase.io import read, write
from molalignlib import assign_atoms

def molalign():

    parser = ArgumentParser()
    parser.add_argument('filelist', nargs='+')
    parser.add_argument('-sort', action='store_true')
    parser.add_argument('-fast', action='store_true')
    parser.add_argument('-mass', action='store_true')
    parser.add_argument('-enan', action='store_true')
    parser.add_argument('-stats', action='store_true')
    parser.add_argument('-trials', type=int)
    parser.add_argument('-tol', type=float)
    parser.add_argument('-scale', type=float)
    parser.add_argument('-count', type=int)
    parser.add_argument('-rec', type=int, default=1)
    parser.add_argument('-out', type=str, default='aligned.xyz')
#    parser.add_argument('-stdin', type=str)
#    parser.add_argument('-stdout', stype=str)
#    parser.add_argument('-live', action='store_true')
#    parser.add_argument('-test', action='store_true')
    args = parser.parse_args()

    if len(args.filelist) == 1:
        atoms0 = read(args.filelist[0], index=0)
        atoms1 = read(args.filelist[0], index=1)
    elif len(args.filelist) == 2:
        atoms0 = read(args.filelist[0], index=0)
        atoms1 = read(args.filelist[1], index=0)
    else:
        print('Error: Too many files')
        raise SystemExit

    if args.enan:
        atoms1.positions[:, 0] = -atoms1.positions[:, 0]

    if args.fast:
        biasing = True
        iteration = True
    else:
        biasing = False
        iteration = False

    if not path.splitext(args.out)[1]:
        print('Error: Output file must have an extension')
        raise SystemExit

    if args.sort:
        assignments = assign_atoms(
            atoms0,
            atoms1,
            biasing = biasing,
            iteration = iteration,
            massweighted = args.mass,
            stats = args.stats,
            records = args.rec,
            tolerance = args.tol,
            scale = args.scale,
            count = args.count,
            trials = args.trials,
        )
        i = assignments.pop(0)
        atoms2 = atoms1[i.order]
        rmsd = atoms2.align_to(atoms0, massweighted=args.mass)
        print('Optimized RMSD = {:.4f}'.format(rmsd))
        write(args.out, atoms0, comment='Reference')
        write(args.out, atoms2, append=True, comment='RMSD {:.4f}'.format(rmsd))
        for i in assignments:
            atoms2 = atoms1[i.order]
            rmsd = atoms2.align_to(atoms0, massweighted=args.mass)
            write(args.out, atoms2, append=True, comment='RMSD {:.4f}'.format(rmsd))
    else:
        rmsd = atoms1.align_to(atoms0, massweighted=args.mass)
        print('RMSD = {:.4f}'.format(rmsd))
        write(args.out, atoms0, comment='Reference')
        write(args.out, atoms1, append=True, comment='RMSD {:.4f}'.format(rmsd))

if __name__ == '__main__':
    molalign()
