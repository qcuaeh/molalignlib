# MolAlign
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

from ase import io
from molalignlib import Alignable, Assignment
from argparse import ArgumentParser

def molalign():

    # Parse arguments
    parser = ArgumentParser()
    parser.add_argument('-sort', action='store_true')
    parser.add_argument('-fast', action='store_true')
    parser.add_argument('-test', action='store_true')
    parser.add_argument('-mass', action='store_true')
    parser.add_argument('-mirror', action='store_true')
    parser.add_argument('-trials', type=int)
    parser.add_argument('-tol', type=float, default=0.35)
    parser.add_argument('-count', type=int, default=10)
    parser.add_argument('-rec', type=int, default=1)
    parser.add_argument('-scale', type=float, default=1.e3)
    parser.add_argument('-out', type=str, default='xyz')
    #parser.add_argument('-live', action='store_true')
    #parser.add_argument('-stdin', action='store_true')
    parser.add_argument('filelist', nargs='+')
    args = parser.parse_args()

    # Read clusters coordinates
    if len(args.filelist) == 1:
        atoms0 = io.read(args.filelist[0], index=0)
        atoms1 = io.read(args.filelist[0], index=1)
    elif len(args.filelist) == 2:
        atoms0 = io.read(args.filelist[0], index=0)
        atoms1 = io.read(args.filelist[1], index=0)
    else:
        print('Error: Too many files')

    if args.mirror:
        atoms1.positions[:, 0] = -atoms1.positions[:, 0]

    if args.fast:
        biasing = True
        iteration = True
    else:
        biasing = False
        iteration = False

    if args.sort:
        assignments = Assignment(
            atoms0,
            atoms1,
            biasing = biasing,
            bias_tol = args.tol,
            bias_scale = args.scale,
            iteration = iteration,
            testing = args.test,
            records = args.rec,
            count = args.count,
            trials = args.trials,
            mass_weighted = args.mass,
        )
        # Align atoms1 to atoms0 for each calculated mapping and write coordinates to file
        for i, order in enumerate(assignments.mapind, start=1):
            alignable1 = Alignable(atoms1[order], mass_weighted=args.mass)
            alignable1.alignto(atoms0)
            io.write('aligned_{}.{ext}'.format(i, ext=args.out), atoms0)
            io.write('aligned_{}.{ext}'.format(i, ext=args.out), alignable1, append=True)
    else:
        alignable1 = Alignable(atoms1, mass_weighted=args.mass)
        alignable1.alignto(atoms0)
        io.write('aligned.xyz', atoms0)
        io.write('aligned.xyz', alignable1, append=True)
        print('RMSD: {:.4f} (only alignment performed)'.format(alignable1.disto(atoms0)))

if __name__ == '__main__':
    molalign()
