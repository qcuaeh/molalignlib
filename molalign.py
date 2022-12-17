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
from argparse import ArgumentParser
from molalignlib import Assignment, Alignment

def molalign():

    # Parse arguments
    parser = ArgumentParser()
    parser.add_argument('-sort', action='store_true')
    parser.add_argument('-fast', action='store_true')
    parser.add_argument('-test', action='store_true')
    parser.add_argument('-mass', action='store_true')
    parser.add_argument('-enan', action='store_true')
    parser.add_argument('-trials', type=int)
    parser.add_argument('-tol', type=float, default=0.35)
    parser.add_argument('-count', type=int, default=10)
    parser.add_argument('-rec', type=int, default=1)
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

    if args.enan:
        atoms1.positions[:, 0] = -atoms1.positions[:, 0]

    if args.fast:
        biasing = True
        iteration = True
    else:
        biasing = False
        iteration = False

    if args.test:
        reproducible = True
    else:
        reproducible = False

    if args.sort:
        assignments = Assignment(
            atoms0,
            atoms1,
            biasing = biasing,
            iteration = iteration,
            reproducible = reproducible,
            records = args.rec,
            biastol = args.tol,
            maxcount = args.count,
            maxtrials = args.trials,
            mw = args.mass,
        )
        # Align atoms1 to atoms0 for each calculated mapping and write coordinates to file
        for i, mapping in enumerate(assignments, start=1):
            alignment = Alignment(atoms0, atoms1[mapping], mw=args.mass)
            if not args.test:
                io.write('aligned_{}.{ext}'.format(i, ext=args.out), atoms0)
                io.write('aligned_{}.{ext}'.format(i, ext=args.out), alignment.align(atoms1[mapping]), append=True)
    else:
        alignment = Alignment(atoms0, atoms1, mw=args.mass)
        print('RMSD = {:.4f} (only alignment performed)'.format(alignment.rmsd))
        if not args.test:
            io.write('aligned.xyz', atoms0)
            io.write('aligned.xyz', alignment.align(atoms1), append=True)

if __name__ == '__main__':
    molalign()
