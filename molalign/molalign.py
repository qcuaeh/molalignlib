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
from molalignlib import Align, Assign

def main():

    # Parse arguments
    parser = ArgumentParser()
    parser.add_argument('-sort', action='store_true')
    parser.add_argument('-fast', action='store_true')
    parser.add_argument('-test', action='store_true')
    parser.add_argument('-mass', action='store_true')
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

    if args.fast:
        biasing = True
        iteration = True
    else:
        biasing = False
        iteration = False

    if args.mass:
        weights0 = atoms0.get_masses()
    else:
        weights0 = None

    align0 = Align(
        atoms = atoms0,
        weights = weights0,
    )

    if args.sort:
        assign0 = Assign(
            atoms = atoms0,
            weights = weights0,
            testing = args.test,
            biasing = biasing,
            iteration = iteration,
            records = args.rec,
            count = args.count,
            trials = args.trials,
            bias_tol = args.tol,
            bias_scale = args.scale,
        )
        assignments = assign0(atoms1)
        # Align atoms1 to atoms0 for each calculated mapping and write coordinates to file
        for i, a in enumerate(assignments, start=1):
            io.write('aligned_{}.{ext}'.format(i, ext=args.out), atoms0)
            io.write('aligned_{}.{ext}'.format(i, ext=args.out), align0(atoms1[a.map]).atoms, append=True)
    else:
        alignment = align0(atoms1)
        print('RMSD: {:.4f} (only alignment performed)'.format(alignment.rmsd))
        io.write('aligned_0.xyz', atoms0)
        io.write('aligned_0.xyz', alignment.atoms, append=True)
