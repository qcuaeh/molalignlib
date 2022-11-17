from ase import io
from argparse import ArgumentParser
from molalignlib import Alignment

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
    #parser.add_argument('-live', action='store_true')
    #parser.add_argument('-stdin', action='store_true')
    #parser.add_argument('-out', type=str, choices=['xyz', 'mol2'], default='xyz')
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

    # Create an alignment object with atoms0 as reference
    alignment0 = Alignment(
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

    if args.sort:
        # Sort atoms1 to minimize the RMSD respect to atoms0
        maplist, mapcount, mindist = alignment0.sorted(atoms1)
        # Align atoms1 to atoms0 for each calculated mapping and write coordinates to file
        for i, mapping in enumerate(maplist, start=1):
            alignedatoms1 = alignment0.aligned(atoms1[mapping])
            io.write('aligned_{}.xyz'.format(i), atoms0)
            io.write('aligned_{}.xyz'.format(i), alignedatoms1, append=True)
    else:
        alignedatoms1 = alignment0.aligned(atoms1)
        #rmsd = np.sqrt((weights0*((alignedatoms1 - atoms0)**2).sum(axis=0)).mean())
        #print('RMSD: {} (only alignment performed)'.format(rmsd))
        io.write('aligned_0.xyz', atoms0)
        io.write('aligned_0.xyz', alignedatoms1, append=True)
