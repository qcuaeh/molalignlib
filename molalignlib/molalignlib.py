import numpy as np
from ase import Atoms

try:
    from .f2py_molalignlib import molalignlib
except ModuleNotFoundError:
    from os import path, symlink
    from subprocess import Popen, PIPE, STDOUT
    import molalignutil
    moduledir = path.dirname(path.abspath(__file__))
    scriptdir = path.dirname(path.abspath(molalignutil.__file__))
    symlink(path.join(scriptdir, 'gnu.env'), path.join(moduledir, 'build.env'))
    command = [path.join(scriptdir, 'compile.sh'), '-all', '-pic', '.', '.']
    with Popen(command, cwd=moduledir, stdout=PIPE, stderr=STDOUT, bufsize=0) as p:
        for line in p.stdout:
            print(line.decode('utf-8').rstrip())
    command = [path.join(scriptdir, 'link-pylib.sh'), '.', 'f2py_molalignlib']
    with Popen(command, cwd=moduledir, stdout=PIPE, stderr=STDOUT, bufsize=0) as p:
        for line in p.stdout:
            print(line.decode('utf-8').rstrip())

class Alignment:
    def __init__(self, rmsd, atoms):
        self.rmsd = rmsd
        self.atoms = atoms

class Assignment:
    def __init__(self, map, count, rmsd):
        self.map = map
        self.count = count
        self.rmsd = rmsd

class Align(Atoms):
    def __init__(
        self,
        atoms,
        weights = None,
    ):
        if isinstance(atoms, Atoms):
            self.__dict__.update(atoms.__dict__)
        else:
            raise TypeError('An Atoms object was expected as argument')
        if weights is None:
            self.weights = np.ones(len(atoms), dtype=float)/len(atoms)
        elif type(weights) is np.ndarray and weights.dtype is float:
            self.weights = weights/sum(weights)
        else:
            raise TypeError('"weights" must be a real numpy array')
    def __call__(self, other):
        if not isinstance(other, Atoms):
            raise TypeError('An Atoms object was expected as first argument')
        if len(other) != len(self):
            raise ValueError('First argument does no have the right length')
        znums0 = self.get_atomic_numbers()
        znums1 = other.get_atomic_numbers()
        types0 = np.ones(len(self), dtype=int)
        types1 = np.ones(len(other), dtype=int)
        coords0 = self.get_positions().transpose() # Convert to column-major order
        coords1 = other.get_positions().transpose() # Convert to column-major order
        rmsd, coords1 = molalignlib.align_atoms(znums0, znums1, types0, types1, \
            coords0, coords1, self.weights)
        # Convert back to row-major order
        atoms1 = Atoms(numbers=znums1, positions=coords1.transpose())
        return Alignment(rmsd, atoms1)

class Assign(Atoms):
    def __init__(
        self,
        atoms,
        weights = None,
        biasing = False,
        bias_tol = 0.2,
        bias_scale = 1.e3,
        iteration = False,
        testing = False,
        records = 1,
        count = 10,
        trials = None,
    ):
        if isinstance(atoms, Atoms):
            self.__dict__.update(atoms.__dict__)
        else:
            raise TypeError('An Atoms object was expected as argument')
        if type(biasing) is not bool:
            raise TypeError('"biasing" must be a boolean')
        if type(iteration) is not bool:
            raise TypeError('"iteration" must be a boolean')
        if type(testing) is not bool:
            raise TypeError('"testing" must be a boolean')
        if type(records) is not int:
            raise TypeError('"records" must be an integer')
        if type(count) is not int:
            raise TypeError('"count" must be an integer')
        if type(bias_tol) is not float:
            raise TypeError('"bias_tol" must be a real')
        if type(bias_scale) is not float:
            raise TypeError('"bias_scale" must be a real')
        if trials is None:
            molalignlib.set_abort_flag(False)
        elif type(trials) is int:
            molalignlib.set_abort_flag(True)
            molalignlib.set_max_trials(trials)
        else:
            raise TypeError('"trials" must be an integer')
        if weights is None:
            self.weights = np.ones(len(atoms), dtype=float)/len(atoms)
        elif type(weights) is np.ndarray and weights.dtype is float:
            self.weights = weights/sum(weights)
        else:
            raise TypeError('"weights" must be a real numpy array')
        self.records = records
        molalignlib.set_bias_flag(biasing)
        molalignlib.set_bias_scale(bias_scale)
        molalignlib.set_bias_tol(bias_tol)
        molalignlib.set_conv_flag(iteration)
        molalignlib.set_test_flag(testing)
        molalignlib.set_max_count(count)
    def __call__(self, other):
        if not isinstance(other, Atoms):
            raise TypeError('An Atoms object was expected as argument')
        if len(other) != len(self):
            raise ValueError('Argument does no have the right length')
        znums0 = self.get_atomic_numbers()
        znums1 = other.get_atomic_numbers()
        types0 = np.ones(len(self), dtype=int)
        types1 = np.ones(len(other), dtype=int)
        coords0 = self.get_positions().transpose() # Convert to column-major order
        coords1 = other.get_positions().transpose() # Convert to column-major order
        nmap, maplist, countlist, rmsdlist = molalignlib.assign_atoms(znums0, znums1, \
            types0, types1, coords0, coords1, self.weights, self.records)
        maplist = maplist - 1 # Convert to 0-based indexing
        return [Assignment(maplist[:, i], countlist[i], rmsdlist[i]) for i in range(nmap)]
