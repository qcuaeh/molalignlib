#!/usr/bin/python3

from sys import argv
from ase import io
import numpy as np

def spin(atoms):
    x = np.random.uniform(0, 1, size=(3,1))
    #hacemos descomposicion QR para generar una matriz de rotacion aleatoria 
    q, r = np.linalg.qr(x, mode='complete')
    q = q * np.linalg.det(q)
    atoms.positions = np.dot(atoms.positions, q)

def distort(atoms):
    tol = float(argv[2])
    x = np.random.uniform(-tol, tol, size=(len(atoms), 3))
    atoms.positions = atoms.positions + x

def shuffle(atoms):
    indices = np.arange(len(atoms))
    np.random.shuffle(indices)
    atoms.__dict__.update(atoms[indices].__dict__)

#np.random.seed(0)
atoms = io.read(argv[1], index=0)
io.write('test.xyz', atoms)
spin(atoms)
#distort(atoms)
shuffle(atoms)
io.write('test.xyz', atoms, append=True)
