! MolAlignLib
! Copyright (C) 2022 José M. Vásquez

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module chemdata
! Purpose: Definition of physical constants

use settings

! elsym: Element symbols
! covradii: Atomic covalent radii (Angstrom)
! vdwradii: Atomic van der Waals radii (Angstrom)
! stdmasses: Standard atomic masses
! valencies: Element valenc

implicit none

integer, parameter :: nelem = 103

character(2), parameter :: elsym(nelem) = [ &
'H ',                                                                                                 'He', &
'Li', 'Be',                                                             'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
'Na', 'Mg',                                                             'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', &
'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', &
'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', &
'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',       &
                  'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', &
'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'        &
]

! Lit.: R.T. Sanderson, Inorganic Chemistry, Reinhold 1967
real(wp), parameter :: covradii(nelem) = [ &
0.31,                                                                                                 0.28, &
1.28, 0.96,                                                             0.84, 0.76, 0.71, 0.66, 0.57, 0.58, &
1.66, 1.41,                                                             1.21, 1.11, 1.07, 1.05, 1.02, 1.06, &
2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16, &
2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.40, &
2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 1.75,       &
                  1.87, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48, 1.40, 1.50, 1.50, &
2.60, 2.21, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.69, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00        &
]

! Lit.: A. Bondi, J. Phys. Chem. 68, 441 (1964)         
real(wp), parameter :: vdwradii(nelem) = [ &
1.20,                                                                                                 1.40, &
1.82, 2.00,                                                             2.00, 1.70, 1.55, 1.52, 1.47, 1.54, &
2.27, 1.73,                                                             2.00, 2.10, 1.80, 1.80, 1.75, 1.88, &
2.75, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 1.63, 1.40, 2.00, 1.87, 2.00, 1.85, 1.90, 1.85, 2.02, &
2.50, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 1.63, 1.72, 1.58, 1.93, 2.17, 2.00, 2.06, 1.98, 2.16, &
3.00, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50,       &
                  2.50, 2.50, 2.50, 2.50, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, &
3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00        &
]

! Lit.: CRC Handbook of Chemistry and Physics, 1989
real(wp), parameter :: stdmasses(nelem) = [ &
1.0,                                                                                                                     4.0, &
6.9,     9.0,                                                                        10.8,  12.0,  14.0,  16.0,  19.0,  20.2, &
23.0,   24.3,                                                                        27.0,  28.1,  31.0,  32.1,  35.5,  39.9, &
39.1,   40.1,  45.0,  47.9,  50.9,  52.0,  54.9,  55.8,  58.9,  58.7,  63.5,  65.4,  69.7,  72.6,  74.9,  79.0,  79.9,  83.8, &
85.5,   87.6,  88.9,  91.2,  92.9,  95.9,  98.0, 101.1, 102.9, 106.4, 107.9, 112.4, 114.8, 118.7, 121.8, 127.6, 126.9, 131.3, &
132.9, 137.3, 138.9, 140.1, 140.9, 144.2, 145.0, 150.4, 152.0, 157.2, 158.9, 162.5, 164.9, 167.3, 168.9, 173.0, 175.0,        &
                     178.5, 180.9, 183.8, 186.2, 190.2, 192.2, 195.1, 197.0, 200.6, 204.4, 207.2, 209.0, 209.0, 210.0, 222.0, &
223.0, 226.0, 227.0, 232.0, 231.0, 238.0, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0, 258.0, 259.0, 262.0         &
]

! s-block: min(v, 2 - v), p-block: min(v, 8 - v), d,f-block: 5
integer, parameter :: valencies(nelem) = [ &
1,                                                 0, &
1, 2,                               3, 4, 3, 2, 1, 0, &
1, 2,                               3, 4, 3, 2, 1, 0, &
1, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 4, 3, 2, 1, 0, &
1, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 4, 3, 2, 1, 0, &
1, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,    &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 4, 3, 2, 1, 0, &
1, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5     &
]

end module
