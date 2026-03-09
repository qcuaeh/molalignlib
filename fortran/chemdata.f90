! MolAlignLib
! Copyright (C) 2025 José M. Vásquez

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
use parameters
use str_utils
implicit none

integer(ik), parameter :: symlen = 3 ! Symbol length
integer(ik), parameter :: num_elems = 103

! Element symbols
character(symlen), parameter :: atomic_symbols(num_elems) = [ &
'h ', &
'he', &
'li', &
'be', &
'b ', &
'c ', &
'n ', &
'o ', &
'f ', &
'ne', &
'na', &
'mg', &
'al', &
'si', &
'p ', &
's ', &
'cl', &
'ar', &
'k ', &
'ca', &
'sc', &
'ti', &
'v ', &
'cr', &
'mn', &
'fe', &
'co', &
'ni', &
'cu', &
'zn', &
'ga', &
'ge', &
'as', &
'se', &
'br', &
'kr', &
'rb', &
'sr', &
'y ', &
'zr', &
'nb', &
'mo', &
'tc', &
'ru', &
'rh', &
'pd', &
'ag', &
'cd', &
'in', &
'sn', &
'sb', &
'te', &
'i ', &
'xe', &
'cs', &
'ba', &
'la', &
'ce', &
'pr', &
'nd', &
'pm', &
'sm', &
'eu', &
'gd', &
'tb', &
'dy', &
'ho', &
'er', &
'tm', &
'yb', &
'lu', &
'hf', &
'ta', &
'w ', &
're', &
'os', &
'ir', &
'pt', &
'au', &
'hg', &
'tl', &
'pb', &
'bi', &
'po', &
'at', &
'rn', &
'fr', &
'ra', &
'ac', &
'th', &
'pa', &
'u ', &
'np', &
'pu', &
'am', &
'cm', &
'bk', &
'cf', &
'es', &
'fm', &
'md', &
'no', &
'lr'  &
]

! Standard atomic masses
! Source: mendeleev Python library
real(rk), target :: atomic_masses(num_elems) = [ &
1.0, &  ! H
4.0, &  ! He
6.9, &  ! Li
9.0, &  ! Be
10.8, &  ! B
12.0, &  ! C
14.0, &  ! N
16.0, &  ! O
19.0, &  ! F
20.2, &  ! Ne
23.0, &  ! Na
24.3, &  ! Mg
27.0, &  ! Al
28.1, &  ! Si
31.0, &  ! P
32.1, &  ! S
35.5, &  ! Cl
39.9, &  ! Ar
39.1, &  ! K
40.1, &  ! Ca
45.0, &  ! Sc
47.9, &  ! Ti
50.9, &  ! V
52.0, &  ! Cr
54.9, &  ! Mn
55.8, &  ! Fe
58.9, &  ! Co
58.7, &  ! Ni
63.5, &  ! Cu
65.4, &  ! Zn
69.7, &  ! Ga
72.6, &  ! Ge
74.9, &  ! As
79.0, &  ! Se
79.9, &  ! Br
83.8, &  ! Kr
85.5, &  ! Rb
87.6, &  ! Sr
88.9, &  ! Y
91.2, &  ! Zr
92.9, &  ! Nb
96.0, &  ! Mo
97.9, &  ! Tc
101.1, &  ! Ru
102.9, &  ! Rh
106.4, &  ! Pd
107.9, &  ! Ag
112.4, &  ! Cd
114.8, &  ! In
118.7, &  ! Sn
121.8, &  ! Sb
127.6, &  ! Te
126.9, &  ! I
131.3, &  ! Xe
132.9, &  ! Cs
137.3, &  ! Ba
138.9, &  ! La
140.1, &  ! Ce
140.9, &  ! Pr
144.2, &  ! Nd
144.9, &  ! Pm
150.4, &  ! Sm
152.0, &  ! Eu
157.2, &  ! Gd
158.9, &  ! Tb
162.5, &  ! Dy
164.9, &  ! Ho
167.3, &  ! Er
168.9, &  ! Tm
173.0, &  ! Yb
175.0, &  ! Lu
178.5, &  ! Hf
180.9, &  ! Ta
183.8, &  ! W
186.2, &  ! Re
190.2, &  ! Os
192.2, &  ! Ir
195.1, &  ! Pt
197.0, &  ! Au
200.6, &  ! Hg
204.4, &  ! Tl
207.2, &  ! Pb
209.0, &  ! Bi
209.0, &  ! Po
210.0, &  ! At
222.0, &  ! Rn
223.0, &  ! Fr
226.0, &  ! Ra
227.0, &  ! Ac
232.0, &  ! Th
231.0, &  ! Pa
238.0, &  ! U
237.0, &  ! Np
244.0, &  ! Pu
243.0, &  ! Am
247.0, &  ! Cm
247.0, &  ! Bk
251.0, &  ! Cf
252.0, &  ! Es
257.0, &  ! Fm
258.0, &  ! Md
259.0, &  ! No
262.0  &  ! Lr
]

! Atomic covalent radii (Angstrom)
! Source: mendeleev Python library
real(rk), parameter :: covalent_radii(num_elems) = [ &
0.32, &  ! H
0.46, &  ! He
1.33, &  ! Li
1.02, &  ! Be
0.85, &  ! B
0.75, &  ! C
0.71, &  ! N
0.63, &  ! O
0.64, &  ! F
0.67, &  ! Ne
1.55, &  ! Na
1.39, &  ! Mg
1.26, &  ! Al
1.16, &  ! Si
1.11, &  ! P
1.03, &  ! S
0.99, &  ! Cl
0.96, &  ! Ar
1.96, &  ! K
1.71, &  ! Ca
1.48, &  ! Sc
1.36, &  ! Ti
1.34, &  ! V
1.22, &  ! Cr
1.19, &  ! Mn
1.16, &  ! Fe
1.11, &  ! Co
1.10, &  ! Ni
1.12, &  ! Cu
1.18, &  ! Zn
1.24, &  ! Ga
1.21, &  ! Ge
1.21, &  ! As
1.16, &  ! Se
1.14, &  ! Br
1.17, &  ! Kr
2.10, &  ! Rb
1.85, &  ! Sr
1.63, &  ! Y
1.54, &  ! Zr
1.47, &  ! Nb
1.38, &  ! Mo
1.28, &  ! Tc
1.25, &  ! Ru
1.25, &  ! Rh
1.20, &  ! Pd
1.28, &  ! Ag
1.36, &  ! Cd
1.42, &  ! In
1.40, &  ! Sn
1.40, &  ! Sb
1.36, &  ! Te
1.33, &  ! I
1.31, &  ! Xe
2.32, &  ! Cs
1.96, &  ! Ba
1.80, &  ! La
1.63, &  ! Ce
1.76, &  ! Pr
1.74, &  ! Nd
1.73, &  ! Pm
1.72, &  ! Sm
1.68, &  ! Eu
1.69, &  ! Gd
1.68, &  ! Tb
1.67, &  ! Dy
1.66, &  ! Ho
1.65, &  ! Er
1.64, &  ! Tm
1.70, &  ! Yb
1.62, &  ! Lu
1.52, &  ! Hf
1.46, &  ! Ta
1.37, &  ! W
1.31, &  ! Re
1.29, &  ! Os
1.22, &  ! Ir
1.23, &  ! Pt
1.24, &  ! Au
1.33, &  ! Hg
1.44, &  ! Tl
1.44, &  ! Pb
1.51, &  ! Bi
1.45, &  ! Po
1.47, &  ! At
1.42, &  ! Rn
2.23, &  ! Fr
2.01, &  ! Ra
1.86, &  ! Ac
1.75, &  ! Th
1.69, &  ! Pa
1.70, &  ! U
1.71, &  ! Np
1.72, &  ! Pu
1.66, &  ! Am
1.66, &  ! Cm
1.68, &  ! Bk
1.68, &  ! Cf
1.65, &  ! Es
1.67, &  ! Fm
1.73, &  ! Md
1.76, &  ! No
1.61  &  ! Lr
]

! Atomic Van der Waals radii (Angstrom)
! Source: mendeleev Python library
real(rk), parameter :: vdw_radii(num_elems) = [ &
1.10, &  ! H
1.40, &  ! He
1.82, &  ! Li
1.53, &  ! Be
1.92, &  ! B
1.70, &  ! C
1.55, &  ! N
1.52, &  ! O
1.47, &  ! F
1.54, &  ! Ne
2.27, &  ! Na
1.73, &  ! Mg
1.84, &  ! Al
2.10, &  ! Si
1.80, &  ! P
1.80, &  ! S
1.75, &  ! Cl
1.88, &  ! Ar
2.75, &  ! K
2.31, &  ! Ca
2.15, &  ! Sc
2.11, &  ! Ti
2.07, &  ! V
2.06, &  ! Cr
2.05, &  ! Mn
2.04, &  ! Fe
2.00, &  ! Co
1.97, &  ! Ni
1.96, &  ! Cu
2.01, &  ! Zn
1.87, &  ! Ga
2.11, &  ! Ge
1.85, &  ! As
1.90, &  ! Se
1.85, &  ! Br
2.02, &  ! Kr
3.03, &  ! Rb
2.49, &  ! Sr
2.32, &  ! Y
2.23, &  ! Zr
2.18, &  ! Nb
2.17, &  ! Mo
2.16, &  ! Tc
2.13, &  ! Ru
2.10, &  ! Rh
2.10, &  ! Pd
2.11, &  ! Ag
2.18, &  ! Cd
1.93, &  ! In
2.17, &  ! Sn
2.06, &  ! Sb
2.06, &  ! Te
1.98, &  ! I
2.16, &  ! Xe
3.43, &  ! Cs
2.68, &  ! Ba
2.43, &  ! La
2.42, &  ! Ce
2.40, &  ! Pr
2.39, &  ! Nd
2.38, &  ! Pm
2.36, &  ! Sm
2.35, &  ! Eu
2.34, &  ! Gd
2.33, &  ! Tb
2.31, &  ! Dy
2.30, &  ! Ho
2.29, &  ! Er
2.27, &  ! Tm
2.26, &  ! Yb
2.24, &  ! Lu
2.23, &  ! Hf
2.22, &  ! Ta
2.18, &  ! W
2.16, &  ! Re
2.16, &  ! Os
2.13, &  ! Ir
2.13, &  ! Pt
2.14, &  ! Au
2.23, &  ! Hg
1.96, &  ! Tl
2.02, &  ! Pb
2.07, &  ! Bi
1.97, &  ! Po
2.02, &  ! At
2.20, &  ! Rn
3.48, &  ! Fr
2.83, &  ! Ra
2.47, &  ! Ac
2.45, &  ! Th
2.43, &  ! Pa
2.41, &  ! U
2.39, &  ! Np
2.43, &  ! Pu
2.44, &  ! Am
2.45, &  ! Cm
2.44, &  ! Bk
2.45, &  ! Cf
2.45, &  ! Es
2.45, &  ! Fm
2.46, &  ! Md
2.46, &  ! No
2.46  &  ! Lr
]

end module
