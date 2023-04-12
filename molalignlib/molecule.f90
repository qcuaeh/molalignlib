!|-module molecule.f95 :
!|  -fortran module with definitions of the derived type molecule .
!|  -intended to declare basic manging methods in pseudo object-oriented style .
!|  -molecule type definitions :
!|    -parameters/fields/vars :
!|      -mtype :-molecule type according the molecule hierarchy (see below) ;
!|      -natom :-number of atoms/particles ;
!|      -nbond :-number of bonds ;
!|      -title :-text line ;
!|      -coord :-Cartesian coordinates matrix in Angstroms ;
!|      -label :-atom labels array ;;;
!|  -molecule hierarchy .
!|    -cluster (1) :
!|      -homoatomic cluster .
!|      -molecule type parameter :-mtype_cluster=1 ;
!|      -stored information :
!|        -natom, coord, title ;;
!|    -xyz :
!|      -heteroatomic molecule .
!|      - ;;
!|  - ;
module molecule

use kinds

implicit none

! parameters
integer, parameter, private :: max_natom=1000, max_ncoord=10, max_namelen=5
integer, parameter, private :: mtype_cluster=1, &   ! homoatomic cluster 
                               mtype_xyz=2, &   ! heteroatomic molecule
                               mtype_mol2=3, &   ! heteroatomic molecule with mol2 information

type molecule
    integer :: mtype, natom, nbond
    character(len=ll) :: title
    real(wp), dimension(:,:), pointer, allocatable :: coord
    character(len=max_namelen), dimension(:), pointer, allocatable :: label
    integer, dimension(:,:), pointer, allocatable :: bonds
end type molecule

contains

end module 
