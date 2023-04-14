!|-module molecule.f90 :| (provisional comments)
!|  -fortran module with definitions of the derived type molecule .
!|  -intended to declare basic managing methods in pseudo object-oriented style .
!|  -atom type definitions :
!|    -parameters :
!|      -id .
!|      -pos .
!|      -atype .
!|      -label .
!|      -coord .
!|      -elsym .
!|      -covrad .
!|      -vdwrad .
!|      -stdmass .
!|      -valence .
!|      -nadj .
!|      -adjlist .
!|      - ;;
!|  -molecule type definitions :
!|    -parameters/fields/vars :
!|      -mtype :
!|        -molecule type according the molecule hierarchy (see below) ;
!|      -natom :
!|        -number of atoms/particles ;
!|      -nbond :
!|        -number of bonds ;
!|      -title :
!|        -text line ;
!|      -coord :
!|        -Cartesian coordinates matrix in Angstroms ;
!|      -label :
!|        -atom labels array ;
!|      -atom :
!|        -list of atom structures (see atom type definition above) .
!|        -this provisional structure duplicates information ;;;
!|  -molecule hierarchy .
!|    -cluster (1) :
!|      -homoatomic cluster .
!|      -molecule type parameter :
!|        -mtype_cluster=1 ;
!|      -stored information :
!|        -mtype, natom, coord, title ;;
!|    -xyz :
!|      -heteroatomic molecule read from a .xyz file .
!|      -molecule type parameter :
!|        -mtype_xyz=2 ;
!|      -stored information :
!|        -mtype, natom, coord, title .
!|        -label ;;
!|    -mol2 :
!|      -heteroatomic molecule read from a .mol2 file .
!|      -molecule type parameter :
!|        -mtype_xyz=3 ;
!|      -stored information :
!|        -mtype, natom, coord, title .
!|        -label .
!|        -nbond, bonds ;;;
!|  -molecule methods :
!|    -constructor new_molecule :
!|      -creates a new instance with a specified type ;;
!|  - ;
module molecule

use kinds

implicit none

! parameters
integer, parameter, private :: max_natom=1000, max_ncoord=10, max_namelen=5
integer, parameter, private :: mtype_cluster=1, &   ! homoatomic cluster 
                               mtype_xyz=2, &   ! heteroatomic molecule
                               mtype_mol2=3, &   ! molecule with mol2 info

type atom
    integer :: id, pos, atype
    character(len=wl) :: label
    real(wp), dimension(3) :: coord
    character(2) :: elsym
    real(wp) :: covrad, vdwrad, stdmass
    integer :: valence, nadj
    integer, dimension(max_ncoord) :: adjlist
end type atom

type molecule
    integer :: mtype, natom, nbond
    character(len=ll) :: title
    real(wp), dimension(:,:), pointer, allocatable :: coord
    character(len=wl), dimension(:), pointer, allocatable :: label
    integer, dimension(:,:), pointer, allocatable :: bonds
    type(atom), dimension(:), pointer, allocatable :: atom
end type molecule

contains

subroutine new_molecule (mol, mtype)
    type(molecule), intent(inout) :: mol
    integer, intent(in) :: mtype
    mol%mtype=mtype
end subroutine new_molecule

subroutine read_xyz (mol, filename)

end subroutine read_xyz

end module


