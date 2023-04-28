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
!|        -Cartesian coordinates matrix in Angstroms .
!|        -matrix size :
!|          -(3 x natom) ;;
!|      -label :
!|        -atom labels array .
!|        -vector size :
!|          -natom ;;
!|      -atom :
!|        -list of atom structures (see atom type definition above) .
!|        -this provisional structure duplicates information .
!|        -vector size :
!|          -natom ;;;;
!|  -molecule hierarchy .
!|    -undefined (0) :
!|      -non specified molecule .
!|      -this is how molecules are initialized .
!|      -molecule type parameter :
!|        -mtype_cluster=0 ;
!|      -defined information :
!|        -mtype ;;
!|    -cluster (1) :
!|      -homoatomic cluster .
!|      -molecule type parameter :
!|        -mtype_cluster=1 ;
!|      -defined information :
!|        -mtype, natom, coord, title ;;
!|    -xyz :
!|      -heteroatomic molecule read from a .xyz file .
!|      -molecule type parameter :
!|        -mtype_xyz=2 ;
!|      -defined information :
!|        -mtype, natom, coord, title .
!|        -label ;;
!|    -mol2 :
!|      -heteroatomic molecule read from a .mol2 file .
!|      -molecule type parameter :
!|        -mtype_xyz=3 ;
!|      -defined information :
!|        -mtype, natom, coord, title .
!|        -label .
!|        -nbond, bonds ;;;
!|  -module definitions (types, variables and parameters) :
!|    - ;
module molecule

use kinds

implicit none

! parameters
integer, parameter, private :: max_natom=1000, max_ncoord=10, max_namelen=5
integer, parameter, private :: mtype_undef=0, &   ! not specified mol type
                               mtype_cluster=1, &   ! homoatomic cluster 
                               mtype_xyz=2, &   ! heteroatomic molecule
                               mtype_mol2=3,   ! molecule with mol2 info

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

!|  -module contains section :
!|    -molecule methods :

!|-subroutine new_molecule (mol, prev_mtype) :
!|  -initializes a new molecule instance with unspecified type .
!|  -arguments :
!|    -mol :
!|      -molecule-type data structure .
!|      -to be initialized as undefined molecule type .
!|      -an already used molecule could be reinitialized ;
!|    -prev_mtype :
!|      -(optional*) prior molecule type of mol ... ;;;;

subroutine new_molecule (mol, prev_mtype)
   type(molecule), intent(inout) :: mol
   integer, intent(in) :: prev_mtype   ! could be set as optional

!*** code for handling deallocation of previously stored data ... ***

   mol%mtype=0
   mol%natom=0
   mol%nbond=0
   mol%title=""
   nullify(coord, label, bonds, atom)
end subroutine new_molecule

subroutine read_xyz (unit, mol)
   integer, intent(in) :: unit
   type(molecule), intent(inout) :: mol

   if (mol%mtype /= 0) then
      call new_molecule (mol)
   end if   

! ****

end subroutine read_xyz

end module

!|  - ;


