subroutine ralign(mol0, mol1, label0, label1, bonds0, bonds1, mapcount, atomaplist, rotquatlist)
! Purpose: Find best product to reactant map

    use iso_fortran_env, only: output_unit
    use iso_fortran_env, only: error_unit

    use globals
    use random
    use decoding
    use utilities
    use alignment
    use sorting
    use remapping
    use translation
    use assortment
    use chemistry
    use messages

    implicit none

    integer, intent(out) :: mapcount
    integer, dimension(:, :), allocatable, intent(out) :: atomaplist
    real, dimension(:, :), allocatable, intent(out) :: rotquatlist
    real, dimension(:, :), allocatable, intent(inout) :: mol0, mol1
    character(32), dimension(:), allocatable, intent(inout) :: label0, label1
    integer, dimension(:, :), allocatable, intent(inout) :: bonds0, bonds1

    integer i
    integer natom0, natom1, nbond0, nbond1
    integer nblock0, nblock1, nequiv0, nequiv1
    integer, dimension(:), allocatable :: seed
    integer, dimension(:), allocatable :: znum0, znum1
    integer, dimension(:), allocatable :: ncoord0, ncoord1, nadjeq0, nadjeq1
    integer, dimension(:), allocatable :: blocksize0, blocksize1, blocktype0, blocktype1
    integer, dimension(:), allocatable :: equivtype0, equivtype1, equivsize0, equivsize1
    integer, dimension(:), allocatable :: identitymap, order0, order1, reorder0
    logical, dimension(:, :), allocatable :: adjmat0, adjmat1
    integer, dimension(:, :), allocatable :: adjlist0, adjlist1
    integer, dimension(:, :), allocatable :: adjeqsize0, adjeqsize1
    real, dimension(:), allocatable :: weights0
    real, dimension(:, :), allocatable :: bias
    real, dimension(nelem) :: atomweight, adjrad

    if (maxcount < 0) then
        call error('Argument is missing!')
    end if

    select case (weighting)
    case ('none')
        atomweight = [(1., i=1, nelem)]
    case ('mass')
        atomweight = stdmatom
    case ('valency')
        atomweight = valency
    case default
        call error('Invalid weighting option: '//trim(weighting))
    end select

    ! Set adjacency radii and check them

    adjrad = covrad + 0.25*(vdwrad - covrad)

    if (any(adjrad < covrad) .or. any(adjrad > vdwrad)) then
        call error('There are unphysical atomic radii!')
    end if

    ! Get number of atoms

    natom0 = size(mol0, dim=1)
    natom1 = size(mol1, dim=1)

    ! Get number of bonds

    nbond0 = size(bonds0, dim=1)
    nbond1 = size(bonds1, dim=1)

    ! Abort if number of atoms differ

    if (natom0 /= natom1) then
        call error('The molecules are not isomers!')
    end if

    ! Allocate arrays

    allocate(weights0(natom0))
    allocate(identitymap(natom0))
    allocate(znum0(natom0), znum1(natom1))
    allocate(ncoord0(natom0), ncoord1(natom1))
    allocate(adjmat0(natom0, natom0), adjmat1(natom1, natom1))
    allocate(adjlist0(maxcoord, natom0), adjlist1(maxcoord, natom1))

    ! Get atomic numbers

    do i = 1, natom0
        label0(i) = upper(label0(i))
        znum0(i) = znum(label0(i))
    end do

    do i = 1, natom1
        label1(i) = upper(label1(i))
        znum1(i) = znum(label1(i))
    end do

    ! Print stats and return if maxcount is zero

    if (maxcount == 0) then
        if (any(label0 /= label1)) then
            call error('Can not align because the atomic labels do not match!')
        end if
        identitymap = [(i, i=1, natom0)]
        weights0 = atomweight(znum0)/sum(atomweight(znum0))
        call print_header()
        call print_stats(0, 0, 0, 0., 0., 0., leastsquaredist(natom0, weights0, mol0, mol1, identitymap))
        call print_footer(.false., .false., 0, 0)
        mapcount = 1
        atomaplist(:, 1) = identitymap
        return
    end if

    allocate(order0(natom0), order1(natom1))
    allocate(adjeqsize0(maxcoord, natom0), adjeqsize1(maxcoord, natom1))
    allocate(blocktype0(natom0), blocktype1(natom1), blocksize0(natom0), blocksize1(natom1))
    allocate(equivtype0(natom0), equivtype1(natom1), equivsize0(natom0), equivsize1(natom1))
    allocate(nadjeq0(natom0), nadjeq1(natom1))
    allocate(atomaplist(natom0, maxrecord))
    allocate(rotquatlist(3, maxrecord))
    allocate(bias(natom0, natom1))

    ! Group atoms by label

    call grouplabels(natom0, label0, nblock0, blocksize0, blocktype0)
    call grouplabels(natom1, label1, nblock1, blocksize1, blocktype1)

    ! Sort equivalent atoms in contiguous blocks

    order0 = sortorder(equivtype0, natom0)
    order1 = sortorder(equivtype1, natom1)

    blocktype0 = blocktype0(order0)
    blocktype1 = blocktype1(order1)

    label0 = label0(order0)
    label1 = label1(order1)

    mol0 = mol0(:, order0)
    mol1 = mol1(:, order1)

    znum0 = znum0(order0)
    znum1 = znum1(order1)

    reorder0 = inversemap(order0)

    ! Abort if there are incompatible atomic symbols

    if (any(label0 /= label1)) then
        call error('The molecules do not have the same number of atoms!')
    end if

    ! Abort if there are incompatible atomic labels

    if (any(label0 /= label1)) then
        call error('There are incompatible atomic labels!')
    end if

    ! Calculate normalized weights

    weights0 = atomweight(znum0)/sum(atomweight(znum0))

    ! Translate molecules to their centroid

    call translate(natom0, centroid(natom0, weights0, mol0), mol0)
    call translate(natom1, centroid(natom1, weights0, mol1), mol1)

    ! Set bias for non equivalent atoms 

    call setadjbias(natom0, nblock0, blocksize0, blocktype0, mol0, ncoord0, adjlist0, mol1, &
        ncoord1, adjlist1, weights0, bias)

    ! Initialize random number generator

    call init_random_seed(seed)
    call random_seed(put=seed)

    ! Remap atoms to minimize distance and difference

    call remapatoms(natom0, nblock0, blocksize0, blocktype0, mol0, ncoord0, adjlist0, adjmat0, &
        nequiv0, equivsize0, equivtype0, nadjeq0, adjeqsize0, mol1, ncoord1, adjlist1, adjmat1, &
        nequiv1, equivsize1, equivtype1, nadjeq1, adjeqsize1, weights0, mapcount, atomaplist, bias)

    ! Calculate optimal rotation matrices

    do i = 1, mapcount
        rotquatlist(:, i) = leastrotquat(natom0, weights0, mol0, mol1, atomaplist(:, i))
    end do

    ! Restore original atom order

    label0 = label0(reorder0)
    label1 = label1(reorder0)

    mol0 = mol0(:, reorder0)
    mol1 = mol1(:, reorder0)

    do i = 1, mapcount
        atomaplist(:, i) = atomaplist(reorder0, i)
    end do

end subroutine
