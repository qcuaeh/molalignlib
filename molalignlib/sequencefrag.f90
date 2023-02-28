! purpose: determines the most convenient atom to start runs over the structure
!          of each fragment in a molecule. (returns nfrag, fragrt)
module sequencefrag

use stdio
use sorting
implicit none
private
public runsequence

contains

subroutine runsequence (natom, nadj, adjlist, blksz, blkid, &
                  neqv, eqvsz, nfrag, fragrt, fragid)

   integer, intent(in) :: natom, neqv
   integer, dimension(:), intent(in) :: nadj, blksz, blkid, eqvsz
   integer, dimension(:, :), intent(in) :: adjlist
   integer, intent(out) :: nfrag, fragrt(natom), fragid(natom)

   integer :: h, n, f, firstn, prevn, offset
   integer, dimension(natom) :: eqvid
   integer, dimension(natom) :: order
   integer, dimension(natom) :: fragsize
   integer, dimension(natom,natom) :: fragind, fragblock, fragequiv, fragcoord
   logical :: track(natom), printInfo = .false.
   character(80) :: fmtstr
   character(3) :: numat

   ! assing equivalence group

   offset = 0
   do h = 1, neqv
      eqvid(offset+1:offset+eqvsz(h)) = h
      offset = offset + eqvsz(h)
   end do

   ! initialization

   nfrag = 0
   fragrt(:) = 1

   order = [ (n, n = 1, natom) ]
   track(:) = .false.
   nfrag = 0
   fragid(:) = 0
   fragsize(:) = 0
   fragind(:,:) = 0
   fragblock(:,:) = 0
   fragequiv(:,:) = 0
   fragcoord(:,:) = 0

   ! detect fragments and populate frag arrays
   ! 'order' could be shuffled at this point to randomize the fragments order
   n = 1
   do while ( n <= natom )
      if ( track(order(n)) ) then
         n = n + 1
      else
         nfrag = nfrag + 1
         call recrun (order(n), track)
         n = 1
      end if
   end do

   ! print initial information
   if ( printInfo ) then
      write (numat, '(i0)') natom
      fmtstr = '(i2,x,a,'//numat//'i3)'
      print *, "fragments detected: ", nfrag
      print fmtstr, 0, "sequence ini:   ", (n, n = 1, natom)
      print fmtstr, 0, "fragids ini:    ", fragid(:natom)
      do f = 1, nfrag
         print fmtstr, f, "frag inds ini:  ", fragind(:fragsize(f),f)
         print fmtstr, f, "frag block size:", fragblock(:fragsize(f),f)
         print fmtstr, f, "frag equiv size:", fragequiv(:fragsize(f),f)
         print fmtstr, f, "frag nadj:    ", fragcoord(:fragsize(f),f)
      end do
   end if

   ! sorts sequentially by blksz, eqvsz, and nadj for each fragment
   order = [ (n, n = 1, natom) ]
   do f = 1, nfrag
      firstn = fragsize(f)
      call groupminsort (fragblock(:firstn,f), order(:firstn), firstn)
      fragind(:firstn,f) = fragind(order(:firstn),f)
      fragblock(:firstn,f) = fragblock(order(:firstn),f)
      fragequiv(:firstn,f) = fragequiv(order(:firstn),f)
      fragcoord(:firstn,f) = fragcoord(order(:firstn),f)
      call groupminsort (fragequiv(:firstn,f), order(:firstn), firstn)
      fragind(:firstn,f) = fragind(order(:firstn),f)
      fragblock(:firstn,f) = fragblock(order(:firstn),f)
      fragequiv(:firstn,f) = fragequiv(order(:firstn),f)
      fragcoord(:firstn,f) = fragcoord(order(:firstn),f)
      prevn = firstn
      call groupminsort (fragcoord(:firstn,f), order(:firstn), firstn)
      if ( firstn /= prevn ) order(:prevn) = invsort (order(:prevn), prevn)
      fragind(:prevn,f) = fragind(order(:prevn),f)
      fragblock(:firstn,f) = fragblock(order(:firstn),f)
      fragequiv(:firstn,f) = fragequiv(order(:firstn),f)
      fragcoord(:prevn,f) = fragcoord(order(:prevn),f)
      fragrt(f) = fragind(1,f)
   end do

   ! print final vectors
   if ( printInfo ) then
      write (numat, '(i0)') natom
      fmtstr = '(i2,x,a,'//numat//'i3)'
      do f = 1, nfrag
         print fmtstr, f, "frag inds fin:  ", fragind(:fragsize(f),f)
         print fmtstr, f, "frag block size:", fragblock(:fragsize(f),f)
         print fmtstr, f, "frag equiv size:", fragequiv(:fragsize(f),f)
         print fmtstr, f, "frag nadj:    ", fragcoord(:fragsize(f),f)
         print '(i2,1x,a,1x,i10)', f, "fragment start ind:", fragrt(f)
      end do
   end if

contains

   ! runs recursivelly over the structure and populates arrays
   recursive subroutine recrun (n, track)
      integer, intent(in) :: n
      logical, dimension(:), intent(inout) :: track
      integer i

      if ( track(n) ) return
      track(n) = .true.

      fragsize(nfrag) = fragsize(nfrag) + 1
      fragid(n) = nfrag
      fragind(fragsize(nfrag),nfrag) = n
      fragblock(fragsize(nfrag),nfrag) = blksz(blkid(n))
      fragequiv(fragsize(nfrag),nfrag) = eqvsz(eqvid(n))
      fragcoord(fragsize(nfrag),nfrag) = nadj(n)
      
      do i = 1, nadj(n)
         call recrun (adjList(i, n), track)
      end do
   end subroutine recrun

   ! sorts vec and returns the order and number of elements in the first group
   subroutine groupminsort (vec, order, firstn)
!        integer, intent(in) :: n
      integer, dimension(:), intent(in) :: vec
      integer, dimension(:), intent(out) :: order
      integer, intent(inout) :: firstn
      integer i, n, firstitem

      order = sortorder (vec, firstn)
      n = firstn
      firstitem = vec(order(1))
      do i = 2, n
         if ( vec(order(i)) /= firstitem ) then
            firstn = i - 1
            return
         end if
      end do
   end subroutine groupminsort

   ! inverts the order of a vector
   function invsort (vec, n)
      integer, intent(in) :: vec(:), n
      integer :: i, p, invsort(n)
      p = n
      do i = 1, n
         invsort(p) = vec(i)
         p = p - 1
      end do
   end function invsort

end subroutine runsequence

end module
