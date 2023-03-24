subroutine mapatoms(natom, coords0, coords1, nadj0, adjlist0, nadj1, adjlist1, equivmat, biasmat, mapping, dist)

   hoffset = 0
   do h = 1, nblk
      do i = hoffset + 1, hoffset + blklen(h)
         do j = hoffset + 1, hoffset + blklen(h)
            goffset = 0
            do g = 1, nblknei0(i)
               do k = goffset + 1, goffset + blkneilen0(i)
                  do l = goffset + 1, goffset + blkneilen1(j)
                     distmat(k, l) = sum((coords0(:, adjlist0(k)) - coords1(:, adjlist1(l)))**2)
                  end do
               end do
               call minperm(distmat, goffset, blkneilen0(k), mapping, dist)
               goffset = goffset + blkneilen0(k)
            end do
         end do
      end do
      hoffset = hoffset + 1
   end do

end subroutine


