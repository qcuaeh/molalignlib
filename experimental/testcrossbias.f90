subroutine vdwcrossbias(natom, nblk, blksz, blkz, nadj0, adjlist0, &
   nadj1, adjlist1, coords0, coords1, biasmat)
! Purpose: Generate the adjacency matrix

    integer, intent(in) :: natom, nblk
    integer, dimension(:), intent(in) :: nadj0, nadj1
    integer, dimension(:), intent(in) :: blksz, blkz
    integer, dimension(:, :), intent(in) :: adjlist0, adjlist1
    real(wp), dimension(:, :), intent(in) :: coords0, coords1
    real(wp), dimension(:, :), intent(out) :: biasmat

    real(wp) :: atomdist
    integer :: h, i, j, k, n, offset2, offset
    integer, dimension(natom) :: blkid
    integer, dimension(nblk, natom) :: n0, n1
    real(wp), dimension(natom) :: radii
    real(wp), dimension(natom, natom) :: d0, d1

    real(wp), parameter :: tolerance = 1.0_wp

   ! set atoms block indices

   offset = 0
   do h = 1, nblk
      blkid(offset+1:offset+blksz(h)) = h
      radii(offset+1:offset+blksz(h)) = covrad(blkz(h)) + 1.0*(vdwrad(blkz(h)) - covrad(blkz(h)))
      offset = offset + blksz(h)
   end do

    do i = 1, natom
        offset = 0
        do h = 1, nblk
            k = 0
            do j = offset + 1, offset + blksz(h)
                if (j /= i) then
                    k = k + 1
                    atomdist = sqrt(sum((coords0(:, i) - coords0(:, j))**2))
                    d0(offset+k, i) = max(0._wp, atomdist - radii(i) - radii(j))
                end if
            end do
            if (k > 1) then
                call sort(d0(:, i), offset+1, offset+k)
            end if
            offset = offset + blksz(h)
        end do
    end do

    do i = 1, natom
        offset = 0
        do h = 1, nblk
            k = 0
            do j = offset + 1, offset + blksz(h)
                if (j /= i) then
                    k = k + 1
                    atomdist = sqrt(sum((coords1(:, i) - coords1(:, j))**2))
                    d1(offset+k, i) = max(0._wp, atomdist - radii(i) - radii(j))
                end if
            end do
            if (k > 1) then
                call sort(d1(:, i), offset+1, offset+k)
            end if
            offset = offset + blksz(h)
        end do
    end do

    offset = 0
    biasmat(:, :) = 0
    n0(:, :) = 0
    n1(:, :) = 0

    do i = 1, natom
       do k = 1, nadj0(i)
          n0(blkid(adjlist0(k, i)), i) = n0(blkid(adjlist0(k, i)), i) + 1
       end do
    end do

    do i = 1, natom
       do k = 1, nadj1(i)
          n1(blkid(adjlist1(k, i)), i) = n1(blkid(adjlist1(k, i)), i) + 1
       end do
    end do

    do h = 1, nblk
        do i = offset + 1, offset + blksz(h)
            do j = offset + 1, offset + blksz(h)

                offset2 = 0

                do k = 1, nblk
                    if (n0(k, i) > n1(k, j)) then
                        print *, i, j, ':', d1(offset2+1:offset2+n1(k, j), j), '/', d1(offset2+n1(k, j)+1:offset2+n0(k, i), j)
                        if (sum(d1(offset2+1:offset2+n0(k, i), j)) > tolerance) then
                            biasmat(i, j) = bias_scale**2
                        end if
                    else if (n0(k, i) < n1(k, j)) then
                        print *, i, j, ':', d0(offset2+1:offset2+n0(k, i), i), '/', d0(offset2+n0(k, i)+1:offset2+n1(k, j), i)
                        if (sum(d0(offset2+1:offset2+n1(k, j), i)) > tolerance) then
                            biasmat(i, j) = bias_scale**2
                        end if
                    end if
                    offset2 = offset2 + blksz(k)
                end do

            end do
        end do

        offset = offset + blksz(h)

    end do

!    print *, biasmat

end subroutine

subroutine mixedcrossbias(natom, nblk, blksz, nadj0, adjlist0, &
   nadj1, adjlist1, coords0, coords1, biasmat)
! Purpose: Generate the adjacency matrix

    integer, intent(in) :: natom, nblk
    integer, dimension(:), intent(in) :: nadj0, nadj1
    integer, dimension(:), intent(in) :: blksz
    integer, dimension(:, :), intent(in) :: adjlist0, adjlist1
    real(wp), dimension(:, :), intent(in) :: coords0, coords1
    real(wp), dimension(:, :), intent(out) :: biasmat

    integer :: h, i, j, k, n, offset2, offset
    integer, dimension(natom) :: blkid
    integer, dimension(nblk, natom) :: n0, n1
    real(wp), dimension(natom, natom) :: d0, d1

    real(wp), parameter :: tolerance = 1.0_wp

   ! set atoms block indices

   offset = 0
   do h = 1, nblk
      blkid(offset+1:offset+blksz(h)) = h
      offset = offset + blksz(h)
   end do

    do i = 1, natom
        offset = 0
        do h = 1, nblk
            k = 0
            do j = offset + 1, offset + blksz(h)
                if (j /= i) then
                    k = k + 1
                    d0(offset+k, i) = sum((coords0(:, i) - coords0(:, j))**2)
                end if
            end do
            if (k > 1) then
                call sort(d0(:, i), offset+1, offset+k)
            end if
            offset = offset + blksz(h)
        end do
    end do

    do i = 1, natom
        offset = 0
        do h = 1, nblk
            k = 0
            do j = offset + 1, offset + blksz(h)
                if (j /= i) then
                    k = k + 1
                    d1(offset+k, i) = sum((coords1(:, i) - coords1(:, j))**2)
                end if
            end do
            if (k > 1) then
                call sort(d1(:, i), offset+1, offset+k)
            end if
            offset = offset + blksz(h)
        end do
    end do

    offset = 0
    biasmat(:, :) = 0
    n0(:, :) = 0
    n1(:, :) = 0

    do i = 1, natom
       do k = 1, nadj0(i)
          n0(blkid(adjlist0(k, i)), i) = n0(blkid(adjlist0(k, i)), i) + 1
       end do
    end do

    do i = 1, natom
       do k = 1, nadj1(i)
          n1(blkid(adjlist1(k, i)), i) = n1(blkid(adjlist1(k, i)), i) + 1
       end do
    end do

    do h = 1, nblk
        do i = offset + 1, offset + blksz(h)
            do j = offset + 1, offset + blksz(h)

                offset2 = 0

                do k = 1, nblk
                    if (n0(k, i) > n1(k, j)) then
!                        print *, i, j, ':', d1(offset2+1:offset2+n1(k, j), j), '/', d1(offset2+n1(k, j)+1:offset2+n0(k, i), j)
                        biasmat(i, j) = bias_scale*sum(d1(offset2+1:offset2+n0(k, i), j))/n0(k, i)
                    else if (n0(k, i) < n1(k, j)) then
!                        print *, i, j, ':', d0(offset2+1:offset2+n0(k, i), i), '/', d0(offset2+n0(k, i)+1:offset2+n1(k, j), i)
                        biasmat(i, j) = bias_scale*sum(d0(offset2+1:offset2+n1(k, j), i))/n1(k, j)
                    end if
                    offset2 = offset2 + blksz(k)
                end do

            end do
        end do

        offset = offset + blksz(h)

    end do

!    print *, biasmat

end subroutine
