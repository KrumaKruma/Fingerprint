! Copyright (C) 2020 Jonas A. Finkler
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Implementation of the algorithm given in section 4 of:
!   Carpaneto, Giorgio, Silvano Martello, and Paolo Toth.
!   "Algorithms and codes for the assignment problem."
!   Annals of operations research 13.1 (1988): 191-223.

module assignmentProblem
    use precision
    implicit none

    real(dp), parameter :: eps =  epsilon(eps)
    private
    public :: solveAP

contains

    subroutine solveAP(n, A, f, cost)
        integer, intent(in) :: n
        real(dp), intent(in) :: A(n,n)
        integer, intent(out) :: f(n)
        real(dp), intent(out) :: cost

        real(dp) :: u(n), v(n)
        integer :: c(n)
        integer :: i
        integer :: j
        integer :: ff(n)

        call initialize(n, A, f, ff, u, v)
        ! After init sum(u) + sum(v) is a lower bound to the total cost

        do i=1,n
            if (f(i)==0) then
                call path(n, A, ff, u, v, c, i, j)
                call increase(n, f, ff, c, i, j)
            end if
        end do

        cost = sum(u) + sum(v)

    end subroutine

    subroutine initialize(n, A, f, ff, u, v)
        integer, intent(in) :: n
        real(dp), intent(in) :: A(n,n)
        integer, intent(out) :: f(n)
        integer, intent(out) :: ff(n)
        real(dp), intent(out) :: u(n)
        real(dp), intent(out) :: v(n)
        integer :: p(n)
        integer :: i, j, k, r
        real(dp) :: tmp
        logical :: assign

        f(:) = 0
        ff(:) = 0

        ! phase 1
        do j=1,n
            ! find r: a(r,j) = min(a(:,j))
            tmp = huge(A)
            do k=1,n
                if (A(k,j) < tmp) then
                    if (A(k,j) == tmp) then ! break ties such that f(r) == 0
                        if (f(r)/=0 .and. f(k)==0) then
                            r = k
                            tmp = A(k,j)
                        end if
                    else
                        r = k
                        tmp = A(k,j)
                    end if
                end if
            end do

            v(j) = A(r,j)

            if (f(r) == 0) then
                ff(j) = r
                f(r) = j
                u(r) = 0._dp
                p(r) = j + 1
            end if
        end do

        ! phase 2
        do i=1,n
            if (f(i) == 0) then
                ! find j: a(i,j) = min(a(i,:) - v(:))
                tmp = huge(A)
                j = 1 ! initialize j, in case first iteration goes to tie braking (can happen if 'huge' is in A)
                do k=1,n
                    if (A(i,k)-v(k) <= tmp) then
                        if (A(i,k)-v(k) == tmp) then ! break ties such that ff(r) == 0
                            if (ff(j)/=0 .and. ff(k)==0) then
                                j = k
                                tmp = A(i,k) - v(k)
                            end if
                        else
                            j = k
                            tmp = A(i,k) - v(k)
                        end if
                    end if
                end do

                u(i) = A(i,j) - v(j)

                assign = ff(j) == 0
                do while((.not. assign) .and. (j<=n))
                    if (abs(A(i,j) - u(i) - v(j)) < eps) then ! todo: should be doable without use of eps?
                        r = ff(j)
                        k = p(r)
                        do while((.not.assign) .and. (k<=n))
                            if (ff(k) == 0 .and. abs(a(r,k) - u(r) - v(k)) < eps) then
                                assign = .true.
                                f(r) = k
                                ff(k) = r
                            else
                                k = k + 1
                            end if
                        end do
                        p(r) = k + 1
                    end if
                    if (.not. assign) then
                        j = j + 1
                    end if
                end do

                if (assign) then
                    f(i) = j
                    ff(j) = i
                    p(i) = j + 1
                end if

            end if
        end do

    end subroutine

    subroutine path(n, A, ff, u, v, c, i, j)
        integer, intent(in) :: n
        real(dp), intent(in) :: A(n,n) ! cost matrix
        integer, intent(inout) :: ff(n) ! ff(i) = row assigned to column i
        real(dp), intent(inout) :: u(n) ! substracted from rows
        real(dp), intent(inout) :: v(n) ! substracted from columns
        integer, intent(out) :: c(n) ! row preceding column j in current path
        integer, intent(in) :: i
        integer, intent(out) :: j
        integer :: lr(n), nlr ! labeled rows
        integer :: uc(n), nuc ! unlabeled columns
        real(dp) :: p(n) ! pi in paper = min(a(:,j)-u(:)-v(j)) for i in LR, i/=ff(i)
        integer :: x

        integer :: r
        logical :: nz
        real(dp) :: d
        integer :: ucj
        real(dp) :: tmp

        lr(1) = i
        nlr = 1
        do x=1,n
            uc(x) = x
        end do
        nuc = n
        p(:) = huge(p)

        do
            r = lr(nlr) ! last element in current path
            ! go through unassigned columns and determine the best predecessor (out of labeled rows)
            ! we update the list with the current row r (newest one)
            do x=1,nuc
                tmp = a(r,uc(x)) - u(r) - v(uc(x))
                if (tmp < p(uc(x))) then ! < p(uc(x)) + eps ?
                    p(uc(x)) = tmp ! cost
                    c(uc(x)) = r ! predecessor
                end if
            end do

            ! check if no p(j) is zero for j in uc
            nz = .true.
            d = huge(d)
            do x=1,nuc
                if (abs(p(uc(x))) < eps) then ! p(uc(x)) is zer
                    ! if (p(uc(x)) /= 0._dp)print*, 'a', p(uc(x))
                    j = uc(x) ! take j to be the first col in UC for which p(j)==0
                    ucj = x ! remember x for removing j from uc
                    nz = .false.
                    exit
                end if
                ! d = min(d, p(uc(x)))
                if (p(uc(x)) < d) then ! d is the min of p(i) for i in uc
                    d = p(uc(x))
                    ! in case no row with zero cost will be found, we use d
                    ! to create a zero entry in p
                    ! we remember where so that we know what our j will be
                    ucj = x
                    j = uc(x)
                end if
            end do

            ! no row with zero cost found
            ! we go and create one by substracting d (the smallest possible number that can create on)
            if (nz) then
                do x=1,nlr
                    u(lr(x)) = u(lr(x)) + d
                end do
                do x=1,n
                    if (abs(p(x)) < eps) then
                        ! if (p(x) /= 0._dp)print*, 'aa', p(x)
                        v(x) = v(x) - d
                    else
                        p(x) = p(x) - d
                    end if
                end do
            end if

            ! j is first col in UC for which p(j) == 0
            if (ff(j) > 0) then
                ! add j to labeled rows
                nlr = nlr + 1
                !
                lr(nlr) = ff(j)
                ! remove j from uc
                uc(ucj) = uc(nuc)
                nuc = nuc - 1
            else
                return
            end if
        end do

    end subroutine

    subroutine increase(n, f, ff, c, i, j)
        integer, intent(in) :: n
        integer, intent(inout) :: f(n)
        integer, intent(inout) :: ff(n)
        integer, intent(in) :: c(n)
        integer, intent(in) :: i
        integer, intent(inout) :: j

        integer :: l, k

        l = i+1
        do while(l/=i)
            l = c(j)
            ff(j) = l
            k = f(l)
            f(l) = j
            j = k
        end do

    end subroutine

end module assignmentProblem
