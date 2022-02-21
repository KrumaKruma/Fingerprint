! Copyright (C) 2021 Marco Krummenacher, Behnam Parsaeifard
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



!!! TODO THERE IS STILL A GOTO STATEMENT...

!------------------------------------------------------------------------------
! OM Fingerprint Package
!------------------------------------------------------------------------------
!
! MODULE: Fingerprint
!
!> Marco Krummenacher
!> University of Basel
!
! DESCRIPTION:
!> In this module a maximum volume simplex compression is applied to the
!> overlap matrix fingerprint. The compressed fingerprint as well as its
!> derivative can be calculated.
!------------------------------------------------------------------------------




module simplex
  IMPLICIT NONE

CONTAINS

  ! Main subroutines
  !---------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  !> This routine calculates the simplex corner and the compressed fingerprint
  !> vector.
  !> @param[in] len_fp length of the fingerprint vector
  !> @param[in] natx number of atoms (environments)
  !> @param[in] nconf number of configurations
  !> @param[in] nsimplex number of corners for the compression
  !> @param[in] fpall fingerprint vecotr of all configurations and environments
  !> @param[out] fp compressed fingerprint
  !> @param[out] fpcorner fingerprint of the simplex corners
  !---------------------------------------------------------------------------

  subroutine SimplexSparse(len_fp, natx, nconf, nsimplex, fpall, fp, fpcorner)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: len_fp
    INTEGER, INTENT(IN) :: natx
    INTEGER, INTENT(IN) :: nconf
    INTEGER, INTENT(IN) :: nsimplex
    REAL(8), DIMENSION(len_fp, natx, nconf), INTENT(IN) :: fpall
    REAL(8), DIMENSION(nsimplex, natx, nconf), INTENT(OUT) :: fp
    REAL(8), DIMENSION(len_fp, 0:nsimplex), INTENT(OUT) :: fpcorner

    !dummy
    INTEGER, DIMENSION(0:nsimplex) :: indat
    INTEGER, DIMENSION(0:nsimplex) :: indconf
    REAL(8), ALLOCATABLE :: xcart(:,:)
    REAL(8), ALLOCATABLE :: distsimplex(:,:)
    REAL(8), ALLOCATABLE :: xcartw(:,:)
    REAL(8), ALLOCATABLE :: distsimplexw(:,:)
    LOGICAL, ALLOCATABLE :: mask(:,:)
    INTEGER :: i, j, iat, iiat, jjat, iconf, iiconf, jjconf
    INTEGER :: isim, nsim, msim
    !REAL(8) :: fpdf
    REAL(8) :: dmax
    REAL(8) :: distmax

    REAL(8) :: local_dmax
    INTEGER :: local_iiat
    INTEGER :: local_iiconf

    REAL(8) :: t1
    REAL(8) :: t2
    REAL(8) :: time
    INTEGER :: count1
    INTEGER :: count2
    INTEGER :: count_rate
    INTEGER :: count_max






    ALLOCATE(xcart(nsimplex,nsimplex))
    ALLOCATE(distsimplex(0:nsimplex,0:nsimplex))
    ALLOCATE(mask(natx,nconf))

    DO j = 0, nsimplex
      DO i = 0, nsimplex
        distsimplex(i,j)=0.d0
      ENDDO
    ENDDO

    DO j = 1, nsimplex
      DO i = 1, nsimplex
        xcart(i,j)=0.d0
      ENDDO
    ENDDO

    DO iconf = 1, nconf
      DO iat = 1, natx
        mask(iat,iconf)=.true.
      ENDDO
    ENDDO
    CALL cpu_time(t1)
    CALL system_clock(count1,count_rate,count_max)
    CALL largest_pair_dist(len_fp,natx,nconf, fpall,iiat,iiconf,jjat,jjconf,distmax)
    !WRITE(*,*) iiat, iiconf, jjat, jjconf, distmax
    CALL cpu_time(t2)
    CALL system_clock(count2,count_rate,count_max)
    WRITE(*,*) 'CPU time 1:    ', t2-t1, " s"
    time=(count2-count1)/float(count_rate)
    WRITE(*,*) 'Elapsed time1 ', time ," s"



    indat(0)=iiat
    indconf(0)=iiconf
    indat(1)=jjat
    indconf(1)=jjconf
    mask(iiat,iiconf)=.false.
    mask(jjat,jjconf)=.false.
    xcart(1,1)=distmax
    distsimplex(0,0)=0.d0
    distsimplex(1,1)=0.d0
    distsimplex(0,1)=distmax
    distsimplex(1,0)=distmax


    CALL cpu_time(t1)
    CALL system_clock(count1,count_rate,count_max)
    DO msim = 2, nsimplex

      dmax=0.d0
      iiconf=-1
      iiat=-1

      !$omp parallel private(iconf, iat, isim,&
      !$omp     &local_dmax,local_iiat,local_iiconf) firstprivate(distsimplex, xcart)

      local_dmax=0.d0
      local_iiconf=-1
      local_iiat=-1

      !$omp do schedule(guided) collapse(2)

      DO iconf = 1, nconf
        DO iat = 1, natx
          IF(mask(iat,iconf)) THEN

            distsimplex(msim,msim)=0.d0
            DO isim = 0, msim-1
              distsimplex(msim,isim)=fpdf(len_fp,fpall(1,iat,iconf),fpall(1,indat(isim),indconf(isim)))
              distsimplex(isim,msim)=distsimplex(msim,isim)
            ENDDO

            CALL add_point_simplex(nsimplex,msim,distsimplex,xcart)
            IF(xcart(msim,msim).gt.local_dmax) THEN
              local_iiconf=iconf
              local_iiat=iat
              local_dmax=xcart(msim,msim)
            ENDIF

          ENDIF
        ENDDO
      ENDDO
      !$omp end do
      !$omp critical
      IF(local_dmax .ge. dmax) THEN
        dmax = local_dmax
        iiat = local_iiat
        iiconf = local_iiconf
      ENDIF
      !$omp end critical
      !$omp end parallel



      CALL cpu_time(t2)
      CALL system_clock(count2,count_rate,count_max)
      WRITE(*,*) 'CPU time 2:    ', t2-t1, " s"
      time=(count2-count1)/float(count_rate)
      WRITE(*,*) 'Elapsed time 2', time ," s"


      IF(dmax.lt.1.d-6) THEN
        nsim=msim-1
        goto 2000
      ENDIF
      CALL cpu_time(t1)
      CALL system_clock(count1,count_rate,count_max)
      distsimplex(msim,msim)=0.d0
      DO isim = 0, msim-1
        distsimplex(msim,isim)=fpdf(len_fp,fpall(1,iiat,iiconf),fpall(1,indat(isim),indconf(isim)))
        distsimplex(isim,msim)=distsimplex(msim,isim)
      ENDDO
      CALL cpu_time(t2)
      CALL system_clock(count2,count_rate,count_max)
      WRITE(*,*) 'CPU time 3:    ', t2-t1, " s"


      CALL cpu_time(t1)
      CALL system_clock(count1,count_rate,count_max)


      CALL add_point_simplex(nsimplex,msim,distsimplex,xcart)


      CALL cpu_time(t2)
      CALL system_clock(count2,count_rate,count_max)
      WRITE(*,*) 'CPU time 4:    ', t2-t1, " s"


      indat(msim)=iiat
      indconf(msim)=iiconf
      mask(iiat,iiconf)=.false.

    ENDDO
    nsim=nsimplex

    2000 continue

    IF(nsim .ne. nsimplex )  WRITE(*,*) "nsim.ne.nsimplex, ", nsim, nsimplex
    DO i = 0, nsim
      DO j = 1, len_fp
        fpcorner(j, i) = fpall(j, indat(i), indconf(i))
      ENDDO
    ENDDO


    ALLOCATE(xcartw(nsim+1,nsim+1))
    ALLOCATE(distsimplexw(0:nsim+1,0:nsim+1))
    DO j = 1, nsim
      DO i = 1, nsim
        xcartw(i,j)=xcart(i,j)
      ENDDO
    ENDDO
    DO j = 0, nsim
      DO i = 0, nsim
        distsimplexw(i,j)=distsimplex(i,j)
      ENDDO
    ENDDO

    !contracted fingerprint
    CALL cpu_time(t1)
    CALL system_clock(count1,count_rate,count_max)

    !$omp parallel private(iconf, iat, isim) firstprivate(distsimplexw, xcartw)
    !$omp do schedule(guided) collapse(2)
    DO iconf = 1, nconf
      DO iat = 1, natx
        distsimplexw(nsim+1,nsim+1)=0.d0
        DO isim = 0, nsim
          distsimplexw(nsim+1,isim)=fpdf(len_fp,fpall(1,iat,iconf),fpall(1,indat(isim), indconf(isim)) )
          distsimplexw(isim,nsim+1)=distsimplexw(nsim+1,isim)
        ENDDO
        CALL add_point_simplex(nsim+1,nsim+1,distsimplexw,xcartw)
        DO isim = 1, nsim
          fp(isim,iat,iconf) = xcartw(isim,nsim+1)
          !write(122,*) fp(isim,iat,iconf)
        ENDDO
      ENDDO
    ENDDO
    !$omp end do
    !$omp end parallel


    CALL cpu_time(t2)
    CALL system_clock(count2,count_rate,count_max)
    WRITE(*,*) 'CPU time 5:    ', t2-t1, " s"


    DEALLOCATE(xcart)
    DEALLOCATE(distsimplex)
    DEALLOCATE(xcartw)
    DEALLOCATE(distsimplexw)
    DEALLOCATE(mask)
  end subroutine

  ! Main subroutines
  !---------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  !> This routine calculates the simplex corner and the compressed fingerprint
  !> vector.
  !> @param[in] len_fp length of the fingerprint vector
  !> @param[in] natx number of atoms (environments)
  !> @param[in] nconf number of configurations
  !> @param[in] nsimplex number of corners for the compression
  !> @param[in] fpall fingerprint vecotor of all configurations and environments
  !> @param[in] fpallder derivative fingerprint vecotor of all configurations and environments
  !> @param[out] fp compressed fingerprint
  !> @param[out] fpgrad compressed fingerprint derivatives
  !> @param[out] fpcorner fingerprint of the simplex corners
  !> @param[in] ws if true, write simpex information to files
  !> @param[in] rs if true, read simplex information from files
  !---------------------------------------------------------------------------



subroutine SimplexSparse_derivatve(len_fp, natx, nconf, nsimplex, fpall,&
   fpallder, fp, fpgrad, fpcorner, ws, rs)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: len_fp
  INTEGER, INTENT(IN) :: natx
  INTEGER, INTENT(IN) :: nconf
  INTEGER, INTENT(IN) :: nsimplex
  REAL(8), DIMENSION(len_fp, natx, nconf), INTENT(IN) :: fpall
  REAL(8), DIMENSION(len_fp, 3, natx, natx, nconf), INTENT(IN) :: fpallder
  REAL(8), DIMENSION(nsimplex, natx, nconf), INTENT(OUT) :: fp
  REAL(8), DIMENSION(nsimplex, 3, natx, natx, nconf), INTENT(OUT) :: fpgrad
  REAL(8), DIMENSION(len_fp, 0:nsimplex), INTENT(OUT) :: fpcorner

  !dummy
  INTEGER, DIMENSION(0:nsimplex) :: indat
  INTEGER, DIMENSION(0:nsimplex) :: indconf
  REAL(8), ALLOCATABLE :: xcart(:,:)
  REAL(8), ALLOCATABLE :: distsimplex(:,:)
  REAL(8), ALLOCATABLE :: xcartw(:,:)
  REAL(8), ALLOCATABLE :: distsimplexw(:,:)
  LOGICAL, ALLOCATABLE :: mask(:,:)
  INTEGER :: i, j, l, iat, jat, iiat, jjat, iconf, iiconf, jjconf
  INTEGER :: isim, jsim, nsim, msim
  !REAL(8) :: fpdf
  REAL(8) :: dmax
  REAL(8) :: distmax
  LOGICAL :: ws
  LOGICAL :: rs

  REAL(8) :: local_dmax
  INTEGER :: local_iiconf
  INTEGER :: local_iiat


  ALLOCATE(xcart(nsimplex,nsimplex))
  ALLOCATE(distsimplex(0:nsimplex,0:nsimplex))
  ALLOCATE(mask(natx,nconf))
  IF(rs) THEN
    nsim = nsimplex
    ALLOCATE(xcartw(nsim+1,nsim+1))
    ALLOCATE(distsimplexw(0:nsim+1,0:nsim+1))

    OPEN(UNIT=12, FILE="../data/fp_corner.dat")
    DO isim = 0,nsimplex
      DO l = 1, len_fp
        READ(12,*) fpcorner(l,isim)
      ENDDO
    ENDDO
    CLOSE(12)

    OPEN(UNIT=11, FILE="../data/distsimplexw.dat")
    DO j = 0, nsimplex
      DO i = 0, nsimplex
        READ(11,*) distsimplexw(i,j)
      ENDDO
    ENDDO
    CLOSE(11)

    OPEN(UNIT=56, FILE="../data/xcartw.dat")
    DO j = 1, nsimplex
      DO i = 1, nsimplex
        READ(56,*) xcartw(i,j)
      ENDDO
    ENDDO
    CLOSE(56)

  ELSE

    DO j = 0, nsimplex
      DO i = 0, nsimplex
        distsimplex(i,j)=0.d0
      ENDDO
    ENDDO

    DO j = 1, nsimplex
      DO i = 1, nsimplex
        xcart(i,j)=0.d0
      ENDDO
    ENDDO

    DO iconf = 1, nconf
      DO iat = 1, natx
        mask(iat,iconf)=.true.
      ENDDO
    ENDDO

    CALL largest_pair_dist(len_fp,natx,nconf, fpall,iiat,iiconf,jjat,jjconf,distmax)
    !WRITE(*,*) iiat, iiconf, jjat, jjconf, distmax

    indat(0)=iiat
    indconf(0)=iiconf
    indat(1)=jjat
    indconf(1)=jjconf
    mask(iiat,iiconf)=.false.
    mask(jjat,jjconf)=.false.
    xcart(1,1)=distmax
    distsimplex(0,0)=0.d0
    distsimplex(1,1)=0.d0
    distsimplex(0,1)=distmax
    distsimplex(1,0)=distmax

    DO msim = 2, nsimplex

      dmax=0.d0
      iiconf=-1
      iiat=-1

      !$omp parallel private(iconf, iat,&
      !$omp     &local_dmax,local_iiat,local_iiconf) firstprivate(distsimplex, xcart)

      local_dmax=0.d0
      local_iiconf=-1
      local_iiat=-1

      !$omp do schedule(guided) collapse(2)
      DO iconf = 1, nconf
        DO iat = 1, natx
          IF(mask(iat,iconf)) THEN

            distsimplex(msim,msim)=0.d0
            DO isim = 0, msim-1
              distsimplex(msim,isim)=fpdf(len_fp,fpall(1,iat,iconf),fpall(1,indat(isim),indconf(isim)))
              distsimplex(isim,msim)=distsimplex(msim,isim)
            ENDDO

            CALL add_point_simplex(nsimplex,msim,distsimplex,xcart)
            IF(xcart(msim,msim).gt.dmax) THEN
              iiconf=iconf
              iiat=iat
              dmax=xcart(msim,msim)
            ENDIF

          ENDIF
        ENDDO
      ENDDO
      !$omp end do
      !$omp critical
      IF(local_dmax .ge. dmax) THEN
        dmax = local_dmax
        iiat = local_iiat
        iiconf = local_iiconf
      ENDIF
      !$omp end critical
      !$omp end parallel





      IF(dmax.lt.1.d-6) THEN
        nsim=msim-1
        goto 2000
      ENDIF

      distsimplex(msim,msim)=0.d0
      DO isim = 0, msim-1
        distsimplex(msim,isim)=fpdf(len_fp,fpall(1,iiat,iiconf),fpall(1,indat(isim),indconf(isim)))
        distsimplex(isim,msim)=distsimplex(msim,isim)
      ENDDO

      CALL add_point_simplex(nsimplex,msim,distsimplex,xcart)

      indat(msim)=iiat
      indconf(msim)=iiconf
      mask(iiat,iiconf)=.false.

    ENDDO
    nsim=nsimplex

    2000 continue

    IF(nsim .ne. nsimplex )  WRITE(*,*) "nsim.ne.nsimplex, ", nsim, nsimplex
    DO i = 0, nsim
      DO j = 1, len_fp
        fpcorner(j, i) = fpall(j, indat(i), indconf(i))
      ENDDO
    ENDDO


    ALLOCATE(xcartw(nsim+1,nsim+1))
    ALLOCATE(distsimplexw(0:nsim+1,0:nsim+1))
    DO j = 1, nsim
      DO i = 1, nsim
        xcartw(i,j)=xcart(i,j)
      ENDDO
    ENDDO
    DO j = 0, nsim
      DO i = 0, nsim
        distsimplexw(i,j)=distsimplex(i,j)
      ENDDO
    ENDDO

    IF(ws) THEN

      OPEN(UNIT=10, FILE="../data/xcartw.dat")
      DO j = 1, nsim
        DO i = 1, nsim
          WRITE(10,*) xcartw(i,j)
        ENDDO
      ENDDO
      CLOSE(10)

      OPEN(UNIT=11, FILE="../data/distsimplexw.dat")
      DO j = 0, nsim
        DO i = 0, nsim
          WRITE(11,*) distsimplexw(i,j)
        ENDDO
      ENDDO
      CLOSE(11)

      OPEN(UNIT=12, FILE="../data/fp_corner.dat")
      DO isim = 0,nsim
        DO l = 1, len_fp
          WRITE(12,*) fpcorner(l,isim)
        ENDDO
      ENDDO
      CLOSE(12)

    ENDIF
  ENDIF



  !contracted fingerprint
  !$omp parallel private(iconf, iat, isim) firstprivate(distsimplexw, xcartw)
  !$omp do schedule(guided) collapse(2)
  DO iconf = 1, nconf
    DO iat = 1, natx
      distsimplexw(nsim+1,nsim+1)=0.d0
      DO isim = 0, nsim
        distsimplexw(nsim+1,isim)=fpdf(len_fp,fpall(1,iat,iconf),fpcorner(1,isim) )
        distsimplexw(isim,nsim+1)=distsimplexw(nsim+1,isim)
      ENDDO
      CALL add_point_simplex(nsim+1,nsim+1,distsimplexw,xcartw)
      DO isim = 1, nsim
        fp(isim,iat,iconf) = xcartw(isim,nsim+1)
      ENDDO
    ENDDO
  ENDDO
  !$omp end do
  !$omp end parallel


  !contracted derivatives
  !derivative of fingerprint of iat with respect to cartesian coordinates of jat


  !$omp parallel private(iconf, iat, jat, isim, jsim, l) firstprivate(xcartw, xcart, fpcorner)
  !$omp do schedule(guided) collapse(3)
  DO iconf = 1, nconf
    DO jat = 1, natx
      DO iat = 1, natx
        DO isim = 1, nsim
          DO l = 1, 3
            fpgrad(isim,l,jat,iat,iconf) = (1.d0/xcartw(isim,isim))*&
                sum( (fpcorner(:,isim)-fpcorner(:,0))*fpallder(:,l,jat,iat,iconf) )
            DO jsim=1, isim-1
              fpgrad(isim,l,jat,iat,iconf) = fpgrad(isim,l,jat,iat,iconf) &
                 -(xcart(jsim,isim)/xcartw(isim,isim))*fpgrad(jsim,l,jat,iat,iconf)
            ENDDO !jsim
          ENDDO !l
        ENDDO !isim
      ENDDO !iat
    ENDDO !jat
  ENDDO !iconf
  !$omp end do
  !$omp end parallel

  DEALLOCATE(xcart)
  DEALLOCATE(distsimplex)
  DEALLOCATE(xcartw)
  DEALLOCATE(distsimplexw)
  DEALLOCATE(mask)
end subroutine


subroutine add_point_simplex(nsimplex,nsim,dist,xcart)
  IMPLICIT NONE
  INTEGER :: nsim
  INTEGER :: nsimplex
  REAL(8), DIMENSION(nsimplex+1,nsimplex+1) :: dist
  REAL(8), DIMENSION(nsimplex+1,nsimplex+1) :: d2
  REAL(8), DIMENSION(nsimplex,nsimplex) :: xcart

  INTEGER :: i, j, k

  IF(nsim.le.1) STOP 'add_point'

  DO j = 1, nsimplex+1
     DO i = 1, nsimplex+1
       d2(i,j)=dist(i,j)**2
     ENDDO
  ENDDO

  j = nsim
  DO i = 1, j-1
    xcart(i,j)=d2(1,j+1) + d2(1,i+1) - d2(i+1,j+1)
    DO k = i-1, 1, -1
      xcart(i,j)=xcart(i,j)-2.d0*xcart(k,i)*xcart(k,j)
    ENDDO
    xcart(i,j)=xcart(i,j)/(2.d0*xcart(i,i))
        !write(*,*) i,j,xcart(i,j)
  ENDDO
  xcart(j,j)= d2(1,j+1)
  DO k = 1, j-1
    xcart(j,j)= xcart(j,j) - xcart(k,j)*xcart(k,j)
  ENDDO
  xcart(j,j)=sqrt(max(0.d0,xcart(j,j)))

end subroutine add_point_simplex

! Main subroutines
!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> Find largest fingerprint pair distance
!> @param[in] len_fp length of the fingerprint vector
!> @param[in] natx number of atoms (environments)
!> @param[in] nconf number of configurations
!> @param[in] fpall fingerprint vecotor of all configurations and environments
!> @param[out] iiat ith environment which has largest dist to jth environment
!> @param[out] iiconf ith configuration which has largest dist to jth configuration
!> @param[out] jjat jth environment which has largest dist to ith environment
!> @param[out] jjconf jth configuration which has largest dist to ith configuration
!> @param[out] distmax maximal fingerprint distance

!---------------------------------------------------------------------------

subroutine largest_pair_dist(len_fp,natx,nconf,fpall,iiat,iiconf,jjat,jjconf,distmax)
!  Find points with the largest distances
  IMPLICIT NONE
  INTEGER :: len_fp
  INTEGER :: natx
  INTEGER :: nconf
  INTEGER :: iiat
  INTEGER :: jjat
  INTEGER :: iiconf
  INTEGER :: jjconf
  REAL(8) :: distmax
  REAL(8) :: fpd
  REAL(8), DIMENSION(len_fp, natx, nconf) :: fpall

  INTEGER :: iat, jat, iconf, jconf


  REAL(8) :: local_distmax
  INTEGER :: local_iiat
  INTEGER :: local_iiconf
  INTEGER :: local_jjat
  INTEGER :: local_jjconf

  distmax=0.d0

  !$omp parallel private(iconf, iat, jconf, jat, fpd,&
    !$omp     &local_distmax,local_iiat,local_jjat,local_iiconf,local_jjconf)
    local_distmax=-1.d99
    local_iiat=0
    local_iiconf=0
    local_jjat=0
    local_jjconf=0
    !$omp do schedule(guided) collapse(2)
    DO iconf = 1, nconf
      DO iat = 1, natx
        DO jconf = iconf+1, nconf
          DO jat = 1, natx
              fpd=fpdf(len_fp,fpall(1,iat,iconf),fpall(1,jat,jconf))
              IF (fpd .gt. local_distmax) THEN
                local_distmax= fpd
                local_iiat=iat
                local_iiconf=iconf
                local_jjat=jat
                local_jjconf=jconf
              ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !$omp end do
    !$omp critical
    if(local_distmax.gt.distmax) THEN
      distmax=local_distmax
      iiat=local_iiat
      jjat=local_jjat
      iiconf=local_iiconf
      jjconf=local_jjconf
    endif
    !$omp end critical
    !$omp end parallel

end subroutine largest_pair_dist






!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> fingerprint distance of two fingerprint vectors
!> @param[in] len length of the fingerprint vector
!> @param[in] fp1 fingerprint vecotor 1
!> @param[in] fp2 fingerprint vector 2
!> @param[out] fpd fingerprint distance
!---------------------------------------------------------------------------


function fpdf(len,fp1,fp2) result(fpd)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: len
 INTEGER :: l
 REAL(8) :: fpd
 REAL(8) :: tt
 REAL(8), DIMENSION(len), INTENT(IN) :: fp1
 REAL(8), DIMENSION(len), INTENT(IN) :: fp2

 tt=0.d0
 DO l = 1, len
   tt=tt+(fp1(l)-fp2(l))**2
 ENDDO
 fpd=sqrt(tt)

end
!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> citation of the maximum volume simplex method
!---------------------------------------------------------------------------
subroutine cite_simplex()
  WRITE(*,*) "--> B. Parsaeifard et al. Machine Learning: Science and Technology 2, 015018 (2021)"
end subroutine

end module simplex
