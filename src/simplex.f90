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

module simplex
  IMPLICIT NONE

CONTAINS

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
    INTEGER :: i, j, k, l, iat, jat, iiat, jjat, iconf, jconf, iiconf, jjconf
    INTEGER :: isim, jsim, nsim, msim
    REAL(8) :: fpd
    !REAL(8) :: fpdf
    REAL(8) :: dmax
    REAL(8) :: distmax

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

    !contracted fingerprint
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
          write(122,*) fp(isim,iat,iconf)
        ENDDO
      ENDDO
    ENDDO


    DEALLOCATE(xcart)
    DEALLOCATE(distsimplex)
    DEALLOCATE(xcartw)
    DEALLOCATE(distsimplexw)
    DEALLOCATE(mask)
  end subroutine




subroutine SimplexSparse_derivatve(len_fp, natx, nconf, nsimplex, fpall, fpallder, fp, fpgrad, fpcorner)
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
  INTEGER :: i, j, k, l, iat, jat, iiat, jjat, iconf, jconf, iiconf, jjconf
  INTEGER :: isim, jsim, nsim, msim
  REAL(8) :: fpd
  !REAL(8) :: fpdf
  REAL(8) :: dmax
  REAL(8) :: distmax

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

  !contracted fingerprint
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
        write(123,*) fp(isim,iat,iconf)
      ENDDO
    ENDDO
  ENDDO
  !contracted derivatives
  !derivative of fingerprint of iat with respect to cartesian coordinates of jat
  DO iconf = 1, nconf
    DO jat = 1, natx
      DO iat = 1, natx
        DO isim = 1, nsim
          DO l = 1, 3
            fpgrad(isim,l,jat,iat,iconf) = (1.d0/xcart(isim,isim))*&
                sum( (fpcorner(:,isim)-fpcorner(:,0))*fpallder(:,l,jat,iat,iconf) )
            DO jsim=1, isim-1
              fpgrad(isim,l,jat,iat,iconf) = fpgrad(isim,l,jat,iat,iconf) &
                 -(xcart(jsim,isim)/xcart(isim,isim))*fpgrad(jsim,l,jat,iat,iconf)
            ENDDO !jsim
          ENDDO !l
        ENDDO !isim
      ENDDO !iat
    ENDDO !jat
  ENDDO !iconf

!  DO iconf = 1, nconf
!    DO iat = 1, natx
!      DO isim = 1, nsim
!        fp(isim,iat,iconf) = fp(isim,iat,iconf)/xcart(isim,isim)
!        !if(.true.) write(*,'(3i7,e24.17)') iconf,iat,isim,fp(isim,iat,iconf)
!        DO jat = 1, natx
!          DO l = 1, 3
!            fpgrad(isim,l,jat,iat,iconf) = fpgrad(isim,l,jat,iat,iconf)/xcart(isim,isim)
!          ENDDO
!        ENDDO
!      ENDDO
!    ENDDO
!  ENDDO

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

  INTEGER :: ic, jc, iat, jat, iconf, jconf

  distmax=0.d0
  ic=0
  DO iconf = 1, nconf
    DO iat = 1, natx
      ic=ic+1
      jc=0
      DO jconf = 1, nconf
        DO jat = 1, natx
          jc=jc+1
          IF(jc.lt.ic) THEN
            fpd=fpdf(len_fp,fpall(1,iat,iconf),fpall(1,jat,jconf))
            IF (fpd.gt.distmax) THEN
              distmax=fpd
              iiat=iat
              iiconf=iconf
              jjat=jat
              jjconf=jconf
              !        ii=ic
              !        jj=jc
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO
end subroutine largest_pair_dist

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

subroutine cite_simplex()
  WRITE(*,*) "--> B. Parsaeifard et al. Machine Learning: Science and Technology 2, 015018 (2021)"
end subroutine

end module simplex
