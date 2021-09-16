! Copyright (C) 2021 Marco Krummenacher
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



module fingerprint
  IMPLICIT NONE
CONTAINS

! Main subroutines

  SUBROUTINE only_fingerprint(nat, nat_sphere_max, ns, np, width_cutoff, alat,rxyz, symb, fp)
    INTEGER, INTENT(IN) :: nat      !number of atoms in system
    INTEGER, INTENT(IN) :: ns       !s-orbitals
    INTEGER, INTENT(IN) :: np       !p-orbitals
    INTEGER, INTENT(IN) :: nat_sphere_max  !max. number of atoms in sphere
    REAL(8), INTENT(IN) :: width_cutoff  !cutoff radius
    REAL(8), DIMENSION(3,nat), INTENT(IN) :: rxyz  !coorinates
    REAL(8), DIMENSION(3,3), INTENT(IN) :: alat  !unit cell vectors, if a cluster then alat=0.d0
    REAL(8), DIMENSION(nat_sphere_max*(ns+3*np),nat), INTENT(OUT) :: fp  !fingerprint vectors

    INTEGER, PARAMETER :: nex_cutoff = 2
    INTEGER, PARAMETER :: nwork = 100
    INTEGER :: nenv
    INTEGER :: len_fp  !length of fingerprint
    INTEGER :: ixyzmax
    INTEGER :: ienv
    INTEGER :: llat
    INTEGER :: iorb
    INTEGER :: i
    INTEGER :: j
    INTEGER :: norb
    INTEGER :: nat_sphere
    INTEGER :: iat
    INTEGER :: info
    INTEGER, DIMENSION(nat_sphere_max) :: indat

    REAL(8) :: radius_cutoff
    REAL(8) :: t1
    REAL(8) :: t2

    REAL(8), DIMENSION(3,nat_sphere_max) :: rxyz_sphere
    REAL(8), DIMENSION(nat) :: rcov
    REAL(8), DIMENSION(nat_sphere_max) :: rcov_sphere
    REAL(8), DIMENSION(3,3) :: alatalat
    REAL(8), DIMENSION(3) :: eigalat
    REAL(8), DIMENSION(nwork) :: workalat
    REAL(8), DIMENSION(nat_sphere_max) :: amplitude
    REAL(8), DIMENSION(nat_sphere_max) :: deramplitude
    REAL(8), DIMENSION(nat_sphere_max) :: alpha
    REAL(8), DIMENSION(10) :: cs
    REAL(8), DIMENSION(10) :: cp

    REAL(8), ALLOCATABLE :: ovrlp(:,:)
    REAL(8), ALLOCATABLE :: ampovrlp(:,:)
    REAL(8), ALLOCATABLE :: eval(:)

    CHARACTER(len=2), DIMENSION(nat) :: symb

    len_fp = nat_sphere_max*(ns+3*np)
    nenv = nat

    !calculate covalent radii corresponding to atom type
    DO iat = 1, nat
      CALL sym2rcov(symb(iat), rcov(iat))
    ENDDO
    fp = 0.d0

    !loop over all atoms to get a fingerprint vector for each environment
    DO ienv=1, nenv


      radius_cutoff=sqrt(2.d0*nex_cutoff)*width_cutoff

      IF (alat(1,1)*alat(2,2)*alat(3,3).eq.0.d0 ) THEN ! no periodic boundary condition
          ixyzmax=0
      ELSE  ! periodic boundary conditions
        DO i=1,3
          DO j=1,3
            alatalat(i,j)=alat(1,i)*alat(1,j)+alat(2,i)*alat(2,j)+alat(3,i)*alat(3,j)
          ENDDO
        ENDDO
        CALL dsyev('N', 'L', 3, alatalat, 3, eigalat, workalat, nwork, info)
        ! ixyzmax determines over how many periodiv images one has to search to fill the sphere with atoms
        ixyzmax= int(sqrt(1.d0/eigalat(1))*radius_cutoff) + 1
      ENDIF

      amplitude = 0.d0
      deramplitude = 0.d0
      !determine the atoms which are within the sphere
      CALL atoms_sphere(width_cutoff,nex_cutoff,ienv,llat,ixyzmax,nat,nat_sphere_max,nat_sphere,alat,rxyz,rxyz_sphere, &
                              rcov,rcov_sphere,indat,amplitude,deramplitude)

      DO iat=1,nat_sphere
        alpha(iat)=.5d0/rcov_sphere(iat)**2
      ENDDO
      ! Specify the width of the Gaussians if several Gaussians per l-channel are used
      DO i=1,10
        cs(i)=sqrt(2.d0)**(i-1)
        cp(i)=sqrt(2.d0)**(i-1)
      ENDDO
      norb = (ns+3*np)*nat_sphere


      ALLOCATE(ovrlp(norb,norb))
      ALLOCATE(ampovrlp(norb,norb))
      ALLOCATE(eval(norb))
      ovrlp = 0.d0
      ampovrlp = 0.d0
      eval = 0.d0

      ! calculate overlapmatrix
      CALL crtovrlp(nat_sphere,rxyz_sphere,alpha,cs,cp,ns,np,ovrlp)
      ! multiply cutoff function to overlapmatrix
      CALL multamp(nat_sphere,ovrlp,amplitude,norb,ns,np,ampovrlp)
      ! generate eigenvalues
      CALL diagonalizeMatrix(norb, ampovrlp, eval)

      ! eigenvalues in descending order
      DO i=1,norb/2
         t1=eval(i)
         t2=eval(norb-i+1)
         eval(i)=t2
         eval(norb-i+1)=t1
      ENDDO

      DO iorb=1,norb
        fp(iorb,ienv) = eval(iorb)
      ENDDO

      DEALLOCATE(ovrlp)
      DEALLOCATE(ampovrlp)
      DEALLOCATE(eval)
      !write(*,*) "Fingerprint of environment", ienv, " : DONE"
    ENDDO

  END SUBROUTINE


  SUBROUTINE fingerprint_and_derivatives(nat, nat_sphere_max, ns, np, width_cutoff, alat,rxyz,symb, fp, dfp)
    INTEGER, INTENT(IN) :: nat      !number of atoms in system
    INTEGER, INTENT(IN) :: ns       !s-orbitals
    INTEGER, INTENT(IN) :: np       !p-orbitals
    INTEGER, INTENT(IN) :: nat_sphere_max  !max. number of atoms in sphere
    REAL(8), INTENT(IN) :: width_cutoff  !cutoff radius
    REAL(8), DIMENSION(3,nat), INTENT(IN) :: rxyz  !coorinates
    REAL(8), DIMENSION(3,3), INTENT(IN) :: alat  !unit cell vectors, if a cluster then alat=0.d0
    REAL(8), DIMENSION(nat_sphere_max*(ns+3*np),nat), INTENT(OUT) :: fp  !fingerprint vectors
    REAL(8), DIMENSION(nat_sphere_max*(ns+3*np),3,nat,nat),INTENT(OUT) :: dfp !fingerrprint derivatives

    INTEGER, PARAMETER :: nex_cutoff = 2
    INTEGER, PARAMETER :: nwork = 100
    INTEGER :: nenv
    INTEGER :: len_fp !length of fingerprint
    INTEGER :: ixyzmax
    INTEGER :: ienv
    INTEGER :: llat
    INTEGER :: iorb
    INTEGER :: i
    INTEGER :: iat
    INTEGER :: iiat
    INTEGER :: iiorb
    INTEGER :: info
    INTEGER :: iref
    INTEGER :: j
    INTEGER :: l
    INTEGER :: nat_sphere
    INTEGER :: norb
    INTEGER, DIMENSION(nat_sphere_max) :: indat

    REAL(8) :: radius_cutoff
    REAL(8) :: t1
    REAL(8) :: t2
    REAL(8) :: xl
    REAL(8) :: yl
    REAL(8) :: zl

    REAL(8), DIMENSION(3,nat_sphere_max) :: rxyz_sphere
    REAL(8), DIMENSION(nat) :: rcov
    REAL(8), DIMENSION(nat_sphere_max) :: rcov_sphere
    REAL(8), DIMENSION(3,3) :: alatalat
    REAL(8), DIMENSION(3) :: eigalat
    REAL(8), DIMENSION(nwork) :: workalat
    REAL(8), DIMENSION(nat_sphere_max) :: amplitude
    REAL(8), DIMENSION(nat_sphere_max) :: deramplitude
    REAL(8), DIMENSION(nat_sphere_max) :: alpha
    REAL(8), DIMENSION(10) :: cs
    REAL(8), DIMENSION(10) :: cp

    REAL(8), ALLOCATABLE :: ovrlp(:,:)
    REAL(8), ALLOCATABLE :: ampovrlp(:,:)
    REAL(8), ALLOCATABLE :: eval(:)
    REAL(8), ALLOCATABLE :: evecn(:,:)
    REAL(8), ALLOCATABLE :: dovrlpdr(:,:,:,:)
    REAL(8), ALLOCATABLE :: tmpA(:,:)
    REAL(8), ALLOCATABLE :: devaldr(:,:,:)

    REAL(8), EXTERNAL :: DDOT

    CHARACTER(len=2), DIMENSION(nat) :: symb

    dfp = 0.d0
    len_fp = nat_sphere_max*(ns+3*np)
    nenv = nat

    !calculate covalent radii corresponding to atom type
    DO iat = 1, nat
      CALL sym2rcov(symb(iat), rcov(iat))
    ENDDO
    fp = 0.d0

    !loop over all atoms to get a fingerprint vector for each environment
    DO ienv=1, nenv

      radius_cutoff=sqrt(2.d0*nex_cutoff)*width_cutoff

      IF (alat(1,1)*alat(2,2)*alat(3,3).eq.0.d0 ) THEN ! no periodic boundary condition
          ixyzmax=0
      ELSE  ! periodic boundary conditions
        DO i=1,3
          DO j=1,3
            alatalat(i,j)=alat(1,i)*alat(1,j)+alat(2,i)*alat(2,j)+alat(3,i)*alat(3,j)
          ENDDO
        ENDDO
        CALL dsyev('N', 'L', 3, alatalat, 3, eigalat, workalat, nwork, info)
        ! ixyzmax determines over how many periodiv images one has to search to fill the sphere with atoms
        ixyzmax= int(sqrt(1.d0/eigalat(1))*radius_cutoff) + 1
      ENDIF

      amplitude = 0.d0
      deramplitude = 0.d0
      !determine the atoms which are within the sphere
      CALL atoms_sphere(width_cutoff,nex_cutoff,ienv,llat,ixyzmax,nat,nat_sphere_max,nat_sphere,alat,rxyz,rxyz_sphere, &
                              rcov,rcov_sphere,indat,amplitude,deramplitude)

      DO iat=1,nat_sphere
        alpha(iat)=.5d0/rcov_sphere(iat)**2
      ENDDO
      ! Specify the width of the Gaussians if several Gaussians per l-channel are used
      DO i=1,10
        cs(i)=sqrt(2.d0)**(i-1)
        cp(i)=sqrt(2.d0)**(i-1)
      ENDDO
      norb = (ns+3*np)*nat_sphere


      ALLOCATE(ovrlp(norb,norb))
      ALLOCATE(ampovrlp(norb,norb))
      ALLOCATE(eval(norb))
      ALLOCATE(evecn(norb,norb))
      ALLOCATE(dovrlpdr(norb,norb,3,nat_sphere))
      ALLOCATE(tmpA(norb,norb))
      ALLOCATE(devaldr(3,nat_sphere,norb))
      ovrlp = 0.d0
      ampovrlp = 0.d0
      eval = 0.d0
      evecn = 0.d0
      dovrlpdr = 0.d0
      tmpA = 0.d0
      devaldr = 0.d0
      ! calculate overlapmatrix
      CALL crtovrlp(nat_sphere,rxyz_sphere,alpha,cs,cp,ns,np,ovrlp)
      ! multiply cutoff function to overlapmatrix
      CALL multamp(nat_sphere,ovrlp,amplitude,norb,ns,np,ampovrlp)
      ! generate eigenvalues
      evecn = ampovrlp
      CALL diagonalizeMatrix(norb, evecn, eval)

      ! eigenvalues in descending order
      DO i=1,norb/2
         t1=eval(i)
         t2=eval(norb-i+1)
         eval(i)=t2
         eval(norb-i+1)=t1
      ENDDO

      DO iorb=1,norb
        fp(iorb,ienv) = eval(iorb)
      ENDDO

      xl = rxyz_sphere(1,llat)
      yl = rxyz_sphere(2,llat)
      zl = rxyz_sphere(3,llat)

      CALL xyz2ovrlpdr(norb,nat_sphere,rxyz_sphere,rcov_sphere,amplitude,deramplitude,llat,xl,yl,zl,ns,np,dovrlpdr,ovrlp)
      ! calculate fingerprint derivatives
      DO iat=1,nat_sphere
        DO i=1,3
           CALL DGEMM('N','N',norb,norb,norb,1.d0,dovrlpdr(1,1,i,iat),norb,evecn,norb,0.d0,tmpA,norb)
           DO iorb=1,norb
             iiorb = norb-iorb+1
             devaldr(i,iat,iorb) = DDOT(norb,evecn(1,iiorb),1,tmpA(1,iiorb),1)
           ENDDO
         ENDDO
      ENDDO

      DO l = 1, norb
        DO iat=1,nat_sphere
          iiat=indat(iat)
          dfp(l,1,iiat,ienv)=dfp(l,1,iiat,ienv) + devaldr(1,iat,l)
          dfp(l,2,iiat,ienv)=dfp(l,2,iiat,ienv) + devaldr(2,iat,l)
          dfp(l,3,iiat,ienv)=dfp(l,3,iiat,ienv) + devaldr(3,iat,l)
        ENDDO
      ENDDO
!      do l = 1, norb
!        do iat=1,nat_sphere
!          iiat=indat(iat)
!          print*, l, iat, dfp(ienv,1,iiat,l)
!        enddo
!      enddo
!stop

      DEALLOCATE(ovrlp)
      DEALLOCATE(ampovrlp)
      DEALLOCATE(eval)
      DEALLOCATE(evecn)
      DEALLOCATE(dovrlpdr)
      DEALLOCATE(tmpA)
      DEALLOCATE(devaldr)
      !write(*,*) "Fingerprint of environment", ienv, " : DONE"
    ENDDO



  END SUBROUTINE
! subroutines used for calculations

  subroutine atoms_sphere(width_cutoff,nex_cutoff,lat,llat,ixyzmax,nat,natx_sphere,nat_sphere,alat,rxyz,rxyz_sphere, &
                          rcov,rcov_sphere,indat,amplitude,deramplitude)
    IMPLICIT NONE
    INTEGER :: nex_cutoff
    INTEGER :: lat
    INTEGER :: llat
    INTEGER :: ixyzmax
    INTEGER :: nat
    INTEGER :: natx_sphere
    INTEGER :: nat_sphere
    INTEGER :: ix
    INTEGER :: iy
    INTEGER :: iz
    INTEGER :: jat
    REAL(8) :: xj
    REAL(8) :: yj
    REAL(8) :: zj
    REAL(8) :: radius_cutoff
    REAL(8) :: radius_cutoff2
    REAL(8) :: temp
    REAL(8) :: width_cutoff
    REAL(8) :: dist2
    REAL(8) :: factor_cutoff
    REAL(8), DIMENSION(3,natx_sphere) :: rxyz_sphere
    REAL(8), DIMENSION(natx_sphere) :: rcov_sphere
    REAL(8), DIMENSION(natx_sphere) :: amplitude
    REAL(8), DIMENSION(natx_sphere) :: deramplitude
    REAL(8), DIMENSION(3,nat) :: rxyz
    REAL(8), DIMENSION(nat) :: rcov
    REAL(8), DIMENSION(3,3) :: alat
    INTEGER, DIMENSION(natx_sphere) :: indat


    radius_cutoff=sqrt(2.d0*nex_cutoff)*width_cutoff
    radius_cutoff2=radius_cutoff**2
    factor_cutoff=1.d0/(2.d0*nex_cutoff*width_cutoff**2)
  !!  write(*,*) 'width_cutoff,radius_cutoff',width_cutoff,radius_cutoff

       nat_sphere=0
       DO jat = 1, nat
           DO ix = -ixyzmax,ixyzmax
             DO iy = -ixyzmax,ixyzmax
               DO iz = -ixyzmax,ixyzmax
                  xj = rxyz(1, jat) + ix*alat(1,1)+iy*alat(1,2)+iz*alat(1,3)
                  yj = rxyz(2, jat) + ix*alat(2,1)+iy*alat(2,2)+iz*alat(2,3)
                  zj = rxyz(3, jat) + ix*alat(3,1)+iy*alat(3,2)+iz*alat(3,3)
                  dist2 = (xj-rxyz(1, lat))**2+(yj-rxyz(2, lat))**2+(zj-rxyz(3, lat))**2
                  !write(*,*) xj,rxyz(1, lat),yj,rxyz(2, lat),zj,rxyz(3,lat)


                  IF (dist2.le.radius_cutoff2) THEN
                      nat_sphere=nat_sphere+1
                      IF (nat_sphere.gt.natx_sphere) STOP 'enlarge natx_sphere'
                      !amplitude(nat_sphere)=(1.d0 - dist2*factor_cutoff)**nex_cutoff
                      temp = (1.d0 - dist2*factor_cutoff)**(nex_cutoff-1)
                      amplitude(nat_sphere) = temp*(1.d0 - dist2*factor_cutoff)
                      deramplitude(nat_sphere) = -2.d0*factor_cutoff*nex_cutoff*temp

                      rxyz_sphere(1,nat_sphere)=xj
                      rxyz_sphere(2,nat_sphere)=yj
                      rxyz_sphere(3,nat_sphere)=zj
                      rcov_sphere(nat_sphere)=rcov(jat)
                      indat(nat_sphere)=jat
                      IF (dist2.eq.0.d0) llat=nat_sphere
                  ENDIF
               ENDDO
             ENDDO
           ENDDO
       ENDDO
  end subroutine atoms_sphere


  subroutine crtovrlp(nat,rxyz,alpha,cs,cp,ns,np,ovrlp)
    IMPLICIT NONE
    INTEGER :: nat
    INTEGER :: ns
    INTEGER :: np
    INTEGER :: iat
    INTEGER :: jat
    INTEGER :: iorb
    INTEGER :: jorb
    INTEGER :: is
    INTEGER :: ip
    INTEGER :: js
    INTEGER :: jp
    REAL(8) :: r2
    REAL(8) :: sij
    REAL(8) :: ai
    REAL(8) :: aj
    REAL(8) :: t1
    REAL(8) :: t2
    REAL(8) :: t3
    REAL(8) :: t4
    REAL(8) :: t5
    REAL(8) :: xi
    REAL(8) :: yi
    REAL(8) :: zi
    REAL(8) :: xj
    REAL(8) :: yj
    REAL(8) :: zj
    REAL(8) :: xij
    REAL(8) :: yij
    REAL(8) :: zij
    REAL(8), DIMENSION(3,nat) :: rxyz
    REAL(8), DIMENSION(nat*(ns+3*np),nat*(ns+3*np)) :: ovrlp
    REAL(8), DIMENSION(nat) :: alpha
    REAL(8), DIMENSION(10) :: cs
    REAL(8), DIMENSION(10) :: cp


    IF(ns>10 .or. np > 10) STOP 'ns > 10   .or.  np > 10  !'


   ! 1- setup the overlap matrix

    !  <s|s>
    DO jat=1,nat
      DO js=1,ns
        jorb=(jat-1)*(ns+3*np)+js
        aj=alpha(jat)/cs(js)
        xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

        DO iat=1,nat
          DO is=1,ns
            !!iorb=iat+(is-1)*nat
            iorb=(iat-1)*(ns+3*np)+is
            ai= alpha(iat)/cs(is)
            xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

            xij=xi-xj; yij=yi-yj; zij=zi-zj
            r2=xij**2 + yij**2 + zij**2
            t1=ai*aj
            t2=ai+aj

            ! normalized GTOs:
            sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
            ovrlp(iorb,jorb)=sij

          ENDDO
        ENDDO
      ENDDO
    ENDDO


    !  <pi|sj>
    DO jat=1,nat
      DO js=1,ns

        jorb=(jat-1)*(ns+3*np)+js
        aj=alpha(jat)/cs(js)
        xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

        DO iat=1,nat
          DO ip=1,np
            !!iorb=1+(iat-1)*3+ns*nat + (ip-1)*3*nat
            iorb=(iat-1)*(ns+3*np)+ns+ip
            ai= alpha(iat)/cp(ip)
            xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

            xij=xi-xj; yij=yi-yj; zij=zi-zj
            r2=xij**2 + yij**2 + zij**2

            t1=ai*aj
            t2=ai+aj

            ! normalized GTOs:
            sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
            t3=-2.d0*sqrt(ai)*aj/t2
            ovrlp(iorb  ,jorb  )= t3 * xij *sij
            ovrlp(iorb+1,jorb  )= t3 * yij *sij
            ovrlp(iorb+2,jorb  )= t3 * zij *sij

          ENDDO
        ENDDO
      ENDDO
    ENDDO


    !  <si|pj>
    DO jat=1,nat
      DO jp=1,np

        jorb=(jat-1)*(ns+3*np)+ns+jp
        aj=alpha(jat)/cp(jp)
        xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

        DO iat=1,nat
          DO is=1,ns
            !!iorb=iat+(is-1)*nat
            iorb=(iat-1)*(ns+3*np)+is
            ai= alpha(iat)/cs(is)
            xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

            xij=xi-xj; yij=yi-yj; zij=zi-zj
            r2=xij**2 + yij**2 + zij**2

            t1=ai*aj
            t2=ai+aj

            ! normalized GTOs:
            sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
            t3=+2.d0*sqrt(aj)*ai/t2
            ovrlp(iorb,jorb  )= t3 * xij *sij
            ovrlp(iorb,jorb+1)= t3 * yij *sij
            ovrlp(iorb,jorb+2)= t3 * zij *sij

          ENDDO
        ENDDO
      ENDDO
    ENDDO


    !  <p|p>
    DO jat=1,nat
      DO jp=1,np

        jorb=(jat-1)*(ns+3*np)+ns+jp
        aj=alpha(jat)/cp(jp)
        xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

        DO iat=1,nat
          DO ip=1,np
            iorb=(iat-1)*(ns+3*np)+ns+ip
            ai= alpha(iat)/cp(ip)
            xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

            xij=xi-xj; yij=yi-yj; zij=zi-zj
            r2=xij**2 + yij**2 + zij**2
            t1=ai*aj
            t2=ai+aj

            ! normalized GTOs:
            sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
            t4= 2.d0*sqrt(t1)/t2
            t5=-2.d0*t1/t2

            ovrlp(iorb  ,jorb  )= t4 *(1.d0 + t5* xij* xij)  * sij
            ovrlp(iorb+1,jorb  )= t4 *(       t5* yij* xij)  * sij
            ovrlp(iorb+2,jorb  )= t4 *(       t5* zij* xij)  * sij
            ovrlp(iorb  ,jorb+1)= t4 *(       t5* xij* yij)  * sij
            ovrlp(iorb+1,jorb+1)= t4 *(1.d0+  t5* yij* yij)  * sij
            ovrlp(iorb+2,jorb+1)= t4 *(       t5* zij* yij)  * sij
            ovrlp(iorb  ,jorb+2)= t4 *(       t5* xij* zij)  * sij
            ovrlp(iorb+1,jorb+2)= t4 *(       t5* yij* zij)  * sij
            ovrlp(iorb+2,jorb+2)= t4 *(1.d0+  t5* zij* zij)  * sij

          ENDDO
        ENDDO
      ENDDO
    ENDDO

  end subroutine crtovrlp


  subroutine multamp(nat,ovrlp,amplitude,norb,ns,np,ovrla)
    IMPLICIT NONE
    INTEGER :: nat
    INTEGER :: norb
    INTEGER :: ns
    INTEGER :: np
    INTEGER :: i
    INTEGER :: j
    INTEGER :: iat
    INTEGER :: jat
    INTEGER :: iorb
    INTEGER :: jorb
    REAL(8), DIMENSION(nat) :: amplitude
    REAL(8), DIMENSION(norb, norb) :: ovrlp
    REAL(8), DIMENSION(norb, norb) :: ovrla

    DO jat = 1,nat
      DO j = 1,(ns+3*np)
        jorb = (jat -1)*(ns+3*np) + j
        DO iat = 1,nat
          DO i = 1,(ns+3*np)
            iorb = (iat-1)*(ns+3*np) +i
            ovrla(iorb,jorb) = ovrlp(iorb,jorb)*amplitude(iat)*amplitude(jat)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  end subroutine multamp

  subroutine diagonalizeMatrix(n, aa, eval)
    IMPLICIT NONE
    ! Calling arguments
    INTEGER :: n
    REAL(8),DIMENSION(n,n) :: aa
    REAL(8),DIMENSION(n) :: eval
    ! Local variables
    INTEGER :: info
    INTEGER :: lwork
    REAL(8),DIMENSION(:),ALLOCATABLE:: work

    lwork=100*n
    ALLOCATE(work(lwork))
    CALL dsyev('v','l', n, aa, n, eval, work, lwork, info)
    IF(info/=0) STOP ' ERROR in dsyev'
    DEALLOCATE(work)

  end subroutine diagonalizeMatrix


  subroutine sym2rcov(sym,rcov)
  ! returns the covalent radius of atom with chemical symbol sym
    REAL(8)  :: rcov
    CHARACTER(len=2) :: sym  ! chemical symbol
    SELECT CASE (adjustl(trim(sym)))
       ! covalet radius in Angstrom taken from WebElements: http://www.webelements.com/periodicity/covalent_radius/
    CASE('H')
       rcov= 0.37d0
    CASE('He')
       rcov= 0.32d0
    CASE('Li')
       rcov= 1.34d0
    CASE('Be')
       rcov= 0.90d0
    CASE('B')
       rcov= 0.82d0
    CASE('C')
       rcov= 0.77d0
    CASE('N')
       rcov= 0.75d0
    CASE('O')
       rcov= 0.73d0
    CASE('F')
       rcov= 0.71d0
    CASE('Ne')
       rcov= 0.69d0
    CASE('Na')
       rcov= 1.54d0
    CASE('Mg')
       rcov= 1.30d0
    CASE('Al')
       rcov= 1.18d0
    CASE('Si')
       rcov= 1.11d0
    CASE('P')
       rcov= 1.06d0
    CASE('S')
       rcov= 1.02d0
    CASE('Cl')
       rcov= 0.99d0
    CASE('Ar')
       rcov= 0.97d0
    CASE('K')
       rcov= 1.96d0
    CASE('Ca')
       rcov= 1.74d0
    CASE('Sc')
       rcov= 1.44d0
    CASE('Ti')
       rcov= 1.36d0
    CASE('V')
       rcov= 1.25d0
    CASE('Cr')
       rcov= 1.27d0
    CASE('Mn')
       rcov= 1.39d0
    CASE('Fe')
       rcov= 1.25d0
    CASE('Co')
       rcov= 1.26d0
    CASE('Ni')
       rcov= 1.21d0
    CASE('Cu')
       rcov= 1.38d0
    CASE('Zn')
       rcov= 1.31d0
    CASE('Ga')
       rcov= 1.26d0
    CASE('Ge')
       rcov= 1.22d0
    CASE('As')
       rcov= 1.19d0
    CASE('Se')
       rcov= 1.16d0
    CASE('Br')
       rcov= 1.14d0
    CASE('Kr')
       rcov= 1.10d0
    CASE('Rb')
       rcov= 2.11d0
    CASE('Sr')
       rcov= 1.92d0
    CASE('Y')
       rcov= 1.62d0
    CASE('Zr')
       rcov= 1.48d0
    CASE('Nb')
       rcov= 1.37d0
    CASE('Mo')
       rcov= 1.45d0
    CASE('Tc')
       rcov= 1.56d0
    CASE('Ru')
       rcov= 1.26d0
    CASE('Rh')
       rcov= 1.35d0
    CASE('Pd')
       rcov= 1.31d0
    CASE('Ag')
       rcov= 1.53d0
    CASE('Cd')
       rcov= 1.48d0
    CASE('In')
       rcov= 1.44d0
    CASE('Sn')
       rcov= 1.41d0
    CASE('Sb')
       rcov= 1.38d0
    CASE('Te')
       rcov= 1.35d0
    CASE('I')
       rcov= 1.33d0
    CASE('Xe')
       rcov= 1.30d0
    CASE('Cs')
       rcov= 2.25d0
    CASE('Ba')
       rcov= 1.98d0
    CASE('La')
       rcov= 1.69d0
       !     case('Ce')
       !     case('Pr')
       !     case('Nd')
       !     case('Pm')
       !     case('Sm')
       !     case('Eu')
       !     case('Gd')
       !     case('Tb')
       !     case('Dy')
       !     case('Ho')
       !     case('Er')
       !     case('Tm')
       !     case('Yb')
    CASE('Lu')
       rcov= 1.60d0
    CASE('Hf')
       rcov= 1.50d0
    CASE('Ta')
       rcov= 1.38d0
    CASE('W')
       rcov= 1.46d0
    CASE('Re')
       rcov= 1.59d0
    CASE('Os')
       rcov= 1.28d0
    CASE('Ir')
       rcov= 1.37d0
    CASE('Pt')
       rcov= 1.28d0
    CASE('Au')
       rcov= 1.44d0
    CASE('Hg')
       rcov= 1.49d0
    CASE('Tl')
       rcov= 1.48d0
    CASE('Pb')
       rcov= 1.47d0
    CASE('Bi')
       rcov= 1.46d0
       !     case('Po')
       !     case('At')
    CASE('Rn')
       rcov= 1.45d0
    CASE('LA')   ! Lennard Jones atom
       rcov= 1.122462048309373d0
    CASE('LB')  ! Lennard Jones atom
       rcov= 0.9877666025122482d0
       !     case('Fr')
       !     case('Ra')
       !     case('Ac')
       !     case('Th')
       !     case('Pa')
       !     case('U')
       !     case('Np')
       !     case('Pu')
       !     case('Am')
       !     case('Cm')
    CASE DEFAULT
       PRINT*, " Not recognized atomic type "//sym ; STOP
    ENDSELECT

    rcov = rcov /  0.52917720859d0   ! convert to atomic units
    !rcov = 1.25d0 * rcov

  endsubroutine sym2rcov



  subroutine xyz2ovrlpdr(norb,nat,rxyz,rcov,amplitude,deramplitude,lat,xl,yl,zl,ns,np,dovrlpdr,ovrlp)
  ! Calculqates the derivative of all eigenvalues of an atomic fingerprint with
  ! respect to the atomic positions
  ! Nat=NatSphere da wir eval von sphere davor hatten bzw. für mich sollte es
  ! Eval = FP
    IMPLICIT none
    INTEGER :: norb, nmol, nat, ns, np
    INTEGER :: i, iat, imol, iorb, ip, is
    INTEGER :: jat, jmol, jorb, jp, js
    INTEGER :: lat

    REAL(8) :: derj1, derj2, derj3, derl1, derl2, derl3
    REAL(8) :: derx00i, derx00j, derx00l
    REAL(8) :: derx01i, derx01j, derx01l
    REAL(8) :: derx02i, derx02j, derx02l
    REAL(8) :: derx10i, derx10j, derx10l
    REAL(8) :: deri1, deri2, deri3
    REAL(8) :: derx11i, derx11j, derx11l
    REAL(8) :: derx12i, derx12j, derx12l
    REAL(8) :: derx20i, derx20j, derx20l
    REAL(8) :: derx21i, derx21j, derx21l
    REAL(8) :: derx22i, derx22j, derx22l
    REAL(8) :: dery00i, dery00j, dery00l
    REAL(8) :: dery01i, dery01j, dery01l
    REAL(8) :: dery02i, dery02j, dery02l
    REAL(8) :: dery10i, dery10j, dery10l
    REAL(8) :: dery11i, dery11j, dery11l
    REAL(8) :: dery12i, dery12j, dery12l
    REAL(8) :: dery20i, dery20j, dery20l
    REAL(8) :: dery21i, dery21j, dery21l
    REAL(8) :: dery22i, dery22j, dery22l
    REAL(8) :: derz00i, derz00j, derz00l
    REAL(8) :: derz01i, derz01j, derz01l
    REAL(8) :: derz02i, derz02j, derz02l
    REAL(8) :: derz10i, derz10j, derz10l
    REAL(8) :: derz11i, derz11j, derz11l
    REAL(8) :: derz12i, derz12j, derz12l
    REAL(8) :: derz20i, derz20j, derz20l
    REAL(8) :: derz21i, derz21j, derz21l
    REAL(8) :: derz22i, derz22j, derz22l
    REAL(8) :: dipjampl
    REAL(8) :: djpiampl
    REAL(8) :: pijampl
    REAL(8) :: r2, sij
    REAL(8) :: t1, t2, t3, t4, t5, tt
    REAL(8) :: xi, yi, zi, xj, yj, zj, xl, yl, zl
    REAL(8) :: xij, yij, zij
    REAL(8) :: xil, yil, zil
    REAL(8) :: xjl, yjl, zjl
    REAL(8) :: ss1, ss2, ss3


    REAL(8) :: ai, aj

    INTEGER, DIMENSION(nat) :: ind_molecule
    REAL(8), DIMENSION(norb,norb) :: ovrlp
    REAL(8), DIMENSION(3,nat) :: rxyz
    REAL(8), DIMENSION(nat) :: rcov
    REAL(8), DIMENSION(nat) :: amplitude
    REAL(8), DIMENSION(nat) :: deramplitude
    REAL(8), DIMENSION(norb,norb,3,nat) :: dovrlpdr
    REAL(8), DIMENSION(nat) :: alpha
    REAL(8), DIMENSION(10) :: cs
    REAL(8), DIMENSION(10) :: cp

      DO iat=1,nat
         alpha(iat)=.5d0/rcov(iat)**2
      ENDDO
      ! Specify the width of the Gaussians if several Gaussians per l-channel are
      ! used
      DO i=1,10
        cs(i)=sqrt(2.d0)**(i-1)
        cp(i)=sqrt(2.d0)**(i-1)
      ENDDO



  ! Now calculate derivatives
    !  <s|s>
    DO jat=1,nat
      DO js=1,ns
        jorb=(jat-1)*(ns+3*np)+js
        aj=alpha(jat)/cs(js)

        xj=rxyz(1,jat)
        yj=rxyz(2,jat)
        zj=rxyz(3,jat)

        DO iat= 1,nat
          DO is=1,ns
            iorb=(iat-1)*(ns+3*np)+is
            ai= alpha(iat)/cs(is)

            xi=rxyz(1,iat)
            yi=rxyz(2,iat)
            zi=rxyz(3,iat)

            xij=xi-xj
            yij=yi-yj
            zij=zi-zj

            r2=xij**2 + yij**2 + zij**2
            t1=ai*aj
            t2=ai+aj

            ! derivatives
            tt= -2.d0*t1/t2

            xil = xi - xl
            yil = yi - yl
            zil = zi - zl

            xjl = xj - xl
            yjl = yj - yl
            zjl = zj - zl

              pijampl=amplitude(iat)*amplitude(jat)
              dipjampl=deramplitude(iat)*amplitude(jat)
              djpiampl=deramplitude(jat)*amplitude(iat)

              deri1 =  pijampl*(tt*ovrlp(iorb,jorb)*xij) + dipjampl*ovrlp(iorb,jorb)*xil
              derj1 = -pijampl*(tt*ovrlp(iorb,jorb)*xij) + djpiampl*ovrlp(iorb,jorb)*xjl
              derl1 = -xil*dipjampl*ovrlp(iorb,jorb) - djpiampl*ovrlp(iorb,jorb)*xjl

              deri2 =  pijampl*(tt*ovrlp(iorb,jorb)*yij) + dipjampl*ovrlp(iorb,jorb)*yil
              derj2 = -pijampl*(tt*ovrlp(iorb,jorb)*yij) + djpiampl*ovrlp(iorb,jorb)*yjl
              derl2 = -yil*dipjampl*ovrlp(iorb,jorb) - djpiampl*ovrlp(iorb,jorb)*yjl

              deri3 =  pijampl*(tt*ovrlp(iorb,jorb)*zij) + dipjampl*ovrlp(iorb,jorb)*zil
              derj3 = -pijampl*(tt*ovrlp(iorb,jorb)*zij) + djpiampl*ovrlp(iorb,jorb)*zjl
              derl3 = -zil*dipjampl*ovrlp(iorb,jorb) - djpiampl*ovrlp(iorb,jorb)*zjl


              dovrlpdr(iorb,jorb,1,iat)=dovrlpdr(iorb,jorb,1,iat)+deri1
              dovrlpdr(iorb,jorb,1,jat)=dovrlpdr(iorb,jorb,1,jat)+derj1
              dovrlpdr(iorb,jorb,1,lat)=dovrlpdr(iorb,jorb,1,lat)+derl1


              dovrlpdr(iorb,jorb,2,iat)=dovrlpdr(iorb,jorb,2,iat)+deri2
              dovrlpdr(iorb,jorb,2,jat)=dovrlpdr(iorb,jorb,2,jat)+derj2
              dovrlpdr(iorb,jorb,2,lat)=dovrlpdr(iorb,jorb,2,lat)+derl2


              dovrlpdr(iorb,jorb,3,iat)=dovrlpdr(iorb,jorb,3,iat)+deri3
              dovrlpdr(iorb,jorb,3,jat)=dovrlpdr(iorb,jorb,3,jat)+derj3
              dovrlpdr(iorb,jorb,3,lat)=dovrlpdr(iorb,jorb,3,lat)+derl3

          ENDDO
        ENDDO
      ENDDO
    ENDDO


    IF (np.eq.0) GOTO 1111

    !  <pi|sj>
    DO jat=1,nat ! kat, kat ! 1,nat
      DO js=1,ns

        jorb=(jat-1)*(ns+3*np)+js
        aj=alpha(jat)/cs(js)

        xj=rxyz(1,jat)
        yj=rxyz(2,jat)
        zj=rxyz(3,jat)

        DO iat=1,nat
          DO ip=1,np

            iorb=(iat-1)*(ns+3*np)+ns+ip
            ai= alpha(iat)/cp(ip)

            xi=rxyz(1,iat)
            yi=rxyz(2,iat)
            zi=rxyz(3,iat)

            xij=xi-xj
            yij=yi-yj
            zij=zi-zj

            xil = xi - xl
            yil = yi - yl
            zil = zi - zl

            xjl = xj - xl
            yjl = yj - yl
            zjl = zj - zl

            r2=xij**2 + yij**2 + zij**2

            t1=ai*aj
            t2=ai+aj

            ! normalized GTOs:
            sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
            t3=-2.d0*sqrt(ai)*aj/t2

            ! derivatives
            t5=-2.d0*t1/t2

              pijampl=amplitude(iat)*amplitude(jat)

              derx00i = pijampl*(ovrlp(iorb,jorb)*t5*xij+t3*sij) +  &
                     xil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
              derx00j = -pijampl*(ovrlp(iorb,jorb)*t5*xij+t3*sij) + &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*xjl
              derx00l = -xil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -          &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*xjl

              dery00i = pijampl*(ovrlp(iorb,jorb)*t5*yij) +  &
                     yil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
              dery00j = -pijampl*(ovrlp(iorb,jorb)*t5*yij) + &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*yjl
              dery00l = -yil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*yjl

              derz00i = pijampl*(ovrlp(iorb,jorb)*t5*zij) +  &
                     zil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
              derz00j = -pijampl*(ovrlp(iorb,jorb)*t5*zij) + &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*zjl
              derz00l = -zil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*zjl

              derx10i = pijampl*(ovrlp(iorb+1,jorb)*t5*xij) +  &
                     xil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat)
              derx10j = -pijampl*(ovrlp(iorb+1,jorb)*t5*xij) + &
                      amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*xjl
              derx10l = -xil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*xjl

              dery10i = pijampl*(ovrlp(iorb+1,jorb)*t5*yij+t3*sij) +  &
                     yil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat)
              dery10j = -pijampl*(ovrlp(iorb+1,jorb)*t5*yij+t3*sij) + &
                      amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*yjl
              dery10l = -yil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat) -          &
                      amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*yjl

              derz10i = pijampl*(ovrlp(iorb+1,jorb)*t5*zij) +  &
                     zil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat)
              derz10j = -pijampl*(ovrlp(iorb+1,jorb)*t5*zij) + &
                      amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*zjl
              derz10l = -zil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*zjl

              derx20i = pijampl*(ovrlp(iorb+2,jorb)*t5*xij) +  &
                     xil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat)
              derx20j = -pijampl*(ovrlp(iorb+2,jorb)*t5*xij) + &
                      amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*xjl
              derx20l = -xil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*xjl

              dery20i = pijampl*(ovrlp(iorb+2,jorb)*t5*yij) +  &
                     yil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat)
              dery20j = -pijampl*(ovrlp(iorb+2,jorb)*t5*yij) + &
                      amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*yjl
              dery20l = -yil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*yjl

              derz20i = pijampl*(ovrlp(iorb+2,jorb)*t5*zij+t3*sij) +  &
                     zil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat)
              derz20j = -pijampl*(ovrlp(iorb+2,jorb)*t5*zij+t3*sij) + &
                      amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*zjl
              derz20l = -zil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat) -          &
                      amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*zjl


              dovrlpdr(iorb,jorb,1,iat)=dovrlpdr(iorb,jorb,1,iat)+derx00i
              dovrlpdr(iorb+1,jorb,1,iat)=dovrlpdr(iorb+1,jorb,1,iat)+derx10i
              dovrlpdr(iorb+2,jorb,1,iat)=dovrlpdr(iorb+2,jorb,1,iat)+derx20i

              dovrlpdr(iorb,jorb,1,jat)=dovrlpdr(iorb,jorb,1,jat)+derx00j
              dovrlpdr(iorb+1,jorb,1,jat)=dovrlpdr(iorb+1,jorb,1,jat)+derx10j
              dovrlpdr(iorb+2,jorb,1,jat)=dovrlpdr(iorb+2,jorb,1,jat)+derx20j

              dovrlpdr(iorb,jorb,1,lat)=dovrlpdr(iorb,jorb,1,lat)+derx00l
              dovrlpdr(iorb+1,jorb,1,lat)=dovrlpdr(iorb+1,jorb,1,lat)+derx10l
              dovrlpdr(iorb+2,jorb,1,lat)=dovrlpdr(iorb+2,jorb,1,lat)+derx20l


              dovrlpdr(iorb,jorb,2,iat)=dovrlpdr(iorb,jorb,2,iat)+dery00i
              dovrlpdr(iorb+1,jorb,2,iat)=dovrlpdr(iorb+1,jorb,2,iat)+dery10i
              dovrlpdr(iorb+2,jorb,2,iat)=dovrlpdr(iorb+2,jorb,2,iat)+dery20i

              dovrlpdr(iorb,jorb,2,jat)=dovrlpdr(iorb,jorb,2,jat)+dery00j
              dovrlpdr(iorb+1,jorb,2,jat)=dovrlpdr(iorb+1,jorb,2,jat)+dery10j
              dovrlpdr(iorb+2,jorb,2,jat)=dovrlpdr(iorb+2,jorb,2,jat)+dery20j

              dovrlpdr(iorb,jorb,2,lat)=dovrlpdr(iorb,jorb,2,lat)+dery00l
              dovrlpdr(iorb+1,jorb,2,lat)=dovrlpdr(iorb+1,jorb,2,lat)+dery10l
              dovrlpdr(iorb+2,jorb,2,lat)=dovrlpdr(iorb+2,jorb,2,lat)+dery20l


              dovrlpdr(iorb,jorb,3,iat)=dovrlpdr(iorb,jorb,3,iat)+derz00i
              dovrlpdr(iorb+1,jorb,3,iat)=dovrlpdr(iorb+1,jorb,3,iat)+derz10i
              dovrlpdr(iorb+2,jorb,3,iat)=dovrlpdr(iorb+2,jorb,3,iat)+derz20i

              dovrlpdr(iorb,jorb,3,jat)=dovrlpdr(iorb,jorb,3,jat)+derz00j
              dovrlpdr(iorb+1,jorb,3,jat)=dovrlpdr(iorb+1,jorb,3,jat)+derz10j
              dovrlpdr(iorb+2,jorb,3,jat)=dovrlpdr(iorb+2,jorb,3,jat)+derz20j

              dovrlpdr(iorb,jorb,3,lat)=dovrlpdr(iorb,jorb,3,lat)+derz00l
              dovrlpdr(iorb+1,jorb,3,lat)=dovrlpdr(iorb+1,jorb,3,lat)+derz10l
              dovrlpdr(iorb+2,jorb,3,lat)=dovrlpdr(iorb+2,jorb,3,lat)+derz20l

          ENDDO
        ENDDO
      ENDDO
    ENDDO


    !  <si|pj>
    DO jat=1,nat
      DO jp=1,np

        jorb=(jat-1)*(ns+3*np)+ns+jp
        aj=alpha(jat)/cp(jp)

        xj=rxyz(1,jat)
        yj=rxyz(2,jat)
        zj=rxyz(3,jat)

        DO iat=1,nat
          DO is=1,ns
            iorb=(iat-1)*(ns+3*np)+is
            ai= alpha(iat)/cs(is)

            xi=rxyz(1,iat)
            yi=rxyz(2,iat)
            zi=rxyz(3,iat)

            xij=xi-xj
            yij=yi-yj
            zij=zi-zj

            xil = xi - xl
            yil = yi - yl
            zil = zi - zl

            xjl = xj - xl
            yjl = yj - yl
            zjl = zj - zl

            r2=xij**2 + yij**2 + zij**2

            t1=ai*aj
            t2=ai+aj

            ! normalized GTOs:
            sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
            t3=+2.d0*sqrt(aj)*ai/t2

            ! derivatives
            !tt= -2.d0*t1/t2 * sij
            t5=-2.d0*t1/t2


              pijampl=amplitude(iat)*amplitude(jat)

              derx00i = pijampl*(ovrlp(iorb,jorb)*t5*xij+t3*sij) +  &
                     xil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
              derx00j = -pijampl*(ovrlp(iorb,jorb)*t5*xij+t3*sij) + &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*xjl
              derx00l = -xil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -          &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*xjl

              dery00i = pijampl*(ovrlp(iorb,jorb)*t5*yij) +  &
                     yil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
              dery00j = -pijampl*(ovrlp(iorb,jorb)*t5*yij) + &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*yjl
              dery00l = -yil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*yjl

              derz00i = pijampl*(ovrlp(iorb,jorb)*t5*zij) +  &
                     zil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
              derz00j = -pijampl*(ovrlp(iorb,jorb)*t5*zij) + &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*zjl
              derz00l = -zil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*zjl

              derx01i = pijampl*(ovrlp(iorb,jorb+1)*t5*xij) +  &
                     xil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat)
              derx01j = -pijampl*(ovrlp(iorb,jorb+1)*t5*xij) + &
                      amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*xjl
              derx01l = -xil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*xjl

              dery01i = pijampl*(ovrlp(iorb,jorb+1)*t5*yij+t3*sij) +  &
                     yil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat)
              dery01j = -pijampl*(ovrlp(iorb,jorb+1)*t5*yij+t3*sij) + &
                      amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*yjl
              dery01l = -yil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat) -          &
                      amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*yjl

              derz01i = pijampl*(ovrlp(iorb,jorb+1)*t5*zij) +  &
                     zil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat)
              derz01j = -pijampl*(ovrlp(iorb,jorb+1)*t5*zij) + &
                      amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*zjl
              derz01l = -zil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*zjl

              derx02i = pijampl*(ovrlp(iorb,jorb+2)*t5*xij) +  &
                     xil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat)
              derx02j = -pijampl*(ovrlp(iorb,jorb+2)*t5*xij) + &
                      amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*xjl
              derx02l = -xil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*xjl

              dery02i = pijampl*(ovrlp(iorb,jorb+2)*t5*yij) +  &
                     yil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat)
              dery02j = -pijampl*(ovrlp(iorb,jorb+2)*t5*yij) + &
                      amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*yjl
              dery02l = -yil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*yjl

              derz02i = pijampl*(ovrlp(iorb,jorb+2)*t5*zij+t3*sij) +  &
                     zil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat)
              derz02j = -pijampl*(ovrlp(iorb,jorb+2)*t5*zij+t3*sij) + &
                      amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*zjl
              derz02l = -zil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat) -          &
                      amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*zjl

              dovrlpdr(iorb,jorb,1,iat)=dovrlpdr(iorb,jorb,1,iat)+derx00i
              dovrlpdr(iorb,jorb+1,1,iat)=dovrlpdr(iorb,jorb+1,1,iat)+derx01i
              dovrlpdr(iorb,jorb+2,1,iat)=dovrlpdr(iorb,jorb+2,1,iat)+derx02i

              dovrlpdr(iorb,jorb,1,jat)=dovrlpdr(iorb,jorb,1,jat)+derx00j
              dovrlpdr(iorb,jorb+1,1,jat)=dovrlpdr(iorb,jorb+1,1,jat)+derx01j
              dovrlpdr(iorb,jorb+2,1,jat)=dovrlpdr(iorb,jorb+2,1,jat)+derx02j

              dovrlpdr(iorb,jorb,1,lat)=dovrlpdr(iorb,jorb,1,lat)+derx00l
              dovrlpdr(iorb,jorb+1,1,lat)=dovrlpdr(iorb,jorb+1,1,lat)+derx01l
              dovrlpdr(iorb,jorb+2,1,lat)=dovrlpdr(iorb,jorb+2,1,lat)+derx02l


              dovrlpdr(iorb,jorb,2,iat)=dovrlpdr(iorb,jorb,2,iat)+dery00i
              dovrlpdr(iorb,jorb+1,2,iat)=dovrlpdr(iorb,jorb+1,2,iat)+dery01i
              dovrlpdr(iorb,jorb+2,2,iat)=dovrlpdr(iorb,jorb+2,2,iat)+dery02i

              dovrlpdr(iorb,jorb,2,jat)=dovrlpdr(iorb,jorb,2,jat)+dery00j
              dovrlpdr(iorb,jorb+1,2,jat)=dovrlpdr(iorb,jorb+1,2,jat)+dery01j
              dovrlpdr(iorb,jorb+2,2,jat)=dovrlpdr(iorb,jorb+2,2,jat)+dery02j

              dovrlpdr(iorb,jorb,2,lat)=dovrlpdr(iorb,jorb,2,lat)+dery00l
              dovrlpdr(iorb,jorb+1,2,lat)=dovrlpdr(iorb,jorb+1,2,lat)+dery01l
              dovrlpdr(iorb,jorb+2,2,lat)=dovrlpdr(iorb,jorb+2,2,lat)+dery02l


              dovrlpdr(iorb,jorb,3,iat)=dovrlpdr(iorb,jorb,3,iat)+derz00i
              dovrlpdr(iorb,jorb+1,3,iat)=dovrlpdr(iorb,jorb+1,3,iat)+derz01i
              dovrlpdr(iorb,jorb+2,3,iat)=dovrlpdr(iorb,jorb+2,3,iat)+derz02i

              dovrlpdr(iorb,jorb,3,jat)=dovrlpdr(iorb,jorb,3,jat)+derz00j
              dovrlpdr(iorb,jorb+1,3,jat)=dovrlpdr(iorb,jorb+1,3,jat)+derz01j
              dovrlpdr(iorb,jorb+2,3,jat)=dovrlpdr(iorb,jorb+2,3,jat)+derz02j

              dovrlpdr(iorb,jorb,3,lat)=dovrlpdr(iorb,jorb,3,lat)+derz00l
              dovrlpdr(iorb,jorb+1,3,lat)=dovrlpdr(iorb,jorb+1,3,lat)+derz01l
              dovrlpdr(iorb,jorb+2,3,lat)=dovrlpdr(iorb,jorb+2,3,lat)+derz02l

          ENDDO
        ENDDO
      ENDDO
    ENDDO


    !  <p|p>
    DO jat=1,nat
      DO jp=1,np

        jorb=(jat-1)*(ns+3*np)+ns+jp
        aj=alpha(jat)/cp(jp)

        xj=rxyz(1,jat)
        yj=rxyz(2,jat)
        zj=rxyz(3,jat)

        DO iat=1,nat
          DO ip=1,np
            iorb=(iat-1)*(ns+3*np)+ns+ip
            ai= alpha(iat)/cp(ip)

            xi=rxyz(1,iat)
            yi=rxyz(2,iat)
            zi=rxyz(3,iat)

            xij=xi-xj
            yij=yi-yj
            zij=zi-zj

            r2=xij**2 + yij**2 + zij**2
            t1=ai*aj
            t2=ai+aj

            sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
            t4= 2.d0*sqrt(t1)/t2
            t5=-2.d0*t1/t2

            ! derivatives
            xil = xi - xl; yil = yi - yl; zil = zi - zl
            xjl = xj - xl; yjl = yj - yl; zjl = zj - zl
              pijampl=amplitude(iat)*amplitude(jat)
              !dipjampl=deramplitude(iat)*amplitude(jat)
              !djpiampl=deramplitude(jat)*amplitude(iat)

              derx00i = pijampl*(ovrlp(iorb,jorb)*t5*xij+t4*t5*sij*xij*2.d0) +  &
                     xil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
              derx00j = -pijampl*(ovrlp(iorb,jorb)*t5*xij+t4*t5*sij*xij*2.d0) + &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*xjl
              derx00l = -xil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -                      &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*xjl

              dery00i = pijampl*(ovrlp(iorb,jorb)*t5*yij) +  &
                     yil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
              dery00j = -pijampl*(ovrlp(iorb,jorb)*t5*yij) + &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*yjl
              dery00l = -yil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -      &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*yjl

              derz00i = pijampl*(ovrlp(iorb,jorb)*t5*zij) +  &
                     zil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
              derz00j = -pijampl*(ovrlp(iorb,jorb)*t5*zij) + &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*zjl
              derz00l = -zil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*zjl

              derx10i = pijampl*(ovrlp(iorb+1,jorb)*t5*xij+t4*t5*sij*yij) +  &
                     xil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat)
              derx10j = -pijampl*(ovrlp(iorb+1,jorb)*t5*xij+t4*t5*sij*yij) + &
                      amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*xjl
              derx10l = -xil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat) -                   &
                      amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*xjl

              dery10i = pijampl*(ovrlp(iorb+1,jorb)*t5*yij+t4*t5*sij*xij) +  &
                     yil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat)
              dery10j = -pijampl*(ovrlp(iorb+1,jorb)*t5*yij+t4*t5*sij*xij) + &
                      amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*yjl
              dery10l = -yil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat) -                 &
                      amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*yjl

              derz10i = pijampl*(ovrlp(iorb+1,jorb)*t5*zij) +  &
                     zil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat)
              derz10j = -pijampl*(ovrlp(iorb+1,jorb)*t5*zij) + &
                      amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*zjl
              derz10l = -zil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*zjl

              derx20i = pijampl*(ovrlp(iorb+2,jorb)*t5*xij+t4*t5*sij*zij) +  &
                     xil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat)
              derx20j = -pijampl*(ovrlp(iorb+2,jorb)*t5*xij+t4*t5*sij*zij) + &
                      amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*xjl
              derx20l = -xil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat) -                 &
                      amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*xjl

              dery20i = pijampl*(ovrlp(iorb+2,jorb)*t5*yij) +  &
                     yil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat)
              dery20j = -pijampl*(ovrlp(iorb+2,jorb)*t5*yij) + &
                      amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*yjl
              dery20l = -yil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*yjl

              derz20i = pijampl*(ovrlp(iorb+2,jorb)*t5*zij+t4*t5*sij*xij) +  &
                     zil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat)
              derz20j = -pijampl*(ovrlp(iorb+2,jorb)*t5*zij+t4*t5*sij*xij) + &
                      amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*zjl
              derz20l = -zil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat) -                 &
                      amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*zjl

              derx01i = pijampl*(ovrlp(iorb,jorb+1)*t5*xij+t4*t5*sij*yij) +  &
                     xil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat)
              derx01j = -pijampl*(ovrlp(iorb,jorb+1)*t5*xij+t4*t5*sij*yij) + &
                      amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*xjl
              derx01l = -xil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat) -                 &
                      amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*xjl

              dery01i = pijampl*(ovrlp(iorb,jorb+1)*t5*yij+t4*t5*sij*xij) +  &
                     yil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat)
              dery01j = -pijampl*(ovrlp(iorb,jorb+1)*t5*yij+t4*t5*sij*xij) + &
                      amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*yjl
              dery01l = -yil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat) -                 &
                      amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*yjl

              derz01i = pijampl*(ovrlp(iorb,jorb+1)*t5*zij) +  &
                     zil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat)
              derz01j = -pijampl*(ovrlp(iorb,jorb+1)*t5*zij) + &
                      amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*zjl
              derz01l = -zil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*zjl

              derx11i = pijampl*(ovrlp(iorb+1,jorb+1)*t5*xij) +  &
                     xil*deramplitude(iat)*ovrlp(iorb+1,jorb+1)*amplitude(jat)
              derx11j = -pijampl*(ovrlp(iorb+1,jorb+1)*t5*xij) + &
                      amplitude(iat)*ovrlp(iorb+1,jorb+1)*deramplitude(jat)*xjl
              derx11l = -xil*deramplitude(iat)*ovrlp(iorb+1,jorb+1)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb+1,jorb+1)*deramplitude(jat)*xjl

              dery11i = pijampl*(ovrlp(iorb+1,jorb+1)*t5*yij+t4*t5*sij*yij*2.d0) +  &
                     yil*deramplitude(iat)*ovrlp(iorb+1,jorb+1)*amplitude(jat)
              dery11j = -pijampl*(ovrlp(iorb+1,jorb+1)*t5*yij+t4*t5*sij*yij*2.d0) + &
                      amplitude(iat)*ovrlp(iorb+1,jorb+1)*deramplitude(jat)*yjl
              dery11l = -yil*deramplitude(iat)*ovrlp(iorb+1,jorb+1)*amplitude(jat) -                          &
                      amplitude(iat)*ovrlp(iorb+1,jorb+1)*deramplitude(jat)*yjl

              derz11i = pijampl*(ovrlp(iorb+1,jorb+1)*t5*zij) +  &
                     zil*deramplitude(iat)*ovrlp(iorb+1,jorb+1)*amplitude(jat)
              derz11j = -pijampl*(ovrlp(iorb+1,jorb+1)*t5*zij) + &
                      amplitude(iat)*ovrlp(iorb+1,jorb+1)*deramplitude(jat)*zjl
              derz11l = -zil*deramplitude(iat)*ovrlp(iorb+1,jorb+1)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb+1,jorb+1)*deramplitude(jat)*zjl

              derx21i = pijampl*(ovrlp(iorb+2,jorb+1)*t5*xij) +  &
                     xil*deramplitude(iat)*ovrlp(iorb+2,jorb+1)*amplitude(jat)
              derx21j = -pijampl*(ovrlp(iorb+2,jorb+1)*t5*xij) + &
                      amplitude(iat)*ovrlp(iorb+2,jorb+1)*deramplitude(jat)*xjl
              derx21l = -xil*deramplitude(iat)*ovrlp(iorb+2,jorb+1)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb+2,jorb+1)*deramplitude(jat)*xjl

              dery21i = pijampl*(ovrlp(iorb+2,jorb+1)*t5*yij+t4*t5*sij*zij) +  &
                     yil*deramplitude(iat)*ovrlp(iorb+2,jorb+1)*amplitude(jat)
              dery21j = -pijampl*(ovrlp(iorb+2,jorb+1)*t5*yij+t4*t5*sij*zij) + &
                      amplitude(iat)*ovrlp(iorb+2,jorb+1)*deramplitude(jat)*yjl
              dery21l = -yil*deramplitude(iat)*ovrlp(iorb+2,jorb+1)*amplitude(jat)-                  &
                      amplitude(iat)*ovrlp(iorb+2,jorb+1)*deramplitude(jat)*yjl

              derz21i = pijampl*(ovrlp(iorb+2,jorb+1)*t5*zij+t4*t5*sij*yij) +  &
                     zil*deramplitude(iat)*ovrlp(iorb+2,jorb+1)*amplitude(jat)
              derz21j = -pijampl*(ovrlp(iorb+2,jorb+1)*t5*zij+t4*t5*sij*yij) + &
                      amplitude(iat)*ovrlp(iorb+2,jorb+1)*deramplitude(jat)*zjl
              derz21l = -zil*deramplitude(iat)*ovrlp(iorb+2,jorb+1)*amplitude(jat) -                 &
                      amplitude(iat)*ovrlp(iorb+2,jorb+1)*deramplitude(jat)*zjl

              derx02i = pijampl*(ovrlp(iorb,jorb+2)*t5*xij+t4*t5*sij*zij) +  &
                     xil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat)
              derx02j = -pijampl*(ovrlp(iorb,jorb+2)*t5*xij+t4*t5*sij*zij) + &
                      amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*xjl
              derx02l = -xil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat) -                 &
                      amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*xjl

              dery02i = pijampl*(ovrlp(iorb,jorb+2)*t5*yij) +  &
                     yil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat)
              dery02j = -pijampl*(ovrlp(iorb,jorb+2)*t5*yij) + &
                      amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*yjl
              dery02l = -yil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*yjl

              derz02i = pijampl*(ovrlp(iorb,jorb+2)*t5*zij+t4*t5*sij*xij) +  &
                     zil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat)
              derz02j = -pijampl*(ovrlp(iorb,jorb+2)*t5*zij+t4*t5*sij*xij) + &
                      amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*zjl
              derz02l = -zil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat) -                 &
                      amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*zjl

              derx12i = pijampl*(ovrlp(iorb+1,jorb+2)*t5*xij) +  &
                     xil*deramplitude(iat)*ovrlp(iorb+1,jorb+2)*amplitude(jat)
              derx12j = -pijampl*(ovrlp(iorb+1,jorb+2)*t5*xij) + &
                      amplitude(iat)*ovrlp(iorb+1,jorb+2)*deramplitude(jat)*xjl
              derx12l = -xil*deramplitude(iat)*ovrlp(iorb+1,jorb+2)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb+1,jorb+2)*deramplitude(jat)*xjl

              dery12i = pijampl*(ovrlp(iorb+1,jorb+2)*t5*yij+t4*t5*sij*zij) +  &
                     yil*deramplitude(iat)*ovrlp(iorb+1,jorb+2)*amplitude(jat)
              dery12j = -pijampl*(ovrlp(iorb+1,jorb+2)*t5*yij+t4*t5*sij*zij) + &
                      amplitude(iat)*ovrlp(iorb+1,jorb+2)*deramplitude(jat)*yjl
              dery12l = -yil*deramplitude(iat)*ovrlp(iorb+1,jorb+2)*amplitude(jat) -                 &
                      amplitude(iat)*ovrlp(iorb+1,jorb+2)*deramplitude(jat)*yjl

              derz12i = pijampl*(ovrlp(iorb+1,jorb+2)*t5*zij+t4*t5*sij*yij) +  &
                     zil*deramplitude(iat)*ovrlp(iorb+1,jorb+2)*amplitude(jat)
              derz12j = -pijampl*(ovrlp(iorb+1,jorb+2)*t5*zij+t4*t5*sij*yij) + &
                      amplitude(iat)*ovrlp(iorb+1,jorb+2)*deramplitude(jat)*zjl
              derz12l = -zil*deramplitude(iat)*ovrlp(iorb+1,jorb+2)*amplitude(jat) -                 &
                      amplitude(iat)*ovrlp(iorb+1,jorb+2)*deramplitude(jat)*zjl

              derx22i = pijampl*(ovrlp(iorb+2,jorb+2)*t5*xij) +  &
                     xil*deramplitude(iat)*ovrlp(iorb+2,jorb+2)*amplitude(jat)
              derx22j = -pijampl*(ovrlp(iorb+2,jorb+2)*t5*xij) + &
                      amplitude(iat)*ovrlp(iorb+2,jorb+2)*deramplitude(jat)*xjl
              derx22l = -xil*deramplitude(iat)*ovrlp(iorb+2,jorb+2)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb+2,jorb+2)*deramplitude(jat)*xjl

              dery22i = pijampl*(ovrlp(iorb+2,jorb+2)*t5*yij) +  &
                     yil*deramplitude(iat)*ovrlp(iorb+2,jorb+2)*amplitude(jat)
              dery22j = -pijampl*(ovrlp(iorb+2,jorb+2)*t5*yij) + &
                      amplitude(iat)*ovrlp(iorb+2,jorb+2)*deramplitude(jat)*yjl
              dery22l = -yil*deramplitude(iat)*ovrlp(iorb+2,jorb+2)*amplitude(jat) -   &
                      amplitude(iat)*ovrlp(iorb+2,jorb+2)*deramplitude(jat)*yjl

              derz22i = pijampl*(ovrlp(iorb+2,jorb+2)*t5*zij+t4*t5*sij*zij*2.d0) +  &
                     zil*deramplitude(iat)*ovrlp(iorb+2,jorb+2)*amplitude(jat)
              derz22j = -pijampl*(ovrlp(iorb+2,jorb+2)*t5*zij+t4*t5*sij*zij*2.d0) + &
                      amplitude(iat)*ovrlp(iorb+2,jorb+2)*deramplitude(jat)*zjl
              derz22l = -zil*deramplitude(iat)*ovrlp(iorb+2,jorb+2)*amplitude(jat) -                      &
                      amplitude(iat)*ovrlp(iorb+2,jorb+2)*deramplitude(jat)*zjl


              dovrlpdr(iorb,jorb,1,iat)=dovrlpdr(iorb,jorb,1,iat)+derx00i
              dovrlpdr(iorb+1,jorb,1,iat)=dovrlpdr(iorb+1,jorb,1,iat)+derx10i
              dovrlpdr(iorb+2,jorb,1,iat)=dovrlpdr(iorb+2,jorb,1,iat)+derx20i
              dovrlpdr(iorb,jorb+1,1,iat)=dovrlpdr(iorb,jorb+1,1,iat)+derx01i
              dovrlpdr(iorb,jorb+2,1,iat)=dovrlpdr(iorb,jorb+2,1,iat)+derx02i
              dovrlpdr(iorb+1,jorb+1,1,iat)=dovrlpdr(iorb+1,jorb+1,1,iat)+derx11i
              dovrlpdr(iorb+2,jorb+1,1,iat)=dovrlpdr(iorb+2,jorb+1,1,iat)+derx21i
              dovrlpdr(iorb+1,jorb+2,1,iat)=dovrlpdr(iorb+1,jorb+2,1,iat)+derx12i
              dovrlpdr(iorb+2,jorb+2,1,iat)=dovrlpdr(iorb+2,jorb+2,1,iat)+derx22i

              dovrlpdr(iorb,jorb,1,jat)=dovrlpdr(iorb,jorb,1,jat)+derx00j
              dovrlpdr(iorb+1,jorb,1,jat)=dovrlpdr(iorb+1,jorb,1,jat)+derx10j
              dovrlpdr(iorb+2,jorb,1,jat)=dovrlpdr(iorb+2,jorb,1,jat)+derx20j
              dovrlpdr(iorb,jorb+1,1,jat)=dovrlpdr(iorb,jorb+1,1,jat)+derx01j
              dovrlpdr(iorb,jorb+2,1,jat)=dovrlpdr(iorb,jorb+2,1,jat)+derx02j
              dovrlpdr(iorb+1,jorb+1,1,jat)=dovrlpdr(iorb+1,jorb+1,1,jat)+derx11j
              dovrlpdr(iorb+2,jorb+1,1,jat)=dovrlpdr(iorb+2,jorb+1,1,jat)+derx21j
              dovrlpdr(iorb+1,jorb+2,1,jat)=dovrlpdr(iorb+1,jorb+2,1,jat)+derx12j
              dovrlpdr(iorb+2,jorb+2,1,jat)=dovrlpdr(iorb+2,jorb+2,1,jat)+derx22j

              dovrlpdr(iorb,jorb,1,lat)=dovrlpdr(iorb,jorb,1,lat)+derx00l
              dovrlpdr(iorb+1,jorb,1,lat)=dovrlpdr(iorb+1,jorb,1,lat)+derx10l
              dovrlpdr(iorb+2,jorb,1,lat)=dovrlpdr(iorb+2,jorb,1,lat)+derx20l
              dovrlpdr(iorb,jorb+1,1,lat)=dovrlpdr(iorb,jorb+1,1,lat)+derx01l
              dovrlpdr(iorb,jorb+2,1,lat)=dovrlpdr(iorb,jorb+2,1,lat)+derx02l
              dovrlpdr(iorb+1,jorb+1,1,lat)=dovrlpdr(iorb+1,jorb+1,1,lat)+derx11l
              dovrlpdr(iorb+2,jorb+1,1,lat)=dovrlpdr(iorb+2,jorb+1,1,lat)+derx21l
              dovrlpdr(iorb+1,jorb+2,1,lat)=dovrlpdr(iorb+1,jorb+2,1,lat)+derx12l
              dovrlpdr(iorb+2,jorb+2,1,lat)=dovrlpdr(iorb+2,jorb+2,1,lat)+derx22l


              dovrlpdr(iorb,jorb,2,iat)=dovrlpdr(iorb,jorb,2,iat)+dery00i
              dovrlpdr(iorb+1,jorb,2,iat)=dovrlpdr(iorb+1,jorb,2,iat)+dery10i
              dovrlpdr(iorb+2,jorb,2,iat)=dovrlpdr(iorb+2,jorb,2,iat)+dery20i
              dovrlpdr(iorb,jorb+1,2,iat)=dovrlpdr(iorb,jorb+1,2,iat)+dery01i
              dovrlpdr(iorb,jorb+2,2,iat)=dovrlpdr(iorb,jorb+2,2,iat)+dery02i
              dovrlpdr(iorb+1,jorb+1,2,iat)=dovrlpdr(iorb+1,jorb+1,2,iat)+dery11i
              dovrlpdr(iorb+2,jorb+1,2,iat)=dovrlpdr(iorb+2,jorb+1,2,iat)+dery21i
              dovrlpdr(iorb+1,jorb+2,2,iat)=dovrlpdr(iorb+1,jorb+2,2,iat)+dery12i
              dovrlpdr(iorb+2,jorb+2,2,iat)=dovrlpdr(iorb+2,jorb+2,2,iat)+dery22i

              dovrlpdr(iorb,jorb,2,jat)=dovrlpdr(iorb,jorb,2,jat)+dery00j
              dovrlpdr(iorb+1,jorb,2,jat)=dovrlpdr(iorb+1,jorb,2,jat)+dery10j
              dovrlpdr(iorb+2,jorb,2,jat)=dovrlpdr(iorb+2,jorb,2,jat)+dery20j
              dovrlpdr(iorb,jorb+1,2,jat)=dovrlpdr(iorb,jorb+1,2,jat)+dery01j
              dovrlpdr(iorb,jorb+2,2,jat)=dovrlpdr(iorb,jorb+2,2,jat)+dery02j
              dovrlpdr(iorb+1,jorb+1,2,jat)=dovrlpdr(iorb+1,jorb+1,2,jat)+dery11j
              dovrlpdr(iorb+2,jorb+1,2,jat)=dovrlpdr(iorb+2,jorb+1,2,jat)+dery21j
              dovrlpdr(iorb+1,jorb+2,2,jat)=dovrlpdr(iorb+1,jorb+2,2,jat)+dery12j
              dovrlpdr(iorb+2,jorb+2,2,jat)=dovrlpdr(iorb+2,jorb+2,2,jat)+dery22j

              dovrlpdr(iorb,jorb,2,lat)=dovrlpdr(iorb,jorb,2,lat)+dery00l
              dovrlpdr(iorb+1,jorb,2,lat)=dovrlpdr(iorb+1,jorb,2,lat)+dery10l
              dovrlpdr(iorb+2,jorb,2,lat)=dovrlpdr(iorb+2,jorb,2,lat)+dery20l
              dovrlpdr(iorb,jorb+1,2,lat)=dovrlpdr(iorb,jorb+1,2,lat)+dery01l
              dovrlpdr(iorb,jorb+2,2,lat)=dovrlpdr(iorb,jorb+2,2,lat)+dery02l
              dovrlpdr(iorb+1,jorb+1,2,lat)=dovrlpdr(iorb+1,jorb+1,2,lat)+dery11l
              dovrlpdr(iorb+2,jorb+1,2,lat)=dovrlpdr(iorb+2,jorb+1,2,lat)+dery21l
              dovrlpdr(iorb+1,jorb+2,2,lat)=dovrlpdr(iorb+1,jorb+2,2,lat)+dery12l
              dovrlpdr(iorb+2,jorb+2,2,lat)=dovrlpdr(iorb+2,jorb+2,2,lat)+dery22l


              dovrlpdr(iorb,jorb,3,iat)=dovrlpdr(iorb,jorb,3,iat)+derz00i
              dovrlpdr(iorb+1,jorb,3,iat)=dovrlpdr(iorb+1,jorb,3,iat)+derz10i
              dovrlpdr(iorb+2,jorb,3,iat)=dovrlpdr(iorb+2,jorb,3,iat)+derz20i
              dovrlpdr(iorb,jorb+1,3,iat)=dovrlpdr(iorb,jorb+1,3,iat)+derz01i
              dovrlpdr(iorb,jorb+2,3,iat)=dovrlpdr(iorb,jorb+2,3,iat)+derz02i
              dovrlpdr(iorb+1,jorb+1,3,iat)=dovrlpdr(iorb+1,jorb+1,3,iat)+derz11i
              dovrlpdr(iorb+2,jorb+1,3,iat)=dovrlpdr(iorb+2,jorb+1,3,iat)+derz21i
              dovrlpdr(iorb+1,jorb+2,3,iat)=dovrlpdr(iorb+1,jorb+2,3,iat)+derz12i
              dovrlpdr(iorb+2,jorb+2,3,iat)=dovrlpdr(iorb+2,jorb+2,3,iat)+derz22i

              dovrlpdr(iorb,jorb,3,jat)=dovrlpdr(iorb,jorb,3,jat)+derz00j
              dovrlpdr(iorb+1,jorb,3,jat)=dovrlpdr(iorb+1,jorb,3,jat)+derz10j
              dovrlpdr(iorb+2,jorb,3,jat)=dovrlpdr(iorb+2,jorb,3,jat)+derz20j
              dovrlpdr(iorb,jorb+1,3,jat)=dovrlpdr(iorb,jorb+1,3,jat)+derz01j
              dovrlpdr(iorb,jorb+2,3,jat)=dovrlpdr(iorb,jorb+2,3,jat)+derz02j
              dovrlpdr(iorb+1,jorb+1,3,jat)=dovrlpdr(iorb+1,jorb+1,3,jat)+derz11j
              dovrlpdr(iorb+2,jorb+1,3,jat)=dovrlpdr(iorb+2,jorb+1,3,jat)+derz21j
              dovrlpdr(iorb+1,jorb+2,3,jat)=dovrlpdr(iorb+1,jorb+2,3,jat)+derz12j
              dovrlpdr(iorb+2,jorb+2,3,jat)=dovrlpdr(iorb+2,jorb+2,3,jat)+derz22j

              dovrlpdr(iorb,jorb,3,lat)=dovrlpdr(iorb,jorb,3,lat)+derz00l
              dovrlpdr(iorb+1,jorb,3,lat)=dovrlpdr(iorb+1,jorb,3,lat)+derz10l
              dovrlpdr(iorb+2,jorb,3,lat)=dovrlpdr(iorb+2,jorb,3,lat)+derz20l
              dovrlpdr(iorb,jorb+1,3,lat)=dovrlpdr(iorb,jorb+1,3,lat)+derz01l
              dovrlpdr(iorb,jorb+2,3,lat)=dovrlpdr(iorb,jorb+2,3,lat)+derz02l
              dovrlpdr(iorb+1,jorb+1,3,lat)=dovrlpdr(iorb+1,jorb+1,3,lat)+derz11l
              dovrlpdr(iorb+2,jorb+1,3,lat)=dovrlpdr(iorb+2,jorb+1,3,lat)+derz21l
              dovrlpdr(iorb+1,jorb+2,3,lat)=dovrlpdr(iorb+1,jorb+2,3,lat)+derz12l
              dovrlpdr(iorb+2,jorb+2,3,lat)=dovrlpdr(iorb+2,jorb+2,3,lat)+derz22l

          ENDDO
        ENDDO
      ENDDO
    ENDDO


  1111 CONTINUE

  end subroutine xyz2ovrlpdr


  subroutine cite()
    WRITE(*,*) "Author of this test program: Marco Krummenacher"
    WRITE(*,*) "E-Mail: marco.krummenacher@unibas.ch"
    WRITE(*,*) "Please cite the following papers:"
    WRITE(*,*) "--> L. Zhu et al.  J. Chem. Phys. 144, 034203 (2016)"
    WRITE(*,*) "--> A. Sadeghi et al.  J. Chem. Phys. 139, 184118 (2013)"
  end subroutine

end module fingerprint
