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



!------------------------------------------------------------------------------
! OM Fingerprint Package
!------------------------------------------------------------------------------
!
! MAIN: fp_distance
!
!> Marco Krummenacher
!> University of Basel
!
! DESCRIPTION:
!> This program calculates the fingerprint distance based on the overlap matrix
!> fingerprint.
!
!------------------------------------------------------------------------------


program fp_distance
  use fingerprint
  use hung
  IMPLICIT NONE
  INTEGER, PARAMETER :: nat = 40    !number of atoms in the unit cell
  INTEGER, PARAMETER :: nconf = 2   ! number of configurations
  INTEGER, PARAMETER :: ns = 1      !using s-orbitals 1->yes, 0->no
  INTEGER, PARAMETER :: np = 1      !using p-orbitals 1->yes, 0->no
  INTEGER, PARAMETER :: nat_sphere_max = 100 !maximal number of atoms in the sphere
  REAL(8), PARAMETER :: width_cutoff = 5.0 !cutoff radius of the sphere
  REAL(8), DIMENSION(3,3,nconf) :: alat !unit cell vectors
  REAL(8), DIMENSION(3,nat, nconf) :: rxyz !atom coordinates
  REAL(8), DIMENSION(nat_sphere_max*(ns+3*np),nat, nconf) :: fp !fingerprint for each environment
  CHARACTER(len=100) :: filename  !name of the file which contains the structure
  CHARACTER(len=2), DIMENSION(nat,nconf) :: symb !atomic symbols

  INTEGER :: iat
  INTEGER :: jat
  INTEGER :: iorb
  INTEGER :: iconf
  INTEGER :: jconf
  INTEGER :: count1
  INTEGER :: count2
  INTEGER :: count_rate
  INTEGER :: count_max
  INTEGER, DIMENSION(nat) :: iassign

  REAL(8) :: tav
  REAL(8) :: tmax
  REAL(8) :: tmin
  REAL(8) :: tt
  REAL(8) :: fpd
  REAL(8) :: t1
  REAL(8) :: t2
  REAL(8) :: time
  REAL(8), DIMENSION(nat,nat) :: cost


  CALL cpu_time(t1)
  CALL system_clock(count1,count_rate,count_max)
  !reading ascii file...
  !This part can be replaced by any file reader as long as rxyz (coordinates) and
  !symb (atomic symbols) are initialized.
  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "READ FILE"
  filename = "../data/perovskite.ascii"
  CALL read_ascii(filename, nat, nconf, rxyz, alat, symb)
  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"

  !compute fingerprint
  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "COMPUTE FINGEREPRINT"
  DO iconf = 1, nconf
    CALL only_fingerprint(nat, nat_sphere_max, ns, np, width_cutoff, alat(1,1,iconf),rxyz(1,1,iconf), symb(1,iconf), fp(1,1,iconf))
  ENDDO
  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"


  !compute fingerprit distance with hungarian algorithm
  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "COMPUTE FINGEREPRINT DISTANCE"

  DO iconf = 1, nconf
    DO jconf = iconf+1, nconf
      tmin=1.d100
      tmax=-1.d100
      tav=0.d0

      DO iat = 1, nat
        DO jat = 1, nat
          tt=0.d0
          DO iorb=1,(ns+3*np)*nat_sphere_max
            tt=tt+(fp(iorb,iat,iconf)-fp(iorb,jat,jconf))**2
          ENDDO
          tt=sqrt(tt)
          tmin=min(tmin,tt)
          tmax=max(tmax,tt)
          tav=tav+tt
          cost(iat,jat)=tt
        ENDDO !jat
      ENDDO !iat
      tav=tav/nat**2
      CALL solveAP(nat, cost, iassign, fpd)
      tt=-1.d100
      DO iat = 1, nat
        tt=max(tt,cost(iat,iassign(iat)))
      ENDDO
      WRITE(*,*) "  "
      WRITE(*,'(A,I2.2,A,I2.2,A,F10.5)') "FIGNERPRINT DISTANCE OF CONFIGURATION ", iconf," AND ", jconf, ": ", fpd
      WRITE(*,*) "  "
    ENDDO !jconf
  ENDDO !iconf


  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"








  WRITE(*,*) "--------------------------------------"
  CALL cpu_time(t2)
  CALL system_clock(count2,count_rate,count_max)

  WRITE(*,*) 'CPU time:    ', t2-t1, " s"
  time=(count2-count1)/float(count_rate)
  WRITE(*,*) 'Elapsed time ', time ," s"
  WRITE(*,*) "--------------------------------------"


  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "SOFTWARE INFORMATION"
  CALL cite()
  WRITE(*,*) "--------------------------------------"

  WRITE(*,*) "======================================"
  WRITE(*,*) "PROGRAM FINISHED SUCCESSFUL"
  WRITE(*,*) "======================================"
end program

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> Reading ascii file
!
!> @param[in] filename filename of the configurations
!> @param[in] nat number of atoms in the file
!> @param[in] nconf number of configurations in the file
!> @param[out] rxyz xyz-positions of the atoms
!> @param[in] alat unit cell vectors; if cluster alat = 0.0
!> @param[in] symb atomic symbols
!---------------------------------------------------------------------------

subroutine read_ascii(filename, nat, nconf, rxyz, alat, symb)
  IMPLICIT NONE
  INTEGER :: nat
  INTEGER :: natx
  INTEGER :: nconf
  INTEGER :: iconf
  INTEGER :: i
  INTEGER :: j
  INTEGER :: l
  REAL(8) :: dxx
  REAL(8) :: dyx
  REAL(8) :: dyy
  REAL(8) :: dzx
  REAL(8) :: dzy
  REAL(8) :: dzz
  REAL(8) :: convert
  REAL(8), DIMENSION(3,nat,nconf) :: rxyz
  REAL(8), DIMENSION(3,3,nconf) :: alat
  CHARACTER(len=2), DIMENSION(nat,nconf) :: symb
  CHARACTER(len=100) :: filename

  OPEN(UNIT=10, FILE=filename)
  DO iconf = 1, nconf
    READ(10,*) natx
    IF (nat.gt.natx) STOP 'number of atoms do not match'
    READ(10,*) dxx,dyx,dyy
    READ(10,*) dzx,dzy,dzz
    alat(1,1,iconf)=dxx
    alat(2,1,iconf)=0.d0
    alat(3,1,iconf)=0.d0

    alat(1,2,iconf)=dyx
    alat(2,2,iconf)=dyy
    alat(3,2,iconf)=0.d0

    alat(1,3,iconf)=dzx
    alat(2,3,iconf)=dzy
    alat(3,3,iconf)=dzz

    convert=1.d0/0.52917720859d0

    DO j=1,3
      DO i=1,3
        alat(i,j,iconf)=alat(i,j,iconf)*convert
      ENDDO
    ENDDO

    DO i = 1, nat
       READ(10, *) ( rxyz(j, i, iconf), j = 1, 3) ,symb(i,iconf)
       DO l=1,3
         rxyz(l,i,iconf)=rxyz(l,i,iconf)*convert
       ENDDO
    ENDDO
  ENDDO
  CLOSE(10)
end subroutine
