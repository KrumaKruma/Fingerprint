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


program fp_for_crystals
  use fingerprint
  IMPLICIT NONE
  !REAL(8) is used in the present version
  !If single precision is used, the calls to the lapack and blas libriary must be adjusted accordingly
  INTEGER, PARAMETER :: nat = 12    !number of atoms in the unit cell
  INTEGER, PARAMETER :: ns = 1      !using s-orbitals 1->yes, 0->no
  INTEGER, PARAMETER :: np = 1      !using p-orbitals 1->yes, 0->no
  INTEGER, PARAMETER :: nat_sphere_max = 100 !maximal number of atoms in the sphere
  REAL(8), PARAMETER :: width_cutoff = 5.0 !cutoff radius of the sphere
  REAL(8), DIMENSION(3,3) :: alat !unit cell vectors
  REAL(8), DIMENSION(3,nat) :: rxyz !atom coordinates
  REAL(8), DIMENSION(nat, nat_sphere_max*(ns+3*np)) :: fp !fingerprint for each environment
  REAL(8), DIMENSION(nat, 3, nat, nat_sphere_max*(ns+3*np)) :: dfp !fingerprint derivative for each environment
  CHARACTER(len=100) :: filename  !name of the file which contains the structure
  CHARACTER(len=2), DIMENSION(nat) :: symb !atomic symbols

  INTEGER :: i
  INTEGER :: j
  INTEGER :: l
  INTEGER :: iat
  INTEGER :: jat
  INTEGER :: count1
  INTEGER :: count2
  INTEGER :: count_rate
  INTEGER :: count_max

  REAL(8) :: t1
  REAL(8) :: t2
  REAL(8) :: time

  LOGICAL, PARAMETER :: only_fp = .false.  !if true the fingerprint is calculated. if false both fingerprint and derivatives are calculated!

  CALL cpu_time(t1)
  CALL system_clock(count1,count_rate,count_max)
  !reading ascii file...
  !This part can be replaced by any file reader as long as rxyz (coordinates) and
  !symb (atomic symbols) are initialized.
  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "READ FILE"
  filename = "../data/benzene.xyz"
  CALL read_xyz(filename, nat, rxyz, symb)
  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"


  !If the structure is not a crystal then set alat to zero
  alat = 0.d0

  !Either only the fingerprint or the fingerprint and its derivative is computed
  IF (only_fp) THEN
    WRITE(*,*) "--------------------------------------"
    WRITE(*,*) "COMPUTE FINGEREPRINT ONLY"
    CALL only_fingerprint(nat, nat_sphere_max, ns, np, width_cutoff, alat,rxyz, symb, fp)
    WRITE(*,*) "DONE"
    WRITE(*,*) "--------------------------------------"
    WRITE(*,*) "--------------------------------------"
    WRITE(*,*) "WRITE FINGERPRINT TO FILE"
    !Write fingerprint to file
    OPEN(UNIT=10, FILE="fingerprint.dat")
    DO iat = 1, nat
      write(10,*) ( fp(iat, j), j = 1, nat_sphere_max*(ns+3*np))
    ENDDO
    CLOSE(10)
    WRITE(*,*) "DONE"
  ELSE
    WRITE(*,*) "--------------------------------------"
    WRITE(*,*) "COMPUTE FINGERPRINT AND DERIVATIVES"
    CALL fingerprint_and_derivatives(nat, nat_sphere_max, ns, np, width_cutoff, alat,rxyz,symb, fp, dfp)
    WRITE(*,*) "DONE"
    WRITE(*,*) "--------------------------------------"
    WRITE(*,*) "--------------------------------------"
    WRITE(*,*) "WRITE FINGERPRINT AND DERIVATIVES TO FILE"
    !Write fingerprint and its derivative to file
    OPEN(UNIT=10, FILE="fingerprint.dat")
    DO iat = 1, nat
      WRITE(10,*) ( fp(iat, j), j = 1, nat_sphere_max*(ns+3*np))
    ENDDO
    CLOSE(10)
    OPEN(UNIT=10, FILE="fingerprint_derivative.dat")
    DO iat = 1, nat !loop over environments
      DO l = 1, 3   !loop over x,y,z dererivative
        DO jat = 1, nat !loop over all atoms to which it fp is derivated
          WRITE(10,*) (dfp(iat,l,jat,j), j = 1, nat_sphere_max*(ns+3*np))
        ENDDO
      ENDDO
    ENDDO
    CLOSE(10)
    WRITE(*,*) "DONE"
    WRITE(*,*) "--------------------------------------"
  ENDIF

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


subroutine read_xyz(filename, nat, rxyz, symb)
  IMPLICIT NONE
  INTEGER :: nat
  INTEGER :: natx
  INTEGER :: i
  INTEGER :: j
  INTEGER :: l
  REAL(8) :: convert
  REAL(8), DIMENSION(3,nat) :: rxyz
  REAL(8), DIMENSION(3,3) :: alat
  CHARACTER(len=2), DIMENSION(nat) :: symb
  CHARACTER(len=100) :: filename
  CHARACTER(len=100) :: ustring

  OPEN(UNIT=10, FILE=filename)

  READ(10,*) natx
  READ(10,*) ustring
  IF (nat.gt.natx) STOP 'number of atoms do not match'

  convert=1.d0/0.52917720859d0

  DO i = 1, nat
     READ(10, *) symb(i), ( rxyz(j, i), j = 1, 3)
     DO l=1,3
       rxyz(l,i)=rxyz(l,i)*convert
     ENDDO
  ENDDO

  CLOSE(10)
end subroutine
