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

program test_derivative
  use fingerprint
  IMPLICIT NONE
  INTEGER, PARAMETER :: nat = 12    !number of atoms in the unit cell
  INTEGER, PARAMETER :: ifp = 1
  INTEGER, PARAMETER :: kat = 1
  INTEGER, PARAMETER :: nint = 100
  INTEGER, PARAMETER :: ns = 1      !using s-orbitals 1->yes, 0->no
  INTEGER, PARAMETER :: np = 1      !using p-orbitals 1->yes, 0->no
  INTEGER, PARAMETER :: nat_sphere_max = 100 !maximal number of atoms in the sphere
  REAL(8), PARAMETER :: width_cutoff = 5.0 !cutoff radius of the sphere
  INTEGER :: idispl
  INTEGER :: iat
  REAL(8) :: stepsize
  REAL(8) :: pi
  REAL(8) :: fact
  REAL(8) :: fp0
  REAL(8) :: path
  REAL(8) :: sumx
  REAL(8) :: sumy
  REAL(8) :: sumz
  REAL(8) :: t1
  REAL(8) :: t2
  REAL(8) :: t3
  REAL(8), DIMENSION(3,3) :: alat !unit cell vectors
  REAL(8), DIMENSION(3,nat) :: rxyz !atom coordinates
  REAL(8), DIMENSION(3,nat) :: rxyz0 !original atom coordinates
  REAL(8), DIMENSION(3,nat) :: displ
  REAL(8), DIMENSION(3,nat,2) :: dispd
  REAL(8), DIMENSION(nat, nat_sphere_max*(ns+3*np)) :: fp !fingerprint for each environment
  REAL(8), DIMENSION(nat, 3, nat, nat_sphere_max*(ns+3*np)) :: dfp !fingerprint derivative for each environment
  CHARACTER(len=100) :: filename  !name of the file which contains the structure
  CHARACTER(len=2), DIMENSION(nat) :: symb !atomic symbols

  INTEGER :: count1
  INTEGER :: count2
  INTEGER :: count_rate
  INTEGER :: count_max

  REAL(8) :: time1
  REAL(8) :: time2
  REAL(8) :: time

  CALL cpu_time(time1)
  CALL system_clock(count1,count_rate,count_max)
  !reading ascii file...
  !This part can be replaced by any file reader as long as rxyz (coordinates) and
  !symb (atomic symbols) are initialized.
  !Assume that units are in angstroem in the inputfile.
  !After the reading they are converted to atomic units.
  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "READ FILE"
  filename = "../data/LiAl-hydrate.ascii"
  CALL read_ascii(filename, nat, rxyz, alat, symb)
  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"

  pi=4.d0*atan(1.d0)
  stepsize=1.0d0/nint
  rxyz0 = rxyz
  call random_number(dispd)
  dispd=dispd*0.1d0
  path = 0.d0
  fact=2*pi/nint

  OPEN(UNIT=10, FILE="Path_Integration.dat")

  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "STARTING PATH INTEGRATION"
  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "--------------------------------------"

  DO idispl = 0, nint
    IF(mod(idispl,20) .eq. 0) THEN
      WRITE(*,*) "PATH INTEGRATION STEP:", idispl
    ENDIF
    DO iat = 1, nat
      rxyz(1,iat) = rxyz0(1,iat)+dispd(1,iat,1)*sin(2*pi*stepsize*idispl)&
                   +dispd(1,iat,2)*cos(2*pi*stepsize*idispl)
      rxyz(2,iat) = rxyz0(2,iat)+dispd(2,iat,1)*sin(2*pi*stepsize*idispl)&
                   +dispd(2,iat,2)*cos(2*pi*stepsize*idispl)
      rxyz(3,iat) = rxyz0(3,iat)+dispd(3,iat,1)*sin(2*pi*stepsize*idispl)&
                   +dispd(3,iat,2)*cos(2*pi*stepsize*idispl)
      displ(1,iat) = dispd(1,iat,1)*cos(2*pi*stepsize*idispl)&
                   -dispd(1,iat,2)*sin(2*pi*stepsize*idispl)
      displ(2,iat) = dispd(2,iat,1)*cos(2*pi*stepsize*idispl)&
                   -dispd(2,iat,2)*sin(2*pi*stepsize*idispl)
      displ(3,iat) = dispd(3,iat,1)*cos(2*pi*stepsize*idispl)&
                   -dispd(3,iat,2)*sin(2*pi*stepsize*idispl)
    ENDDO

    CALL fingerprint_and_derivatives(nat, nat_sphere_max, ns, np, width_cutoff, alat,rxyz,symb, fp, dfp)

    IF(idispl .eq. 0) THEN
      fp0 = fp(kat,ifp)
    ENDIF

    sumx = 0.d0
    sumy = 0.d0
    sumz = 0.d0
    DO iat = 1, nat
      sumx = sumx + dfp(kat,1,iat,ifp)
      sumy = sumy + dfp(kat,2,iat,ifp)
      sumz = sumz + dfp(kat,3,iat,ifp)
    ENDDO
    IF(sumx**2 .gt. 1.d-12) WRITE(*,*) "sumx", sumx
    IF(sumy**2 .gt. 1.d-12) WRITE(*,*) "sumy", sumy
    IF(sumz**2 .gt. 1.d-12) WRITE(*,*) "sumz", sumz

    IF(idispl .lt. nint) THEN
      t1 = 0.d0
      t2 = 0.d0
      t3 = 0.d0
      DO iat = 1, nat
        t1 = t1 + dfp(kat,1,iat,ifp)*displ(1,iat)
        t2 = t2 + dfp(kat,2,iat,ifp)*displ(2,iat)
        t3 = t3 + dfp(kat,3,iat,ifp)*displ(3,iat)
      ENDDO
    ENDIF

    path = path + (t1+t2+t3)*fact
    WRITE(10,'(I3, e24.7,e24.7)') idispl, fp(kat,ifp), fp0+path

  ENDDO

  CLOSE(10)
  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "RESULTS HAVE BEEN WRITTEN TO FILE: Path_Integration.dat"
  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"

  CALL cpu_time(time2)
  CALL system_clock(count2,count_rate,count_max)

  WRITE(*,*) 'CPU time:    ', time2-time1, " s"
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

end program test_derivative




subroutine read_ascii(filename, nat, rxyz, alat, symb)
  IMPLICIT NONE
  INTEGER :: nat
  INTEGER :: natx
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
  REAL(8), DIMENSION(3,nat) :: rxyz
  REAL(8), DIMENSION(3,3) :: alat
  CHARACTER(len=2), DIMENSION(nat) :: symb
  CHARACTER(len=100) :: filename

  OPEN(UNIT=10, FILE=filename)

  READ(10,*) natx
  IF (nat.gt.natx) STOP 'number of atoms do not match'
  READ(10,*) dxx,dyx,dyy
  READ(10,*) dzx,dzy,dzz
  alat(1,1)=dxx
  alat(2,1)=0.d0
  alat(3,1)=0.d0

  alat(1,2)=dyx
  alat(2,2)=dyy
  alat(3,2)=0.d0

  alat(1,3)=dzx
  alat(2,3)=dzy
  alat(3,3)=dzz
  !here units are converted to atomic units
  convert=1.d0/0.52917720859d0

  DO j=1,3
    DO i=1,3
      alat(i,j)=alat(i,j)*convert
    ENDDO
  ENDDO

  DO i = 1, nat
     READ(10, *) ( rxyz(j, i), j = 1, 3) ,symb(i)
     DO l=1,3
       rxyz(l,i)=rxyz(l,i)*convert
     ENDDO
  ENDDO

  CLOSE(10)
end subroutine
