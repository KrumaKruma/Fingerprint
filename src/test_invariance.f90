program test_invariance
  use rot
  use fingerprint
  use hung
  IMPLICIT NONE
  INTEGER, PARAMETER :: nat = 12    !number of atoms in the unit cell
  INTEGER, PARAMETER :: nconf = 2
  INTEGER, PARAMETER :: ns = 1      !using s-orbitals 1->yes, 0->no
  INTEGER, PARAMETER :: np = 1      !using p-orbitals 1->yes, 0->no
  INTEGER, PARAMETER :: nat_sphere_max = 100 !maximal number of atoms in the sphere
  INTEGER :: iat
  INTEGER :: jat
  INTEGER :: iconf
  INTEGER :: jconf
  INTEGER :: iorb
  INTEGER, DIMENSION(nat) :: iassign
  REAL(8), PARAMETER :: width_cutoff = 5.0 !cutoff radius of the sphere
  REAL(8) :: fpd
  REAL(8) :: tav
  REAL(8) :: tmax
  REAL(8) :: tmin
  REAL(8) :: tt
  REAL(8) :: sum
  REAL(8), DIMENSION(nat,nat) :: cost
  REAL(8), DIMENSION(3,3) :: alat !unit cell vectors
  REAL(8), DIMENSION(3,nat) :: rxyz0 !atom coordinates
  REAL(8), DIMENSION(3,nat) :: rxyz !atom coordinates
  REAL(8), DIMENSION(nat, nat_sphere_max*(ns+3*np),2) :: fp !fingerprint for each environment
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
  filename = "../data/benzene.xyz"
  CALL read_xyz(filename, nat, rxyz, symb)
  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "COMPUTE FINGERPRINT BEFORE ROTATION"

  CALL only_fingerprint(nat, nat_sphere_max, ns, np, width_cutoff, alat,rxyz, symb, fp(1,1,1))
  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"

  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "RANDOM ROTANION OF THE SYSTEM"
  rxyz0 = rxyz
  CALL rotation(nat, rxyz)
  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "COMPUTE FINGERPRINT AFTER ROTATION"
  CALL only_fingerprint(nat, nat_sphere_max, ns, np, width_cutoff, alat,rxyz, symb, fp(1,1,2))
  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "COMPUTE FINGERPRINT DISTANCE AND COORDINATE DISTANCE"

  sum = 0.d0
  DO iat=1,nat
    sum = sum + (rxyz0(1,iat)-rxyz(1,iat))**2 + (rxyz0(2,iat)-rxyz(2,iat))**2 + (rxyz0(3,iat)-rxyz(3,iat))**2
  ENDDO
  WRITE(*,'(A,F10.5)') "COORDINATE DISTANCE OF CONFIGURATION 1 AND 2: ", sum

  DO iconf = 1, nconf
    DO jconf = iconf+1, nconf
      tmin=1.d100
      tmax=-1.d100
      tav=0.d0

      DO iat = 1, nat
        DO jat = 1, nat
          tt=0.d0
          DO iorb=1,(ns+3*np)*nat_sphere_max
            tt=tt+(fp(iat,iorb,iconf)-fp(jat,iorb,jconf))**2
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


  CALL cpu_time(time2)
  CALL system_clock(count2,count_rate,count_max)

  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"
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




END PROGRAM

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
