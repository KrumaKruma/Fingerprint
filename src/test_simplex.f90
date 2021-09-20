program test_simplex
  use fingerprint
  use simplex
  IMPLICIT NONE
  INTEGER, PARAMETER :: nat = 20
  INTEGER, PARAMETER :: nconf = 1000
  INTEGER, PARAMETER :: nsimplex = 8
  INTEGER, PARAMETER :: ifp = 1
  INTEGER, PARAMETER :: kat = 1
  INTEGER, PARAMETER :: nint = 100
  INTEGER, PARAMETER :: ns = 1      !using s-orbitals 1->yes, 0->no
  INTEGER, PARAMETER :: np = 1      !using p-orbitals 1->yes, 0->no
  INTEGER, PARAMETER :: nat_sphere_max = 100 !maximal number of atoms in the sphere
  REAL(8), PARAMETER :: width_cutoff = 5.0 !cutoff radius of the sphere
  INTEGER :: iat
  INTEGER :: iconf
  INTEGER :: idispl
  INTEGER :: len_fp
  INTEGER :: nconf2
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
  REAL(8), DIMENSION(3,nat) :: rxyz_test !atom coordinates
  REAL(8), DIMENSION(3,nat) :: rxyz0 !atom coordinates
  REAL(8), DIMENSION(3,3) :: alat_test !unit cell vectors
  REAL(8), DIMENSION(3,3,nconf) :: alat !unit cell vectors
  REAL(8), DIMENSION(3,nat,nconf) :: rxyz
  REAL(8), DIMENSION(3,nat) :: displ
  REAL(8), DIMENSION(3,nat,2) :: dispd
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: fp
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: fp_test
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: dfp !fingerprint derivative for each environment
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: dfp_test !fingerprint derivative for each environment
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: fp_contracted
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: fp_contracted_test
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: dfp_contracted
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: dfp_contracted_test
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: fpcorner
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: fpcorner_test

  CHARACTER(len=2), DIMENSION(nat,nconf) :: symb
  CHARACTER(len=2), DIMENSION(nat) :: symb_test
  CHARACTER(len=100) :: filename
  LOGICAL :: ws
  LOGICAL :: rs

  INTEGER :: count1
  INTEGER :: count2
  INTEGER :: count_rate
  INTEGER :: count_max

  REAL(8) :: time1
  REAL(8) :: time2
  REAL(8) :: time

  CALL cpu_time(time1)
  CALL system_clock(count1,count_rate,count_max)

  len_fp = nat_sphere_max*(ns+3*np)

  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "READ FILE"
  filename = "../data/paracetamol.xyz"
  CALL read_xyz(filename, nat, nconf, rxyz, symb)
  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"



  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "COMPUTE FINGEREPRINT AND DERIVATIVES"
  ALLOCATE(fp(len_fp,nat, nconf))
  ALLOCATE(dfp(len_fp,3,nat,nat,nconf))
  DO iconf = 1, nconf
    CALL fingerprint_and_derivatives(nat, nat_sphere_max, ns, np, width_cutoff,&
     alat(1,1,iconf),rxyz(1,1,iconf),symb(1,iconf), fp(1,1,iconf), dfp(1,1,1,1,iconf))
  ENDDO
  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"



  WRITE(*,*) "COMPUTE SIMPLEX CONTRACTION AND DERIVATIVES"
  ws = .true.
  rs = .false.
  ALLOCATE(fpcorner(len_fp, 0:nsimplex))
  ALLOCATE(fp_contracted(nsimplex, nat, nconf))
  ALLOCATE(dfp_contracted(nsimplex, 3, nat, nat, nconf))
  CALL SimplexSparse_derivatve(len_fp, nat, nconf, nsimplex, fp, dfp, fp_contracted, dfp_contracted, fpcorner, ws, rs)
  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"
  DEALLOCATE(fp)
  DEALLOCATE(dfp)
  DEALLOCATE(fp_contracted)
  DEALLOCATE(dfp_contracted)

  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "STARTING PATH INTEGRATION"
  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "--------------------------------------"

  ws = .false.
  rs = .true.

  rxyz_test(:,:) = rxyz(:,:,1)
  symb_test(:) = symb(:,1)
  alat_test(:,:) = alat(:,:,1)
  pi=4.d0*atan(1.d0)
  stepsize=1.0d0/nint
  rxyz0 = rxyz_test
  call random_number(dispd)
  dispd=dispd*0.01d0
  path = 0.d0
  fact=2*pi/nint
  nconf2 = 1

  ALLOCATE(fp_test(len_fp,nat,1))
  ALLOCATE(dfp_test(len_fp,3,nat,nat,1))
  ALLOCATE(fp_contracted_test(nsimplex,nat,1))
  ALLOCATE(dfp_contracted_test(nsimplex,3,nat,nat,1))
  OPEN(UNIT=10, FILE="Path_Integration_Simplex.dat")
  DO idispl = 0, nint
    IF(mod(idispl,20) .eq. 0) THEN
      WRITE(*,*) "PATH INTEGRATION STEP:", idispl
    ENDIF
    DO iat = 1, nat
      rxyz_test(1,iat) = rxyz0(1,iat)+dispd(1,iat,1)*sin(2*pi*stepsize*idispl)&
                   +dispd(1,iat,2)*cos(2*pi*stepsize*idispl)
      rxyz_test(2,iat) = rxyz0(2,iat)+dispd(2,iat,1)*sin(2*pi*stepsize*idispl)&
                   +dispd(2,iat,2)*cos(2*pi*stepsize*idispl)
      rxyz_test(3,iat) = rxyz0(3,iat)+dispd(3,iat,1)*sin(2*pi*stepsize*idispl)&
                   +dispd(3,iat,2)*cos(2*pi*stepsize*idispl)
      displ(1,iat) = dispd(1,iat,1)*cos(2*pi*stepsize*idispl)&
                   -dispd(1,iat,2)*sin(2*pi*stepsize*idispl)
      displ(2,iat) = dispd(2,iat,1)*cos(2*pi*stepsize*idispl)&
                   -dispd(2,iat,2)*sin(2*pi*stepsize*idispl)
      displ(3,iat) = dispd(3,iat,1)*cos(2*pi*stepsize*idispl)&
                   -dispd(3,iat,2)*sin(2*pi*stepsize*idispl)
    ENDDO
    fp_test = 0.d0
    dfp_test = 0.d0
    fp_contracted_test = 0.d0
    dfp_contracted_test = 0.d0
    fpcorner = 0.d0

    CALL fingerprint_and_derivatives(nat, nat_sphere_max, ns, np, width_cutoff,&
      alat_test,rxyz_test,symb_test, fp_test(1,1,1), dfp_test(1,1,1,1,1))

    CALL SimplexSparse_derivatve(len_fp, nat, nconf2, nsimplex, fp_test,&
       dfp_test, fp_contracted_test, dfp_contracted_test, fpcorner, ws, rs)

    IF(idispl .eq. 0) THEN
      fp0 = fp_contracted_test(ifp,kat,1)
    ENDIF

    sumx = 0.d0
    sumy = 0.d0
    sumz = 0.d0
    DO iat = 1, nat
      sumx = sumx + dfp_contracted_test(ifp,1,iat,kat,1)
      sumy = sumy + dfp_contracted_test(ifp,2,iat,kat,1)
      sumz = sumz + dfp_contracted_test(ifp,3,iat,kat,1)
    ENDDO
    IF(sumx**2 .gt. 1.d-12) WRITE(*,*) "sumx", sumx
    IF(sumy**2 .gt. 1.d-12) WRITE(*,*) "sumy", sumy
    IF(sumz**2 .gt. 1.d-12) WRITE(*,*) "sumz", sumz

    IF(idispl .lt. nint) THEN
      t1 = 0.d0
      t2 = 0.d0
      t3 = 0.d0
      DO iat = 1, nat
        t1 = t1 + dfp_contracted_test(ifp,1,iat,kat,1)*displ(1,iat)
        t2 = t2 + dfp_contracted_test(ifp,2,iat,kat,1)*displ(2,iat)
        t3 = t3 + dfp_contracted_test(ifp,3,iat,kat,1)*displ(3,iat)
      ENDDO
    ENDIF

    path = path + (t1+t2+t3)*fact
    WRITE(10,'(I3, e24.7,e24.7)') idispl, fp_contracted_test(ifp,kat,1), fp0+path

  ENDDO

  CLOSE(10)

  DEALLOCATE(fpcorner)
  DEALLOCATE(fp_test)
  DEALLOCATE(dfp_test)
  DEALLOCATE(fp_contracted_test)
  DEALLOCATE(dfp_contracted_test)

  CALL cpu_time(time2)
  CALL system_clock(count2,count_rate,count_max)


  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "RESULTS HAVE BEEN WRITTEN TO FILE: Path_Integration_Simplex.dat"
  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"

  WRITE(*,*) 'CPU time:    ', time2-time1, " s"
  time=(count2-count1)/float(count_rate)
  WRITE(*,*) 'Elapsed time ', time ," s"
  WRITE(*,*) "--------------------------------------"


  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "SOFTWARE INFORMATION"
  CALL cite()
  CALL cite_simplex()
  WRITE(*,*) "--------------------------------------"

  WRITE(*,*) "======================================"
  WRITE(*,*) "PROGRAM FINISHED SUCCESSFUL"
  WRITE(*,*) "======================================"

end program


subroutine read_xyz(filename, nat, nconf, rxyz, symb)
  IMPLICIT NONE
  INTEGER :: nat
  INTEGER :: natx
  INTEGER :: nconf
  INTEGER :: iconf
  INTEGER :: i
  INTEGER :: j
  INTEGER :: l
  REAL(8) :: convert
  REAL(8) :: energy
  REAL(8), DIMENSION(3,nat,nconf) :: rxyz
  CHARACTER(len=2), DIMENSION(nat,nconf) :: symb
  CHARACTER(len=100) :: usting
  CHARACTER(len=100) :: filename

  OPEN(UNIT=10, FILE=filename)
  DO iconf = 1, nconf
    READ(10,*) natx, usting, energy
    READ(10,*) usting
    IF (nat.gt.natx) STOP 'number of atoms do not match'

    convert=1.d0/0.52917720859d0

    DO i = 1, nat
       READ(10, *) symb(i,iconf),( rxyz(j, i, iconf), j = 1, 3)
       DO l=1,3
         rxyz(l,i,iconf)=rxyz(l,i,iconf)*convert
       ENDDO
    ENDDO
  ENDDO
  CLOSE(10)
end subroutine
