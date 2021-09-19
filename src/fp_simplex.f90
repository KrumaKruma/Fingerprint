program fp_simplex
  use fingerprint
  use simplex
  IMPLICIT NONE
  INTEGER, PARAMETER :: nat = 20
  INTEGER, PARAMETER :: nconf = 1000
  INTEGER, PARAMETER :: ns = 1      !using s-orbitals 1->yes, 0->no
  INTEGER, PARAMETER :: np = 1      !using p-orbitals 1->yes, 0->no
  INTEGER, PARAMETER :: nat_sphere_max = 100 !maximal number of atoms in the sphere
  INTEGER, PARAMETER :: nsimplex = 8
  REAL(8), PARAMETER :: width_cutoff = 5.0 !cutoff radius of the sphere
  REAL(8), PARAMETER :: noise1 = 1.d-2
  REAL(8), PARAMETER :: noise2 = 1.d-2
  REAL(8), DIMENSION(3,3,nconf) :: alat
  REAL(8), DIMENSION(3,nat, nconf) :: rxyz !atom coordinates
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: fp
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: dfp !fingerprint derivative for each environment
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: fp_contracted
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: fpcorner
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: dfp_contracted
  REAL(8), ALLOCATABLE, DIMENSION(:) :: fpd1_list
  REAL(8), ALLOCATABLE, DIMENSION(:) :: fpd2_list
  !REAL(8), DIMENSION(nat, nat_sphere_max*(ns+3*np),nconf) :: fp !fingerprint for each environment
  CHARACTER(len=100) :: filename  !name of the file which contains the structure
  CHARACTER(len=100) :: out_file  !name of the file which contains the structure
  CHARACTER(len=100) :: out_file2  !name of the file which contains the structure
  CHARACTER(len=100) :: out_file3  !name of the file which contains the structure
  CHARACTER(len=2), DIMENSION(nat,nconf) :: symb !atomic symbols

  LOGICAL, PARAMETER :: only_fp = .false.  !if true the fingerprint is calculated. if false both fingerprint and derivatives are calculated!

  INTEGER :: iconf
  INTEGER :: len_fp

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
  alat = 0.d0
  ALLOCATE(fp(len_fp,nat, nconf))
  ALLOCATE(fp_contracted(nsimplex, nat, nconf))
  ALLOCATE(fpcorner(len_fp, 0:nsimplex))
  ALLOCATE(fpd1_list((nconf*(nconf-1)*nat*nat)/2))
  ALLOCATE(fpd2_list((nconf*(nconf-1)*nat*nat)/2))


  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "READ FILE"
  filename = "../data/paracetamol.xyz"
  CALL read_xyz(filename, nat, nconf, rxyz, symb)
  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"

  !compute fingerprint
  IF(only_fp) THEN
    WRITE(*,*) "--------------------------------------"
    WRITE(*,*) "COMPUTE FINGEREPRINT"
    DO iconf = 1, nconf
      CALL only_fingerprint(nat, nat_sphere_max, ns, np, width_cutoff,&
          alat(1,1,iconf),rxyz(1,1,iconf), symb(1,iconf), fp(1,1,iconf))
    ENDDO
    WRITE(*,*) "DONE"
    WRITE(*,*) "--------------------------------------"
  ELSE
    WRITE(*,*) "--------------------------------------"
    WRITE(*,*) "COMPUTE FINGEREPRINT AND DERIVATIVES"
    ALLOCATE(dfp(len_fp,3,nat,nat,nconf))
    DO iconf = 1, nconf
      CALL fingerprint_and_derivatives(nat, nat_sphere_max, ns, np, width_cutoff,&
       alat(1,1,iconf),rxyz(1,1,iconf),symb(1,iconf), fp(1,1,iconf), dfp(1,1,1,1,iconf))
    ENDDO
    WRITE(*,*) "DONE"
    WRITE(*,*) "--------------------------------------"
  ENDIF

  IF(only_fp) THEN
    WRITE(*,*) "--------------------------------------"
    WRITE(*,*) "COMPUTE SIMPLEX CONTRACTION"
    CALL SimplexSparse(len_fp, nat, nconf, nsimplex, fp, fp_contracted, fpcorner)
    WRITE(*,*) "DONE"
    WRITE(*,*) "--------------------------------------"
  ELSE
    WRITE(*,*) "--------------------------------------"
    WRITE(*,*) "COMPUTE SIMPLEX CONTRACTION AND DERIVATIVES"
    ALLOCATE(dfp_contracted(nsimplex, 3, nat, nat, nconf))
    CALL SimplexSparse_derivatve(len_fp, nat, nconf, nsimplex, fp, dfp, fp_contracted, dfp_contracted, fpcorner)
    WRITE(*,*) "DONE"
    WRITE(*,*) "--------------------------------------"
  ENDIF


  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "COMPUTE FINGERPRINT DISTANCES"
  CALL correlation(nconf, nat, len_fp, nsimplex, fp, fp_contracted, fpd1_list, fpd2_list)
  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"

  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "COMPUTE HISTORGRAM OUTPUT"
  out_file = "Correlation.dat"
  CALL histogram(out_file, nat, nconf, noise1, noise2, fpd1_list, fpd2_list)
  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"

  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "WIRTE CONTRACTED FINGERPRINT AND DERIVATIVES TO A FILE"
  out_file2 = "fp_contracted.dat"
  out_file3 = "dfp_contracted.dat"
  CALL write_fp(only_fp,out_file2, out_file3, nat, nconf, nsimplex, fp_contracted, dfp_contracted)
  WRITE(*,*) "DONE"
  WRITE(*,*) "--------------------------------------"

  DEALLOCATE(fp)
  DEALLOCATE(fp_contracted)
  DEALLOCATE(fpcorner)
  DEALLOCATE(fpd1_list)
  DEALLOCATE(fpd2_list)

  IF(.NOT. only_fp) THEN
    DEALLOCATE(dfp)
    DEALLOCATE(dfp_contracted)
  ENDIF


  CALL cpu_time(time2)
  CALL system_clock(count2,count_rate,count_max)

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

subroutine correlation(nconf, nat, len_fp, nsimplex, fp, fp_contracted, fpd1_list, fpd2_list)
  IMPLICIT NONE
  INTEGER :: nconf
  INTEGER :: nat
  INTEGER :: len_fp
  INTEGER :: nsimplex
  INTEGER :: l, iat, jat, iconf, jconf
  INTEGER :: counter
  REAL(8) :: dist1
  REAL(8) :: dist2
  REAL(8), DIMENSION(len_fp,nat, nconf) :: fp
  REAL(8), DIMENSION(nsimplex, nat, nconf) :: fp_contracted
  REAL(8), DIMENSION((nconf*(nconf-1)*nat*nat)/2) :: fpd1_list
  REAL(8), DIMENSION((nconf*(nconf-1)*nat*nat)/2) :: fpd2_list

  counter = 1

  DO iconf = 1, nconf
    DO iat = 1, nat
      DO jconf = iconf+1, nconf
        DO jat = 1, nat

          dist1=0.d0
          DO l = 1, len_fp
            dist1=dist1+(fp(l,iat, iconf)-fp(l,jat, jconf))**2
          ENDDO
          dist1=sqrt(dist1)
          fpd1_list(counter) = dist1
          !!print*, iconf, iat, jconf, jat
          dist2=0.d0
          DO l = 1, nsimplex
            dist2=dist2+(fp_contracted(l,iat, iconf)-fp_contracted(l,jat, jconf))**2
          ENDDO
          dist2=sqrt(dist2)
          fpd2_list(counter) = dist2
          counter = counter + 1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !!print*, "COUNTER:   ", counter
end subroutine


subroutine histogram(out_file, nat, nconf, noise1, noise2, fpd1_list, fpd2_list)
  IMPLICIT NONE
  INTEGER :: nat
  INTEGER :: nconf
  INTEGER :: i, j
  INTEGER :: ifp
  INTEGER :: nfp
  INTEGER :: nd1
  INTEGER :: nd2
  REAL(8) :: noise1
  REAL(8) :: noise2
  REAL(8) :: fp1max
  REAL(8) :: fp2max
  REAL(8) :: fp1
  REAL(8) :: fp2
  REAL(8), DIMENSION((nconf*(nconf-1)*nat*nat)/2) :: fpd1_list
  REAL(8), DIMENSION((nconf*(nconf-1)*nat*nat)/2) :: fpd2_list
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: den
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nden
  CHARACTER(len=100) :: out_file


  nfp = (nconf*(nconf-1)*nat*(nat-1))/2
  fp1max = maxval(fpd1_list)/noise1
  fp2max = maxval(fpd2_list)/noise2

  allocate(den(0:int(fp1max)+1,0:int(fp2max)+1))
  allocate(nden(0:int(fp1max)+1,0:int(fp2max)+1))

  DO i = 0, (int(fp1max)+1)
    DO j = 0, (int(fp2max)+1)
      den(i,j) = .false.
      nden(i,j) = 0
    ENDDO
  ENDDO

  DO ifp = 1, nfp
    fp1=fpd1_list(ifp)/noise1
    fp2=fpd2_list(ifp)/noise2

    nd1=int(fp1)
    nd2=int(fp2)

    den(nd1,nd2) = .true.
    nden(nd1,nd2) = nden(nd1,nd2) + 1

  ENDDO
  print*, out_file
  OPEN(UNIT=20,FILE=out_file)

  DO i = 0, (int(fp1max)+1)
    DO j = 0, (int(fp2max)+1)
      IF(den(i,j)) THEN
        WRITE(20,'(3(1x,e24.17))') i+0.0,j+0.0,nden(i,j)+0.0
      ENDIF
    ENDDO
    WRITE(20,*) " "
  ENDDO

  CLOSE(20)


end subroutine


subroutine write_fp(only_fp, out_file2, out_file3, nat, nconf, nsimplex, fp_contracted, dfp_contracted)
  IMPLICIT NONE
  INTEGER :: nat
  INTEGER :: nconf
  INTEGER :: nsimplex
  INTEGER :: i, l, iat, jat, iconf
  REAL(8), DIMENSION(nsimplex, nat, nconf) :: fp_contracted
  REAL(8), DIMENSION(nsimplex,3, nat, nat, nconf) :: dfp_contracted
  CHARACTER(len=100) :: out_file2
  CHARACTER(len=100) :: out_file3
  LOGICAL :: only_fp

  OPEN(UNIT=23, FILE=out_file2)
  DO iconf = 1, nconf
    DO iat = 1, nat
      WRITE(23,'(50(e17.7))') (fp_contracted(l, iat, iconf), l = 1, nsimplex)
    ENDDO
  ENDDO
  CLOSE(23)

  IF(.NOT. only_fp) THEN
    OPEN(UNIT=45, FILE=out_file3)
    DO iconf = 1, nconf
      DO iat = 1, nat
        DO i = 1, 3
          DO jat = 1, nat
            WRITE(45,'(50(e17.7))') (dfp_contracted(l,i, iat, jat,iconf), l = 1, nsimplex)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    CLOSE(45)
  ENDIF

end subroutine
