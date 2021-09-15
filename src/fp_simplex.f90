program fp_simplex
  use fingerprint
  use simplex
  IMPLICIT NONE
  INTEGER, PARAMETER :: nat = 20
  INTEGER, PARAMETER :: nconf = 1000
  INTEGER, PARAMETER :: ns = 1      !using s-orbitals 1->yes, 0->no
  INTEGER, PARAMETER :: np = 1      !using p-orbitals 1->yes, 0->no
  INTEGER, PARAMETER :: nat_sphere_max = 100 !maximal number of atoms in the sphere
  REAL(8), PARAMETER :: width_cutoff = 5.0 !cutoff radius of the sphere
  REAL(8), DIMENSION(3,3,nconf) :: alat
  REAL(8), DIMENSION(3,nat, nconf) :: rxyz !atom coordinates
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: fp
  !REAL(8), DIMENSION(nat, nat_sphere_max*(ns+3*np),nconf) :: fp !fingerprint for each environment
  CHARACTER(len=100) :: filename  !name of the file which contains the structure
  CHARACTER(len=2), DIMENSION(nat,nconf) :: symb !atomic symbols

  LOGICAL, PARAMETER :: only_fp = .false.  !if true the fingerprint is calculated. if false both fingerprint and derivatives are calculated!

  INTEGER :: iconf




  alat = 0.d0
  ALLOCATE(fp(nat, nat_sphere_max*(ns+3*np),nconf))

  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "READ FILE"
  filename = "../data/paracetamol.xyz"
  CALL read_xyz(filename, nat, nconf, rxyz, symb)
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
  DEALLOCATE(fp)







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
