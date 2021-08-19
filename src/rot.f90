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



module rot
  IMPLICIT NONE
CONTAINS


SUBROUTINE torque(nat,rxyz,vxyz,tx,ty,tz)
    ! checks whether the vector vxyz has a rotational (torque) component
  IMPLICIT NONE
  INTEGER :: nat
  INTEGER :: iat
  REAL(8) :: tx
  REAL(8) :: ty
  REAL(8) :: tz
  REAL(8) :: cmx
  REAL(8) :: cmy
  REAL(8) :: cmz
  REAL(8),DIMENSION(3,nat) :: rxyz
  REAL(8),DIMENSION(3,nat) :: vxyz

  ! center of mass
  cmx=0.d0
  cmy=0.d0
  cmz=0.d0
  DO iat=1,nat
     cmx=cmx+rxyz(1,iat)
     cmy=cmy+rxyz(2,iat)
     cmz=cmz+rxyz(3,iat)
  ENDDO
  cmx=cmx/nat
  cmy=cmy/nat
  cmz=cmz/nat

  ! torque
  tx=0.d0 ; ty=0.d0 ; tz=0.d0
  DO iat=1,nat
     tx=tx+(rxyz(2,iat)-cmy)*vxyz(3,iat)-(rxyz(3,iat)-cmz)*vxyz(2,iat)
     ty=ty+(rxyz(3,iat)-cmz)*vxyz(1,iat)-(rxyz(1,iat)-cmx)*vxyz(3,iat)
     tz=tz+(rxyz(1,iat)-cmx)*vxyz(2,iat)-(rxyz(2,iat)-cmy)*vxyz(1,iat)
  ENDDO

END SUBROUTINE torque




SUBROUTINE moment(nat,vxyz,sx,sy,sz)
    ! checks whether the vector vxyz has a translational component
  IMPLICIT NONE
  INTEGER :: nat
  INTEGER :: iat
  REAL(8) :: sx
  REAL(8) :: sy
  REAL(8) :: sz
  REAL(8), DIMENSION(3,nat) :: vxyz

  sx=0.d0
  sy=0.d0
  sz=0.d0
  DO iat=1,nat
     sx=sx+vxyz(1,iat)
     sy=sy+vxyz(2,iat)
     sz=sz+vxyz(3,iat)
  ENDDO

END SUBROUTINE moment


SUBROUTINE rotation(nat, rxyz)
  implicit none

  INTEGER :: i
  INTEGER :: iat
  INTEGER :: nat

  REAL(8) :: factor
  REAL(8) :: pi
  REAL(8) :: alpha
  REAL(8) :: a
  REAL(8) :: b
  REAL(8) :: c

  REAL(8), dimension(3) :: n
  REAL(8), dimension(3,nat) :: rxyz
  REAL(8), dimension(3,nat) :: rxyzr
  REAL(8), dimension(3,3) :: rotmat

  n(1) = 1.d0
  n(2) = 1.d0
  n(3) = 1.d0

  factor = 0.d0
  DO i=1,3
     factor = factor+n(i)**2
  ENDDO
  factor = 1/sqrt(factor)
  DO i=1,3
     n(i) = factor * n(i)
  ENDDO
!  write(*,*) sqrt(n(1)**2+n(2)**2+n(3)**2)
  ! Construct rotation matrix
  pi=4.d0*atan(1.d0)
  call random_number(alpha)
  alpha = alpha*pi

  a = 1.d0-cos(alpha)
  b = cos(alpha)
  c = sin(alpha)

  rotmat(1,1) = n(1)*n(1)*a + b
  rotmat(1,2) = n(1)*n(2)*a - n(3)*c
  rotmat(1,3) = n(1)*n(3)*a + n(2)*c
  rotmat(2,1) = n(2)*n(1)*a + n(3)*c
  rotmat(2,2) = n(2)*n(2)*a + b
  rotmat(2,3) = n(2)*n(3)*a - n(1)*c
  rotmat(3,1) = n(3)*n(1)*a - n(2)*c
  rotmat(3,2) = n(3)*n(2)*a + n(1)*c
  rotmat(3,3) = n(3)*n(3)*a + b


  ! Rotate coordinate

  DO iat=1,nat
    rxyzr(1,iat) = rotmat(1,1)*rxyz(1,iat)+rotmat(1,2)*rxyz(2,iat)+rotmat(1,3)*rxyz(3,iat)
    rxyzr(2,iat) = rotmat(2,1)*rxyz(1,iat)+rotmat(2,2)*rxyz(2,iat)+rotmat(2,3)*rxyz(3,iat)
    rxyzr(3,iat) = rotmat(3,1)*rxyz(1,iat)+rotmat(3,2)*rxyz(2,iat)+rotmat(3,3)*rxyz(3,iat)
  ENDDO

  rxyz(:,:) = rxyzr(:,:)

end subroutine

end module
