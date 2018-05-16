module fft_tools

!----------------------------------------------------------------------------
! NOTES: 
! The commands fft,ifft and ifft2,fft2 defined here use the FFTW library and 
! are meant to do exactly what their MATLAB equivalent's do
!
!----------------------------------------------------------------------------

use, intrinsic  :: iso_c_binding

implicit none
private
!---------------------------------------------

public :: fft,ifft,fft2,ifft2 

include 'fftw3.f03'

!======================================================================
contains

function fft(f) result(fk)
! calculate 1d complex to complex fft

integer*8                               :: plan 
complex,dimension(:)                    :: f
complex,dimension(size(f))              :: fk
integer                                 :: N

N = size(f)
call dfftw_plan_dft_1d(plan,N,f,fk,FFTW_FORWARD,FFTW_ESTIMATE)
call dfftw_execute_dft(plan, f, fk)
call dfftw_destroy_plan(plan)

end function fft


function ifft(fk) result(f)
! calculate 1d complex to complex ifft

integer*8                               :: plan
complex,dimension(:)                    :: fk
complex,dimension(size(fk))             :: f
integer                                 :: N
real                                    :: scale

N = size(fk)
call dfftw_plan_dft_1d(plan,N,fk,f,FFTW_BACKWARD,FFTW_ESTIMATE)
call dfftw_execute_dft(plan, fk, f)
call dfftw_destroy_plan(plan)

! scale
scale = 1.0/N
f = scale*f

end function ifft


function fft2(f) result(fk)
! calculate 2d complex to complex fft

integer*8                               :: plan 
complex,dimension(:,:)                  :: f
complex,dimension(size(f,1),size(f,2))  :: fk
integer                                 :: Nx,Ny

Nx = size(f,1)
Ny = size(f,2)

call dfftw_plan_dft_2d(plan,Ny,Nx,f,fk,FFTW_FORWARD,FFTW_ESTIMATE)
call dfftw_execute_dft(plan, f, fk)
call dfftw_destroy_plan(plan)

end function fft2


function ifft2(fk) result(f)
! calculate 2d complex to complex ifft

integer*8                                 :: plan
complex,dimension(:,:)                    :: fk
complex,dimension(size(fk,1),size(fk,2))  :: f
integer                                   :: Nx,Ny
real                                      :: scale

Nx = size(fk,1)
Ny = size(fk,2)
call dfftw_plan_dft_2d(plan,Ny,Nx,fk,f,FFTW_BACKWARD,FFTW_ESTIMATE)
call dfftw_execute_dft(plan, fk, f)
call dfftw_destroy_plan(plan)

! scale
scale = 1.0/(Nx*Ny)
f = scale*f

end function ifft2


end module fft_tools
