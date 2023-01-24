!#################################################
! Directionally-split 3D Fast Fourier Transform  #
! Reference: Numerical Recipes, Press et' al'    #
!#################################################

PROGRAM DFT_3d


IMPLICIT NONE
INCLUDE 'mpif.h'



INTEGER, PARAMETER :: nx = 32
INTEGER, PARAMETER :: ny = nx
INTEGER, PARAMETER :: nz = nx
INTEGER, PARAMETER :: nranks_x = 1
INTEGER, PARAMETER :: nranks_y = 1
INTEGER, PARAMETER :: nranks_z = 1
REAL(8), PARAMETER :: TWOPI = 8.d0*ATAN(1.d0)

REAL*8, ALLOCATABLE :: fx(:,:,:,:), fk_re(:,:,:,:), fk_im(:,:,:,:)
REAL*8 :: Lx, Ly, Lz, dx
REAL*8 :: x, y, z
INTEGER :: nmax, ix, iy, iz, kx, ky, kz
REAL*8 :: t1, t2, t3, t4


ALLOCATE(fk_re(1:nx,1:ny,1:nz,3), fk_im(1:nx,1:ny,1:nz,3))
ALLOCATE(fx(1:nx,1:ny,1:nz,3))

! grid size
Lx = 1.d0
Ly = 1.d0
Lz = 1.d0
dx = Lx/DBLE(nx)
nmax = MAX(nx,ny,nz)


! test values
DO iz = 1, nz
    z = (iz-0.5)*dx
    DO iy = 1, ny
        y = (iy-0.5)*dx
        DO ix = 1, nx
            x = (ix-0.5)*dx

            fx(ix,iy,iz,:) = SIN(4*TWOPI*y)

        END DO
    END DO
END DO

PRINT*,'Performing FFT.'

t1 = MPI_WTIME()
CALL fft_3d(fx, fk_re, fk_im)
t2 =  MPI_WTIME()



OPEN(UNIT=12, FILE='output_3d.dat', FORM = 'UNFORMATTED', ACCESS = 'STREAM')

DO iz = 1, nz
DO iy = 1, ny
DO ix = 1, nx
    WRITE(12) fx(ix,iy,iz,1), fk_re(ix,iy,iz,1), fk_im(ix,iy,iz,1)
END DO
END DO
END DO

CLOSE(UNIT=12)

DEALLOCATE(fx, fk_re, fk_im)



PRINT*,''
PRINT*,'FFT time (sec) = ', t2-t1
PRINT*,''
PRINT*,'Done.'
PRINT*,''


CONTAINS


SUBROUTINE fft_3d(fx, fk_re, fk_im)

    REAL*8, INTENT(IN) :: fx(1:nx,1:ny,1:nz,3) 
    REAL*8, INTENT(INOUT) :: fk_re(1:nx,1:ny,1:nz,3), fk_im(1:nx,1:ny,1:nz,3) 
    REAL*8 :: dft_buffer_in_re(1:nx), dft_buffer_in_im(1:nx), dft_buffer_out_re(1:nx),dft_buffer_out_im(1:nx)   
    INTEGER :: kx, ky, kz, i, ix, iy, iz
  
    
    ! clear output array    
    fk_re = 0.d0
    fk_im = 0.d0
        
    PRINT*,' Doing x pass..'    
        
    ! FFT in x-direction
     DO i = 1, 3
        DO iz = 1, nx
            DO iy = 1, nx
                            
                ! copy strips-along x into 1d buffer    
                DO ix = 1, nx
                    dft_buffer_in_re(ix) = fx(ix,iy,iz,i)
                    dft_buffer_in_im(ix) = 0.d0
                END DO
    
                ! perform 1D inverse DFT 
                !CALL direct_idft_1d(dft_buffer_in_re,dft_buffer_in_im,dft_buffer_out_re,dft_buffer_out_im)
                CALL cfft_1d(dft_buffer_in_re,dft_buffer_in_im,dft_buffer_out_re,dft_buffer_out_im)
    
                ! copy back into delv array 
                DO kx = 1, nx
                    fk_re(kx,iy,iz,i) = dft_buffer_out_re(kx)  
                    fk_im(kx,iy,iz,i) = dft_buffer_out_im(kx)  
                END DO

            END DO       
        END DO
    END DO
    

    PRINT*,' Doing y pass..'    


    ! DFT in y-direction
    DO i = 1, 3
        DO iz = 1, nx
            DO kx = 1, nx
                            
                ! copy strips-along y into 1d buffer    
                DO iy = 1, nx
                    dft_buffer_in_re(iy) = fk_re(kx,iy,iz,i)
                    dft_buffer_in_im(iy) = fk_im(kx,iy,iz,i)
                END DO
    
                ! perform 1D inverse DFT 
                !CALL direct_idft_1d(dft_buffer_in_re,dft_buffer_in_im,dft_buffer_out_re,dft_buffer_out_im)
                CALL cfft_1d(dft_buffer_in_re,dft_buffer_in_im,dft_buffer_out_re,dft_buffer_out_im)
                
                ! copy back into delvk array 
                DO ky = 1, nx
                    fk_re(kx,ky,iz,i) = dft_buffer_out_re(ky)  
                    fk_im(kx,ky,iz,i) = dft_buffer_out_im(ky)  
                END DO
                
            
            END DO
        END DO
    END DO
    
    
    PRINT*,' Doing z pass..'    


    ! DFT in z-direction. (also apply the IDFT prefactor of 1/nx*ny*nz)
    DO i = 1, 3  
        DO ky = 1, nx
            DO kx = 1, nx
                            
                ! copy strips-along z into 1d buffer    
                DO iz = 1,nx
                    dft_buffer_in_re(iz) = fk_re(kx,ky,iz,i)
                    dft_buffer_in_im(iz) = fk_im(kx,ky,iz,i)
                END DO
    
                ! perform 1D inverse DFT 
                !CALL direct_idft_1d(dft_buffer_in_re,dft_buffer_in_im,dft_buffer_out_re,dft_buffer_out_im)
                CALL cfft_1d(dft_buffer_in_re,dft_buffer_in_im,dft_buffer_out_re,dft_buffer_out_im)
                
                DO kz = 1, nx
                    fk_re(kx,ky,kz,i) = dft_buffer_out_re(kz)  
                    fk_im(kx,ky,kz,i) = dft_buffer_out_im(kz)               
                END DO
                
            END DO        
        END DO
    END DO


END SUBROUTINE fft_3d



! This subroutine performs a 1D discrete fourier transform via direct summation
SUBROUTINE direct_dft_1d(in_re, in_im, out_re, out_im)


    REAL(8), INTENT(IN) :: in_re(0:nx-1), in_im(0:nx-1)  
    REAL(8), INTENT(INOUT)  :: out_re(-nx/2:-1+nx/2), out_im(-nx/2:-1+nx/2)    
    INTEGER :: sgn = 1   ! DFT for sgn = 1 and Inverse DFT for sgn = -1
    INTEGER :: ix, kx
    REAL(8) :: theta, theta0
    
    theta0 =  TWOPI / DBLE(nx)
    
    ! clear output arrays
    out_re = 0.d0
    out_im = 0.d0
        
    DO kx =-nx/2,-1+nx/2
        DO ix = 0, nx-1 
            theta = theta0 * DBLE(ix * kx)
            out_re(kx) = out_re(kx) + (in_re(ix) * DCOS(theta) + DBLE(sgn) * in_im(ix) * DSIN(theta))
            out_im(kx) = out_im(kx) + (in_im(ix) * DCOS(theta) - DBLE(sgn) * in_re(ix) * DSIN(theta))
        END DO
    END DO   
  

END SUBROUTINE direct_dft_1d


! 1D Daniel-Lanczos FFT algorithm for complex  input
SUBROUTINE cfft_1d(in_re, in_im, out_re, out_im)

    REAL(4), INTENT(IN) :: in_re(1:nx_tot), in_im(1:nx_tot)  
    REAL(4), INTENT(INOUT)  :: out_re(-nx_tot/2:-1+nx_tot/2), out_im(-nx_tot/2:-1+nx_tot/2)    
    INTEGER :: sgn = -1   ! DFT for sgn = 1 and Inverse DFT for sgn = -1
    INTEGER :: ix, kx
 

    REAL(4) :: buffer(1:2*nx_tot)   ! input array gets replaced by output
    REAL(4) :: tempr, tempi, theta, wi, wr, wpi, wpr, wtemp
    INTEGER :: i, j, n, m, mmax, istep
    INTEGER :: i1, i2, i3, i4, n2p3
    REAL*4 :: c1, c2, h1i, h1r, h2i, h2r, wis, wrs
    
 
    ! clear output arrays
    out_re = 0.d0
    out_im = 0.d0
    
    
    ! load input array into work buffer
    ix = 1
    DO i = 1, nx_tot
        buffer(ix)   = in_re(i)   
        buffer(ix+1) = in_im(i) 
        ix = ix+2        
    END DO
  
  
    !***************************************
    ! Sort input array in bit-reversed order  
    !***************************************    
    n = 2*nx_tot
    j = 1
    
    DO i = 1, n, 2
    
        IF(j .GT. i) THEN  ! swap the two complex numbers
            tempr = buffer(j)
            tempi = buffer(j+1)
            buffer(j) = buffer(i)
            buffer(j+1) = buffer(i+1)
            buffer(i) = tempr
            buffer(i+1) = tempi                    
        END IF
    
        m = n/2
        DO WHILE ((m .GT. 2) .AND. (j .GT. m))
            j = j - m 
            m = m / 2
        END DO
        j = j + m
    END DO
      
    !********************************************************************************
    ! Using Danielson-Laczos lemma, compute the DFT by summing up the 1-pt base DFT's
    !********************************************************************************     
    mmax = 2
    
    DO WHILE(n .GT. mmax) 
    
        ! initialize for trigonometric recurrence
        istep = 2 * mmax
        theta = TWOPI / (sgn*mmax)  
        wpr = -2.d0 * SIN(0.5d0 * theta)**2 
        wpi =  SIN(theta)
        wr = 1.d0
        wi = 0.d0
        
        DO m = 1, mmax, 2 
            DO i = m, n, istep
       
                j = i + mmax
                
                ! apply Danielson-Lanczos lemma
                tempr = wr*buffer(j) - wi*buffer(j+1)
                tempi = wr*buffer(j+1) + wi*buffer(j)
                buffer(j) = buffer(i) - tempr 
                buffer(j+1) = buffer(i+1) - tempi 
                buffer(i) = buffer(i) + tempr 
                buffer(i+1) = buffer(i+1) + tempi 
                
            END DO
            
            wtemp = wr
            wr = wr*wpr - wi*wpi + wr
            wi = wi*wpr + wtemp*wpi + wi
            
        END DO
        mmax = istep
        
    END DO
      
    ix = 1
    DO kx = 0, -1+nx_tot/2
        out_re(kx) =  buffer(ix)
        out_im(kx) =  -buffer(ix+1)
        out_re(-kx) = buffer(ix)
        out_im(-kx) = buffer(ix+1)
        ix = ix + 2
    END DO



END SUBROUTINE cfft_1d




END PROGRAM DFT_3d