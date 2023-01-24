PROGRAM DFT_1d

IMPLICIT NONE

INCLUDE 'mpif.h'


INTEGER, PARAMETER :: nx = 32
INTEGER, PARAMETER :: ny = nx
INTEGER, PARAMETER :: nz = nx
INTEGER, PARAMETER :: k_max = nx/2 
REAL(8), PARAMETER :: TWOPI = 8.d0*ATAN(1.d0)

REAL*8, ALLOCATABLE :: fx_re(:), fx_im(:)
REAL*8, ALLOCATABLE :: fk_re(:), fk_im(:)
REAL*8, ALLOCATABLE :: fk2_re(:), fk2_im(:)
REAL*8 :: Lx, Ly, Lz, dx
REAL*8 :: k, x
REAL*8 :: t1, t2, t3, t4
INTEGER :: ix


ALLOCATE(fx_re(1:nx),fx_im(1:nx))
ALLOCATE(fk_re(-nx/2:-1+nx/2),fk_im(-nx/2:-1+nx/2))
ALLOCATE(fk2_re(-nx/2:-1+nx/2),fk2_im(-nx/2:-1+nx/2))

! grid size
Lx = 1.d0
Ly = 1.d0
Lz = 1.d0
dx = Lx/DBLE(nx)


fx_re = 0.d0
fx_im = 0.d0

! prepare input arrays
DO ix = 1, nx
    
    x = (ix-0.5)*dx
    !IF(x .GT. 0.25d0 .AND. x .LE. 0.76d0) THEN
    !    fx_re(ix-1) = 1.d0
    !ELSE
    !    fx_re(ix-1) = 0.d0
    !END IF

    !fx_re(ix) = 2.0*COS(TWOPI*x) 
    fx_re(ix) = 2.0*SIN(4*TWOPI*x) 
    
END DO

!fx_re = (/ 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0/)

!PRINT*,''
!PRINT*,'fx = ',fx_re
!PRINT*,''


! perform dft
t1 = MPI_WTIME()
CALL direct_dft_1d(fx_re, fk_re, fk_im)
t2 = MPI_WTIME()


t3 = MPI_WTIME()
CALL cfft_1d(fx_re, fx_im, fk2_re, fk2_im)
t4 = MPI_WTIME()

PRINT*,''
PRINT*,'DFT time (sec) = ',t2-t1
PRINT*,'FFT time (sec) = ',t4-t3
PRINT*,''


!PRINT*,'Direct DFT: '
!PRINT*,''
!PRINT*,'fk_re = ',fk_re
!PRINT*,''
!PRINT*,'fk_im = ',fk_im
!PRINT*,''



!PRINT*,'FFT: '
!PRINT*,''
!PRINT*,'fk_re = ',fk_re
!PRINT*,''
!PRINT*,'fk_im = ',fk_im
!PRINT*,''



OPEN(UNIT=12, FILE='output.dat', FORM = 'UNFORMATTED', ACCESS = 'STREAM')

DO ix = 1, nx
    k = ix-1-nx/2 !* twopi / Lx
    x = ix * dx
    WRITE(12) x, fx_re(ix), k,fk_re(k),fk_im(k),fk2_re(k),fk2_im(k) 
END DO

CLOSE(UNIT=12)

DEALLOCATE(fx_re, fx_im, fk_re, fk_im)


PRINT*,'Done.'


CONTAINS


! This subroutine performs a 1D discrete fourier transform via direct summation
SUBROUTINE direct_dft_1d(in_re, out_re, out_im)

    REAL(8), INTENT(IN) :: in_re(0:nx-1)  
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
            out_re(kx) = out_re(kx) + (in_re(ix) * DCOS(theta))
            out_im(kx) = out_im(kx) + (- DBLE(sgn) * in_re(ix) * DSIN(theta))
        END DO
    END DO   
  

END SUBROUTINE direct_dft_1d




! 1D Daniel-Lanczos FFT algorithm for complex  input
SUBROUTINE cfft_1d(in_re, in_im, out_re, out_im)

    REAL(8), INTENT(IN) :: in_re(1:nx), in_im(1:nx)  
    REAL(8), INTENT(INOUT)  :: out_re(-nx/2:-1+nx/2), out_im(-nx/2:-1+nx/2)    
    INTEGER :: sgn = 1   ! DFT for sgn = 1 and Inverse DFT for sgn = -1
    INTEGER :: ix, kx
 

    REAL(8) :: buffer(1:nx)   ! input array gets replaced by output
    REAL(8) :: tempr, tempi, theta, wi, wr, wpi, wpr, wtemp
    INTEGER :: i, j, n, m, mmax, istep
    INTEGER :: i1, i2, i3, i4, n2p3
    REAL*8 :: c1, c2, h1i, h1r, h2i, h2r, wis, wrs
    
 
    ! clear output arrays
    out_re = 0.d0
    out_im = 0.d0
    
    
    ! load input array into work buffer
    ix = 1
    DO i = 1, nx
        buffer(ix)   = in_re(i)   
        buffer(ix+1) = in_im(i) 
        ix = ix+2        
    END DO
  
  
    !***************************************
    ! Sort input array in bit-reversed order  
    !***************************************    
    n = 2*nx
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
    DO kx = 0, -1+nx/2
        out_re(kx) =  buffer(ix)
        out_im(kx) =  -buffer(ix+1)
        out_re(-kx) = buffer(ix)
        out_im(-kx) = buffer(ix+1)
        ix = ix + 2
    END DO

    !out_im(0) = 0.d0


END SUBROUTINE cfft_1d



END PROGRAM DFT_1d