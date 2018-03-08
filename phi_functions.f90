module phi_functions
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none

    interface expm1
        module procedure dexpm1
        module procedure cexpm1
    end interface 

    interface phi
        module procedure dphi
        module procedure cphi
    end interface phi

contains

    elemental function dexpm1(x)
        real(c_double), intent(in) :: x
        real(c_double) :: dexpm1

        interface
            pure function dexpm1_c(x) bind(c, name='expm1')
                import c_double
                real(c_double), intent(in), value :: x
                real(c_double) ::  dexpm1_c
            end function dexpm1_c
        end interface

        dexpm1 = dexpm1_c(x)
    end function dexpm1


    elemental subroutine sincos(x, s, c)
        real(c_double), intent(in) :: x
        real(c_double), intent(out) :: c
        real(c_double), intent(out) :: s

        interface
            pure subroutine sincos_c(x, s, c) bind(c, name='sincos')
                import c_double
                real(c_double), intent(in), value :: x
                real(c_double), intent(out) :: c
                real(c_double), intent(out) :: s
            end subroutine sincos_c
        end interface

        call sincos_c(x, s, c)
    end subroutine sincos


    elemental function cexpm1(z)
        complex(c_double), intent(in) :: z
        complex(c_double) :: cexpm1
        real(c_double) :: a, b, em1, s, c, h

        a = real(z, kind=c_double)
        b = imag(z)

        call sincos(b, s, c)
        em1 = dexpm1(a)
        h = c-1.0_c_double
        if (abs(h)<0.01_c_double) then ! more stable evaluation
            h = -s**2/(1.0_c_double+c)
        end if
        cexpm1 = cmplx(em1*c + h, (em1+1.0_c_double)*s, kind=c_double)
    end function cexpm1


    elemental function dphi(x, n)
        real(c_double), intent(in) :: x
        integer, intent(in) :: n
        real(c_double) :: dphi

        integer :: m, k, fac_k
    
        if (n==0) then
            dphi = exp(x)
        else if (abs(x)<=0.5) then ! use taylor series
            if (n==1) then
                m=15
            else if (n==2) then
                m=14
            else if (n<=5) then
                m=13
            else if (n<=8) then
                m=12
            else
                m=11
            end if
            dphi = 1.0
            do k=1,m
                dphi = dphi*x/real(n+m+1-k, kind=c_double) + 1.0_c_double
            end do
            fac_k = 1
            do k=1,n
                fac_k = k*fac_k
            end do
            dphi = dphi/real(fac_k, kind=c_double)
        else ! use recursion formula
            dphi = dexpm1(x)/x
            fac_k = 1            
            do k=1,n-1
                fac_k = k*fac_k
                dphi = (dphi - 1.0_c_double/real(fac_k, kind=c_double))/x
            end do
        end if         
    end function


    elemental function cphi(x, n)
        complex(c_double), intent(in) :: x
        integer, intent(in) :: n
        complex(c_double) :: cphi

        integer :: m, k, fac_k
        
        if (n==0) then
            cphi = exp(x)
        else if (abs(x)<=0.5) then ! use taylor series
            if (n==1) then
                m=15
            else if (n==2) then
                m=14
            else if (n<=5) then
                m=13
            else if (n<=8) then
                m=12
            else
                m=11
            end if
            cphi = 1.0
            do k=1,m
                cphi = cphi*x/real(n+m+1-k, kind=c_double) + 1.0_c_double
            end do
            fac_k = 1
            do k=1,n
                fac_k = k*fac_k
            end do
            cphi = cphi/real(fac_k, kind=c_double)
        else ! use recursion formula        
            cphi = cexpm1(x)/x
            fac_k = 1            
            do k=1,n-1
                fac_k = k*fac_k
                cphi = (cphi - 1.0_c_double/real(fac_k, kind=c_double))/x
            end do   
        end if                 
    end function

end module phi_functions



