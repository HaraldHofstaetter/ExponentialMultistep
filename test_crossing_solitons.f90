program test_soliton
    use wavefunctions_fourier1d
    use exponential_multistep
    implicit none

    type(fourier1d) :: m
    type(wf_fourier1d) :: psi
    type(adaptive_adams_lawson_time_stepper) :: time_stepper1
    type(adaptive_adams_exponential_time_stepper) :: time_stepper2

    real(kind=prec) :: err, t0, tend, tol, time_old
    integer :: p, step

    m = fourier1d(1024, -16.0_prec, +16.0_prec)
    call m%set_B(B)
    
    t0 = 0.0_prec
    tend = 5.0_prec

    psi = wf_fourier1d(m)
    call psi%set(psi_0, t0)

    p = 4
    tol = 1e-5_prec

    step = 0
    time_old = psi%time
    time_stepper1 = adaptive_adams_lawson_time_stepper(psi, t0, tend, tol, p)
    do while (.not. time_stepper1%done() )
         call time_stepper1%next
         print *, step, psi%time, psi%time-time_old
         time_old = psi%time
         step = step + 1
    end do


  ! clean up
    call psi%finalize
    call m%finalize

contains   

    function psi_0(x, t)
        complex(kind=prec) :: psi_0
        real(kind=prec), intent(in) :: x
        real(kind=prec), intent(in) :: t
        
        real(prec), parameter :: a1 = 2_prec
        real(prec), parameter :: b1 = 1_prec
        real(prec), parameter :: c1 = 5_prec

        real(prec), parameter :: a2 =  2_prec
        real(prec), parameter :: b2 = -3_prec
        real(prec), parameter :: c2 = -5_prec

        psi_0 = (a1/cosh(a1*(x-c1))) * exp(cmplx(0.0_prec, -b1*x, prec)) &
               +(a2/cosh(a2*(x-c2))) * exp(cmplx(0.0_prec, -b2*x, prec)) 
    end function psi_0

    function B(u, x, t)
        complex(kind=prec) :: B 
        complex(kind=prec), intent(in) :: u
        real(kind=prec), intent(in) :: x
        real(kind=prec), intent(in) :: t

        B = -abs(u)**2*u
    end function B


end program test_soliton
