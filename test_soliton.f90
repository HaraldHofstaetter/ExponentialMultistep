program test_soliton
    use wavefunctions_fourier1d
    use exponential_multistep
    implicit none

    type(fourier1d) :: m
    type(wf_fourier1d) :: psi, psi_ex
    type(adaptive_adams_lawson_time_stepper) :: time_stepper1
    type(adaptive_adams_exponential_time_stepper) :: time_stepper2

    real(kind=prec) :: err, t0, tend, tol, time_old
    integer :: p, step

    m = fourier1d(1024, -16.0_prec, +16.0_prec)
    call m%set_B(B)
    
    t0 = 0.0_prec
    tend = 1.0_prec

    psi_ex = wf_fourier1d(m)
    call psi_ex%set(soliton, tend)


    psi = wf_fourier1d(m)
    call psi%set(soliton, t0)

    p = 4
    tol = 1e-5_prec

!  ! *** Lawson ***
!
!    call psi%set(soliton, t0)
!    step = 0
!    time_old = psi%time
!    time_stepper1 = adaptive_adams_lawson_time_stepper(psi, t0, tend, tol, p, tol)
!    do while (.not. time_stepper1%done() )
!         call time_stepper1%next
!
!         print *, "step=",step," t=", psi%time, " dt=", psi%time-time_old
!         time_old = psi%time
!         step = step + 1
!    end do
!
!  ! compute the error of the (numerically) propagated solution
!  ! compared to the exact solution
!    err = psi%distance(psi_ex)
!
!    print *, "err = ", err
!
  ! *** Exponential (with phi functions ***

    call psi%set(soliton, t0)
    step = 0
    time_old = psi%time
    time_stepper2 = adaptive_adams_exponential_time_stepper(psi, t0, tend, tol, p, tol)
    do while (.not. time_stepper2%done() )
         call time_stepper2%next

         print *, "step=",step," t=", psi%time, " dt=", psi%time-time_old
         time_old = psi%time
         step = step + 1
    end do

  ! compute the error of the (numerically) propagated solution
  ! compared to the exact solution
    err = psi%distance(psi_ex)

    print *, "err = ", err

  ! clean up
    call psi%finalize
    call psi_ex%finalize
    call m%finalize

contains
    function soliton(x, t)
        complex(kind=prec) :: soliton
        real(kind=prec), intent(in) :: x
        real(kind=prec), intent(in) :: t

        real(prec), parameter :: a = 2_prec
        real(prec), parameter :: b = 1_prec
        real(prec), parameter :: c = 0_prec
        real(kind=prec) :: h
        h = .5_prec*(a**2 - b**2)*t - b*x
        soliton = (a/cosh(a*(b*t+x-c))) * exp(cmplx(0_prec, h, prec)) 
    end function soliton

    function B(u, x, t)
        complex(kind=prec) :: B 
        complex(kind=prec), intent(in) :: u
        real(kind=prec), intent(in) :: x
        real(kind=prec), intent(in) :: t

        real(prec), parameter :: alpha = 0.5_prec

        complex(kind=prec) :: u_ex

        u_ex = soliton(x, t)

        B = -(1.0_prec-alpha)*abs(u_ex)**2*u -alpha*abs(u)**2*u
    end function B

end program test_soliton
