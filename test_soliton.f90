program test_soliton
    use wavefunctions_fourier1d
    use exponential_multistep
    implicit none

    type(fourier1d) :: m
    type(wf_fourier1d) :: psi, psi_ex
    type(adaptive_adams_lawson) :: time_stepper

    real(kind=prec) :: err, t0, tend, tol, time_old
    integer :: p, step

    m = fourier1d(256, -16.0_prec, +16.0_prec)
    call m%set_B(B)
    
    t0 = 0.0_prec
    tend = 1.0_prec

    psi_ex = wf_fourier1d(m)
    call psi_ex%set(soliton, tend)


    psi = wf_fourier1d(m)
    call psi%set(soliton, t0)

!    call psi%eval_B(psi_ex)
!    print *, "time=", psi%time
!    print *, "soliton(-.125,0)=", soliton(-.125_prec, 0.0_prec)
!    print *, "B(..)=", B((1.9239578542199571_prec,0.24175518693680464_prec), -.125_prec, .0_prec)
!    call psi_ex%print

    p = 4
    tol = 1e-8_prec

    time_stepper = adaptive_adams_lawson(psi, t0, tend, tol, p)
    step = 0
    time_old = psi%time
    do while (.not. time_stepper%done() )
         print *, "step=",step," t=", psi%time, " dt=", psi%time-time_old
         
         time_old = psi%time
         step = step + 1
       call time_stepper%next
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

        real(prec), parameter :: alpha = 0.0_prec

        complex(kind=prec) :: u_ex

        u_ex = soliton(x, t)

        B = -(1.0_prec-alpha)*abs(u_ex)**2*u -alpha*abs(u)**2*u
    end function B

end program test_soliton
