program test_soliton
    use wavefunctions_fourier1d
    use exponential_multistep
    implicit none
!
! Solve the nonlinear Schr\"{o}dinger equation  
!    $$i\partial_t\psi(x,t) = A\psi(x,t) + B(\psi(x,t),x,t))$$    
! where $A=-\frac{1}{2}\Delta$ and $B(u,x,t) = V(x,t)*u -\frac{1}{2}|u|^2\u$ 
! is defined by the Fortran function B(u,x,t) below.
! Here the time-dependent potential is given by $V(x,t) = \frac{1}{2}|\psi_{ex}(x,t)|^2$
! where $psi_{ex}(x,t)=\frac{a e^{\frac{1}{2}i((a^2-b^2)t-bx)}{\cosh(a(bt+x-c))}$
!
! $\psi_{ex}(x,t)$ defined by the Fortran function soliton(x,t) is an exact 
! to this problem.
!

    type(fourier1d) :: m
    type(wf_fourier1d) :: psi, psi_ref
    type(adaptive_adams_lawson_time_stepper) :: time_stepper
    !type(adaptive_adams_exponential_time_stepper) :: time_stepper

    real(kind=prec) :: err, t0, tend, tol, time_old
    integer :: p, step

  ! define the spectral method on a grid of 1024 points on
  ! the computational domain [-16, +16]
    m = fourier1d(1024, -16.0_prec, +16.0_prec)

  ! register the function B    
    call m%set_B(B)
    
    t0 = 0.0_prec
    tend = 1.0_prec

  ! allocate a wavefunction psi for the spectral method m
    psi = wf_fourier1d(m)

  ! define the initial value for psi at t=t0 by the function soliton 
    call psi%set(soliton, t0)

  ! define the order and tolerance for the adaptive time stepper
    p = 4            ! order
    tol = 1e-5_prec  ! tolerance

 ! run the adaptive time stepper
    call psi%set(soliton, t0)
    step = 0
    time_old = psi%time
    time_stepper = adaptive_adams_lawson_time_stepper(psi, t0, tend, tol, p, tol)
    !time_stepper = adaptive_adams_exponential_time_stepper(psi, t0, tend, tol, p, tol)
    do while (.not. time_stepper%done() )
         call time_stepper%next
       ! print number of step, current time, and current stepsize
         print *, "step=",step," t=", psi%time, " dt=", psi%time-time_old
         time_old = psi%time
         step = step + 1
    end do

  ! allocate a wavefunction psi_ref for the reference solution 
    psi_ref = wf_fourier1d(m)

  ! define the reference solution at t=tend by the function soliton
    call psi_ref%set(soliton, tend)


  ! compute the error of the (numerically) propagated solution
  ! compared to the reference solution.
  ! Should be small (if not something went wrong :( )
    err = psi%distance(psi_ref)

    print *, "err = ", err

  ! clean up
    call psi%finalize
    call psi_ref%finalize
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
