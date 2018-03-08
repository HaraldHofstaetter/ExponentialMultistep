program test_crossing_solitons
    use wavefunctions_fourier1d
    use exponential_multistep
    implicit none
!    
! Solve the nonlinear Schr\"{o}dinger equation  
!    $$i\partial_t\psi(x,t) = A\psi(x,t) + B(\psi(x,t),x,t))$$    
! where $A=-\frac{1}{2}\Delta$ and $B(u,x,t) = -1/2|u|^2\u$.
! is defined by the Fortran function B(u,x,t) below.
!
! The initial value defined by the Fortran function psi_0(x,t) 
! below gives a solution  with two crossing solitons.
!
    
    type(fourier1d) :: m
    type(wf_fourier1d) :: psi
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
    tend = 5.0_prec

  ! allocate a wavefunction for the spectral method m
    psi = wf_fourier1d(m)
   
  ! define the initial value for psi at t=t0 by the function psi_0
    call psi%set(psi_0, t0)

  ! define the order and tolerance for the adaptive time stepper
    p = 4            ! order
    tol = 1e-5_prec  ! tolerance

  ! run the adaptive time stepper  
    step = 0
    time_old = psi%time
    time_stepper = adaptive_adams_lawson_time_stepper(psi, t0, tend, tol, p)
    !time_stepper = adaptive_adams_exponential_time_stepper(psi, t0, tend, tol, p)
    do while (.not. time_stepper%done() )
         call time_stepper%next
       ! print number of step, current time, and current stepsize  
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


end program test_crossing_solitons
