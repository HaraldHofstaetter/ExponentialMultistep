program test_gaussian
    use wavefunctions_fourier1d
    implicit none
!
! Solve the free Schr\"{o}dinger equation
!   $\partial_t\psi(x,t) = A\psi(x,t)$ where $A=-\frac{1}{2}\Delta$
!
! The initial value is such that the solution is a moving and 
! spreading Gaussian wave packet.
!
    type(fourier1d) :: m
    type(wf_fourier1d) :: psi, psi_ex
    real(kind=prec) :: err, t0, tend

  ! define the spectral method on a grid of 256 points on
  ! the computational domain [-32, +32]
    m = fourier1d(256, -32.0_prec, +32.0_prec)

    t0 = 0.0_prec
    tend = 5.0_prec

  ! allocate a wavefunction psi for the spectral method m
    psi = wf_fourier1d(m)

  ! define the initial value for psi at t=t0 by the function gaussian   
    call psi%set(gaussian, t0) 

  ! propagate psi from t=t0 to t=tend
    call psi%propagate_A(tend-t0)

  ! allocate a wavefunction psi_ref for the reference solution
    psi_ref = wf_fourier1d(m)
  
  ! define the reference solution at t=tend by the function gaussian
    call psi_ref%set(gaussian, tend)

  ! compute the error of the (numerically) propagated solution
  ! compared to the exact solution
    err = psi%distance(psi_ref)

  ! Should be small (if not something went wrong :( )
    print *, "err = ", err

  ! clean up
    call psi%finalize
    call psi_ex%finalize
    call m%finalize

contains
    function gaussian(x, t)
        complex(kind=prec) :: gaussian
        real(kind=prec), intent(in) :: x
        real(kind=prec), intent(in) :: t

        real(kind=prec), parameter :: hbar = 1_prec
        real(kind=prec), parameter :: mass = 1_prec
        real(kind=prec), parameter :: a  = 3.0_prec
        real(kind=prec), parameter :: kx = +0.9_prec
        real(kind=prec), parameter :: pi = 4.0_prec*atan(1.0_prec)

        real(prec) :: phi, theta, a1, cx
        complex(prec) :: b
        !cf. C. Cohen-Tannoudji, B. Diu, F. Laloe: Quantenmechanik 2.Aufl., eq. (1.165)
        theta = 0.5_prec*atan(2.0_prec*hbar*t/(mass*a**2))
        phi = -theta - 0.5_prec*hbar*(kx**2)*t/mass
        a1 = (2.0_prec*a**2/pi/(a**4+4.0_prec*(hbar*t/mass)**2))**0.25_prec
        b = cmplx(a**2, 2.0_prec*hbar*t/mass, kind=prec)
        cx = hbar*kx*t/mass
        gaussian = a1*exp(cmplx(0_prec, phi + kx*x, kind=prec) - ((x-cx)**2)/b)  
    end function gaussian
end program test_gaussian
