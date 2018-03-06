module wavefunctions_fourier1d
    use, intrinsic :: iso_c_binding

    implicit none

  ! FFTW stuff
    include 'fftw3.f03'
    integer(C_INT) :: fftw_planning_rigor = FFTW_PATIENT
    character(len=64, kind=c_char) :: fftw_wisdom_file = c_char_'fftw_wisdom' // c_null_char    
    
    integer, parameter :: prec=selected_real_kind(p=15)

    public :: fourier1d, wf_fourier1d


    type fourier1d
        real(kind=prec) :: xmin
        real(kind=prec) :: xmax
        integer         :: nx

        real(kind=prec) :: dx
        real(kind=prec), allocatable :: nodes_x(:)
        real(kind=prec), allocatable :: eigvals_A(:)
    contains
        procedure :: finalize => finalize_method
    end type fourier1d
   
    interface  fourier1d ! constructor
        module procedure new_method
    end interface fourier1d


    type wf_fourier1d
      ! high level data
        class(fourier1d), pointer   :: m 
        real(kind=prec)             :: time = 0.0_prec
        complex(kind=prec), pointer :: u(:)

      ! lower level data
        type(c_ptr)  :: plan_forward 
        type(c_ptr)  :: plan_backward        
    contains
        procedure :: set
        procedure :: print
        procedure :: propagate_A
        procedure :: copy
        procedure :: scale
        procedure :: axpy
        procedure :: norm
        procedure :: distance
        procedure :: finalize => finalize_wf
    end type wf_fourier1d

    interface wf_fourier1d ! constructor
        module procedure new_wf
    end interface wf_fourier1d

contains

    function new_method(nx, xmin, xmax) result(this)
        implicit none
        type(fourier1d)          :: this
        real(kind=prec), intent(in) :: xmin 
        real(kind=prec), intent(in) :: xmax
        integer, intent(in)      :: nx

        integer :: ix, k, m
        real(kind=prec) :: d
        logical, save :: fftw_already_initialized = .false.
        integer(c_int) :: ret
        real(kind=prec), parameter :: pi = 4.0_prec*atan(1.0_prec)

        this%xmin = xmin
        this%xmax = xmax
        this%nx = nx
        this%dx = (xmax-xmin)/nx

        allocate( this%nodes_x(1:nx) )
            this%nodes_x = (/ ( xmin+this%dx*real(ix, prec), ix = 1, nx )  /)

        allocate( this%eigvals_A(1:nx) )
            m = nx/2
            d = xmax-xmin
            this%eigvals_A(1:m)  = (/ (real(k, kind=prec)**2, k=1-1, m-1) /)
            this%eigvals_A(m:nx) = (/ (real(k, kind=prec)**2, k=m-nx-1, -1) /)
            this%eigvals_A = (-(2.0_prec*pi/d)**2/2.0_prec) * this%eigvals_A

      ! initialize if not already been done
        if (.not. fftw_already_initialized) then
            ret = fftw_import_wisdom_from_filename(fftw_wisdom_file)
            fftw_already_initialized = .true.
        end if
    end function new_method


    subroutine finalize_method(this)
        class(fourier1d), intent(inout) :: this

        deallocate( this%nodes_x )
        deallocate( this%eigvals_A )
    end subroutine finalize_method


    function new_wf(m) result(this) ! constructor for wave functions
        implicit none
        type(wf_fourier1d) :: this
        class(fourier1d), target, intent(inout) :: m
        
        type(c_ptr) :: p
        integer(c_int) :: ret

        this%m => m
        this%time = 0.0_prec

      ! let allocation of storage for wavefunction
      ! be done by FFTW (because of alignment issues)
        p = fftw_alloc_complex(int(this%m%nx, kind=C_SIZE_T))        
        call c_f_pointer(p, this%u, [this%m%nx])

      ! create plans
        this%plan_forward  = fftw_plan_dft_1d(this%m%nx, this%u, this%u, FFTW_FORWARD, fftw_planning_rigor)
        this%plan_backward = fftw_plan_dft_1d(this%m%nx, this%u, this%u, FFTW_BACKWARD, fftw_planning_rigor)

      ! save plans  
        ret = fftw_export_wisdom_to_filename(fftw_wisdom_file)
    end

    subroutine finalize_wf(this)
        use, intrinsic :: iso_c_binding, only: c_loc, c_associated
        class(wf_fourier1d), intent(inout) :: this

        call fftw_destroy_plan(this%plan_forward)
        call fftw_destroy_plan(this%plan_backward)
        call fftw_free(c_loc(this%u(1)))
    end subroutine finalize_wf


    subroutine set(this, f, t)
        implicit none
        class(wf_fourier1d) :: this
        complex(kind=prec), external :: f
        real(kind=prec), intent(in), optional :: t

        integer :: ix
        real(kind=prec) :: x
        if(present(t)) then
            this%time = t
        else
            this%time = 0.0_prec 
        end if
        do ix = 1,this%m%nx
            this%u(ix) = f(this%m%nodes_x(ix), this%time)
        end do
    end subroutine set


    subroutine print(this)
        implicit none
        class(wf_fourier1d) :: this
        integer :: ix

        do ix = 1,this%m%nx
            print *, this%m%nodes_x(ix), real(this%u(ix)), aimag(this%u(ix))
        end do
    end subroutine print


    subroutine propagate_A(this, dt) ! this = exp(dt*A)*this (where A = -i*Laplacian)
        implicit none
        class(wf_fourier1d), intent(inout) :: this
        real(kind=prec), intent(in) :: dt

      ! transform to frequency space 
        call fftw_execute_dft(this%plan_forward, this%u, this%u) 
      
      ! apply exponentiated operator (diagonal in frequency space)
        this%u = exp(cmplx(0.0_prec, dt, kind=prec)*this%m%eigvals_A) * this%u

      ! back-transform to real space 
        call fftw_execute_dft(this%plan_backward, this%u, this%u)

      ! scale by 1/nx (not already done by fftw)  
        this%u = cmplx(1.0_prec/this%m%nx, kind=prec)*this%u

      ! propagate time:
        this%time = this%time + dt        
    end subroutine propagate_A


    subroutine scale(this, fac) ! this = fac*this
        implicit none
        class(wf_fourier1d), intent(inout) :: this
        complex(kind=prec), intent(in) :: fac

        this%u = fac*this%u
    end subroutine scale


    subroutine copy(this, source) ! this = source
        class(wf_fourier1d), intent(inout) :: this
        class(wf_fourier1d), intent(inout) :: source 
 
        if (.not.associated(source%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    
        this%u = source%u
        this%time = source%time
    end subroutine copy


    subroutine axpy(this, other, fac) ! this = this + fac*other
        class(wf_fourier1d), intent(inout) :: this
        class(wf_fourier1d), intent(inout) :: other 
        complex(kind=prec), intent(in) :: fac
 
        if (.not.associated(other%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    
        this%u = this%u + fac*other%u
    end subroutine axpy

    function norm(this) result(n)
        class(wf_fourier1d), intent(inout) :: this
        real(kind=prec) :: n

        n = norm2([norm2(real(this%u)), norm2(aimag(this%u))])*sqrt(this%m%dx)
    end function norm

    function distance(this, other) result(d)
        class(wf_fourier1d), intent(inout) :: this
        class(wf_fourier1d), intent(inout) :: other 
        real(kind=prec) :: d

        if (.not.associated(other%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    
        d = norm2([norm2(real(this%u)-real(other%u)), norm2(aimag(this%u)-aimag(other%u))])*sqrt(this%m%dx)
    end function distance


end module wavefunctions_fourier1d
