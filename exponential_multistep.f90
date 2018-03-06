module exponential_multistep
!    time_stepper = adaptive_adams_lawson(psi, t0, tend, tol, dt0, p)
!    do while (.not. time_stepper%done() )
!
!         ! do something with psi
!
!       call time_stepper%next
!    end do

    use wavefunctions_fourier1d

    implicit none 
    
    type adaptive_adams_lawson ! iterator
        private
        
        real(kind=prec) :: t
        type(wf_fourier1d), pointer :: psi
        real(kind=prec) :: t0
        real(kind=prec) :: tend
        real(kind=prec) :: tol
        real(kind=prec) :: dt

        integer :: n
        integer :: n1
        integer :: ptr
        real(kind=prec), allocatable :: t_back(:)
        type(wf_fourier1d), allocatable :: rhs_back(:)
        type(wf_fourier1d) :: psi0 
        type(wf_fourier1d) :: psi1 

        real(kind=prec), allocatable :: tt(:)
        real(kind=prec), allocatable :: c(:)

        logical :: bootstrap_mode
        integer :: n1_final

    contains
        procedure :: done
        procedure :: next
    end type adaptive_adams_lawson

    interface adaptive_adams_lawson ! constructor
        module procedure new_adaptive_adams_lawson
    end interface adaptive_adams_lawson

contains
    function new_adaptive_adams_lawson(psi, t0, tend, tol, dt0, p) result(this)
        implicit none
        type(adaptive_adams_lawson) :: this
        type(wf_fourier1d), target, intent(in) :: psi
        real(kind=prec), intent(in) :: t0
        real(kind=prec), intent(in) :: tend
        real(kind=prec), intent(in) :: tol
        real(kind=prec), intent(in), optional :: dt0
        integer, intent(in), optional :: p

        integer :: n_final, j

        this%psi => psi
        this%t = t0
        this%tend = tend
        this%tol = tol 

        if (present(dt0)) then
            this%dt = dt0
        else
            this%dt = tol
        end if
        if (present(p)) then
            this%n1_final = p
        else
            this%n1_final = 4 
        end if

        this%n   = 2 
        this%n1  = 1
        this%bootstrap_mode = .true.
        n_final = this%n1_final + 1

        allocate( this%t_back(n_final) )
        allocate( this%rhs_back(n_final) )
        do j =1,n_final
            this%rhs_back(j) = wf_fourier1d(psi%m)
        end do
        this%psi0 = wf_fourier1d(psi%m)
        this%psi1 = wf_fourier1d(psi%m)

        allocate( this%tt(n_final) )
        allocate( this%c(n_final) )

        this%ptr = this%n1

    end function new_adaptive_adams_lawson

    function done(this) result(f)
        implicit none
        class(adaptive_adams_lawson), intent(inout) :: this
        logical :: f

        integer :: j

        f = (this%t .ge. this%tend)
        if (f) then
            ! destructor
            deallocate( this%t_back )
            do j = 1, this%n1_final + 1
                call this%rhs_back(j)%finalize
            end do
            deallocate( this%rhs_back )
            call this%psi0%finalize
            call this%psi1%finalize
            deallocate( this%tt )
            deallocate( this%c )
        end if
    end function done


    subroutine next(this)
        implicit none
        class(adaptive_adams_lawson), intent(inout) :: this

        real(kind=prec), parameter :: facmin = 0.25_prec
        real(kind=prec), parameter :: facmax = 4.0_prec
        real(kind=prec), parameter :: fac    = 0.9_prec

        real(kind=prec) :: dt0, err
        integer :: ptr0
        integer :: j, k, k1 

        call this%psi0%copy(this%psi) ! psi0 = psi
        call this%psi1%copy(this%psi) ! psi1 = psi

        ptr0 = this%ptr
        err = 2.0_prec

        do while(err .ge. 1.0_prec)
            this%dt = min(this%dt, this%tend-this%t)
            dt0 = this%dt

          ! *** predictor

          ! recompute coefficients
            do k= 1, this%n1
                this%tt(k) = this%t_back(mod(k-this%n1+this%ptr-1, this%n)+1)
            end do    
            do k= 1, this%n1
                this%tt(k) = this%tt(k)-this%tt(this%n1)
            end do
            do k= 0, this%n1-1
                this%c(k+1) =  this%dt**k/(k+1.0_prec)
            end do
            call solve_vander_trans(this%tt, this%c, this%n1)  !! TODO

            do k= 1, this%n1
                k1 = mod(k-this%n1+this%ptr-1, this%n) + 1
                call this%psi1%axpy(this%rhs_back(k1), cmplx(this%c(k)*this%dt, kind=prec))
            end do
            call this%psi1%propagate_A(this%dt)   

            if (this%bootstrap_mode) then
                this%ptr = this%ptr + 1
            else
                this%ptr = mod(this%ptr, this%n) + 1
            end if

            call this%psi1%eval_B(this.rhs_back(this%ptr)) !! TODO
            this%t_back(this%ptr) = this%t+this%dt

            call this%rhs_back(this%ptr)%propagate_A(-this%dt)
          
          ! *** corrector

          ! recompute coefficients (one more than for predictor)
            do k = 1, this%n1+1
                this%tt(k) = this%t_back(mod(k-this%n1+this%ptr-2, this%n)+1)
            end do    
            do k = 1, this%n1+1
                this%tt(k) = this%tt(k)-this%tt(this%n1)
            end do
            do k = 0, this%n1
                this%c(k+1) =  this%dt**k/(k+1.0_prec)
            end do
            call solve_vander_trans(this%tt, this%c, this%n1+1)  !! TODO

            do k = 1, this%n1+1
                k1 = mod(k-this%n1+this%ptr-2, this%n) + 1
                call this%psi%axpy(this%rhs_back(k1), cmplx(this%c(k)*this%dt, kind=prec))
            end do
            call this%psi%propagate_A(this%dt)   

          ! compute error estimate and new stepsize
            err = this%psi%distance(this%psi1)/this%tol
            this%dt = this%dt*min(facmax, max(facmin, fac*(1.0_prec/err)**(1.0_prec/(this%n1+1))))

            if (err .ge. 1.0_prec) then ! reject step
                call this%psi1%copy(this%psi0)
                call this%psi%copy(this%psi0)
                this%ptr = ptr0
                print *, "t=", this%t, " err=", err, " dt=", this%dt, " rejected..."
            end if
        end do

        call this%psi%eval_B(this.rhs_back(this%ptr)) !! TODO
        this%t = this%t + dt0
        this%t_back(this%ptr) = this%t
    
        do k = 1, this%n
            if (k .ne. this%ptr) then 
                call this%rhs_back(k)%propagate_A(dt0)
            end if
        end do  
        if (this%bootstrap_mode) then
            if (this%n1_final .gt. this%n1) then
                this%n1 = this%n1+1
                this%n = this%n+1
            else
                this%bootstrap_mode = .false.
            end if
        end if
    end subroutine next



end module exponential_multistep
