module exponential_multistep
! Usage:
!   time_stepper = adaptive_adams_lawson_time_stepper(psi, t0, tend, tol, dt0, p)
!   do while (.not. time_stepper%done() )
!
!        ! do something ...
!
!      call time_stepper%next
!
!         ! do something ...
!
!   end do
 
    use wavefunctions_fourier1d

    implicit none 
    
    type adaptive_adams_lawson_time_stepper 
        private
        
        real(kind=prec) :: t
        type(wf_fourier1d), pointer :: psi !>>>
        real(kind=prec) :: t0
        real(kind=prec) :: tend
        real(kind=prec) :: tol
        real(kind=prec) :: dt

        integer :: p
        integer :: ptr
        real(kind=prec), allocatable :: t_back(:)
        type(wf_fourier1d), allocatable :: rhs_back(:) !>>>
        type(wf_fourier1d) :: psi0 !>>> 
        type(wf_fourier1d) :: psi1 !>>>

        real(kind=prec), allocatable :: tt(:)
        real(kind=prec), allocatable :: c(:)

        logical :: bootstrap_mode
        integer :: p_final

    contains
        procedure :: done => done_adams_lawson
        procedure :: next => next_adams_lawson
    end type adaptive_adams_lawson_time_stepper

    interface adaptive_adams_lawson_time_stepper ! constructor
        module procedure new_adaptive_adams_lawson_time_stepper
    end interface adaptive_adams_lawson_time_stepper



    type adaptive_adams_exponential_time_stepper 
        private
        
        real(kind=prec) :: t
        type(wf_fourier1d), pointer :: psi !>>>
        real(kind=prec) :: t0
        real(kind=prec) :: tend
        real(kind=prec) :: tol
        real(kind=prec) :: dt

        integer :: p
        integer :: ptr
        real(kind=prec), allocatable :: t_back(:)
        type(wf_fourier1d), allocatable :: rhs_back(:) !>>>
        type(wf_fourier1d) :: psi0 !>>> 
        type(wf_fourier1d) :: psi1 !>>>
        type(wf_fourier1d) :: acc  !>>>

        real(kind=prec), allocatable :: tt(:)
        real(kind=prec), allocatable :: c(:,:)

        logical :: bootstrap_mode
        integer :: p_final

    contains
        procedure :: done => done_adams_exponential
        procedure :: next => next_adams_exponential
    end type adaptive_adams_exponential_time_stepper

    interface adaptive_adams_exponential_time_stepper ! constructor
        module procedure new_adaptive_adams_exponential_time_stepper
    end interface adaptive_adams_exponential_time_stepper




contains
    function new_adaptive_adams_lawson_time_stepper(psi, t0, tend, tol, p, dt0) result(this)
        implicit none
        type(adaptive_adams_lawson_time_stepper) :: this
        type(wf_fourier1d), target, intent(in) :: psi
        real(kind=prec), intent(in) :: t0
        real(kind=prec), intent(in) :: tend
        real(kind=prec), intent(in) :: tol
        integer, intent(in), optional :: p
        real(kind=prec), intent(in), optional :: dt0

        integer :: j

        this%psi => psi !>>>
        this%t = t0
        this%tend = tend
        this%tol = tol 

        if (present(dt0)) then
            this%dt = dt0
        else
            this%dt = tol
        end if
        if (present(p)) then
            this%p_final = p
        else
            this%p_final = 4 !default value
        end if

        this%p  = 1
        this%bootstrap_mode = .true.

        allocate( this%t_back(this%p_final+1) )
        this%t_back = 0.0_prec
        allocate( this%rhs_back(this%p_final+1) ) !>>>
        do j =1,this%p_final+1
            this%rhs_back(j) = wf_fourier1d(psi%m) !>>>
        end do

        !>>> rhs_back(1) = B(psi)
        call psi%eval_B(this%rhs_back(1))
        
        this%psi0 = wf_fourier1d(psi%m) !>>>
        this%psi1 = wf_fourier1d(psi%m) !>>>

        allocate( this%tt(this%p_final+1) )
        allocate( this%c(this%p_final+1) )

        this%ptr = this%p
    end function new_adaptive_adams_lawson_time_stepper


    function done_adams_lawson(this) result(f)
        implicit none
        class(adaptive_adams_lawson_time_stepper), intent(inout) :: this
        logical :: f

        integer :: j

        f = (this%t .ge. this%tend)
        if (f) then
            ! destructor
            deallocate( this%t_back )
            do j = 1, this%p_final + 1
                call this%rhs_back(j)%finalize !>>>
            end do
            deallocate( this%rhs_back ) !>>>
            call this%psi0%finalize !>>>
            call this%psi1%finalize !>>>
            deallocate( this%tt )
            deallocate( this%c )
        end if
    end function done_adams_lawson


    subroutine next_adams_lawson(this)
        implicit none
        class(adaptive_adams_lawson_time_stepper), intent(inout) :: this

      ! parameters for stepsize selection    
        real(kind=prec), parameter :: facmin = 0.25_prec
        real(kind=prec), parameter :: facmax = 4.0_prec
        real(kind=prec), parameter :: fac    = 0.9_prec

        real(kind=prec) :: dt0, err, h
        integer :: ptr0
        integer :: j, k, k1 

        call this%psi0%copy(this%psi) !>>> psi0 = psi
        call this%psi1%copy(this%psi) !>>> psi1 = psi

        ptr0 = this%ptr
        err = 2.0_prec

        do while(err .ge. 1.0_prec)
            this%dt = min(this%dt, this%tend-this%t)
            dt0 = this%dt

          ! *** predictor

          ! recompute coefficients
            do k= 1, this%p
                this%tt(k) = this%t_back(mod(k-this%p+this%ptr-1 +this%p+1, this%p+1)+1)
            end do    
            do k= 1, this%p
                this%tt(k) = this%tt(k)-this%tt(this%p)
            end do
            do k= 0, this%p-1
                this%c(k+1) =  this%dt**k/(k+1.0_prec)
            end do
            call solve_vander(this%tt, this%c, this%p)  

            do k= 1, this%p
                k1 = mod(k-this%p+this%ptr-1 +this%p+1, this%p+1) + 1
                !>>> psi1 = psi1 + (c(k)*dt)*rhs_back(k1)
                call this%psi1%axpy(this%rhs_back(k1), cmplx(this%c(k)*this%dt, kind=prec))
            end do
            !>>> psi1 = exp(dt*A)*psi1
            call this%psi1%propagate_A(this%dt) 

            if (this%bootstrap_mode) then
                this%ptr = this%ptr + 1
            else
                this%ptr = mod(this%ptr, this%p+1) + 1
            end if

            !>>> rhs_back(ptr) = B(psi1)
            call this%psi1%eval_B(this%rhs_back(this%ptr)) 
            this%t_back(this%ptr) = this%t+this%dt

            call this%rhs_back(this%ptr)%propagate_A(-this%dt)
          ! *** corrector

          ! recompute coefficients (one more than for predictor)
            do k = 1, this%p+1
                this%tt(k) = this%t_back(mod(k-this%p+this%ptr-2 +this%p+1, this%p+1)+1)
            end do    
            h = this%tt(this%p)
            do k = 1, this%p+1
                this%tt(k) = this%tt(k)-h
            end do
            do k = 0, this%p
                this%c(k+1) =  this%dt**k/(k+1.0_prec)
            end do
            call solve_vander(this%tt, this%c, this%p+1) 

            do k = 1, this%p+1
                k1 = mod(k-this%p+this%ptr-2 +this%p+1, this%p+1) + 1
                !>>> psi = psi + (c(k)*dt)*rhs_back(k1)
                call this%psi%axpy(this%rhs_back(k1), cmplx(this%c(k)*this%dt, kind=prec))
            end do
            !>>> psi = exp(dt*A)*psi
            call this%psi%propagate_A(this%dt)   

          ! compute error estimate and new stepsize
            !>>> err=||psi-psi1||/tol
            err = this%psi%distance(this%psi1)/this%tol
          ! compute new stepize using eq. (II.4.13) in   
          ! Hairer/Norsett/Wanner, Solving Ordinary Differential Equations I, 2nd ed.
            this%dt = this%dt*min(facmax, max(facmin, fac*(1.0_prec/err)**(1.0_prec/(real(this%p, kind=prec)+1_prec))))

            if (err .ge. 1.0_prec) then ! reject step
                call this%psi1%copy(this%psi0) !>>> psi1 = psi0
                call this%psi%copy(this%psi0)  !>>> psi  = psi0
                this%ptr = ptr0
                print *, "t=", this%t, " err=", err, " dt=", this%dt, " rejected..."
            end if
        end do

        !>>> rhs_back(ptr) = B(psi)
        call this%psi%eval_B(this%rhs_back(this%ptr)) 
        this%t = this%t + dt0
        this%t_back(this%ptr) = this%t
    
        do k = 1, this%p+1
            if (k .ne. this%ptr) then 
                !>> rhs_back(k) = exp(dt0*A)*rhs_back(k)
                call this%rhs_back(k)%propagate_A(dt0)
            end if
        end do  
        if (this%bootstrap_mode) then
            if (this%p_final .gt. this%p) then
                this%p = this%p+1
            else
                this%bootstrap_mode = .false.
            end if
        end if
    end subroutine next_adams_lawson



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    function new_adaptive_adams_exponential_time_stepper(psi, t0, tend, tol, p, dt0) result(this)
        implicit none
        type(adaptive_adams_exponential_time_stepper) :: this
        type(wf_fourier1d), target, intent(in) :: psi
        real(kind=prec), intent(in) :: t0
        real(kind=prec), intent(in) :: tend
        real(kind=prec), intent(in) :: tol
        integer, intent(in), optional :: p
        real(kind=prec), intent(in), optional :: dt0

        integer :: j

        this%psi => psi !>>>
        this%t = t0
        this%tend = tend
        this%tol = tol 

        if (present(dt0)) then
            this%dt = dt0
        else
            this%dt = tol
        end if
        if (present(p)) then
            this%p_final = p
        else
            this%p_final = 4 !default value
        end if

        this%p  = 1
        this%bootstrap_mode = .true.

        allocate( this%t_back(this%p_final+1) )
        this%t_back = 0.0_prec
        allocate( this%rhs_back(this%p_final+1) ) !>>>
        do j =1,this%p_final+1
            this%rhs_back(j) = wf_fourier1d(psi%m) !>>>
        end do

        !>>> rhs_back(1) = B(psi)
        call psi%eval_B(this%rhs_back(1))
        
        this%psi0 = wf_fourier1d(psi%m) !>>>
        this%psi1 = wf_fourier1d(psi%m) !>>>
        this%acc =  wf_fourier1d(psi%m) !>>>

        allocate( this%tt(this%p_final+1) )
        allocate( this%c(this%p_final+1, this%p_final+1) )

        this%ptr = this%p
    end function new_adaptive_adams_exponential_time_stepper


    function done_adams_exponential(this) result(f)
        implicit none
        class(adaptive_adams_exponential_time_stepper), intent(inout) :: this
        logical :: f

        integer :: j

        f = (this%t .ge. this%tend)
        if (f) then
            ! destructor
            deallocate( this%t_back )
            do j = 1, this%p_final + 1
                call this%rhs_back(j)%finalize !>>>
            end do
            deallocate( this%rhs_back ) !>>>
            call this%psi0%finalize !>>>
            call this%psi1%finalize !>>>
            call this%acc%finalize  !>>>
            deallocate( this%tt )
            deallocate( this%c )
        end if
    end function done_adams_exponential


    subroutine next_adams_exponential(this)
        implicit none
        class(adaptive_adams_exponential_time_stepper), intent(inout) :: this

      ! parameters for stepsize selection    
        real(kind=prec), parameter :: facmin = 0.25_prec
        real(kind=prec), parameter :: facmax = 4.0_prec
        real(kind=prec), parameter :: fac    = 0.9_prec

        real(kind=prec) :: dt0, err, h
        integer :: ptr0
        integer :: j, k, n, k1 

        call this%psi0%copy(this%psi) !>>> psi0 = psi
        call this%psi1%copy(this%psi) !>>> psi1 = psi

        ptr0 = this%ptr
        err = 2.0_prec

        do while(err .ge. 1.0_prec)
            this%dt = min(this%dt, this%tend-this%t)
            dt0 = this%dt

          ! *** predictor

          ! recompute coefficients
            do k = 1, this%p
                this%tt(k) = this%t_back(mod(k-this%p+this%ptr-1 +this%p+1, this%p+1)+1)
            end do    
            do k = 1, this%p
                this%tt(k) = this%tt(k)-this%tt(this%p)
            end do
            call adams_multistep_coefficients(this%tt, this%c, this%dt, this%p_final+1, this%p)

            !>>> psi1 = exp(dt*A)*psi1
            call this%psi1%propagate_A(this%dt) 

            do n = 1, this%p
                this%acc%u = 0.0_prec !>>> acc=0
                do k = 1, this%p
                    k1 = mod(k-this%p+this%ptr-1 +this%p+1, this%p+1) + 1
                    !>>> acc = acc + c(k,n)*rhs_back(k1)
                    call this%acc%axpy(this%rhs_back(k1), cmplx(this%c(k,n), kind=prec))
                end do
                !>>> psi1 = psi1 + dt*phi_n(dt*A)*acc
                call this%acc%add_phi_A(this%psi1, this%dt, n, this%dt)
            end do


            if (this%bootstrap_mode) then
                this%ptr = this%ptr + 1
            else
                this%ptr = mod(this%ptr, this%p+1) + 1
            end if

            !>>> rhs_back(ptr) = B(psi1)
            call this%psi1%eval_B(this%rhs_back(this%ptr)) 
            this%t_back(this%ptr) = this%t+this%dt

          ! *** corrector

          ! recompute coefficients (one more than for predictor)
            do k = 1, this%p+1
                this%tt(k) = this%t_back(mod(k-this%p+this%ptr-2 +this%p+1, this%p+1)+1)
            end do    
            h = this%tt(this%p)
            do k = 1, this%p+1
                this%tt(k) = this%tt(k)-h
            end do
            call adams_multistep_coefficients(this%tt, this%c, this%dt, this%p_final+1, this%p+1)

            !>>> psi = exp(dt*A)*psi
            call this%psi%propagate_A(this%dt) 

            do n = 1, this%p+1
                this%acc%u = 0.0_prec !>>> acc=0
                do k = 1, this%p+1
                    k1 = mod(k-this%p+this%ptr-2 +this%p+1, this%p+1) + 1
                    !>>> acc = acc + c(k,n)*rhs_back(k1)
                    call this%acc%axpy(this%rhs_back(k1), cmplx(this%c(k,n), kind=prec))
                end do
                !>>> psi = psi + dt*phi_n(dt*A)*acc
                call this%acc%add_phi_A(this%psi, this%dt, n, this%dt)
            end do
            
          ! compute error estimate and new stepsize
            !>>> err=||psi-psi1||/tol
            err = this%psi%distance(this%psi1)/this%tol
          ! compute new stepize using eq. (II.4.13) in   
          ! Hairer/Norsett/Wanner, Solving Ordinary Differential Equations I, 2nd ed.
            this%dt = this%dt*min(facmax, max(facmin, fac*(1.0_prec/err)**(1.0_prec/(real(this%p, kind=prec)+1_prec))))

            if (err .ge. 1.0_prec) then ! reject step
                call this%psi1%copy(this%psi0) !>>> psi1 = psi0
                call this%psi%copy(this%psi0)  !>>> psi  = psi0
                this%ptr = ptr0
                print *, "t=", this%t, " err=", err, " dt=", this%dt, " rejected..."
            end if
        end do

        !>>> rhs_back(ptr) = B(psi)
        call this%psi%eval_B(this%rhs_back(this%ptr)) 
        this%t = this%t + dt0
        this%t_back(this%ptr) = this%t
    
        if (this%bootstrap_mode) then
            if (this%p_final .gt. this%p) then
                this%p = this%p+1
            else
                this%bootstrap_mode = .false.
            end if
        end if
    end subroutine next_adams_exponential



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    subroutine solve_vander(x, b, n)
        integer, intent(in) :: n
        real(kind=prec), intent(in) :: x(n)
        real(kind=prec), intent(inout) :: b(n)

      ! Algorithm 4.6.2 from Golub/van Loan, Matrix Computations.
      ! Solve Vz=b, where V is the Vandermonde matrix with respect 
      ! to x=(x(1),...,x(n)); b is overwritten by the solution z.

        integer i, k

        do k = 1, n-1
            do i = n, k+1, -1
                b(i) = b(i) - x(k)*b(i-1)
            end do
        end do
        do k = n-1, 1, -1
            do i = k+1, n
                b(i) = b(i)/(x(i)-x(i-k))
            end do
            do i = k, n-1
                b(i) = b(i) - b(i+1)
            end do
       end do
    end subroutine solve_vander


    subroutine inverse_vander(x, Vi, ld, n)
        integer, intent(in) :: ld 
        integer, intent(in) :: n
        real(kind=prec), intent(in) :: x(n)
        real(kind=prec), intent(inout) :: Vi(ld, *)

        integer :: i, j

        do i = 1, n
            do j=1,n
                Vi(j, i) = 0.0_prec
                end do 
            Vi(i, i) = 1.0_prec
            call solve_vander(x, Vi(1,i), n)
        end do
    end subroutine inverse_vander


    subroutine adams_multistep_coefficients(x, C, dt, ld, n)
        integer, intent(in) :: ld 
        integer, intent(in) :: n
        real(kind=prec), intent(in) :: x(n)
        real(kind=prec), intent(inout) :: C(ld, *)
        real(kind=prec), intent(in) :: dt
        
        real(kind=prec) :: f
        integer :: i,j
   
        call inverse_vander(x, C, ld, n)

        !multiply from the left by diag(dt^i*factorial(i))
        f = 1.0_prec
        do i = 1, n
            do j=1, n
                C(j,i) = C(j,i)*f
            end do
            f = f*dt*i
        end do
    end subroutine adams_multistep_coefficients
    
end module exponential_multistep
