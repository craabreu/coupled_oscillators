module mThermostat

implicit none

type :: NoseHooverChain
    integer :: nchain                ! Number of thermostats in the chain
    integer :: nloops                ! Number of integration loops
    real(8) :: kT                    ! Set-point temperature
    real(8) :: invQ_eta
    real(8), allocatable :: v_eta(:,:)
  contains
    procedure :: integrate_v_eta => NoseHooverChain_integrate_v_eta
    procedure :: allocate => NoseHooverChain_allocate
    procedure :: initialize => NoseHooverChain_initialize
    procedure :: transform => NoseHooverChain_transform
    procedure :: force => NoseHooverChain_force
    procedure :: integrate => NoseHooverChain_integrate
end type

type, extends(NoseHooverChain) :: IsokineticNoseHooverChain
    real(8) :: factor
    real(8) :: sqrt_n
  contains
    procedure :: initialize => IsokineticNoseHooverChain_initialize
    procedure :: transform => IsokineticNoseHooverChain_transform
    procedure :: force => IsokineticNoseHooverChain_force
end type

type, extends(NoseHooverChain) :: SubkineticNoseHooverChain
    real(8) :: sqrt_n
  contains
    procedure :: initialize => SubkineticNoseHooverChain_initialize
    procedure :: transform => SubkineticNoseHooverChain_transform
    procedure :: force => SubkineticNoseHooverChain_force
end type

contains

  pure subroutine NoseHooverChain_allocate(self, Nf, nchain, nloops, kT, Q_eta)
    class(NoseHooverChain), intent(inout)        :: self
    integer,                intent(in)           :: Nf, nchain, nloops
    real(8),                intent(in)           :: kT, Q_eta
    self%nchain = nchain
    self%nloops = nloops
    self%kT = kT
    self%invQ_eta = 1.d0/Q_eta
    allocate( self%v_eta(Nf, nchain), source = 0d0 )
  end subroutine NoseHooverChain_allocate

  pure subroutine NoseHooverChain_initialize(self, Nf, nchain, nloops, kT, Q_eta, n)
    class(NoseHooverChain), intent(inout)        :: self
    integer,                intent(in)           :: Nf, nchain, nloops
    real(8),                intent(in)           :: kT, Q_eta
    real(8),                intent(in), optional :: n
    call self%allocate( Nf, nchain, nloops, kT, Q_eta )
  end subroutine NoseHooverChain_initialize

  pure subroutine NoseHooverChain_integrate_v_eta(self, v_eta, force, damp, dt)
    class(NoseHooverChain), intent(in)    :: self
    real(8),                intent(inout) :: v_eta(:)
    real(8),                intent(in)    :: force(:), damp(:), dt
    v_eta = v_eta + 0.5d0*dt*force*self%invQ_eta
    v_eta = v_eta*exp(-damp*dt)
    v_eta = v_eta + 0.5d0*dt*force*self%invQ_eta
  end subroutine NoseHooverChain_integrate_v_eta

  pure subroutine NoseHooverChain_transform(self, p, dt)
    class(NoseHooverChain), intent(inout) :: self
    real(8),                intent(inout) :: p(:)
    real(8),                intent(in)    :: dt
    p = p*exp(-self%v_eta(:,1)*dt)
  end subroutine NoseHooverChain_transform

  pure function NoseHooverChain_force(self, p) result( force )
    class(NoseHooverChain), intent(in) :: self
    real(8),                intent(in) :: p(:)
    real(8)                            :: force(size(p))
    force = p*p - self%kT
  end function NoseHooverChain_force

  pure subroutine NoseHooverChain_integrate(self, p, timestep)
    class(NoseHooverChain), intent(inout) :: self
    real(8),                intent(inout) :: p(:)
    real(8),                intent(in)    :: timestep

    integer :: loop, j
    real(8) :: dt, half_dt

    associate(v => self%v_eta, n => self%nchain)
      dt = timestep/self%nloops
      half_dt = 0.5d0*dt
      do loop = 1, self%nloops
        v(:,n) = v(:,n) + half_dt*(v(:,n-1)**2 - self%kT)*self%invQ_eta
        if (n > 1) then
          do j = n - 1, 2, -1
            call self%integrate_v_eta(v(:,j), v(:,j-1)**2 - self%kT, v(:,j+1), half_dt)
          end do
          call self%integrate_v_eta(v(:,1), self%force(p), v(:,2), half_dt)
        end if
        call self%transform(p, dt)
        if (n > 1) then
          call self%integrate_v_eta(v(:,1), self%force(p), v(:,2), half_dt)
          do j = 2, n - 1
            call self%integrate_v_eta(v(:,j), v(:,j-1)**2 - self%kT, v(:,j+1), half_dt)
          end do
        end if
        v(:,n) = v(:,n) + half_dt*(v(:,n-1)**2 - self%kT)*self%invQ_eta
      end do
    end associate
  end subroutine NoseHooverChain_integrate

  pure subroutine IsokineticNoseHooverChain_initialize(self, Nf, nchain, nloops, kT, Q_eta, n)
    class(IsokineticNoseHooverChain), intent(inout)        :: self
    integer,                          intent(in)           :: Nf, nchain, nloops
    real(8),                          intent(in)           :: kT, Q_eta
    real(8),                          intent(in), optional :: n
    real(8) :: nn
    call self%allocate( Nf, nchain, nloops, kT, Q_eta )
    if (present(n)) then
      self%factor = (n + 1d0)/n
      self%sqrt_n = sqrt(n)
    else
      self%factor = 2d0
      self%sqrt_n = 1.d0
    end if
  end subroutine IsokineticNoseHooverChain_initialize

  pure subroutine IsokineticNoseHooverChain_transform(self, p, dt)
    class(IsokineticNoseHooverChain), intent(inout) :: self
    real(8),                          intent(inout) :: p(:)
    real(8),                          intent(in)    :: dt
    p = self%sqrt_n*asinh(exp(-self%v_eta(:,1)*dt)*sinh(p/self%sqrt_n))
  end subroutine IsokineticNoseHooverChain_transform

  pure function IsokineticNoseHooverChain_force(self, p) result( force )
    class(IsokineticNoseHooverChain), intent(in) :: self
    real(8),                          intent(in) :: p(:)
    real(8)                                      :: force(size(p))
    real(8) :: v(size(p))
    v = self%sqrt_n*tanh(p/self%sqrt_n)
    force = self%factor*v*v - self%kT
  end function IsokineticNoseHooverChain_force

  pure subroutine SubkineticNoseHooverChain_initialize(self, Nf, nchain, nloops, kT, Q_eta, n)
    class(SubkineticNoseHooverChain), intent(inout)        :: self
    integer,                          intent(in)           :: Nf, nchain, nloops
    real(8),                          intent(in)           :: kT, Q_eta
    real(8),                          intent(in), optional :: n
    real(8) :: nn
    call self%allocate( Nf, nchain, nloops, kT, Q_eta )
    if (present(n)) then
      self%sqrt_n = sqrt(n)
    else
      self%sqrt_n = 1.d0
    end if
  end subroutine SubkineticNoseHooverChain_initialize

  pure subroutine SubkineticNoseHooverChain_transform(self, p, dt)
    class(SubkineticNoseHooverChain), intent(inout) :: self
    real(8),                          intent(inout) :: p(:)
    real(8),                          intent(in)    :: dt
    p = exp(-self%v_eta(:,1)*dt)*p
  end subroutine SubkineticNoseHooverChain_transform

  pure function SubkineticNoseHooverChain_force(self, p) result( force )
    class(SubkineticNoseHooverChain), intent(in) :: self
    real(8),                          intent(in) :: p(:)
    real(8)                                      :: force(size(p))
    force = self%sqrt_n*p*tanh(p/self%sqrt_n) - self%kT
  end function SubkineticNoseHooverChain_force

end module mThermostat
