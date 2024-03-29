program harmonic

use mRandom
use mThermostat

implicit none

integer, parameter :: nrespa = NRESPA
integer, parameter :: nsamples = 1000000
integer, parameter :: isample = 100
integer, parameter :: nsteps = nsamples*isample
real(8), parameter :: n = NN
real(8), parameter :: gamma = 1d0
real(8), parameter :: L = 10d0
real(8), parameter :: dt = 0.01d0
real(8), parameter :: Q_eta = 1d0

real(8), parameter :: pi = 4d0*datan(1d0)
real(8), parameter :: mksoft = -(2d0*pi/L)**2
real(8), parameter :: kstiff = 2d0*pi**2*(1d0 - 1d0/L**2)
real(8), parameter :: mKfast(2, 2) = reshape([-kstiff, kstiff, kstiff, -kstiff], [2, 2])
real(8), parameter :: K(2, 2) = reshape([-mksoft+kstiff, -kstiff, -kstiff, -mksoft+kstiff], [2, 2])
real(8), parameter :: dt_2 = 0.5d0*dt
real(8), parameter :: O1 = exp(-gamma*dt)
real(8), parameter :: O2 = sqrt(1d0 - O1**2)
real(8), parameter :: sqrt_n = sqrt(n)

integer :: nseeds
integer, allocatable :: seed(:)
real(8) :: q(2), p(2), v_eta(2)
type(xrsr128) :: random
type(NoseHooverChain) :: nhc, xo_nhc
type(IsokineticNoseHooverChain) :: isok_nhc, sinr_nhc
type(SubkineticNoseHooverChain) :: subk_nhc

#define normalvec [random%normal(),random%normal()]
#define F_fast(q) matmul(mKfast, q)
#define F_slow(q) (mksoft*(q))
#define v(p) (sqrt_n*tanh((p)/sqrt_n))

#define SAMPLE if (mod(step, isample) == 0) write(*,'(4F13.8)') q, p

call random_seed(size=nseeds)
allocate( seed(nseeds) )
call random_seed(get=seed)
call random % setup( seed(1) )

call xo_nhc % initialize( 2, 3, nrespa, 1d0, Q_eta, n )
call nhc % initialize( 2, 3, 1, 1d0, Q_eta, n )
call isok_nhc % initialize( 2, 3, 1, 1d0, Q_eta, n )
call sinr_nhc % initialize( 2, 1, 1, 1d0, Q_eta, n )
call subk_nhc % initialize( 2, 3, 1, 1d0, Q_eta, n )

q = normalvec
p = normalvec
v_eta = normalvec

print*, 'q1 q2 p1 p2'
call METHOD(nsteps)

contains
  !---------------------------
  ! subroutine Verlet(nsteps)
  !  integer, intent(in) :: nsteps
  !  integer :: step, substep
  !  do step = 1, nsteps
  !    p = p + nrespa*dt_2*F_slow(q)
  !    do substep = 1, nrespa
  !      p = p + dt_2*F_fast(q)
  !      q = q + dt*p
  !      p = p + dt_2*F_fast(q)
  !    end do
  !    p = p + nrespa*dt_2*F_slow(q)
  !    SAMPLE
  !  end do
  ! end subroutine
  !---------------------------
  subroutine BAOAB(nsteps)
    integer, intent(in) :: nsteps
    integer :: step, substep
    print*, '# BAOAB: n, nrespa, seed = ', int(n), nrespa, seed(1)
    do step = 1, nsteps
      p = p + nrespa*dt_2*F_slow(q)
      do substep = 1, nrespa
        p = p + dt_2*F_fast(q)
        q = q + dt_2*p
        p = O1*p + O2*normalvec
        q = q + dt_2*p
        p = p + dt_2*F_fast(q)
      end do
      p = p + nrespa*dt_2*F_slow(q)
      SAMPLE
    end do
  end subroutine
  !---------------------------
  subroutine memory_BAOAB(nsteps)
    integer, intent(in) :: nsteps
    integer :: step, substep
    real(8) :: F(2)
    print*, '# memory_BAOAB: n, nrespa, seed = ', int(n), nrespa, seed(1)
    do step = 1, nsteps
      F = F_slow(q)
      do substep = 1, nrespa
        p = p + dt_2*(F_fast(q) + F)
        q = q + dt_2*p
        p = O1*p + O2*normalvec
        q = q + dt_2*p
        p = p + dt_2*(F_fast(q) + F)
      end do
      p = p + nrespa*dt_2*(F_slow(q) - F)
      SAMPLE
    end do
  end subroutine
  !---------------------------
  subroutine middle_NHC(nsteps)
    integer, intent(in) :: nsteps
    integer :: step, substep
    print*, '# middle_NHC: n, nrespa, seed = ', int(n), nrespa, seed(1)
    do step = 1, nsteps
      p = p + nrespa*dt_2*F_slow(q)
      do substep = 1, nrespa
        p = p + dt_2*F_fast(q)
        q = q + dt_2*p
        call nhc % integrate( p, dt )
        q = q + dt_2*p
        p = p + dt_2*F_fast(q)
      end do
      p = p + nrespa*dt_2*F_slow(q)
      SAMPLE
    end do
  end subroutine
  !---------------------------
  subroutine xi_respa_NHC(nsteps)
    integer, intent(in) :: nsteps
    integer :: step, substep
    print*, '# xi_respa_NHC: n, nrespa, seed = ', int(n), nrespa, seed(1)
    do step = 1, nsteps
      call nhc % integrate( p, dt_2 )
      p = p + nrespa*dt_2*F_slow(q)
      do substep = 1, nrespa
        if (substep > 1) call nhc % integrate( p, dt_2 )
        p = p + dt_2*F_fast(q)
        q = q + dt_2*p
        q = q + dt_2*p
        p = p + dt_2*F_fast(q)
        if (substep < nrespa) call nhc % integrate( p, dt_2 )
      end do
      p = p + nrespa*dt_2*F_slow(q)
      call nhc % integrate( p, dt_2 )
      SAMPLE
    end do
  end subroutine
  !---------------------------
  subroutine xo_respa_NHC(nsteps)
    integer, intent(in) :: nsteps
    integer :: step, substep
    print*, '# xo_respa_NHC: n, nrespa, seed = ', int(n), nrespa, seed(1)
    do step = 1, nsteps
      call xo_nhc % integrate( p, nrespa*dt_2 )
      p = p + nrespa*dt_2*F_slow(q)
      do substep = 1, nrespa
        p = p + dt_2*F_fast(q)
        q = q + dt*p
        p = p + dt_2*F_fast(q)
      end do
      p = p + nrespa*dt_2*F_slow(q)
      call xo_nhc % integrate( p, nrespa*dt_2 )
      SAMPLE
    end do
  end subroutine
  !---------------------------
  subroutine middle_IsoK_NHC(nsteps)
    integer, intent(in) :: nsteps
    integer :: step, substep
    print*, '# middle_IsoK_NHC: n, nrespa, seed = ', int(n), nrespa, seed(1)
    do step = 1, nsteps
      p = p + nrespa*dt_2*F_slow(q)
      do substep = 1, nrespa
        p = p + dt_2*F_fast(q)
        q = q + dt_2*v(p)
        call isok_nhc % integrate( p, dt )
        q = q + dt_2*v(p)
        p = p + dt_2*F_fast(q)
      end do
      p = p + nrespa*dt_2*F_slow(q)
      SAMPLE
    end do
  end subroutine
  !---------------------------
  subroutine middle_SubK_NHC(nsteps)
    integer, intent(in) :: nsteps
    integer :: step, substep
    print*, '# middle_SubK_NHC: n, nrespa, seed = ', int(n), nrespa, seed(1)
    do step = 1, nsteps
      p = p + nrespa*dt_2*F_slow(q)
      do substep = 1, nrespa
        p = p + dt_2*F_fast(q)
        q = q + dt_2*v(p)
        call subk_nhc % integrate( p, dt )
        q = q + dt_2*v(p)
        p = p + dt_2*F_fast(q)
      end do
      p = p + nrespa*dt_2*F_slow(q)
      SAMPLE
    end do
  end subroutine
  !---------------------------
  subroutine middle_NHL(nsteps)
    integer, intent(in) :: nsteps
    integer :: step, substep
    real(8) :: vs(2)
    print*, '# middle_NHL: n, nrespa, seed = ', int(n), nrespa, seed(1)
    do step = 1, nsteps
      p = p + nrespa*dt_2*F_slow(q)
      do substep = 1, nrespa
        p = p + dt_2*F_fast(q)
        vs = v(p)
        q = q + dt_2*vs
        v_eta = v_eta + (dt_2/Q_eta)*(p*vs - 1d0)
        p = exp(-dt_2*v_eta)*p
        v_eta = O1*v_eta + (O2/sqrt(Q_eta))*normalvec
        p = exp(-dt_2*v_eta)*p
        vs = v(p)
        v_eta = v_eta + (dt_2/Q_eta)*(p*vs - 1d0)
        q = q + dt_2*vs
        p = p + dt_2*F_fast(q)
      end do
      p = p + nrespa*dt_2*F_slow(q)
      SAMPLE
    end do
  end subroutine
  !---------------------------
  subroutine middle_SINR(nsteps)
    integer, intent(in) :: nsteps
    integer :: step, substep
    real(8) :: vs(2)
    print*, '# middle_SINR: n, nrespa, seed = ', int(n), nrespa, seed(1)
    do step = 1, nsteps
      p = p + nrespa*dt_2*F_slow(q)
      do substep = 1, nrespa
        p = p + dt_2*F_fast(q)
        vs = v(p)
        q = q + dt_2*vs
        v_eta = v_eta + (dt_2/Q_eta)*((n+1)/n*vs*vs - 1d0)
        p = sqrt_n*asinh(exp(-dt_2*v_eta)*sinh(p/sqrt_n))
        v_eta = O1*v_eta + (O2/sqrt(Q_eta))*normalvec
        p = sqrt_n*asinh(exp(-dt_2*v_eta)*sinh(p/sqrt_n))
        vs = v(p)
        v_eta = v_eta + (dt_2/Q_eta)*((n+1)/n*vs*vs - 1d0)
        q = q + dt_2*vs
        p = p + dt_2*F_fast(q)
      end do
      p = p + nrespa*dt_2*F_slow(q)
      SAMPLE
    end do
  end subroutine
  !---------------------------
  subroutine middle_SINR_alt(nsteps)
    integer, intent(in) :: nsteps
    integer :: step, substep
    real(8) :: vs(2)
    print*, '# middle_SINR_alt: n, nrespa, seed = ', int(n), nrespa, seed(1)
    do step = 1, nsteps
      p = p + nrespa*dt_2*F_slow(q)
      do substep = 1, nrespa
        p = p + dt_2*F_fast(q)
        q = q + dt_2*v(p)
        p = sqrt_n*asinh(exp(-dt_2*v_eta)*sinh(p/sqrt_n))
        vs = v(p)
        v_eta = O1*v_eta + (1 - O1)*((n+1)/n*vs*vs - 1d0)/(Q_eta*gamma) + (O2/sqrt(Q_eta))*normalvec
        p = sqrt_n*asinh(exp(-dt_2*v_eta)*sinh(p/sqrt_n))
        q = q + dt_2*v(p)
        p = p + dt_2*F_fast(q)
      end do
      p = p + nrespa*dt_2*F_slow(q)
      SAMPLE
    end do
  end subroutine
  !---------------------------
end program harmonic
