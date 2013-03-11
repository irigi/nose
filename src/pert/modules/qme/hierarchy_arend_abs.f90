program main
 implicit none
   integer, parameter:: NSB = 12 ! number of SB coupling terms
   integer, parameter:: tmax = 6 ! was 8 for N=10
   integer, parameter:: Ntimestep = 200
   real*8, parameter:: dt = 0.01

   integer, parameter:: Nsys = 12 ! system size
   real*8, parameter:: omega_0 = 5.24
   real*8, parameter:: JAT = 1.0
   real*8, parameter:: JAA = 0.875
   real*8, parameter:: JTT = 0.685
   real*8, parameter:: pa = 3.1415927 ! pulse angle

   integer:: Nind ! number of indices in hierarchy
   complex*16, allocatable:: rho(:,:), prhox(:,:), prhodx(:,:), mu(:,:)
     complex*16, allocatable:: HS(:,:)
   complex*16, allocatable:: V(:, :,:)
   complex*16, allocatable:: opLeft(:,:,:), opPlusLeft(:,:,:,:), opMinLeft(:,:,:,:)
   real*8, allocatable:: lambda(:), beta(:), gamma(:), Dtrans(:), Dlong(:)
   real*8, parameter:: pi = 3.1415927

   complex*16, parameter:: izero = dcmplx(0.0, 0.0)
   complex*16, parameter:: iconst = dcmplx(0.0, 1.0)

   complex*16, allocatable:: signal(:), signalx(:), signaly(:), signalz(:)

   integer:: Nhier = 1! number of density matrices in hierarchy
   integer:: tt, s

   integer, allocatable :: perm(:,:)
   integer, allocatable :: permplus(:,:), permmin(:,:)

integer:: tierstart, tier, kk1, kk2, currentindex, nn
logical:: permexists
integer, allocatable:: currentperm(:)


   Nind = NSB 

do tt = 1, tmax
  Nhier = Nhier + numpermt(tt)
end do

write(*,*) 'number of elements in hierarchy = ', Nhier

 allocate(rho(Nhier+1, Nsys),stat=s) ! +1: store zeros in the last density matrix
 if (s .ne. 0) stop 'memory allocation FAILED'
! allocate(prhodt(Nhier+1, Nsys, Nsys),stat=s)
! if (s .ne. 0) stop 'memory allocation FAILED'
 allocate(prhodx(Nhier+1, Nsys),stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
 allocate(prhox(Nhier+1, Nsys),stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
 rho = 0.0
 prhodx = 0.0
 prhox = 0.0
 !allocate(prhodm(Nhier+1, Nsys, Nsys),stat=s)
 !if (s .ne. 0) stop 'memory allocation FAILED'
 allocate(HS(Nsys, Nsys), stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
 Hs = 0.0
 allocate(V(NSB, Nsys, Nsys), stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
V = 0.0

 allocate(beta(NSB),stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
 allocate(gamma(NSB),stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
 allocate(lambda(NSB),stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
 allocate(Dlong(NSB),stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
 allocate(Dtrans(NSB),stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'

 allocate(perm(Nhier, Nind),stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
allocate(currentperm(Nind),stat=s)
if (s .ne. 0) stop 'memory allocation FAILED'

  allocate(opLeft(Nhier, Nsys, Nsys),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
  allocate(opPlusLeft(Nhier, NSB, Nsys, Nsys),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
  allocate(opMinLeft(Nhier, NSB, Nsys, Nsys),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'

  allocate(mu(Nsys, 3),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'

allocate(permplus(Nhier, NSB),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'

allocate(permmin(Nhier, NSB),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'

allocate(signal(Ntimestep),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
allocate(signalx(Ntimestep),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
allocate(signaly(Ntimestep),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
allocate(signalz(Ntimestep),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'



! parameters for system-bath coupling
lambda(1) = 10.0 ! 21.5
gamma(1) = 2.69
beta(1) = 0.19
Dlong(1) = 1.0
Dtrans(1) = 0.0

do s = 2, NSB
  lambda(s) = lambda(1)
  gamma(s) = gamma(1)
  beta(s) = beta(1)
  Dlong(s) = Dlong(1)
  Dtrans(s) = Dtrans(1)
end do


! system Hamiltonian
do s = 1, Nsys/2
  HS(2*s-1, 2*s-1) = omega_0/2
  HS(2*s, 2*s) = -omega_0/2
  HS(2*s-1, 2*s) = JAT
  HS(2*s, 2*s-1) = JAT
end do

do s = 1, (Nsys/2 - 1)
  HS(2*s-1, 2*s+1) = dcmplx(JAA)
  HS(2*s+1, 2*s-1) = dcmplx(JAA)
  HS(2*s, 2*s+2) = dcmplx(JTT)
  HS(2*s+2, 2*s) = dcmplx(JTT)
end do

! transition dipoles
mu(1,1) = 1.0
mu(1,2) = 0.0
mu(1,3) = 0.0
mu(2,1) = cos(pi*117/180)
mu(2,2) = sin(pi*117/180)
mu(2,3) = 0.0
do s=2, Nsys/2
  mu(2*s-1,1) = cos((s-1)*pi*36/180)
  mu(2*s-1,2) = sin((s-1)*pi*36/180)
  mu(2*s-1,3) = 0.0
  mu(2*s,1) = cos(((s-1)*36+117)*pi/180)
  mu(2*s,2) = sin(((s-1)*36+117)*pi/180)
  mu(2*s,3) = 0.0
end do

! system-bath coupling, uncorrelated
do s=1, Nsys
  V(s,s,s) = dcmplx(Dlong(1))
end do



! build index
tierstart = 1 ! first element of the current tier
currentindex = 2

currentperm = 0.0
perm = 0.0

do tier = 1, tmax
  tierstart = tierstart + numpermt(tier-1)
  !write(*,*) "tier = ", tier, ", tierstart = ", tierstart
  do kk1 = currentindex-numpermt(tier-1), currentindex-1 ! loop over all elements in the previous tier
    do kk2 = 1, Nind ! try to add a ball at position kk2
        currentperm(:) = perm(kk1, :)
        currentperm(kk2) = currentperm(kk2) + 1
        permexists = .False.
        do nn = tierstart, currentindex-1 ! check if it is already in the list
           if ( ALL(perm(nn, :) == currentperm) ) then
              permexists = .True.
!              write(*,*) "permutation ", currentperm, " exists"
           end if
        end do
        if (permexists == .False.) then ! if not, add it to the list
          perm(currentindex, :) = currentperm(:)
          currentindex = currentindex + 1                       
!          write(*,*) "new permutation ", currentperm
        end if
        
     end do
   end do
end do

! cache plus and minus
do nn = 1, Nhier
  do s = 1, NSB
    permplus(nn, s) = nplus(nn, s)
    permmin(nn, s) = nmin(nn, s)
  end do
end do

write(*,*) "index complete"



call initmult

! Initial condition, x polarization
do nn = 1, Nhier
do s=1, Nsys
  rho(nn, s) = mu(s, 1)
end do
end do

do nn = 1, Ntimestep
   signalx(nn) = 0
   do s = 1, Nsys
     signalx(nn) = signalx(nn) + mu(s, 1) * rho(1, s)
   end do
   call propagate(dt)
end do

! y polarization
rho(:,:) = 0.0
do nn = 1, Nhier
do s=1, Nsys
  rho(nn, s) = mu(s, 2)
end do
end do

do nn = 1, Ntimestep
   signaly(nn) = 0
   do s = 1, Nsys
     signaly(nn) = signaly(nn) + mu(s, 2) * rho(1, s)
   end do
   call propagate(dt)
end do

! z polarization
rho(:,:) = 0.0
do nn = 1, Nhier
do s=1, Nsys
  rho(nn, s) = mu(s, 3)
end do
end do

do nn = 1, Ntimestep
   signalz(nn) = 0
   do s = 1, Nsys
     signalz(nn) = signalz(nn) + mu(s, 3) * rho(1, s)
   end do
   call propagate(dt)
end do



open (55, File='out.dat' )
  do nn = 1, Ntimestep
     signal(nn) = (signalx(nn) + signaly(nn) + signalz(nn))/3.0
     write(55,*), (nn-1)*dt, real(signal(nn)), imag(signal(nn))
  end do

close(55)



deallocate(currentperm)
deallocate(perm)
deallocate(lambda)
deallocate(gamma)
deallocate(beta)
deallocate(V)
deallocate(HS)
deallocate(prhodx)
deallocate(prhox)
deallocate(rho)



contains

function fact(x)
  real*8:: j, fact, result
  integer*4, intent(in):: x
  result = 1
  do j=1, x
   result = result * j
  end do
  fact = result
end function fact

function numpermt(tier) ! number of permutations in tier
  integer, intent(in):: tier
  integer*4 :: numpermt, tmp
      numpermt = int(fact(Nind+tier-1)/(fact(Nind-1)*fact(tier)))
end function numpermt


function permindex (permutation)
   integer:: permindex
   integer, intent(in):: permutation(Nind)
   logical:: isfound
   integer:: findex, tt, fend

        isfound = .False.
        findex = 1
        do tt = 1, sum(permutation)-1 ! find the tier where this permutation should be
                findex = findex + numpermt (tt)
        end do
        fend = findex + numpermt(sum(permutation))+1
        if (findex >= Nhier) then ! the permutation does not exist
                isfound = .True.
                findex = 0
        end if
        do while ((isfound == .false.) .and. (findex < fend))  ! actually, findex >= fend should never happen if all permutations are present in perm
                                                 ! but we include it for use when perm can be pruned.
          if (ALL(perm(findex,:) == permutation(:))) then
             isfound = .True.
          end if
          findex = findex + 1
        end do

        if (isfound == .false.) then
                findex = 0
        end if
        
        permindex = findex - 1

end function permindex

function nplus (nin, jin)
  integer, intent(in):: nin, jin
  integer:: nplus
     
  currentperm = perm(nin, :)
  currentperm(jin) = currentperm(jin) + 1
     
  nplus = permindex(currentperm)

  if (nplus == -1) then
    nplus = Nhier + 1 ! the Nhier + 1 'th density matrix contains only zeros
  end if

end function nplus


function nmin (nin, jin)
  integer, intent(in):: nin, jin
  integer:: nmin
  
  currentperm = perm(nin, :)

  nmin = -1
  if (currentperm(jin) > 0) then
    currentperm(jin) = currentperm(jin) - 1
    nmin = permindex(currentperm)
  end if

  if (nmin == -1) then
    nmin = Nhier + 1 ! the Nhier + 1 'th density matrix contains only zeros
  end if

end function nmin


function nu(j, m)
  integer, intent(in):: j, m
  complex*16:: nu
  if (m == 0) then
    nu = gamma(j)
  else 
    nu = 2*pi*m/beta(j)
  end if
end function nu

function cot(x)
  real*8, intent(in):: x
  real*8:: cot
  cot = cos(x) / sin(x)
end function cot

function cconst(j)
  integer, intent(in):: j
  complex*16:: cconst

  cconst = 2*lambda(j)/beta(j) - dcmplx(0.0,1.0) * lambda(j) * gamma(j)

end function cconst

subroutine initmult
integer:: n, j
complex*16, parameter:: iconst = dcmplx(0.0, 1.0)
complex*16:: musum, nu1, jsum
complex*16, allocatable:: identity(:,:)

allocate(identity(Nsys, Nsys))

do n=1, Nsys
  identity(n, n) = 1.0
end do

opLeft = 0.0

do n = 1, Nhier

  musum = 0
  do j=1, NSB
    musum = musum + perm(n, j) * gamma(j)
  end do

  opLeft(n,:,:) = opLeft(n,:,:) - iconst * HS
  opLeft(n,:,:) = opLeft(n,:,:) - musum * identity

  do j=1, NSB
    ! first low temperature correction term, see Ishizaki PNAS 2009 
    nu1 = 2*pi/beta(j)
    jsum = 2*(lambda(j) / beta(j)) * 2*gamma(j)/(nu1*nu1 - gamma(j)*gamma(j))
    opLeft(n,:,:) = opLeft(n,:,:) - jsum * MATMUL(V(j,:,:), V(j,:,:))
    
   opPlusLeft(n,j,:,:) =  opPlusLeft(n,j,:,:) - iconst * V(j,:,:)
 
   opMinLeft(n,j,:,:) =  opMinLeft(n,j,:,:) - iconst*perm(n, j)*cconst(j) * V(j,:,:)
 
   ! first low temperature correction term, see Ishizaki PNAS 2009 
   opMinLeft(n,j,:,:) = opMinLeft(n,j,:,:) - iconst * perm(n,j)* jsum * gamma(j) * V(j,:,:)
   
  end do

end do

deallocate(identity)

end subroutine initmult


subroutine Lmult (rhoin, result)
complex*16, intent(in):: rhoin(:,:)
complex*16:: result(:,:)
integer:: n,j, np
complex*16, parameter:: iconst = dcmplx(0.0, 1.0)


result(:,:) = 0.0

!$OMP PARALLEL
do n = 1, Nhier

  result(n,:) = result(n,:) + MATMUL(opLeft(n,:,:), rhoin(n,:))

  do j=1, NSB
    result(n,:) = result(n,:) + MATMUL(opPlusLeft(n,j, :,:), rhoin(permplus(n,j),:))
    result(n,:) = result(n,:) + MATMUL(opMinLeft(n,j, :,:), rhoin(permmin(n,j),:))
  end do

end do
!$OMP END PARALLEL

end subroutine Lmult

subroutine propagate (dt)
  real*8, intent(in):: dt
  

!  second order
call Lmult(rho, prhox)
rho = rho + dt * 0.25 * prhox
prhodx = rho + 0.66666666666666666667 * dt * prhox
call Lmult(prhodx, prhox) ! prhox -> prhodt for speed instead of memory
!rho = rho + dt * (0.25*prhox + 0.75*prhodt) for speed instead of memory
rho = rho + dt * 0.75*prhox


! fourth order
!  call Lmult(rho, prhodt)
!  prhox = rho + (dt/2.0) * prhodt
!  call Lmult(prhox, prhodx)
!  prhox = rho + (dt/2.0) * prhodx
!  call Lmult(prhox, prhodm)
!  prhox = rho + dt * prhodm
!  prhodm = prhodm + prhodx
!  call Lmult(prhox, prhodx)
!  rho = rho + (dt/6.0)*(prhodt + prhodx + 2.0*prhodm)
end subroutine propagate


end program main
