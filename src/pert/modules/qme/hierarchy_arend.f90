! at the moment, there is no loop over polarizations
! this is OK, as long as the two dipoles are parallel
program main
 implicit none
   integer, parameter:: Nsys = 2 ! system size
   integer, parameter:: NSB = 2 ! number of SB coupling terms
   integer, parameter:: tmax = 8 ! depth of the hierarchy
   integer, parameter:: Ntimestept1 = 200 ! number of time steps during t1 in outer loop
   integer, parameter:: Ntimestept1in = 25 ! number of time steps during t1 in inner loop
   integer, parameter:: Ntimestept2 = 0 ! number of time steps during t2 in outer loop
   integer, parameter:: Ntimestept2in = 25 ! number of time steps during t2 in inner loop
   integer, parameter:: Ntimestept3 = 200 ! number of time steps during t3 in outer loop
   integer, parameter:: Ntimestept3in = 25 ! number of time steps during t3 in inner loop

   real*8, parameter:: dt = 0.002

   logical, parameter:: calculatereph = .true.
   logical, parameter:: calculatenonreph = .true.


   integer:: Nind ! number of indices in hierarchy
   complex*16, allocatable:: rho1(:,:), prhox1(:,:), prhodx1(:,:)
   complex*16, allocatable:: rho2(:,:,:), prhox2(:,:,:), prhodx2(:,:,:)
   complex*16, allocatable:: rho3s(:,:), prhox3s(:,:), prhodx3s(:,:)
   complex*16, allocatable:: rho3(:,:,:), prhox3(:,:,:), prhodx3(:,:,:)
   complex*16, allocatable:: HS(:,:), HS2(:,:)
   complex*16, allocatable:: V(:, :,:), V2(:,:,:), mu(:,:)
   complex*16, allocatable:: opLeft2(:,:,:), opRight2(:,:,:), opLRLeft2(:,:,:,:), opLRRight2(:,:,:,:)
   complex*16, allocatable:: opPlusLeft2(:,:,:,:), opPlusRight2(:,:,:,:), opMinLeft2(:,:,:,:), opMinRight2(:,:,:,:)
   complex*16, allocatable:: opLeft3(:,:,:), opRight3(:,:,:), opLRLeft3(:,:,:,:), opLRRight3(:,:,:,:)
   complex*16, allocatable:: opPlusLeft3(:,:,:,:), opPlusRight3(:,:,:,:), opMinLeft3(:,:,:,:), opMinRight3(:,:,:,:)
   complex*16, allocatable:: opLeft1(:,:,:), opPlusLeft1(:,:,:,:), opMinLeft1(:,:,:,:)
   complex*16, allocatable:: signal(:,:), signalpar(:,:), signalper(:,:)
   real*8, allocatable:: lambda(:), beta(:), gamma(:), Dtrans(:), Dlong(:)
   real*8, parameter:: pi = 3.1415927
   complex*16:: pc, ps, rhosum, mem, orcoeffpar
   complex*16, parameter:: izero = dcmplx(0.0, 0.0)
   complex*16, parameter:: iconst = dcmplx(0.0, 1.0)

   integer:: Nhier = 1! number of density matrices in hierarchy
   integer:: tt, s, nin, s2, nt1, nt1in, nt3, nt2, w1, w2

   integer, allocatable :: perm(:,:)
  integer, allocatable :: permplus(:,:), permmin(:,:)

integer:: tierstart, tier, kk1, kk2, currentindex, nn, n, m, np, mp, w, wp
integer:: dir1, dir2, dir3, dir4, pol
logical:: permexists
integer, allocatable:: currentperm(:)

! polarization components
! during t1, we need delta = x (1), y (2) and z (3)
! during t2, we have delta gamma = xx (1), xy (2), xz (3), yx (4), yy (5),
! yz (6), zx (7), zy (8), zz (9)
! during t2, we then need to multiply these to form
! delta gamma beta = xxx (1), xxy (2), xxz(3), xyx (4), xyy(5), xzx(6),
! xzz (7), yxx (8), yxy (9), yyx (10), yyy (11), yyz (12), yzy (13), 
! yzz (14), zxx (15), zxz (16), zyy (17), zyz (18), zzx (19), zzy (20), 
! zzz (21)

integer, parameter:: dirinda(21) = (/ 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3 /)
integer, parameter:: dirindb(21) = (/ 1, 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3, 3 /)
integer, parameter:: dirindc(21) = (/ 1, 2, 3, 1, 2, 1, 3, 1, 2, 1, 2, 3, 2, 3, 1, 3, 2, 3, 1, 2, 3 /)
integer, parameter:: dirindd(21) = (/ 1, 2, 3, 2, 1, 3, 1, 2, 1, 1, 2, 3, 3, 2, 3, 1, 3, 2, 1, 2, 3 /)

! xxxx, xxyy, xxzz, xyxy, xyyx, xzxz, xzzx, yxxy, yxyx, yyxx, 
! yyyy, yyzz, yzyz, yzzy, zxxz, zxzx, zyyz, zyzy, zzxx, zzyy, zzzz

real*8, parameter:: orcoeffZZYY(21) = (/ 2.0/30, 4.0/30, 4.0/30, -1.0/30, -1.0/30, -1.0/30, -1.0/30, -1.0/30, -1.0/30, 4.0/30, 2.0/30, 4.0/30, -1.0/30, -1.0/30,-1.0/30, -1.0/30, -1.0/30, -1.0/30, 4.0/30, 4.0/30, 2.0/30 /)

real*8, parameter:: orcoeffZZZZ(21) = (/ 6.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 6.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 6.0/30 /)

real*8, parameter:: orcoeffZYZY(21) = (/ 2.0/30, -1.0/30, -1.0/30, 4.0/30, -1.0/30, 4.0/30, -1.0/30, -1.0/30, 4.0/30, -1.0/30, 2.0/30, -1.0/30, 4.0/30, -1.0/30, -1.0/30, 4.0/30, -1.0/30, 4.0/30, -1.0/30, -1.0/30, 2.0/30 /)


Nind = NSB 

do tt = 1, tmax
  Nhier = Nhier + numpermt(tt)
end do

write(*,*) 'number of elements in hierarchy = ', Nhier

 allocate(rho1(Nhier+1, Nsys),stat=s) ! +1: store zeros in the last density matrix
 if (s .ne. 0) stop 'memory allocation FAILED'
 allocate(prhodx1(Nhier+1, Nsys),stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
 allocate(prhox1(Nhier+1, Nsys),stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
 rho1 = 0.0
 prhodx1 = 0.0
 prhox1 = 0.0

 allocate(rho2(Nhier+1, Nsys, Nsys),stat=s) ! +1: store zeros in the last density matrix
 if (s .ne. 0) stop 'memory allocation FAILED'
 allocate(prhodx2(Nhier+1, Nsys, Nsys),stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
 allocate(prhox2(Nhier+1, Nsys, Nsys),stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
 rho2 = 0.0
 prhodx2 = 0.0
 prhox2 = 0.0

 allocate(rho3(Nhier+1, (Nsys*(Nsys-1))/2, Nsys),stat=s) ! +1: store zeros in the last density matrix
 if (s .ne. 0) stop 'memory allocation FAILED'
 allocate(prhodx3(Nhier+1, (Nsys*(Nsys-1))/2, Nsys),stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
 allocate(prhox3(Nhier+1, (Nsys*(Nsys-1))/2, Nsys),stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
 rho3 = 0.0
 prhodx3 = 0.0
 prhox3 = 0.0

 allocate(rho3s(Nhier+1, Nsys),stat=s) ! +1: store zeros in the last density matrix
 if (s .ne. 0) stop 'memory allocation FAILED'
 allocate(prhodx3s(Nhier+1, Nsys),stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
 allocate(prhox3s(Nhier+1, Nsys),stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
 rho3s = 0.0
 prhodx3s = 0.0
 prhox3s = 0.0


 allocate(HS(Nsys, Nsys), stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
 Hs = 0.0
 allocate(V(NSB, Nsys, Nsys), stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
V = 0.0

 allocate(HS2((Nsys*(Nsys-1))/2, (Nsys*(Nsys-1))/2), stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
 HS2 = 0.0

 allocate(V2(NSB, (Nsys*(Nsys-1))/2, (Nsys*(Nsys-1))/2), stat=s)
 if (s .ne. 0) stop 'memory allocation FAILED'
V2 = 0.0


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

allocate(opLeft1(Nhier, Nsys, Nsys),stat=s)
if (s .ne. 0) stop 'memory allocation FAILED'
allocate(opPlusLeft1(Nhier, NSB, Nsys, Nsys),stat=s)
if (s .ne. 0) stop 'memory allocation FAILED'
allocate(opMinLeft1(Nhier, NSB, Nsys, Nsys),stat=s)
if (s .ne. 0) stop 'memory allocation FAILED'

  allocate(opLeft2(Nhier, Nsys, Nsys),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
  allocate(opRight2(Nhier, Nsys, Nsys),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
  allocate(opLRLeft2(Nhier, NSB, Nsys, Nsys),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
  allocate(opLRRight2(Nhier, NSB, Nsys, Nsys),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
  allocate(opPlusLeft2(Nhier, NSB, Nsys, Nsys),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
  allocate(opPlusRight2(Nhier, NSB, Nsys, Nsys),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
  allocate(opMinLeft2(Nhier, NSB, Nsys, Nsys),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
  allocate(opMinRight2(Nhier, NSB, Nsys, Nsys),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'

  allocate(opLeft3(Nhier, (Nsys*(Nsys-1))/2, (Nsys*(Nsys-1))/2),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
  allocate(opRight3(Nhier, Nsys, Nsys),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
  allocate(opLRLeft3(Nhier, NSB, (Nsys*(Nsys-1))/2, (Nsys*(Nsys-1))/2),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
  allocate(opLRRight3(Nhier, NSB, Nsys, Nsys),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
  allocate(opPlusLeft3(Nhier, NSB, (Nsys*(Nsys-1))/2, (Nsys*(Nsys-1))/2),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
  allocate(opPlusRight3(Nhier, NSB, Nsys, Nsys),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
  allocate(opMinLeft3(Nhier, NSB, (Nsys*(Nsys-1))/2, (Nsys*(Nsys-1))/2),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
  allocate(opMinRight3(Nhier, NSB, Nsys, Nsys),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'


allocate(permplus(Nhier, NSB),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'

allocate(permmin(Nhier, NSB),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'

allocate(signal(Ntimestept1, Ntimestept3),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
allocate(signalpar(Ntimestept1, Ntimestept3),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'
allocate(signalper(Ntimestept1, Ntimestept3),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'

  allocate(mu(Nsys, 3),stat=s)
  if (s .ne. 0) stop 'memory allocation FAILED'

! parameters for system-bath coupling
lambda(1) = 0.5 ! 10.0
gamma(1) = 2.0
beta(1) = 0.13
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
HS(1,1) = 3.0
HS(2,2) = -3.0
!HS(3,3) = 0.0
HS(1,2) = -5.0
HS(2,1) = -5.0
!HS(2,3) = 1.0
!HS(3,2) = 1.0

! transition dipoles
mu(1,1) = 1.0
mu(1,2) = 0.0
mu(1,3) = 0.0
mu(2,1) = 1.0
mu(2,2) = 0.0
mu(2,3) = 0.0
!mu(3,1) = 1.0
!mu(3,2) = 0.0
!mu(3,3) = 0.0

! system-bath coupling, no correlation
do s=1, NSB
  V(s,s,s) = dcmplx(Dlong(s))
end do

! 2-quantum Hamiltonian and system-bath coupling operators
do n = 1, Nsys
  do m = (n+1), Nsys
    w = (n-1)*Nsys + m - (n*(n+1))/2

    do np = 1, Nsys
      do mp = (np+1), Nsys
        wp = (np-1)*Nsys + mp - (np*(np+1))/2

	if (n == np) then
           HS2(w, wp) = HS2(w, wp) + HS(m, mp)
           do s = 1, NSB
	     V2(s, w, wp) = V2(s, w, wp) + V(s, m, mp)      
           end do
        end if
	if (n == mp) then
           HS2(w, wp) = HS2(w, wp) + HS(m, np)
           do s = 1, NSB
	     V2(s, w, wp) = V2(s, w, wp) + V(s, m, np)      
           end do
        end if
	if (m == np) then
           HS2(w, wp) = HS2(w, wp) + HS(n, mp)
           do s = 1, NSB
	     V2(s, w, wp) = V2(s, w, wp) + V(s, n, mp)      
           end do
        end if
	if (m == mp) then
           HS2(w, wp) = HS2(w, wp) + HS(n, np)
           do s = 1, NSB
	     V2(s, w, wp) = V2(s, w, wp) + V(s, n, np)      
           end do
        end if

      end do
     end do
  end do
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


call initmult1
call initmult2
call initmult3


! REPHASING

if (calculatereph) then

signalpar(:,:) = 0.0
signalper(:,:) = 0.0

! loop over polarizations
do pol = 1, 1 ! 21

dir1 = dirinda(pol)
dir2 = dirindb(pol)
dir3 = dirindc(pol)
dir4 = dirindd(pol)
orcoeffpar = orcoeffZZZZ(pol)


! Initial condition
do nn = 1, Nhier
  do s=1, Nsys
    rho1(nn, s) = mu(s, dir1)
  end do
end do


do nt1 = 1, Ntimestept1
  write(*,'(A6,I2,A12,I5,A3,I5)') "pol = ", pol, " / 21; t1 = ", nt1, " / " , Ntimestept1

! GB

! initialize t3
  do nn = 1, Nhier
    mem = 0.0
    do s = 1, Nsys
      mem = mem + mu(s, dir2) * rho1(nn, s)
    end do
    do s = 1, Nsys
      rho3s(nn, s) = mu(s, dir3) * mem
    end do
  end do

 do nt3 = 1, Ntimestept3

     do s = 1, Nsys
        signalpar(nt1, nt3) = signalpar(nt1, nt3) + orcoeffpar * mu(s, dir4) * rho3s(1, s)
     end do

     ! propagate t3
     do nin = 1, Ntimestept3in
       call propagate3s(dt)
     end do


 end do



! SE and IA

  ! initialize t2
  do nn = 1, Nhier
    do s = 1, Nsys
      do s2 = 1, Nsys
        rho2(nn, s, s2) = mu(s, dir2) * rho1(nn, s2)
      end do
    end do
  end do

 ! propagate t2
 do nt2 = 1, Ntimestept2
   do nin = 1, Ntimestept2in
     call propagate2(dt)
   end do
 end do

! SE

 ! initialize t3
  rho3s = 0.0
  do nn = 1, Nhier
    do s = 1, Nsys
      do s2 = 1, Nsys
        rho3s(nn, s) = rho3s(nn, s) + mu(s2, dir3) * rho2(nn, s, s2)
      end do
    end do
  end do

 do nt3 = 1, Ntimestept3

     ! calculate signal
     do s = 1, Nsys
        signalpar(nt1, nt3) = signalpar(nt1, nt3) + orcoeffpar * mu(s, dir4) * rho3s(1, s)
      end do

     ! propagate t3
     do nin = 1, Ntimestept3in
       call propagate3s(dt)
     end do


 end do




! IA
  ! initialize t3
  rho3 = 0.0
  do nn = 1, Nhier
    do s = 1, Nsys
     do w1 = 1, Nsys
       do w2 = (w1+1), Nsys
          w = (w1-1)*Nsys - (w1*(w1+1))/2 + w2
          rho3(nn, w, s) = rho3(nn, w, s) + mu(w2, dir3) * rho2(nn, w1, s)
          rho3(nn, w, s) = rho3(nn, w, s) + mu(w1, dir3) * rho2(nn, w2, s)
        end do
      end do
    end do
  end do


 do nt3 = 1, Ntimestept3

     ! calculate signal
     do w1 = 1, Nsys
       do w2 = (w1+1), Nsys
          w = (w1-1)*Nsys - (w1*(w1+1))/2 + w2
          signalpar(nt1, nt3) = signalpar(nt1, nt3) - orcoeffpar * mu(w2, dir4) * rho3(1, w, w1)
          signalpar(nt1, nt3) = signalpar(nt1, nt3) - orcoeffpar * mu(w1, dir4) * rho3(1, w, w2)
        end do
      end do

     ! propagate t3
     do nin = 1, Ntimestept3in
       call propagate3(dt)
     end do


 end do

 ! propagate t1
 do nin = 1, Ntimestept1in
   call propagate1reph(dt)
 end do


end do

end do ! pol

write(*,*) "rephasing calculation complete"

!output
open (55, File='reph.par.re' )
do nt1 = 1, Ntimestept1
  do nn = 1, Ntimestept3
  write(55,'(F12.6)', advance='no'), real(signalpar(nt1, nn))
  end do
  write(55,*)
end do
close(55)

open (55, File='reph.per.re' )
do nt1 = 1, Ntimestept1
  do nn = 1, Ntimestept3
  write(55,'(F12.6)', advance='no'), real(signalper(nt1, nn))
  end do
  write(55,*)
end do
close(55)

open (55, File='reph.par.im' )
do nt1 = 1, Ntimestept1
  do nn = 1, Ntimestept3
  write(55,'(F12.6)', advance='no'), imag(signalpar(nt1, nn))
  end do
  write(55,*)
end do
close(55)

open (55, File='reph.per.im' )
do nt1 = 1, Ntimestept1
  do nn = 1, Ntimestept3
  write(55,'(F12.6)', advance='no'), imag(signalper(nt1, nn))
  end do
  write(55,*)
end do
close(55)



end if ! calculatereph



! NON-REPHASING

if (calculatenonreph) then

signalpar(:,:) = 0.0
signalper(:,:) = 0.0

! loop over polarizations
do pol = 1,1
dir1 = dirinda(pol)
dir2 = dirindb(pol)
dir3 = dirindc(pol)
dir4 = dirindd(pol)
orcoeffpar = orcoeffZZZZ(pol)



! Initial condition
do nn = 1, Nhier
  do s=1, Nsys
    rho1(nn, s) = mu(s, dir1)
  end do
end do


do nt1 = 1, Ntimestept1
  write(*,'(A6,I2,A12,I5,A3,I5)') "pol = ", pol, " / 21; t1 = ", nt1, " / " , Ntimestept1
 
! GB

! initialize t3
  do nn = 1, Nhier
    mem = 0.0
    do s = 1, Nsys
      mem = mem + mu(s, dir2) * rho1(nn, s)
    end do
    do s = 1, Nsys
      rho3s(nn, s) = mu(s, dir3) * mem
    end do
  end do

 do nt3 = 1, Ntimestept3

     do s = 1, Nsys
        signalpar(nt1, nt3) = signalpar(nt1, nt3) + orcoeffpar * mu(s, dir4) * rho3s(1, s)
     end do

     ! propagate t3
     do nin = 1, Ntimestept3in
       call propagate3s(dt)
     end do


 end do



! SE and IA

  ! initialize t2
  do nn = 1, Nhier
    do s = 1, Nsys
      do s2 = 1, Nsys
        rho2(nn, s, s2) = mu(s2, dir2) * rho1(nn, s)
      end do
    end do
  end do

 ! propagate t2
 do nt2 = 1, Ntimestept2
   do nin = 1, Ntimestept2in
     call propagate2(dt)
   end do
 end do

! SE

 ! initialize t3
  rho3s = 0.0
  do nn = 1, Nhier
    do s = 1, Nsys
      do s2 = 1, Nsys
        rho3s(nn, s) = rho3s(nn, s) + mu(s2, dir3) * rho2(nn, s, s2)
      end do
    end do
  end do

 do nt3 = 1, Ntimestept3

     ! calculate signal
     do s = 1, Nsys
        signalpar(nt1, nt3) = signalpar(nt1, nt3) + orcoeffpar * mu(s, dir4) * rho3s(1, s)
      end do

     ! propagate t3
     do nin = 1, Ntimestept3in
       call propagate3s(dt)
     end do


 end do




! IA
  ! initialize t3
  rho3 = 0.0
  do nn = 1, Nhier
    do s = 1, Nsys
     do w1 = 1, Nsys
       do w2 = (w1+1), Nsys
          w = (w1-1)*Nsys - (w1*(w1+1))/2 + w2
          rho3(nn, w, s) = rho3(nn, w, s) + mu(w2, dir3) * rho2(nn, w1, s)
          rho3(nn, w, s) = rho3(nn, w, s) + mu(w1, dir3) * rho2(nn, w2, s)
        end do
      end do
    end do
  end do


 do nt3 = 1, Ntimestept3

     ! calculate signal
     do w1 = 1, Nsys
       do w2 = (w1+1), Nsys
          w = (w1-1)*Nsys - (w1*(w1+1))/2 + w2
          signalpar(nt1, nt3) = signalpar(nt1, nt3) - orcoeffpar * mu(w2, dir4) * rho3(1, w, w1)
          signalpar(nt1, nt3) = signalpar(nt1, nt3) - orcoeffpar * mu(w1, dir4) * rho3(1, w, w2)
        end do
      end do

     ! propagate t3
     do nin = 1, Ntimestept3in
       call propagate3(dt)
     end do


 end do

 ! propagate t1
 do nin = 1, Ntimestept1in
   call propagate1nonreph(dt)
 end do


end do

end do ! pol

write(*,*) "non-rephasing calculation complete"

!output
open (55, File='nonreph.par.re' )
do nt1 = 1, Ntimestept1
  do nn = 1, Ntimestept3
  write(55,'(F12.6)', advance='no'), real(signalpar(nt1, nn))
  end do
  write(55,*)
end do
close(55)

open (55, File='nonreph.per.re' )
do nt1 = 1, Ntimestept1
  do nn = 1, Ntimestept3
  write(55,'(F12.6)', advance='no'), real(signalper(nt1, nn))
  end do
  write(55,*)
end do
close(55)

open (55, File='nonreph.par.im' )
do nt1 = 1, Ntimestept1
  do nn = 1, Ntimestept3
  write(55,'(F12.6)', advance='no'), imag(signalpar(nt1, nn))
  end do
  write(55,*)
end do
close(55)

open (55, File='nonreph.per.im' )
do nt1 = 1, Ntimestept1
  do nn = 1, Ntimestept3
  write(55,'(F12.6)', advance='no'), imag(signalper(nt1, nn))
  end do
  write(55,*)
end do
close(55)



end if ! calculatenonreph














! deallocate


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



subroutine initmult1
! initialize the propagator for a coherence |1><0|
integer:: n, j
complex*16, parameter:: iconst = dcmplx(0.0, 1.0)
complex*16:: musum, nu1, jsum
complex*16, allocatable:: identity(:,:)

allocate(identity(Nsys, Nsys))

do n=1, Nsys
  identity(n, n) = 1.0
end do

opLeft1 = 0.0
opPlusLeft1 = 0.0
opMinLeft1 = 0.0

do n = 1, Nhier

  musum = 0
  do j=1, NSB
    musum = musum + perm(n, j) * gamma(j)
  end do


  opLeft1(n,:,:) = opLeft1(n,:,:) - iconst * HS

  opLeft1(n,:,:) = opLeft1(n,:,:) - musum * identity

  do j=1, NSB
    ! first low temperature correction term, see Ishizaki PNAS 2009 
    nu1 = 2*pi/beta(j)
    jsum = 2*(lambda(j) / beta(j)) * 2*gamma(j)/(nu1*nu1 - gamma(j)*gamma(j))
    opLeft1(n,:,:) = opLeft1(n,:,:) - jsum * MATMUL(V(j,:,:), V(j,:,:))

   opPlusLeft1(n,j,:,:) =  opPlusLeft1(n,j,:,:) - iconst * V(j,:,:)

   opMinLeft1(n,j,:,:) =  opMinLeft1(n,j,:,:) - iconst*perm(n, j)*cconst(j) * V(j,:,:)

   ! first low temperature correction term, see Ishizaki PNAS 2009 
   opMinLeft1(n,j,:,:) = opMinLeft1(n,j,:,:) - iconst * perm(n,j)*jsum * gamma(j) * V(j,:,:)

  end do

end do

deallocate(identity)

end subroutine initmult1





subroutine initmult2
! initialize the propagator for populations and coherences |1><1'|
integer:: n, j
complex*16, parameter:: iconst = dcmplx(0.0, 1.0)
complex*16:: musum, nu1, jsum
complex*16, allocatable:: identity(:,:)

allocate(identity(Nsys, Nsys))

do n=1, Nsys
  identity(n, n) = 1.0
end do

opLeft2 = 0.0
opRight2 = 0.0
opLRLeft2 = 0.0
opLRRight2 = 0.0
opPlusLeft2 = 0.0
opPlusRight2 = 0.0
opMinLeft2 = 0.0
opMinRight2 = 0.0

do n = 1, Nhier

  musum = 0
  do j=1, NSB
    musum = musum + perm(n, j) * gamma(j)
  end do


  opLeft2(n,:,:) = opLeft2(n,:,:) - iconst * HS
  opRight2(n,:,:) = opRight2(n,:,:) + iconst * HS

  opLeft2(n,:,:) = opLeft2(n,:,:) - musum * identity

  do j=1, NSB

    ! first low temperature correction term, see Ishizaki PNAS 2009 
    nu1 = 2*pi/beta(j)
    jsum = 2*(lambda(j) / beta(j)) * 2*gamma(j)/(nu1*nu1 - gamma(j)*gamma(j))
    opLeft2(n,:,:) = opLeft2(n,:,:) - jsum * MATMUL(V(j,:,:), V(j,:,:))
    opRight2(n,:,:)= opRight2(n,:,:) -  jsum * MATMUL(V(j,:,:), V(j,:,:))
    opLRLeft2(n,j,:,:) =  opLRLeft2(n,j,:,:) + 2*jsum * V(j,:,:)
    opLRRight2(n,j,:,:) =  opLRRight2(n,j,:,:) + V(j,:,:)

   opPlusLeft2(n,j,:,:) =  opPlusLeft2(n,j,:,:) - iconst * V(j,:,:)
   opPlusRight2(n,j,:,:) =  opPlusRight2(n,j,:,:) + iconst * V(j,:,:)

   opMinLeft2(n,j,:,:) =  opMinLeft2(n,j,:,:) - iconst*perm(n, j)*cconst(j) * V(j,:,:)
   opMinRight2(n,j,:,:) = opMinRight2(n,j,:,:) + iconst*perm(n, j)*conjg(cconst(j)) * V(j,:,:)

   ! first low temperature correction term, see Ishizaki PNAS 2009 
   opMinLeft2(n,j,:,:) = opMinLeft2(n,j,:,:) - iconst * perm(n,j)*jsum * gamma(j) * V(j,:,:)
   opMinRight2(n,j,:,:) = opMinRight2(n,j,:,:) + iconst * perm(n,j)*jsum * gamma(j) * V(j,:,:)

  end do

end do

deallocate(identity)

end subroutine initmult2


subroutine initmult3
! initialize the propagator for 2-1 exciton coherences |2><1|
integer:: n, j
complex*16, parameter:: iconst = dcmplx(0.0, 1.0)
complex*16:: musum, nu1, jsum
complex*16, allocatable:: identity(:,:)

allocate(identity((Nsys*(Nsys-1))/2, (Nsys*(Nsys-1))/2))

identity = 0.0
do n=1, (Nsys*(Nsys-1))/2
  identity(n, n) = 1.0
end do

opLeft3 = 0.0
opRight3 = 0.0
opLRLeft3 = 0.0
opLRRight3 = 0.0
opPlusLeft3 = 0.0
opPlusRight3 = 0.0
opMinLeft3 = 0.0
opMinRight3 = 0.0

do n = 1, Nhier

  musum = 0
  do j=1, NSB
    musum = musum + perm(n, j) * gamma(j)
  end do


  opLeft3(n,:,:) = opLeft3(n,:,:) - iconst * HS2
  opRight3(n,:,:) = opRight3(n,:,:) + iconst * HS

  opLeft3(n,:,:) = opLeft3(n,:,:) - musum * identity

  do j=1, NSB
    ! first low temperature correction term, see Ishizaki PNAS 2009 
    nu1 = 2*pi/beta(j)
    jsum = 2*(lambda(j) / beta(j)) * 2*gamma(j)/(nu1*nu1 - gamma(j)*gamma(j))
    opLeft3(n,:,:) = opLeft3(n,:,:) - jsum * MATMUL(V2(j,:,:), V2(j,:,:))
    opRight3(n,:,:)= opRight3(n,:,:) -  jsum * MATMUL(V(j,:,:), V(j,:,:))
    opLRLeft3(n,j,:,:) =  opLRLeft3(n,j,:,:) + 2*jsum * V2(j,:,:)
    opLRRight3(n,j,:,:) =  opLRRight3(n,j,:,:) + V(j,:,:)
   opPlusLeft3(n,j,:,:) =  opPlusLeft3(n,j,:,:) - iconst * V2(j,:,:)
   opPlusRight3(n,j,:,:) =  opPlusRight3(n,j,:,:) + iconst * V(j,:,:)
   opMinLeft3(n,j,:,:) =  opMinLeft3(n,j,:,:) - iconst*perm(n, j)*cconst(j) * V2(j,:,:)
   opMinRight3(n,j,:,:) = opMinRight3(n,j,:,:) + iconst*perm(n, j)*conjg(cconst(j)) * V(j,:,:)
   ! first low temperature correction term, see Ishizaki PNAS 2009 
   opMinLeft3(n,j,:,:) = opMinLeft3(n,j,:,:) - iconst * perm(n,j)*jsum * gamma(j) * V2(j,:,:)
   opMinRight3(n,j,:,:) = opMinRight3(n,j,:,:) + iconst * perm(n,j)*jsum * gamma(j) * V(j,:,:)

  end do

end do

deallocate(identity)

end subroutine initmult3





subroutine Lmult1 (rhoin, result)
complex*16, intent(in):: rhoin(:,:)
complex*16:: result(:,:)
integer:: n,j, np
complex*16, parameter:: iconst = dcmplx(0.0, 1.0)


result(:,:) = 0.0


do n = 1, Nhier

  result(n,:) = result(n,:) + MATMUL(opLeft1(n,:,:), rhoin(n,:))

  do j=1, NSB
    result(n,:) = result(n,:) + MATMUL(opPlusLeft1(n, j, :, :), rhoin(permplus(n,j),:))
    result(n,:) = result(n,:) + MATMUL(opMinLeft1(n,j,:,:), rhoin(permmin(n,j),:))
  end do

end do

end subroutine Lmult1

subroutine propagate1reph (dt)
  real*8, intent(in):: dt
  

!  second order
rho1 = conjg(rho1) ! Lmult1 propates the ket, want the bra for rephasing response.
call Lmult1(rho1, prhox1)
rho1 = rho1 + dt * 0.25 * prhox1
prhodx1 = rho1 + 0.66666666666666666667 * dt * prhox1
call Lmult1(prhodx1, prhox1) ! prhox1 -> prhodt1 for speed instead of memory
!rho1 = rho1 + dt * (0.25*prhox1 + 0.75*prhodt1) for speed instead of memory
rho1 = rho1 + dt * 0.75*prhox1
rho1 = conjg(rho1)


! fourth order
!  call Lmult(rho1, prhodt1)
!  prhox1 = rho1 + (dt/2.0) * prhodt1
!  call Lmult1(prhox1, prhodx1)
!  prhox1 = rho1 + (dt/2.0) * prhodx1
!  call Lmult1(prhox1, prhodm1)
!  prhox1 = rho1 + dt * prhodm1
!  prhodm1 = prhodm1 + prhodx1
!  call Lmult1(prhox1, prhodx1)
!  rho1 = rho1 + (dt/6.0)*(prhodt1 + prhodx1 + 2.0*prhodm1)
end subroutine propagate1reph


subroutine propagate1nonreph (dt)
  real*8, intent(in):: dt
  

!  second order

call Lmult1(rho1, prhox1)
rho1 = rho1 + dt * 0.25 * prhox1
prhodx1 = rho1 + 0.66666666666666666667 * dt * prhox1
call Lmult1(prhodx1, prhox1) 

rho1 = rho1 + dt * 0.75*prhox1

end subroutine propagate1nonreph


subroutine Lmult2 (rhoin, result)
complex*16, intent(in):: rhoin(:,:,:)
complex*16:: result(:,:,:)
integer:: n,j, np
complex*16, parameter:: iconst = dcmplx(0.0, 1.0)


result(:,:,:) = 0.0


do n = 1, Nhier

  result(n,:,:) = result(n,:,:) + MATMUL(opLeft2(n,:,:), rhoin(n,:,:))
  result(n,:,:) = result(n,:,:) + MATMUL(rhoin(n,:,:), opRight2(n,:,:))

  do j=1, NSB
    result(n,:,:) = result(n,:,:) + MATMUL(MATMUL(opLRLeft2(n,j,:,:), rhoin(n,:,:)), opLRRight2(n,j,:,:))
    result(n,:,:) = result(n,:,:) + MATMUL(opPlusLeft2(n,j, :,:), rhoin(permplus(n,j),:,:))
    result(n,:,:) = result(n,:,:) + MATMUL(rhoin(permplus(n,j),:,:), opPlusRight2(n,j,:,:))
    result(n,:,:) = result(n,:,:) + MATMUL(opMinLeft2(n,j, :,:), rhoin(permmin(n,j),:,:))
    result(n,:,:) = result(n,:,:) + MATMUL(rhoin(permmin(n,j),:,:), opMinRight2(n,j,:,:))
  end do

end do

end subroutine Lmult2

subroutine propagate2 (dt)
  real*8, intent(in):: dt
  

!  second order
call Lmult2(rho2, prhox2)
rho2 = rho2 + dt * 0.25 * prhox2
prhodx2 = rho2 + 0.66666666666666666667 * dt * prhox2
call Lmult2(prhodx2, prhox2) ! prhox -> prhodt for speed instead of memory
rho2 = rho2 + dt * 0.75*prhox2
end subroutine propagate2


subroutine Lmult3 (rhoin, result)
complex*16, intent(in):: rhoin(:,:,:)
complex*16:: result(:,:,:)
integer:: n,j, np
complex*16, parameter:: iconst = dcmplx(0.0, 1.0)


result(:,:,:) = 0.0


do n = 1, Nhier

  result(n,:,:) = result(n,:,:) + MATMUL(opLeft3(n,:,:), rhoin(n,:,:))
  result(n,:,:) = result(n,:,:) + MATMUL(rhoin(n,:,:), opRight3(n,:,:))

  do j=1, NSB
    result(n,:,:) = result(n,:,:) + MATMUL(MATMUL(opLRLeft3(n,j,:,:), rhoin(n,:,:)), opLRRight3(n,j,:,:))
    result(n,:,:) = result(n,:,:) + MATMUL(opPlusLeft3(n,j, :,:), rhoin(permplus(n,j),:,:))
    result(n,:,:) = result(n,:,:) + MATMUL(rhoin(permplus(n,j),:,:), opPlusRight3(n,j,:,:))
    result(n,:,:) = result(n,:,:) + MATMUL(opMinLeft3(n,j, :,:), rhoin(permmin(n,j),:,:))
    result(n,:,:) = result(n,:,:) + MATMUL(rhoin(permmin(n,j),:,:), opMinRight3(n,j,:,:))
  end do

end do

end subroutine Lmult3


subroutine propagate3 (dt)
  real*8, intent(in):: dt
  

!  second order
call Lmult3(rho3, prhox3)
rho3 = rho3 + dt * 0.25 * prhox3
prhodx3 = rho3 + 0.66666666666666666667 * dt * prhox3
call Lmult3(prhodx3, prhox3) ! prhox -> prhodt for speed instead of memory
rho3 = rho3 + dt * 0.75*prhox3

end subroutine propagate3


subroutine propagate3s (dt)
  real*8, intent(in):: dt

!  second order
call Lmult1(rho3s, prhox3s)
rho3s = rho3s + dt * 0.25 * prhox3s
prhodx3s = rho3s + 0.66666666666666666667 * dt * prhox3s
call Lmult1(prhodx3s, prhox3s) 
rho3s = rho3s + dt * 0.75*prhox3s

end subroutine propagate3s



end program main
