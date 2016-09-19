module init

implicit none
include 'elem_config90.inc'
integer :: ns !number of species
real :: eps !required accuracy
integer :: nreg !parameter  IN THE CALCULATION OF THE VISCOSITY.
!                        0     default. THE BANANA AND
!                              PFIRSCH SCHLUETER CONTRIBUTIONS
!                              ARE WEIGHTED
!                        1     THE BANANA REGIME IS FORCED.
!                        2     THE PFIRSCH SCHLUETER REGIME IS
!                              FORCED.
!                        NOTE : FORCING THE BANANA REGIME IS NOT 
!                               THE SAME THING AS CALCULATING THE
!                               BANANA PLATEAU CONTRIBUTION. FOR
!                               HIGH COLLISION FREQ. THIS CONTRI-
!                               BUTION USUALLY DECREASES. WHEN
!                               NREG = 1 THERE IS NO SUCH DECR.

integer, save, allocatable :: sigma(:)  !THE SIGMA'S DETERMINE WHICH OF THE COUPLING 
!                        TERMS IN THE EQUATION FOR THE PFIRSCH 
!                        SCHLUETER REGIME ARE TAKEN INTO ACOUNT. 
!                        SIGMA(1)= 0,1  IS THE CROSS COUPLING BETWEEN 
!                                    THE HEAT EQUATION AND THE EQUATION
!                                    FOR THE THIRD LAGUERRE HARM.
!                        SIGMA(2)= 0.,1 IS THE CROSS COUPLING BETWEEN 
!                                    THE EQUATION FOR THE THIRD HARMONIC
!                                    AND THE TEMPERATURE PERTURB.
!                        SIGMA(3)= 0,1 IS THE ENERGY EXCHANGE BETWEEN 
!                                    THE DIFFERENT SPECIES
!                        SIGMA(4)= 0,1 IS THE COUPLING TO THE "DENSITY 
!                                    PERTURB." IN THE EQUATION FOR THE 
!                                    THRIRD LAGUERRE HARMONIC
!                        IN THE APPROXIMATIONS OF HIRSHMANN AND SIGMAR
!                        SIGMA(1)=0 SIGMA(2)=0 SIGMA(3)=1 SIGMA(4)=1
!                        THAT IS THE FIRST TWO TERMS ARE NEGLECTED.
!                        DEFAULT OF ALL: 1.
integer :: nleg ! order of legendre polynomials in expasnsion up to 3(=default).    
integer :: nenergy !parameter: ENERGY SCATTERING IS TAKEN INTO ACCOUNT IN THE 
!                        CALCULATION OF THE VISCOSITY 
!                        0 NOT ACCOUNTED FOR 
!                        1 ACCOUNTED FOR (default)    
integer :: ncof ! PARAMETERS THAT DETERMINES WHETHER ION-
!                        ELECTRON COLLISIONS ARE TAKEN INTO ACCOUNT
!                        0 NOT TAKEN INTO ACCOUNT
!                        1 TAKEN INTO ACCOUNT (NORMAL VALUE, 1) 
! logical :: neogeo !      LOGICAL IF TRUE THEN THE GEOMETRY DEP.
! !                        PARAMETERS ARE RECALCULATED.
! logical :: neofrc !      LOGICAL IF TRUE THEN THE FRICTION AND 
! !                       VISCOSITY MATRIX ARE !NOT! NEWLY CALCULATED
integer :: ic     !       Control: THE CONTRIBUTION FOR WHICH THE COEFFICIENTS
!                        ARE CALCULATED. 
!                        0 THE CLASSICAL PARTICLE FLUX
!                        1 BANANA PLATEAU CONTRIBUTION
!                        2 PFIRSCH SCHLUETER CONTRIBUTION
!                        3 BOTH BANANA-PLATEAU, CLASSICAL, AND P.S.
integer :: ncm           ! maximum number of charge states

integer, save ,allocatable :: nc(:) ! DIM:(ns) number of charge states
real, save ,allocatable :: zsp(:,:) !DIM:(ns,nc) charge number
real, save, allocatable :: m(:) !DIM:(ns) mass in units of proton mass 
real, save, allocatable :: T(:) !DIM:(ns) temperautre in units of kEV
real, save, allocatable :: den(:,:) !DIM:(ns,nc) density in units of 10^19 m^-3
real, save, allocatable :: ds(:,:,:) !DIM:(ns,nc,2); 1 Pressure gradient; 2 Temperature gradient
!                                    !  divide by R^2 to get gkw gradients
!
!                                          D LN P_IJ 
!                              DS(I,J,1) = ---------   
!                                            D rho
!
!                                          D LN T_I
!                              DS(I,J,2) = --------
!                                           D rho


integer :: ISEL !geometry switch: 1 ASSUME HAMADA COORDINATES AND USE THE
!                                   ROUTINE VISGEOM
!                                 2 ASSUME CIRCULAR GEOMETRY AND USE CIRC-
!                                   GEOM TO DETERMINE THE PARAMETERS
!                                 3 READ FROM FILE
!                                 4 READ FROM CHEASE OUTPUT, IN THAT CASE
!                                   THE FLUX SURFACE LABEL IS 
!                                   RHO=(RMAX-RMIN)/2/R0EXP as in GKW
integer :: ISHOT !shot number, only used when ISEL=3
real :: rho !flux surface label
real :: e   ! inverse aspect ration epsilon
real :: q           ! safety factor
real :: Rn          ! major radius in units of R_ref= ?
real :: Bn          ! magnetic field strength in midplane, i.e. poloidal angle= pi/2
real :: eparr     ! THE PARALLEL ELECTRIC FIELD TIMES THE 
!                        MAJOR RADIUS, I.E. THE LOOP VOLTAGE. 
real, save, allocatable :: CFF(:,:,:) !used in main for saving results


integer :: vnlin_order ! which orders are considered 3(default). 0= source ignored 
integer :: vnlin_contribution ! 0(default) for all, 1: ExB only; 2: Curv only; 3: Coriolis only 
real, save, allocatable :: vnlin_drive(:,:,:)
logical :: l_vnlin
logical :: l_test

contains

subroutine get_ns_and_ncm()

use error, only : neo_abort
implicit none

namelist /control/ ns, eps, nreg, sigma1, sigma2,&
& sigma3, sigma4, nleg, nenergy, ncof, ic
!neogeo, neofrc
namelist /geometry/ isel, ishot, rho, e, q, Rn, Bn, eparr
namelist /species/ mas, temp, ncharge, Z, n, dp, dt 

integer :: io_stat, i,j 
integer :: ierr
integer :: sigma1, sigma2, sigma3, sigma4 !for input file
real, dimension(nionmax) :: z !for input file
real, dimension(nionmax) :: n, dp, dt
integer:: ncharge 
real :: mas, temp 

open(30,file='input.dat',FORM='formatted',STATUS='old', &
         POSITION='rewind', ACTION='read', IOSTAT=io_stat)
 if(io_stat /= 0)call neo_abort('input.dat not found!')
ns=0 
ncm=0
read(30,NML=control,iostat=io_stat)
do i=1, ns

 read(30,NML=species,iostat=io_stat) 
 if(ncharge>ncm) then
  ncm=ncharge
 end if

end do
close(unit=30)

if(ns ==0)call neo_abort('Set number of species and fill Species namelist')
if(ncm ==0)call neo_abort('At least one charge state per species')


end subroutine

subroutine set_defaults()
 
 implicit none
 !default
 l_test=.false.
 
 eps=1E-5
 nreg=1
 sigma=1
 nleg=3
 nenergy=1
 ncof=1
 ic=1
 l_vnlin = .false.
 
 isel=2
 rho=0.05
 e=0.19
 q=1.4
 Rn=6.
 Bn=5.3
 Eparr = 0.

 m=0.
 T=0.
 zsp=0.
 den=0.
 ds=0.
 
 vnlin_order = 0
 vnlin_contribution = 0

end subroutine


subroutine read_input()

use error, only : neo_abort
implicit none

namelist /control/ ns, eps, nreg, sigma1, sigma2,&
& sigma3, sigma4, nleg, nenergy, ncof, ic, l_vnlin, l_test

!neogeo, neofrc
namelist /geometry/ isel, ishot, rho, e, q, Rn, Bn, eparr
namelist /species/ mas, temp, ncharge, Z, n, dp, dt 
namelist /vnlin/ vnlin_order, vnlin_contribution
integer :: io_stat, i,j 
integer :: ierr
integer :: sigma1, sigma2, sigma3, sigma4 !for input file
real, dimension(nionmax) :: z !for input file
real, dimension(nionmax) :: n, dp, dt
integer:: ncharge 
real :: mas, temp 


open(30,file='input.dat',FORM='formatted',STATUS='old', &
         POSITION='rewind', ACTION='read', IOSTAT=io_stat)
 if(io_stat /= 0)call neo_abort('input.dat not found!')
open(31, file='input.out', iostat=io_stat)
 if(io_stat /= 0) call neo_abort('Failed to open input.out!')

!CONTROL INPUT
read(30,NML=control,iostat=io_stat) 

sigma(1)=sigma1
sigma(2)=sigma2
sigma(3)=sigma3
sigma(4)=sigma4
write(31,*) '------------------------------------------------'
write(31,*) '&Control'
write(31,*) 'ns = ', ns
write(31,*) 'eps = ', eps   
write(31,*) 'nreg = ', nreg 
write(31,*) 'sigma1 = ', sigma(1)   
write(31,*) 'sigma2 = ', sigma(2)  
write(31,*) 'sigma3 = ', sigma(3)  
write(31,*) 'sigma4 = ', sigma(4)
write(31,*) 'nleg = ', nleg  
write(31,*) 'nenergy = ', nenergy
write(31,*) 'ncof = ', ncof
! write(31,*) 'neogeo = ', neogeo
! write(31,*) 'neofrc = ', neofrc
write(31,*) 'ic = ', ic
write(31,*) 'l_vnlin = ', l_vnlin
write(31,*) 'l_vnlin = ', l_test

 
!GEOM INPUT
read(30,NML=geometry, iostat=io_stat)

write(31,*) '------------------------------------------------'
write(31,*) '&Geometry'
write(31,*) 'isel = ', isel
write(31,*) 'ishot = ', ishot   
write(31,*) 'rho = ', rho 
write(31,*) 'e = ', e   
write(31,*) 'q = ', q  
write(31,*) 'Rn = ', Rn  
write(31,*) 'Bn = ', Bn
write(31,*) 'Eparr = ', Eparr   


!SPECIES INPUT
do i=1, ns

read(30,NML=species, iostat=io_stat)
m(i)=mas*1.6726219E-27
T(i)=temp
nc(i)=ncharge
write(31,*) '------------------------------------------------'
write(31,*) '&Species'
write(31,*) 'nc= ', nc(i)
write(31,*) 'mas= ', mas
write(31,*) 'temp = ', temp

write(31, '(A)', advance='NO') ' z = '
do j=1, ncharge
zsp(i,j)=z(j)
write(31, '(F5.2)', advance='NO') zsp(i,j)
end do
write(31,*) ' '

write(31, '(A)', advance='NO') ' n = '
do j=1, ncharge
den(i,j)=n(j)
write(31, '(F5.2)', advance='NO') den(i,j)
end do
write(31,*) ' '

write(31, '(A)', advance='NO') ' dp = '
do j=1, ncharge
ds(i,j,1)=dp(j)
write(31, '(F5.2)', advance='NO') ds(i,j,1)
end do
write(31,*) ' '

write(31, '(A)', advance='NO') ' dt = '
do j=1, ncharge
ds(i,j,2)=dt(j)
write(31, '(F5.2)', advance='NO') ds(i,j,2)
end do
write(31,*) ' '

end do

!vnlin_control
read(30,NML=vnlin, iostat=io_stat)
write(31,*) '------------------------------------------------'
write(31,*) '&vnlin'
write(31,*) 'vnlin_order= ', vnlin_order
write(31,*) 'vnlin_contribution= ', vnlin_contribution

close(unit=30)
close(unit=31)



end subroutine

subroutine vnlin()

real :: vnlin_in(3,3)
real :: dum, dum1
integer :: k,l, io_stat
real :: rs, norm, vth


open(31, file='vnlin_moments_AV.dat', iostat=io_stat)
 if(io_stat /= 0) then
   write(*,*) ' supply vnlin_moments_AV.dat   &
                &  or set vnlin_Source=.false. '
   stop 1
 end if 
 
 do k=1,3
   do l=1,3
     read(31,*) vnlin_in(k,l), dum, dum1
   end  do
 end do
 

 vnlin_drive=.0
 
 if(l_vnlin) then
 ds=ds/(rn)     !! dirty solution to get same gradients as in gkw

 do k=1,3
   if(vnlin_order.ge.k) then
     if(vnlin_contribution.eq.0 .or. vnlin_contribution.eq.1) &
       & vnlin_drive(2,1,k) = vnlin_drive(2,1,k) + vnlin_in(1,k)
     if(vnlin_contribution.eq.0 .or. vnlin_contribution.eq.2) &
       & vnlin_drive(2,1,k) = vnlin_drive(2,1,k) + vnlin_in(2,k)
     if(vnlin_contribution.eq.0 .or. vnlin_contribution.eq.3) &
       & vnlin_drive(2,1,k) = vnlin_drive(2,1,k) + vnlin_in(3,k)
     end if
 end do
 
 if(vnlin_order .eq. -1) then
     if(vnlin_contribution.eq.0 .or. vnlin_contribution.eq.1) &
       & vnlin_drive(2,1,1) = vnlin_drive(2,1,1) + vnlin_in(1,1)
     if(vnlin_contribution.eq.0 .or. vnlin_contribution.eq.2) &
       & vnlin_drive(2,1,1) = vnlin_drive(2,1,1) + vnlin_in(2,1)
     if(vnlin_contribution.eq.0 .or. vnlin_contribution.eq.3) &
       & vnlin_drive(2,1,1) = vnlin_drive(2,1,1) + vnlin_in(3,1)
   end if
 if(vnlin_order .eq. -2) then
     if(vnlin_contribution.eq.0 .or. vnlin_contribution.eq.1) &
       & vnlin_drive(2,1,2) = vnlin_drive(2,1,2) + vnlin_in(1,2)
     if(vnlin_contribution.eq.0 .or. vnlin_contribution.eq.2) &
       & vnlin_drive(2,1,2) = vnlin_drive(2,1,2) + vnlin_in(2,2)
     if(vnlin_contribution.eq.0 .or. vnlin_contribution.eq.3) &
       & vnlin_drive(2,1,2) = vnlin_drive(2,1,2) + vnlin_in(3,2)
   end if
 if(vnlin_order .eq. -3) then
     if(vnlin_contribution.eq.0 .or. vnlin_contribution.eq.1) &
       & vnlin_drive(2,1,3) = vnlin_drive(2,1,3) + vnlin_in(1,3)
     if(vnlin_contribution.eq.0 .or. vnlin_contribution.eq.2) &
       & vnlin_drive(2,1,3) = vnlin_drive(2,1,3) + vnlin_in(2,3)
     if(vnlin_contribution.eq.0 .or. vnlin_contribution.eq.3) &
       & vnlin_drive(2,1,3) = vnlin_drive(2,1,3) + vnlin_in(3,3)
   end if
!rhostar
rs=sqrt(2*m(2)*T(2)*1.6022E-16)/(zsp(2,1)*1.6022E-19*Bn*Rn)
norm= rs*rs*Bn*den(2,1)*1E19*T(2)*1.6022E-16/Rn   ! 
!write(*,*) rs, norm
do k=1,3
vnlin_drive(2,1,k)=vnlin_drive(2,1,k)* norm
end do
!write(*,*) vnlin_drive
end if

end subroutine


subroutine check()

 use error, only : neo_abort, neo_warn
 implicit none
 
 integer :: i,j,k
 real :: soSmall
 
 soSmall= 1.E-15
 
 if ( nreg .ne. 1) call neo_warn('only nreg=1 (banana regime) tested, or switched off')
 if ( ic .ne. 1) call neo_warn('only ic=1 (banana contribution) tested, or switched off')
 if ( isel .ne. 2) call neo_warn('only isel=2 (circular geometry) tested, or switched off')
 
 if ( nc(1).ne.1 .or. zsp(1,1).ne.-1 ) call neo_abort('species 1 are the electrons')
 
 do i= 1,ns
   if ( m(i) == 0.) call neo_abort('set mass of all species')
   if ( T(i) == 0.) call neo_abort('set temperautre of all species')
   do j=1,nc(i)
     if ( abs(zsp(i,j)) < soSmall) call neo_abort('set charge number for all charge states')
     if ( abs(den(i,j)) < soSmall) call neo_abort('set density for all charge states')
     if ((.not.l_test).and. abs(ds(i,j,1)) < soSmall) call neo_abort('set pressure gradient for all charge states')
     if ((.not.l_test).and. abs(ds(i,j,2)) < soSmall) call neo_abort('set temperautre gradient for all charge states')
   end do
 end do

end subroutine

end module init