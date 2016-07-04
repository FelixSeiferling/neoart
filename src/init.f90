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

integer, dimension(4) :: sigma  !THE SIGMA'S DETERMINE WHICH OF THE COUPLING 
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
logical :: neogeo !      LOGICAL IF TRUE THEN THE GEOMETRY DEP.
!                        PARAMETERS ARE RECALCULATED.
logical :: neofrc !      LOGICAL IF TRUE THEN THE FRICTION AND 
!                       VISCOSITY MATRIX ARE !NOT! NEWLY CALCULATED
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
real, save, allocatable :: T(:) !DIM:(ns) temperautre in units of kEV (OR THERMAL VEL???)
real, save, allocatable :: den(:,:) !DIM:(ns,nc) density in units of 10^19 m^-3 (OR THERMAL VEL???)
real, save, allocatable :: ds(:,:,:) !DIM:(ns,nc,2); 1 Pressure gradient; 2 Temperature gradient
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

!integer :: NAR. NZM

contains

subroutine read_input()

implicit none

namelist /control/ ns, eps, nreg, sigma1, sigma2,&
& sigma3, sigma4, nleg, nenergy, ncof, neogeo, neofrc, ic
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
ncm=0
read(30,NML=control,iostat=io_stat)
do i=1, ns

 read(30,NML=species,iostat=io_stat) 
 if(ncharge>ncm) then
  ncm=ncharge
 end if

end do
close(unit=30)

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
write(31,*) 'neogeo = ', neogeo
write(31,*) 'neofrc = ', neofrc
write(31,*) 'ic = ', ic

allocate(nc(ns),stat=ierr)
 if (ierr /= 0) call neo_abort('could not allocate charge state array')
allocate(zsp(ns,ncm),stat=ierr)
 if (ierr /= 0) call neo_abort('could not allocate charge number array')
allocate(m(ns),stat=ierr)
 if (ierr /= 0) call neo_abort('could not allocate mass array')
allocate(T(ns),stat=ierr)
 if (ierr /= 0) call neo_abort('could not allocate temperature array')
allocate(den(ns,ncm),stat=ierr)
 if (ierr /= 0) call neo_abort('could not allocate density array')
allocate(ds(ns,ncm,2),stat=ierr)
 if (ierr /= 0) call neo_abort('could not allocate mass array')
 
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
m(i)=mas
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
ds(i,j,1)=dt(j)
write(31, '(F5.2)', advance='NO') ds(i,j,2)
end do
write(31,*) ' '

end do
close(unit=30)
close(unit=31)



end subroutine

 subroutine neo_abort(message)

  character (len=*),intent(in) :: message

  write(*,*) 'NEOART STOPPED: ', message
  stop 1

 end subroutine


end module init