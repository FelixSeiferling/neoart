!-----------------------------------------------------------
!>    this subroutine calculates the classical contribution
!>    to the transport
!>
!>    input  nsm  : maximum number of species
!>           ncm  : maximum number of charge states
!>           ns   : actual number of species
!>           nc   : array(nsm) actual number of charge 
!>                  states per species
!>           t    : array(nsm) temperature of every species
!>           zsp  : array(nsm,ncm) the charge in units of e
!>                  of every charge state
!>           xi   : array(nsm,ncm) weigting factor for the 
!>                  reduced charge state method
!>           la   : array(3,3,nsm) test particle friction 
!>                  coefficients
!>           lab  : array(3,3,nsm,nsm) friction coefficients
!>           ds   : array(nsm,ncm,2) the thermodynamic
!>                  forces. for definition see routine 
!>                  neoart
!>
!>    output coeffc : the tranport coefficients (see neoart
!>                    for definition
!-----------------------------------------------------------
subroutine class(coeffc)
 
 
  use init, only :  ns,nc, ncm, zsp, m, t, den, ds
  use collision, only : xi
  use friction, only : la,lab
                       
  implicit none

  integer :: i,j,k,l
  real :: coeffc,v

  dimension coeffc(ns,ncm,4)

  do i = 1, ns ;  do j = 1, nc(i)
          
    coeffc(i,j,1) = 0.
    coeffc(i,j,2) = 0.
 
    do k = 1, ns ; do l = 1, nc(k)

      v = t(k)/zsp(k,l)
     coeffc(i,j,1) = coeffc(i,j,1) + v*xi(i,j)*               &
          &  xi(k,l)*(lab(1,1,i,k)*ds(k,l,1)+lab(1,2,i,k)*    &
          &  ds(k,l,2)) 
     coeffc(i,j,2) = coeffc(i,j,2) + v*xi(i,j)*               &
          &  xi(k,l)*(lab(2,1,i,k)*ds(k,l,1)+lab(2,2,i,k)*    &
          &  ds(k,l,2)) 
     if ((k.eq.i).and.(l.eq.j)) then 
       coeffc(i,j,1) = coeffc(i,j,1) + v*xi(i,j)*(            &
          &   la(1,1,i)*ds(i,j,1)+la(1,2,i)*ds(i,j,2))
       coeffc(i,j,2) = coeffc(i,j,2) + v*xi(i,j)*(            &
          &   la(2,1,i)*ds(i,j,1)+la(2,2,i)*ds(i,j,2))
     endif

    end do ; end do 
  end do ; end do

return
end subroutine class
