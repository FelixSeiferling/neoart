! C--------------------------------------------------------------------
! C     THE ROUTINE THAT CALCULATES THE PFIRSCH SCHLUETER CONTRIBUTION.
! C
! C     INPUT   NSC    : MAXIMUM NUMBER OF SPECIES (USED TO CHECK 
! C                      CONSISTENCY
! C             NCC    : MAXIMUM NUMBER OF CHARGE STATES (USED TO 
! C                      CHECK CONSISTENCY
! C             NS     : NUMBER OF SPECIES
! C             NC     : ARRAY(NS) NUMBER OF CHARGE STATES PER SPECIES
! C             NEOFRC : LOGICAL IF TRUE THEN THE FRICTION MATRICES
! C                      ARE NOT NEWLY CALCULATED
! C             TAU    : ARRAY(NSC,NSC) NORMALIZED COLLISION FREQUEN-
! C                      CIES (M_A N_A / TAU_AB)
! C             XI     : ARRAY(NSC,NCC) RELATIVE WEIGHT OF THE CHARGE
! C                      STATES
! C             LA     : FRICTION COEFFICIENTS
! C             LAB    : FRICTION COEFFICIENTS
! C             m    : ARRAY(NSC) MASS OF THE SPECIES (KG)
! C             T   : ARRAY(NSC) TEMPERATURE IN KEV
! C             DEN    : ARRAY(NSC,NCC) DENSITY IN 10^19 M^-3
! C             ZSP    : CHARGE NORMALIZED TO E, ARRAY(NSM,NCM) 
! C             SIGMA  : CONTROL PARAMETERS THAT DETERMINE WHICH OF 
! C                      THE COUPLING TERMS ARE TAKEN INTO ACCOUNT. 
! C                      FOR ALL COUPLING TERMS ALL SIGMAS ARE 1
! C                      SEE ROUTINE NEOART FOR DEFINITION
! C             ds  : THE THERMODYNAMIC FORCES. ARRAY(NSM,NCM,3)
! C                      SEE ROUTINE NEOART FOR DEFINITION
! C             RNQ    : LENGTH OF THE FIELD LINE. IN CIRCULAR GEOMETRY
! C                      THIS CAN BE APPROXIMATED BY MAJOR RADIUS 
! C                      TIMES SAFETY FACTOR. 
! C     OUTPUT  UAI    : THE VELOCITIES
! C             COEFF  : NORMALIZED TRANSPORT COEFFICIENTS. 
! C
! C     THE ROUTINE CALLS THE FOLLOWING SUBROUTINES
! C     PERR   : ERROR HANDLING
! C     PENQ   : TO CALCULATE THE COUPLING COEFFICIENTS
! C     LUDCMP : TO CALCULATE THE LU DECOMPOSITION
! C     LUBKSB : TO CALCULATE THE SOLUTION OF THE MATRIX EQUATION   
! C--------------------------------------------------------------------
subroutine PS(uai,coeff)
      
      
      use init, only : ns,nc,ncm,ds
      use init, only : m,t,den,zsp,sigma
      use collision, only : tau, xi
      use friction, only : la,lab
      use geometry, only : rnq
      IMPLICIT NONE
 

       integer :: i,j,k,l,ii,ij,ik,indx
       integer :: il, im, in
       real :: a1,a2,a3
       real :: lg, ta,tb,ma,mb
       real :: pp, qq, pab, qab
       real :: dum1, dum2
       real :: oa, kab, ha, da, rm
       real :: dnom, lga,gk, rmoi
       real :: rmo, rmoo, rmooo
       real :: dv, ev, fv
       real :: aaa, akeep, d
       real :: uo, uoo
       real :: sor, soro, soroo
       real :: b, bkeep, btus
       real :: sol, usolo
       real :: f1,f2
      logical :: nliter, nlchan
      real, dimension(ns,ncm,3) :: uai, uaio
      real, dimension(ns,ncm,4) :: coeff
      dimension pab(2,2,ns,ns),qab(2,2,ns,ns), &
       &          kab(2,2,ns,ns), oa(2,2,ns), &
       &	  pp(2,2),qq(2,2),ha(2,3,ns), &
       &          rm(2,2,ns,ncm),gk(ns), rmoi(2,2,ns), & 
       &          rmo(2,2,ns),rmoo(2,2,ns),rmooo(2,2,ns), &
       &          dv(2,3,ns,ns),ev(2,3,ns,ns),fv(2,3,ns,ns), &
       &          aaa(2*ns,2*ns), akeep(2*ns,2*ns), &
       &          uo(ns,3),uoo(ns,3), &
       &          sor(2,ns,ncm),soro(2,ns), soroo(2,ns), &
       &          b(2*ns),bkeep(2*ns),btus(2*ns), &
       &          sol(2*ns), usolo(ns,3)


      coeff=0.
      uai=0.
       nliter = .FALSE.
       nlchan = .false.
            
!     SOME COEFFICIENTS THAT APPEAR IN THE EQUATIONS
      a1 = 4./25.
      a2 = 16./175.
      a3 = (8./35.)**2



!     SINCE THE THERMODYNAMIC FORCES ARE PROPORTIONAL TO T/ Z 
!     THE UAIO QUANTITIES ARE MULTIPLIED BY THIS NUMBER.
      do i = 1, ns
        do j = 1, nc(i) 
          do k = 1, 2
            UAIO(I,J,K) = T(I)*ds(I,J,K) /(ZSP(I,J))
          end do
          UAIO(I,J,3) = 0.
      end do; end do


!     THE COEFFICIENT THAT DETERMINES THE INFLUENCE OF THE FIELD
!     LINE LENGTH
      LG = -1./(RNQ)**2

!     CALCULATED THE FRICTION COEFFICIENTS
      do i = 1, ns
        do j = 1, ns 
          TA = T(I)
          TB = T(J)
          MA = m(I)
          MB = m(J)       
          CALL PENQ(TA,TB,MA,MB,PP,QQ)
          do k = 1, 2
            do l = 1, 2
              PAB(K,L,I,J) = PP(K,L)
              QAB(K,L,I,J) = QQ(K,L)
            end do
          end do
        end do
      end do

!     calculate the kab and oa matrices
      do i = 1, ns
        dum1 = 1./(t(i))
        oa(1,1,i) = 0.
        oa(1,2,i) = 0.
        oa(2,1,i) = 0.
        oa(2,2,i) = 0.
        do j = 1, ns
          dum2 = 1./(t(j))
          kab(1,1,i,j) = -a1*dum1*sigma(3)*tau(i,j)* &
     &                   pab(1,1,i,j)
          kab(1,2,i,j) = - a2*dum2*sigma(1)*tau(i,j)* &
     &                   qab(1,2,i,j)
          kab(2,1,i,j) = a2*dum1*sigma(2)*tau(i,j)* &
     &                   (pab(2,1,i,j)+pab(1,1,i,j))
          kab(2,2,i,j) = a3*dum2*tau(i,j)*sigma(4)* &
     &                   (qab(2,2,i,j)+sigma(1)*qab(1,2,i,j))
          oa(1,1,i) = oa(1,1,i) + a1*dum1*sigma(3)* &
     &                tau(i,j)*pab(1,1,i,j)
          oa(1,2,i) = oa(1,2,i) - a2*dum1*sigma(1)* &
     &                tau(i,j)*pab(1,2,i,j)
          oa(2,1,i) = oa(2,1,i) - a2*dum1*sigma(2)* &
     &                tau(i,j)*(pab(2,1,i,j) + &
     &                pab(1,1,i,j))
          oa(2,2,i) = oa(2,2,i) + a3*dum1*tau(i,j)*sigma(4)* &
     &                (pab(2,2,i,j)+ sigma(1)*pab(1,2,i,j))
      end do; end do

!    now calculate the ha coefficients
      do i = 1, ns
        do k = 1, 2
          do l = 1, 3
            ha(k,l,i) = 0.
            do j = 1, 2
              ha(k,l,i) = ha(k,l,i) - oa(k,j,i)* la(j+1,l,i)
            end do
          end do
        end do
      end do

      do i = 1, ns
        da = ha(1,2,i)*ha(2,3,i)-ha(1,3,i)*ha(2,2,i)
        ta = ha(1,2,i)+ha(2,3,i)
        lga = LG*m(i)*1.6e22
        do j = 1, nc(i)
          dnom = (xi(i,j)/den(i,j))/ (LGA**2 + (xi(i,j)/ &
     &       (den(i,j)))**2*LGA*ta+(xi(i,j)/(den(i,j)))**4*da)
          dum1 = (xi(i,j)/den(i,j))**2
          rm(1,1,i,j) =  dnom*(lga+dum1*ha(2,3,i))
          rm(1,2,i,j) = -dnom*dum1*ha(1,3,i)
          rm(2,1,i,j) = -dnom*dum1*ha(2,2,i)
          rm(2,2,i,j) =  dnom*(lga+dum1*ha(1,2,i))
        end do
      end do
 
      do i = 1, ns
        gk(i) = 0.
        do k = 1, 2
          do l = 1, 2
          rmo(k,l,i) = 0.
          rmoo(k,l,i) = 0.
          rmooo(k,l,i) = 0.
          end do
        end do
        do j = 1, nc(i)
          dum1 = xi(i,j)**2/den(i,j)
          dum2 = xi(i,j)*dum1/den(i,j)
          gk(i) = gk(i) + dum1
          do k = 1, 2
            do l = 1, 2
            rmo(k,l,i) = rmo(k,l,i)+ xi(i,j)*rm(k,l,i,j)
            rmoo(k,l,i) = rmoo(k,l,i)+dum1*rm(k,l,i,j)
            rmooo(k,l,i) = rmooo(k,l,i)+dum2*rm(k,l,i,j)
            end do
          end do
        end do
      end do
 
      do i = 1, ns
        dnom = 1./(rmo(1,1,i)*rmo(2,2,i)- &
     &         rmo(1,2,i)*rmo(2,1,i)) 
        rmoi(1,1,i) =  dnom*rmo(2,2,i) 
        rmoi(1,2,i) = -dnom*rmo(1,2,i)
        rmoi(2,1,i) = -dnom*rmo(2,1,i)
        rmoi(2,2,i) =  dnom*rmo(1,1,i)
      end do

      do i = 1, ns
        do j = 1, ns

          do k = 1, 2
            do l = 1, 3
              dv(k,l,i,j) = 0.
              ev(k,l,i,j) = 0.
              fv(k,l,i,j) = 0.
              do ii = 1, 2
                dv(k,l,i,j) = dv(k,l,i,j) +  &
     &                        kab(k,ii,i,j)*la(ii+1,l,j)
                ev(k,l,i,j) = ev(k,l,i,j) +  &
     &                        oa(k,ii,i)*lab(ii+1,l,i,j)
                do ij = 1, ns
                  fv(k,l,i,j) = fv(k,l,i,j) +  &
     &               gk(ij)*kab(k,ii,i,ij)*lab(ii+1,l,ij,j)
                end do
              end do
            end do
          end do
        end do
      end do
 
!    build the matrix
       do i = 1, ns; do j = 1, ns; do k = 1, 2; do l = 1, 2
              dum1 = 0.
              dum2 = 0.
              do ii = 1, 2
                dum1 = dum1 + rmo(k,ii,i)*fv(ii,l+1,i,j) &
     &                 +rmoo(k,ii,i)*ev(ii,l+1,i,j)
                do ij = 1, 2
                  do ik = 1, 2
                     dum2 = dum2 + rmo(k,ii,i)*dv(ii,ij+1,i,j)* &
     &                 rmoo(ij,ik,j)*rmoi(ik,l,j)
                  end do
                end do
              end do

              aaa(2*(i-1)+k,2*(j-1)+l) = - dum1 - dum2

              dum1 = 0.
              do ii = 1, ns
                do ij = 1, 2
                  do ik = 1, 2 
                    do il = 1, 2
                      dum2 = 0.
                      do im = 1, 2
                        do in = 1, 2
                          dum2 = dum2 + rmoo(ik,im,ii)* &
     &                      rmoi(im,in,ii)*rmoo(in,il,ii)
                        end do
                      end do
                      dum2 = rmooo(ik,il,ii)-dum2
                      dum1 = dum1 + rmo(k,ij,i)*dv(ij,ik+1,i,ii)* &
     &                    dum2*ev(il,l+1,ii,j)     
                    end do
                  end do 
                end do
              end do
              aaa(2*(i-1)+k,2*(j-1)+l) = aaa(2*(i-1)+k,2*(j-1)+l) &
     &           - dum1
   end do; end do; end do; end do

      do i = 1, 2*ns
        aaa(i,i) = aaa(i,i) + 1.
      end do

      if (nlchan) then
      end if

      if (nliter) then
        do i = 1, 2*ns
          do j = 1, 2*ns
            akeep(i,j) = aaa(i,j)
          end do
        end do
      end if

!     now do the lU decomposition 
      call ludcmp(aaa,2*ns,2*ns,indx,d)

      do i = 1, ns
        do j = 1, 3
          uo(i,j) = 0.
          uoo(i,j) = 0.
          do k = 1, nc(i)
            uo(i,j) = uo(i,j) + xi(i,k)*uaio(i,k,j) 
            uoo(i,j)=uoo(i,j)+xi(i,k)**2*uaio(i,k,j)/ &
     &         (den(i,k))
          end do
        end do
      end do

!     calculate the right hand side solutions
      do i = 1, ns
        do j = 1, nc(i)
          do k = 1, 2
            sor(k,i,j) = 0.
            do l = 1, 3
              do ii = 1, 2
                sor(k,i,j) = sor(k,i,j) - (xi(i,j)/den(i,j))* &
     &            rm(k,ii,i,j)*ha(ii,l,i)*uaio(i,j,l)
                do ij = 1, ns
                  sor(k,i,j) = sor(k,i,j) + rm(k,ii,i,j)* &
     &            (fv(ii,l,i,ij)*uo(ij,l)+dv(ii,l,i,ij)*uoo(ij,l)+ &
     &            (xi(i,j)/den(i,j))*ev(ii,l,i,ij)*uo(ij,l))
                end do
              end do
            end do
          end do
        end do
      end do


      do i= 1, ns
        do k = 1, 2
          soro(k,i) = 0.
          soroo(k,i) = 0.
          do  j = 1, nc(i)
            soro(k,i) = soro(k,i) + xi(i,j)*sor(k,i,j)
            soroo(k,i) = soroo(k,i) + xi(i,j)**2*sor(k,i,j)/ &
     &                   den(i,j)
          end do
        end do
      end do

      do i = 1, ns
        do  k = 1, 2
          b(2*(i-1)+k) = soro(k,i)
          do j = 1, ns
            do l = 1, 2
              do ii = 1, 2
                b(2*(i-1)+k) = b(2*(i-1)+k) + rmo(k,l,i)* &
     &            dv(l,ii+1,i,j)*soroo(ii,j)
                do ij = 1, 2
                  do ik = 1, 2
                    b(2*(i-1)+k) = b(2*(i-1)+k) - rmo(k,l,i)* &
     &                dv(l,ii+1,i,j)*rmoo(ii,ij,j)* &
     &                rmoi(ij,ik,j)*soro(ik,j)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do

      if (nlchan) then
      end if

      if (nliter) then
        do i = 1, 2*ns
          bkeep(i) = b(i)
        end do
      end if
      
!      calculate the solution through backsubstition   
      call lubksb(aaa,2*ns,2*ns,indx,b)
 
      if (nliter) then
        do i = 1, 2*ns
          sol(i) = b(i)
        end do
        do k=1, 5
          do i = 1, 2*ns
            btus(i) = bkeep(i)
            do j = 1, 2*ns
              btus(i) = btus(i)-akeep(i,j)*sol(j)
            end do
          end do
          call lubksb(aaa,2*ns,2*ns,indx,btus)
          do i = 1, 2*ns
            sol(i) = sol(i)+btus(i)
          end do 
        end do
        do i = 1, 2*ns
          b(i) = sol(i)
        end do
      end if

!      now calculate the velocities
      do i = 1, ns
        do j = 1, nc(i)
          do k = 1, 2
            uai(i,j,k+1) = sor(k,i,j)
            do l = 1, 2
              do ii = 1, 2
                uai(i,j,k+1) = uai(i,j,k+1) + rm(k,l,i,j)* &
     &             rmoi(l,ii,i) * (b(2*(i-1)+ii) - soro(ii,i))
              end do
            end do
            do l = 1, 2
              do ii = 1, 2
                dum1 = 0.
                if (l.eq.ii) dum1 = dum1 + xi(i,j)/den(i,j)
                do ij = 1, 2
                  dum1 = dum1 - rmoi(l,ij,i)*rmoo(ij,ii,i)
                end do
                do ij = 1, 2 
                  do ik = 1, ns
                    uai(i,j,k+1) = uai(i,j,k+1) + rm(k,l,i,j)* &
     &                dum1*ev(ii,ij+1,i,ik)*b(2*(ik-1)+ij)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
                  
  
      do i = 1, ns
        do j = 1, nc(i)
          uai(i,j,1) = uaio(i,j,1)
          uai(i,j,2) = uai(i,j,2)+uaio(i,j,2)
          usolo(i,1) = uo(i,1)
          usolo(i,2) = uo(i,2) + b(2*(i-1)+1)
          usolo(i,3) = b(2*(i-1)+2)
        end do
      end do
 
!     calculate the friction 
      do i = 1, ns
        f1 = 0.
        f2 = 0.
        do j = 1, ns
          do k = 1, 3
            f1 = f1 + lab(1,k,i,j)*usolo(j,k)
            f2 = f2 + lab(2,k,i,j)*usolo(j,k)
          end do
        end do
        do j = 1, nc(i)
          coeff(i,j,1) = 0.
          coeff(i,j,2) = 0.
          do k = 1, 3
            coeff(i,j,1) = coeff(i,j,1) + xi(i,j)*(la(1,k,i)* &
     &                     uai(i,j,k))
            coeff(i,j,2) = coeff(i,j,2) + xi(i,j)*(la(2,k,i)* &
     &                     uai(i,j,k))
          end do
          coeff(i,j,1) = coeff(i,j,1) + xi(i,j)*f1
          coeff(i,j,2) = coeff(i,j,2) + xi(i,j)*f2
        end do
      end do
 
   return
end subroutine
