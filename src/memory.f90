module memory

contains
 subroutine alloc_all()
 
 use init, only : ns,ncm,cff
 use init, only : nc,zsp,m,t,den,ds,sigma
 use collision, only : xi,tau
 use friction, only : la,lab
 use geometry, only : fm, nmaxgr
 use error, only : neo_abort
 implicit none
 integer :: ierr
 
  allocate(nc(ns),stat=ierr)
    if (ierr /= 0) call neo_abort('could not allocate charge state array')
    nc=0
  allocate(zsp(ns,ncm),stat=ierr)
    if (ierr /= 0) call neo_abort('could not allocate charge number array')
    zsp=0.
  allocate(m(ns),stat=ierr)
    if (ierr /= 0) call neo_abort('could not allocate mass array')
    m=0.
  allocate(T(ns),stat=ierr)
    if (ierr /= 0) call neo_abort('could not allocate temperature array')
    T=0.
  allocate(den(ns,ncm),stat=ierr)
    if (ierr /= 0) call neo_abort('could not allocate density array')
    den=0.
  allocate(ds(ns,ncm,2),stat=ierr)
    if (ierr /= 0) call neo_abort('could not allocate mass array')
    ds=0.
  allocate(sigma(4),stat=ierr)
    if (ierr /= 0) call neo_abort('could not allocate mass array')
    sigma=0
 
  allocate(XI(ns,ncm),stat=ierr)
    if (ierr /= 0) call neo_abort('could not allocate array')
    xi=0.
  allocate(TAU(ns,ncm),stat=ierr)
    if (ierr /= 0) call neo_abort('could not allocate array')
    tau=0.
    
  allocate(la(3,3,ns),stat=ierr)
    if (ierr /= 0) call neo_abort('could not allocate array')
    la=0.
  allocate(lab(3,3,ns,ns),stat=ierr)
    if (ierr /= 0) call neo_abort('could not allocate array')
    lab=0.
    
  allocate(fm(nmaxgr),stat=ierr)
    if (ierr /= 0) call neo_abort('could not allocate array')
    fm=0.
 
  allocate(cff(ns,ncm,4),stat=ierr)
    if (ierr /= 0) call neo_abort('could not allocate array')
    cff=0.
    
 end subroutine

 subroutine dealloc_all()
 
 use init, only : nc,zsp,m,t,den,ds,sigma,cff
 use collision, only : xi,tau
 use friction, only : la,lab
 use geometry, only : fm, nmaxgr
 implicit none
 
 
  deallocate(nc)
  deallocate(zsp)
  deallocate(m)
  deallocate(T)
  deallocate(den)
  deallocate(ds)
  deallocate(sigma)
 
  deallocate(xi)
  deallocate(tau)
    
  deallocate(la)
  deallocate(lab)
  
  deallocate(fm)
  
  deallocate(cff)
    
end subroutine
end module