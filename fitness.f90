    function fitness(source,rx,ry,rz,nst,nt,dt,p,vmodel)
    use raytracing
    implicit none
    
    integer :: nst,nt
    real*8 :: source(3), rx(nst), ry(nst), rz(nst), p(nst,nt)
    real*8 :: dt
    type(Velocity) :: vmodel
    
    real*8 :: p0(nt), sump(nt)
    
    real*8 :: dk,dd,daze,dazs
    real*8 :: fp,ftime
    
    integer :: ist
    
    real*8 :: fitness
    
    sump = 0.0
    do ist=1,nst
        !call distaz(ry(ist),rx(ist),source(2),source(1),dk,dd,daze,dazs)
        call spheredist(ry(ist),rx(ist),source(2),source(1),dk)
        call ray_out(vmodel,dk,source(3),'P',fp,ftime)
        !print *,dk, source(3), fp, ftime
        call vshift(p(ist,:),nt,int(ftime/dt),p0)
        sump = sump+p0
    end do
    
    fitness = -1*maxval(sump)
    
    end function fitness
    
    subroutine vshift(p,nt,ishift,p0)
    implicit none
    
    integer :: nt, ishift
    real*8 :: p(nt),p0(nt)
    
    p0 = 0.0
    p0(1:nt-ishift) = p(ishift+1:nt)
    
    end subroutine vshift
