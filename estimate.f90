    subroutine estimate(ddir,grid,vm,np,nd,dist)
    ! Input
    ! grid : an array with size np * 3 contains the possible source locations
    ! vm : the brightness
    ! np : the number of the source locations
    ! Output
    ! dist : the distance from each source location to result
    
    use m_npy
    implicit none
    include 'omp_lib.h'
    
    character(len=40) :: ddir
    real*8 :: grid(np,nd), vm(np),dist(np)
    integer :: np,nd
    integer :: idx0(1), idx
    
    real*8 :: dk,dd,daze,dazs
    real*8 :: dist_min, dist_max, ddist, psum
    real*8, allocatable :: drange(:),p(:),v(:)
    real*8 :: p0(np)
    integer :: ndrange
    integer :: i,j
    
    real*8 :: t1, t2
    
    idx0 = maxloc(vm)
    idx = idx0(1)
    !print *, grid(idx,:),vm(idx)
    
    !psum = 0.0
    t1 = omp_get_wtime()
    !$omp parallel default(shared) private(i,dk,dd,daze,dazs)
    !$omp do
    do i=1,np
        !call distaz(grid(i,2),grid(i,1),grid(idx,2),grid(idx,1),dk,dd,daze,dazs)
        call spheredist(grid(i,2),grid(i,1),grid(idx,2),grid(idx,1),dk)
        dist(i) = sqrt(dk*dk+(grid(i,3)-grid(idx,3))*(grid(i,3)-grid(idx,3)))
        !psum = psum+exp(vm(i))
        p0(i) = exp(vm(i))
    end do
    !$omp end do nowait
    !$omp end parallel
    t2 = omp_get_wtime()
    print *,"time for distance: ",t2-t1 
    
    psum = sum(p0)
    dist_min = minval(dist)
    dist_max = maxval(dist)
    ndrange = 101
    ddist = (dist_max-dist_min)/(ndrange-1)
    allocate(drange(ndrange))
    allocate(p(ndrange))
    allocate(v(ndrange))
    do i=1,ndrange
        drange(i) = dist_min + (i-1)*ddist
    end do
    do i=1,ndrange
        v(i) = 0.0
        p(i) = 0.0
        do j=1,np
            if (dist(j) < drange(i)) then
                v(i) = v(i) + 1
                !p(i) = p(i) + exp(vm(j))
                p(i) = p(i) + p0(j)
            end if
        end do
        !print *, i,p(i),psum
        v(i) = v(i)/np
        p(i) = p(i)/psum
    end do
    
    open(10,file=trim(ddir)//'pv.txt')
    do i=1,ndrange
        write(10,*)drange(i),p(i),v(i)
    end do
    close(10)
! save brightness in numpy format
    call add_npz(trim(ddir)//'bright.npz','grid',grid)
    call add_npz(trim(ddir)//'bright.npz','vm',vm)
    call add_npz(trim(ddir)//'bright.npz','dist',dist)
    
    idx0 = maxloc(dist)
    !call distaz(grid(idx0(1),2),grid(idx0(1),1),grid(idx,2),grid(idx,1),dk,dd,daze,dazs)
    call spheredist(grid(idx0(1),2),grid(idx0(1),1),grid(idx,2),grid(idx,1),dk)
    print *, "max. dist.", grid(idx0(1),:),grid(idx,:),dk,dist(idx0(1))
    
    deallocate(drange)
    deallocate(p)
    deallocate(v)
    
    end subroutine estimate
