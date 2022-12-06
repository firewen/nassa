     program main
    
    use raytracing
    use omp_lib
    implicit none
    
    integer :: nst
    character*80 :: buffer
    integer :: status
    
    character(len=40) :: filename,stfile,velfile,ddir
    real*8 :: lb(3),hb(3),dh(3)
    real*8,allocatable :: rx(:),ry(:),rz(:)
    character*30,allocatable :: stnm(:)
    type(Velocity) :: vmodel
    
    integer :: iflag, nx,ny,nz,np
    integer :: i,j
    character*80 :: fnm
    
    integer :: ns,nr,niter,nd
    real*8,allocatable :: vm(:),vf(:,:)
    real*8,allocatable :: dist(:)
    
    real*8, allocatable :: pall(:,:)
    integer :: nt
    real*8 :: dt
    integer :: idx0(1), idx
    
    real*8 :: twin1, twin2
    integer :: itw1, itw2, llen
    real*8, allocatable :: pp(:,:)
    integer :: evid
    
    character*20 :: tmp
 
    interface
        subroutine load_st(stfile,rx,ry,rz,stnm,nst)
            character(len=*) :: stfile
            real*8, allocatable :: rx(:),ry(:),rz(:)
            character*30,allocatable :: stnm(:)
            integer :: nst
        end subroutine load_st
        
    end interface
    
    call getarg(1,ddir)
    call getarg(2,tmp)
    read(tmp,'(f)') twin1 
    call getarg(3,tmp)
    read(tmp,'(f)') twin2 
    call getarg(4,tmp)
    read(tmp,'(i)') nt 
    call getarg(5,tmp)
    read(tmp,'(f)') dt 
    call getarg(6,tmp)
    read(tmp,'(i)') evid 
!    print *, ddir,"twin1",twin1,nt,dt

    !ddir = '2014.08.19.12.41.34/'
    
!    nt = 50001
!    dt = 0.01
    
    !twin1 = 97.59
    !twin2 = 123.95
    
    filename=trim(ddir) // 'input.txt'
    call get_par(filename,lb,hb,dh,stfile,velfile,iflag)
    
    call load_st(trim(ddir) // stfile,rx,ry,rz,stnm,nst)

    !print *, lb, stfile, rx(1),stnm(1)

    itw1 = int(twin1/dt)
    itw2 = int(twin2/dt)
    llen = itw2-itw1+1
    allocate(pp(nst,llen))
    pp = 0.0
    do i=1,nst
        fnm = trim(ddir)//'probs/'//trim(stnm(i))//'.prob'
        open(10,file=fnm,recl=25,form='formatted',access='direct',status='old')
        do j=1,llen
            read(10,fmt="(e24.17)",rec=itw1+j)pp(i,j)
        end do
        close(10)
    end do
    print *, "probs have been loaded"
    print *,shape(pp)
    nt = llen
        
!    pp = pall(:,itw1:itw2)
!    print *,shape(pp)
!    nt = llen
!
!    allocate(pall(nst,nt))
!    pall = 0.0
!    do i=1,nst
!        fnm = trim(ddir)//'probs/'//trim(stnm(i))//'.prob'
!        open(10,file=fnm)
!        do j=1,nt
!            read(10,*)pall(i,j)
!        end do
!        close(10)
!    end do
    
    !call omp_set_num_threads(16)
    print *, 'procs',omp_get_num_procs()
    print *, 'threads',omp_get_num_threads() 
    
    ns = 100
    nr = 100
    niter = 50
    nd = 3
    np = ns+ns*niter
    call model_info(trim(ddir) // velfile,vmodel,0d0)
    
    allocate(vf(np,nd))
    allocate(vm(np))
    call NA(rx,ry,rz,nst,vmodel,pp,nt,ns,nr,niter,nd,lb,hb,dt,vm,vf)
    
    idx0=minloc(vm)
    print *, idx0(1),vf(idx0(1),:),vm(idx0(1))
    
    !open(10,file='out.txt')
    !do i=1,ns+ns*niter
    !    write(10,'(4(f8.4,1x))')(vf(i,j),j=1,nd),vm(i)
    !    !write(*,*)(vf(i,j),j=1,nd),vm(i)
    !end do
    !close(10)
    
    write(fnm,'(i)')evid
    fnm = trim(ddir)//'outfig/na'//trim(adjustl(fnm))
    print *, fnm
    allocate(dist(np))
    call estimate(fnm,vf,-vm,np,nd,dist)
    
    open(10,file=trim(fnm)//'bright.txt')
    do i=1,np
        write(10,'(5(f8.4,1x))')(vf(i,j),j=1,3),-vm(i),dist(i)
    end do
    close(10)
    
    
!    deallocate(pall)
    deallocate(pp)
    deallocate(vf)
    deallocate(vm)
    deallocate(dist)
    deallocate(rx)    
    deallocate(ry)    
    deallocate(rz)    
    deallocate(stnm)    

    end program main
