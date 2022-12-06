    subroutine load_st(stfile,rx,ry,rz,stnm,nst)
    ! Input
    ! stfile : a file name which contain the station information
    ! Output
    ! rx : x coordinate of stations
    ! ry : y coordinate of stations
    ! rz : z coordinate of stations
    ! vmodel : velocity model
    
    implicit none
    
    character(len=*) :: stfile
    real*8, allocatable :: rx(:),ry(:),rz(:)
    character*30,allocatable :: stnm(:)
    
    integer :: nst, ist
    
    character*80 :: buffer
    integer :: status
    
    open(10,file=stfile)
    nst = 0
    do while(.true.)
        read(10,fmt='(a79)',iostat=status) buffer
        if(status /= 0) exit
        nst = nst+1
    end do
    close(10)
    
    allocate(rx(nst))
    allocate(ry(nst))
    allocate(rz(nst))
    allocate(stnm(nst))
    open(10,file=stfile)
    do ist=1,nst
        read(10,*) stnm(ist),rx(ist),ry(ist),rz(ist)
    end do
    close(10)
    
    end subroutine load_st