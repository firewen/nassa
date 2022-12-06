    subroutine get_par(filename,lb,hb,dh,stfile,velfile,iflag)
    ! Input: 
    ! filename  : the file name wich includes all inversion parameters
    ! Output:
    ! lb : the lower range of (x,y,z)
    ! hb : the upper range of (x,y,z)
    ! dh : the step of (x,y,z)
    ! stfile : a file name which contain the station information
    ! velfile : a file name which contain the velocity model

    implicit none
    
    character(len=*) :: filename,stfile,velfile
    real*8 :: lb(3),hb(3),dh(3)
    integer :: iflag
    
    character*80 :: list,tmp
    
    open(10,file=filename,form='formatted',status='old')
    read(10,'(a80)') list
    read(10,'(a80)') list
    read(10,*) lb(1),hb(1),dh(1)
    read(10,'(a80)') list
    read(10,*) lb(2),hb(2),dh(2)
    read(10,'(a80)') list
    read(10,*) lb(3),hb(3),dh(3)
    read(10,'(a80)') list
    read(10,*) stfile
    read(10,'(a80)') list
    read(10,*) velfile
    !read(10,'(a80)') list
    !read(10,*) iflag
    
    close(10)
    
    return
    end subroutine get_par
