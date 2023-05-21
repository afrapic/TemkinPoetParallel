module debug_module
    logical, parameter :: debug = .true.
    integer, parameter :: stdout = 6
    integer :: iprint = 0 
    integer :: nout = 6

contains

    subroutine set_debuglevel(level,iam)
        implicit none
        character level
        integer iam
        character*2 itxt
        !
        !   Default: processor 0 writes on standard out
        !
        nout = 6
        !
        !   Each processor writes in a separate output
        !
        if(level.eq.'a'.or.level.eq.'A') then
            iprint = iam
            nout = 3000+iam
            write(itxt,'(i2.2)') iam
            open(unit=nout,file='tp'//itxt//'.dat')
            !
            !   Null printing, except for debug
            !
        else if(level.eq.'n'.or.level.eq.'N') then
            iprint = -1
        end if

        return
    end subroutine set_debuglevel

    subroutine timing(timei,iam,nproc,comment)
        !
        ! Prints timing data, returns actual time in timei
        !
        implicit none
        integer iam,nproc
        real timef,timei,time        
        character*(*) comment
        timef = 0.0
        if (iam.eq.iprint) then
            call cpu_time(timef)
            time=timef-timei
            write(nout,'(a30,f9.3,a,i4)') comment,time/60.,"min., processors = ",nproc
        endif
        timei = timef
        return
    end subroutine timing

end module debug_module
