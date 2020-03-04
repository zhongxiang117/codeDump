      program testREWIND

      character*80, envalue, tmp
      integer*4     i

      call getenv('TESTFILE',envalue)
      print*, envalue

      if (envalue.eq.'HEAD  ')  then

        open(UNIT=999,FILE=envalue,ACCESS='APPEND')
        do i=1,10
            write(999,'(I4)') i
        enddo

        close(999)
      endif


      open(UNIT=999,FILE=envalue)

      read(999,'(A20)') tmp
      print*, tmp

      read(999,'(A20)') tmp
      print*, tmp

      read(999,'(A20)') tmp
      print*, tmp


      rewind 999
      read(999,'(A20)') tmp
      print*, tmp

      read(999,'(A20)') tmp
      print*, tmp

      read(999,'(A20)') tmp
      print*, tmp


 
      rewind 999
      read(999,'(A20)') tmp
      print*, tmp

      end
