      program rmod

      character*80 line
c	open files for writing
      open(unit=3,file='mod.7')
      open(unit=5,file='mod.5')
      open(unit=15,file='flag')
      open(unit=9,file='old')

c	open file for temp and log g (must be in executed directory
      open(unit=2,file='new')

c	set the modellist
      open(unit=1,file='./models/ap00k2odfnew.dat')

c	set the extra line info file
      open(unit=4,file='./ref/xtra_line_info')


c	get temp and logg from file new (unit 2)
      read(2,*)xtt,xgg
       if(xtt.le.11000.) then
        write(15,*)'IOPADD=1 IRSCT=1 IOPHMI=1'
       else
        write(15,*)'IOPADD=0'
       endif

      iflag=0
      jflag=0
      i=0
 1     i=i+1
       read(1,1000,end=999)teff,xlogg
c      write(*,*)teff,xlogg,i
       if(teff.eq.xtt.and.xgg.eq.xlogg) then
        jflag=1
       endif
       if(jflag.eq.1) write(3,1000)teff,xlogg
       if(iflag.eq.0) then
        do 12 k=1,21
         read(1,1001)line
         if(jflag.eq.1) write(3,1001)line
 12     continue

c        read(1,1001)line
         read(1,1002)nd
         ne=69
        if(jflag.eq.1) write(3,1002)ne
         ne=74
c       if(jflag.eq.1) write(3,1001)line
        do 2 k=1,74
         read(1,1001)line
         if(jflag.eq.1) write(3,1001)line
 2      continue
       endif
       if(iflag.eq.1) then
        do 3 k=1,88
         read(1,1001)line
         if(jflag.eq.1) write(3,1001)line
 3      continue
       endif
       if(jflag.eq.2) then
        rewind(2)
        write(2,1000)teff,xlogg
        write(9,1000)xtt,xgg
        goto 999
       endif
       if(jflag.eq.1) jflag=2
      goto 1
 999  continue

      read(4,*)
      write(5,1000)xtt,xgg
      l=0
  5    l=l+1
       read(4,1001,end=888)line
       write(5,1001)line
      goto 5
 888  continue

 1000 format(6x,F6.0,10x,f4.2)
 1001 format(A80)
 1002 format(10x,i3)

      end
