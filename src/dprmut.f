      subroutine dprmut (x,npar,jpvt,job)
      integer npar,jpvt(npar),job
      double precision x(npar)
c
c Purpose: permute the elements of the array x according to the index 
c	vector jpvt (either forward or backward permutation).
c
c On Entry:
c   x(npar)		array to be permuted
c   npar		size of x (and jpvt)
c   jpvt		indices of the permutation
c   job			indicator of forward or backward permutation
c			if job = 0 forward permutation  
c				x(jpvt(i)) moved to x(i)
c			if job is nonzero backward permutation 
c				x(i) moved to x(jpvt(i))
c On Exit:
c   x(npar)		array with permuted entries
c
c   Written:	Yin Ling	U. of Maryland, August,1978
c
c $Header: dprmut.f,v 2.1 86/04/08 14:05:53 lindstrom Exp $
c
      integer i,j,k
      double precision t
c
      if (npar .le. 1) then
         return
      endif
      do 10 j = 1,npar
         jpvt(j) = -jpvt(j)
   10 continue
      if (job .eq. 0) then
c		forward permutation
         do 30 i = 1,npar 
            if (jpvt(i) .gt. 0) then
               goto 30
            endif
            j = i
            jpvt(j) = -jpvt(j)
            k = jpvt(j)
c           while
   20       if (jpvt(k) .lt. 0) then
               t = x(j)
               x(j) = x(k)
               x(k) = t
               jpvt(k) = -jpvt(k)
               j = k
               k = jpvt(k)
               goto 20
c           endwhile
            endif
   30    continue
      endif
      if (job .ne. 0 ) then
c			backward permutation
         do 50 i = 1,npar 
            if (jpvt(i) .gt. 0) then
               goto 50
            endif
            jpvt(i) = -jpvt(i)
            j = jpvt(i)
c           while
   40       if (j .ne. i) then
               t = x(i)
               x(i) = x(j)
               x(j) = t
               jpvt(j) = -jpvt(j)
               j = jpvt(j)
               goto 40
c           endwhile
            endif
   50    continue
      endif
      return
      end
