C Output from Public domain Ratfor, version 1.0
      subroutine dqrslm (x, ldx, n, k, qraux, a, lda, job, info, work)
      integer ldx, n, k, lda, job, info
      double precision x(ldx,*), qraux(*), a(lda,*), work(*)
      double precision tmp, alph, ddot
      integer i, j, step
      info = 0
      if( lda .lt. n .or. n .lt. k .or. k .lt. 1 )then
      info = -1
      return
      endif
      if( job .ne. 0 .and. job .ne. 1 )then
      info = 1
      return
      endif
      if( job .eq. 0 )then
      j = 1
      step = 1
      else
      j = k
      step = -1
      endif
23006 if( j .ge. 1 .and. j .le. k )then
      if( qraux(j) .eq. 0.0d0 )then
      j = j + step
      goto 23006
      endif
      tmp = x(j,j)
      x(j,j) = qraux(j)
      i=1
23010 if(.not.(i.lt.j))goto 23012
      alph = - ddot (n-j+1, x(j,j), 1, a(j,i), 1) / x(j,j)
      call daxpy (n-j+1, alph, x(j,j), 1, a(j,i), 1)
23011 i=i+1
      goto 23010
23012 continue
      alph = 1.d0 / x(j,j)
      call dsymv ('l', n-j+1, alph, a(j,j), lda, x(j,j), 1, 0.d0, work(j
     *), 1)
      alph = - ddot (n-j+1, work(j), 1, x(j,j), 1) / 2.d0 / x(j,j)
      call daxpy (n-j+1, alph, x(j,j), 1, work(j), 1)
      call dsyr2 ('l', n-j+1, -1.d0, x(j,j), 1, work(j), 1, a(j,j), lda)
      x(j,j) = tmp
      j = j + step
      goto 23006
      endif
23007 continue
      return
      end
