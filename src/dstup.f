      subroutine dstup (s, lds, nobs, nnull, qraux, jpvt, y, q, ldqr, 
     *ldqc, nq, info, work)
      integer lds, nobs, nnull, jpvt(*), ldqr, ldqc, nq, info
      double precision s(lds,*), y(*), qraux(*), q(ldqr,ldqc,*), work(*)
      double precision dum
      integer j
      info = 0
      if(.not.( nobs .lt. 1 .or. nobs .gt. lds .or. nobs .gt. ldqr .or. 
     *nobs .gt. ldqc ))goto 23000
      info = -1
      return
23000 continue
      j=1
23002 if(.not.(j.le.nnull))goto 23004
      jpvt(j) = 0
      j=j+1
      goto 23002
23004 continue
      call dqrdc (s, lds, nobs, nnull, qraux, jpvt, work, 1)
      call dqrsl (s, lds, nobs, nnull, qraux, y, dum, y, work, dum, dum,
     * 01100, info)
      if(.not.( info .ne. 0 ))goto 23005
      return
23005 continue
      j=1
23007 if(.not.(j.le.nq))goto 23009
      call dqrslm (s, lds, nobs, nnull, qraux, q(1,1,j), ldqr, 0, info, 
     *work)
      j=j+1
      goto 23007
23009 continue
      return
      end
