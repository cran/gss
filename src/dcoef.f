      subroutine dcoef (s, lds, nobs, nnull, qraux, jpvt, z, q, ldq, 
     *nlaht, c, d, info, twk)
      integer lds, nobs, nnull, jpvt(*), ldq, info
      double precision s(lds,*), qraux(*), z(*), q(ldq,*), nlaht, c(*), 
     *d(*), twk(2,*)
      double precision dum, ddot
      integer n, n0
      info = 0
      if(.not.( nnull .lt. 1 .or. nnull .ge. nobs .or. nobs .gt. lds 
     *.or. nobs .gt. ldq ))goto 23000
      info = -1
      return
23000 continue
      n0 = nnull
      n = nobs - nnull
      call dset (n, 10.d0 ** nlaht, twk(2,1), 2)
      call daxpy (n, 1.d0, q(n0+1,n0+1), ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(n0+1,n0+2), ldq+1, twk(1,2), 2)
      call dpbfa (twk, 2, n, 1, info)
      if(.not.( info .ne. 0 ))goto 23002
      info = -2
      return
23002 continue
      call dpbsl (twk, 2, n, 1, z(n0+1))
      call dcopy (n-2, q(n0+2,n0+1), ldq+1, twk, 1)
      call dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, twk, z(n0+2), z(n0+2), 
     *dum, dum, dum, dum, 10000, info)
      call dset (n0, 0.d0, c, 1)
      call dcopy (n, z(n0+1), 1, c(n0+1), 1)
      call dqrsl (s, lds, nobs, nnull, qraux, c, c, dum, dum, dum, dum, 
     *10000, info)
      j=1
23004 if(.not.(j.le.n0))goto 23006
      d(j) = z(j) - ddot (n, z(n0+1), 1, q(n0+1,j), 1)
      j=j+1
      goto 23004
23006 continue
      call dtrsl (s, lds, n0, d, 01, info)
      call dprmut (d, n0, jpvt, 1)
      return
      end
