      subroutine dcrdr (s, lds, nobs, nnull, qraux, jpvt, q, ldq, nlaht,
     & r, ldr, nr, cr, ldcr, dr, lddr, wk, info)
      integer lds, nobs, nnull, jpvt(*), ldq, ldr, nr, ldcr, lddr, info
      double precision s(lds,*), qraux(*), q(ldq,*), nlaht, r(ldr,*), 
     &cr(ldcr,*), dr(lddr,*), wk(2,*)
      double precision dum, ddot
      integer i, j, n, n0
      info = 0
      if(.not.( nnull .lt. 1 .or. nnull .ge. nobs .or. nobs .gt. lds 
     &.or. nobs .gt. ldq .or. ldr .lt. nobs .or. nr .lt. 1 .or. ldcr 
     &.lt. nobs .or. lddr .lt. nnull ))goto 23000
      info = -1
      return
23000 continue
      n0 = nnull
      n = nobs - nnull
      j=1
23002 if(.not.(j.le.nr))goto 23004
      call dcopy (nobs, r(1,j), 1, cr(1,j), 1)
      j=j+1
      goto 23002
23004 continue
      call dcopy (n-2, q(n0+2,n0+1), ldq+1, wk, 1)
      j=1
23005 if(.not.(j.le.nr))goto 23007
      call dqrsl (s, lds, nobs, nnull, qraux, cr(1,j), dum, cr(1,j), 
     &dum, dum, dum, 01000, info)
      call dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, wk, cr(n0+2,j), dum, cr(
     &n0+2,j), dum, dum, dum, 01000, info)
      j=j+1
      goto 23005
23007 continue
      call dset (n, 10.d0 ** nlaht, wk(2,1), 2)
      call daxpy (n, 1.d0, q(n0+1,n0+1), ldq+1, wk(2,1), 2)
      call dcopy (n-1, q(n0+1,n0+2), ldq+1, wk(1,2), 2)
      call dpbfa (wk, 2, n, 1, info)
      if(.not.( info .ne. 0 ))goto 23008
      info = -2
      return
23008 continue
      j=1
23010 if(.not.(j.le.nr))goto 23012
      call dpbsl (wk, 2, n, 1, cr(n0+1,j))
      j=j+1
      goto 23010
23012 continue
      call dcopy (n-2, q(n0+2,n0+1), ldq+1, wk, 1)
      j=1
23013 if(.not.(j.le.nr))goto 23015
      call dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, wk, cr(n0+2,j), cr(n0+2,
     &j), dum, dum, dum, dum, 10000, info)
      j=j+1
      goto 23013
23015 continue
      j=1
23016 if(.not.(j.le.nr))goto 23018
      i=1
23019 if(.not.(i.le.n0))goto 23021
      dr(i,j) = cr(i,j) - ddot (n, cr(n0+1,j), 1, q(n0+1,i), 1)
      i=i+1
      goto 23019
23021 continue
      call dtrsl (s, lds, n0, dr(1,j), 01, info)
      call dprmut (dr(1,j), n0, jpvt, 1)
      j=j+1
      goto 23016
23018 continue
      j=1
23022 if(.not.(j.le.nr))goto 23024
      call dset (n0, 0.d0, cr(1,j), 1)
      call dqrsl (s, lds, nobs, nnull, qraux, cr(1,j), cr(1,j), dum, 
     &dum, dum, dum, 10000, info)
      j=j+1
      goto 23022
23024 continue
      return
      end
