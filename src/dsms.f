C Output from Public domain Ratfor, version 1.0
      subroutine dsms (s, lds, nobs, nnull, jpvt, q, ldq, nlaht, sms, ld
     *sms, wk, info)
      integer lds, nobs, nnull, jpvt(*), ldq, ldsms, info
      double precision s(lds,*), q(ldq,*), nlaht, sms(ldsms,*), wk(2,*)
      double precision dum, ddot
      integer i, j, n, n0
      info = 0
      if( nnull .lt. 1 .or. nnull .ge. nobs .or. nobs .gt. lds .or. nobs
     * .gt. ldq .or. ldsms .lt. nnull )then
      info = -1
      return
      endif
      n0 = nnull
      n = nobs - nnull
      call dcopy (n-2, q(n0+2,n0+1), ldq+1, wk, 1)
      j=1
23002 if(.not.(j.le.n0))goto 23004
      call dcopy (n, q(n0+1,j), 1, q(j,n0+1), ldq)
      call dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, wk, q(n0+2,j), dum, q(n0+
     *2,j), dum, dum, dum, 01000, info)
23003 j=j+1
      goto 23002
23004 continue
      call dset (n, 10.d0 ** nlaht, wk(2,1), 2)
      call daxpy (n, 1.d0, q(n0+1,n0+1), ldq+1, wk(2,1), 2)
      call dcopy (n-1, q(n0+1,n0+2), ldq+1, wk(1,2), 2)
      call dpbfa (wk, 2, n, 1, info)
      if( info .ne. 0 )then
      info = -2
      return
      endif
      j=1
23007 if(.not.(j.le.n0))goto 23009
      call dpbsl (wk, 2, n, 1, q(n0+1,j))
23008 j=j+1
      goto 23007
23009 continue
      call dcopy (n-2, q(n0+2,n0+1), ldq+1, wk, 1)
      j=1
23010 if(.not.(j.le.n0))goto 23012
      call dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, wk, q(n0+2,j), q(n0+2,j),
     * dum, dum, dum, dum, 10000, info)
23011 j=j+1
      goto 23010
23012 continue
      i=1
23013 if(.not.(i.le.n0))goto 23015
      j=1
23016 if(.not.(j.lt.i))goto 23018
      sms(i,j) = sms(j,i)
23017 j=j+1
      goto 23016
23018 continue
      j=i
23019 if(.not.(j.le.n0))goto 23021
      sms(i,j) = q(j,i) - ddot (n, q(n0+1,j), 1, q(i,n0+1), ldq)
23020 j=j+1
      goto 23019
23021 continue
      sms(i,i) = sms(i,i) + 10.d0**nlaht
23014 i=i+1
      goto 23013
23015 continue
      j=1
23022 if(.not.(j.le.n0))goto 23024
      call dtrsl (s, lds, n0, sms(1,j), 01, info)
23023 j=j+1
      goto 23022
23024 continue
      i=1
23025 if(.not.(i.le.n0))goto 23027
      call dcopy (n0, sms(i,1), ldsms, wk, 1)
      call dtrsl (s, lds, n0, wk, 01, info)
      call dprmut (wk, n0, jpvt, 1)
      call dcopy (n0, wk, 1, sms(i,1), ldsms)
23026 i=i+1
      goto 23025
23027 continue
      j=1
23028 if(.not.(j.le.n0))goto 23030
      call dprmut (sms(1,j), n0, jpvt, 1)
23029 j=j+1
      goto 23028
23030 continue
      j=1
23031 if(.not.(j.le.n0))goto 23033
      call dcopy (n, q(j,n0+1), ldq, q(n0+1,j), 1)
23032 j=j+1
      goto 23031
23033 continue
      return
      end
