      subroutine dcore (vmu, q, ldq, nobs, nnull, tol, z, job, limnla, 
     *nlaht, score, varht, info, twk, work)
      character*1 vmu
      integer ldq, nobs, nnull, job, info
      double precision q(ldq,*), tol, z(*), limnla(2), nlaht, score(*), 
     *varht, twk(2,*), work(*)
      double precision dum, low, upp, dasum, mchpr
      integer n0, n, j
      info = 0
      if(.not.( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' ))
     *goto 23000
      info = -3
      return
23000 continue
      if(.not.( nnull .lt. 1 .or. nobs .le. nnull .or. nobs .gt. ldq ))
     *goto 23002
      info = -1
      return
23002 continue
      n0 = nnull
      n = nobs - nnull
      call dsytr (q(n0+1,n0+1), ldq, n, tol, info, work)
      if(.not.( info .ne. 0 ))goto 23004
      return
23004 continue
      call dcopy (n-2, q(n0+2,n0+1), ldq+1, work, 1)
      call dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, work, z(n0+2), dum, z(n0+
     *2), dum, dum, dum, 01000, info)
      if(.not.( job .eq. 0 ))goto 23006
      mchpr = 1.d0
23008 if(.not.( 1.d0 + mchpr .gt. 1.d0 ))goto 23009
      mchpr = mchpr / 2.d0
      goto 23008
23009 continue
      mchpr = mchpr * 2.d0
      limnla(2) = dmax1 (dasum (n, q(n0+1,n0+1), ldq+1) * 1.d2, mchpr)
      limnla(1) = limnla(2) * mchpr
      limnla(2) = dlog10 (limnla(2))
      limnla(1) = dlog10 (limnla(1))
23006 continue
      low = limnla(1)
      upp = limnla(2)
      if(.not.( job .le. 0 ))goto 23010
      call dgold (vmu, q(n0+1,n0+1), ldq, n, z(n0+1), low, upp, nlaht, 
     *score(1), varht, info, twk, work)
      if(.not.( vmu .eq. 'v' ))goto 23012
      score(1) = score(1) * dfloat (nobs) / dfloat (n)
23012 continue
      if(.not.( vmu .eq. 'm' ))goto 23014
      score(1) = score(1) * dfloat (n) / dfloat (nobs)
23014 continue
      if(.not.( vmu .eq. 'u' ))goto 23016
      score(1) = score(1) * dfloat (n) / dfloat (nobs) + 2.d0 * varht
23016 continue
      goto 23011
23010 continue
      call deval (vmu, q(n0+1,n0+1), ldq, n, z(n0+1), job, low, upp, 
     *nlaht, score, varht, info, twk, work)
      dum = dfloat (nobs) / dfloat (n)
      j=1
23018 if(.not.(j.le.job+1))goto 23020
      if(.not.( vmu .eq. 'v' ))goto 23021
      score(j) = score(j) * dum
23021 continue
      if(.not.( vmu .eq. 'm' ))goto 23023
      score(j) = score(j) / dum
23023 continue
      if(.not.( vmu .eq. 'u' ))goto 23025
      score(j) = score(j) / dum + 2.d0 * varht
23025 continue
      j=j+1
      goto 23018
23020 continue
23011 continue
      return
      end
