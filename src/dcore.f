C Output from Public domain Ratfor, version 1.0
      subroutine dcore (vmu, q, ldq, nobs, nnull, tol, z, job, limnla, n
     *laht, score, varht, info, twk, work)
      character vmu
      integer ldq, nobs, nnull, job, info
      double precision q(ldq,*), tol, z(*), limnla(2), nlaht, score(*), 
     *varht, twk(2,*), work(*)
      double precision dum, low, upp, dasum, mchpr
      integer n0, n, j
      info = 0
      if( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' )then
      info = -3
      return
      endif
      if( nnull .lt. 1 .or. nobs .le. nnull .or. nobs .gt. ldq )then
      info = -1
      return
      endif
      n0 = nnull
      n = nobs - nnull
      call dsytr (q(n0+1,n0+1), ldq, n, tol, info, work)
      if( info .ne. 0 )then
      return
      endif
      call dcopy (n-2, q(n0+2,n0+1), ldq+1, work, 1)
      call dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, work, z(n0+2), dum, z(n0+
     *2), dum, dum, dum, 01000, info)
      if( job .eq. 0 )then
      mchpr = 1.d0
23008 if( 1.d0 + mchpr .gt. 1.d0 )then
      mchpr = mchpr / 2.d0
      goto 23008
      endif
23009 continue
      mchpr = mchpr * 2.d0
      limnla(2) = dmax1 (dasum (n, q(n0+1,n0+1), ldq+1) * 1.d2, mchpr)
      limnla(1) = limnla(2) * mchpr
      limnla(2) = dlog10 (limnla(2))
      limnla(1) = dlog10 (limnla(1))
      endif
      low = limnla(1)
      upp = limnla(2)
      if( job .le. 0 )then
      call dgold (vmu, q(n0+1,n0+1), ldq, n, z(n0+1), low, upp, nlaht, s
     *core(1), varht, info, twk, work)
      if( vmu .eq. 'v' )then
      score(1) = score(1) * dfloat (nobs) / dfloat (n)
      endif
      if( vmu .eq. 'm' )then
      score(1) = score(1) * dfloat (n) / dfloat (nobs)
      endif
      if( vmu .eq. 'u' )then
      score(1) = score(1) * dfloat (n) / dfloat (nobs) + 2.d0 * varht
      endif
      else
      call deval (vmu, q(n0+1,n0+1), ldq, n, z(n0+1), job, low, upp, nla
     *ht, score, varht, info, twk, work)
      dum = dfloat (nobs) / dfloat (n)
      j=1
23018 if(.not.(j.le.job+1))goto 23020
      if( vmu .eq. 'v' )then
      score(j) = score(j) * dum
      endif
      if( vmu .eq. 'm' )then
      score(j) = score(j) / dum
      endif
      if( vmu .eq. 'u' )then
      score(j) = score(j) / dum + 2.d0 * varht
      endif
23019 j=j+1
      goto 23018
23020 continue
      endif
      return
      end
