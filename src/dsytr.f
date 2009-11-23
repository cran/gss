C Output from Public domain Ratfor, version 1.0
      subroutine dsytr (x, ldx, n, tol, info, work)
      integer ldx, n, info
      double precision x(ldx,*), tol, work(*)
      double precision nrmtot, nrmxj, alph, toltot, tolcum, toluni, dn, 
     *ddot
      integer j
      info = 0
      if( ldx .lt. n .or. n .le. 2 )then
      info = -1
      return
      endif
      nrmtot = ddot (n, x, ldx+1, x, ldx+1)
      j=1 
23002 if(.not.(j.lt.n ))goto 23004
      nrmtot = nrmtot + 2.d0 * ddot (n-j, x(j+1,j), 1, x(j+1,j), 1)
23003 j=j+1 
      goto 23002
23004 continue
      toltot = 1.d0
23005 if( 1.d0 + toltot .gt. 1.d0 )then
      toltot = toltot / 2.d0
      goto 23005
      endif
23006 continue
      toltot = 4.d0 * toltot ** 2
      if( toltot .lt. tol )then
      toltot = tol
      endif
      toltot = toltot * nrmtot
      dn = dfloat (n)
      toluni = toltot * 6.d0 / dn / ( dn - 1.d0 ) / ( 2.d0 * dn - 1.d0 )
      tolcum = 0.d0
      j=1 
23009 if(.not.(j.lt.n-1 ))goto 23011
      nrmtot = nrmtot - x(j,j) * x(j,j)
      nrmxj = ddot (n-j, x(j+1,j), 1, x(j+1,j), 1)
      dn = dfloat (n-j)
      tolcum = tolcum + toluni * dn * dn
      if( 2.d0 * nrmxj .le. tolcum )then
      x(j,j+1) = 0.d0
      call dscal (n-j, 0.d0, x(j+1,j), 1)
      tolcum = tolcum - 2.d0 * nrmxj
      toltot = toltot - 2.d0 * nrmxj
      goto 23010
      endif
      if( x(j+1,j) .lt. 0.d0 )then
      x(j,j+1) = dsqrt (nrmxj)
      else
      x(j,j+1) = - dsqrt (nrmxj)
      endif
      nrmtot = nrmtot - 2.d0 * nrmxj
      call dscal (n-j, -1.d0/x(j,j+1), x(j+1,j), 1)
      x(j+1,j) = 1.d0 + x(j+1,j)
      alph = 1.d0 / x(j+1,j)
      call dsymv ('l', n-j, alph, x(j+1,j+1), ldx, x(j+1,j), 1, 0.d0, wo
     *rk(j+1), 1)
      alph = - ddot (n-j, work(j+1), 1, x(j+1,j), 1) / 2.d0 / x(j+1,j)
      call daxpy (n-j, alph, x(j+1,j), 1, work(j+1), 1)
      call dsyr2 ('l', n-j, -1.d0, x(j+1,j), 1, work(j+1), 1, x(j+1,j+1)
     *, ldx)
23010 j=j+1 
      goto 23009
23011 continue
      x(n-1,n) = x(n,n-1)
      return
      end
