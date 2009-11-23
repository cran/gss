C Output from Public domain Ratfor, version 1.0
      subroutine dmcdc (a, lda, p, e, jpvt, info)
      integer lda, p, jpvt(*), info
      double precision a(lda,*), e(*)
      double precision beta, delta, theta, tmp, dasum, ddot
      integer i, j, jmax, jtmp, idamax
      info = 0
      if( lda .lt. p .or. p .lt. 1 )then
      info = -1
      return
      endif
      tmp = 1.d0
23002 if( 1.d0 + tmp .gt. 1.d0 )then
      tmp = tmp / 2.d0
      goto 23002
      endif
23003 continue
      jmax = idamax (p, a, lda+1)
      beta = dmax1 (2.d0 * tmp, dabs (a(jmax,jmax)))
      tmp = dsqrt (dfloat (p*p-1))
      if( tmp .lt. 1.d0 )then
      tmp = 1.d0
      endif
      j=2
23006 if(.not.(j.le.p))goto 23008
      jmax = idamax (j-1, a(1,j), 1)
      beta = dmax1 (beta, dabs (a(jmax,j)) / tmp)
23007 j=j+1
      goto 23006
23008 continue
      delta = dasum (p, a, lda+1) / dfloat (p) * 1.d-7
      delta = dmax1 (delta, 1.d-10)
      j=1
23009 if(.not.(j.le.p))goto 23011
      jpvt(j) = j
23010 j=j+1
      goto 23009
23011 continue
      j=1
23012 if(.not.(j.le.p))goto 23014
      jmax = idamax (p-j+1, a(j,j), lda+1) + j - 1
      if( jmax .ne. j )then
      call dswap (j-1, a(1,j), 1, a(1,jmax), 1)
      call dswap (jmax-j-1, a(j,j+1), lda, a(j+1,jmax), 1)
      call dswap (p-jmax, a(j,jmax+1), lda, a(jmax,jmax+1), lda)
      tmp = a(j,j)
      a(j,j) = a(jmax,jmax)
      a(jmax,jmax) = tmp
      jtmp = jpvt(j)
      jpvt(j) = jpvt(jmax)
      jpvt(jmax) = jtmp
      endif
      i=1
23017 if(.not.(i.lt.j))goto 23019
      a(i,j) = a(i,j) / a(i,i)
23018 i=i+1
      goto 23017
23019 continue
      i=j+1
23020 if(.not.(i.le.p))goto 23022
      a(j,i) = a(j,i) - ddot (j-1, a(1,j), 1, a(1,i), 1)
23021 i=i+1
      goto 23020
23022 continue
      if( j .eq. p )then
      theta = 0.d0
      else
      jmax = idamax (p-j, a(j,j+1), lda) + j
      theta = dabs (a(j,jmax))
      endif
      tmp = dmax1 (delta, dabs (a(j,j)), theta ** 2 / beta)
      e(j) = tmp - a(j,j)
      a(j,j) = tmp
      i=j+1
23025 if(.not.(i.le.p))goto 23027
      a(i,i) = a(i,i) - a(j,i) ** 2 / a(j,j)
23026 i=i+1
      goto 23025
23027 continue
23013 j=j+1
      goto 23012
23014 continue
      j=1
23028 if(.not.(j.le.p))goto 23030
      a(j,j) = dsqrt (a(j,j))
      call dscal (p-j, a(j,j), a(j,j+1), lda)
23029 j=j+1
      goto 23028
23030 continue
      return
      end
