C Output from Public domain Ratfor, version 1.01
      subroutine reg (sr, nobs, nnull, q, nxi, y, method, alpha, varht, 
     *score, dc, mchpr, v, mu, jpvt, wk, rkv, info)
      double precision sr(nobs,*), q(nxi,*), y(*), alpha, varht, score, 
     *dc(*), mchpr, v(nnull+nxi,*), mu(*), wk(*)
      integer nobs, nnull, nxi, method, jpvt(*), rkv, info
      double precision ddot, dasum, rss, trc, dum
      integer i, j, nn, idamax, infowk, idum
      info = 0
      nn = nnull + nxi
      i=1
23000 if(.not.(i.le.nn))goto 23002
      mu(i) = ddot (nobs, sr(1,i), 1, y, 1)
      j=i
23003 if(.not.(j.le.nn))goto 23005
      v(i,j) = ddot (nobs, sr(1,i), 1, sr(1,j), 1)
      if(i.gt.nnull)then
      v(i,j) = v(i,j) + q(i-nnull,j-nnull)
      endif
23004 j=j+1
      goto 23003
23005 continue
23001 i=i+1
      goto 23000
23002 continue
      infowk = 0
      i=1
23008 if(.not.(i.le.nn))goto 23010
      infowk = infowk + jpvt(i)
23009 i=i+1
      goto 23008
23010 continue
      call dchdc (v, nn, nn, wk, jpvt, 1, rkv)
      j = idamax (rkv-infowk, v(infowk+1,infowk+1), nn+1)
23011 if(v(rkv,rkv).lt.v(infowk+j,infowk+j)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23011
      endif
23012 continue
      i=rkv+1
23013 if(.not.(i.le.nn))goto 23015
      v(i,i) = v(j,j)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23014 i=i+1
      goto 23013
23015 continue
      call dcopy (nn, mu, 1, dc, 1)
      call dprmut (dc, nn, jpvt, 0)
      call dtrsl (v, nn, nn, dc, 11, infowk)
      call dset (nn-rkv, 0.d0, dc(rkv+1), 1)
      call dtrsl (v, nn, nn, dc, 01, infowk)
      call dprmut (dc, nn, jpvt, 1)
      if(method.eq.4)then
      return
      endif
      i=1
23018 if(.not.(i.le.nobs))goto 23020
      wk(i) = y(i) - ddot (nn, sr(i,1), nobs, dc, 1)
23019 i=i+1
      goto 23018
23020 continue
      if(method.eq.5)then
      wk(nobs+1) = ddot (nobs, wk, 1, wk, 1) / dfloat (nobs)
      i=1
23023 if(.not.(i.le.nobs))goto 23025
      call dcopy (nn, sr(i,1), nobs, mu, 1)
      call dprmut (mu, nn, jpvt, 0)
      call dtrsl (v, nn, nn, mu, 11, infowk)
      wk(i) = ddot (nn, mu, 1, mu, 1)
23024 i=i+1
      goto 23023
23025 continue
      return
      endif
      if(method.eq.3)then
      rss = ddot (nobs, y, 1, wk, 1)
      if(nnull.gt.0)then
      call dqrdc (sr, nobs, nobs, nnull, wk, idum, dum, 0)
      i=1
23030 if(.not.(i.le.nxi))goto 23032
      call dqrsl (sr, nobs, nobs, nnull, wk, sr(1,nnull+i), dum, sr(1,nn
     *ull+i), dum, dum, dum, 01000, infowk)
23031 i=i+1
      goto 23030
23032 continue
      endif
      call dcopy (nxi, q, nxi+1, wk, 1)
      i=1
23033 if(.not.(i.le.nxi))goto 23035
      j=i
23036 if(.not.(j.le.nxi))goto 23038
      q(i,j) = q(i,j) + ddot (nobs-nnull, sr(nnull+1,nnull+i), 1, sr(nnu
     *ll+1,nnull+j), 1)
23037 j=j+1
      goto 23036
23038 continue
23034 i=i+1
      goto 23033
23035 continue
      i=1
23039 if(.not.(i.le.nxi))goto 23041
      j=i
23042 if(.not.(j.le.nxi))goto 23044
      sr(i,j) = q(i,j)
      sr(j,i) = q(i,j)
      q(i,j) = q(j,i)
23043 j=j+1
      goto 23042
23044 continue
23040 i=i+1
      goto 23039
23041 continue
      call dcopy (nxi, wk, 1, q, nxi+1)
      call dsyev ('n', 'u', nxi, sr, nobs, mu, wk, 3*nxi, info)
      trc = 0.d0
      i=1
23045 if(.not.(i.le.rkv-nnull))goto 23047
      trc = trc + dlog (mu(nxi-i+1))
23046 i=i+1
      goto 23045
23047 continue
      call dsyev ('n', 'u', nxi, q, nxi, mu, wk, 3*nxi, info)
      i=1
23048 if(.not.(i.le.rkv-nnull))goto 23050
      trc = trc - dlog (mu(nxi-i+1))
23049 i=i+1
      goto 23048
23050 continue
      score = rss / dfloat (nobs) * dexp (trc/dfloat(nobs-nnull))
      varht = rss / dfloat (nobs-nnull)
      else
      rss = ddot (nobs, wk, 1, wk, 1) / dfloat (nobs)
      i=1
23051 if(.not.(i.le.nobs))goto 23053
      call dcopy (nn, sr(i,1), nobs, mu, 1)
      call dprmut (mu, nn, jpvt, 0)
      call dtrsl (v, nn, nn, mu, 11, infowk)
      wk(i) = ddot (nn, mu, 1, mu, 1)
23052 i=i+1
      goto 23051
23053 continue
      trc = dasum (nobs, wk, 1) / dfloat (nobs)
      if(method.eq.2)then
      score = rss / (1.d0-alpha*trc)**2
      varht = rss / (1.d0-trc)
      else
      score = rss + 2.d0 * varht * alpha * trc
      endif
      endif
      wk(1) = rss
      wk(2) = trc
      return
      end
      subroutine regaux (v, nn, jpvt, rkv, r, nr, sms, nnull, wk)
      double precision v(nn,*), r(nn,*), sms(nnull,*), wk(nn,*)
      integer nn, jpvt(*), rkv, nr, nnull
      double precision ddot
      integer i, j, infowk
      i=1
23056 if(.not.(i.le.nr))goto 23058
      call dprmut (r(1,i), nn, jpvt, 0)
      call dtrsl (v, nn, nn, r(1,i), 11, infowk)
      if(nn-rkv.gt.0)then
      call dset (nn-rkv, 0.d0, r(rkv+1,i), 1)
      endif
      call dtrsl (v, nn, nn, r(1,i), 01, infowk)
      call dprmut (r(1,i), nn, jpvt, 1)
23057 i=i+1
      goto 23056
23058 continue
      call dset (nn*nnull, 0.d0, wk, 1)
      call dset (nnull, 1.d0, wk, nn+1)
      i=1
23061 if(.not.(i.le.nnull))goto 23063
      call dtrsl (v, nn, nn, wk(1,i), 11, infowk)
23062 i=i+1
      goto 23061
23063 continue
      i=1
23064 if(.not.(i.le.nnull))goto 23066
      j=i
23067 if(.not.(j.le.nnull))goto 23069
      sms(i,j) = ddot (nn, wk(1,i), 1, wk(1,j), 1)
      sms(j,i) = sms(i,j)
23068 j=j+1
      goto 23067
23069 continue
23065 i=i+1
      goto 23064
23066 continue
      return
      end
