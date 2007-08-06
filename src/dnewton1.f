C Output from Public domain Ratfor, version 1.0
      subroutine dnewton1 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qdrs
     *, nqd, qdwt, prec, maxiter, mchpr, mrs, wt, fit, mu, v, jpvt, cdne
     *w, wtnew, fitnew, wk, info)
      integer nxis, nxi, nobs, cntsum, cnt(*), nqd, maxiter, jpvt(*), in
     *fo
      double precision cd(*), q(nxi,*), rs(nxis,*), qdrs(nqd,*), qdwt(*)
     *, prec, mchpr, mrs(*), wt(*), fit(*), mu(*), v(nxis,*), cdnew(*), 
     *wtnew(*), fitnew(*), wk(*)
      integer i, j, k, iter, flag, rkv, idamax, infowk
      double precision wtsum, tmp, ddot, fitmean, lkhd, mumax, wtsumnew,
     * lkhdnew, disc, disc0, trc
      info = 0
      i=1
23000 if(.not.(i.le.nxis))goto 23002
      mrs(i) = 0.d0
      if(cntsum.eq.0)then
      j=1
23005 if(.not.(j.le.nobs))goto 23007
      mrs(i) = mrs(i) + rs(i,j)
23006 j=j+1
      goto 23005
23007 continue
      mrs(i) = mrs(i) / dfloat (nobs)
      else
      j=1
23008 if(.not.(j.le.nobs))goto 23010
      mrs(i) = mrs(i) + rs(i,j) * dfloat (cnt(j))
23009 j=j+1
      goto 23008
23010 continue
      mrs(i) = mrs(i) / dfloat (cntsum)
      endif
23001 i=i+1
      goto 23000
23002 continue
      wtsum = 0.d0
      i=1
23011 if(.not.(i.le.nqd))goto 23013
      wt(i) = qdwt(i) * dexp (ddot (nxis, qdrs(i,1), nqd, cd, 1))
      wtsum = wtsum + wt(i)
23012 i=i+1
      goto 23011
23013 continue
      fitmean = 0.d0
      i=1
23014 if(.not.(i.le.nobs))goto 23016
      tmp = ddot (nxis, rs(1,i), 1, cd, 1) - dlog (wtsum)
      fit(i) = dexp (tmp)
      if(cntsum.ne.0)then
      tmp = tmp * dfloat (cnt(i))
      endif
      fitmean = fitmean + tmp
23015 i=i+1
      goto 23014
23016 continue
      if(cntsum.eq.0)then
      fitmean = fitmean / dfloat (nobs)
      else
      fitmean = fitmean / dfloat (cntsum)
      endif
      call dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, wk, 1)
      lkhd = ddot (nxi, cd, 1, wk, 1) / 2.d0 - fitmean
      iter = 0
      flag = 0
23021 continue
      iter = iter + 1
      i=1
23024 if(.not.(i.le.nxis))goto 23026
      mu(i) = - ddot (nqd, wt, 1, qdrs(1,i), 1) / wtsum
23025 i=i+1
      goto 23024
23026 continue
      i=1
23027 if(.not.(i.le.nxis))goto 23029
      j=i
23030 if(.not.(j.le.nxis))goto 23032
      v(i,j) = 0.d0
      k=1
23033 if(.not.(k.le.nqd))goto 23035
      v(i,j) = v(i,j) + wt(k) * qdrs(k,i) * qdrs(k,j)
23034 k=k+1
      goto 23033
23035 continue
      v(i,j) = v(i,j) / wtsum - mu(i) * mu(j)
      if(j.le.nxi)then
      v(i,j) = v(i,j) + q(i,j)
      endif
23031 j=j+1
      goto 23030
23032 continue
23028 i=i+1
      goto 23027
23029 continue
      call daxpy (nxis, 1.d0, mrs, 1, mu, 1)
      call dsymv ('u', nxi, -1.d0, q, nxi, cd, 1, 1.d0, mu, 1)
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23038 if(.not.(i.le.nxis))goto 23040
      jpvt(i) = 0
23039 i=i+1
      goto 23038
23040 continue
      call dchdc (v, nxis, nxis, wk, jpvt, 1, rkv)
23041 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23041
      endif
23042 continue
      i=rkv+1
23043 if(.not.(i.le.nxis))goto 23045
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23044 i=i+1
      goto 23043
23045 continue
23046 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dset (nxis-rkv, 0.d0, cdnew(rkv+1), 1)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      wtsumnew = 0.d0
      i=1
23049 if(.not.(i.le.nqd))goto 23051
      wtnew(i) = qdwt(i) * dexp (ddot (nxis, qdrs(i,1), nqd, cdnew, 1))
      wtsumnew = wtsumnew + wtnew(i)
23050 i=i+1
      goto 23049
23051 continue
      fitmean = 0.d0
      i=1
23052 if(.not.(i.le.nobs))goto 23054
      tmp = ddot (nxis, rs(1,i), 1, cdnew, 1) - dlog (wtsumnew)
      if(tmp.gt.3.d2)then
      flag = flag + 1
      goto 23054
      endif
      fitnew(i) = dexp (tmp)
      if(cntsum.ne.0)then
      tmp = tmp * dfloat (cnt(i))
      endif
      fitmean = fitmean + tmp
23053 i=i+1
      goto 23052
23054 continue
      if(cntsum.eq.0)then
      fitmean = fitmean / dfloat (nobs)
      else
      fitmean = fitmean / dfloat (cntsum)
      endif
      call dsymv ('u', nxi, 1.d0, q, nxi, cdnew, 1, 0.d0, wk, 1)
      lkhdnew = ddot (nxi, cdnew, 1, wk, 1) / 2.d0 - fitmean
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      wtsum = 0.d0
      i=1
23063 if(.not.(i.le.nqd))goto 23065
      wt(i) = qdwt(i)
      wtsum = wtsum + wt(i)
23064 i=i+1
      goto 23063
23065 continue
      call dset (nobs, 1.d0/wtsum, fit, 1)
      fitmean = - dlog (wtsum)
      lkhd = - fitmean
      iter = 0
      goto 23048
      endif
      if(flag.eq.3)then
      goto 23048
      endif
      if(lkhdnew-lkhd.lt.1.d1*(1.d0+dabs(lkhd))*mchpr)then
      goto 23048
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23048
      endif
23047 goto 23046
23048 continue
      if(flag.eq.1)then
      flag = 2
      goto 23022
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      i=1
23076 if(.not.(i.le.nqd))goto 23078
      disc = dmax1 (disc, dabs(wt(i)-wtnew(i))/(1.d0+dabs(wt(i))))
23077 i=i+1
      goto 23076
23078 continue
      i=1
23079 if(.not.(i.le.nobs))goto 23081
      disc = dmax1 (disc, dabs(fit(i)-fitnew(i))/(1.d0+dabs(fit(i))))
23080 i=i+1
      goto 23079
23081 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
      disc0 = dmax1 ((mumax/(1.d0+lkhd))**2, dabs(lkhd-lkhdnew)/(1.d0+da
     *bs(lkhd)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nqd, wtnew, 1, wt, 1)
      wtsum = wtsumnew
      call dcopy (nobs, fitnew, 1, fit, 1)
      lkhd = lkhdnew
      if(disc0.lt.prec)then
      goto 23023
      endif
      if(disc.lt.prec)then
      goto 23023
      endif
      if(iter.lt.maxiter)then
      goto 23022
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      wtsum = 0.d0
      i=1
23090 if(.not.(i.le.nqd))goto 23092
      wt(i) = qdwt(i)
      wtsum = wtsum + wt(i)
23091 i=i+1
      goto 23090
23092 continue
      call dset (nobs, 1.d0/wtsum, fit, 1)
      fitmean = - dlog (wtsum)
      lkhd = - fitmean
      iter = 0
      flag = 2
      else
      info = 2
      goto 23023
      endif
23022 goto 23021
23023 continue
      i=1
23093 if(.not.(i.le.nobs))goto 23095
      call daxpy (nxis, -1.d0, mrs, 1, rs(1,i), 1)
      call dprmut (rs(1,i), nxis, jpvt, 0)
      if(cntsum.ne.0)then
      call dscal (nxis, dsqrt(dfloat(cnt(i))), rs(1,i), 1)
      endif
      call dtrsl (v, nxis, nxis, rs(1,i), 11, infowk)
23094 i=i+1
      goto 23093
23095 continue
      trc = ddot (nobs*nxis, rs, 1, rs, 1)
      if(cntsum.eq.0)then
      trc = trc / dfloat(nobs) / (dfloat(nobs)-1.d0)
      lkhd = 0.d0
      i=1
23100 if(.not.(i.le.nobs))goto 23102
      lkhd = lkhd + dlog (fit(i))
23101 i=i+1
      goto 23100
23102 continue
      lkhd = lkhd / dfloat (nobs)
      else
      trc = trc / dfloat(cntsum) / (dfloat(cntsum)-1.d0)
      lkhd = 0.d0
      i=1
23103 if(.not.(i.le.nobs))goto 23105
      lkhd = lkhd + dfloat (cnt(i)) * dlog (fit(i))
23104 i=i+1
      goto 23103
23105 continue
      lkhd = lkhd / dfloat (cntsum)
      endif
      mrs(1) = lkhd
      mrs(2) = trc
      mrs(3) = wtsum
      return
      end
