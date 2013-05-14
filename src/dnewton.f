C Output from Public domain Ratfor, version 1.01
      subroutine dnewton (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qdrs,
     * nqd, nt, bwt, qdwt, prec, maxiter, mchpr, jpvt, wk, info)
      integer nxis, nxi, nobs, cntsum, cnt(*), nqd, nt, maxiter, jpvt(*)
     *, info
      double precision cd(*), q(nxi,*), rs(nxis,*), qdrs(nqd,*), bwt(*),
     * qdwt(nt,*), prec, mchpr, wk(*)
      integer imrs, iwt, iwtsum, ifit, imu, imuwk, iv, ivwk, icdnew, iwt
     *new, iwtsumnew, ifitnew, iwk
      imrs = 1
      iwt = imrs + max0 (nxis, 3)
      iwtsum = iwt + nqd*nt
      ifit = iwtsum + nt
      imu = ifit + nobs
      imuwk = imu + nxis
      iv = imuwk + nxis
      ivwk = iv + nxis*nxis
      icdnew = ivwk + nxis*nxis
      iwtnew = icdnew + nxis
      iwtsumnew = iwtnew + nqd*nt
      ifitnew = iwtsumnew + nt
      iwk = ifitnew + nobs
      call dnewton1 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qdrs, nqd,
     * nt, bwt, qdwt, prec, maxiter, mchpr, wk(imrs), wk(iwt), wk(iwtsum
     *), wk(ifit), wk(imu), wk(imuwk), wk(iv), wk(ivwk), jpvt, wk(icdnew
     *), wk(iwtnew), wk(iwtsumnew), wk(ifitnew), wk(iwk), info)
      return
      end
      subroutine dnewton1 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qdrs
     *, nqd, nt, bwt, qdwt, prec, maxiter, mchpr, mrs, wt, wtsum, fit, m
     *u, muwk, v, vwk, jpvt, cdnew, wtnew, wtsumnew, fitnew, wk, info)
      integer nxis, nxi, nobs, cntsum, cnt(*), nqd, nt, maxiter, jpvt(*)
     *, info
      double precision cd(*), q(nxi,*), rs(nxis,*), qdrs(nqd,*), bwt(*),
     * qdwt(nt,*), prec, mchpr, mrs(*), wt(nt,*), wtsum(*), fit(*), mu(*
     *), muwk(*), v(nxis,*), vwk(nxis,*), cdnew(*), wtnew(nt,*), wtsumne
     *w(*), fitnew(*), wk(*)
      integer i, j, k, m, iter, flag, rkv, idamax, infowk
      double precision norm, tmp, ddot, fitmean, lkhd, mumax, lkhdnew, d
     *isc, disc0, trc
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
      m=1
23011 if(.not.(m.le.nt))goto 23013
      wtsum(m) = 0.d0
23012 m=m+1
      goto 23011
23013 continue
      i=1
23014 if(.not.(i.le.nqd))goto 23016
      tmp = dexp (ddot (nxis, qdrs(i,1), nqd, cd, 1))
      m=1
23017 if(.not.(m.le.nt))goto 23019
      wt(m,i) = qdwt(m,i) * tmp
      wtsum(m) = wtsum(m) + wt(m,i)
23018 m=m+1
      goto 23017
23019 continue
23015 i=i+1
      goto 23014
23016 continue
      norm = 0.d0
      m=1
23020 if(.not.(m.le.nt))goto 23022
      norm = norm + bwt(m) * dlog (wtsum(m))
23021 m=m+1
      goto 23020
23022 continue
      fitmean = 0.d0
      i=1
23023 if(.not.(i.le.nobs))goto 23025
      tmp = ddot (nxis, rs(1,i), 1, cd, 1)
      fit(i) = dexp (tmp)
      if(cntsum.ne.0)then
      tmp = tmp * dfloat (cnt(i))
      endif
      fitmean = fitmean + tmp
23024 i=i+1
      goto 23023
23025 continue
      if(cntsum.eq.0)then
      fitmean = fitmean / dfloat (nobs)
      else
      fitmean = fitmean / dfloat (cntsum)
      endif
      call dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, wk, 1)
      lkhd = ddot (nxi, cd, 1, wk, 1) / 2.d0 - fitmean + norm
      iter = 0
      flag = 0
23030 continue
      iter = iter + 1
      call dset(nxis, 0.d0, mu, 1)
      call dset(nxis*nxis, 0.d0, v, 1)
      m=1
23033 if(.not.(m.le.nt))goto 23035
      i=1
23036 if(.not.(i.le.nxis))goto 23038
      muwk(i) = - ddot (nqd, wt(m,1), nt, qdrs(1,i), 1) / wtsum(m)
23037 i=i+1
      goto 23036
23038 continue
      i=1
23039 if(.not.(i.le.nxis))goto 23041
      j=i
23042 if(.not.(j.le.nxis))goto 23044
      vwk(i,j) = 0.d0
      k=1
23045 if(.not.(k.le.nqd))goto 23047
      vwk(i,j) = vwk(i,j) + wt(m,k) * qdrs(k,i) * qdrs(k,j)
23046 k=k+1
      goto 23045
23047 continue
      vwk(i,j) = vwk(i,j) / wtsum(m) - muwk(i) * muwk(j)
23043 j=j+1
      goto 23042
23044 continue
23040 i=i+1
      goto 23039
23041 continue
      call daxpy (nxis, bwt(m), muwk, 1, mu, 1)
      call daxpy (nxis*nxis, bwt(m), vwk, 1, v, 1)
23034 m=m+1
      goto 23033
23035 continue
      i=1
23048 if(.not.(i.le.nxi))goto 23050
      j=i
23051 if(.not.(j.le.nxi))goto 23053
      v(i,j) = v(i,j) + q(i,j)
23052 j=j+1
      goto 23051
23053 continue
23049 i=i+1
      goto 23048
23050 continue
      call daxpy (nxis, 1.d0, mrs, 1, mu, 1)
      call dsymv ('u', nxi, -1.d0, q, nxi, cd, 1, 1.d0, mu, 1)
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23054 if(.not.(i.le.nxis))goto 23056
      jpvt(i) = 0
23055 i=i+1
      goto 23054
23056 continue
      call dchdc (v, nxis, nxis, wk, jpvt, 1, rkv)
23057 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23057
      endif
23058 continue
      i=rkv+1
23059 if(.not.(i.le.nxis))goto 23061
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23060 i=i+1
      goto 23059
23061 continue
23062 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dset (nxis-rkv, 0.d0, cdnew(rkv+1), 1)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      m=1
23065 if(.not.(m.le.nt))goto 23067
      wtsumnew(m) = 0.d0
23066 m=m+1
      goto 23065
23067 continue
      i=1
23068 if(.not.(i.le.nqd))goto 23070
      tmp = ddot (nxis, qdrs(i,1), nqd, cdnew, 1)
      if(tmp.gt.3.d2)then
      flag = flag + 1
      goto 23070
      endif
      tmp = dexp (tmp)
      m=1
23073 if(.not.(m.le.nt))goto 23075
      wtnew(m,i) = qdwt(m,i) * tmp
      wtsumnew(m) = wtsumnew(m) + wtnew(m,i)
23074 m=m+1
      goto 23073
23075 continue
23069 i=i+1
      goto 23068
23070 continue
      norm = 0.d0
      m=1
23076 if(.not.(m.le.nt))goto 23078
      norm = norm + bwt(m) * dlog (wtsumnew(m))
23077 m=m+1
      goto 23076
23078 continue
      if((flag.eq.0).or.(flag.eq.2))then
      fitmean = 0.d0
      i=1
23081 if(.not.(i.le.nobs))goto 23083
      tmp = ddot (nxis, rs(1,i), 1, cdnew, 1)
      if(tmp.gt.3.d2)then
      flag = flag + 1
      goto 23083
      endif
      fitnew(i) = dexp (tmp)
      if(cntsum.ne.0)then
      tmp = tmp * dfloat (cnt(i))
      endif
      fitmean = fitmean + tmp
23082 i=i+1
      goto 23081
23083 continue
      if(cntsum.eq.0)then
      fitmean = fitmean / dfloat (nobs)
      else
      fitmean = fitmean / dfloat (cntsum)
      endif
      call dsymv ('u', nxi, 1.d0, q, nxi, cdnew, 1, 0.d0, wk, 1)
      lkhdnew = ddot (nxi, cdnew, 1, wk, 1) / 2.d0 - fitmean + norm
      endif
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      call dcopy (nt*nqd, qdwt, 1, wt, 1)
      lkhd = 0.d0
      m=1
23092 if(.not.(m.le.nt))goto 23094
      wtsum(m) = 0.d0
      i=1
23095 if(.not.(i.le.nqd))goto 23097
      wtsum(m) = wtsum(m) + wt(m,i)
23096 i=i+1
      goto 23095
23097 continue
      lkhd = lkhd + bwt(m) * dlog (wtsum(m))
23093 m=m+1
      goto 23092
23094 continue
      call dset (nobs, 1.d0, fit, 1)
      iter = 0
      goto 23064
      endif
      if(flag.eq.3)then
      goto 23064
      endif
      if(lkhdnew-lkhd.lt.1.d1*(1.d0+dabs(lkhd))*mchpr)then
      goto 23064
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23064
      endif
23063 goto 23062
23064 continue
      if(flag.eq.1)then
      flag = 2
      goto 23031
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      i=1
23108 if(.not.(i.le.nqd))goto 23110
      m=1
23111 if(.not.(m.le.nt))goto 23113
      disc = dmax1 (disc, dabs(wt(m,i)-wtnew(m,i))/(1.d0+dabs(wt(m,i))))
23112 m=m+1
      goto 23111
23113 continue
23109 i=i+1
      goto 23108
23110 continue
      i=1
23114 if(.not.(i.le.nobs))goto 23116
      disc = dmax1 (disc, dabs(fit(i)-fitnew(i))/(1.d0+dabs(fit(i))))
23115 i=i+1
      goto 23114
23116 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
      disc0 = dmax1 ((mumax/(1.d0+lkhd))**2, dabs(lkhd-lkhdnew)/(1.d0+da
     *bs(lkhd)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nqd*nt, wtnew, 1, wt, 1)
      call dcopy (nt, wtsumnew, 1, wtsum, 1)
      call dcopy (nobs, fitnew, 1, fit, 1)
      lkhd = lkhdnew
      if(disc0.lt.prec)then
      goto 23032
      endif
      if(disc.lt.prec)then
      goto 23032
      endif
      if(iter.lt.maxiter)then
      goto 23031
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      call dcopy (nt*nqd, qdwt, 1, wt, 1)
      lkhd = 0.d0
      m=1
23125 if(.not.(m.le.nt))goto 23127
      wtsum(m) = 0.d0
      i=1
23128 if(.not.(i.le.nqd))goto 23130
      wtsum(m) = wtsum(m) + wt(m,i)
23129 i=i+1
      goto 23128
23130 continue
      lkhd = lkhd + bwt(m) * dlog (wtsum(m))
23126 m=m+1
      goto 23125
23127 continue
      call dset (nobs, 1.d0, fit, 1)
      iter = 0
      flag = 2
      else
      info = 2
      goto 23032
      endif
23031 goto 23030
23032 continue
      i=1
23131 if(.not.(i.le.nobs))goto 23133
      call daxpy (nxis, -1.d0, mrs, 1, rs(1,i), 1)
      call dprmut (rs(1,i), nxis, jpvt, 0)
      if(cntsum.ne.0)then
      call dscal (nxis, dsqrt(dfloat(cnt(i))), rs(1,i), 1)
      endif
      call dtrsl (v, nxis, nxis, rs(1,i), 11, infowk)
      if(nxis-rkv.gt.0)then
      call dset (nxis-rkv, 0.d0, rs(rkv+1,i), 1)
      endif
23132 i=i+1
      goto 23131
23133 continue
      trc = ddot (nobs*nxis, rs, 1, rs, 1)
      if(cntsum.eq.0)then
      trc = trc / dfloat(nobs) / (dfloat(nobs)-1.d0)
      lkhd = 0.d0
      i=1
23140 if(.not.(i.le.nobs))goto 23142
      lkhd = lkhd + dlog (fit(i))
23141 i=i+1
      goto 23140
23142 continue
      lkhd = lkhd / dfloat (nobs)
      else
      trc = trc / dfloat(cntsum) / (dfloat(cntsum)-1.d0)
      lkhd = 0.d0
      i=1
23143 if(.not.(i.le.nobs))goto 23145
      lkhd = lkhd + dfloat (cnt(i)) * dlog (fit(i))
23144 i=i+1
      goto 23143
23145 continue
      lkhd = lkhd / dfloat (cntsum)
      endif
      m=1
23146 if(.not.(m.le.nt))goto 23148
      lkhd = lkhd - bwt(m) * dlog (wtsum(m))
23147 m=m+1
      goto 23146
23148 continue
      mrs(1) = lkhd
      mrs(2) = trc
      return
      end
