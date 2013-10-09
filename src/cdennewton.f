C Output from Public domain Ratfor, version 1.01
      subroutine cdennewton (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qd
     *rs, nqd, nx, xxwt, qdwt, prec, maxiter, mchpr, jpvt, wk, info)
      integer nxis, nxi, nobs, cntsum, cnt(*), nqd, nx, maxiter, jpvt(*)
     *, info
      double precision cd(*), q(nxi,*), rs(nxis,*), qdrs(nqd,nxis,*), xx
     *wt(*), qdwt(*), prec, mchpr, wk(*)
      integer iwt, iwtsum, imrs, ifit, imu, imuwk, iv, ivwk, icdnew, iwt
     *new, iwtnewsum, ifitnew, iwk
      iwt = 1
      iwtsum = iwt + nqd*nx
      imrs = iwtsum + nx
      ifit = imrs + nxis
      imu = ifit + nobs
      imuwk = imu + nxis
      iv = imuwk + nxis
      ivwk = iv + nxis*nxis
      icdnew = ivwk + nxis*nxis
      iwtnew = icdnew + nxis
      iwtnewsum = iwtnew + nqd*nx
      ifitnew = iwtnewsum + nx
      iwk = ifitnew + nobs
      call cdennewton1 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qdrs, n
     *qd, nx, xxwt, qdwt, prec, maxiter, mchpr, wk(iwt), wk(iwtsum), wk(
     *imrs), wk(ifit), wk(imu), wk(imuwk), wk(iv), wk(ivwk), jpvt, wk(ic
     *dnew), wk(iwtnew), wk(iwtnewsum), wk(ifitnew), wk(iwk), info)
      return
      end
      subroutine cdennewton1 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, q
     *drs, nqd, nx, xxwt, qdwt, prec, maxiter, mchpr, wt, wtsum, mrs, fi
     *t, mu, muwk, v, vwk, jpvt, cdnew, wtnew, wtnewsum, fitnew, wk, inf
     *o)
      integer nxis, nxi, nobs, cntsum, cnt(*), nqd, nx, maxiter, jpvt(*)
     *, info
      double precision cd(*), q(nxi,*), rs(nxis,*), qdrs(nqd,nxis,*), xx
     *wt(*), qdwt(*), prec, mchpr, wt(nqd,*), wtsum(*), mrs(*), fit(*), 
     *mu(*), muwk(*), v(nxis,*), vwk(nxis,*), cdnew(*), wtnew(nqd,*), wt
     *newsum(*), fitnew(*), wk(*)
      integer i, j, k, kk, iter, flag, rkv, idamax, infowk
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
      if(cntsum.eq.0)then
      trc = 1.d0 / dfloat (nobs)
      else
      trc = 1.d0 / dfloat (cntsum)
      endif
      norm = 0.d0
      kk=1
23013 if(.not.(kk.le.nx))goto 23015
      wtsum(kk) = 0.d0
      i=1
23016 if(.not.(i.le.nqd))goto 23018
      wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1)) * qdwt(i)
      wtsum(kk) = wtsum(kk) + wt(i,kk)
23017 i=i+1
      goto 23016
23018 continue
      norm = norm + xxwt(kk) * dlog (wtsum(kk))
23014 kk=kk+1
      goto 23013
23015 continue
      fitmean = 0.d0
      i=1
23019 if(.not.(i.le.nobs))goto 23021
      tmp = ddot (nxis, rs(1,i), 1, cd, 1)
      fit(i) = dexp (tmp)
      if(cntsum.ne.0)then
      tmp = tmp * dfloat (cnt(i))
      endif
      fitmean = fitmean + tmp
23020 i=i+1
      goto 23019
23021 continue
      call dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, wk, 1)
      lkhd = ddot (nxi, cd, 1, wk, 1) / 2.d0 - fitmean * trc + norm
      iter = 0
      flag = 0
23024 continue
      iter = iter + 1
      call dset (nxis, 0.d0, mu, 1)
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23027 if(.not.(kk.le.nx))goto 23029
      i=1
23030 if(.not.(i.le.nxis))goto 23032
      muwk(i) = - ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1) / wtsum(kk)
23031 i=i+1
      goto 23030
23032 continue
      i=1
23033 if(.not.(i.le.nxis))goto 23035
      j=i
23036 if(.not.(j.le.nxis))goto 23038
      vwk(i,j) = 0.d0
      k=1
23039 if(.not.(k.le.nqd))goto 23041
      vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23040 k=k+1
      goto 23039
23041 continue
      vwk(i,j) = vwk(i,j) / wtsum(kk) - muwk(i) * muwk(j)
23037 j=j+1
      goto 23036
23038 continue
23034 i=i+1
      goto 23033
23035 continue
      call daxpy (nxis, xxwt(kk), muwk, 1, mu, 1)
      call daxpy (nxis*nxis, xxwt(kk), vwk, 1, v, 1)
23028 kk=kk+1
      goto 23027
23029 continue
      i=1
23042 if(.not.(i.le.nxi))goto 23044
      j=i
23045 if(.not.(j.le.nxi))goto 23047
      v(i,j) = v(i,j) + q(i,j)
23046 j=j+1
      goto 23045
23047 continue
23043 i=i+1
      goto 23042
23044 continue
      call daxpy (nxis, 1.d0, mrs, 1, mu, 1)
      call dsymv ('u', nxi, -1.d0, q, nxi, cd, 1, 1.d0, mu, 1)
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23048 if(.not.(i.le.nxis))goto 23050
      jpvt(i) = 0
23049 i=i+1
      goto 23048
23050 continue
      call dchdc (v, nxis, nxis, wk, jpvt, 1, rkv)
23051 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23051
      endif
23052 continue
      i=rkv+1
23053 if(.not.(i.le.nxis))goto 23055
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23054 i=i+1
      goto 23053
23055 continue
23056 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dset (nxis-rkv, 0.d0, cdnew(rkv+1), 1)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      norm = 0.d0
      kk=1
23059 if(.not.(kk.le.nx))goto 23061
      wtnewsum(kk) = 0.d0
      i=1
23062 if(.not.(i.le.nqd))goto 23064
      wtnew(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cdnew, 1)) * qd
     *wt(i)
      wtnewsum(kk) = wtnewsum(kk) + wtnew(i,kk)
23063 i=i+1
      goto 23062
23064 continue
      norm = norm + xxwt(kk) * dlog (wtnewsum(kk))
23060 kk=kk+1
      goto 23059
23061 continue
      if((flag.eq.0).or.(flag.eq.2))then
      fitmean = 0.d0
      i=1
23067 if(.not.(i.le.nobs))goto 23069
      tmp = ddot (nxis, rs(1,i), 1, cdnew, 1)
      if(tmp.gt.3.d2)then
      flag = flag + 1
      goto 23069
      endif
      fitnew(i) = dexp (tmp)
      if(cntsum.ne.0)then
      tmp = tmp * dfloat (cnt(i))
      endif
      fitmean = fitmean + tmp
23068 i=i+1
      goto 23067
23069 continue
      call dsymv ('u', nxi, 1.d0, q, nxi, cdnew, 1, 0.d0, wk, 1)
      lkhdnew = ddot (nxi, cdnew, 1, wk, 1) / 2.d0 - fitmean * trc + nor
     *m
      endif
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      kk=1
23076 if(.not.(kk.le.nx))goto 23078
      call dcopy (nqd, qdwt, 1, wt(1,kk), 1)
23077 kk=kk+1
      goto 23076
23078 continue
      call dset (nx, 1.d0, wtsum, 1)
      call dset (nobs, 1.d0, fit, 1)
      lkhd = 0.d0
      iter = 0
      goto 23058
      endif
      if(flag.eq.3)then
      goto 23058
      endif
      if(lkhdnew-lkhd.lt.1.d1*(1.d0+dabs(lkhd))*mchpr)then
      goto 23058
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23058
      endif
23057 goto 23056
23058 continue
      if(flag.eq.1)then
      flag = 2
      goto 23025
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      kk=1
23089 if(.not.(kk.le.nx))goto 23091
      i=1
23092 if(.not.(i.le.nqd))goto 23094
      disc = dmax1 (disc, dabs(wt(i,kk)-wtnew(i,kk))/(1.d0+dabs(wt(i,kk)
     *)))
23093 i=i+1
      goto 23092
23094 continue
23090 kk=kk+1
      goto 23089
23091 continue
      i=1
23095 if(.not.(i.le.nobs))goto 23097
      disc = dmax1 (disc, dabs(fit(i)-fitnew(i))/(1.d0+dabs(fit(i))))
23096 i=i+1
      goto 23095
23097 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
      disc0 = dmax1 ((mumax/(1.d0+lkhd))**2, dabs(lkhd-lkhdnew)/(1+dabs(
     *lkhd)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nqd*nx, wtnew, 1, wt, 1)
      call dcopy (nx, wtnewsum, 1, wtsum, 1)
      call dcopy (nobs, fitnew, 1, fit, 1)
      lkhd = lkhdnew
      if(disc0.lt.prec)then
      goto 23026
      endif
      if(disc.lt.prec)then
      goto 23026
      endif
      if(iter.lt.maxiter)then
      goto 23025
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      kk=1
23106 if(.not.(kk.le.nx))goto 23108
      call dcopy (nqd, qdwt, 1, wt(1,kk), 1)
23107 kk=kk+1
      goto 23106
23108 continue
      call dset (nx, 1.d0, wtsum, 1)
      call dset (nobs, 1.d0, fit, 1)
      lkhd = 0.d0
      iter = 0
      flag = 2
      else
      info = 2
      goto 23026
      endif
23025 goto 23024
23026 continue
      i=1
23109 if(.not.(i.le.nobs))goto 23111
      call daxpy (nxis, -1.d0, mrs, 1, rs(1,i), 1)
      call dprmut (rs(1,i), nxis, jpvt, 0)
      if(cntsum.ne.0)then
      call dscal (nxis, dsqrt(dfloat(cnt(i))), rs(1,i), 1)
      endif
      call dtrsl (v, nxis, nxis, rs(1,i), 11, infowk)
23110 i=i+1
      goto 23109
23111 continue
      trc = ddot (nobs*nxis, rs, 1, rs, 1)
      if(cntsum.eq.0)then
      trc = trc / dfloat(nobs) / (dfloat(nobs)-1.d0)
      lkhd = 0.d0
      i=1
23116 if(.not.(i.le.nobs))goto 23118
      lkhd = lkhd + dlog (fit(i))
23117 i=i+1
      goto 23116
23118 continue
      lkhd = lkhd / dfloat (nobs)
      else
      trc = trc / dfloat(cntsum) / (dfloat(cntsum)-1.d0)
      lkhd = 0.d0
      i=1
23119 if(.not.(i.le.nobs))goto 23121
      lkhd = lkhd + dfloat (cnt(i)) * dlog (fit(i))
23120 i=i+1
      goto 23119
23121 continue
      lkhd = lkhd / dfloat (cntsum)
      endif
      kk=1
23122 if(.not.(kk.le.nx))goto 23124
      lkhd = lkhd - xxwt(kk) * dlog (wtsum(kk))
23123 kk=kk+1
      goto 23122
23124 continue
      wtsum(1) = lkhd
      wtsum(2) = trc
      return
      end
      subroutine cdenrkl (cd, nxis, qdrs, nqd, nx, xxwt, qdwt, wt0, mchp
     *r, wt, wtnew, mu, muwk, v, vwk, jpvt, cdnew, prec, maxiter, info)
      integer nxis, nqd, nx, jpvt(*), maxiter, info
      double precision cd(*), qdrs(nqd,nxis,*), xxwt(*), qdwt(*), wt0(nq
     *d,*), mchpr, wt(nqd,*), wtnew(nqd,*), mu(*), muwk(*), v(nxis,*), v
     *wk(nxis,*), cdnew(*), prec
      integer i, j, k, kk, iter, flag, idamax, infowk
      double precision ddot, dasum, rkl, tmp, mumax, rklnew, disc, disc0
      kk=1
23125 if(.not.(kk.le.nx))goto 23127
      i=1
23128 if(.not.(i.le.nqd))goto 23130
      wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1)) * qdwt(i)
23129 i=i+1
      goto 23128
23130 continue
      call dscal (nqd, 1.d0/dasum(nqd,wt(1,kk),1), wt(1,kk), 1)
23126 kk=kk+1
      goto 23125
23127 continue
      rkl = 0.d0
      kk=1
23131 if(.not.(kk.le.nx))goto 23133
      tmp = 0.d0
      i=1
23134 if(.not.(i.le.nqd))goto 23136
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23135 i=i+1
      goto 23134
23136 continue
      rkl = rkl + xxwt(kk) * tmp
23132 kk=kk+1
      goto 23131
23133 continue
      iter = 0
      flag = 0
23137 continue
      iter = iter + 1
      call dset (nxis, 0.d0, mu, 1)
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23140 if(.not.(kk.le.nx))goto 23142
      i=1
23143 if(.not.(i.le.nxis))goto 23145
      muwk(i) = ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1)
23144 i=i+1
      goto 23143
23145 continue
      i=1
23146 if(.not.(i.le.nxis))goto 23148
      j=i
23149 if(.not.(j.le.nxis))goto 23151
      vwk(i,j) = 0.d0
      k=1
23152 if(.not.(k.le.nqd))goto 23154
      vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23153 k=k+1
      goto 23152
23154 continue
      vwk(i,j) = vwk(i,j) - muwk(i) * muwk(j)
23150 j=j+1
      goto 23149
23151 continue
      muwk(i) = ddot (nqd, wt0(1,kk), 1, qdrs(1,i,kk), 1) - muwk(i)
23147 i=i+1
      goto 23146
23148 continue
      call daxpy (nxis, xxwt(kk), muwk, 1, mu, 1)
      call daxpy (nxis*nxis, xxwt(kk), vwk, 1, v, 1)
23141 kk=kk+1
      goto 23140
23142 continue
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23155 if(.not.(i.le.nxis))goto 23157
      jpvt(i) = 0
23156 i=i+1
      goto 23155
23157 continue
      call dmcdc (v, nxis, nxis, cdnew, jpvt, infowk)
23158 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      kk=1
23161 if(.not.(kk.le.nx))goto 23163
      i=1
23164 if(.not.(i.le.nqd))goto 23166
      wtnew(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cdnew, 1)) * qd
     *wt(i)
23165 i=i+1
      goto 23164
23166 continue
      call dscal (nqd, 1.d0/dasum(nqd,wtnew(1,kk),1), wtnew(1,kk), 1)
23162 kk=kk+1
      goto 23161
23163 continue
      if((flag.eq.0).or.(flag.eq.2))then
      rklnew = 0.d0
      kk=1
23169 if(.not.(kk.le.nx))goto 23171
      tmp = 0.d0
      i=1
23172 if(.not.(i.le.nqd))goto 23174
      tmp = tmp + dlog(wt0(i,kk)/wtnew(i,kk)) * wt0(i,kk)
23173 i=i+1
      goto 23172
23174 continue
      rklnew = rklnew + xxwt(kk) * tmp
23170 kk=kk+1
      goto 23169
23171 continue
      endif
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      kk=1
23177 if(.not.(kk.le.nx))goto 23179
      call dcopy (nqd, qdwt, 1, wt(1,kk), 1)
      call dscal (nqd, 1.d0/dasum(nqd,wt(1,kk),1), wt(1,kk), 1)
23178 kk=kk+1
      goto 23177
23179 continue
      rkl = 0.d0
      kk=1
23180 if(.not.(kk.le.nx))goto 23182
      tmp = 0.d0
      i=1
23183 if(.not.(i.le.nqd))goto 23185
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23184 i=i+1
      goto 23183
23185 continue
      rkl = rkl + xxwt(kk) * tmp
23181 kk=kk+1
      goto 23180
23182 continue
      iter = 0
      goto 23160
      endif
      if(flag.eq.3)then
      goto 23160
      endif
      if(rklnew-rkl.lt.1.d1*(1.d0+dabs(rkl))*mchpr)then
      goto 23160
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23160
      endif
23159 goto 23158
23160 continue
      if(flag.eq.1)then
      flag = 2
      goto 23138
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      kk=1
23196 if(.not.(kk.le.nx))goto 23198
      i=1
23199 if(.not.(i.le.nqd))goto 23201
      disc = dmax1 (disc, dabs(wt(i,kk)-wtnew(i,kk))/(1.d0+dabs(wt(i,kk)
     *)))
23200 i=i+1
      goto 23199
23201 continue
23197 kk=kk+1
      goto 23196
23198 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(rkl)))**2)
      disc0 = dmax1 ((mumax/(1.d0+rkl))**2, dabs(rkl-rklnew)/(1+dabs(rkl
     *)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nqd*nx, wtnew, 1, wt, 1)
      rkl = rklnew
      if(disc0.lt.prec)then
      goto 23139
      endif
      if(disc.lt.prec)then
      goto 23139
      endif
      if(iter.lt.maxiter)then
      goto 23138
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      call dset (nqd*nx, 1.d0/dfloat(nqd), wt, 1)
      rkl = 0.d0
      kk=1
23210 if(.not.(kk.le.nx))goto 23212
      tmp = 0.d0
      i=1
23213 if(.not.(i.le.nqd))goto 23215
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23214 i=i+1
      goto 23213
23215 continue
      rkl = rkl + xxwt(kk) * tmp
23211 kk=kk+1
      goto 23210
23212 continue
      iter = 0
      flag = 2
      else
      info = 2
      goto 23139
      endif
23138 goto 23137
23139 continue
      rkl = 0.d0
      kk=1
23216 if(.not.(kk.le.nx))goto 23218
      tmp = 0.d0
      i=1
23219 if(.not.(i.le.nqd))goto 23221
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23220 i=i+1
      goto 23219
23221 continue
      rkl = rkl + xxwt(kk) * tmp
23217 kk=kk+1
      goto 23216
23218 continue
      wt(1,1) = rkl
      return
      end
