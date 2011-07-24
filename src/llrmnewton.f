C Output from Public domain Ratfor, version 1.01
      subroutine llrmnewton (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qd
     *rs, nqd, nx, xxwt, prec, maxiter, mchpr, jpvt, wk, info)
      integer nxis, nxi, nobs, cntsum, cnt(*), nqd, nx, maxiter, jpvt(*)
     *, info
      double precision cd(*), q(nxi,*), rs(nxis,*), qdrs(nqd,nxis,*), xx
     *wt(*), prec, mchpr, wk(*)
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
      call llrmnewton1 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qdrs, n
     *qd, nx, xxwt, prec, maxiter, mchpr, wk(iwt), wk(iwtsum), wk(imrs),
     * wk(ifit), wk(imu), wk(imuwk), wk(iv), wk(ivwk), jpvt, wk(icdnew),
     * wk(iwtnew), wk(iwtnewsum), wk(ifitnew), wk(iwk), info)
      return
      end
      subroutine llrmnewton1 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, q
     *drs, nqd, nx, xxwt, prec, maxiter, mchpr, wt, wtsum, mrs, fit, mu,
     * muwk, v, vwk, jpvt, cdnew, wtnew, wtnewsum, fitnew, wk, info)
      integer nxis, nxi, nobs, cntsum, cnt(*), nqd, nx, maxiter, jpvt(*)
     *, info
      double precision cd(*), q(nxi,*), rs(nxis,*), qdrs(nqd,nxis,*), xx
     *wt(*), prec, mchpr, wt(nqd,*), wtsum(*), mrs(*), fit(*), mu(*), mu
     *wk(*), v(nxis,*), vwk(nxis,*), cdnew(*), wtnew(nqd,*), wtnewsum(*)
     *, fitnew(*), wk(*)
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
      wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1))
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
      wtnew(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cdnew, 1))
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
      tmp = dfloat (nqd)
      call dset (nqd*nx, 1.d0, wt, 1)
      call dset (nx, tmp, wtsum, 1)
      call dset (nobs, 1.d0/tmp, fit, 1)
      fitmean = - dlog (tmp)
      lkhd = - fitmean
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
23086 if(.not.(kk.le.nx))goto 23088
      i=1
23089 if(.not.(i.le.nqd))goto 23091
      disc = dmax1 (disc, dabs(wt(i,kk)-wtnew(i,kk))/(1.d0+dabs(wt(i,kk)
     *)))
23090 i=i+1
      goto 23089
23091 continue
23087 kk=kk+1
      goto 23086
23088 continue
      i=1
23092 if(.not.(i.le.nobs))goto 23094
      disc = dmax1 (disc, dabs(fit(i)-fitnew(i))/(1.d0+dabs(fit(i))))
23093 i=i+1
      goto 23092
23094 continue
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
      tmp = dfloat (nqd)
      call dset (nqd*nx, 1.d0, wt, 1)
      call dset (nx, tmp, wtsum, 1)
      call dset (nobs, 1.d0/tmp, fit, 1)
      fitmean = - dlog (tmp)
      lkhd = - fitmean
      iter = 0
      flag = 2
      else
      info = 2
      goto 23026
      endif
23025 goto 23024
23026 continue
      i=1
23103 if(.not.(i.le.nobs))goto 23105
      call daxpy (nxis, -1.d0, mrs, 1, rs(1,i), 1)
      call dprmut (rs(1,i), nxis, jpvt, 0)
      if(cntsum.ne.0)then
      call dscal (nxis, dsqrt(dfloat(cnt(i))), rs(1,i), 1)
      endif
      call dtrsl (v, nxis, nxis, rs(1,i), 11, infowk)
23104 i=i+1
      goto 23103
23105 continue
      trc = ddot (nobs*nxis, rs, 1, rs, 1)
      if(cntsum.eq.0)then
      trc = trc / dfloat(nobs) / (dfloat(nobs)-1.d0)
      lkhd = 0.d0
      i=1
23110 if(.not.(i.le.nobs))goto 23112
      lkhd = lkhd + dlog (fit(i))
23111 i=i+1
      goto 23110
23112 continue
      lkhd = lkhd / dfloat (nobs)
      else
      trc = trc / dfloat(cntsum) / (dfloat(cntsum)-1.d0)
      lkhd = 0.d0
      i=1
23113 if(.not.(i.le.nobs))goto 23115
      lkhd = lkhd + dfloat (cnt(i)) * dlog (fit(i))
23114 i=i+1
      goto 23113
23115 continue
      lkhd = lkhd / dfloat (cntsum)
      endif
      kk=1
23116 if(.not.(kk.le.nx))goto 23118
      lkhd = lkhd - xxwt(kk) * dlog (wtsum(kk))
23117 kk=kk+1
      goto 23116
23118 continue
      wtsum(1) = lkhd
      wtsum(2) = trc
      return
      end
      subroutine llrmaux (cd, nxis, q, nxi, qdrs, nqd, nx, xxwt, mchpr, 
     *wt, wtsum, mu, v, vwk, jpvt)
      integer nxis, nxi, nqd, nx, jpvt(*)
      double precision cd(*), q(nxi,*), qdrs(nqd,nxis,*), xxwt(*), mchpr
     *, wt(nqd,*), wtsum(*), mu(*), v(nxis,*), vwk(nxis,*)
      integer i, j, k, kk, rkv
      double precision ddot
      kk=1
23119 if(.not.(kk.le.nx))goto 23121
      wtsum(kk) = 0.d0
      i=1
23122 if(.not.(i.le.nqd))goto 23124
      wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1))
      wtsum(kk) = wtsum(kk) + wt(i,kk)
23123 i=i+1
      goto 23122
23124 continue
23120 kk=kk+1
      goto 23119
23121 continue
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23125 if(.not.(kk.le.nx))goto 23127
      i=1
23128 if(.not.(i.le.nxis))goto 23130
      mu(i) = ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1) / wtsum(kk)
23129 i=i+1
      goto 23128
23130 continue
      i=1
23131 if(.not.(i.le.nxis))goto 23133
      j=i
23134 if(.not.(j.le.nxis))goto 23136
      vwk(i,j) = 0.d0
      k=1
23137 if(.not.(k.le.nqd))goto 23139
      vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23138 k=k+1
      goto 23137
23139 continue
      vwk(i,j) = vwk(i,j) / wtsum(kk) - mu(i) * mu(j)
23135 j=j+1
      goto 23134
23136 continue
23132 i=i+1
      goto 23131
23133 continue
      call daxpy (nxis*nxis, xxwt(kk), vwk, 1, v, 1)
23126 kk=kk+1
      goto 23125
23127 continue
      i=1
23140 if(.not.(i.le.nxi))goto 23142
      j=i
23143 if(.not.(j.le.nxi))goto 23145
      v(i,j) = v(i,j) + q(i,j)
23144 j=j+1
      goto 23143
23145 continue
23141 i=i+1
      goto 23140
23142 continue
      i=1
23146 if(.not.(i.le.nxis))goto 23148
      jpvt(i) = 0
23147 i=i+1
      goto 23146
23148 continue
      call dchdc (v, nxis, nxis, vwk, jpvt, 1, rkv)
23149 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23149
      endif
23150 continue
      i=rkv+1
23151 if(.not.(i.le.nxis))goto 23153
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23152 i=i+1
      goto 23151
23153 continue
      return
      end
      subroutine llrmrkl (cd, nxis, qdrs, nqd, nx, xxwt, wt0, offset, mc
     *hpr, wt, wtnew, mu, muwk, v, vwk, jpvt, cdnew, prec, maxiter, info
     *)
      integer nxis, nqd, nx, jpvt(*), maxiter, info
      double precision cd(*), qdrs(nqd,nxis,*), xxwt(*), wt0(nqd,*), off
     *set(nqd,*), mchpr, wt(nqd,*), wtnew(nqd,*), mu(*), muwk(*), v(nxis
     *,*), vwk(nxis,*), cdnew(*), prec
      integer i, j, k, kk, iter, flag, idamax, infowk
      double precision ddot, dasum, rkl, tmp, mumax, rklnew, disc, disc0
      kk=1
23154 if(.not.(kk.le.nx))goto 23156
      i=1
23157 if(.not.(i.le.nqd))goto 23159
      wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1) + offset(i,
     *kk))
23158 i=i+1
      goto 23157
23159 continue
      call dscal (nqd, 1.d0/dasum(nqd,wt(1,kk),1), wt(1,kk), 1)
23155 kk=kk+1
      goto 23154
23156 continue
      rkl = 0.d0
      kk=1
23160 if(.not.(kk.le.nx))goto 23162
      tmp = 0.d0
      i=1
23163 if(.not.(i.le.nqd))goto 23165
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23164 i=i+1
      goto 23163
23165 continue
      rkl = rkl + xxwt(kk) * tmp
23161 kk=kk+1
      goto 23160
23162 continue
      iter = 0
      flag = 0
23166 continue
      iter = iter + 1
      call dset (nxis, 0.d0, mu, 1)
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23169 if(.not.(kk.le.nx))goto 23171
      i=1
23172 if(.not.(i.le.nxis))goto 23174
      muwk(i) = ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1)
23173 i=i+1
      goto 23172
23174 continue
      i=1
23175 if(.not.(i.le.nxis))goto 23177
      j=i
23178 if(.not.(j.le.nxis))goto 23180
      vwk(i,j) = 0.d0
      k=1
23181 if(.not.(k.le.nqd))goto 23183
      vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23182 k=k+1
      goto 23181
23183 continue
      vwk(i,j) = vwk(i,j) - muwk(i) * muwk(j)
23179 j=j+1
      goto 23178
23180 continue
      muwk(i) = ddot (nqd, wt0(1,kk), 1, qdrs(1,i,kk), 1) - muwk(i)
23176 i=i+1
      goto 23175
23177 continue
      call daxpy (nxis, xxwt(kk), muwk, 1, mu, 1)
      call daxpy (nxis*nxis, xxwt(kk), vwk, 1, v, 1)
23170 kk=kk+1
      goto 23169
23171 continue
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23184 if(.not.(i.le.nxis))goto 23186
      jpvt(i) = 0
23185 i=i+1
      goto 23184
23186 continue
      call dmcdc (v, nxis, nxis, cdnew, jpvt, infowk)
23187 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      kk=1
23190 if(.not.(kk.le.nx))goto 23192
      i=1
23193 if(.not.(i.le.nqd))goto 23195
      wtnew(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cdnew, 1) + off
     *set(i,kk))
23194 i=i+1
      goto 23193
23195 continue
      call dscal (nqd, 1.d0/dasum(nqd,wtnew(1,kk),1), wtnew(1,kk), 1)
23191 kk=kk+1
      goto 23190
23192 continue
      if((flag.eq.0).or.(flag.eq.2))then
      rklnew = 0.d0
      kk=1
23198 if(.not.(kk.le.nx))goto 23200
      tmp = 0.d0
      i=1
23201 if(.not.(i.le.nqd))goto 23203
      tmp = tmp + dlog(wt0(i,kk)/wtnew(i,kk)) * wt0(i,kk)
23202 i=i+1
      goto 23201
23203 continue
      rklnew = rklnew + xxwt(kk) * tmp
23199 kk=kk+1
      goto 23198
23200 continue
      endif
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      kk=1
23206 if(.not.(kk.le.nx))goto 23208
      i=1
23209 if(.not.(i.le.nqd))goto 23211
      wt(i,kk) = dexp (offset(i,kk))
23210 i=i+1
      goto 23209
23211 continue
      call dscal (nqd, 1.d0/dasum(nqd,wt(1,kk),1), wt(1,kk), 1)
23207 kk=kk+1
      goto 23206
23208 continue
      rkl = 0.d0
      kk=1
23212 if(.not.(kk.le.nx))goto 23214
      tmp = 0.d0
      i=1
23215 if(.not.(i.le.nqd))goto 23217
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23216 i=i+1
      goto 23215
23217 continue
      rkl = rkl + xxwt(kk) * tmp
23213 kk=kk+1
      goto 23212
23214 continue
      iter = 0
      goto 23189
      endif
      if(flag.eq.3)then
      goto 23189
      endif
      if(rklnew-rkl.lt.1.d1*(1.d0+dabs(rkl))*mchpr)then
      goto 23189
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23189
      endif
23188 goto 23187
23189 continue
      if(flag.eq.1)then
      flag = 2
      goto 23167
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      kk=1
23228 if(.not.(kk.le.nx))goto 23230
      i=1
23231 if(.not.(i.le.nqd))goto 23233
      disc = dmax1 (disc, dabs(wt(i,kk)-wtnew(i,kk))/(1.d0+dabs(wt(i,kk)
     *)))
23232 i=i+1
      goto 23231
23233 continue
23229 kk=kk+1
      goto 23228
23230 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(rkl)))**2)
      disc0 = dmax1 ((mumax/(1.d0+rkl))**2, dabs(rkl-rklnew)/(1+dabs(rkl
     *)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nqd*nx, wtnew, 1, wt, 1)
      rkl = rklnew
      if(disc0.lt.prec)then
      goto 23168
      endif
      if(disc.lt.prec)then
      goto 23168
      endif
      if(iter.lt.maxiter)then
      goto 23167
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      call dset (nqd*nx, 1.d0/dfloat(nqd), wt, 1)
      rkl = 0.d0
      kk=1
23242 if(.not.(kk.le.nx))goto 23244
      tmp = 0.d0
      i=1
23245 if(.not.(i.le.nqd))goto 23247
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23246 i=i+1
      goto 23245
23247 continue
      rkl = rkl + xxwt(kk) * tmp
23243 kk=kk+1
      goto 23242
23244 continue
      iter = 0
      flag = 2
      else
      info = 2
      goto 23168
      endif
23167 goto 23166
23168 continue
      rkl = 0.d0
      kk=1
23248 if(.not.(kk.le.nx))goto 23250
      tmp = 0.d0
      i=1
23251 if(.not.(i.le.nqd))goto 23253
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23252 i=i+1
      goto 23251
23253 continue
      rkl = rkl + xxwt(kk) * tmp
23249 kk=kk+1
      goto 23248
23250 continue
      wt(1,1) = rkl
      return
      end
