C Output from Public domain Ratfor, version 1.01
      subroutine llrmnewton (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qd
     *rs, nqd, nx, xxwt, prec, maxiter, mchpr, jpvt, wk, info)
      integer nxis, nxi, nobs, cntsum, nqd, nx, maxiter, jpvt(*), info
      double precision cd(*), q(nxi,*), rs(nxis,*), cnt(*), qdrs(nqd,nxi
     *s,*), xxwt(*), prec, mchpr, wk(*)
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
      integer nxis, nxi, nobs, cntsum, nqd, nx, maxiter, jpvt(*), info
      double precision cd(*), q(nxi,*), rs(nxis,*), cnt(*), qdrs(nqd,nxi
     *s,*), xxwt(*), prec, mchpr, wt(nqd,*), wtsum(*), mrs(*), fit(*), m
     *u(*), muwk(*), v(nxis,*), vwk(nxis,*), cdnew(*), wtnew(nqd,*), wtn
     *ewsum(*), fitnew(*), wk(*)
      integer i, j, k, kk, iter, flag, rkv, idamax, infowk
      double precision cnt1, norm, tmp, ddot, fitmean, lkhd, mumax, lkhd
     *new, disc, disc0, trc
      info = 0
      cnt1 = 0.d0
      j=1
23000 if(.not.(j.le.nobs))goto 23002
      cnt1 = cnt1 + cnt(j)
23001 j=j+1
      goto 23000
23002 continue
      i=1
23003 if(.not.(i.le.nxis))goto 23005
      mrs(i) = 0.d0
      if(cntsum.eq.0)then
      j=1
23008 if(.not.(j.le.nobs))goto 23010
      mrs(i) = mrs(i) + rs(i,j)
23009 j=j+1
      goto 23008
23010 continue
      mrs(i) = mrs(i) / dfloat (nobs)
      else
      j=1
23011 if(.not.(j.le.nobs))goto 23013
      mrs(i) = mrs(i) + rs(i,j) * cnt(j)
23012 j=j+1
      goto 23011
23013 continue
      mrs(i) = mrs(i) / cnt1
      endif
23004 i=i+1
      goto 23003
23005 continue
      if(cntsum.eq.0)then
      trc = 1.d0 / dfloat (nobs)
      else
      trc = 1.d0 / cnt1
      endif
      norm = 0.d0
      kk=1
23016 if(.not.(kk.le.nx))goto 23018
      wtsum(kk) = 0.d0
      i=1
23019 if(.not.(i.le.nqd))goto 23021
      wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1))
      wtsum(kk) = wtsum(kk) + wt(i,kk)
23020 i=i+1
      goto 23019
23021 continue
      norm = norm + xxwt(kk) * dlog (wtsum(kk))
23017 kk=kk+1
      goto 23016
23018 continue
      fitmean = 0.d0
      i=1
23022 if(.not.(i.le.nobs))goto 23024
      tmp = ddot (nxis, rs(1,i), 1, cd, 1)
      fit(i) = dexp (tmp)
      if(cntsum.ne.0)then
      tmp = tmp * cnt(i)
      endif
      fitmean = fitmean + tmp
23023 i=i+1
      goto 23022
23024 continue
      call dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, wk, 1)
      lkhd = ddot (nxi, cd, 1, wk, 1) / 2.d0 - fitmean * trc + norm
      iter = 0
      flag = 0
23027 continue
      iter = iter + 1
      call dset (nxis, 0.d0, mu, 1)
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23030 if(.not.(kk.le.nx))goto 23032
      i=1
23033 if(.not.(i.le.nxis))goto 23035
      muwk(i) = - ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1) / wtsum(kk)
23034 i=i+1
      goto 23033
23035 continue
      i=1
23036 if(.not.(i.le.nxis))goto 23038
      j=i
23039 if(.not.(j.le.nxis))goto 23041
      vwk(i,j) = 0.d0
      k=1
23042 if(.not.(k.le.nqd))goto 23044
      vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23043 k=k+1
      goto 23042
23044 continue
      vwk(i,j) = vwk(i,j) / wtsum(kk) - muwk(i) * muwk(j)
23040 j=j+1
      goto 23039
23041 continue
23037 i=i+1
      goto 23036
23038 continue
      call daxpy (nxis, xxwt(kk), muwk, 1, mu, 1)
      call daxpy (nxis*nxis, xxwt(kk), vwk, 1, v, 1)
23031 kk=kk+1
      goto 23030
23032 continue
      i=1
23045 if(.not.(i.le.nxi))goto 23047
      j=i
23048 if(.not.(j.le.nxi))goto 23050
      v(i,j) = v(i,j) + q(i,j)
23049 j=j+1
      goto 23048
23050 continue
23046 i=i+1
      goto 23045
23047 continue
      call daxpy (nxis, 1.d0, mrs, 1, mu, 1)
      call dsymv ('u', nxi, -1.d0, q, nxi, cd, 1, 1.d0, mu, 1)
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23051 if(.not.(i.le.nxis))goto 23053
      jpvt(i) = 0
23052 i=i+1
      goto 23051
23053 continue
      call dchdc (v, nxis, nxis, wk, jpvt, 1, rkv)
23054 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23054
      endif
23055 continue
      i=rkv+1
23056 if(.not.(i.le.nxis))goto 23058
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23057 i=i+1
      goto 23056
23058 continue
23059 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dset (nxis-rkv, 0.d0, cdnew(rkv+1), 1)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      norm = 0.d0
      kk=1
23062 if(.not.(kk.le.nx))goto 23064
      wtnewsum(kk) = 0.d0
      i=1
23065 if(.not.(i.le.nqd))goto 23067
      wtnew(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cdnew, 1))
      wtnewsum(kk) = wtnewsum(kk) + wtnew(i,kk)
23066 i=i+1
      goto 23065
23067 continue
      norm = norm + xxwt(kk) * dlog (wtnewsum(kk))
23063 kk=kk+1
      goto 23062
23064 continue
      if((flag.eq.0).or.(flag.eq.2))then
      fitmean = 0.d0
      i=1
23070 if(.not.(i.le.nobs))goto 23072
      tmp = ddot (nxis, rs(1,i), 1, cdnew, 1)
      if(tmp.gt.3.d2)then
      flag = flag + 1
      goto 23072
      endif
      fitnew(i) = dexp (tmp)
      if(cntsum.ne.0)then
      tmp = tmp * cnt(i)
      endif
      fitmean = fitmean + tmp
23071 i=i+1
      goto 23070
23072 continue
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
      goto 23061
      endif
      if(flag.eq.3)then
      goto 23061
      endif
      if(lkhdnew-lkhd.lt.1.d1*(1.d0+dabs(lkhd))*mchpr)then
      goto 23061
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23061
      endif
23060 goto 23059
23061 continue
      if(flag.eq.1)then
      flag = 2
      goto 23028
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
      goto 23029
      endif
      if(disc.lt.prec)then
      goto 23029
      endif
      if(iter.lt.maxiter)then
      goto 23028
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
      goto 23029
      endif
23028 goto 23027
23029 continue
      i=1
23106 if(.not.(i.le.nobs))goto 23108
      call daxpy (nxis, -1.d0, mrs, 1, rs(1,i), 1)
      call dprmut (rs(1,i), nxis, jpvt, 0)
      if(cntsum.ne.0)then
      call dscal (nxis, dsqrt(cnt(i)), rs(1,i), 1)
      endif
      call dtrsl (v, nxis, nxis, rs(1,i), 11, infowk)
23107 i=i+1
      goto 23106
23108 continue
      trc = ddot (nobs*nxis, rs, 1, rs, 1)
      if(cntsum.eq.0)then
      trc = trc / dfloat(nobs) / (dfloat(nobs)-1.d0)
      lkhd = 0.d0
      i=1
23113 if(.not.(i.le.nobs))goto 23115
      lkhd = lkhd + dlog (fit(i))
23114 i=i+1
      goto 23113
23115 continue
      lkhd = lkhd / dfloat (nobs)
      else
      trc = trc / cnt1 / (cnt1-1.d0)
      lkhd = 0.d0
      i=1
23116 if(.not.(i.le.nobs))goto 23118
      lkhd = lkhd + cnt(i) * dlog (fit(i))
23117 i=i+1
      goto 23116
23118 continue
      lkhd = lkhd / cnt1
      endif
      kk=1
23119 if(.not.(kk.le.nx))goto 23121
      lkhd = lkhd - xxwt(kk) * dlog (wtsum(kk))
23120 kk=kk+1
      goto 23119
23121 continue
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
23122 if(.not.(kk.le.nx))goto 23124
      wtsum(kk) = 0.d0
      i=1
23125 if(.not.(i.le.nqd))goto 23127
      wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1))
      wtsum(kk) = wtsum(kk) + wt(i,kk)
23126 i=i+1
      goto 23125
23127 continue
23123 kk=kk+1
      goto 23122
23124 continue
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23128 if(.not.(kk.le.nx))goto 23130
      i=1
23131 if(.not.(i.le.nxis))goto 23133
      mu(i) = ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1) / wtsum(kk)
23132 i=i+1
      goto 23131
23133 continue
      i=1
23134 if(.not.(i.le.nxis))goto 23136
      j=i
23137 if(.not.(j.le.nxis))goto 23139
      vwk(i,j) = 0.d0
      k=1
23140 if(.not.(k.le.nqd))goto 23142
      vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23141 k=k+1
      goto 23140
23142 continue
      vwk(i,j) = vwk(i,j) / wtsum(kk) - mu(i) * mu(j)
23138 j=j+1
      goto 23137
23139 continue
23135 i=i+1
      goto 23134
23136 continue
      call daxpy (nxis*nxis, xxwt(kk), vwk, 1, v, 1)
23129 kk=kk+1
      goto 23128
23130 continue
      i=1
23143 if(.not.(i.le.nxi))goto 23145
      j=i
23146 if(.not.(j.le.nxi))goto 23148
      v(i,j) = v(i,j) + q(i,j)
23147 j=j+1
      goto 23146
23148 continue
23144 i=i+1
      goto 23143
23145 continue
      i=1
23149 if(.not.(i.le.nxis))goto 23151
      jpvt(i) = 0
23150 i=i+1
      goto 23149
23151 continue
      call dchdc (v, nxis, nxis, vwk, jpvt, 1, rkv)
23152 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23152
      endif
23153 continue
      i=rkv+1
23154 if(.not.(i.le.nxis))goto 23156
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23155 i=i+1
      goto 23154
23156 continue
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
23157 if(.not.(kk.le.nx))goto 23159
      i=1
23160 if(.not.(i.le.nqd))goto 23162
      wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1) + offset(i,
     *kk))
23161 i=i+1
      goto 23160
23162 continue
      call dscal (nqd, 1.d0/dasum(nqd,wt(1,kk),1), wt(1,kk), 1)
23158 kk=kk+1
      goto 23157
23159 continue
      rkl = 0.d0
      kk=1
23163 if(.not.(kk.le.nx))goto 23165
      tmp = 0.d0
      i=1
23166 if(.not.(i.le.nqd))goto 23168
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23167 i=i+1
      goto 23166
23168 continue
      rkl = rkl + xxwt(kk) * tmp
23164 kk=kk+1
      goto 23163
23165 continue
      iter = 0
      flag = 0
23169 continue
      iter = iter + 1
      call dset (nxis, 0.d0, mu, 1)
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23172 if(.not.(kk.le.nx))goto 23174
      i=1
23175 if(.not.(i.le.nxis))goto 23177
      muwk(i) = ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1)
23176 i=i+1
      goto 23175
23177 continue
      i=1
23178 if(.not.(i.le.nxis))goto 23180
      j=i
23181 if(.not.(j.le.nxis))goto 23183
      vwk(i,j) = 0.d0
      k=1
23184 if(.not.(k.le.nqd))goto 23186
      vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23185 k=k+1
      goto 23184
23186 continue
      vwk(i,j) = vwk(i,j) - muwk(i) * muwk(j)
23182 j=j+1
      goto 23181
23183 continue
      muwk(i) = ddot (nqd, wt0(1,kk), 1, qdrs(1,i,kk), 1) - muwk(i)
23179 i=i+1
      goto 23178
23180 continue
      call daxpy (nxis, xxwt(kk), muwk, 1, mu, 1)
      call daxpy (nxis*nxis, xxwt(kk), vwk, 1, v, 1)
23173 kk=kk+1
      goto 23172
23174 continue
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23187 if(.not.(i.le.nxis))goto 23189
      jpvt(i) = 0
23188 i=i+1
      goto 23187
23189 continue
      call dmcdc (v, nxis, nxis, cdnew, jpvt, infowk)
23190 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      kk=1
23193 if(.not.(kk.le.nx))goto 23195
      i=1
23196 if(.not.(i.le.nqd))goto 23198
      wtnew(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cdnew, 1) + off
     *set(i,kk))
23197 i=i+1
      goto 23196
23198 continue
      call dscal (nqd, 1.d0/dasum(nqd,wtnew(1,kk),1), wtnew(1,kk), 1)
23194 kk=kk+1
      goto 23193
23195 continue
      if((flag.eq.0).or.(flag.eq.2))then
      rklnew = 0.d0
      kk=1
23201 if(.not.(kk.le.nx))goto 23203
      tmp = 0.d0
      i=1
23204 if(.not.(i.le.nqd))goto 23206
      tmp = tmp + dlog(wt0(i,kk)/wtnew(i,kk)) * wt0(i,kk)
23205 i=i+1
      goto 23204
23206 continue
      rklnew = rklnew + xxwt(kk) * tmp
23202 kk=kk+1
      goto 23201
23203 continue
      endif
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      kk=1
23209 if(.not.(kk.le.nx))goto 23211
      i=1
23212 if(.not.(i.le.nqd))goto 23214
      wt(i,kk) = dexp (offset(i,kk))
23213 i=i+1
      goto 23212
23214 continue
      call dscal (nqd, 1.d0/dasum(nqd,wt(1,kk),1), wt(1,kk), 1)
23210 kk=kk+1
      goto 23209
23211 continue
      rkl = 0.d0
      kk=1
23215 if(.not.(kk.le.nx))goto 23217
      tmp = 0.d0
      i=1
23218 if(.not.(i.le.nqd))goto 23220
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23219 i=i+1
      goto 23218
23220 continue
      rkl = rkl + xxwt(kk) * tmp
23216 kk=kk+1
      goto 23215
23217 continue
      iter = 0
      goto 23192
      endif
      if(flag.eq.3)then
      goto 23192
      endif
      if(rklnew-rkl.lt.1.d1*(1.d0+dabs(rkl))*mchpr)then
      goto 23192
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23192
      endif
23191 goto 23190
23192 continue
      if(flag.eq.1)then
      flag = 2
      goto 23170
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      kk=1
23231 if(.not.(kk.le.nx))goto 23233
      i=1
23234 if(.not.(i.le.nqd))goto 23236
      disc = dmax1 (disc, dabs(wt(i,kk)-wtnew(i,kk))/(1.d0+dabs(wt(i,kk)
     *)))
23235 i=i+1
      goto 23234
23236 continue
23232 kk=kk+1
      goto 23231
23233 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(rkl)))**2)
      disc0 = dmax1 ((mumax/(1.d0+rkl))**2, dabs(rkl-rklnew)/(1+dabs(rkl
     *)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nqd*nx, wtnew, 1, wt, 1)
      rkl = rklnew
      if(disc0.lt.prec)then
      goto 23171
      endif
      if(disc.lt.prec)then
      goto 23171
      endif
      if(iter.lt.maxiter)then
      goto 23170
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      call dset (nqd*nx, 1.d0/dfloat(nqd), wt, 1)
      rkl = 0.d0
      kk=1
23245 if(.not.(kk.le.nx))goto 23247
      tmp = 0.d0
      i=1
23248 if(.not.(i.le.nqd))goto 23250
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23249 i=i+1
      goto 23248
23250 continue
      rkl = rkl + xxwt(kk) * tmp
23246 kk=kk+1
      goto 23245
23247 continue
      iter = 0
      flag = 2
      else
      info = 2
      goto 23171
      endif
23170 goto 23169
23171 continue
      rkl = 0.d0
      kk=1
23251 if(.not.(kk.le.nx))goto 23253
      tmp = 0.d0
      i=1
23254 if(.not.(i.le.nqd))goto 23256
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23255 i=i+1
      goto 23254
23256 continue
      rkl = rkl + xxwt(kk) * tmp
23252 kk=kk+1
      goto 23251
23253 continue
      wt(1,1) = rkl
      return
      end
