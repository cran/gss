C Output from Public domain Ratfor, version 1.01
      subroutine llrmnewton (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qd
     *rs, nqd, nx, xxwt, idx, prec, maxiter, mchpr, wk, info)
      integer nxis, nxi, nobs, cntsum, cnt(*), nqd, nx, idx(*), maxiter,
     * info
      double precision cd(*), q(nxi,*), rs(nxis,*), qdrs(nqd,nxis,*), xx
     *wt(*), prec, mchpr, wk(*)
      integer iwt, iwtsum, imrs, ifit, imu, imuwk, iv, ivwk, ijpvt, icdn
     *ew, iwtnew, iwtnewsum, ifitnew, iwk
      iwt = 1
      iwtsum = iwt + nqd*nx
      imrs = iwtsum + nx
      ifit = imrs + nxis
      imu = ifit + nobs
      imuwk = imu + nxis
      iv = imuwk + nxis
      ivwk = iv + nxis*nxis
      ijpvt = ivwk + nxis*nxis
      icdnew = ijpvt + nxis
      iwtnew = icdnew + nxis
      iwtnewsum = iwtnew + nqd*nx
      ifitnew = iwtnewsum + nx
      iwk = ifitnew + nobs
      call llrmnewton1 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qdrs, n
     *qd, nx, xxwt, idx, prec, maxiter, mchpr, wk(iwt), wk(iwtsum), wk(i
     *mrs), wk(ifit), wk(imu), wk(imuwk), wk(iv), wk(ivwk), wk(ijpvt), w
     *k(icdnew), wk(iwtnew), wk(iwtnewsum), wk(ifitnew), wk(iwk), info)
      return
      end
      subroutine llrmnewton1 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, q
     *drs, nqd, nx, xxwt, idx, prec, maxiter, mchpr, wt, wtsum, mrs, fit
     *, mu, muwk, v, vwk, jpvt, cdnew, wtnew, wtnewsum, fitnew, wk, info
     *)
      integer nxis, nxi, nobs, cntsum, cnt(*), nqd, nx, idx(*), maxiter,
     * jpvt(*), info
      double precision cd(*), q(nxi,*), rs(nxis,*), qdrs(nqd,nxis,*), xx
     *wt(*), prec, mchpr, wt(nqd,*), wtsum(*), mrs(*), fit(*), mu(*), mu
     *wk(*), v(nxis,*), vwk(nxis,*), cdnew(*), wtnew(nqd,*), wtnewsum(*)
     *, fitnew(*), wk(*)
      integer i, j, k, kk, iter, flag, rkv, idamax, infowk
      double precision tmp, ddot, fitmean, lkhd, mumax, lkhdnew, disc, d
     *isc0, trc
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
23014 kk=kk+1
      goto 23013
23015 continue
      fitmean = 0.d0
      i=1
23019 if(.not.(i.le.nobs))goto 23021
      tmp = ddot (nxis, rs(1,i), 1, cd, 1) - dlog (wtsum(idx(i)))
      fit(i) = dexp (tmp)
      if(cntsum.ne.0)then
      tmp = tmp * dfloat (cnt(i))
      endif
      fitmean = fitmean + tmp
23020 i=i+1
      goto 23019
23021 continue
      call dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, wk, 1)
      lkhd = lkhd + ddot (nxi, cd, 1, wk, 1) / 2.d0 - fitmean * trc
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
23060 kk=kk+1
      goto 23059
23061 continue
      if((flag.eq.0).or.(flag.eq.2))then
      fitmean = 0.d0
      i=1
23067 if(.not.(i.le.nobs))goto 23069
      tmp = ddot (nxis, rs(1,i), 1, cdnew, 1) - dlog (wtnewsum(idx(i)))
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
      lkhdnew = ddot (nxi, cdnew, 1, wk, 1) / 2.d0 - fitmean * trc
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
23116 if(.not.(kk.le.nx))goto 23118
      wtsum(kk) = 0.d0
      i=1
23119 if(.not.(i.le.nqd))goto 23121
      wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1))
      wtsum(kk) = wtsum(kk) + wt(i,kk)
23120 i=i+1
      goto 23119
23121 continue
23117 kk=kk+1
      goto 23116
23118 continue
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23122 if(.not.(kk.le.nx))goto 23124
      i=1
23125 if(.not.(i.le.nxis))goto 23127
      mu(i) = ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1) / wtsum(kk)
23126 i=i+1
      goto 23125
23127 continue
      i=1
23128 if(.not.(i.le.nxis))goto 23130
      j=i
23131 if(.not.(j.le.nxis))goto 23133
      vwk(i,j) = 0.d0
      k=1
23134 if(.not.(k.le.nqd))goto 23136
      vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23135 k=k+1
      goto 23134
23136 continue
      vwk(i,j) = vwk(i,j) / wtsum(kk) - mu(i) * mu(j)
23132 j=j+1
      goto 23131
23133 continue
23129 i=i+1
      goto 23128
23130 continue
      call daxpy (nxis*nxis, xxwt(kk), vwk, 1, v, 1)
23123 kk=kk+1
      goto 23122
23124 continue
      i=1
23137 if(.not.(i.le.nxi))goto 23139
      j=i
23140 if(.not.(j.le.nxi))goto 23142
      v(i,j) = v(i,j) + q(i,j)
23141 j=j+1
      goto 23140
23142 continue
23138 i=i+1
      goto 23137
23139 continue
      i=1
23143 if(.not.(i.le.nxis))goto 23145
      jpvt(i) = 0
23144 i=i+1
      goto 23143
23145 continue
      call dchdc (v, nxis, nxis, vwk, jpvt, 1, rkv)
23146 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23146
      endif
23147 continue
      i=rkv+1
23148 if(.not.(i.le.nxis))goto 23150
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23149 i=i+1
      goto 23148
23150 continue
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
23151 if(.not.(kk.le.nx))goto 23153
      i=1
23154 if(.not.(i.le.nqd))goto 23156
      wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1) + offset(i,
     *kk))
23155 i=i+1
      goto 23154
23156 continue
      call dscal (nqd, 1.d0/dasum(nqd,wt(1,kk),1), wt(1,kk), 1)
23152 kk=kk+1
      goto 23151
23153 continue
      rkl = 0.d0
      kk=1
23157 if(.not.(kk.le.nx))goto 23159
      tmp = 0.d0
      i=1
23160 if(.not.(i.le.nqd))goto 23162
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23161 i=i+1
      goto 23160
23162 continue
      rkl = rkl + xxwt(kk) * tmp
23158 kk=kk+1
      goto 23157
23159 continue
      iter = 0
      flag = 0
23163 continue
      iter = iter + 1
      call dset (nxis, 0.d0, mu, 1)
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23166 if(.not.(kk.le.nx))goto 23168
      i=1
23169 if(.not.(i.le.nxis))goto 23171
      muwk(i) = ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1)
23170 i=i+1
      goto 23169
23171 continue
      i=1
23172 if(.not.(i.le.nxis))goto 23174
      j=i
23175 if(.not.(j.le.nxis))goto 23177
      vwk(i,j) = 0.d0
      k=1
23178 if(.not.(k.le.nqd))goto 23180
      vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23179 k=k+1
      goto 23178
23180 continue
      vwk(i,j) = vwk(i,j) - muwk(i) * muwk(j)
23176 j=j+1
      goto 23175
23177 continue
      muwk(i) = ddot (nqd, wt0(1,kk), 1, qdrs(1,i,kk), 1) - muwk(i)
23173 i=i+1
      goto 23172
23174 continue
      call daxpy (nxis, xxwt(kk), muwk, 1, mu, 1)
      call daxpy (nxis*nxis, xxwt(kk), vwk, 1, v, 1)
23167 kk=kk+1
      goto 23166
23168 continue
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23181 if(.not.(i.le.nxis))goto 23183
      jpvt(i) = 0
23182 i=i+1
      goto 23181
23183 continue
      call dmcdc (v, nxis, nxis, cdnew, jpvt, infowk)
23184 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      kk=1
23187 if(.not.(kk.le.nx))goto 23189
      i=1
23190 if(.not.(i.le.nqd))goto 23192
      wtnew(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cdnew, 1) + off
     *set(i,kk))
23191 i=i+1
      goto 23190
23192 continue
      call dscal (nqd, 1.d0/dasum(nqd,wtnew(1,kk),1), wtnew(1,kk), 1)
23188 kk=kk+1
      goto 23187
23189 continue
      if((flag.eq.0).or.(flag.eq.2))then
      rklnew = 0.d0
      kk=1
23195 if(.not.(kk.le.nx))goto 23197
      tmp = 0.d0
      i=1
23198 if(.not.(i.le.nqd))goto 23200
      tmp = tmp + dlog(wt0(i,kk)/wtnew(i,kk)) * wt0(i,kk)
23199 i=i+1
      goto 23198
23200 continue
      rklnew = rklnew + xxwt(kk) * tmp
23196 kk=kk+1
      goto 23195
23197 continue
      endif
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      call dset (nqd*nx, 1.d0/dfloat(nqd), wt, 1)
      rkl = 0.d0
      kk=1
23203 if(.not.(kk.le.nx))goto 23205
      tmp = 0.d0
      i=1
23206 if(.not.(i.le.nqd))goto 23208
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23207 i=i+1
      goto 23206
23208 continue
      rkl = rkl + xxwt(kk) * tmp
23204 kk=kk+1
      goto 23203
23205 continue
      iter = 0
      goto 23186
      endif
      if(flag.eq.3)then
      goto 23186
      endif
      if(rklnew-rkl.lt.1.d1*(1.d0+dabs(rkl))*mchpr)then
      goto 23186
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23186
      endif
23185 goto 23184
23186 continue
      if(flag.eq.1)then
      flag = 2
      goto 23164
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      kk=1
23219 if(.not.(kk.le.nx))goto 23221
      i=1
23222 if(.not.(i.le.nqd))goto 23224
      disc = dmax1 (disc, dabs(wt(i,kk)-wtnew(i,kk))/(1.d0+dabs(wt(i,kk)
     *)))
23223 i=i+1
      goto 23222
23224 continue
23220 kk=kk+1
      goto 23219
23221 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(rkl)))**2)
      disc0 = dmax1 ((mumax/(1.d0+rkl))**2, dabs(rkl-rklnew)/(1+dabs(rkl
     *)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nqd*nx, wtnew, 1, wt, 1)
      rkl = rklnew
      if(disc0.lt.prec)then
      goto 23165
      endif
      if(disc.lt.prec)then
      goto 23165
      endif
      if(iter.lt.maxiter)then
      goto 23164
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      call dset (nqd*nx, 1.d0/dfloat(nqd), wt, 1)
      rkl = 0.d0
      kk=1
23233 if(.not.(kk.le.nx))goto 23235
      tmp = 0.d0
      i=1
23236 if(.not.(i.le.nqd))goto 23238
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23237 i=i+1
      goto 23236
23238 continue
      rkl = rkl + xxwt(kk) * tmp
23234 kk=kk+1
      goto 23233
23235 continue
      iter = 0
      flag = 2
      else
      info = 2
      goto 23165
      endif
23164 goto 23163
23165 continue
      rkl = 0.d0
      kk=1
23239 if(.not.(kk.le.nx))goto 23241
      tmp = 0.d0
      i=1
23242 if(.not.(i.le.nqd))goto 23244
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23243 i=i+1
      goto 23242
23244 continue
      rkl = rkl + xxwt(kk) * tmp
23240 kk=kk+1
      goto 23239
23241 continue
      wt(1,1) = rkl
      return
      end
