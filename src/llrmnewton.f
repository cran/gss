C Output from Public domain Ratfor, version 1.04
      subroutine llrmnewton (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qd
     *rs, nqd, nx, xxwt, qdwt, prec, maxiter, mchpr, jpvt, wk, info)
      integer nxis, nxi, nobs, nqd, nx, maxiter, jpvt(*), info
      double precision cd(*), q(nxi,*), rs(nxis,*), cntsum, cnt(*), qdrs
     *(nqd,nxis,*), xxwt(*), qdwt(*), prec, mchpr, wk(*)
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
     *qd, nx, xxwt, qdwt, prec, maxiter, mchpr, wk(iwt), wk(iwtsum), wk(
     *imrs), wk(ifit), wk(imu), wk(imuwk), wk(iv), wk(ivwk), jpvt, wk(ic
     *dnew), wk(iwtnew), wk(iwtnewsum), wk(ifitnew), wk(iwk), info)
      return
      end
      subroutine llrmnewton1 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, q
     *drs, nqd, nx, xxwt, qdwt, prec, maxiter, mchpr, wt, wtsum, mrs, fi
     *t, mu, muwk, v, vwk, jpvt, cdnew, wtnew, wtnewsum, fitnew, wk, inf
     *o)
      integer nxis, nxi, nobs, nqd, nx, maxiter, jpvt(*), info
      double precision cd(*), q(nxi,*), rs(nxis,*), cntsum, cnt(*), qdrs
     *(nqd,nxis,*), xxwt(*), qdwt(*), prec, mchpr, wt(nqd,*), wtsum(*), 
     *mrs(*), fit(*), mu(*), muwk(*), v(nxis,*), vwk(nxis,*), cdnew(*), 
     *wtnew(nqd,*), wtnewsum(*), fitnew(*), wk(*)
      integer i, j, k, kk, iter, flag, rkv, idamax, infowk
      double precision norm, tmp, ddot, fitmean, lkhd, mumax, lkhdnew, d
     *isc, disc0, trc
      info = 0
      i=1
23000 if(.not.(i.le.nxis))goto 23002
      mrs(i) = 0.d0
      if(.not.(cntsum.gt.0.d0))then
      j=1
23005 if(.not.(j.le.nobs))goto 23007
      mrs(i) = mrs(i) + rs(i,j)
23006 j=j+1
      goto 23005
23007 continue
      mrs(i) = mrs(i) / dble (nobs)
      else
      j=1
23008 if(.not.(j.le.nobs))goto 23010
      mrs(i) = mrs(i) + rs(i,j) * cnt(j)
23009 j=j+1
      goto 23008
23010 continue
      mrs(i) = mrs(i) / cntsum
      endif
23001 i=i+1
      goto 23000
23002 continue
      if(.not.(cntsum.gt.0.d0))then
      trc = 1.d0 / dble (nobs)
      else
      trc = 1.d0 / cntsum
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
      if(cntsum.gt.0.d0)then
      tmp = tmp * cnt(i)
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
      if(cntsum.gt.0.d0)then
      tmp = tmp * cnt(i)
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
      tmp = 0.d0
      i=1
23079 if(.not.(i.le.nqd))goto 23081
      tmp = tmp + qdwt(i)
23080 i=i+1
      goto 23079
23081 continue
      call dset (nx, tmp, wtsum, 1)
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
23092 if(.not.(kk.le.nx))goto 23094
      i=1
23095 if(.not.(i.le.nqd))goto 23097
      disc = dmax1 (disc, dabs(wt(i,kk)-wtnew(i,kk))/(1.d0+dabs(wt(i,kk)
     *)))
23096 i=i+1
      goto 23095
23097 continue
23093 kk=kk+1
      goto 23092
23094 continue
      i=1
23098 if(.not.(i.le.nobs))goto 23100
      disc = dmax1 (disc, dabs(fit(i)-fitnew(i))/(1.d0+dabs(fit(i))))
23099 i=i+1
      goto 23098
23100 continue
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
23109 if(.not.(kk.le.nx))goto 23111
      call dcopy (nqd, qdwt, 1, wt(1,kk), 1)
23110 kk=kk+1
      goto 23109
23111 continue
      tmp = 0.d0
      i=1
23112 if(.not.(i.le.nqd))goto 23114
      tmp = tmp + qdwt(i)
23113 i=i+1
      goto 23112
23114 continue
      call dset (nx, tmp, wtsum, 1)
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
23115 if(.not.(i.le.nobs))goto 23117
      call daxpy (nxis, -1.d0, mrs, 1, rs(1,i), 1)
      call dprmut (rs(1,i), nxis, jpvt, 0)
      if(cntsum.gt.0.d0)then
      call dscal (nxis, dsqrt(cnt(i)), rs(1,i), 1)
      endif
      call dtrsl (v, nxis, nxis, rs(1,i), 11, infowk)
23116 i=i+1
      goto 23115
23117 continue
      trc = ddot (nobs*nxis, rs, 1, rs, 1)
      if(.not.(cntsum.gt.0.d0))then
      trc = trc / dble(nobs) / (dble(nobs)-1.d0)
      lkhd = 0.d0
      i=1
23122 if(.not.(i.le.nobs))goto 23124
      lkhd = lkhd + dlog (fit(i))
23123 i=i+1
      goto 23122
23124 continue
      lkhd = lkhd / dble (nobs)
      else
      trc = trc / cntsum / (cntsum-1.d0)
      lkhd = 0.d0
      i=1
23125 if(.not.(i.le.nobs))goto 23127
      lkhd = lkhd + cnt(i) * dlog (fit(i))
23126 i=i+1
      goto 23125
23127 continue
      lkhd = lkhd / cntsum
      endif
      kk=1
23128 if(.not.(kk.le.nx))goto 23130
      lkhd = lkhd - xxwt(kk) * dlog (wtsum(kk))
23129 kk=kk+1
      goto 23128
23130 continue
      wtsum(1) = lkhd
      wtsum(2) = trc
      return
      end
      subroutine llrmaux (cd, nxis, q, nxi, qdrs, nqd, nx, xxwt, qdwt, m
     *chpr, wt, wtsum, mu, v, vwk, jpvt)
      integer nxis, nxi, nqd, nx, jpvt(*)
      double precision cd(*), q(nxi,*), qdrs(nqd,nxis,*), xxwt(*), qdwt(
     **), mchpr, wt(nqd,*), wtsum(*), mu(*), v(nxis,*), vwk(nxis,*)
      integer i, j, k, kk, rkv
      double precision ddot
      kk=1
23131 if(.not.(kk.le.nx))goto 23133
      wtsum(kk) = 0.d0
      i=1
23134 if(.not.(i.le.nqd))goto 23136
      wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1)) * qdwt(i)
      wtsum(kk) = wtsum(kk) + wt(i,kk)
23135 i=i+1
      goto 23134
23136 continue
23132 kk=kk+1
      goto 23131
23133 continue
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23137 if(.not.(kk.le.nx))goto 23139
      i=1
23140 if(.not.(i.le.nxis))goto 23142
      mu(i) = ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1) / wtsum(kk)
23141 i=i+1
      goto 23140
23142 continue
      i=1
23143 if(.not.(i.le.nxis))goto 23145
      j=i
23146 if(.not.(j.le.nxis))goto 23148
      vwk(i,j) = 0.d0
      k=1
23149 if(.not.(k.le.nqd))goto 23151
      vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23150 k=k+1
      goto 23149
23151 continue
      vwk(i,j) = vwk(i,j) / wtsum(kk) - mu(i) * mu(j)
23147 j=j+1
      goto 23146
23148 continue
23144 i=i+1
      goto 23143
23145 continue
      call daxpy (nxis*nxis, xxwt(kk), vwk, 1, v, 1)
23138 kk=kk+1
      goto 23137
23139 continue
      i=1
23152 if(.not.(i.le.nxi))goto 23154
      j=i
23155 if(.not.(j.le.nxi))goto 23157
      v(i,j) = v(i,j) + q(i,j)
23156 j=j+1
      goto 23155
23157 continue
23153 i=i+1
      goto 23152
23154 continue
      i=1
23158 if(.not.(i.le.nxis))goto 23160
      jpvt(i) = 0
23159 i=i+1
      goto 23158
23160 continue
      call dchdc (v, nxis, nxis, vwk, jpvt, 1, rkv)
23161 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23161
      endif
23162 continue
      i=rkv+1
23163 if(.not.(i.le.nxis))goto 23165
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23164 i=i+1
      goto 23163
23165 continue
      return
      end
      subroutine llrmrkl (cd, nxis, qdrs, nqd, nx, xxwt, qdwt, wt0, offs
     *et, mchpr, wt, wtnew, mu, muwk, v, vwk, jpvt, cdnew, prec, maxiter
     *, info)
      integer nxis, nqd, nx, jpvt(*), maxiter, info
      double precision cd(*), qdrs(nqd,nxis,*), xxwt(*), qdwt(*), wt0(nq
     *d,*), offset(nqd,*), mchpr, wt(nqd,*), wtnew(nqd,*), mu(*), muwk(*
     *), v(nxis,*), vwk(nxis,*), cdnew(*), prec
      integer i, j, k, kk, iter, flag, idamax, infowk
      double precision ddot, dasum, rkl, tmp, mumax, rklnew, disc, disc0
      kk=1
23166 if(.not.(kk.le.nx))goto 23168
      i=1
23169 if(.not.(i.le.nqd))goto 23171
      wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1) + offset(i,
     *kk)) * qdwt(i)
23170 i=i+1
      goto 23169
23171 continue
      call dscal (nqd, 1.d0/dasum(nqd,wt(1,kk),1), wt(1,kk), 1)
23167 kk=kk+1
      goto 23166
23168 continue
      rkl = 0.d0
      kk=1
23172 if(.not.(kk.le.nx))goto 23174
      tmp = 0.d0
      i=1
23175 if(.not.(i.le.nqd))goto 23177
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23176 i=i+1
      goto 23175
23177 continue
      rkl = rkl + xxwt(kk) * tmp
23173 kk=kk+1
      goto 23172
23174 continue
      iter = 0
      flag = 0
23178 continue
      iter = iter + 1
      call dset (nxis, 0.d0, mu, 1)
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23181 if(.not.(kk.le.nx))goto 23183
      i=1
23184 if(.not.(i.le.nxis))goto 23186
      muwk(i) = ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1)
23185 i=i+1
      goto 23184
23186 continue
      i=1
23187 if(.not.(i.le.nxis))goto 23189
      j=i
23190 if(.not.(j.le.nxis))goto 23192
      vwk(i,j) = 0.d0
      k=1
23193 if(.not.(k.le.nqd))goto 23195
      vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23194 k=k+1
      goto 23193
23195 continue
      vwk(i,j) = vwk(i,j) - muwk(i) * muwk(j)
23191 j=j+1
      goto 23190
23192 continue
      muwk(i) = ddot (nqd, wt0(1,kk), 1, qdrs(1,i,kk), 1) - muwk(i)
23188 i=i+1
      goto 23187
23189 continue
      call daxpy (nxis, xxwt(kk), muwk, 1, mu, 1)
      call daxpy (nxis*nxis, xxwt(kk), vwk, 1, v, 1)
23182 kk=kk+1
      goto 23181
23183 continue
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23196 if(.not.(i.le.nxis))goto 23198
      jpvt(i) = 0
23197 i=i+1
      goto 23196
23198 continue
      call dmcdc (v, nxis, nxis, cdnew, jpvt, infowk)
23199 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      kk=1
23202 if(.not.(kk.le.nx))goto 23204
      i=1
23205 if(.not.(i.le.nqd))goto 23207
      wtnew(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cdnew, 1) + off
     *set(i,kk)) * qdwt(i)
23206 i=i+1
      goto 23205
23207 continue
      call dscal (nqd, 1.d0/dasum(nqd,wtnew(1,kk),1), wtnew(1,kk), 1)
23203 kk=kk+1
      goto 23202
23204 continue
      if((flag.eq.0).or.(flag.eq.2))then
      rklnew = 0.d0
      kk=1
23210 if(.not.(kk.le.nx))goto 23212
      tmp = 0.d0
      i=1
23213 if(.not.(i.le.nqd))goto 23215
      tmp = tmp + dlog(wt0(i,kk)/wtnew(i,kk)) * wt0(i,kk)
23214 i=i+1
      goto 23213
23215 continue
      rklnew = rklnew + xxwt(kk) * tmp
23211 kk=kk+1
      goto 23210
23212 continue
      endif
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      kk=1
23218 if(.not.(kk.le.nx))goto 23220
      i=1
23221 if(.not.(i.le.nqd))goto 23223
      wt(i,kk) = dexp (offset(i,kk)) * qdwt(i)
23222 i=i+1
      goto 23221
23223 continue
      call dscal (nqd, 1.d0/dasum(nqd,wt(1,kk),1), wt(1,kk), 1)
23219 kk=kk+1
      goto 23218
23220 continue
      rkl = 0.d0
      kk=1
23224 if(.not.(kk.le.nx))goto 23226
      tmp = 0.d0
      i=1
23227 if(.not.(i.le.nqd))goto 23229
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23228 i=i+1
      goto 23227
23229 continue
      rkl = rkl + xxwt(kk) * tmp
23225 kk=kk+1
      goto 23224
23226 continue
      iter = 0
      goto 23201
      endif
      if(flag.eq.3)then
      goto 23201
      endif
      if(rklnew-rkl.lt.1.d1*(1.d0+dabs(rkl))*mchpr)then
      goto 23201
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23201
      endif
23200 goto 23199
23201 continue
      if(flag.eq.1)then
      flag = 2
      goto 23179
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      kk=1
23240 if(.not.(kk.le.nx))goto 23242
      i=1
23243 if(.not.(i.le.nqd))goto 23245
      disc = dmax1 (disc, dabs(wt(i,kk)-wtnew(i,kk))/(1.d0+dabs(wt(i,kk)
     *)))
23244 i=i+1
      goto 23243
23245 continue
23241 kk=kk+1
      goto 23240
23242 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(rkl)))**2)
      disc0 = dmax1 ((mumax/(1.d0+rkl))**2, dabs(rkl-rklnew)/(1+dabs(rkl
     *)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nqd*nx, wtnew, 1, wt, 1)
      rkl = rklnew
      if(disc0.lt.prec)then
      goto 23180
      endif
      if(disc.lt.prec)then
      goto 23180
      endif
      if(iter.lt.maxiter)then
      goto 23179
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      kk=1
23254 if(.not.(kk.le.nx))goto 23256
      i=1
23257 if(.not.(i.le.nqd))goto 23259
      wt(i,kk) = dexp (offset(i,kk)) * qdwt(i)
23258 i=i+1
      goto 23257
23259 continue
      call dscal (nqd, 1.d0/dasum(nqd,wt(1,kk),1), wt(1,kk), 1)
23255 kk=kk+1
      goto 23254
23256 continue
      rkl = 0.d0
      kk=1
23260 if(.not.(kk.le.nx))goto 23262
      tmp = 0.d0
      i=1
23263 if(.not.(i.le.nqd))goto 23265
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23264 i=i+1
      goto 23263
23265 continue
      rkl = rkl + xxwt(kk) * tmp
23261 kk=kk+1
      goto 23260
23262 continue
      iter = 0
      flag = 2
      else
      info = 2
      goto 23180
      endif
23179 goto 23178
23180 continue
      rkl = 0.d0
      kk=1
23266 if(.not.(kk.le.nx))goto 23268
      tmp = 0.d0
      i=1
23269 if(.not.(i.le.nqd))goto 23271
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23270 i=i+1
      goto 23269
23271 continue
      rkl = rkl + xxwt(kk) * tmp
23267 kk=kk+1
      goto 23266
23268 continue
      wt(1,1) = rkl
      return
      end
