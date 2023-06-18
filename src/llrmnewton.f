C Output from Public domain Ratfor, version 1.01
      subroutine llrmnewton (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qd
     *rs, nqd, nx, xxwt, qdwt, prec, maxiter, mchpr, jpvt, wk, info)
      integer nxis, nxi, nobs, cntsum, nqd, nx, maxiter, jpvt(*), info
      double precision cd(*), q(nxi,*), rs(nxis,*), cnt(*), qdrs(nqd,nxi
     *s,*), xxwt(*), qdwt(*), prec, mchpr, wk(*)
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
      integer nxis, nxi, nobs, cntsum, nqd, nx, maxiter, jpvt(*), info
      double precision cd(*), q(nxi,*), rs(nxis,*), cnt(*), qdrs(nqd,nxi
     *s,*), xxwt(*), qdwt(*), prec, mchpr, wt(nqd,*), wtsum(*), mrs(*), 
     *fit(*), mu(*), muwk(*), v(nxis,*), vwk(nxis,*), cdnew(*), wtnew(nq
     *d,*), wtnewsum(*), fitnew(*), wk(*)
      integer i, j, k, kk, iter, flag, rkv, idamax, infowk
      double precision cnt1, norm, tmp, ddot, fitmean, lkhd, mumax, lkhd
     *new, disc, disc0, trc
      info = 0
      if(cntsum.ne.0)then
      cnt1 = 0.d0
      j=1
23002 if(.not.(j.le.nobs))goto 23004
      cnt1 = cnt1 + cnt(j)
23003 j=j+1
      goto 23002
23004 continue
      endif
      i=1
23005 if(.not.(i.le.nxis))goto 23007
      mrs(i) = 0.d0
      if(cntsum.eq.0)then
      j=1
23010 if(.not.(j.le.nobs))goto 23012
      mrs(i) = mrs(i) + rs(i,j)
23011 j=j+1
      goto 23010
23012 continue
      mrs(i) = mrs(i) / dble (nobs)
      else
      j=1
23013 if(.not.(j.le.nobs))goto 23015
      mrs(i) = mrs(i) + rs(i,j) * cnt(j)
23014 j=j+1
      goto 23013
23015 continue
      mrs(i) = mrs(i) / cnt1
      endif
23006 i=i+1
      goto 23005
23007 continue
      if(cntsum.eq.0)then
      trc = 1.d0 / dble (nobs)
      else
      trc = 1.d0 / cnt1
      endif
      norm = 0.d0
      kk=1
23018 if(.not.(kk.le.nx))goto 23020
      wtsum(kk) = 0.d0
      i=1
23021 if(.not.(i.le.nqd))goto 23023
      wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1)) * qdwt(i)
      wtsum(kk) = wtsum(kk) + wt(i,kk)
23022 i=i+1
      goto 23021
23023 continue
      norm = norm + xxwt(kk) * dlog (wtsum(kk))
23019 kk=kk+1
      goto 23018
23020 continue
      fitmean = 0.d0
      i=1
23024 if(.not.(i.le.nobs))goto 23026
      tmp = ddot (nxis, rs(1,i), 1, cd, 1)
      fit(i) = dexp (tmp)
      if(cntsum.ne.0)then
      tmp = tmp * cnt(i)
      endif
      fitmean = fitmean + tmp
23025 i=i+1
      goto 23024
23026 continue
      call dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, wk, 1)
      lkhd = ddot (nxi, cd, 1, wk, 1) / 2.d0 - fitmean * trc + norm
      iter = 0
      flag = 0
23029 continue
      iter = iter + 1
      call dset (nxis, 0.d0, mu, 1)
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23032 if(.not.(kk.le.nx))goto 23034
      i=1
23035 if(.not.(i.le.nxis))goto 23037
      muwk(i) = - ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1) / wtsum(kk)
23036 i=i+1
      goto 23035
23037 continue
      i=1
23038 if(.not.(i.le.nxis))goto 23040
      j=i
23041 if(.not.(j.le.nxis))goto 23043
      vwk(i,j) = 0.d0
      k=1
23044 if(.not.(k.le.nqd))goto 23046
      vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23045 k=k+1
      goto 23044
23046 continue
      vwk(i,j) = vwk(i,j) / wtsum(kk) - muwk(i) * muwk(j)
23042 j=j+1
      goto 23041
23043 continue
23039 i=i+1
      goto 23038
23040 continue
      call daxpy (nxis, xxwt(kk), muwk, 1, mu, 1)
      call daxpy (nxis*nxis, xxwt(kk), vwk, 1, v, 1)
23033 kk=kk+1
      goto 23032
23034 continue
      i=1
23047 if(.not.(i.le.nxi))goto 23049
      j=i
23050 if(.not.(j.le.nxi))goto 23052
      v(i,j) = v(i,j) + q(i,j)
23051 j=j+1
      goto 23050
23052 continue
23048 i=i+1
      goto 23047
23049 continue
      call daxpy (nxis, 1.d0, mrs, 1, mu, 1)
      call dsymv ('u', nxi, -1.d0, q, nxi, cd, 1, 1.d0, mu, 1)
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23053 if(.not.(i.le.nxis))goto 23055
      jpvt(i) = 0
23054 i=i+1
      goto 23053
23055 continue
      call dchdc (v, nxis, nxis, wk, jpvt, 1, rkv)
23056 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23056
      endif
23057 continue
      i=rkv+1
23058 if(.not.(i.le.nxis))goto 23060
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23059 i=i+1
      goto 23058
23060 continue
23061 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dset (nxis-rkv, 0.d0, cdnew(rkv+1), 1)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      norm = 0.d0
      kk=1
23064 if(.not.(kk.le.nx))goto 23066
      wtnewsum(kk) = 0.d0
      i=1
23067 if(.not.(i.le.nqd))goto 23069
      wtnew(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cdnew, 1)) * qd
     *wt(i)
      wtnewsum(kk) = wtnewsum(kk) + wtnew(i,kk)
23068 i=i+1
      goto 23067
23069 continue
      norm = norm + xxwt(kk) * dlog (wtnewsum(kk))
23065 kk=kk+1
      goto 23064
23066 continue
      if((flag.eq.0).or.(flag.eq.2))then
      fitmean = 0.d0
      i=1
23072 if(.not.(i.le.nobs))goto 23074
      tmp = ddot (nxis, rs(1,i), 1, cdnew, 1)
      if(tmp.gt.3.d2)then
      flag = flag + 1
      goto 23074
      endif
      fitnew(i) = dexp (tmp)
      if(cntsum.ne.0)then
      tmp = tmp * cnt(i)
      endif
      fitmean = fitmean + tmp
23073 i=i+1
      goto 23072
23074 continue
      call dsymv ('u', nxi, 1.d0, q, nxi, cdnew, 1, 0.d0, wk, 1)
      lkhdnew = ddot (nxi, cdnew, 1, wk, 1) / 2.d0 - fitmean * trc + nor
     *m
      endif
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      kk=1
23081 if(.not.(kk.le.nx))goto 23083
      call dcopy (nqd, qdwt, 1, wt(1,kk), 1)
23082 kk=kk+1
      goto 23081
23083 continue
      tmp = 0.d0
      i=1
23084 if(.not.(i.le.nqd))goto 23086
      tmp = tmp + qdwt(i)
23085 i=i+1
      goto 23084
23086 continue
      call dset (nx, tmp, wtsum, 1)
      call dset (nobs, 1.d0, fit, 1)
      lkhd = 0.d0
      iter = 0
      goto 23063
      endif
      if(flag.eq.3)then
      goto 23063
      endif
      if(lkhdnew-lkhd.lt.1.d1*(1.d0+dabs(lkhd))*mchpr)then
      goto 23063
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23063
      endif
23062 goto 23061
23063 continue
      if(flag.eq.1)then
      flag = 2
      goto 23030
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      kk=1
23097 if(.not.(kk.le.nx))goto 23099
      i=1
23100 if(.not.(i.le.nqd))goto 23102
      disc = dmax1 (disc, dabs(wt(i,kk)-wtnew(i,kk))/(1.d0+dabs(wt(i,kk)
     *)))
23101 i=i+1
      goto 23100
23102 continue
23098 kk=kk+1
      goto 23097
23099 continue
      i=1
23103 if(.not.(i.le.nobs))goto 23105
      disc = dmax1 (disc, dabs(fit(i)-fitnew(i))/(1.d0+dabs(fit(i))))
23104 i=i+1
      goto 23103
23105 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
      disc0 = dmax1 ((mumax/(1.d0+lkhd))**2, dabs(lkhd-lkhdnew)/(1+dabs(
     *lkhd)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nqd*nx, wtnew, 1, wt, 1)
      call dcopy (nx, wtnewsum, 1, wtsum, 1)
      call dcopy (nobs, fitnew, 1, fit, 1)
      lkhd = lkhdnew
      if(disc0.lt.prec)then
      goto 23031
      endif
      if(disc.lt.prec)then
      goto 23031
      endif
      if(iter.lt.maxiter)then
      goto 23030
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      kk=1
23114 if(.not.(kk.le.nx))goto 23116
      call dcopy (nqd, qdwt, 1, wt(1,kk), 1)
23115 kk=kk+1
      goto 23114
23116 continue
      tmp = 0.d0
      i=1
23117 if(.not.(i.le.nqd))goto 23119
      tmp = tmp + qdwt(i)
23118 i=i+1
      goto 23117
23119 continue
      call dset (nx, tmp, wtsum, 1)
      call dset (nobs, 1.d0, fit, 1)
      lkhd = 0.d0
      iter = 0
      flag = 2
      else
      info = 2
      goto 23031
      endif
23030 goto 23029
23031 continue
      i=1
23120 if(.not.(i.le.nobs))goto 23122
      call daxpy (nxis, -1.d0, mrs, 1, rs(1,i), 1)
      call dprmut (rs(1,i), nxis, jpvt, 0)
      if(cntsum.ne.0)then
      call dscal (nxis, dsqrt(cnt(i)), rs(1,i), 1)
      endif
      call dtrsl (v, nxis, nxis, rs(1,i), 11, infowk)
23121 i=i+1
      goto 23120
23122 continue
      trc = ddot (nobs*nxis, rs, 1, rs, 1)
      if(cntsum.eq.0)then
      trc = trc / dble(nobs) / (dble(nobs)-1.d0)
      lkhd = 0.d0
      i=1
23127 if(.not.(i.le.nobs))goto 23129
      lkhd = lkhd + dlog (fit(i))
23128 i=i+1
      goto 23127
23129 continue
      lkhd = lkhd / dble (nobs)
      else
      trc = trc / cnt1 / (cnt1-1.d0)
      lkhd = 0.d0
      i=1
23130 if(.not.(i.le.nobs))goto 23132
      lkhd = lkhd + cnt(i) * dlog (fit(i))
23131 i=i+1
      goto 23130
23132 continue
      lkhd = lkhd / cnt1
      endif
      kk=1
23133 if(.not.(kk.le.nx))goto 23135
      lkhd = lkhd - xxwt(kk) * dlog (wtsum(kk))
23134 kk=kk+1
      goto 23133
23135 continue
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
23136 if(.not.(kk.le.nx))goto 23138
      wtsum(kk) = 0.d0
      i=1
23139 if(.not.(i.le.nqd))goto 23141
      wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1)) * qdwt(i)
      wtsum(kk) = wtsum(kk) + wt(i,kk)
23140 i=i+1
      goto 23139
23141 continue
23137 kk=kk+1
      goto 23136
23138 continue
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23142 if(.not.(kk.le.nx))goto 23144
      i=1
23145 if(.not.(i.le.nxis))goto 23147
      mu(i) = ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1) / wtsum(kk)
23146 i=i+1
      goto 23145
23147 continue
      i=1
23148 if(.not.(i.le.nxis))goto 23150
      j=i
23151 if(.not.(j.le.nxis))goto 23153
      vwk(i,j) = 0.d0
      k=1
23154 if(.not.(k.le.nqd))goto 23156
      vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23155 k=k+1
      goto 23154
23156 continue
      vwk(i,j) = vwk(i,j) / wtsum(kk) - mu(i) * mu(j)
23152 j=j+1
      goto 23151
23153 continue
23149 i=i+1
      goto 23148
23150 continue
      call daxpy (nxis*nxis, xxwt(kk), vwk, 1, v, 1)
23143 kk=kk+1
      goto 23142
23144 continue
      i=1
23157 if(.not.(i.le.nxi))goto 23159
      j=i
23160 if(.not.(j.le.nxi))goto 23162
      v(i,j) = v(i,j) + q(i,j)
23161 j=j+1
      goto 23160
23162 continue
23158 i=i+1
      goto 23157
23159 continue
      i=1
23163 if(.not.(i.le.nxis))goto 23165
      jpvt(i) = 0
23164 i=i+1
      goto 23163
23165 continue
      call dchdc (v, nxis, nxis, vwk, jpvt, 1, rkv)
23166 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23166
      endif
23167 continue
      i=rkv+1
23168 if(.not.(i.le.nxis))goto 23170
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23169 i=i+1
      goto 23168
23170 continue
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
23171 if(.not.(kk.le.nx))goto 23173
      i=1
23174 if(.not.(i.le.nqd))goto 23176
      wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1) + offset(i,
     *kk)) * qdwt(i)
23175 i=i+1
      goto 23174
23176 continue
      call dscal (nqd, 1.d0/dasum(nqd,wt(1,kk),1), wt(1,kk), 1)
23172 kk=kk+1
      goto 23171
23173 continue
      rkl = 0.d0
      kk=1
23177 if(.not.(kk.le.nx))goto 23179
      tmp = 0.d0
      i=1
23180 if(.not.(i.le.nqd))goto 23182
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23181 i=i+1
      goto 23180
23182 continue
      rkl = rkl + xxwt(kk) * tmp
23178 kk=kk+1
      goto 23177
23179 continue
      iter = 0
      flag = 0
23183 continue
      iter = iter + 1
      call dset (nxis, 0.d0, mu, 1)
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23186 if(.not.(kk.le.nx))goto 23188
      i=1
23189 if(.not.(i.le.nxis))goto 23191
      muwk(i) = ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1)
23190 i=i+1
      goto 23189
23191 continue
      i=1
23192 if(.not.(i.le.nxis))goto 23194
      j=i
23195 if(.not.(j.le.nxis))goto 23197
      vwk(i,j) = 0.d0
      k=1
23198 if(.not.(k.le.nqd))goto 23200
      vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23199 k=k+1
      goto 23198
23200 continue
      vwk(i,j) = vwk(i,j) - muwk(i) * muwk(j)
23196 j=j+1
      goto 23195
23197 continue
      muwk(i) = ddot (nqd, wt0(1,kk), 1, qdrs(1,i,kk), 1) - muwk(i)
23193 i=i+1
      goto 23192
23194 continue
      call daxpy (nxis, xxwt(kk), muwk, 1, mu, 1)
      call daxpy (nxis*nxis, xxwt(kk), vwk, 1, v, 1)
23187 kk=kk+1
      goto 23186
23188 continue
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23201 if(.not.(i.le.nxis))goto 23203
      jpvt(i) = 0
23202 i=i+1
      goto 23201
23203 continue
      call dmcdc (v, nxis, nxis, cdnew, jpvt, infowk)
23204 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      kk=1
23207 if(.not.(kk.le.nx))goto 23209
      i=1
23210 if(.not.(i.le.nqd))goto 23212
      wtnew(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cdnew, 1) + off
     *set(i,kk)) * qdwt(i)
23211 i=i+1
      goto 23210
23212 continue
      call dscal (nqd, 1.d0/dasum(nqd,wtnew(1,kk),1), wtnew(1,kk), 1)
23208 kk=kk+1
      goto 23207
23209 continue
      if((flag.eq.0).or.(flag.eq.2))then
      rklnew = 0.d0
      kk=1
23215 if(.not.(kk.le.nx))goto 23217
      tmp = 0.d0
      i=1
23218 if(.not.(i.le.nqd))goto 23220
      tmp = tmp + dlog(wt0(i,kk)/wtnew(i,kk)) * wt0(i,kk)
23219 i=i+1
      goto 23218
23220 continue
      rklnew = rklnew + xxwt(kk) * tmp
23216 kk=kk+1
      goto 23215
23217 continue
      endif
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      kk=1
23223 if(.not.(kk.le.nx))goto 23225
      i=1
23226 if(.not.(i.le.nqd))goto 23228
      wt(i,kk) = dexp (offset(i,kk)) * qdwt(i)
23227 i=i+1
      goto 23226
23228 continue
      call dscal (nqd, 1.d0/dasum(nqd,wt(1,kk),1), wt(1,kk), 1)
23224 kk=kk+1
      goto 23223
23225 continue
      rkl = 0.d0
      kk=1
23229 if(.not.(kk.le.nx))goto 23231
      tmp = 0.d0
      i=1
23232 if(.not.(i.le.nqd))goto 23234
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23233 i=i+1
      goto 23232
23234 continue
      rkl = rkl + xxwt(kk) * tmp
23230 kk=kk+1
      goto 23229
23231 continue
      iter = 0
      goto 23206
      endif
      if(flag.eq.3)then
      goto 23206
      endif
      if(rklnew-rkl.lt.1.d1*(1.d0+dabs(rkl))*mchpr)then
      goto 23206
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23206
      endif
23205 goto 23204
23206 continue
      if(flag.eq.1)then
      flag = 2
      goto 23184
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      kk=1
23245 if(.not.(kk.le.nx))goto 23247
      i=1
23248 if(.not.(i.le.nqd))goto 23250
      disc = dmax1 (disc, dabs(wt(i,kk)-wtnew(i,kk))/(1.d0+dabs(wt(i,kk)
     *)))
23249 i=i+1
      goto 23248
23250 continue
23246 kk=kk+1
      goto 23245
23247 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(rkl)))**2)
      disc0 = dmax1 ((mumax/(1.d0+rkl))**2, dabs(rkl-rklnew)/(1+dabs(rkl
     *)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nqd*nx, wtnew, 1, wt, 1)
      rkl = rklnew
      if(disc0.lt.prec)then
      goto 23185
      endif
      if(disc.lt.prec)then
      goto 23185
      endif
      if(iter.lt.maxiter)then
      goto 23184
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      kk=1
23259 if(.not.(kk.le.nx))goto 23261
      i=1
23262 if(.not.(i.le.nqd))goto 23264
      wt(i,kk) = dexp (offset(i,kk)) * qdwt(i)
23263 i=i+1
      goto 23262
23264 continue
      call dscal (nqd, 1.d0/dasum(nqd,wt(1,kk),1), wt(1,kk), 1)
23260 kk=kk+1
      goto 23259
23261 continue
      rkl = 0.d0
      kk=1
23265 if(.not.(kk.le.nx))goto 23267
      tmp = 0.d0
      i=1
23268 if(.not.(i.le.nqd))goto 23270
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23269 i=i+1
      goto 23268
23270 continue
      rkl = rkl + xxwt(kk) * tmp
23266 kk=kk+1
      goto 23265
23267 continue
      iter = 0
      flag = 2
      else
      info = 2
      goto 23185
      endif
23184 goto 23183
23185 continue
      rkl = 0.d0
      kk=1
23271 if(.not.(kk.le.nx))goto 23273
      tmp = 0.d0
      i=1
23274 if(.not.(i.le.nqd))goto 23276
      tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
23275 i=i+1
      goto 23274
23276 continue
      rkl = rkl + xxwt(kk) * tmp
23272 kk=kk+1
      goto 23271
23273 continue
      wt(1,1) = rkl
      return
      end
