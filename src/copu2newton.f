C Output from Public domain Ratfor, version 1.01
      subroutine copu2newton (cd, nxis, q, nxi, rs0, n0, cntsum0, cnt0, 
     *qdrs, nqd, qdrs1, wt1, n1, cntsum1, cnt1, qdrs2, wt2, n2, cntsum2,
     * cnt2, wt3, n3, cntsum3, cnt3, nt, twt, qdwt, tind, prec, maxiter,
     * mchpr, jpvt, wk, info)
      integer nxis, nxi, n0, cntsum0, cnt0(*), nqd, n1, cntsum1, cnt1(*)
     *, n2, cntsum2, cnt2(*), n3, cntsum3, cnt3(*), nt, tind(*), maxiter
     *, jpvt(*), info
      double precision cd(*), q(nxi,*), rs0(nxis,*), qdrs(nxis,*), qdrs1
     *(nqd,nxis,*), wt1(nqd,*), qdrs2(nqd,nxis,*), wt2(nqd,*), wt3(nqd,2
     *,*), twt(*), qdwt(nqd,2,*), prec, mchpr, wk(*)
      integer imrs, imrs2, ieta, ieta1, ieta2, imu, imuwk, iv, ivwk, icd
     *new, imut, iwk
      imrs = 1
      imrs2 = imrs + nxis
      ieta = imrs2 + nxis
      ieta1 = ieta + nqd*nqd
      ieta2 = ieta1 + nqd*n1
      imu = ieta2 + nqd*n2
      imuwk = imu + nxis
      iv = imuwk + nxis
      ivwk = iv + nxis*nxis
      icdnew = ivwk + nxis*nxis
      imut = icdnew + nxis
      iwk = imut + nxis*nt
      call copu2newton1 (cd, nxis, q, nxi, rs0, n0, cntsum0, cnt0, qdrs,
     * nqd, qdrs1, wt1, n1, cntsum1, cnt1, qdrs2, wt2, n2, cntsum2, cnt2
     *, wt3, n3, cntsum3, cnt3, nt, twt, qdwt, tind, prec, maxiter, mchp
     *r, wk(imrs), wk(imrs2), wk(ieta), wk(ieta1), wk(ieta2), wk(imu), w
     *k(imuwk), wk(iv), wk(ivwk), jpvt, wk(icdnew), wk(imut), wk(iwk), i
     *nfo)
      return
      end
      subroutine copu2newton1 (cd, nxis, q, nxi, rs0, n0, cntsum0, cnt0,
     * qdrs, nqd, qdrs1, wt1, n1, cntsum1, cnt1, qdrs2, wt2, n2, cntsum2
     *, cnt2, wt3, n3, cntsum3, cnt3, nt, twt, qdwt, tind, prec, maxiter
     *, mchpr, mrs, mrs2, eta, eta1, eta2, mu, muwk, v, vwk, jpvt, cdnew
     *, mut, wk, info)
      integer nxis, nxi, n0, cntsum0, cnt0(*), nqd, n1, cntsum1, cnt1(*)
     *, n2, cntsum2, cnt2(*), n3, cntsum3, cnt3(*), nt, tind(*), maxiter
     *, jpvt(*), info
      double precision cd(*), q(nxi,*), rs0(nxis,*), qdrs(nxis,*), qdrs1
     *(nqd,nxis,*), wt1(nqd,*), qdrs2(nqd,nxis,*), wt2(nqd,*), wt3(nqd,2
     *,*), twt(*), qdwt(nqd,2,*), prec, mchpr, mrs(*), mrs2(*), eta(*), 
     *eta1(nqd,*), eta2(nqd,*), mu(*), muwk(*), v(nxis,*), vwk(nxis,*), 
     *cdnew(*), mut(nxis,*), wk(*)
      integer nobs, i, j, k, kk, iter, flag, rkv, idamax, infowk
      double precision tmp, ddot, dasum, lkhd, mumax, lkhdnew, disc, dis
     *c0, trc
      info = 0
      if(cntsum0.eq.0)then
      nobs = n0
      else
      nobs = cntsum0
      endif
      if(cntsum1.eq.0)then
      nobs = nobs + n1
      else
      nobs = nobs + cntsum1
      endif
      if(cntsum2.eq.0)then
      nobs = nobs + n2
      else
      nobs = nobs + cntsum2
      endif
      if(cntsum3.eq.0)then
      nobs = nobs + n3
      else
      nobs = nobs + cntsum3
      endif
      i=1
23008 if(.not.(i.le.nxis))goto 23010
      mrs(i) = 0.d0
      j=1
23011 if(.not.(j.le.n0))goto 23013
      if(cntsum0.eq.0)then
      mrs(i) = mrs(i) + rs0(i,j)
      else
      mrs(i) = mrs(i) + rs0(i,j) * dfloat (cnt0(j))
      endif
23012 j=j+1
      goto 23011
23013 continue
23009 i=i+1
      goto 23008
23010 continue
      i=1
23016 if(.not.(i.le.nqd*nqd))goto 23018
      eta(i) = dexp (ddot (nxis, qdrs(1,i), 1, cd, 1))
23017 i=i+1
      goto 23016
23018 continue
      lkhd = 0.d0
      i=1
23019 if(.not.(i.le.n0))goto 23021
      tmp = ddot (nxis, rs0(1,i), 1, cd, 1)
      if(cntsum0.ne.0)then
      tmp = tmp * dfloat (cnt0(i))
      endif
      lkhd = lkhd - tmp
23020 i=i+1
      goto 23019
23021 continue
      i=1
23024 if(.not.(i.le.n1))goto 23026
      j=1
23027 if(.not.(j.le.nqd))goto 23029
      eta1(j,i) = dexp (ddot (nxis, qdrs1(j,1,i), nqd, cd, 1)) * wt1(j,i
     *)
23028 j=j+1
      goto 23027
23029 continue
      if(cntsum1.eq.0)then
      lkhd = lkhd - dlog (dasum (nqd, eta1(1,i), 1))
      else
      lkhd = lkhd - dlog (dasum (nqd, eta1(1,i), 1)) * dfloat (cnt1(i))
      endif
23025 i=i+1
      goto 23024
23026 continue
      i=1
23032 if(.not.(i.le.n2))goto 23034
      j=1
23035 if(.not.(j.le.nqd))goto 23037
      eta2(j,i) = dexp (ddot (nxis, qdrs2(j,1,i), nqd, cd, 1)) * wt2(j,i
     *)
23036 j=j+1
      goto 23035
23037 continue
      if(cntsum2.eq.0)then
      lkhd = lkhd - dlog (dasum (nqd, eta2(1,i), 1))
      else
      lkhd = lkhd - dlog (dasum (nqd, eta2(1,i), 1)) * dfloat (cnt2(i))
      endif
23033 i=i+1
      goto 23032
23034 continue
      i=1
23040 if(.not.(i.le.n3))goto 23042
      tmp = 0.d0
      j=1
23043 if(.not.(j.le.nqd))goto 23045
      tmp = tmp + ddot (nqd, eta((j-1)*nqd+1), 1, wt3(1,1,i), 1) * wt3(j
     *,2,i)
23044 j=j+1
      goto 23043
23045 continue
      if(cntsum3.eq.0)then
      lkhd = lkhd - dlog (tmp)
      else
      lkhd = lkhd - dlog (tmp) * dfloat (cnt3(i))
      endif
23041 i=i+1
      goto 23040
23042 continue
      lkhd = lkhd / dfloat (nobs)
      i=1
23048 if(.not.(i.le.nt))goto 23050
      tmp = 0.d0
      j=1
23051 if(.not.(j.le.nqd))goto 23053
      tmp = tmp + ddot (nqd, eta((j-1)*nqd+1), 1, qdwt(1,1,i), 1) * qdwt
     *(j,2,i)
23052 j=j+1
      goto 23051
23053 continue
      lkhd = lkhd + dlog (tmp) * twt(i)
23049 i=i+1
      goto 23048
23050 continue
      call dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, muwk, 1)
      lkhd = lkhd + ddot (nxi, cd, 1, muwk, 1) / 2.d0
      iter = 0
      flag = 0
23054 continue
      iter = iter + 1
      call dcopy (nxis, mrs, 1, mu, 1)
      call dset(nxis*nxis, 0.d0, v, 1)
      i=1
23057 if(.not.(i.le.n1))goto 23059
      tmp = dasum (nqd, eta1(1,i), 1)
      j=1
23060 if(.not.(j.le.nxis))goto 23062
      muwk(j) = ddot (nqd, eta1(1,i), 1, qdrs1(1,j,i), 1) / tmp
23061 j=j+1
      goto 23060
23062 continue
      j=1
23063 if(.not.(j.le.nxis))goto 23065
      k=j
23066 if(.not.(k.le.nxis))goto 23068
      vwk(j,k) = 0.d0
      kk=1
23069 if(.not.(kk.le.nqd))goto 23071
      vwk(j,k) = vwk(j,k) + eta1(kk,i)*qdrs1(kk,j,i)*qdrs1(kk,k,i)
23070 kk=kk+1
      goto 23069
23071 continue
      vwk(j,k) = vwk(j,k) / tmp - muwk(j) * muwk(k)
23067 k=k+1
      goto 23066
23068 continue
23064 j=j+1
      goto 23063
23065 continue
      if(cntsum1.eq.0)then
      call daxpy (nxis, 1.d0, muwk, 1, mu, 1)
      call daxpy (nxis*nxis, -1.d0, vwk, 1, v, 1)
      else
      call daxpy (nxis, dfloat (cnt1(i)), muwk, 1, mu, 1)
      call daxpy (nxis*nxis, -dfloat (cnt1(i)), vwk, 1, v, 1)
      endif
23058 i=i+1
      goto 23057
23059 continue
      i=1
23074 if(.not.(i.le.n2))goto 23076
      tmp = dasum (nqd, eta2(1,i), 1)
      j=1
23077 if(.not.(j.le.nxis))goto 23079
      muwk(j) = ddot (nqd, eta2(1,i), 1, qdrs2(1,j,i), 1) / tmp
23078 j=j+1
      goto 23077
23079 continue
      j=1
23080 if(.not.(j.le.nxis))goto 23082
      k=j
23083 if(.not.(k.le.nxis))goto 23085
      vwk(j,k) = 0.d0
      kk=1
23086 if(.not.(kk.le.nqd))goto 23088
      vwk(j,k) = vwk(j,k) + eta2(kk,i)*qdrs2(kk,j,i)*qdrs2(kk,k,i)
23087 kk=kk+1
      goto 23086
23088 continue
      vwk(j,k) = vwk(j,k) / tmp - muwk(j) * muwk(k)
23084 k=k+1
      goto 23083
23085 continue
23081 j=j+1
      goto 23080
23082 continue
      if(cntsum2.eq.0)then
      call daxpy (nxis, 1.d0, muwk, 1, mu, 1)
      call daxpy (nxis*nxis, -1.d0, vwk, 1, v, 1)
      else
      call daxpy (nxis, dfloat (cnt2(i)), muwk, 1, mu, 1)
      call daxpy (nxis*nxis, -dfloat (cnt2(i)), vwk, 1, v, 1)
      endif
23075 i=i+1
      goto 23074
23076 continue
      i=1
23091 if(.not.(i.le.n3))goto 23093
      j=1
23094 if(.not.(j.le.nqd))goto 23096
      k=1
23097 if(.not.(k.le.nqd))goto 23099
      wk((k-1)*nqd+j) = eta((k-1)*nqd+j) * wt3(j,1,i) * wt3(k,2,i)
23098 k=k+1
      goto 23097
23099 continue
23095 j=j+1
      goto 23094
23096 continue
      tmp = dasum (nqd*nqd, wk, 1)
      j=1
23100 if(.not.(j.le.nxis))goto 23102
      muwk(j) = ddot (nqd*nqd, wk, 1, qdrs(j,1), nxis) / tmp
23101 j=j+1
      goto 23100
23102 continue
      j=1
23103 if(.not.(j.le.nxis))goto 23105
      k=j
23106 if(.not.(k.le.nxis))goto 23108
      vwk(j,k) = 0.d0
      kk=1
23109 if(.not.(kk.le.nqd*nqd))goto 23111
      vwk(j,k) = vwk(j,k) + wk(kk)*qdrs(j,kk)*qdrs(k,kk)
23110 kk=kk+1
      goto 23109
23111 continue
      vwk(j,k) = vwk(j,k) / tmp - muwk(j) * muwk(k)
23107 k=k+1
      goto 23106
23108 continue
23104 j=j+1
      goto 23103
23105 continue
      if(cntsum3.eq.0)then
      call daxpy (nxis, 1.d0, muwk, 1, mu, 1)
      call daxpy (nxis*nxis, -1.d0, vwk, 1, v, 1)
      else
      call daxpy (nxis, dfloat (cnt3(i)), muwk, 1, mu, 1)
      call daxpy (nxis*nxis, -dfloat (cnt3(i)), vwk, 1, v, 1)
      endif
23092 i=i+1
      goto 23091
23093 continue
      call dscal (nxis, 1.d0/dfloat(nobs), mu, 1)
      call dscal (nxis*nxis, 1.d0/dfloat(nobs), v, 1)
      i=1
23114 if(.not.(i.le.nt))goto 23116
      j=1
23117 if(.not.(j.le.nqd))goto 23119
      k=1
23120 if(.not.(k.le.nqd))goto 23122
      wk((k-1)*nqd+j) = eta((k-1)*nqd+j) * qdwt(j,1,i) * qdwt(k,2,i)
23121 k=k+1
      goto 23120
23122 continue
23118 j=j+1
      goto 23117
23119 continue
      tmp = dasum (nqd*nqd, wk, 1)
      j=1
23123 if(.not.(j.le.nxis))goto 23125
      muwk(j) = ddot (nqd*nqd, wk, 1, qdrs(j,1), nxis) / tmp
23124 j=j+1
      goto 23123
23125 continue
      j=1
23126 if(.not.(j.le.nxis))goto 23128
      k=j
23129 if(.not.(k.le.nxis))goto 23131
      vwk(j,k) = 0.d0
      kk=1
23132 if(.not.(kk.le.nqd*nqd))goto 23134
      vwk(j,k) = vwk(j,k) + wk(kk)*qdrs(j,kk)*qdrs(k,kk)
23133 kk=kk+1
      goto 23132
23134 continue
      vwk(j,k) = vwk(j,k) / tmp - muwk(j) * muwk(k)
23130 k=k+1
      goto 23129
23131 continue
23127 j=j+1
      goto 23126
23128 continue
      call dcopy (nxis, muwk, 1, mut(1,i), 1)
      call daxpy (nxis, -twt(i), muwk, 1, mu, 1)
      call daxpy (nxis*nxis, twt(i), vwk, 1, v, 1)
23115 i=i+1
      goto 23114
23116 continue
      call dcopy (nxis, mu, 1, mrs2, 1)
      call dsymv ('u', nxi, -1.d0, q, nxi, cd, 1, 1.d0, mu, 1)
      i=1
23135 if(.not.(i.le.nxi))goto 23137
      j=i
23138 if(.not.(j.le.nxi))goto 23140
      v(i,j) = v(i,j) + q(i,j)
23139 j=j+1
      goto 23138
23140 continue
23136 i=i+1
      goto 23135
23137 continue
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23141 if(.not.(i.le.nxis))goto 23143
      jpvt(i) = 0
23142 i=i+1
      goto 23141
23143 continue
      call dchdc (v, nxis, nxis, muwk, jpvt, 1, rkv)
23144 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23144
      endif
23145 continue
      i=rkv+1
23146 if(.not.(i.le.nxis))goto 23148
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23147 i=i+1
      goto 23146
23148 continue
23149 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dset (nxis-rkv, 0.d0, cdnew(rkv+1), 1)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      i=1
23152 if(.not.(i.le.nqd*nqd))goto 23154
      tmp = ddot (nxis, qdrs(1,i), 1, cdnew, 1)
      if(tmp.gt.3.d2)then
      flag = flag + 1
      goto 23154
      endif
      wk(i) = dexp (tmp)
23153 i=i+1
      goto 23152
23154 continue
      lkhdnew = 0.d0
      i=1
23157 if(.not.(i.le.n0))goto 23159
      tmp = ddot (nxis, rs0(1,i), 1, cdnew, 1)
      if(cntsum0.ne.0)then
      tmp = tmp * dfloat (cnt0(i))
      endif
      lkhdnew = lkhdnew - tmp
23158 i=i+1
      goto 23157
23159 continue
      i=1
23162 if(.not.(i.le.n1))goto 23164
      j=1
23165 if(.not.(j.le.nqd))goto 23167
      eta1(j,i) = dexp (ddot (nxis, qdrs1(j,1,i), nqd, cdnew, 1)) * wt1(
     *j,i)
23166 j=j+1
      goto 23165
23167 continue
      if(cntsum1.eq.0)then
      lkhdnew = lkhdnew - dlog (dasum (nqd, eta1(1,i), 1))
      else
      lkhdnew = lkhdnew - dlog (dasum (nqd, eta1(1,i), 1)) * dfloat (cnt
     *1(i))
      endif
23163 i=i+1
      goto 23162
23164 continue
      i=1
23170 if(.not.(i.le.n2))goto 23172
      j=1
23173 if(.not.(j.le.nqd))goto 23175
      eta2(j,i) = dexp (ddot (nxis, qdrs2(j,1,i), nqd, cdnew, 1)) * wt2(
     *j,i)
23174 j=j+1
      goto 23173
23175 continue
      if(cntsum2.eq.0)then
      lkhdnew = lkhdnew - dlog (dasum (nqd, eta2(1,i), 1))
      else
      lkhdnew = lkhdnew - dlog (dasum (nqd, eta2(1,i), 1)) * dfloat (cnt
     *2(i))
      endif
23171 i=i+1
      goto 23170
23172 continue
      i=1
23178 if(.not.(i.le.n3))goto 23180
      tmp = 0.d0
      j=1
23181 if(.not.(j.le.nqd))goto 23183
      tmp = tmp + ddot (nqd, wk((j-1)*nqd+1), 1, wt3(1,1,i), 1) * wt3(j,
     *2,i)
23182 j=j+1
      goto 23181
23183 continue
      if(cntsum3.eq.0)then
      lkhdnew = lkhdnew - dlog (tmp)
      else
      lkhdnew = lkhdnew - dlog (tmp) * dfloat (cnt3(i))
      endif
23179 i=i+1
      goto 23178
23180 continue
      lkhdnew = lkhdnew / dfloat (nobs)
      i=1
23186 if(.not.(i.le.nt))goto 23188
      tmp = 0.d0
      j=1
23189 if(.not.(j.le.nqd))goto 23191
      tmp = tmp + ddot (nqd, wk((j-1)*nqd+1), 1, qdwt(1,1,i), 1) * qdwt(
     *j,2,i)
23190 j=j+1
      goto 23189
23191 continue
      lkhdnew = lkhdnew + dlog (tmp) * twt(i)
23187 i=i+1
      goto 23186
23188 continue
      call dsymv ('u', nxi, 1.d0, q, nxi, cdnew, 1, 0.d0, muwk, 1)
      lkhdnew = lkhdnew + ddot (nxi, cdnew, 1, muwk, 1) / 2.d0
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      call dset (nqd*nqd, 1.d0, eta, 1)
      call dcopy (nqd*n1, wt1, 1, eta1, 1)
      call dcopy (nqd*n2, wt2, 1, eta2, 1)
      lkhd = 0.d0
      i=1
23194 if(.not.(i.le.n1))goto 23196
      if(cntsum1.eq.0)then
      lkhd = lkhd - dlog (dasum (nqd, eta1(1,i), 1))
      else
      lkhd = lkhd - dlog (dasum (nqd, eta1(1,i), 1)) * dfloat (cnt1(i))
      endif
23195 i=i+1
      goto 23194
23196 continue
      i=1
23199 if(.not.(i.le.n2))goto 23201
      if(cntsum2.eq.0)then
      lkhd = lkhd - dlog (dasum (nqd, eta2(1,i), 1))
      else
      lkhd = lkhd - dlog (dasum (nqd, eta2(1,i), 1)) * dfloat (cnt2(i))
      endif
23200 i=i+1
      goto 23199
23201 continue
      i=1
23204 if(.not.(i.le.n3))goto 23206
      tmp = dasum (nqd, wt3(1,1,i), 1) * dasum (nqd, wt3(1,2,i), 1)
      if(cntsum3.eq.0)then
      lkhdnew = lkhdnew - dlog (tmp)
      else
      lkhdnew = lkhdnew - dlog (tmp) * dfloat (cnt3(i))
      endif
23205 i=i+1
      goto 23204
23206 continue
      lkhd = lkhd / dfloat (nobs)
      i=1
23209 if(.not.(i.le.nt))goto 23211
      tmp = dasum (nqd, qdwt(1,1,i), 1) * dasum (nqd, qdwt(1,2,i), 1)
      lkhd = lkhd + dlog (tmp) * twt(i)
23210 i=i+1
      goto 23209
23211 continue
      iter = 0
      goto 23151
      endif
      if(flag.eq.3)then
      goto 23151
      endif
      if(lkhdnew-lkhd.lt.1.d1*(1.d0+dabs(lkhd))*mchpr)then
      goto 23151
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23151
      endif
23150 goto 23149
23151 continue
      if(flag.eq.1)then
      flag = 2
      goto 23055
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      i=1
23222 if(.not.(i.le.nqd*nqd))goto 23224
      disc = dmax1 (disc, dabs(eta(i)-wk(i))/(1.d0+dabs(eta(i))))
23223 i=i+1
      goto 23222
23224 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
      disc0 = dmax1 ((mumax/(1.d0+lkhd))**2, dabs(lkhd-lkhdnew)/(1.d0+da
     *bs(lkhd)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nqd*nqd, wk, 1, eta, 1)
      lkhd = lkhdnew
      if(disc0.lt.prec)then
      goto 23056
      endif
      if(disc.lt.prec)then
      goto 23056
      endif
      if(iter.lt.maxiter)then
      goto 23055
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      call dset (nqd*nqd, 1.d0, eta, 1)
      call dcopy (nqd*n1, wt1, 1, eta1, 1)
      call dcopy (nqd*n2, wt2, 1, eta2, 1)
      lkhd = 0.d0
      i=1
23233 if(.not.(i.le.n1))goto 23235
      if(cntsum1.eq.0)then
      lkhd = lkhd - dlog (dasum (nqd, eta1(1,i), 1))
      else
      lkhd = lkhd - dlog (dasum (nqd, eta1(1,i), 1)) * dfloat (cnt1(i))
      endif
23234 i=i+1
      goto 23233
23235 continue
      i=1
23238 if(.not.(i.le.n2))goto 23240
      if(cntsum2.eq.0)then
      lkhd = lkhd - dlog (dasum (nqd, eta2(1,i), 1))
      else
      lkhd = lkhd - dlog (dasum (nqd, eta2(1,i), 1)) * dfloat (cnt2(i))
      endif
23239 i=i+1
      goto 23238
23240 continue
      i=1
23243 if(.not.(i.le.n3))goto 23245
      tmp = dasum (nqd, wt3(1,1,i), 1) * dasum (nqd, wt3(1,2,i), 1)
      if(cntsum3.eq.0)then
      lkhdnew = lkhdnew - dlog (tmp)
      else
      lkhdnew = lkhdnew - dlog (tmp) * dfloat (cnt3(i))
      endif
23244 i=i+1
      goto 23243
23245 continue
      lkhd = lkhd / dfloat (nobs)
      i=1
23248 if(.not.(i.le.nt))goto 23250
      tmp = dasum (nqd, qdwt(1,1,i), 1) * dasum (nqd, qdwt(1,2,i), 1)
      lkhd = lkhd + dlog (tmp) * twt(i)
23249 i=i+1
      goto 23248
23250 continue
      iter = 0
      flag = 2
      else
      info = 2
      goto 23056
      endif
23055 goto 23054
23056 continue
      trc = 0.d0
      i=1
23251 if(.not.(i.le.n0))goto 23253
      call dcopy (nxis, rs0(1,i), 1, muwk, 1)
      if(nt.gt.1)then
      call daxpy (nxis, -1.d0, mut(1,tind(i)), 1, muwk, 1)
      else
      call daxpy (nxis, -1.d0, mut, 1, muwk, 1)
      endif
      call daxpy (nxis, -1.d0, mrs2, 1, muwk, 1)
      call dprmut (muwk, nxis, jpvt, 0)
      if(cntsum0.ne.0)then
      call dscal (nxis, dsqrt(dfloat(cnt0(i))), muwk, 1)
      endif
      call dtrsl (v, nxis, nxis, muwk, 11, infowk)
      if(nxis-rkv.gt.0)then
      call dset (nxis-rkv, 0.d0, muwk(rkv+1), 1)
      endif
      trc = trc + ddot (nxis, muwk, 1, muwk, 1)
23252 i=i+1
      goto 23251
23253 continue
      i=1
23260 if(.not.(i.le.n1))goto 23262
      tmp = dasum (nqd, eta1(1,i), 1)
      j=1
23263 if(.not.(j.le.nxis))goto 23265
      muwk(j) = ddot (nqd, eta1(1,i), 1, qdrs1(1,j,i), 1) / tmp
23264 j=j+1
      goto 23263
23265 continue
      if(nt.gt.1)then
      call daxpy (nxis, -1.d0, mut(1,tind(n0+i)), 1, muwk, 1)
      else
      call daxpy (nxis, -1.d0, mut, 1, muwk, 1)
      endif
      call daxpy (nxis, -1.d0, msr2, 1, muwk, 1)
      call dprmut (muwk, nxis, jpvt, 0)
      if(cntsum1.ne.0)then
      call dscal (nxis, dsqrt(dfloat(cnt1(i))), muwk, 1)
      endif
      call dtrsl (v, nxis, nxis, muwk, 11, infowk)
      if(nxis-rkv.gt.0)then
      call dset (nxis-rkv, 0.d0, muwk(rkv+1), 1)
      endif
      trc = trc + ddot (nxis, muwk, 1, muwk, 1)
23261 i=i+1
      goto 23260
23262 continue
      i=1
23272 if(.not.(i.le.n2))goto 23274
      tmp = dasum (nqd, eta2(1,i), 1)
      j=1
23275 if(.not.(j.le.nxis))goto 23277
      muwk(j) = ddot (nqd, eta2(1,i), 1, qdrs2(1,j,i), 1) / tmp
23276 j=j+1
      goto 23275
23277 continue
      if(nt.gt.1)then
      call daxpy (nxis, -1.d0, mut(1,tind(n0+n1+i)), 1, muwk, 1)
      else
      call daxpy (nxis, -1.d0, mut, 1, muwk, 1)
      endif
      call daxpy (nxis, -1.d0, mrs2, 1, muwk, 1)
      call dprmut (muwk, nxis, jpvt, 0)
      if(cntsum2.ne.0)then
      call dscal (nxis, dsqrt(dfloat(cnt2(i))), muwk, 1)
      endif
      call dtrsl (v, nxis, nxis, muwk, 11, infowk)
      if(nxis-rkv.gt.0)then
      call dset (nxis-rkv, 0.d0, muwk(rkv+1), 1)
      endif
      trc = trc + ddot (nxis, muwk, 1, muwk, 1)
23273 i=i+1
      goto 23272
23274 continue
      i=1
23284 if(.not.(i.le.n3))goto 23286
      j=1
23287 if(.not.(j.le.nqd))goto 23289
      k=1
23290 if(.not.(k.le.nqd))goto 23292
      wk((k-1)*nqd+j) = eta((k-1)*nqd+j) * wt3(j,1,i) * wt3(k,2,i)
23291 k=k+1
      goto 23290
23292 continue
23288 j=j+1
      goto 23287
23289 continue
      tmp = dasum (nqd*nqd, wk, 1)
      j=1
23293 if(.not.(j.le.nxis))goto 23295
      muwk(j) = ddot (nqd*nqd, wk, 1, qdrs(j,1), nxis) / tmp
23294 j=j+1
      goto 23293
23295 continue
      if(nt.gt.1)then
      call daxpy (nxis, -1.d0, mut(1,tind(n0+n1+n2+i)), 1, muwk, 1)
      else
      call daxpy (nxis, -1.d0, mut, 1, muwk, 1)
      endif
      call daxpy (nxis, -1.d0, mrs2, 1, muwk, 1)
      call dprmut (muwk, nxis, jpvt, 0)
      if(cntsum3.ne.0)then
      call dscal (nxis, dsqrt(dfloat(cnt3(i))), muwk, 1)
      endif
      call dtrsl (v, nxis, nxis, muwk, 11, infowk)
      if(nxis-rkv.gt.0)then
      call dset (nxis-rkv, 0.d0, muwk(rkv+1), 1)
      endif
      trc = trc + ddot (nxis, muwk, 1, muwk, 1)
23285 i=i+1
      goto 23284
23286 continue
      trc = trc / dfloat(nobs) / (dfloat(nobs)-1.d0)
      call dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, muwk, 1)
      lkhd = lkhd - ddot (nxi, cdnew, 1, muwk, 1) / 2.d0
      mrs(1) = lkhd
      mrs(2) = trc
      return
      end
