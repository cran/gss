C Output from Public domain Ratfor, version 1.01
      subroutine hzdnewton10 (cd, nxis, q, nxi, rs, nt, nobs, cntsum, cn
     *t, intrs, rho, prec, maxiter, mchpr, jpvt, wk, info)
      integer nxis, nxi, nt, nobs, cntsum, cnt(*), maxiter, jpvt(*), inf
     *o
      double precision cd(*), q(nxi,*), rs(nt,*), intrs(*), rho(*), prec
     *, mchpr, wk(*)
      integer iwt, imu, iv, icdnew, iwtnew, iwk
      iwt = 1
      imu = iwt + nt
      iv = imu + nxis
      icdnew = iv + nxis*nxis
      iwtnew = icdnew + nxis
      iwk = iwtnew + nt
      call hzdnewton101 (cd, nxis, q, nxi, rs, nt, nobs, cntsum, cnt, in
     *trs, rho, prec, maxiter, mchpr, wk(iwt), wk(imu), wk(iv), jpvt, wk
     *(icdnew), wk(iwtnew), wk(iwk), info)
      return
      end
      subroutine hzdnewton101 (cd, nxis, q, nxi, rs, nt, nobs, cntsum, c
     *nt, intrs, rho, prec, maxiter, mchpr, wt, mu, v, jpvt, cdnew, wtne
     *w, wk, info)
      integer nxis, nxi, nt, nobs, cntsum, cnt(*), maxiter, jpvt(*), inf
     *o
      double precision cd(*), q(nxi,*), rs(nt,*), intrs(*), rho(*), prec
     *, mchpr, wt(*), mu(*), v(nxis,*), cdnew(*), wtnew(*), wk(*)
      integer i, j, k, iter, flag, rkv, idamax, infowk
      double precision tmp, ddot, dasum, lkhd, mumax, lkhdnew, disc, dis
     *c0
      info = 0
      i=1
23000 if(.not.(i.le.nt))goto 23002
      tmp = ddot (nxis, rs(i,1), nt, cd, 1)
      wt(i) = dexp (-tmp) * rho(i)
      if(cntsum.ne.0)then
      wt(i) = wt(i) * dfloat (cnt(i))
      endif
23001 i=i+1
      goto 23000
23002 continue
      call dscal (nt, 1/dfloat(nobs), wt, 1)
      lkhd = dasum(nt, wt, 1) + ddot (nxis, intrs, 1, cd, 1)
      call dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, wk, 1)
      lkhd = lkhd + ddot (nxi, cd, 1, wk, 1) / 2.d0
      iter = 0
      flag = 0
23005 continue
      iter = iter + 1
      i=1
23008 if(.not.(i.le.nxis))goto 23010
      mu(i) = ddot (nt, wt, 1, rs(1,i), 1)
      j=i
23011 if(.not.(j.le.nxis))goto 23013
      v(i,j) = 0.d0
      k=1
23014 if(.not.(k.le.nt))goto 23016
      v(i,j) = v(i,j) + wt(k) * rs(k,i) * rs(k,j)
23015 k=k+1
      goto 23014
23016 continue
      if(j.le.nxi)then
      v(i,j) = v(i,j) + q(i,j)
      endif
23012 j=j+1
      goto 23011
23013 continue
23009 i=i+1
      goto 23008
23010 continue
      call daxpy (nxis, -1.d0, intrs, 1, mu, 1)
      call dsymv ('u', nxi, -1.d0, q, nxi, cd, 1, 1.d0, mu, 1)
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23019 if(.not.(i.le.nxis))goto 23021
      jpvt(i) = 0
23020 i=i+1
      goto 23019
23021 continue
      call dchdc (v, nxis, nxis, wk, jpvt, 1, rkv)
23022 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23022
      endif
23023 continue
      i=rkv+1
23024 if(.not.(i.le.nxis))goto 23026
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23025 i=i+1
      goto 23024
23026 continue
23027 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dset (nxis-rkv, 0.d0, cdnew(rkv+1), 1)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      i=1
23030 if(.not.(i.le.nt))goto 23032
      tmp = ddot (nxis, rs(i,1), nt, cdnew, 1)
      if(-tmp.gt.3.d2)then
      flag = flag + 1
      goto 23032
      endif
      wtnew(i) = dexp (-tmp) * rho(i)
      if(cntsum.ne.0)then
      wtnew(i) = wtnew(i) * dfloat (cnt(i))
      endif
23031 i=i+1
      goto 23030
23032 continue
      call dscal (nt, 1/dfloat(nobs), wtnew, 1)
      lkhdnew = dasum(nt, wtnew, 1) + ddot (nxis, intrs, 1, cdnew, 1)
      call dsymv ('u', nxi, 1.d0, q, nxi, cdnew, 1, 0.d0, wk, 1)
      lkhdnew = lkhdnew + ddot (nxi, cdnew, 1, wk, 1) / 2.d0
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      i=1
23039 if(.not.(i.le.nt))goto 23041
      wt(i) = rho(i)
      if(cntsum.ne.0)then
      wt(i) = wt(i) * dfloat (cnt(i))
      endif
23040 i=i+1
      goto 23039
23041 continue
      call dscal (nt, 1/dfloat(nobs), wt, 1)
      lkhd = dasum (nt, wt, 1)
      iter = 0
      goto 23029
      endif
      if(flag.eq.3)then
      goto 23029
      endif
      if(lkhdnew-lkhd.lt.1.d1*(1.d0+dabs(lkhd))*mchpr)then
      goto 23029
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/(1.d0+mumax)).lt.1.d1*mchpr)then
      goto 23029
      endif
23028 goto 23027
23029 continue
      if(flag.eq.1)then
      flag = 2
      goto 23006
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      i=1
23054 if(.not.(i.le.nt))goto 23056
      disc = dmax1 (disc, dabs(wt(i)-wtnew(i))/(1.d0+dabs(wt(i))))
23055 i=i+1
      goto 23054
23056 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
      disc0 = dmax1 ((mumax/(1.d0+dabs(lkhd)))**2, dabs(lkhd-lkhdnew)/(1
     *.d0+dabs(lkhd)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nt, wtnew, 1, wt, 1)
      lkhd = lkhdnew
      if(disc0.lt.prec)then
      goto 23007
      endif
      if(disc.lt.prec)then
      goto 23007
      endif
      if(iter.lt.maxiter)then
      goto 23006
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      i=1
23065 if(.not.(i.le.nt))goto 23067
      wt(i) = rho(i)
      if(cntsum.ne.0)then
      wt(i) = wt(i) * dfloat (cnt(i))
      endif
23066 i=i+1
      goto 23065
23067 continue
      call dscal (nt, 1/dfloat(nobs), wt, 1)
      lkhd = dasum (nt, wt, 1)
      iter = 0
      flag = 2
      else
      info = 2
      goto 23007
      endif
23006 goto 23005
23007 continue
      lkhd = dasum (nt, wt, 1) + ddot (nxis, intrs, 1, cd, 1)
      tmp = 0.d0
      disc = 0.d0
      i=1
23070 if(.not.(i.le.nt))goto 23072
      call dcopy (nxis, rs(i,1), nt, wk, 1)
      call dprmut (wk, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, wk, 11, infowk)
      call dset (nxis-rkv, 0.d0, wk(rkv+1), 1)
      wtnew(i) = wt(i) * ddot (nxis, wk, 1, wk, 1)
      if(cntsum.ne.0)then
      wtnew(i) = wtnew(i) / dfloat (cnt(i))
      endif
      tmp = tmp + wt(i) * (dexp (wtnew(i)/(1.d0-wtnew(i))) - 1.d0)
      if(cntsum.ne.0)then
      disc = disc + dfloat(cnt(i)) * wtnew(i)/(1.d0-wtnew(i))
      else
      disc = disc + wtnew(i)/(1.d0-wtnew(i))
      endif
23071 i=i+1
      goto 23070
23072 continue
      wt(1) = lkhd
      wt(2) = tmp
      wt(3) = disc/dfloat(nobs)
      return
      end
      subroutine hzdaux101 (cd, nxis, q, nxi, rs, nt, rho, mchpr, v, jpv
     *t)
      integer nxis, nxi, nt, jpvt(*)
      double precision cd(*), q(nxi,*), rs(nt,*), rho(*), mchpr, v(nxis,
     **)
      integer i, j, k, rkv
      double precision tmp, ddot
      i=1
23077 if(.not.(i.le.nt))goto 23079
      tmp = ddot (nxis, rs(i,1), nt, cd, 1)
      rho(i) = dexp (-tmp) * rho(i)
23078 i=i+1
      goto 23077
23079 continue
      i=1
23080 if(.not.(i.le.nxis))goto 23082
      j=i
23083 if(.not.(j.le.nxis))goto 23085
      v(i,j) = 0.d0
      k=1
23086 if(.not.(k.le.nt))goto 23088
      v(i,j) = v(i,j) + rho(k) * rs(k,i) * rs(k,j)
23087 k=k+1
      goto 23086
23088 continue
      if(j.le.nxi)then
      v(i,j) = v(i,j) + q(i,j)
      endif
23084 j=j+1
      goto 23083
23085 continue
23081 i=i+1
      goto 23080
23082 continue
      i=1
23091 if(.not.(i.le.nxis))goto 23093
      jpvt(i) = 0
23092 i=i+1
      goto 23091
23093 continue
      call dchdc (v, nxis, nxis, cd, jpvt, 1, rkv)
23094 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23094
      endif
23095 continue
      i=rkv+1
23096 if(.not.(i.le.nxis))goto 23098
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23097 i=i+1
      goto 23096
23098 continue
      return
      end
