C Output from Public domain Ratfor, version 1.01
      subroutine cdennewton10 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, 
     *intrs, prec, maxiter, mchpr, jpvt, wk, info)
      integer nxis, nxi, nobs, cntsum, cnt(*), maxiter, jpvt(*), info
      double precision cd(*), q(nxi,*), rs(nobs,*), intrs(*), prec, mchp
     *r, wk(*)
      integer iwt, imu, iv, icdnew, iwtnew, iwk
      iwt = 1
      imu = iwt + nobs
      iv = imu + nxis
      icdnew = iv + nxis*nxis
      iwtnew = icdnew + nxis
      iwk = iwtnew + nobs
      call cdennewton101 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, intrs
     *, prec, maxiter, mchpr, wk(iwt), wk(imu), wk(iv), jpvt, wk(icdnew)
     *, wk(iwtnew), wk(iwk), info)
      return
      end
      subroutine cdennewton101 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt,
     * intrs, prec, maxiter, mchpr, wt, mu, v, jpvt, cdnew, wtnew, wk, i
     *nfo)
      integer nxis, nxi, nobs, cntsum, cnt(*), maxiter, jpvt(*), info
      double precision cd(*), q(nxi,*), rs(nobs,*), intrs(*), prec, mchp
     *r, wt(*), mu(*), v(nxis,*), cdnew(*), wtnew(*), wk(*)
      integer i, j, k, iter, flag, rkv, idamax, infowk
      double precision tmp, ddot, dasum, wtsum, lkhd, mumax, wtsumnew, l
     *khdnew, disc, disc0
      info = 0
      i=1
23000 if(.not.(i.le.nobs))goto 23002
      tmp = ddot (nxis, rs(i,1), nobs, cd, 1)
      wt(i) = dexp (-tmp)
      if(cntsum.ne.0)then
      wt(i) = wt(i) * dfloat (cnt(i))
      endif
23001 i=i+1
      goto 23000
23002 continue
      wtsum = dasum (nobs, wt, 1)
      lkhd = dlog (wtsum) + ddot (nxis, intrs, 1, cd, 1)
      call dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, wk, 1)
      lkhd = lkhd + ddot (nxi, cd, 1, wk, 1) / 2.d0
      call dscal (nobs, 1.d0/wtsum, wt, 1)
      iter = 0
      flag = 0
23005 continue
      iter = iter + 1
      i=1
23008 if(.not.(i.le.nxis))goto 23010
      mu(i) = ddot (nobs, wt, 1, rs(1,i), 1)
23009 i=i+1
      goto 23008
23010 continue
      i=1
23011 if(.not.(i.le.nxis))goto 23013
      j=i
23014 if(.not.(j.le.nxis))goto 23016
      v(i,j) = 0.d0
      k=1
23017 if(.not.(k.le.nobs))goto 23019
      v(i,j) = v(i,j) + wt(k) * rs(k,i) * rs(k,j)
23018 k=k+1
      goto 23017
23019 continue
      v(i,j) = v(i,j) - mu(i) * mu(j)
      if(j.le.nxi)then
      v(i,j) = v(i,j) + q(i,j)
      endif
23015 j=j+1
      goto 23014
23016 continue
23012 i=i+1
      goto 23011
23013 continue
      call daxpy (nxis, -1.d0, intrs, 1, mu, 1)
      call dsymv ('u', nxi, -1.d0, q, nxi, cd, 1, 1.d0, mu, 1)
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23022 if(.not.(i.le.nxis))goto 23024
      jpvt(i) = 0
23023 i=i+1
      goto 23022
23024 continue
      call dchdc (v, nxis, nxis, wk, jpvt, 1, rkv)
23025 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23025
      endif
23026 continue
      i=rkv+1
23027 if(.not.(i.le.nxis))goto 23029
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23028 i=i+1
      goto 23027
23029 continue
23030 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dset (nxis-rkv, 0.d0, cdnew(rkv+1), 1)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      i=1
23033 if(.not.(i.le.nobs))goto 23035
      tmp = ddot (nxis, rs(i,1), nobs, cdnew, 1)
      if(-tmp.gt.3.d2)then
      flag = flag + 1
      goto 23035
      endif
      wtnew(i) = dexp (-tmp)
      if(cntsum.ne.0)then
      wtnew(i) = wtnew(i) * dfloat (cnt(i))
      endif
23034 i=i+1
      goto 23033
23035 continue
      wtsumnew = dasum (nobs, wtnew, 1)
      lkhdnew = dlog (wtsumnew) + ddot (nxis, intrs, 1, cdnew, 1)
      call dsymv ('u', nxi, 1.d0, q, nxi, cdnew, 1, 0.d0, wk, 1)
      lkhdnew = lkhdnew + ddot (nxi, cdnew, 1, wk, 1) / 2.d0
      call dscal (nobs, 1.d0/wtsumnew, wtnew, 1)
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      if(cntsum.ne.0)then
      i=1
23044 if(.not.(i.le.nobs))goto 23046
      wt(i) = dfloat (cnt(i))
23045 i=i+1
      goto 23044
23046 continue
      else
      call dset (nobs, 1.d0, wt, 1)
      endif
      wtsum = dasum (nobs, wt, 1)
      lkhd = dlog (wtsum)
      call dscal (nobs, 1.d0/wtsum, wt, 1)
      iter = 0
      goto 23032
      endif
      if(flag.eq.3)then
      goto 23032
      endif
      if(lkhdnew-lkhd.lt.1.d1*(1.d0+dabs(lkhd))*mchpr)then
      goto 23032
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/(1.d0+mumax)).lt.1.d1*mchpr)then
      goto 23032
      endif
23031 goto 23030
23032 continue
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
23057 if(.not.(i.le.nobs))goto 23059
      disc = dmax1 (disc, dabs(wt(i)-wtnew(i))/(1.d0+dabs(wt(i))))
23058 i=i+1
      goto 23057
23059 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
      disc0 = dmax1 ((mumax/(1.d0+dabs(lkhd)))**2, dabs(lkhd-lkhdnew)/(1
     *.d0+dabs(lkhd)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nobs, wtnew, 1, wt, 1)
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
      if(cntsum.ne.0)then
      i=1
23070 if(.not.(i.le.nobs))goto 23072
      wt(i) = dfloat (cnt(i))
23071 i=i+1
      goto 23070
23072 continue
      else
      call dset (nobs, 1.d0, wt, 1)
      endif
      wtsum = dasum (nobs, wt, 1)
      lkhd = dlog (wtsum)
      call dscal (nobs, 1.d0/wtsum, wt, 1)
      iter = 0
      flag = 2
      else
      info = 2
      goto 23007
      endif
23006 goto 23005
23007 continue
      i=1
23073 if(.not.(i.le.nxis))goto 23075
      j=i
23076 if(.not.(j.le.nxis))goto 23078
      v(i,j) = 0.d0
      k=1
23079 if(.not.(k.le.nobs))goto 23081
      v(i,j) = v(i,j) + wt(k) * rs(k,i) * rs(k,j)
23080 k=k+1
      goto 23079
23081 continue
      if(j.le.nxi)then
      v(i,j) = v(i,j) + q(i,j)
      endif
23077 j=j+1
      goto 23076
23078 continue
23074 i=i+1
      goto 23073
23075 continue
      i=1
23084 if(.not.(i.le.nxis))goto 23086
      jpvt(i) = 0
23085 i=i+1
      goto 23084
23086 continue
      call dchdc (v, nxis, nxis, wk, jpvt, 1, rkv)
23087 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23087
      endif
23088 continue
      i=rkv+1
23089 if(.not.(i.le.nxis))goto 23091
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23090 i=i+1
      goto 23089
23091 continue
      i=1
23092 if(.not.(i.le.nobs))goto 23094
      call dcopy (nxis, rs(i,1), nobs, wk, 1)
      call dprmut (wk, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, wk, 11, infowk)
      call dset (nxis-rkv, 0.d0, wk(rkv+1), 1)
      wtnew(i) = wt(i) * ddot (nxis, wk, 1, wk, 1)
      if(cntsum.ne.0)then
      wtnew(i) = wtnew(i) / dfloat (cnt(i))
      endif
23093 i=i+1
      goto 23092
23094 continue
      call dcopy (nobs, wtnew, 1, wt, 1)
      return
      end
