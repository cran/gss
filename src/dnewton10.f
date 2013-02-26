C Output from Public domain Ratfor, version 1.01
      subroutine dnewton10 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, int
     *rs, prec, maxiter, mchpr, jpvt, wk, info)
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
      call dnewton101 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, intrs, p
     *rec, maxiter, mchpr, wk(iwt), wk(imu), wk(iv), jpvt, wk(icdnew), w
     *k(iwtnew), wk(iwk), info)
      return
      end
      subroutine dnewton101 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, in
     *trs, prec, maxiter, mchpr, wt, mu, v, jpvt, cdnew, wtnew, wk, info
     *)
      integer nxis, nxi, nobs, cntsum, cnt(*), maxiter, jpvt(*), info
      double precision cd(*), q(nxi,*), rs(nobs,*), intrs(*), prec, mchp
     *r, wt(*), mu(*), v(nxis,*), cdnew(*), wtnew(*), wk(*)
      integer i, j, k, iter, flag, rkv, idamax, infowk
      double precision wtsum, tmp, ddot, lkhd, mumax, wtsumnew, lkhdnew,
     * disc, disc0
      info = 0
      wtsum = 0.d0
      i=1
23000 if(.not.(i.le.nobs))goto 23002
      tmp = ddot (nxis, rs(i,1), nobs, cd, 1)
      wt(i) = dexp (-tmp)
      if(cntsum.ne.0)then
      wt(i) = wt(i) * dfloat (cnt(i))
      endif
      wtsum = wtsum + wt(i)
23001 i=i+1
      goto 23000
23002 continue
      if(cntsum.eq.0)then
      lkhd = wtsum / dfloat (nobs)
      else
      lkhd = wtsum / dfloat (cntsum)
      endif
      lkhd = dlog (lkhd) + ddot (nxis, intrs, 1, cd, 1)
      call dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, wk, 1)
      lkhd = lkhd + ddot (nxi, cd, 1, wk, 1) / 2.d0
      iter = 0
      flag = 0
23007 continue
      iter = iter + 1
      i=1
23010 if(.not.(i.le.nxis))goto 23012
      mu(i) = ddot (nobs, wt, 1, rs(1,i), 1) / wtsum
23011 i=i+1
      goto 23010
23012 continue
      i=1
23013 if(.not.(i.le.nxis))goto 23015
      j=i
23016 if(.not.(j.le.nxis))goto 23018
      v(i,j) = 0.d0
      k=1
23019 if(.not.(k.le.nobs))goto 23021
      v(i,j) = v(i,j) + wt(k) * rs(k,i) * rs(k,j)
23020 k=k+1
      goto 23019
23021 continue
      v(i,j) = v(i,j) / wtsum - mu(i) * mu(j)
      if(j.le.nxi)then
      v(i,j) = v(i,j) + q(i,j)
      endif
23017 j=j+1
      goto 23016
23018 continue
23014 i=i+1
      goto 23013
23015 continue
      call daxpy (nxis, -1.d0, intrs, 1, mu, 1)
      call dsymv ('u', nxi, -1.d0, q, nxi, cd, 1, 1.d0, mu, 1)
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23024 if(.not.(i.le.nxis))goto 23026
      jpvt(i) = 0
23025 i=i+1
      goto 23024
23026 continue
      call dchdc (v, nxis, nxis, wk, jpvt, 1, rkv)
23027 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23027
      endif
23028 continue
      i=rkv+1
23029 if(.not.(i.le.nxis))goto 23031
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23030 i=i+1
      goto 23029
23031 continue
23032 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dset (nxis-rkv, 0.d0, cdnew(rkv+1), 1)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      wtsumnew = 0.d0
      i=1
23035 if(.not.(i.le.nobs))goto 23037
      tmp = ddot (nxis, rs(i,1), nobs, cdnew, 1)
      if(-tmp.gt.3.d2)then
      flag = flag + 1
      goto 23037
      endif
      wtnew(i) = dexp (-tmp)
      if(cntsum.ne.0)then
      wtnew(i) = wtnew(i) * dfloat (cnt(i))
      endif
      wtsumnew = wtsumnew + wtnew(i)
23036 i=i+1
      goto 23035
23037 continue
      if(cntsum.eq.0)then
      lkhdnew = wtsumnew / dfloat (nobs)
      else
      lkhdnew = wtsumnew / dfloat (cntsum)
      endif
      lkhdnew = dlog (lkhdnew) + ddot (nxis, intrs, 1, cdnew, 1)
      call dsymv ('u', nxi, 1.d0, q, nxi, cdnew, 1, 0.d0, wk, 1)
      lkhdnew = lkhdnew + ddot (nxi, cdnew, 1, wk, 1) / 2.d0
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      wtsum = 0.d0
      i=1
23046 if(.not.(i.le.nobs))goto 23048
      if(cntsum.ne.0)then
      wt(i) = dfloat (cnt(i))
      else
      wt(i) = 1.d0
      endif
      wtsum = wtsum + wt(i)
23047 i=i+1
      goto 23046
23048 continue
      lkhd = 0.d0
      iter = 0
      goto 23034
      endif
      if(flag.eq.3)then
      goto 23034
      endif
      if(lkhdnew-lkhd.lt.1.d1*(1.d0+dabs(lkhd))*mchpr)then
      goto 23034
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/(1.d0+mumax)).lt.1.d1*mchpr)then
      goto 23034
      endif
23033 goto 23032
23034 continue
      if(flag.eq.1)then
      flag = 2
      goto 23008
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      i=1
23061 if(.not.(i.le.nobs))goto 23063
      disc = dmax1 (disc, dabs(wt(i)-wtnew(i))/(1.d0+dabs(wt(i))))
23062 i=i+1
      goto 23061
23063 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
      disc0 = dmax1 ((mumax/(1.d0+dabs(lkhd)))**2, dabs(lkhd-lkhdnew)/(1
     *.d0+dabs(lkhd)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nobs, wtnew, 1, wt, 1)
      wtsum = wtsumnew
      lkhd = lkhdnew
      if(disc0.lt.prec)then
      goto 23009
      endif
      if(disc.lt.prec)then
      goto 23009
      endif
      if(iter.lt.maxiter)then
      goto 23008
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      wtsum = 0.d0
      i=1
23072 if(.not.(i.le.nobs))goto 23074
      if(cntsum.ne.0)then
      wt(i) = dfloat (cnt(i))
      else
      wt(i) = 1.d0
      endif
      wtsum = wtsum + wt(i)
23073 i=i+1
      goto 23072
23074 continue
      lkhd = 0.d0
      iter = 0
      flag = 2
      else
      info = 2
      goto 23009
      endif
23008 goto 23007
23009 continue
      call dscal (nobs, 1.d0/wtsum, wt, 1)
      i=1
23077 if(.not.(i.le.nxis))goto 23079
      j=i
23080 if(.not.(j.le.nxis))goto 23082
      v(i,j) = 0.d0
      k=1
23083 if(.not.(k.le.nobs))goto 23085
      v(i,j) = v(i,j) + wt(k) * rs(k,i) * rs(k,j)
23084 k=k+1
      goto 23083
23085 continue
      if(j.le.nxi)then
      v(i,j) = v(i,j) + q(i,j)
      endif
23081 j=j+1
      goto 23080
23082 continue
23078 i=i+1
      goto 23077
23079 continue
      i=1
23088 if(.not.(i.le.nxis))goto 23090
      jpvt(i) = 0
23089 i=i+1
      goto 23088
23090 continue
      call dchdc (v, nxis, nxis, wk, jpvt, 1, rkv)
23091 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23091
      endif
23092 continue
      i=rkv+1
23093 if(.not.(i.le.nxis))goto 23095
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23094 i=i+1
      goto 23093
23095 continue
      i=1
23096 if(.not.(i.le.nobs))goto 23098
      call dcopy (nxis, rs(i,1), nobs, wk, 1)
      call dprmut (wk, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, wk, 11, infowk)
      call dset (nxis-rkv, 0.d0, wk(rkv+1), 1)
      wtnew(i) = wt(i) * ddot (nxis, wk, 1, wk, 1)
      if(cntsum.ne.0)then
      wtnew(i) = wtnew(i) / dfloat (cnt(i))
      endif
23097 i=i+1
      goto 23096
23098 continue
      call dcopy (nobs, wtnew, 1, wt, 1)
      return
      end
