C Output from Public domain Ratfor, version 1.01
      subroutine dnewton9 (dc, nsr, sr, q, nobs, ez, qdsr, nqd, qdwt, pr
     *ec, maxiter, mchpr, jpvt, wk, info)
      integer nsr, nobs, nqd, maxiter, jpvt(*), info
      double precision dc(*), q(nsr-1,*), sr(nsr,*), ez(*), qdsr(nqd,*),
     * qdwt(*), prec, mchpr, wk(*)
      integer imsr, iwt, ifit, imu, iv, idcnew, iwtnew, ifitnew, iwk
      imsr = 1
      iwt = imsr + nsr
      ifit = iwt + nqd
      imu = ifit + nobs
      iv = imu + nsr
      idcnew = iv + nsr*nsr
      iwtnew = idcnew + nsr
      ifitnew = iwtnew + nqd
      iwk = ifitnew + nobs
      call dnewton91 (dc, nsr, sr, q, nobs, ez, qdsr, nqd, qdwt, prec, m
     *axiter, mchpr, wk(imsr), wk(iwt), wk(ifit), wk(imu), wk(iv), jpvt,
     * wk(idcnew), wk(iwtnew), wk(ifitnew), wk(iwk), info)
      return
      end
      subroutine dnewton91 (dc, nsr, sr, q, nobs, ez, qdsr, nqd, qdwt, p
     *rec, maxiter, mchpr, msr, wt, fit, mu, v, jpvt, dcnew, wtnew, fitn
     *ew, wk, info)
      integer nsr, nobs, nqd, maxiter, jpvt(*), info
      double precision dc(*), q(nsr-1,*), sr(nsr,*), ez(*), qdsr(nqd,*),
     * qdwt(*), prec, mchpr, msr(*), wt(*), fit(*), mu(*), v(nsr,*), dcn
     *ew(*), wtnew(*), fitnew(*), wk(*)
      integer nr, i, j, k, iter, flag, rkv, idamax, infowk
      double precision dasum, tmp, ddot, ezsum, wtsum, fitmean, lkhd, mu
     *max, wtsumnew, lkhdnew, disc, disc0, trc
      info = 0
      nr = nsr -1
      ezsum = dasum (nobs, ez, 1)
      call dscal (nobs, 1.d0/ezsum, ez, 1)
      call dset (nsr, 0.d0, msr, 1)
      i=1
23000 if(.not.(i.le.nobs))goto 23002
      call daxpy (nsr, ez(i), sr(1,i), 1, msr, 1)
23001 i=i+1
      goto 23000
23002 continue
      i=1
23003 if(.not.(i.le.nqd))goto 23005
      wt(i) = dexp (ddot (nsr, qdsr(i,1), nqd, dc, 1)) * qdwt(i)
23004 i=i+1
      goto 23003
23005 continue
      wtsum = dasum (nqd, wt, 1)
      fitmean = 0.d0
      i=1
23006 if(.not.(i.le.nobs))goto 23008
      tmp = ddot (nsr, sr(1,i), 1, dc, 1)
      fit(i) = dexp (tmp)
      fitmean = fitmean + tmp * ez(i)
23007 i=i+1
      goto 23006
23008 continue
      call dsymv ('u', nr, 1.d0, q, nr, dc(2), 1, 0.d0, wk, 1)
      lkhd = ddot (nr, dc(2), 1, wk, 1) / 2.d0 - fitmean + dlog (wtsum)
      iter = 0
      flag = 0
23009 continue
      iter = iter + 1
      i=1
23012 if(.not.(i.le.nsr))goto 23014
      mu(i) = - ddot (nqd, wt, 1, qdsr(1,i), 1) / wtsum
23013 i=i+1
      goto 23012
23014 continue
      i=1
23015 if(.not.(i.le.nsr))goto 23017
      j=i
23018 if(.not.(j.le.nsr))goto 23020
      v(i,j) = 0.d0
      k=1
23021 if(.not.(k.le.nqd))goto 23023
      v(i,j) = v(i,j) + wt(k) * qdsr(k,i) * qdsr(k,j)
23022 k=k+1
      goto 23021
23023 continue
      v(i,j) = v(i,j) / wtsum - mu(i) * mu(j)
23019 j=j+1
      goto 23018
23020 continue
23016 i=i+1
      goto 23015
23017 continue
      i=2
23024 if(.not.(i.le.nsr))goto 23026
      j=i
23027 if(.not.(j.le.nsr))goto 23029
      v(i,j) = v(i,j) + q(i-1,j-1)
23028 j=j+1
      goto 23027
23029 continue
23025 i=i+1
      goto 23024
23026 continue
      call daxpy (nsr, 1.d0, msr, 1, mu, 1)
      call dsymv ('u', nr, -1.d0, q, nr, dc(2), 1, 1.d0, mu(2), 1)
      mumax = dabs(mu(idamax(nsr, mu, 1)))
      i=1
23030 if(.not.(i.le.nsr))goto 23032
      jpvt(i) = 0
23031 i=i+1
      goto 23030
23032 continue
      call dchdc (v, nsr, nsr, wk, jpvt, 1, rkv)
23033 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23033
      endif
23034 continue
      i=rkv+1
23035 if(.not.(i.le.nsr))goto 23037
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23036 i=i+1
      goto 23035
23037 continue
23038 continue
      call dcopy (nsr, mu, 1, dcnew, 1)
      call dprmut (dcnew, nsr, jpvt, 0)
      call dtrsl (v, nsr, nsr, dcnew, 11, infowk)
      call dset (nsr-rkv, 0.d0, dcnew(rkv+1), 1)
      call dtrsl (v, nsr, nsr, dcnew, 01, infowk)
      call dprmut (dcnew, nsr, jpvt, 1)
      call daxpy (nsr, 1.d0, dc, 1, dcnew, 1)
      i=1
23041 if(.not.(i.le.nqd))goto 23043
      tmp = ddot (nsr, qdsr(i,1), nqd, dcnew, 1)
      if(tmp.gt.3.d2)then
      flag = flag + 1
      goto 23043
      endif
      wtnew(i) = qdwt(i) * dexp (tmp)
23042 i=i+1
      goto 23041
23043 continue
      wtsumnew = dasum (nqd, wtnew, 1)
      if((flag.eq.0).or.(flag.eq.2))then
      fitmean = 0.d0
      i=1
23048 if(.not.(i.le.nobs))goto 23050
      tmp = ddot (nsr, sr(1,i), 1, dcnew, 1)
      if(tmp.gt.3.d2)then
      flag = flag + 1
      goto 23050
      endif
      fitnew(i) = dexp (tmp)
      fitmean = fitmean + tmp * ez(i)
23049 i=i+1
      goto 23048
23050 continue
      call dsymv ('u', nr, 1.d0, q, nr, dcnew(2), 1, 0.d0, wk, 1)
      lkhdnew = ddot (nr, dcnew(2), 1, wk, 1) / 2.d0 - fitmean + dlog (w
     *tsumnew)
      endif
      if(flag.eq.1)then
      call dset (nsr, 0.d0, dc, 1)
      call dcopy (nqd, qdwt, 1, wt, 1)
      wtsum = dasum (nqd, wt, 1)
      lkhd = dlog (wtsum)
      call dset (nobs, 1.d0, fit, 1)
      iter = 0
      goto 23040
      endif
      if(flag.eq.3)then
      goto 23040
      endif
      if(lkhdnew-lkhd.lt.1.d1*(1.d0+dabs(lkhd))*mchpr)then
      goto 23040
      endif
      call dscal (nsr, .5d0, mu, 1)
      if(dabs(mu(idamax(nsr, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23040
      endif
23039 goto 23038
23040 continue
      if(flag.eq.1)then
      flag = 2
      goto 23010
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      i=1
23065 if(.not.(i.le.nqd))goto 23067
      disc = dmax1 (disc, dabs(wt(i)-wtnew(i))/(1.d0+dabs(wt(i))))
23066 i=i+1
      goto 23065
23067 continue
      i=1
23068 if(.not.(i.le.nobs))goto 23070
      disc = dmax1 (disc, dabs(fit(i)-fitnew(i))/(1.d0+dabs(fit(i))))
23069 i=i+1
      goto 23068
23070 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
      disc0 = dmax1 ((mumax/(1.d0+lkhd))**2, dabs(lkhd-lkhdnew)/(1.d0+da
     *bs(lkhd)))
      call dcopy (nsr, dcnew, 1, dc, 1)
      call dcopy (nqd, wtnew, 1, wt, 1)
      wtsum = wtsumnew
      call dcopy (nobs, fitnew, 1, fit, 1)
      lkhd = lkhdnew
      if(disc0.lt.prec)then
      goto 23011
      endif
      if(disc.lt.prec)then
      goto 23011
      endif
      if(iter.lt.maxiter)then
      goto 23010
      endif
      if(flag.eq.0)then
      call dset (nsr, 0.d0, dc, 1)
      call dcopy (nqd, qdwt, 1, wt, 1)
      wtsum = dasum (nqd, wt, 1)
      lkhd = dlog (wtsum)
      call dset (nobs, 1.d0, fit, 1)
      iter = 0
      flag = 2
      else
      info = 2
      goto 23011
      endif
23010 goto 23009
23011 continue
      i=1
23079 if(.not.(i.le.nobs))goto 23081
      call daxpy (nsr, -1.d0, msr, 1, sr(1,i), 1)
      call dprmut (sr(1,i), nsr, jpvt, 0)
      call dscal (nsr, dsqrt(ez(i)), sr(1,i), 1)
      call dtrsl (v, nsr, nsr, sr(1,i), 11, infowk)
      if(nsr-rkv.gt.0)then
      call dset (nsr-rkv, 0.d0, sr(rkv+1,i), 1)
      endif
23080 i=i+1
      goto 23079
23081 continue
      trc = ddot (nobs*nsr, sr, 1, sr, 1)
      lkhd = 0.d0
      i=1
23084 if(.not.(i.le.nobs))goto 23086
      lkhd = lkhd + dlog (fit(i)) * ez(i)
23085 i=i+1
      goto 23084
23086 continue
      lkhd = lkhd - dlog (wtsum)
      msr(1) = lkhd
      msr(2) = trc / dmax1(1.d0, ezsum-1.d0)
      return
      end
