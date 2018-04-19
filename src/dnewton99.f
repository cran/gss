C Output from Public domain Ratfor, version 1.01
      subroutine dnewton99 (dc, nsr, sr, q, nobs, cntsum, cnt, ff, qdsr,
     * nqd, qdwt, prec, maxiter, mchpr, jpvt, wk, info)
      integer nsr, nobs, cntsum, cnt(*), nqd, maxiter, jpvt(*), info
      double precision dc(*), q(nsr-1,*), sr(nsr,*), ff(nobs,*), qdsr(nq
     *d,*), qdwt(*), prec, mchpr, wk(*)
      integer iwt, ifit, imu, imuwk, iv, ivwk, idcnew, iwtnew, ifitnew, 
     *iwk
      iwt = 1
      ifit = iwt + nqd
      imu = ifit + nobs
      imuwk = imu + nsr
      iv = imuwk + nsr
      ivwk = iv + nsr*nsr
      idcnew = ivwk + nsr*nsr
      iwtnew = idcnew + nsr
      ifitnew = iwtnew + nqd
      iwk = ifitnew + nobs
      call dnewton991 (dc, nsr, sr, q, nobs, cntsum, cnt, ff, qdsr, nqd,
     * qdwt, prec, maxiter, mchpr, wk(iwt), wk(ifit), wk(imu), wk(imuwk)
     *, wk(iv), wk(ivwk), jpvt, wk(idcnew), wk(iwtnew), wk(ifitnew), wk(
     *iwk), info)
      return
      end
      subroutine dnewton991 (dc, nsr, sr, q, nobs, cntsum, cnt, ff, qdsr
     *, nqd, qdwt, prec, maxiter, mchpr, wt, fit, mu, muwk, v, vwk, jpvt
     *, dcnew, wtnew, fitnew, wk, info)
      integer nsr, nobs, cntsum, cnt(*), nqd, maxiter, jpvt(*), info
      double precision dc(*), q(nsr-1,*), sr(nsr,*), ff(nobs,*), qdsr(nq
     *d,*), qdwt(*), prec, mchpr, wt(*), fit(*), mu(*), muwk(*), v(nsr,*
     *), vwk(nsr,*), dcnew(*), wtnew(*), fitnew(*), wk(*)
      integer nr, i, j, k, iter, flag, idamax, infowk
      double precision dasum, tmp, ddot, wtsum, lkhd, mumax, wtsumnew, l
     *khdnew, disc, disc0
      info = 0
      nr = nsr - 1
      i=1
23000 if(.not.(i.le.nqd))goto 23002
      wt(i) = dexp (ddot (nsr, qdsr(i,1), nqd, dc, 1)) * qdwt(i)
23001 i=i+1
      goto 23000
23002 continue
      wtsum = dasum (nqd, wt, 1)
      lkhd = 0.d0
      i=1
23003 if(.not.(i.le.nobs))goto 23005
      fit(i) = dexp (ddot (nsr, sr(1,i), 1, dc, 1)) / wtsum
      if(cntsum.eq.0)then
      lkhd = lkhd + dlog (ff(i,1)+ff(i,2)*fit(i))
      else
      lkhd = lkhd + dlog (ff(i,1)+ff(i,2)*fit(i)) * dfloat (cnt(i))
      endif
23004 i=i+1
      goto 23003
23005 continue
      if(cntsum.eq.0)then
      lkhd = lkhd / dfloat (nobs)
      else
      lkhd = lkhd / dfloat (cntsum)
      endif
      call dsymv ('u', nr, 1.d0, q, nr, dc(2), 1, 0.d0, wk, 1)
      lkhd = ddot (nr, dc(2), 1, wk, 1) / 2.d0 - lkhd
      iter = 0
      flag = 0
23010 continue
      iter = iter + 1
      i=1
23013 if(.not.(i.le.nsr))goto 23015
      muwk(i) = ddot (nqd, wt, 1, qdsr(1,i), 1) / wtsum
23014 i=i+1
      goto 23013
23015 continue
      i=1
23016 if(.not.(i.le.nsr))goto 23018
      j=i
23019 if(.not.(j.le.nsr))goto 23021
      vwk(i,j) = 0.d0
      k=1
23022 if(.not.(k.le.nqd))goto 23024
      vwk(i,j) = vwk(i,j) + wt(k) * qdsr(k,i) * qdsr(k,j)
23023 k=k+1
      goto 23022
23024 continue
      vwk(i,j) = vwk(i,j) / wtsum - muwk(i) * muwk(j)
23020 j=j+1
      goto 23019
23021 continue
23017 i=i+1
      goto 23016
23018 continue
      call dset (nsr, 0.d0, mu, 1)
      call dset (nsr*nsr, 0.d0, v, 1)
      if(cntsum.eq.0)then
      k=1
23027 if(.not.(k.le.nobs))goto 23029
      tmp = ff(k,2)*fit(k)
      tmp = tmp / (ff(k,1)+tmp)
      call daxpy (nsr, tmp, muwk, 1, mu, 1)
      call daxpy (nsr, -tmp, sr(1,k), 1, mu, 1)
      i=1
23030 if(.not.(i.le.nsr))goto 23032
      j=i
23033 if(.not.(j.le.nsr))goto 23035
      v(i,j) = v(i,j) + tmp*vwk(i,j)
      v(i,j) = v(i,j) - tmp*(1.d0-tmp)*(muwk(i)-sr(i,k))*(muwk(j)-sr(j,k
     *))
23034 j=j+1
      goto 23033
23035 continue
23031 i=i+1
      goto 23030
23032 continue
23028 k=k+1
      goto 23027
23029 continue
      call dscal (nsr, -1.d0/dfloat(nobs), mu, 1)
      call dscal (nsr*nsr, 1.d0/dfloat(nobs), v, 1)
      else
      k=1
23036 if(.not.(k.le.nobs))goto 23038
      tmp = ff(k,2)*fit(k)
      tmp = tmp / (ff(k,1)+tmp)
      call daxpy (nsr, tmp*dfloat(cnt(i)), muwk, 1, mu, 1)
      call daxpy (nsr, -tmp*dfloat(cnt(i)), sr(1,k), 1, mu, 1)
      i=1
23039 if(.not.(i.le.nsr))goto 23041
      j=i
23042 if(.not.(j.le.nsr))goto 23044
      v(i,j) = v(i,j) + tmp*vwk(i,j)*dfloat(cnt(i))
      v(i,j) = v(i,j) - tmp*(1.d0-tmp)*dfloat(cnt(i))* (muwk(i)-sr(i,k))
     **(muwk(j)-sr(j,k))
23043 j=j+1
      goto 23042
23044 continue
23040 i=i+1
      goto 23039
23041 continue
23037 k=k+1
      goto 23036
23038 continue
      call dscal (nsr, -1.d0/dfloat(cntsum), mu, 1)
      call dscal (nsr*nsr, 1.d0/dfloat(cntsum), v, 1)
      endif
      call dsymv ('u', nr, -1.d0, q, nr, dc(2), 1, 1.d0, mu(2), 1)
      i=2
23045 if(.not.(i.le.nsr))goto 23047
      j=i
23048 if(.not.(j.le.nsr))goto 23050
      v(i,j) = v(i,j) + q(i-1,j-1)
23049 j=j+1
      goto 23048
23050 continue
23046 i=i+1
      goto 23045
23047 continue
      mumax = dabs(mu(idamax(nsr, mu, 1)))
      call dmcdc (v, nsr, nsr, wk, jpvt, infowk)
23051 continue
      call dcopy (nsr, mu, 1, dcnew, 1)
      call dprmut (dcnew, nsr, jpvt, 0)
      call dtrsl (v, nsr, nsr, dcnew, 11, infowk)
      call dtrsl (v, nsr, nsr, dcnew, 01, infowk)
      call dprmut (dcnew, nsr, jpvt, 1)
      call daxpy (nsr, 1.d0, dc, 1, dcnew, 1)
      i=1
23054 if(.not.(i.le.nqd))goto 23056
      tmp = ddot (nsr, qdsr(i,1), nqd, dcnew, 1)
      if(tmp.gt.3.d2)then
      flag = flag + 1
      goto 23056
      endif
      wtnew(i) = qdwt(i) * dexp (tmp)
23055 i=i+1
      goto 23054
23056 continue
      wtsumnew = dasum (nqd, wtnew, 1)
      if((flag.eq.0).or.(flag.eq.2))then
      lkhdnew = 0.d0
      i=1
23061 if(.not.(i.le.nobs))goto 23063
      tmp = ddot (nsr, sr(1,i), 1, dcnew, 1)
      if(tmp.gt.3.d2)then
      flag = flag + 1
      goto 23063
      endif
      fitnew(i) = dexp (tmp) / wtsumnew
      if(cntsum.eq.0)then
      lkhdnew = lkhdnew + dlog (ff(i,1)+ff(i,2)*fitnew(i))
      else
      lkhdnew = lkhdnew + dlog (ff(i,1)+ff(i,2)*fitnew(i)) * dfloat (cnt
     *(i))
      endif
23062 i=i+1
      goto 23061
23063 continue
      if(cntsum.eq.0)then
      lkhdnew = lkhdnew / dfloat (nobs)
      else
      lkhdnew = lkhdnew / dfloat (cntsum)
      endif
      call dsymv ('u', nr, 1.d0, q, nr, dcnew(2), 1, 0.d0, wk, 1)
      lkhdnew = ddot (nr, dcnew(2), 1, wk, 1) / 2.d0 - lkhdnew
      endif
      if(flag.eq.1)then
      call dset (nsr, 0.d0, dc, 1)
      call dcopy (nqd, qdwt, 1, wt, 1)
      wtsum = dasum (nqd, wt, 1)
      call dset (nobs, 1.d0/wtsum, fit, 1)
      lkhd = 0.d0
      i=1
23072 if(.not.(i.le.nobs))goto 23074
      if(cntsum.eq.0)then
      lkhd = lkhd + dlog (ff(i,1)+ff(i,2)*fit(i))
      else
      lkhd = lkhd + dlog (ff(i,1)+ff(i,2)*fit(i)) * dfloat (cnt(i))
      endif
23073 i=i+1
      goto 23072
23074 continue
      if(cntsum.eq.0)then
      lkhd = lkhd / dfloat (nobs)
      else
      lkhd = lkhd / dfloat (cntsum)
      endif
      iter = 0
      goto 23053
      endif
      if(flag.eq.3)then
      goto 23053
      endif
      if(lkhdnew-lkhd.lt.1.d1*(1.d0+dabs(lkhd))*mchpr)then
      goto 23053
      endif
      call dscal (nsr, .5d0, mu, 1)
      if(dabs(mu(idamax(nsr, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23053
      endif
23052 goto 23051
23053 continue
      if(flag.eq.1)then
      flag = 2
      goto 23011
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      i=1
23089 if(.not.(i.le.nqd))goto 23091
      disc = dmax1 (disc, dabs(wt(i)-wtnew(i))/(1.d0+dabs(wt(i))))
23090 i=i+1
      goto 23089
23091 continue
      i=1
23092 if(.not.(i.le.nobs))goto 23094
      disc = dmax1 (disc, dabs(fit(i)-fitnew(i))/(1.d0+dabs(fit(i))))
23093 i=i+1
      goto 23092
23094 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
      disc0 = dmax1 ((mumax/(1.d0+lkhd))**2, dabs(lkhd-lkhdnew)/(1.d0+da
     *bs(lkhd)))
      call dcopy (nsr, dcnew, 1, dc, 1)
      call dcopy (nqd, wtnew, 1, wt, 1)
      wtsum = wtsumnew
      call dcopy (nobs, fitnew, 1, fit, 1)
      lkhd = lkhdnew
      if(disc0.lt.prec)then
      goto 23012
      endif
      if(disc.lt.prec)then
      goto 23012
      endif
      if(iter.lt.maxiter)then
      goto 23011
      endif
      if(flag.eq.0)then
      call dset (nsr, 0.d0, dc, 1)
      call dcopy (nqd, qdwt, 1, wt, 1)
      wtsum = dasum (nqd, wt, 1)
      call dset (nobs, 1.d0/wtsum, fit, 1)
      lkhd = 0.d0
      i=1
23103 if(.not.(i.le.nobs))goto 23105
      if(cntsum.eq.0)then
      lkhd = lkhd + dlog (ff(i,1)+ff(i,2)*fit(i))
      else
      lkhd = lkhd + dlog (ff(i,1)+ff(i,2)*fit(i)) * dfloat (cnt(i))
      endif
23104 i=i+1
      goto 23103
23105 continue
      if(cntsum.eq.0)then
      lkhd = lkhd / dfloat (nobs)
      else
      lkhd = lkhd / dfloat (cntsum)
      endif
      iter = 0
      flag = 2
      else
      info = 2
      goto 23012
      endif
23011 goto 23010
23012 continue
      return
      end
