C Output from Public domain Ratfor, version 1.0
      subroutine hzdnewton1 (cd, nxis, q, nxi, rs, nt, nobs, cntsum, cnt
     *, qdrs, nqd, qdwt, nx, prec, maxiter, mchpr, mrs, wt, fit, mu, muw
     *k, v, vwk, jpvt, cdnew, wtnew, fitnew, wk, info)
      integer nxis, nxi, nt, nobs, cntsum, cnt(*), nqd, nx, maxiter, jpv
     *t(*), info
      double precision cd(*), q(nxi,*), rs(nxis,*), qdrs(nqd,nxis,*), qd
     *wt(nqd,*), prec, mchpr, mrs(*), wt(nqd,*), fit(*), mu(*), muwk(*),
     * v(nxis,*), vwk(nxis,*), cdnew(*), wtnew(nqd,*), fitnew(*), wk(*)
      integer i, j, k, kk, iter, flag, rkv, idamax, infowk
      double precision tmp, ddot, fitmean, dasum, lkhd, mumax, lkhdnew, 
     *disc, disc0, trc
      info = 0
      i=1
23000 if(.not.(i.le.nxis))goto 23002
      mrs(i) = 0.d0
      j=1
23003 if(.not.(j.le.nt))goto 23005
      if(cntsum.eq.0)then
      mrs(i) = mrs(i) + rs(i,j)
      else
      mrs(i) = mrs(i) + rs(i,j) * dfloat (cnt(j))
      endif
23004 j=j+1
      goto 23003
23005 continue
      mrs(i) = mrs(i) / dfloat (nobs)
23001 i=i+1
      goto 23000
23002 continue
      kk=1
23008 if(.not.(kk.le.nx))goto 23010
      i=1
23011 if(.not.(i.le.nqd))goto 23013
      wt(i,kk) = qdwt(i,kk) * dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1
     *))
23012 i=i+1
      goto 23011
23013 continue
23009 kk=kk+1
      goto 23008
23010 continue
      fitmean = 0.d0
      i=1
23014 if(.not.(i.le.nt))goto 23016
      tmp = ddot (nxis, rs(1,i), 1, cd, 1)
      fit(i) = dexp (tmp)
      if(cntsum.ne.0)then
      tmp = tmp * dfloat (cnt(i))
      endif
      fitmean = fitmean + tmp
23015 i=i+1
      goto 23014
23016 continue
      fitmean = fitmean / dfloat (nobs) - dasum (nqd*nx, wt, 1)
      call dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, wk, 1)
      lkhd = ddot (nxi, cd, 1, wk, 1) / 2.d0 - fitmean
      iter = 0
      flag = 0
23019 continue
      iter = iter + 1
      call dset (nxis, 0.d0, mu, 1)
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23022 if(.not.(kk.le.nx))goto 23024
      i=1
23025 if(.not.(i.le.nxis))goto 23027
      muwk(i) = - ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1)
      j=i
23028 if(.not.(j.le.nxis))goto 23030
      vwk(i,j) = 0.d0
      k=1
23031 if(.not.(k.le.nqd))goto 23033
      vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23032 k=k+1
      goto 23031
23033 continue
23029 j=j+1
      goto 23028
23030 continue
23026 i=i+1
      goto 23025
23027 continue
      call daxpy (nxis, 1.d0, muwk, 1, mu, 1)
      call daxpy (nxis*nxis, 1.d0, vwk, 1, v, 1)
23023 kk=kk+1
      goto 23022
23024 continue
      i=1
23034 if(.not.(i.le.nxi))goto 23036
      j=i
23037 if(.not.(j.le.nxi))goto 23039
      v(i,j) = v(i,j) + q(i,j)
23038 j=j+1
      goto 23037
23039 continue
23035 i=i+1
      goto 23034
23036 continue
      call daxpy (nxis, 1.d0, mrs, 1, mu, 1)
      call dsymv ('u', nxi, -1.d0, q, nxi, cd, 1, 1.d0, mu, 1)
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23040 if(.not.(i.le.nxis))goto 23042
      jpvt(i) = 0
23041 i=i+1
      goto 23040
23042 continue
      call dchdc (v, nxis, nxis, wk, jpvt, 1, rkv)
23043 if(v(rkv,rkv).lt.v(1,1)*dsqrt(mchpr))then
      rkv = rkv - 1
      goto 23043
      endif
23044 continue
      i=rkv+1
23045 if(.not.(i.le.nxis))goto 23047
      v(i,i) = v(1,1)
      call dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
23046 i=i+1
      goto 23045
23047 continue
23048 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dset (nxis-rkv, 0.d0, cdnew(rkv+1), 1)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      kk=1
23051 if(.not.(kk.le.nx))goto 23053
      i=1
23054 if(.not.(i.le.nqd))goto 23056
      tmp = ddot (nxis, qdrs(i,1,kk), nqd, cdnew, 1)
      if(tmp.gt.3.d2)then
      flag = flag + 1
      goto 23056
      endif
      wtnew(i,kk) = qdwt(i,kk) * dexp (tmp)
23055 i=i+1
      goto 23054
23056 continue
23052 kk=kk+1
      goto 23051
23053 continue
      fitmean = 0.d0
      i=1
23059 if(.not.(i.le.nt))goto 23061
      tmp = ddot (nxis, rs(1,i), 1, cdnew, 1)
      if(tmp.gt.3.d2)then
      flag = flag + 1
      goto 23061
      endif
      fitnew(i) = dexp (tmp)
      if(cntsum.ne.0)then
      tmp = tmp * dfloat (cnt(i))
      endif
      fitmean = fitmean + tmp
23060 i=i+1
      goto 23059
23061 continue
      fitmean = fitmean / dfloat (nobs) - dasum (nqd*nx, wtnew, 1)
      call dsymv ('u', nxi, 1.d0, q, nxi, cdnew, 1, 0.d0, wk, 1)
      lkhdnew = ddot (nxi, cdnew, 1, wk, 1) / 2.d0 - fitmean
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      call dcopy (nqd*nx, qdwt, 1, wt, 1)
      fitmean = - dasum (nqd*nx, wt, 1)
      lkhd = - fitmean
      iter = 0
      goto 23050
      endif
      if(flag.eq.3)then
      goto 23050
      endif
      if(lkhdnew-lkhd.lt.1.d1*(1.d0+dabs(lkhd))*mchpr)then
      goto 23050
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23050
      endif
23049 goto 23048
23050 continue
      if(flag.eq.1)then
      flag = 2
      goto 23020
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      kk=1
23078 if(.not.(kk.le.nx))goto 23080
      i=1
23081 if(.not.(i.le.nqd))goto 23083
      disc = dmax1 (disc, dabs(wt(i,kk)-wtnew(i,kk))/(1.d0+dabs(wt(i,kk)
     *)))
23082 i=i+1
      goto 23081
23083 continue
23079 kk=kk+1
      goto 23078
23080 continue
      i=1
23084 if(.not.(i.le.nt))goto 23086
      disc = dmax1 (disc, dabs(fit(i)-fitnew(i))/(1.d0+dabs(fit(i))))
23085 i=i+1
      goto 23084
23086 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
      disc0 = dmax1 ((mumax/(1.d0+lkhd))**2, dabs(lkhd-lkhdnew)/(1+dabs(
     *lkhd)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nqd*nx, wtnew, 1, wt, 1)
      call dcopy (nt, fitnew, 1, fit, 1)
      lkhd = lkhdnew
      if(disc0.lt.prec)then
      goto 23021
      endif
      if(disc.lt.prec)then
      goto 23021
      endif
      if(iter.lt.maxiter)then
      goto 23020
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      call dcopy (nqd*nx, qdwt, 1, wt, 1)
      fitmean = - dasum (nqd*nx, wt, 1)
      lkhd = - fitmean
      iter = 0
      flag = 2
      else
      info = 2
      goto 23021
      endif
23020 goto 23019
23021 continue
      infowk = 0
      i=1
23095 if(.not.(i.le.nt))goto 23097
      call daxpy (nxis, -1.d0, mrs, 1, rs(1,i), 1)
      call dprmut (rs(1,i), nxis, jpvt, 0)
      if(cntsum.ne.0)then
      call dscal (nxis, dsqrt(dfloat(cnt(i))), rs(1,i), 1)
      infowk = infowk + cnt(i)
      endif
      call dtrsl (v, nxis, nxis, rs(1,i), 11, infowk)
23096 i=i+1
      goto 23095
23097 continue
      call dprmut (mrs, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, mrs, 11, infowk)
      if(cntsum.ne.0)then
      tmp = dfloat (nobs-infowk)
      else
      tmp = dfloat (nobs-nt)
      endif
      trc = ddot (nxis*nt, rs, 1, rs, 1) + tmp * ddot (nxis, mrs, 1, mrs
     *, 1)
      trc = trc / dfloat(nobs) / (dfloat(nobs)-1.d0)
      mrs(1) = fitmean
      mrs(2) = trc
      kk=1
23102 if(.not.(kk.le.nx))goto 23104
      i=1
23105 if(.not.(i.le.nqd))goto 23107
      wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1))
23106 i=i+1
      goto 23105
23107 continue
23103 kk=kk+1
      goto 23102
23104 continue
      return
      end
