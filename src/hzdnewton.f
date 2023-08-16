C Output from Public domain Ratfor, version 1.04
      subroutine hzdnewton (cd, nxis, q, nxi, rs, nt, nobs, cntsum, cnt,
     * qdrs, nqd, qdwt, nx, prec, maxiter, mchpr, jpvt, wk, info)
      integer nxis, nxi, nt, nobs, nqd, nx, maxiter, jpvt(*), info
      double precision cd(*), q(nxi,*), rs(nxis,*), cntsum, cnt(*), qdrs
     *(nqd,nxis,*), qdwt(nqd,*), prec, mchpr, wk(*)
      integer imrs, iwt, ifit, imu, imuwk, iv, ivwk, icdnew, iwtnew, ifi
     *tnew, iwk
      imrs = 1
      iwt = imrs + max0 (nxis, 2)
      ifit = iwt + nqd*nx
      imu = ifit + nt
      imuwk = imu + nxis
      iv = imuwk + nxis
      ivwk = iv + nxis*nxis
      icdnew = ivwk + nxis*nxis
      iwtnew = icdnew + nxis
      ifitnew = iwtnew + nqd*nx
      iwk = ifitnew + nt
      call hzdnewton1 (cd, nxis, q, nxi, rs, nt, nobs, cntsum, cnt, qdrs
     *, nqd, qdwt, nx, prec, maxiter, mchpr, wk(imrs), wk(iwt), wk(ifit)
     *, wk(imu), wk(imuwk), wk(iv), wk(ivwk), jpvt, wk(icdnew), wk(iwtne
     *w), wk(ifitnew), wk(iwk), info)
      return
      end
      subroutine hzdnewton1 (cd, nxis, q, nxi, rs, nt, nobs, cntsum, cnt
     *, qdrs, nqd, qdwt, nx, prec, maxiter, mchpr, mrs, wt, fit, mu, muw
     *k, v, vwk, jpvt, cdnew, wtnew, fitnew, wk, info)
      integer nxis, nxi, nt, nobs, nqd, nx, maxiter, jpvt(*), info
      double precision cd(*), q(nxi,*), rs(nxis,*), cntsum, cnt(*), qdrs
     *(nqd,nxis,*), qdwt(nqd,*), prec, mchpr, mrs(*), wt(nqd,*), fit(*),
     * mu(*), muwk(*), v(nxis,*), vwk(nxis,*), cdnew(*), wtnew(nqd,*), f
     *itnew(*), wk(*)
      integer i, j, k, kk, iter, flag, rkv, idamax, infowk
      double precision tmp, ddot, fitmean, dasum, lkhd, mumax, lkhdnew, 
     *disc, disc0, trc
      info = 0
      i=1
23000 if(.not.(i.le.nxis))goto 23002
      mrs(i) = 0.d0
      j=1
23003 if(.not.(j.le.nt))goto 23005
      if(.not.(cntsum.gt.0.d0))then
      mrs(i) = mrs(i) + rs(i,j)
      else
      mrs(i) = mrs(i) + rs(i,j) * cnt(j)
      endif
23004 j=j+1
      goto 23003
23005 continue
      mrs(i) = mrs(i) / dble (nobs)
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
      if(cntsum.gt.0.d0)then
      tmp = tmp * cnt(i)
      endif
      fitmean = fitmean + tmp
23015 i=i+1
      goto 23014
23016 continue
      fitmean = fitmean / dble (nobs) - dasum (nqd*nx, wt, 1)
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
      if((flag.eq.1).or.(flag.eq.3))then
      goto 23053
      endif
23052 kk=kk+1
      goto 23051
23053 continue
      if((flag.eq.0).or.(flag.eq.2))then
      fitmean = 0.d0
      i=1
23063 if(.not.(i.le.nt))goto 23065
      tmp = ddot (nxis, rs(1,i), 1, cdnew, 1)
      if(tmp.gt.3.d2)then
      flag = flag + 1
      goto 23065
      endif
      fitnew(i) = dexp (tmp)
      if(cntsum.gt.0.d0)then
      tmp = tmp * cnt(i)
      endif
      fitmean = fitmean + tmp
23064 i=i+1
      goto 23063
23065 continue
      fitmean = fitmean / dble (nobs) - dasum (nqd*nx, wtnew, 1)
      call dsymv ('u', nxi, 1.d0, q, nxi, cdnew, 1, 0.d0, wk, 1)
      lkhdnew = ddot (nxi, cdnew, 1, wk, 1) / 2.d0 - fitmean
      endif
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
23082 if(.not.(kk.le.nx))goto 23084
      i=1
23085 if(.not.(i.le.nqd))goto 23087
      disc = dmax1 (disc, dabs(wt(i,kk)-wtnew(i,kk))/(1.d0+dabs(wt(i,kk)
     *)))
23086 i=i+1
      goto 23085
23087 continue
23083 kk=kk+1
      goto 23082
23084 continue
      i=1
23088 if(.not.(i.le.nt))goto 23090
      disc = dmax1 (disc, dabs(fit(i)-fitnew(i))/(1.d0+dabs(fit(i))))
23089 i=i+1
      goto 23088
23090 continue
      disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
      disc0 = dmax1 ((mumax/(1.d0+lkhd))**2, dabs(lkhd-lkhdnew)/(1.d0+da
     *bs(lkhd)))
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
      i=1
23099 if(.not.(i.le.nt))goto 23101
      call dprmut (rs(1,i), nxis, jpvt, 0)
      if(cntsum.gt.0.d0)then
      call dscal (nxis, dsqrt(cnt(i)), rs(1,i), 1)
      endif
      call dtrsl (v, nxis, nxis, rs(1,i), 11, infowk)
23100 i=i+1
      goto 23099
23101 continue
      call dprmut (mrs, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, mrs, 11, infowk)
      trc = ddot (nxis*nt, rs, 1, rs, 1) - dble (nobs) * ddot (nxis, mrs
     *, 1, mrs, 1)
      trc = trc / dble(nobs) / (dble(nobs)-1.d0)
      mrs(1) = fitmean
      mrs(2) = trc
      kk=1
23104 if(.not.(kk.le.nx))goto 23106
      i=1
23107 if(.not.(i.le.nqd))goto 23109
      wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1))
23108 i=i+1
      goto 23107
23109 continue
23105 kk=kk+1
      goto 23104
23106 continue
      return
      end
