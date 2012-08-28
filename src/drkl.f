C Output from Public domain Ratfor, version 1.01
      subroutine drkl (cd, nxis, qdrs, nqd, nt, bwt, qdwt, wt0, mchpr, p
     *rec, maxiter, jpvt, wk, info)
      integer nxis, nqd, nt, maxiter, jpvt(*), info
      double precision cd(*), qdrs(nqd,*), bwt(*), qdwt(nt,*), wt0(*), m
     *chpr, prec, wk(*)
      integer imrs, iwt, iwtsum, imu, imuwk, iv, ivwk, icdnew, iwtnew, i
     *wtnewsum
      imrs = 1
      iwt = imrs + nxis
      iwtsum = iwt + nt*nqd
      imu = iwtsum + nt
      imuwk = imu + nxis
      iv = imuwk + nxis
      ivwk = iv + nxis*nxis
      icdnew = ivwk + nxis*nxis
      iwtnew = icdnew + nxis
      iwtnewsum = iwtnew + nt*nqd
      call drkl1 (cd, nxis, qdrs, nqd, nt, bwt, qdwt, wt0, mchpr, prec, 
     *maxiter, wk(imrs), wk(iwt), wk(iwtsum), wk(imu), wk(imuwk), wk(iv)
     *, wk(ivwk), jpvt, wk(icdnew), wk(iwtnew), wk(iwtnewsum), info)
      return
      end
      subroutine drkl1 (cd, nxis, qdrs, nqd, nt, bwt, qdwt, wt0, mchpr, 
     *prec, maxiter, mrs, wt, wtsum, mu, muwk, v, vwk, jpvt, cdnew, wtne
     *w, wtnewsum, info)
      integer nxis, nqd, nt, maxiter, jpvt(*), info
      double precision cd(*), qdrs(nqd,*), bwt(*), qdwt(nt,*), wt0(*), m
     *chpr, mrs(*), wt(nt,*), wtsum(*), mu(*), muwk(*), v(nxis,*), vwk(n
     *xis,*), cdnew(*), wtnew(nt,*), wtnewsum(*), prec
      integer i, j, k, m, iter, flag, idamax, infowk
      double precision tmp, ddot, rkl, rklnew, mumax, disc, disc0
      info = 0
      call dset (nxis, 0.d0, mrs, 1)
      m=1
23000 if(.not.(m.le.nt))goto 23002
      i=1
23003 if(.not.(i.le.nqd))goto 23005
      wt(m,i) = qdwt(m,i) * wt0(i)
23004 i=i+1
      goto 23003
23005 continue
      i=1
23006 if(.not.(i.le.nxis))goto 23008
      muwk(i) = ddot (nqd, qdrs(1,i), 1, wt(m,1), nt)
23007 i=i+1
      goto 23006
23008 continue
      call daxpy (nxis, bwt(m), muwk, 1, mrs, 1)
      wtsum(m) = 0.d0
23001 m=m+1
      goto 23000
23002 continue
      i=1
23009 if(.not.(i.le.nqd))goto 23011
      tmp = dexp (ddot (nxis, qdrs(i,1), nqd, cd, 1))
      m=1
23012 if(.not.(m.le.nt))goto 23014
      wt(m,i) = qdwt(m,i) * tmp
      wtsum(m) = wtsum(m) + wt(m,i)
23013 m=m+1
      goto 23012
23014 continue
23010 i=i+1
      goto 23009
23011 continue
      rkl = 0.d0
      m=1
23015 if(.not.(m.le.nt))goto 23017
      tmp = 0.d0
      i=1
23018 if(.not.(i.le.nqd))goto 23020
      disc = wt0(i) * qdwt(m,i)
      tmp = tmp + dlog (disc/wt(m,i)) * disc
23019 i=i+1
      goto 23018
23020 continue
      rkl = rkl + bwt(m) * (tmp + dlog (wtsum(m)))
23016 m=m+1
      goto 23015
23017 continue
      iter = 0
      flag = 0
23021 continue
      iter = iter + 1
      call dset(nxis, 0.d0, mu, 1)
      call dset(nxis*nxis, 0.d0, v, 1)
      m=1
23024 if(.not.(m.le.nt))goto 23026
      i=1
23027 if(.not.(i.le.nxis))goto 23029
      muwk(i) = - ddot (nqd, wt(m,1), nt, qdrs(1,i), 1) / wtsum(m)
23028 i=i+1
      goto 23027
23029 continue
      i=1
23030 if(.not.(i.le.nxis))goto 23032
      j=i
23033 if(.not.(j.le.nxis))goto 23035
      vwk(i,j) = 0.d0
      k=1
23036 if(.not.(k.le.nqd))goto 23038
      vwk(i,j) = vwk(i,j) + wt(m,k) * qdrs(k,i) * qdrs(k,j)
23037 k=k+1
      goto 23036
23038 continue
      vwk(i,j) = vwk(i,j) / wtsum(m) - muwk(i) * muwk(j)
23034 j=j+1
      goto 23033
23035 continue
23031 i=i+1
      goto 23030
23032 continue
      call daxpy (nxis, bwt(m), muwk, 1, mu, 1)
      call daxpy (nxis*nxis, bwt(m), vwk, 1, v, 1)
23025 m=m+1
      goto 23024
23026 continue
      call daxpy (nxis, 1.d0, mrs, 1, mu, 1)
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23039 if(.not.(i.le.nxis))goto 23041
      jpvt(i) = 0
23040 i=i+1
      goto 23039
23041 continue
      call dmcdc (v, nxis, nxis, muwk, jpvt, infowk)
23042 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      call dset (nt, 0.d0, wtnewsum, 1)
      i=1
23045 if(.not.(i.le.nqd))goto 23047
      tmp = ddot (nxis, qdrs(i,1), nqd, cdnew, 1)
      if(tmp.gt.3.d2)then
      flag = flag + 1
      goto 23047
      endif
      m=1
23050 if(.not.(m.le.nt))goto 23052
      wtnew(m,i) = qdwt(m,i) * dexp (tmp)
      wtnewsum(m) = wtnewsum(m) + wtnew(m,i)
23051 m=m+1
      goto 23050
23052 continue
23046 i=i+1
      goto 23045
23047 continue
      rklnew = 0.d0
      m=1
23053 if(.not.(m.le.nt))goto 23055
      tmp = 0.d0
      i=1
23056 if(.not.(i.le.nqd))goto 23058
      disc = wt0(i) * qdwt(m,i)
      tmp = tmp + dlog (disc/wtnew(m,i)) * disc
23057 i=i+1
      goto 23056
23058 continue
      rklnew = rklnew + bwt(m) * (tmp + dlog (wtnewsum(m)))
23054 m=m+1
      goto 23053
23055 continue
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      call dcopy (nt*nqd, qdwt, 1, wt, 1)
      call dset (nt, 0.d0, wtsum, 1)
      m=1
23061 if(.not.(m.le.nt))goto 23063
      i=1
23064 if(.not.(i.le.nqd))goto 23066
      wtsum(m) = wtsum(m) + wt(m,i)
23065 i=i+1
      goto 23064
23066 continue
23062 m=m+1
      goto 23061
23063 continue
      rkl = 0.d0
      m=1
23067 if(.not.(m.le.nt))goto 23069
      tmp = 0.d0
      i=1
23070 if(.not.(i.le.nqd))goto 23072
      tmp = tmp + dlog (wt0(i)) * wt0(i) * qdwt(m,i)
23071 i=i+1
      goto 23070
23072 continue
      rkl = rkl + bwt(m) * (tmp + dlog (wtsum(m)))
23068 m=m+1
      goto 23067
23069 continue
      iter = 0
      goto 23044
      endif
      if(rklnew-rkl.lt.1.d1*(1.d0+dabs(rkl))*mchpr)then
      goto 23044
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23044
      endif
23043 goto 23042
23044 continue
      if(flag.eq.1)then
      flag = 2
      goto 23022
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      i=1
23081 if(.not.(i.le.nqd))goto 23083
      m=1
23084 if(.not.(m.le.nt))goto 23086
      disc = dmax1 (disc, dabs(wt(m,i)-wtnew(m,i))/(1.d0+dabs(wt(m,i))))
23085 m=m+1
      goto 23084
23086 continue
23082 i=i+1
      goto 23081
23083 continue
      disc = dmax1 (disc, (mumax/(1.d0+rkl))**2)
      disc = dmax1 (disc, dabs(rkl-rklnew)/(1.d0+dabs(rkl)))
      if(disc.lt.prec)then
      goto 23023
      endif
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nt*nqd, wtnew, 1, wt, 1)
      call dcopy (nt, wtnewsum, 1, wtsum, 1)
      rkl = rklnew
      if(iter.lt.maxiter)then
      goto 23022
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      call dcopy (nt*nqd, qdwt, 1, wt, 1)
      call dset (nt, 0.d0, wtsum, 1)
      m=1
23093 if(.not.(m.le.nt))goto 23095
      i=1
23096 if(.not.(i.le.nqd))goto 23098
      wtsum(m) = wtsum(m) + wt(m,i)
23097 i=i+1
      goto 23096
23098 continue
23094 m=m+1
      goto 23093
23095 continue
      rkl = 0.d0
      m=1
23099 if(.not.(m.le.nt))goto 23101
      tmp = 0.d0
      i=1
23102 if(.not.(i.le.nqd))goto 23104
      tmp = tmp + dlog (wt0(i)) * wt0(i) * qdwt(m,i)
23103 i=i+1
      goto 23102
23104 continue
      rkl = rkl + bwt(m) * (tmp + dlog (wtsum(m)))
23100 m=m+1
      goto 23099
23101 continue
      iter = 0
      flag = 2
      else
      info = 2
      goto 23023
      endif
23022 goto 23021
23023 continue
      i=1
23105 if(.not.(i.le.nqd))goto 23107
      wt0(i) = dexp (ddot (nxis, qdrs(i,1), nqd, cd, 1))
23106 i=i+1
      goto 23105
23107 continue
      return
      end
