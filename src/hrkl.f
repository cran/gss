C Output from Public domain Ratfor, version 1.0
      subroutine hrkl (cd, nxis, qdrs, nqd, nx, qdwt, wt0, mchpr, wt, mu
     *, mu0, v, jpvt, wk, cdnew, wtnew, prec, maxiter, info)
      integer nxis, nqd, nx, jpvt(*), maxiter, info
      double precision cd(*), qdrs(nqd,nxis,*), qdwt(nqd,*), wt0(nqd,*),
     * mchpr, wt(nqd,*), mu(*), mu0(*), v(nxis,*), wk(*), cdnew(*), wtne
     *w(nqd,*), prec
      integer i, j, k, kk, idamax, iter, flag, infowk
      double precision tmp, ddot, dasum, rkl, mumax, rklnew, disc, disc0
      info = 0
      call dset (nxis, 0.d0, mu0, 1)
      kk=1
23000 if(.not.(kk.le.nx))goto 23002
      i=1
23003 if(.not.(i.le.nxis))goto 23005
      mu0(i) = mu0(i) + ddot (nqd, wt0(1,kk), 1, qdrs(1,i,kk), 1)
23004 i=i+1
      goto 23003
23005 continue
23001 kk=kk+1
      goto 23000
23002 continue
      rkl = 0.d0
      kk=1
23006 if(.not.(kk.le.nx))goto 23008
      i=1
23009 if(.not.(i.le.nqd))goto 23011
      tmp = ddot (nxis, qdrs(i,1,kk), nqd, cd, 1)
      wt(i,kk) = qdwt(i,kk) * dexp (tmp)
      rkl = rkl + (wt(i,kk) - wt0(i,kk)*tmp)
23010 i=i+1
      goto 23009
23011 continue
23007 kk=kk+1
      goto 23006
23008 continue
      iter = 0
      flag = 0
23012 continue
      iter = iter + 1
      call dset (nxis, 0.d0, mu, 1)
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23015 if(.not.(kk.le.nx))goto 23017
      i=1
23018 if(.not.(i.le.nxis))goto 23020
      mu(i) = mu(i) - ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1)
      j=i
23021 if(.not.(j.le.nxis))goto 23023
      k=1
23024 if(.not.(k.le.nqd))goto 23026
      v(i,j) = v(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23025 k=k+1
      goto 23024
23026 continue
23022 j=j+1
      goto 23021
23023 continue
23019 i=i+1
      goto 23018
23020 continue
23016 kk=kk+1
      goto 23015
23017 continue
      call daxpy (nxis, 1.d0, mu0, 1, mu, 1)
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23027 if(.not.(i.le.nxis))goto 23029
      jpvt(i) = 0
23028 i=i+1
      goto 23027
23029 continue
      call dmcdc (v, nxis, nxis, wk, jpvt, infowk)
23030 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      rklnew = 0.d0
      kk=1
23033 if(.not.(kk.le.nx))goto 23035
      i=1
23036 if(.not.(i.le.nqd))goto 23038
      tmp = ddot (nxis, qdrs(i,1,kk), nqd, cdnew, 1)
      if(tmp.gt.3.d2)then
      flag = flag + 1
      goto 23038
      endif
      wtnew(i,kk) = qdwt(i,kk) * dexp (tmp)
      rklnew = rklnew + (wtnew(i,kk) - wt0(i,kk)*tmp)
23037 i=i+1
      goto 23036
23038 continue
23034 kk=kk+1
      goto 23033
23035 continue
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      rkl = dasum (nqd*nx, qdwt, 1)
      call dcopy (nqd*nx, qdwt, 1, wt, 1)
      iter = 0
      goto 23032
      endif
      if(rklnew-rkl.lt.1.d1*(1.d0+dabs(rkl))*mchpr)then
      goto 23032
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23032
      endif
23031 goto 23030
23032 continue
      if(flag.eq.1)then
      flag = 2
      goto 23013
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      kk=1
23051 if(.not.(kk.le.nx))goto 23053
      i=1
23054 if(.not.(i.le.nqd))goto 23056
      disc = dmax1 (disc, dabs(wt(i,kk)-wtnew(i,kk))/(1.d0+dabs(wt(i,kk)
     *)))
23055 i=i+1
      goto 23054
23056 continue
23052 kk=kk+1
      goto 23051
23053 continue
      disc = dmax1 (disc, (mumax/(1.d0+rkl))**2)
      disc0 = dmax1 ((mumax/(1.d0+rkl))**2, dabs(rkl-rklnew)/(1+dabs(rkl
     *)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nqd*nx, wtnew, 1, wt, 1)
      rkl = rklnew
      if(disc0.lt.prec)then
      goto 23014
      endif
      if(disc.lt.prec)then
      goto 23014
      endif
      if(iter.lt.maxiter)then
      goto 23013
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      rkl = dasum (nqd*nx, qdwt, 1)
      call dcopy (nqd*nx, qdwt, 1, wt, 1)
      iter = 0
      flag = 2
      else
      info = 2
      goto 23014
      endif
23013 goto 23012
23014 continue
      kk=1
23065 if(.not.(kk.le.nx))goto 23067
      i=1
23068 if(.not.(i.le.nqd))goto 23070
      wt0(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1))
23069 i=i+1
      goto 23068
23070 continue
23066 kk=kk+1
      goto 23065
23067 continue
      return
      end
