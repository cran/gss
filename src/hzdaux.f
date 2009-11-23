C Output from Public domain Ratfor, version 1.01
      subroutine hzdaux1 (cd, nxis, q, nxi, qdrs, nqd, qdwt, nx, mchpr, 
     *wt, v, vwk, jpvt)
      integer nxis, nxi, nqd, nx, jpvt(*)
      double precision cd(*), q(nxi,*), qdrs(nqd,nxis,*), qdwt(nqd,*), m
     *chpr, wt(nqd,*), v(nxis,*), vwk(nxis,*)
      integer i, j, k, kk, rkv
      double precision ddot
      kk=1
23000 if(.not.(kk.le.nx))goto 23002
      i=1
23003 if(.not.(i.le.nqd))goto 23005
      wt(i,kk) = qdwt(i,kk) * dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1
     *))
23004 i=i+1
      goto 23003
23005 continue
23001 kk=kk+1
      goto 23000
23002 continue
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23006 if(.not.(kk.le.nx))goto 23008
      i=1
23009 if(.not.(i.le.nxis))goto 23011
      j=i
23012 if(.not.(j.le.nxis))goto 23014
      vwk(i,j) = 0.d0
      k=1
23015 if(.not.(k.le.nqd))goto 23017
      vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23016 k=k+1
      goto 23015
23017 continue
23013 j=j+1
      goto 23012
23014 continue
23010 i=i+1
      goto 23009
23011 continue
      call daxpy (nxis*nxis, 1.d0, vwk, 1, v, 1)
23007 kk=kk+1
      goto 23006
23008 continue
      i=1
23018 if(.not.(i.le.nxi))goto 23020
      j=i
23021 if(.not.(j.le.nxi))goto 23023
      v(i,j) = v(i,j) + q(i,j)
23022 j=j+1
      goto 23021
23023 continue
23019 i=i+1
      goto 23018
23020 continue
      i=1
23024 if(.not.(i.le.nxis))goto 23026
      jpvt(i) = 0
23025 i=i+1
      goto 23024
23026 continue
      call dchdc (v, nxis, nxis, vwk, jpvt, 1, rkv)
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
      return
      end
      subroutine hzdaux2 (v, nxis, jpvt, r, nr, se)
      double precision v(nxis,*), r(nxis,*), se(*)
      integer nxis, jpvt(*), nr
      double precision ddot
      integer i, infowk
      i=1
23032 if(.not.(i.le.nr))goto 23034
      call dprmut (r(1,i), nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, r(1,i), 11, infowk)
      se(i) = dsqrt (ddot (nxis, r(1,i), 1, r(1,i), 1))
23033 i=i+1
      goto 23032
23034 continue
      return
      end
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
23035 if(.not.(kk.le.nx))goto 23037
      i=1
23038 if(.not.(i.le.nxis))goto 23040
      mu0(i) = mu0(i) + ddot (nqd, wt0(1,kk), 1, qdrs(1,i,kk), 1)
23039 i=i+1
      goto 23038
23040 continue
23036 kk=kk+1
      goto 23035
23037 continue
      rkl = 0.d0
      kk=1
23041 if(.not.(kk.le.nx))goto 23043
      i=1
23044 if(.not.(i.le.nqd))goto 23046
      tmp = ddot (nxis, qdrs(i,1,kk), nqd, cd, 1)
      wt(i,kk) = qdwt(i,kk) * dexp (tmp)
      rkl = rkl + (wt(i,kk) - wt0(i,kk)*tmp)
23045 i=i+1
      goto 23044
23046 continue
23042 kk=kk+1
      goto 23041
23043 continue
      iter = 0
      flag = 0
23047 continue
      iter = iter + 1
      call dset (nxis, 0.d0, mu, 1)
      call dset (nxis*nxis, 0.d0, v, 1)
      kk=1
23050 if(.not.(kk.le.nx))goto 23052
      i=1
23053 if(.not.(i.le.nxis))goto 23055
      mu(i) = mu(i) - ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1)
      j=i
23056 if(.not.(j.le.nxis))goto 23058
      k=1
23059 if(.not.(k.le.nqd))goto 23061
      v(i,j) = v(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
23060 k=k+1
      goto 23059
23061 continue
23057 j=j+1
      goto 23056
23058 continue
23054 i=i+1
      goto 23053
23055 continue
23051 kk=kk+1
      goto 23050
23052 continue
      call daxpy (nxis, 1.d0, mu0, 1, mu, 1)
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23062 if(.not.(i.le.nxis))goto 23064
      jpvt(i) = 0
23063 i=i+1
      goto 23062
23064 continue
      call dmcdc (v, nxis, nxis, wk, jpvt, infowk)
23065 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      rklnew = 0.d0
      kk=1
23068 if(.not.(kk.le.nx))goto 23070
      i=1
23071 if(.not.(i.le.nqd))goto 23073
      tmp = ddot (nxis, qdrs(i,1,kk), nqd, cdnew, 1)
      if(tmp.gt.3.d2)then
      flag = flag + 1
      goto 23073
      endif
      wtnew(i,kk) = qdwt(i,kk) * dexp (tmp)
      rklnew = rklnew + (wtnew(i,kk) - wt0(i,kk)*tmp)
23072 i=i+1
      goto 23071
23073 continue
23069 kk=kk+1
      goto 23068
23070 continue
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      rkl = dasum (nqd*nx, qdwt, 1)
      call dcopy (nqd*nx, qdwt, 1, wt, 1)
      iter = 0
      goto 23067
      endif
      if(rklnew-rkl.lt.1.d1*(1.d0+dabs(rkl))*mchpr)then
      goto 23067
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23067
      endif
23066 goto 23065
23067 continue
      if(flag.eq.1)then
      flag = 2
      goto 23048
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      kk=1
23086 if(.not.(kk.le.nx))goto 23088
      i=1
23089 if(.not.(i.le.nqd))goto 23091
      disc = dmax1 (disc, dabs(wt(i,kk)-wtnew(i,kk))/(1.d0+dabs(wt(i,kk)
     *)))
23090 i=i+1
      goto 23089
23091 continue
23087 kk=kk+1
      goto 23086
23088 continue
      disc = dmax1 (disc, (mumax/(1.d0+rkl))**2)
      disc0 = dmax1 ((mumax/(1.d0+rkl))**2, dabs(rkl-rklnew)/(1+dabs(rkl
     *)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nqd*nx, wtnew, 1, wt, 1)
      rkl = rklnew
      if(disc0.lt.prec)then
      goto 23049
      endif
      if(disc.lt.prec)then
      goto 23049
      endif
      if(iter.lt.maxiter)then
      goto 23048
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      rkl = dasum (nqd*nx, qdwt, 1)
      call dcopy (nqd*nx, qdwt, 1, wt, 1)
      iter = 0
      flag = 2
      else
      info = 2
      goto 23049
      endif
23048 goto 23047
23049 continue
      kk=1
23100 if(.not.(kk.le.nx))goto 23102
      i=1
23103 if(.not.(i.le.nqd))goto 23105
      wt0(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1))
23104 i=i+1
      goto 23103
23105 continue
23101 kk=kk+1
      goto 23100
23102 continue
      return
      end
