C Output from Public domain Ratfor, version 1.0
      subroutine drkl (cd, nxis, qdrs, nqd, qdwt, wt0, mchpr, wt, eta, m
     *u, v, jpvt, wk, cdnew, wtnew, prec, maxiter, info)
      integer nxis, nqd, jpvt(*), maxiter, info
      double precision cd(*), qdrs(nqd,*), qdwt(*), wt0(*), mchpr, wt(*)
     *, eta(*), mu(*), v(nxis,*), wk(*), cdnew(*), wtnew(*), prec
      integer i, j, k, iter, flag, idamax, infowk
      double precision wtsum, ddot, rkl, mumax, wtsumnew, rklnew, disc, 
     *disc0
      info = 0
      wtsum = 0.d0
      rkl = 0.d0
      i=1
23000 if(.not.(i.le.nqd))goto 23002
      wt0(i) = qdwt(i) * wt0(i)
      eta(i) = ddot (nxis, qdrs(i,1), nqd, cd, 1)
      wt(i) = qdwt(i) * dexp (eta(i))
      wtsum = wtsum + wt(i)
      rkl = rkl + dlog(wt0(i)/wt(i)) * wt0(i)
23001 i=i+1
      goto 23000
23002 continue
      rkl = rkl + dlog(wtsum)
      iter = 0
      flag = 0
23003 continue
      iter = iter + 1
      i=1
23006 if(.not.(i.le.nxis))goto 23008
      mu(i) = ddot (nqd, wt, 1, qdrs(1,i), 1) / wtsum
23007 i=i+1
      goto 23006
23008 continue
      i=1
23009 if(.not.(i.le.nxis))goto 23011
      j=i
23012 if(.not.(j.le.nxis))goto 23014
      v(i,j) = 0.d0
      k=1
23015 if(.not.(k.le.nqd))goto 23017
      v(i,j) = v(i,j) + wt(k) * qdrs(k,i) * qdrs(k,j)
23016 k=k+1
      goto 23015
23017 continue
      v(i,j) = v(i,j) / wtsum - mu(i) * mu(j)
23013 j=j+1
      goto 23012
23014 continue
23010 i=i+1
      goto 23009
23011 continue
      i=1
23018 if(.not.(i.le.nxis))goto 23020
      mu(i) = ddot (nqd, wt0, 1, qdrs(1,i), 1) - mu(i)
23019 i=i+1
      goto 23018
23020 continue
      mumax = dabs(mu(idamax(nxis, mu, 1)))
      i=1
23021 if(.not.(i.le.nxis))goto 23023
      jpvt(i) = 0
23022 i=i+1
      goto 23021
23023 continue
      call dmcdc (v, nxis, nxis, wk, jpvt, infowk)
23024 continue
      call dcopy (nxis, mu, 1, cdnew, 1)
      call dprmut (cdnew, nxis, jpvt, 0)
      call dtrsl (v, nxis, nxis, cdnew, 11, infowk)
      call dtrsl (v, nxis, nxis, cdnew, 01, infowk)
      call dprmut (cdnew, nxis, jpvt, 1)
      call daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
      wtsumnew = 0.d0
      rklnew = 0.d0
      i=1
23027 if(.not.(i.le.nqd))goto 23029
      eta(i) = ddot (nxis, qdrs(i,1), nqd, cdnew, 1)
      if(eta(i).gt.3.d2)then
      flag = flag + 1
      goto 23029
      endif
      wtnew(i) = qdwt(i) * dexp (eta(i))
      wtsumnew = wtsumnew + wtnew(i)
      rklnew = rklnew + dlog(wt0(i)/wtnew(i)) * wt0(i)
23028 i=i+1
      goto 23027
23029 continue
      rklnew = rklnew + dlog (wtsumnew)
      if(flag.eq.1)then
      call dset (nxis, 0.d0, cd, 1)
      wtsum = 0.d0
      rkl = 0.d0
      i=1
23034 if(.not.(i.le.nqd))goto 23036
      wt(i) = qdwt(i)
      wtsum = wtsum + wt(i)
      rkl = rkl + dlog(wt0(i)/wt(i)) * wt0(i)
23035 i=i+1
      goto 23034
23036 continue
      rkl = rkl + dlog (wtsum)
      iter = 0
      goto 23026
      endif
      if(rklnew-rkl.lt.1.d1*(1.d0+dabs(rkl))*mchpr)then
      goto 23026
      endif
      call dscal (nxis, .5d0, mu, 1)
      if(dabs(mu(idamax(nxis, mu, 1))/mumax).lt.1.d1*mchpr)then
      goto 23026
      endif
23025 goto 23024
23026 continue
      if(flag.eq.1)then
      flag = 2
      goto 23004
      endif
      if(flag.eq.3)then
      info = 1
      return
      endif
      disc = 0.d0
      i=1
23045 if(.not.(i.le.nqd))goto 23047
      disc = dmax1 (disc, dabs(wt(i)-wtnew(i))/(1.d0+dabs(wt(i))))
23046 i=i+1
      goto 23045
23047 continue
      disc = dmax1 (disc, (mumax/(1.d0+rkl))**2)
      disc0 = dmax1 ((mumax/(1.d0+rkl))**2, dabs(rkl-rklnew)/(1+dabs(rkl
     *)))
      call dcopy (nxis, cdnew, 1, cd, 1)
      call dcopy (nqd, wtnew, 1, wt, 1)
      wtsum = wtsumnew
      rkl = rklnew
      if(disc0.lt.prec)then
      goto 23005
      endif
      if(disc.lt.prec)then
      goto 23005
      endif
      if(iter.lt.maxiter)then
      goto 23004
      endif
      if(flag.eq.0)then
      call dset (nxis, 0.d0, cd, 1)
      wtsum = 0.d0
      rkl = 0.d0
      i=1
23056 if(.not.(i.le.nqd))goto 23058
      wt(i) = qdwt(i)
      wtsum = wtsum + wt(i)
      rkl = rkl + dlog(wt0(i)/wt(i)) * wt0(i)
23057 i=i+1
      goto 23056
23058 continue
      rkl = rkl + dlog (wtsum)
      iter = 0
      flag = 2
      else
      info = 2
      goto 23005
      endif
23004 goto 23003
23005 continue
      wtsum = 0.d0
      i=1
23059 if(.not.(i.le.nqd))goto 23061
      wt0(i) = dexp (ddot (nxis, qdrs(i,1), nqd, cd, 1))
      wtsum = wtsum + qdwt(i) * wt0(i)
23060 i=i+1
      goto 23059
23061 continue
      call dscal (nqd, 1.d0/wtsum, wt0, 1)
      return
      end
