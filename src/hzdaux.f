C Output from Public domain Ratfor, version 1.0
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
