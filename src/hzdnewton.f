C Output from Public domain Ratfor, version 1.0
      subroutine hzdnewton (cd, nxis, q, nxi, rs, nt, nobs, cntsum, cnt,
     * qdrs, nqd, qdwt, nx, prec, maxiter, mchpr, wk, info)
      integer nxis, nxi, nt, nobs, cntsum, cnt(*), nqd, nx, maxiter, inf
     *o
      double precision cd(*), q(nxi,*), rs(nxis,*), qdrs(nqd,nxis,*), qd
     *wt(nqd,*), prec, mchpr, wk(*)
      integer imrs, iwt, ifit, imu, imuwk, iv, ivwk, ijpvt, icdnew, iwtn
     *ew, ifitnew, iwk
      imrs = 1
      iwt = imrs + max0 (nxis, 2)
      ifit = iwt + nqd*nx
      imu = ifit + nt
      imuwk = imu + nxis
      iv = imuwk + nxis
      ivwk = iv + nxis*nxis
      ijpvt = ivwk + nxis*nxis
      icdnew = ijpvt + nxis
      iwtnew = icdnew + nxis
      ifitnew = iwtnew + nqd*nx
      iwk = ifitnew + nt
      call hzdnewton1 (cd, nxis, q, nxi, rs, nt, nobs, cntsum, cnt, qdrs
     *, nqd, qdwt, nx, prec, maxiter, mchpr, wk(imrs), wk(iwt), wk(ifit)
     *, wk(imu), wk(imuwk), wk(iv), wk(ivwk), wk(ijpvt), wk(icdnew), wk(
     *iwtnew), wk(ifitnew), wk(iwk), info)
      return
      end
