C Output from Public domain Ratfor, version 1.0
      subroutine dnewton (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qdrs,
     * nqd, qdwt, prec, maxiter, mchpr, wk, info)
      integer nxis, nxi, nobs, cntsum, cnt(*), nqd, maxiter, info
      double precision cd(*), q(nxi,*), rs(nxis,*), qdrs(nqd,*), qdwt(*)
     *, prec, mchpr, wk(*)
      integer imrs, iwt, ifit, imu, iv, ijpvt, icdnew, iwtnew, ifitnew, 
     *iwk
      imrs = 1
      iwt = imrs + max0 (nxis, 3)
      ifit = iwt + nqd
      imu = ifit + nobs
      iv = imu + nxis
      ijpvt = iv + nxis*nxis
      icdnew = ijpvt + nxis
      iwtnew = icdnew + nxis
      ifitnew = iwtnew + nqd
      iwk = ifitnew + nobs
      call dnewton1 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qdrs, nqd,
     * qdwt, prec, maxiter, mchpr, wk(imrs), wk(iwt), wk(ifit), wk(imu),
     * wk(iv), wk(ijpvt), wk(icdnew), wk(iwtnew), wk(ifitnew), wk(iwk), 
     *info)
      return
      end
