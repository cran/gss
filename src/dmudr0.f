C Output from Public domain Ratfor, version 1.01
      subroutine dmudr0 (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y,
     * tol, init, prec, maxite, theta, nlaht, score, varht, c, d, iwk, w
     *k, info)
      integer vmu
      integer lds, nobs, nnull, ldqr, ldqc, nq, init, maxite, info, iwk(
     **)
      double precision s(lds,*), q(ldqr,ldqc,*), y(*), tol, prec, theta(
     **), nlaht, score, varht, c(*), d(*), wk(*)
      character*1 vmu1
      integer n, n0
      integer iqraux, itraux, itwk, iqwk, iywk, ithewk, ihes, igra, ihwk
     *1, ihwk2, igwk1, igwk2, ikwk, iwork1, iwork2, ijpvt, ipvtwk
      if( vmu .eq. 1 )then
      vmu1 = 'v'
      endif
      if( vmu .eq. 2 )then
      vmu1 = 'm'
      endif
      if( vmu .eq. 3 )then
      vmu1 = 'u'
      endif
      n = nobs
      n0 = nnull
      iqraux = 1
      itraux = iqraux + n0
      itwk = itraux + (n-n0-2)
      iqwk = itwk + 2 * (n-n0)
      iywk = iqwk + n * n
      ithewk = iywk + n
      ihes = ithewk + nq
      igra = ihes + nq * nq
      ihwk1 = igra + nq
      ihwk2 = ihwk1 + nq * nq
      igwk1 = ihwk2 + nq * nq
      igwk2 = igwk1 + nq
      ikwk = igwk2 + nq
      iwork1 = ikwk + (n-n0) * (n-n0) * nq
      iwork2 = iwork1 + n
      ijpvt = 1
      ipvtwk = ijpvt + n0
      call dmudr1 (vmu1, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y, tol,
     * init, prec, maxite, theta, nlaht, score, varht, c, d, wk(iqraux),
     * iwk(ijpvt), wk(itwk), wk(itraux), wk(iqwk), wk(iywk), wk(ithewk),
     * wk(ihes), wk(igra), wk(ihwk1), wk(ihwk2), wk(igwk1), wk(igwk2), i
     *wk(ipvtwk), wk(ikwk), wk(iwork1), wk(iwork2), info)
      return
      end
