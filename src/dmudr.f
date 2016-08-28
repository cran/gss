C Output from Public domain Ratfor, version 1.0
      subroutine dmudr (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y, 
     *tol, init, prec, maxite, theta, nlaht, score, varht, c, d, wk, inf
     *o)
      integer lds, nobs, nnull, ldqr, ldqc, nq, init, maxite, info
      double precision s(lds,*), q(ldqr,ldqc,*), y(*), tol, prec, theta(
     **), nlaht, score, varht, c(*), d(*), wk(*)
      character vmu
      integer n, n0
      integer iqraux, itraux, itwk, iqwk, iywk, ithewk, ihes, igra, ihwk
     *1, ihwk2, igwk1, igwk2, ikwk, iwork1, iwork2, ijpvt, ipvtwk
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
      ijpvt = iwork2 + n
      ipvtwk = ijpvt + n0
      call dmudr1 (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y, tol, 
     *init, prec, maxite, theta, nlaht, score, varht, c, d, wk(iqraux), 
     *wk(ijpvt), wk(itwk), wk(itraux), wk(iqwk), wk(iywk), wk(ithewk), w
     *k(ihes), wk(igra), wk(ihwk1), wk(ihwk2), wk(igwk1), wk(igwk2), wk(
     *ipvtwk), wk(ikwk), wk(iwork1), wk(iwork2), info)
      return
      end
