C Output from Public domain Ratfor, version 1.0
      subroutine dmudr0 (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y,
     * tol, init, prec, maxite, theta, nlaht, score, varht, c, d, wk, in
     *fo)
      integer vmu
      integer lds, nobs, nnull, ldqr, ldqc, nq, init, maxite, info
      double precision s(lds,*), q(ldqr,ldqc,*), y(*), tol, prec, theta(
     **), nlaht, score, varht, c(*), d(*), wk(*)
      character vmu1
      if( vmu .eq. 1 )then
      vmu1 = 'v'
      endif
      if( vmu .eq. 2 )then
      vmu1 = 'm'
      endif
      if( vmu .eq. 3 )then
      vmu1 = 'u'
      endif
      call dmudr (vmu1, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y, tol, 
     *init, prec, maxite, theta, nlaht, score, varht, c, d, wk, info)
      return
      end
