C Output from Public domain Ratfor, version 1.0
      subroutine dsidr0 (vmu, s, lds, nobs, nnull, y, q, ldq, tol, job, 
     *limnla, nlaht, score, varht, c, d, qraux, jpvt, wk, info)
      integer vmu
      integer lds, nobs, nnull, ldq, job, jpvt(*), info
      double precision s(lds,*), y(*), q(ldq,*), tol, limnla(2), nlaht, 
     *score(*), varht, c(*), d(*), qraux(*), wk(*)
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
      call dsidr (vmu1, s, lds, nobs, nnull, y, q, ldq, tol, job, limnla
     *, nlaht, score, varht, c, d, qraux, jpvt, wk, info)
      return
      end
