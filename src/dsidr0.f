C Output from Public domain Ratfor, version 1.05
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
      info = 0
      if( nnull .lt. 1 .or. nnull .ge. nobs .or. nobs .gt. lds .or. nobs
     * .gt. ldq )then
      info = -1
      return
      endif
      call dstup (s, lds, nobs, nnull, qraux, jpvt, y, q, ldq, nobs, 1, 
     *info, wk)
      if( info .ne. 0 )then
      return
      endif
      call dcore (vmu1, q, ldq, nobs, nnull, tol, y, job, limnla, nlaht,
     * score, varht, info, wk, wk(2*nobs+1))
      if( info .ne. 0 )then
      return
      endif
      call dcoef (s, lds, nobs, nnull, qraux, jpvt, y, q, ldq, nlaht, c,
     * d, info, wk)
      return
      end
