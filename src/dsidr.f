      subroutine dsidr (vmu, s, lds, nobs, nnull, y, q, ldq, tol, job, 
     *limnla, nlaht, score, varht, c, d, qraux, jpvt, wk, info)
      character*1 vmu
      integer lds, nobs, nnull, ldq, job, jpvt(*), info
      double precision s(lds,*), y(*), q(ldq,*), tol, limnla(2), nlaht, 
     *score(*), varht, c(*), d(*), qraux(*), wk(*)
      info = 0
      if(.not.( nnull .lt. 1 .or. nnull .ge. nobs .or. nobs .gt. lds 
     *.or. nobs .gt. ldq ))goto 23000
      info = -1
      return
23000 continue
      if(.not.( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' ))
     *goto 23002
      info = -3
      return
23002 continue
      call dstup (s, lds, nobs, nnull, qraux, jpvt, y, q, ldq, nobs, 1, 
     *info, wk)
      if(.not.( info .ne. 0 ))goto 23004
      return
23004 continue
      call dcore (vmu, q, ldq, nobs, nnull, tol, y, job, limnla, nlaht, 
     *score, varht, info, wk, wk(2*nobs+1))
      if(.not.( info .ne. 0 ))goto 23006
      return
23006 continue
      call dcoef (s, lds, nobs, nnull, qraux, jpvt, y, q, ldq, nlaht, c,
     * d, info, wk)
      return
      end
