C Output from Public domain Ratfor, version 1.0
      subroutine dtrev (vmu, t, ldt, n, z, score, varht, info, work)
      character vmu
      integer n, info
      double precision t(ldt,*), z(*), score, varht, work(*)
      double precision nume, deno, tmp, alph, la, dasum, ddot
      integer j
      info = 0
      if( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' )then
      info = -3
      return
      endif
      la = t(1,1)
      alph = dfloat (n) / dasum (n, t(2,1), ldt)
      call dscal (n, alph, t(2,1), ldt)
      call dscal (n-1, alph, t(1,2), ldt)
      call dpbfa (t, ldt, n, 1, info)
      if( info .ne. 0 )then
      return
      endif
      call dcopy (n, z, 1, work, 1)
      call dpbsl (t, ldt, n, 1, work)
      if( vmu .eq. 'v' )then
      tmp = 1.d0 / t(2,n) / t(2,n)
      deno = tmp
      j=n-1
23006 if(.not.(j.gt.0))goto 23008
      tmp = ( 1.d0 + t(1,j+1) * t(1,j+1) * tmp ) / t(2,j) / t(2,j)
      deno = deno + tmp
23007 j=j-1
      goto 23006
23008 continue
      nume = ddot (n, work, 1, work, 1) / dfloat (n)
      deno = deno / dfloat (n)
      varht = alph * la * nume / deno
      score = nume / deno / deno
      endif
      if( vmu .eq. 'm' )then
      deno = dlog (t(2,n))
      j=n-1
23011 if(.not.(j.gt.0))goto 23013
      deno = deno + dlog (t(2,j))
23012 j=j-1
      goto 23011
23013 continue
      nume = ddot (n, z, 1, work, 1) / dfloat (n)
      varht = alph * la * nume
      score = nume * dexp (2.d0 * deno / dfloat (n))
      endif
      if( vmu .eq. 'u' )then
      nume = ddot (n, work, 1, work, 1) / dfloat (n)
      tmp = 1.d0 / t(2,n) / t(2,n)
      deno = tmp
      j=n-1
23016 if(.not.(j.gt.0))goto 23018
      tmp = ( 1.d0 + t(1,j+1) * t(1,j+1) * tmp ) / t(2,j) / t(2,j)
      deno = deno + tmp
23017 j=j-1
      goto 23016
23018 continue
      deno = deno / dfloat (n)
      score = alph * alph * la * la * nume - 2.d0 * varht * alph * la * 
     *deno
      endif
      return
      end
