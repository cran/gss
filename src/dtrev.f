      subroutine dtrev (vmu, t, ldt, n, z, score, varht, info, work)
      character*1 vmu
      integer n, info
      double precision t(ldt,*), z(*), score, varht, work(*)
      double precision nume, deno, tmp, alph, la, dasum, ddot
      integer j
      info = 0
      if(.not.( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' ))
     *goto 23000
      info = -3
      return
23000 continue
      la = t(1,1)
      alph = dfloat (n) / dasum (n, t(2,1), ldt)
      call dscal (n, alph, t(2,1), ldt)
      call dscal (n-1, alph, t(1,2), ldt)
      call dpbfa (t, ldt, n, 1, info)
      if(.not.( info .ne. 0 ))goto 23002
      return
23002 continue
      call dcopy (n, z, 1, work, 1)
      call dpbsl (t, ldt, n, 1, work)
      if(.not.( vmu .eq. 'v' ))goto 23004
      tmp = 1.d0 / t(2,n) / t(2,n)
      deno = tmp
      j=n-1
23006 if(.not.(j.gt.0))goto 23008
      tmp = ( 1.d0 + t(1,j+1) * t(1,j+1) * tmp ) / t(2,j) / t(2,j)
      deno = deno + tmp
      j=j-1
      goto 23006
23008 continue
      nume = ddot (n, work, 1, work, 1) / dfloat (n)
      deno = deno / dfloat (n)
      varht = alph * la * nume / deno
      score = nume / deno / deno
23004 continue
      if(.not.( vmu .eq. 'm' ))goto 23009
      deno = dlog (t(2,n))
      j=n-1
23011 if(.not.(j.gt.0))goto 23013
      deno = deno + dlog (t(2,j))
      j=j-1
      goto 23011
23013 continue
      nume = ddot (n, z, 1, work, 1) / dfloat (n)
      varht = alph * la * nume
      score = nume * dexp (2.d0 * deno / dfloat (n))
23009 continue
      if(.not.( vmu .eq. 'u' ))goto 23014
      nume = ddot (n, work, 1, work, 1) / dfloat (n)
      tmp = 1.d0 / t(2,n) / t(2,n)
      deno = tmp
      j=n-1
23016 if(.not.(j.gt.0))goto 23018
      tmp = ( 1.d0 + t(1,j+1) * t(1,j+1) * tmp ) / t(2,j) / t(2,j)
      deno = deno + tmp
      j=j-1
      goto 23016
23018 continue
      deno = deno / dfloat (n)
      score = alph * alph * la * la * nume - 2.d0 * varht * alph * la * 
     *deno
23014 continue
      return
      end
