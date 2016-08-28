C Output from Public domain Ratfor, version 1.0
      subroutine deval (vmu, q, ldq, n, z, nint, low, upp, nlaht, score,
     * varht, info, twk, work)
      character vmu
      integer ldq, n, nint, info
      double precision q(ldq,*), z(*), low, upp, nlaht, score(*), varht,
     * twk(2,*), work(*)
      double precision tmp, minscr, mlo, varhtwk
      integer j
      info = 0
      if( upp .lt. low )then
      mlo = low
      low = upp
      upp = mlo
      endif
      if( (vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u') .or. nint
     * .lt. 1 )then
      info = -3
      return
      endif
      if( 1 .gt. n .or. n .gt. ldq )then
      info = -1
      return
      endif
      j=1
23006 if(.not.(j.le.nint+1))goto 23008
      tmp = low + dfloat (j-1) * ( upp - low ) / dfloat (nint)
      call dset (n, 10.d0 ** (tmp), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**tmp
      call dtrev (vmu, twk, 2, n, z, score(j), varht, info, work)
      if( info .ne. 0 )then
      info = -2
      return
      endif
      if( score(j) .le. minscr .or. j .eq. 1 )then
      minscr = score(j)
      nlaht = tmp
      varhtwk = varht
      endif
23007 j=j+1
      goto 23006
23008 continue
      varht = varhtwk
      return
      end
