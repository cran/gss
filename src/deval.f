      subroutine deval (vmu, q, ldq, n, z, nint, low, upp, nlaht, score,
     * varht, info, twk, work)
      character*1 vmu
      integer ldq, n, nint, info
      double precision q(ldq,*), z(*), low, upp, nlaht, score(*), varht,
     * twk(2,*), work(*)
      double precision tmp, minscr, mlo, varhtwk
      integer j
      info = 0
      if(.not.( upp .lt. low ))goto 23000
      mlo = low
      low = upp
      upp = mlo
23000 continue
      if(.not.( (vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u') 
     *.or. nint .lt. 1 ))goto 23002
      info = -3
      return
23002 continue
      if(.not.( 1 .gt. n .or. n .gt. ldq ))goto 23004
      info = -1
      return
23004 continue
      j=1
23006 if(.not.(j.le.nint+1))goto 23008
      tmp = low + dfloat (j-1) * ( upp - low ) / dfloat (nint)
      call dset (n, 10.d0 ** (tmp), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**tmp
      call dtrev (vmu, twk, 2, n, z, score(j), varht, info, work)
      if(.not.( info .ne. 0 ))goto 23009
      info = -2
      return
23009 continue
      if(.not.( score(j) .le. minscr .or. j .eq. 1 ))goto 23011
      minscr = score(j)
      nlaht = tmp
      varhtwk = varht
23011 continue
      j=j+1
      goto 23006
23008 continue
      varht = varhtwk
      return
      end
