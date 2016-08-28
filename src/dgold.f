C Output from Public domain Ratfor, version 1.0
      subroutine dgold (vmu, q, ldq, n, z, low, upp, nlaht, score, varht
     *, info, twk, work)
      character vmu
      integer ldq, n, info
      double precision q(ldq,*), z(*), low, upp, nlaht, score, varht, tw
     *k(2,*), work(*)
      double precision ratio, mlo, mup, tmpl, tmpu
      ratio = ( dsqrt (5.d0) - 1.d0 ) / 2.d0
      info = 0
      if( upp .lt. low )then
      mlo = low
      low = upp
      upp = mlo
      endif
      if( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' )then
      info = -3
      return
      endif
      if( n .lt. 1 .or. n .gt. ldq )then
      info = -1
      return
      endif
      mlo = upp - ratio * (upp - low)
      call dset (n, 10.d0 ** (mlo), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**mlo
      call dtrev (vmu, twk, 2, n, z, tmpl, varht, info, work)
      if( info .ne. 0 )then
      info = -2
      return
      endif
      mup = low + ratio * (upp - low)
      call dset (n, 10.d0 ** (mup), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**mup
      call dtrev (vmu, twk, 2, n, z, tmpu, varht, info, work)
      if( info .ne. 0 )then
      info = -2
      return
      endif
23010 continue
      if( mup - mlo .lt. 1.d-7 )then
      goto 23012
      endif
      if( tmpl .lt. tmpu )then
      upp = mup
      mup = mlo
      tmpu = tmpl
      mlo = upp - ratio * (upp - low)
      call dset (n, 10.d0 ** (mlo), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**mlo
      call dtrev (vmu, twk, 2, n, z, tmpl, varht, info, work)
      if( info .ne. 0 )then
      info = -2
      return
      endif
      else
      low = mlo
      mlo = mup
      tmpl = tmpu
      mup = low + ratio * (upp - low)
      call dset (n, 10.d0 ** (mup), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**mup
      call dtrev (vmu, twk, 2, n, z, tmpu, varht, info, work)
      if( info .ne. 0 )then
      info = -2
      return
      endif
      endif
23011 goto 23010
23012 continue
      nlaht = ( mup + mlo ) / 2.d0
      call dset (n, 10.d0 ** (nlaht), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**nlaht
      call dtrev (vmu, twk, 2, n, z, score, varht, info, work)
      if( info .ne. 0 )then
      info = -2
      return
      endif
      return
      end
