      subroutine dgold (vmu, q, ldq, n, z, low, upp, nlaht, score, 
     *varht, info, twk, work)
      character*1 vmu
      integer ldq, n, info
      double precision q(ldq,*), z(*), low, upp, nlaht, score, varht, 
     *twk(2,*), work(*)
      double precision ratio, mlo, mup, tmpl, tmpu
      ratio = ( dsqrt (5.d0) - 1.d0 ) / 2.d0
      info = 0
      if(.not.( upp .lt. low ))goto 23000
      mlo = low
      low = upp
      upp = mlo
23000 continue
      if(.not.( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' ))
     *goto 23002
      info = -3
      return
23002 continue
      if(.not.( n .lt. 1 .or. n .gt. ldq ))goto 23004
      info = -1
      return
23004 continue
      mlo = upp - ratio * (upp - low)
      call dset (n, 10.d0 ** (mlo), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**mlo
      call dtrev (vmu, twk, 2, n, z, tmpl, varht, info, work)
      if(.not.( info .ne. 0 ))goto 23006
      info = -2
      return
23006 continue
      mup = low + ratio * (upp - low)
      call dset (n, 10.d0 ** (mup), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**mup
      call dtrev (vmu, twk, 2, n, z, tmpu, varht, info, work)
      if(.not.( info .ne. 0 ))goto 23008
      info = -2
      return
23008 continue
23010 continue
      if(.not.( mup - mlo .lt. 1.d-7 ))goto 23013
      goto 23012
23013 continue
      if(.not.( tmpl .lt. tmpu ))goto 23015
      upp = mup
      mup = mlo
      tmpu = tmpl
      mlo = upp - ratio * (upp - low)
      call dset (n, 10.d0 ** (mlo), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**mlo
      call dtrev (vmu, twk, 2, n, z, tmpl, varht, info, work)
      if(.not.( info .ne. 0 ))goto 23017
      info = -2
      return
23017 continue
      goto 23016
23015 continue
      low = mlo
      mlo = mup
      tmpl = tmpu
      mup = low + ratio * (upp - low)
      call dset (n, 10.d0 ** (mup), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**mup
      call dtrev (vmu, twk, 2, n, z, tmpu, varht, info, work)
      if(.not.( info .ne. 0 ))goto 23019
      info = -2
      return
23019 continue
23016 continue
23011 goto 23010
23012 continue
      nlaht = ( mup + mlo ) / 2.d0
      call dset (n, 10.d0 ** (nlaht), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**nlaht
      call dtrev (vmu, twk, 2, n, z, score, varht, info, work)
      if(.not.( info .ne. 0 ))goto 23021
      info = -2
      return
23021 continue
      return
      end
