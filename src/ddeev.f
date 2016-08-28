C Output from Public domain Ratfor, version 1.0
      subroutine ddeev (vmu, nobs, q, ldqr, ldqc, n, nq, u, ldu, uaux, t
     *, x, theta, nlaht, score, varht, hes, ldh, gra, hwk1, hwk2, gwk1, 
     *gwk2, kwk, ldk, work1, work2, work3, info)
      character vmu
      integer nobs, ldqr, ldqc, n, nq, ldu, ldh, ldk, info
      double precision q(ldqr,ldqc,*), u(ldu,*), uaux(*), t(2,*), x(*), 
     *theta(*), nlaht, score, varht, hes(ldh,*), gra(*), hwk1(nq,*), hwk
     *2(nq,*), gwk1(*), gwk2(*), kwk(ldk,ldk,*), work1(*), work2(*), wor
     *k3(*)
      double precision trc, det, dum, ddot
      integer i, j, m
      info = 0
      call dset (nq, 0.d0, gra, 1)
      call dset (nq*nq, 0.d0, hes, 1)
      if( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' )then
      info = -3
      return
      endif
      if( nobs .lt. n .or. ldqr .lt. n .or. ldqc .lt. n .or. nq .le. 0 .
     *or. ldu .lt. n-1 .or. ldh .lt. nq .or. ldk .lt. n )then
      info = -1
      return
      endif
      i=2
23004 if(.not.(i.le.nq))goto 23006
      if( theta(i) .le. -25.d0 )then
      goto 23005
      endif
      j=1
23009 if(.not.(j.le.n))goto 23011
      call dcopy (n-j+1, q(j,j,i), 1, kwk(j,j,i), 1)
      call dscal (n-j+1, 10.d0 ** theta(i), kwk(j,j,i), 1)
23010 j=j+1
      goto 23009
23011 continue
      call dqrslm (u, ldu, n-1, n-2, uaux, kwk(2,2,i), n, 0, info, work1
     *)
      call dqrsl (u, ldu, n-1, n-2, uaux, kwk(2,1,i), dum, kwk(2,1,i), d
     *um, dum, dum, 01000, info)
23005 i=i+1
      goto 23004
23006 continue
      call dcopy (n, t(2,1), 2, kwk(1,1,1), n+1)
      call dcopy (n-1, t(1,2), 2, kwk(2,1,1), n+1)
      j=1
23012 if(.not.(j.lt.n-1))goto 23014
      call dset (n-j-1, 0.d0, kwk(j+2,j,1), 1)
23013 j=j+1
      goto 23012
23014 continue
      i=2
23015 if(.not.(i.le.nq))goto 23017
      if( theta(i) .le. -25.d0 )then
      goto 23016
      endif
      j=1
23020 if(.not.(j.le.n))goto 23022
      call daxpy (n-j+1, -1.d0, kwk(j,j,i), 1, kwk(j,j,1), 1)
23021 j=j+1
      goto 23020
23022 continue
23016 i=i+1
      goto 23015
23017 continue
      i=1
23023 if(.not.(i.le.nq))goto 23025
      if( theta(i) .le. -25.d0 )then
      goto 23024
      endif
      j=1
23028 if(.not.(j.lt.n))goto 23030
      call dcopy (n-j, kwk(j+1,j,i), 1, kwk(j,j+1,i), n)
23029 j=j+1
      goto 23028
23030 continue
23024 i=i+1
      goto 23023
23025 continue
      call dset (n, 10.d0 ** nlaht, work1, 1)
      call daxpy (n, 1.d0, work1, 1, t(2,1), 2)
      call dpbfa (t, 2, n, 1, info)
      if( info .ne. 0 )then
      info = -2
      return
      endif
      i=1
23033 if(.not.(i.le.nq))goto 23035
      if( theta(i) .le. -25.d0 )then
      goto 23034
      endif
      j=1
23038 if(.not.(j.le.n))goto 23040
      call dpbsl (t, 2, n, 1, kwk(1,j,i))
23039 j=j+1
      goto 23038
23040 continue
23034 i=i+1
      goto 23033
23035 continue
      call dcopy (n, x, 1, work1, 1)
      call dpbsl (t, 2, n, 1, work1)
      if( vmu .ne. 'm' )then
      call dcopy (n, work1, 1, work2, 1)
      call dscal (n, 2.d0, work2, 1)
      else
      call dcopy (n, x, 1, work2, 1)
      endif
      i=1
23043 if(.not.(i.le.nq))goto 23045
      if( theta(i) .le. -25.d0 )then
      goto 23044
      endif
      call dgemv ('t', n, n, 1.d0, kwk(1,1,i), n, work2, 1, 0.d0, work3,
     * 1)
      gwk1(i) = - ddot (n, work1, 1, work3, 1)
23044 i=i+1
      goto 23043
23045 continue
      i=1
23048 if(.not.(i.le.nq))goto 23050
      gwk2(i) = 0.d0
      if( theta(i) .le. -25.d0 )then
      goto 23049
      endif
      j=1
23053 if(.not.(j.le.n))goto 23055
      if( vmu .ne. 'm' )then
      call dcopy (n, kwk(1,j,i), 1, work1, 1)
      call dpbsl (t, 2, n, 1, work1)
      gwk2(i) = gwk2(i) - work1(j)
      else
      gwk2(i) = gwk2(i) - kwk(j,j,i)
      endif
23054 j=j+1
      goto 23053
23055 continue
23049 i=i+1
      goto 23048
23050 continue
      if( vmu .ne. 'm' )then
      call dcopy (n, x, 1, work1, 1)
      call dpbsl (t, 2, n, 1, work1)
      i=1
23060 if(.not.(i.le.nq))goto 23062
      if( theta(i) .le. -25.d0 )then
      goto 23061
      endif
      call dgemv ('n', n, n, 1.d0, kwk(1,1,i), n, work1, 1, 0.d0, work2,
     * 1)
      j=1
23065 if(.not.(j.le.i))goto 23067
      if( theta(j) .le. -25.d0 )then
      goto 23066
      endif
      call dgemv ('n', n, n, 1.d0, kwk(1,1,j), n, work1, 1, 0.d0, work3,
     * 1)
      hwk1(i,j) = 2.d0 * ddot (n, work2, 1, work3, 1)
      call dgemv ('t', n, n, 1.d0, kwk(1,1,j), n, work1, 1, 0.d0, work3,
     * 1)
      hwk1(i,j) = hwk1(i,j) + 2.d0 * ddot (n, work2, 1, work3, 1)
23066 j=j+1
      goto 23065
23067 continue
      call dgemv ('t', n, n, 1.d0, kwk(1,1,i), n, work1, 1, 0.d0, work2,
     * 1)
      j=1
23070 if(.not.(j.le.i))goto 23072
      if( theta(j) .le. -25.d0 )then
      goto 23071
      endif
      call dgemv ('n', n, n, 1.d0, kwk(1,1,j), n, work1, 1, 0.d0, work3,
     * 1)
      hwk1(i,j) = hwk1(i,j) + 2.d0 * ddot (n, work2, 1, work3, 1)
23071 j=j+1
      goto 23070
23072 continue
23061 i=i+1
      goto 23060
23062 continue
      else
      call dcopy (n, x, 1, work1, 1)
      call dpbsl (t, 2, n, 1, work1)
      i=1
23075 if(.not.(i.le.nq))goto 23077
      if( theta(i) .le. -25.d0 )then
      goto 23076
      endif
      call dgemv ('n', n, n, 1.d0, kwk(1,1,i), n, work1, 1, 0.d0, work2,
     * 1)
      j=1
23080 if(.not.(j.le.i))goto 23082
      if( theta(j) .le. -25.d0 )then
      goto 23081
      endif
      call dgemv ('t', n, n, 1.d0, kwk(1,1,j), n, x, 1, 0.d0, work3, 1)
      hwk1(i,j) = 2.d0 * ddot (n, work2, 1, work3, 1)
23081 j=j+1
      goto 23080
23082 continue
23076 i=i+1
      goto 23075
23077 continue
      endif
      i=1
23085 if(.not.(i.le.nq))goto 23087
      if( theta(i) .le. -25.d0 )then
      goto 23086
      endif
      hwk1(i,i) = hwk1(i,i) + gwk1(i)
23086 i=i+1
      goto 23085
23087 continue
      i=1
23090 if(.not.(i.le.nq))goto 23092
      if( theta(i) .le. -25.d0 )then
      goto 23091
      endif
      m=1
23095 if(.not.(m.le.i))goto 23097
      hwk2(i,m) = 0.d0
      if( theta(m) .le. -25.d0 )then
      goto 23096
      endif
      j=1
23100 if(.not.(j.le.n))goto 23102
      if( vmu .ne. 'm' )then
      call dcopy (n, kwk(1,j,m), 1, work1, 1)
      call dpbsl (t, 2, n, 1, work1)
      hwk2(i,m) = hwk2(i,m) + 2.d0 * ddot (n, kwk(j,1,i), n, work1, 1)
      else
      hwk2(i,m) = hwk2(i,m) + ddot (n, kwk(j,1,i), n, kwk(1,j,m), 1)
      endif
23101 j=j+1
      goto 23100
23102 continue
23096 m=m+1
      goto 23095
23097 continue
23091 i=i+1
      goto 23090
23092 continue
      i=1
23105 if(.not.(i.le.nq))goto 23107
      if( theta(i) .le. -25.d0 )then
      goto 23106
      endif
      hwk2(i,i) = hwk2(i,i) + gwk2(i)
23106 i=i+1
      goto 23105
23107 continue
      if( vmu .eq. 'v' )then
      trc = dfloat (nobs) * 10.d0 ** (-nlaht) * varht / score
      i=1
23112 if(.not.(i.le.nq))goto 23114
      if( theta(i) .le. -25.d0 )then
      goto 23113
      endif
      gra(i) = gwk1(i) / trc / trc - 2.d0 * score * gwk2(i) / trc / dflo
     *at(nobs)
23113 i=i+1
      goto 23112
23114 continue
      call dscal (nq, dfloat (nobs), gra, 1)
      endif
      if( vmu .eq. 'u' )then
      dum = 10.d0 ** nlaht
      i=1
23119 if(.not.(i.le.nq))goto 23121
      if( theta(i) .le. -25.d0 )then
      goto 23120
      endif
      gra(i) = dum * dum * gwk1(i) - 2.d0 * varht * dum * gwk2(i)
23120 i=i+1
      goto 23119
23121 continue
      call dscal (nq, 1.d0/dfloat (n), gra, 1)
      endif
      if( vmu .eq. 'm' )then
      det = 10.d0 ** (-nlaht) * varht / score
      i=1
23126 if(.not.(i.le.nq))goto 23128
      if( theta(i) .le. -25.d0 )then
      goto 23127
      endif
      gra(i) = gwk1(i) / det - dfloat (nobs) / dfloat (n) * score * gwk2
     *(i)
23127 i=i+1
      goto 23126
23128 continue
      call dscal (nq, 1.d0 / dfloat (nobs), gra, 1)
      endif
      if( vmu .eq. 'v' )then
      i=1
23133 if(.not.(i.le.nq))goto 23135
      if( theta(i) .le. -25.d0 )then
      goto 23134
      endif
      j=1
23138 if(.not.(j.le.i))goto 23140
      if( theta(j) .le. -25.d0 )then
      goto 23139
      endif
      hes(i,j) = hwk1(i,j) / trc / trc - 2.d0 * gwk1(i) * gwk2(j) / trc 
     *** 3 - 2.d0 * gwk1(j) * gwk2(i) / trc ** 3 - 2.d0 * score * hwk2(i
     *,j) / trc / dfloat (nobs) + 6.d0 * score * gwk2(i) * gwk2(j) / trc
     * / trc / dfloat (nobs)
23139 j=j+1
      goto 23138
23140 continue
      call dscal (i, dfloat (nobs), hes(i,1), ldh)
23134 i=i+1
      goto 23133
23135 continue
      endif
      if( vmu .eq. 'u' )then
      i=1
23145 if(.not.(i.le.nq))goto 23147
      if( theta(i) .le. -25.d0 )then
      goto 23146
      endif
      j=1
23150 if(.not.(j.le.i))goto 23152
      if( theta(j) .le. -25.d0 )then
      goto 23151
      endif
      hes(i,j) = dum * dum * hwk1(i,j) - 2.d0 * varht * dum * hwk2(i,j)
23151 j=j+1
      goto 23150
23152 continue
      call dscal (i, 1.d0/dfloat (n), hes(i,1), ldh)
23146 i=i+1
      goto 23145
23147 continue
      endif
      if( vmu .eq. 'm' )then
      i=1
23157 if(.not.(i.le.nq))goto 23159
      if( theta(i) .le. -25.d0 )then
      goto 23158
      endif
      j=1
23162 if(.not.(j.le.i))goto 23164
      if( theta(j) .le. -25.d0 )then
      goto 23163
      endif
      hes(i,j) = hwk1(i,j) / det - gwk1(i) * gwk2(j) / det / dfloat (n) 
     *- gwk1(j) * gwk2(i) / det / dfloat (n) - dfloat (nobs) / dfloat (n
     *) * score * hwk2(i,j) + dfloat (nobs) / dfloat (n) ** 2 * score * 
     *gwk2(i) * gwk2(j)
23163 j=j+1
      goto 23162
23164 continue
      call dscal (i, 1.d0 / dfloat (nobs), hes(i,1), ldh)
23158 i=i+1
      goto 23157
23159 continue
      endif
      return
      end
