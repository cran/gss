      subroutine dmudr1 (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y,
     & tol, init, prec, maxite, theta, nlaht, score, varht, c, d, qraux,
     & jpvt, twk, traux, qwk, ywk, thewk, hes, gra, hwk1, hwk2, gwk1, 
     &gwk2, pvtwk, kwk, work1, work2, info)
      integer lds, nobs, nnull, ldqr, ldqc, nq, init, maxite, jpvt(*), 
     &pvtwk(*), info
      double precision s(lds,*), q(ldqr,ldqc,*), y(*), tol, prec, theta(
     &*), nlaht, score, varht, c(*), d(*), qraux(*), traux(*), twk(2,*),
     & qwk(nobs,*), ywk(*), thewk(*), hes(nq,*), gra(*), hwk1(nq,*), 
     &hwk2(nq,*), gwk1(*), gwk2(*), kwk(nobs-nnull,nobs-nnull,*), work1(
     &*), work2(*)
      character*1 vmu
      double precision alph, scrold, scrwk, nlawk, limnla(2), tmp, 
     &dasum, ddot
      integer n, n0, i, j, iwk, maxitwk, idamax, job
      info = 0
      n0 = nnull
      n = nobs - nnull
      maxitwk = maxite
      if(.not.( (vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u') 
     &.or. (init .ne. 0 .and. init .ne. 1) .or. (maxitwk .le.0) .or. (
     &prec .le. 0.d0) ))goto 23000
      info = -3
      return
23000 continue
      if(.not.( lds .lt. nobs .or. nobs .le. n0 .or. n0 .lt. 1 .or. 
     &ldqr .lt. nobs .or. ldqc .lt. nobs .or. nq .le. 0 ))goto 23002
      info = -1
      return
23002 continue
      call dstup (s, lds, nobs, n0, qraux, jpvt, y, q, ldqr, ldqc, nq, 
     &info, work1)
      if(.not.( info .ne. 0 ))goto 23004
      return
23004 continue
      if(.not.( init .eq. 1 ))goto 23006
      call dcopy (nq, theta, 1, thewk, 1)
      goto 23007
23006 continue
      i=1
23008 if(.not.(i.le.nq))goto 23010
      thewk(i) = dasum (n, q(n0+1,n0+1,i), ldqr+1)
      if(.not.( thewk(i) .gt. 0.d0 ))goto 23011
      thewk(i) = 1.d0 / thewk(i)
23011 continue
      i=i+1
      goto 23008
23010 continue
      j=1
23013 if(.not.(j.le.nobs))goto 23015
      call dset (nobs-j+1, 0.d0, qwk(j,j), 1)
      j=j+1
      goto 23013
23015 continue
      i=1
23016 if(.not.(i.le.nq))goto 23018
      j=1
23019 if(.not.(j.le.nobs))goto 23021
      call daxpy (nobs-j+1, thewk(i), q(j,j,i), 1, qwk(j,j), 1)
      j=j+1
      goto 23019
23021 continue
      i=i+1
      goto 23016
23018 continue
      call dcopy (nobs, y, 1, ywk, 1)
      call dcore (vmu, qwk, nobs, nobs, n0, tol, ywk, 0, limnla, nlawk, 
     &scrwk, varht, info, twk, work1)
      if(.not.(info .ne. 0 ))goto 23022
      return
23022 continue
      call dcoef (s, lds, nobs, n0, qraux, jpvt, ywk, qwk, nobs, nlawk, 
     &c, d, info, twk)
      call dqrsl (s, lds, nobs, n0, qraux, c, tmp, c, tmp, tmp, tmp, 
     &01000, info)
      i=1
23024 if(.not.(i.le.nq))goto 23026
      call dsymv('l', n, thewk(i), q(n0+1,n0+1,i), ldqr, c(n0+1), 1, 0.
     &d0, work1, 1)
      thewk(i) = ddot (n, c(n0+1), 1, work1, 1) * thewk(i)
      if(.not.( thewk(i) .gt. 0.d0 ))goto 23027
      thewk(i) = dlog10 (thewk(i))
      goto 23028
23027 continue
      thewk(i) = -25.d0
23028 continue
      i=i+1
      goto 23024
23026 continue
23007 continue
      scrold = 1.d10
      job = 0
23029 continue
      if(.not.( nq .eq. 1 ))goto 23032
      theta(1) = 0.d0
      goto 23031
23032 continue
      j=1
23034 if(.not.(j.le.nobs))goto 23036
      call dset (nobs-j+1, 0.d0, qwk(j,j), 1)
      j=j+1
      goto 23034
23036 continue
      i=1
23037 if(.not.(i.le.nq))goto 23039
      if(.not.( thewk(i) .le. -25.d0 ))goto 23040
      goto 23038
23040 continue
      j=1
23042 if(.not.(j.le.nobs))goto 23044
      call daxpy (nobs-j+1, 10.d0 ** thewk(i), q(j,j,i), 1, qwk(j,j), 1)
      j=j+1
      goto 23042
23044 continue
23038 i=i+1
      goto 23037
23039 continue
      call dcopy (nobs, y, 1, ywk, 1)
      call dcore (vmu, qwk, nobs, nobs, n0, tol, ywk, job, limnla, 
     &nlawk, scrwk, varht, info, twk, work1)
      if(.not.(info .ne. 0 ))goto 23045
      return
23045 continue
      if(.not.( scrold .lt. scrwk ))goto 23047
      tmp = dabs (gwk1(idamax (nq, gwk1, 1)))
      if(.not.( alph * tmp .gt. - prec ))goto 23049
      info = -5
      return
23049 continue
      alph = alph / 2.d0
      i=1
23051 if(.not.(i.le.nq))goto 23053
      thewk(i) = theta(i) + alph * gwk1(i)
      i=i+1
      goto 23051
23053 continue
      goto 23030
23047 continue
      maxitwk = maxitwk - 1
      call dcopy (n-2, qwk(n0+2,n0+1), nobs+1, traux, 1)
      call dcopy (n, qwk(n0+1,n0+1), nobs+1, twk(2,1), 2)
      call dcopy (n-1, qwk(n0+1,n0+2), nobs+1, twk(1,2), 2)
      call ddeev (vmu, nobs, q(n0+1,n0+1,1), ldqr, ldqc, n, nq, qwk(n0+
     &2,n0+1), nobs, traux, twk, ywk(n0+1), thewk, nlawk, scrwk, varht, 
     &hes, nq, gra, hwk1, hwk2, gwk1, gwk2, kwk, n, work1, work2, c, 
     &info)
      iwk = 0
      i=1
23054 if(.not.(i.le.nq))goto 23056
      if(.not.( thewk(i) .le. -25.d0 ))goto 23057
      goto 23055
23057 continue
      iwk = iwk + 1
      call dcopy (nq, hes(1,i), 1, hes(1,iwk), 1)
23055 i=i+1
      goto 23054
23056 continue
      iwk = 0
      i=1
23059 if(.not.(i.le.nq))goto 23061
      if(.not.( thewk(i) .le. -25.d0 ))goto 23062
      goto 23060
23062 continue
      iwk = iwk + 1
      call dcopy (nq, hes(i,1), nq, hes(iwk,1), nq)
      gwk1(iwk) = gra(i)
      work2(iwk) = gra(i)
23060 i=i+1
      goto 23059
23061 continue
      i=1
23064 if(.not.(i.lt.iwk))goto 23066
      call dcopy (iwk-i, hes(i+1,i), 1, hes(i,i+1), nq)
      i=i+1
      goto 23064
23066 continue
      call dmcdc (hes, nq, iwk, gwk2, pvtwk, info)
      call dprmut (gwk1, iwk, pvtwk, 0)
      call dposl (hes, nq, iwk, gwk1)
      call dprmut (gwk1, iwk, pvtwk, 1)
      alph = -1.d0
      j = iwk
      i=nq
23067 if(.not.(i.ge.1))goto 23069
      if(.not.( thewk(i) .le. -25.0 ))goto 23070
      gwk1(i) = 0.d0
      goto 23071
23070 continue
      gwk1(i) = gwk1(iwk)
      iwk = iwk - 1
23071 continue
      i=i-1
      goto 23067
23069 continue
      call dscal (nq, 1.d0/dlog(1.d1), gwk1, 1)
      tmp = dabs (gwk1(idamax (nq, gwk1, 1)))
      if(.not.( tmp .gt. 1.d0 ))goto 23072
      call dscal (nq, 1.d0/tmp, gwk1, 1)
23072 continue
      i=1
23074 if(.not.(i.le.nq))goto 23076
      if(.not.( thewk(i) .le. -25.d0 ))goto 23077
      goto 23075
23077 continue
      thewk(i) = thewk(i) - nlawk
23075 i=i+1
      goto 23074
23076 continue
      call dcopy (nq, thewk, 1, theta, 1)
      tmp = gra(idamax (nq, gra, 1)) ** 2
      if(.not.( tmp .lt. prec ** 2 .or. scrold - scrwk .lt. prec * (
     &scrwk + 1.d0) .and. tmp .lt. prec * (scrwk + 1.d0) ** 2 ))goto 230
     &79
      goto 23031
23079 continue
      if(.not.( maxitwk .lt. 1 ))goto 23081
      info = -4
      return
23081 continue
      scrold = scrwk
      i=1
23083 if(.not.(i.le.nq))goto 23085
      thewk(i) = thewk(i) + alph * gwk1(i)
      i=i+1
      goto 23083
23085 continue
      job = -1
      limnla(1) = -1.d0
      limnla(2) = 1.d0
23030 goto 23029
23031 continue
      j=1
23086 if(.not.(j.le.nobs))goto 23088
      call dset (nobs-j+1, 0.d0, qwk(j,j), 1)
      j=j+1
      goto 23086
23088 continue
      i=1
23089 if(.not.(i.le.nq))goto 23091
      if(.not.( theta(i) .le. -25.d0 ))goto 23092
      goto 23090
23092 continue
      j=1
23094 if(.not.(j.le.nobs))goto 23096
      call daxpy (nobs-j+1, 10.d0 ** theta(i), q(j,j,i), 1, qwk(j,j), 1)
      j=j+1
      goto 23094
23096 continue
23090 i=i+1
      goto 23089
23091 continue
      call dcopy (nobs, y, 1, ywk, 1)
      call dcore (vmu, qwk, nobs, nobs, n0, tol, ywk, job, limnla, 
     &nlaht, score, varht, info, twk, work1)
      if(.not.(info .ne. 0 ))goto 23097
      return
23097 continue
      call dcoef (s, lds, nobs, n0, qraux, jpvt, ywk, qwk, nobs, nlaht, 
     &c, d, info, twk)
      return
      end
