
#::::::::::::::
#   dnewton1
#::::::::::::::

subroutine  dnewton1 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt,
                      qdrs, nqd, qdwt, prec, maxiter, mchpr, 
                      mrs, wt, fit, mu, v, jpvt, cdnew, wtnew,
                      fitnew, wk, info)

integer  nxis, nxi, nobs, cntsum, cnt(*), nqd, maxiter, jpvt(*), info
double precision  cd(*), q(nxi,*), rs(nxis,*), qdrs(nqd,*), qdwt(*), prec, mchpr,
                  mrs(*), wt(*), fit(*), mu(*), v(nxis,*), cdnew(*), wtnew(*),
                  fitnew(*), wk(*)

integer  i, j, k, iter, flag, rkv, idamax, infowk
double precision  wtsum, tmp, ddot, fitmean, lkhd, mumax, wtsumnew, lkhdnew, disc, trc

#   Calculate constants
info = 0
for (i=1;i<=nxis;i=i+1) {
    mrs(i) = 0.d0
    if (cntsum==0) {
        for (j=1;j<=nobs;j=j+1)  mrs(i) = mrs(i) + rs(i,j)
        mrs(i) = mrs(i) / dfloat (nobs)
    }
    else {
        for (j=1;j<=nobs;j=j+1)  mrs(i) = mrs(i) + rs(i,j) * dfloat (cnt(j))
        mrs(i) = mrs(i) / dfloat (cntsum)
    }
}
#   Initialization
wtsum = 0.d0
for (i=1;i<=nqd;i=i+1) {
    wt(i) = qdwt(i) * dexp (ddot (nxis, qdrs(i,1), nqd, cd, 1))
    wtsum = wtsum + wt(i)
}
fitmean = 0.d0
for (i=1;i<=nobs;i=i+1) {
    tmp = ddot (nxis, rs(1,i), 1, cd, 1) - dlog (wtsum)
    fit(i) = dexp (tmp)
    if (cntsum!=0)  tmp = tmp * dfloat (cnt(i))
    fitmean = fitmean + tmp
}
if (cntsum==0)  fitmean = fitmean / dfloat (nobs)
else  fitmean = fitmean / dfloat (cntsum)
call  dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, wk, 1)
lkhd = ddot (nxi, cd, 1, wk, 1) / 2.d0 - fitmean
iter = 0
flag = 0
#   Iteration
repeat {
    iter = iter + 1
    # Calculate hessian and gradient
    for (i=1;i<=nxis;i=i+1)
        mu(i) = - ddot (nqd, wt, 1, qdrs(1,i), 1) / wtsum
    for (i=1;i<=nxis;i=i+1) {
        for (j=i;j<=nxis;j=j+1) {
            v(i,j) = 0.d0
            for (k=1;k<=nqd;k=k+1)
                v(i,j) = v(i,j) + wt(k) * qdrs(k,i) * qdrs(k,j)
            v(i,j) = v(i,j) / wtsum - mu(i) * mu(j)
            if (j<=nxi)  v(i,j) = v(i,j) + q(i,j)
        }
    }
    call  daxpy (nxis, 1.d0, mrs, 1, mu, 1)
    call  dsymv ('u', nxi, -1.d0, q, nxi, cd, 1, 1.d0, mu, 1)
    mumax = mu(idamax(nxis, mu, 1))
    #   Cholesky factorization
    for (i=1;i<=nxis;i=i+1)  jpvt(i) = 0
    call  dchdc (v, nxis, nxis, wk, jpvt, 1, rkv)
    while (v(rkv,rkv)<v(1,1)*dsqrt(mchpr))  rkv = rkv - 1
    for (i=rkv+1;i<=nxis;i=i+1) {
        v(i,i) = v(1,1)
        call  dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
    }
    #   Update coefficients
    repeat {
        call  dcopy (nxis, mu, 1, cdnew, 1)
        call  dprmut (cdnew, nxis, jpvt, 0)
        call  dtrsl (v, nxis, nxis, cdnew, 11, infowk)
        call  dset (nxis-rkv, 0.d0, cdnew(rkv+1), 1)
        call  dtrsl (v, nxis, nxis, cdnew, 01, infowk)
        call  dprmut (cdnew, nxis, jpvt, 1)
        call  daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
        wtsumnew = 0.d0
        for (i=1;i<=nqd;i=i+1) {
            wtnew(i) = qdwt(i) * dexp (ddot (nxis, qdrs(i,1), nqd, cdnew, 1))
            wtsumnew = wtsumnew + wtnew(i)
        }
        fitmean = 0.d0
        for (i=1;i<=nobs;i=i+1) {
            tmp = ddot (nxis, rs(1,i), 1, cdnew, 1) - dlog (wtsumnew)
            if (tmp>3.d2) {
                flag = flag + 1
                break
            }
            fitnew(i) = dexp (tmp)
            if (cntsum!=0)  tmp = tmp * dfloat (cnt(i))
            fitmean = fitmean + tmp
        }
        if (cntsum==0)  fitmean = fitmean / dfloat (nobs)
        else  fitmean = fitmean / dfloat (cntsum)
        if (flag==1) {
            #   Reset iteration with uniform starting value
            call  dset (nxis, 0.d0, cd, 1)
            wtsum = 0.d0
            for (i=1;i<=nqd;i=i+1) {
                wt(i) = qdwt(i)
                wtsum = wtsum + wt(i)
            }
            call  dset (nobs, 1.d0/wtsum, fit, 1)
            fitmean = - dlog (wtsum)
            lkhd = - fitmean
            iter = 0
            break
        }
        call  dsymv ('u', nxi, 1.d0, q, nxi, cdnew, 1, 0.d0, wk, 1)
        lkhdnew = ddot (nxi, cdnew, 1, wk, 1) / 2.d0 - fitmean
        if (lkhdnew-lkhd<2.d0*(1.d0+dabs(lkhd))*mchpr)  break
        call  dscal (nxis, .5d0, mu, 1)
    }
    if (flag==1) {
        flag = 2
        next
    }
    if (flag==3) {
        info = 1
        return
    }
    #   Calculate convergence criterion
    disc = 0.d0
    for (i=1;i<=nqd;i=i+1)
        disc = dmax1 (disc, dabs(wt(i)-wtnew(i))/(1.d0+dabs(wt(i))))
    for (i=1;i<=nobs;i=i+1)
        disc = dmax1 (disc, dabs(fit(i)-fitnew(i))/(1.d0+dabs(fit(i))))
    disc = dmax1 (disc, (mumax/(1.d0+lkhd))**2)
    #   Set to new values
    call  dcopy (nxis, cdnew, 1, cd, 1)
    call  dcopy (nqd, wtnew, 1, wt, 1)
    wtsum = wtsumnew
    call  dcopy (nobs, fitnew, 1, fit, 1)
    lkhd = lkhdnew
    #   Check convergence
    if (disc<prec)  break
    if (iter<maxiter)  next
    if (flag==0) {
        #   Reset iteration with uniform starting value
        call  dset (nxis, 0.d0, cd, 1)
        wtsum = 0.d0
        for (i=1;i<=nqd;i=i+1) {
            wt(i) = qdwt(i)
            wtsum = wtsum + wt(i)
        }
        call  dset (nobs, 1.d0/wtsum, fit, 1)
        fitmean = - dlog (wtsum)
        lkhd = - fitmean
        iter = 0
        flag = 2
    }
    else {
        info = 2
        break
    }
}
#   Calculate proxy loss
for (i=1;i<=nobs;i=i+1) {
    call  daxpy (nxis, -1.d0, mrs, 1, rs(1,i), 1)
    call  dprmut (rs(1,i), nxis, jpvt, 0)
    if (cntsum!=0)  call  dscal (nxis, dsqrt(dfloat(cnt(i))), rs(1,i), 1)
    call  dtrsl (v, nxis, nxis, rs(1,i), 11, infowk)
}
trc = ddot (nobs*nxis, rs, 1, rs, 1)
if (cntsum==0) {
    trc = trc / dfloat(nobs) / (dfloat(nobs)-1.d0)
    lkhd = 0.d0
    for (i=1;i<=nobs;i=i+1)  lkhd = lkhd + dlog (fit(i))
    lkhd = lkhd / dfloat (nobs)
}
else {
    trc = trc / dfloat(cntsum) / (dfloat(cntsum)-1.d0)
    lkhd = 0.d0
    for (i=1;i<=nobs;i=i+1)  lkhd = lkhd + dfloat (cnt(i)) * dlog (fit(i))
    lkhd = lkhd / dfloat (cntsum)
}
mrs(1) = lkhd
mrs(2) = trc
mrs(3) = wtsum

return
end
