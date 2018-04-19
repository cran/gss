
#::::::::::::::
#   dnewton9
#::::::::::::::

subroutine  dnewton9 (dc, nsr, sr, q, nobs, ez, qdsr, nqd, qdwt,
                      prec, maxiter, mchpr, jpvt, wk, info)

integer  nsr, nobs, nqd, maxiter, jpvt(*), info
double precision  dc(*), q(nsr-1,*), sr(nsr,*), ez(*), qdsr(nqd,*), qdwt(*),
                  prec, mchpr, wk(*)

integer  imsr, iwt, ifit, imu, iv, idcnew, iwtnew, ifitnew, iwk

imsr = 1
iwt = imsr + nsr
ifit = iwt + nqd
imu = ifit + nobs
iv = imu + nsr
idcnew = iv + nsr*nsr
iwtnew = idcnew + nsr
ifitnew = iwtnew + nqd
iwk = ifitnew + nobs

call  dnewton91 (dc, nsr, sr, q, nobs, ez, qdsr, nqd, qdwt,
                 prec, maxiter, mchpr, wk(imsr), wk(iwt),
                 wk(ifit), wk(imu), wk(iv), jpvt, wk(idcnew),
                 wk(iwtnew), wk(ifitnew), wk(iwk), info)

return
end


#:::::::::::::::
#   dnewton91
#:::::::::::::::

subroutine  dnewton91 (dc, nsr, sr, q, nobs, ez, qdsr, nqd, qdwt,
                       prec, maxiter, mchpr, msr, wt, fit, mu, v,
                       jpvt, dcnew, wtnew, fitnew, wk, info)

integer  nsr, nobs, nqd, maxiter, jpvt(*), info
double precision  dc(*), q(nsr-1,*), sr(nsr,*), ez(*), qdsr(nqd,*), qdwt(*),
                  prec, mchpr, msr(*), wt(*), fit(*), mu(*), v(nsr,*),
                  dcnew(*), wtnew(*), fitnew(*), wk(*)

integer  nr, i, j, k, iter, flag, rkv, idamax, infowk
double precision  dasum, tmp, ddot, ezsum, wtsum, fitmean, lkhd, mumax,
                  wtsumnew, lkhdnew, disc, disc0, trc

#   Calculate constants
info = 0
nr = nsr -1
ezsum = dasum (nobs, ez, 1)
call  dscal (nobs, 1.d0/ezsum, ez, 1)
call  dset (nsr, 0.d0, msr, 1)
for (i=1;i<=nobs;i=i+1)  call  daxpy (nsr, ez(i), sr(1,i), 1, msr, 1)
#   Initialization
for (i=1;i<=nqd;i=i+1) wt(i) = dexp (ddot (nsr, qdsr(i,1), nqd, dc, 1)) * qdwt(i)
wtsum = dasum (nqd, wt, 1)
fitmean = 0.d0
for (i=1;i<=nobs;i=i+1) {
    tmp = ddot (nsr, sr(1,i), 1, dc, 1)
    fit(i) = dexp (tmp)
    fitmean = fitmean + tmp * ez(i)
}
call  dsymv ('u', nr, 1.d0, q, nr, dc(2), 1, 0.d0, wk, 1)
lkhd = ddot (nr, dc(2), 1, wk, 1) / 2.d0 - fitmean + dlog (wtsum)
iter = 0
flag = 0
#   Iteration
repeat {
    iter = iter + 1
    # Calculate hessian and gradient
    for (i=1;i<=nsr;i=i+1)
        mu(i) = - ddot (nqd, wt, 1, qdsr(1,i), 1) / wtsum
    for (i=1;i<=nsr;i=i+1) {
        for (j=i;j<=nsr;j=j+1) {
            v(i,j) = 0.d0
            for (k=1;k<=nqd;k=k+1)
                v(i,j) = v(i,j) + wt(k) * qdsr(k,i) * qdsr(k,j)
            v(i,j) = v(i,j) / wtsum - mu(i) * mu(j)
        }
    }
    for (i=2;i<=nsr;i=i+1) {
        for (j=i;j<=nsr;j=j+1)  v(i,j) = v(i,j) + q(i-1,j-1)
    }
    call  daxpy (nsr, 1.d0, msr, 1, mu, 1)
    call  dsymv ('u', nr, -1.d0, q, nr, dc(2), 1, 1.d0, mu(2), 1)
    mumax = dabs(mu(idamax(nsr, mu, 1)))
    #   Cholesky factorization
    for (i=1;i<=nsr;i=i+1)  jpvt(i) = 0
    call  dchdc (v, nsr, nsr, wk, jpvt, 1, rkv)
    while (v(rkv,rkv)<v(1,1)*dsqrt(mchpr))  rkv = rkv - 1
    for (i=rkv+1;i<=nsr;i=i+1) {
        v(i,i) = v(1,1)
        call  dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
    }
    #   Update coefficients
    repeat {
        call  dcopy (nsr, mu, 1, dcnew, 1)
        call  dprmut (dcnew, nsr, jpvt, 0)
        call  dtrsl (v, nsr, nsr, dcnew, 11, infowk)
        call  dset (nsr-rkv, 0.d0, dcnew(rkv+1), 1)
        call  dtrsl (v, nsr, nsr, dcnew, 01, infowk)
        call  dprmut (dcnew, nsr, jpvt, 1)
        call  daxpy (nsr, 1.d0, dc, 1, dcnew, 1)
        for (i=1;i<=nqd;i=i+1) {
            tmp = ddot (nsr, qdsr(i,1), nqd, dcnew, 1)
            if (tmp>3.d2) {
                flag = flag + 1
                break
            }
            wtnew(i) = qdwt(i) * dexp (tmp)
        }
        wtsumnew = dasum (nqd, wtnew, 1)
        if ((flag==0)|(flag==2)) {
            fitmean = 0.d0
            for (i=1;i<=nobs;i=i+1) {
                tmp = ddot (nsr, sr(1,i), 1, dcnew, 1)
                if (tmp>3.d2) {
                    flag = flag + 1
                    break
                }
                fitnew(i) = dexp (tmp)
                fitmean = fitmean + tmp * ez(i)
            }
            call  dsymv ('u', nr, 1.d0, q, nr, dcnew(2), 1, 0.d0, wk, 1)
            lkhdnew = ddot (nr, dcnew(2), 1, wk, 1) / 2.d0 - fitmean + dlog (wtsumnew)
        }
        #   Reset iteration with uniform starting value
        if (flag==1) {
            call  dset (nsr, 0.d0, dc, 1)
            call  dcopy (nqd, qdwt, 1, wt, 1)
            wtsum = dasum (nqd, wt, 1)
            lkhd = dlog (wtsum)
            call  dset (nobs, 1.d0, fit, 1)
            iter = 0
            break
        }
        if (flag==3)  break
        if (lkhdnew-lkhd<1.d1*(1.d0+dabs(lkhd))*mchpr)  break
        call  dscal (nsr, .5d0, mu, 1)
        if (dabs(mu(idamax(nsr, mu, 1))/mumax)<1.d1*mchpr)  break
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
    disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
    disc0 = dmax1 ((mumax/(1.d0+lkhd))**2, dabs(lkhd-lkhdnew)/(1.d0+dabs(lkhd)))
    #   Set to new values
    call  dcopy (nsr, dcnew, 1, dc, 1)
    call  dcopy (nqd, wtnew, 1, wt, 1)
    wtsum = wtsumnew
    call  dcopy (nobs, fitnew, 1, fit, 1)
    lkhd = lkhdnew
    #   Check convergence
    if (disc0<prec)  break
    if (disc<prec)  break
    if (iter<maxiter)  next
    if (flag==0) {
        #   Reset iteration with uniform starting value
        call  dset (nsr, 0.d0, dc, 1)
        call  dcopy (nqd, qdwt, 1, wt, 1)
        wtsum = dasum (nqd, wt, 1)
        lkhd = dlog (wtsum)
        call  dset (nobs, 1.d0, fit, 1)
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
    call  daxpy (nsr, -1.d0, msr, 1, sr(1,i), 1)
    call  dprmut (sr(1,i), nsr, jpvt, 0)
    call  dscal (nsr, dsqrt(ez(i)), sr(1,i), 1)
    call  dtrsl (v, nsr, nsr, sr(1,i), 11, infowk)
    if (nsr-rkv>0)  call  dset (nsr-rkv, 0.d0, sr(rkv+1,i), 1)
}
trc = ddot (nobs*nsr, sr, 1, sr, 1)
lkhd = 0.d0
for (i=1;i<=nobs;i=i+1)  lkhd = lkhd + dlog (fit(i)) * ez(i)
lkhd = lkhd - dlog (wtsum)
msr(1) = lkhd
msr(2) = trc / dmax1(1.d0, ezsum-1.d0)

return
end
