
#:::::::::::::::
#   dnewton99
#:::::::::::::::

subroutine  dnewton99 (dc, nsr, sr, q, nobs, cntsum, cnt, ff, qdsr, nqd,
                       qdwt, prec, maxiter, mchpr, jpvt, wk, info)

integer  nsr, nobs, cntsum, cnt(*), nqd, maxiter, jpvt(*), info
double precision  dc(*), q(nsr-1,*), sr(nsr,*), ff(nobs,*), qdsr(nqd,*), qdwt(*),
                  prec, mchpr, wk(*)

integer  iwt, ifit, imu, imuwk, iv, ivwk, idcnew, iwtnew, ifitnew, iwk

iwt = 1
ifit = iwt + nqd
imu = ifit + nobs
imuwk = imu + nsr
iv = imuwk + nsr
ivwk = iv + nsr*nsr
idcnew = ivwk + nsr*nsr
iwtnew = idcnew + nsr
ifitnew = iwtnew + nqd
iwk = ifitnew + nobs

call  dnewton991 (dc, nsr, sr, q, nobs, cntsum, cnt, ff, qdsr, nqd,
                  qdwt, prec, maxiter, mchpr, wk(iwt), wk(ifit), wk(imu),
                  wk(imuwk), wk(iv), wk(ivwk), jpvt, wk(idcnew), wk(iwtnew),
                  wk(ifitnew), wk(iwk), info)

return
end


#::::::::::::::::
#   dnewton991
#::::::::::::::::

subroutine  dnewton991 (dc, nsr, sr, q, nobs, cntsum, cnt, ff, qdsr, nqd,
                        qdwt, prec, maxiter, mchpr, wt, fit, mu, muwk,
                        v, vwk, jpvt, dcnew, wtnew, fitnew, wk, info)

integer  nsr, nobs, cntsum, cnt(*), nqd, maxiter, jpvt(*), info
double precision  dc(*), q(nsr-1,*), sr(nsr,*), ff(nobs,*), qdsr(nqd,*), qdwt(*),
                  prec, mchpr, wt(*), fit(*), mu(*), muwk(*), v(nsr,*),
                  vwk(nsr,*), dcnew(*), wtnew(*), fitnew(*), wk(*)

integer  nr, i, j, k, iter, flag, idamax, infowk
double precision  dasum, tmp, ddot, wtsum, lkhd, mumax, wtsumnew, lkhdnew, disc, disc0

#   Calculate constants
info = 0
nr = nsr - 1
#   Initialization
for (i=1;i<=nqd;i=i+1) wt(i) = dexp (ddot (nsr, qdsr(i,1), nqd, dc, 1)) * qdwt(i)
wtsum = dasum (nqd, wt, 1)
lkhd = 0.d0
for (i=1;i<=nobs;i=i+1) {
    fit(i) = dexp (ddot (nsr, sr(1,i), 1, dc, 1)) / wtsum
    if (cntsum==0)  lkhd = lkhd + dlog (ff(i,1)+ff(i,2)*fit(i))
    else lkhd = lkhd + dlog (ff(i,1)+ff(i,2)*fit(i)) * dfloat (cnt(i))
}
if (cntsum==0)  lkhd = lkhd / dfloat (nobs)
else  lkhd = lkhd / dfloat (cntsum)
call  dsymv ('u', nr, 1.d0, q, nr, dc(2), 1, 0.d0, wk, 1)
lkhd = ddot (nr, dc(2), 1, wk, 1) / 2.d0 - lkhd
iter = 0
flag = 0
#   Iteration
repeat {
    iter = iter + 1
    # Calculate hessian and gradient
    for (i=1;i<=nsr;i=i+1)
        muwk(i) = ddot (nqd, wt, 1, qdsr(1,i), 1) / wtsum
    for (i=1;i<=nsr;i=i+1) {
        for (j=i;j<=nsr;j=j+1) {
            vwk(i,j) = 0.d0
            for (k=1;k<=nqd;k=k+1)
                vwk(i,j) = vwk(i,j) + wt(k) * qdsr(k,i) * qdsr(k,j)
            vwk(i,j) = vwk(i,j) / wtsum - muwk(i) * muwk(j)
        }
    }
    call  dset (nsr, 0.d0, mu, 1)
    call  dset (nsr*nsr, 0.d0, v, 1)
    if (cntsum==0) {
        for (k=1;k<=nobs;k=k+1) {
            tmp = ff(k,2)*fit(k)
            tmp = tmp / (ff(k,1)+tmp)
            call  daxpy (nsr, tmp, muwk, 1, mu, 1)
            call  daxpy (nsr, -tmp, sr(1,k), 1, mu, 1)
            for (i=1;i<=nsr;i=i+1) {
                for (j=i;j<=nsr;j=j+1) {
                    v(i,j) = v(i,j) + tmp*vwk(i,j)
                    v(i,j) = v(i,j) - tmp*(1.d0-tmp)*(muwk(i)-sr(i,k))*(muwk(j)-sr(j,k))
                }
            }
        }
        call  dscal (nsr, -1.d0/dfloat(nobs), mu, 1)
        call  dscal (nsr*nsr, 1.d0/dfloat(nobs), v, 1)
    }
    else {
        for (k=1;k<=nobs;k=k+1) {
            tmp = ff(k,2)*fit(k)
            tmp = tmp / (ff(k,1)+tmp)
            call  daxpy (nsr, tmp*dfloat(cnt(i)), muwk, 1, mu, 1)
            call  daxpy (nsr, -tmp*dfloat(cnt(i)), sr(1,k), 1, mu, 1)
            for (i=1;i<=nsr;i=i+1) {
                for (j=i;j<=nsr;j=j+1) {
                    v(i,j) = v(i,j) + tmp*vwk(i,j)*dfloat(cnt(i))
                    v(i,j) = v(i,j) - tmp*(1.d0-tmp)*dfloat(cnt(i))*
                      (muwk(i)-sr(i,k))*(muwk(j)-sr(j,k))
                }
            }
        }
        call  dscal (nsr, -1.d0/dfloat(cntsum), mu, 1)
        call  dscal (nsr*nsr, 1.d0/dfloat(cntsum), v, 1)            
    }
    call  dsymv ('u', nr, -1.d0, q, nr, dc(2), 1, 1.d0, mu(2), 1)
    for (i=2;i<=nsr;i=i+1) {
        for (j=i;j<=nsr;j=j+1)  v(i,j) = v(i,j) + q(i-1,j-1)
    }
    mumax = dabs(mu(idamax(nsr, mu, 1)))
    #   Cholesky factorization
    call  dmcdc (v, nsr, nsr, wk, jpvt, infowk)
    #   Update coefficients
    repeat {
        call  dcopy (nsr, mu, 1, dcnew, 1)
        call  dprmut (dcnew, nsr, jpvt, 0)
        call  dtrsl (v, nsr, nsr, dcnew, 11, infowk)
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
            lkhdnew = 0.d0
            for (i=1;i<=nobs;i=i+1) {
                tmp = ddot (nsr, sr(1,i), 1, dcnew, 1)
                if (tmp>3.d2) {
                    flag = flag + 1
                    break
                }
                fitnew(i) = dexp (tmp) / wtsumnew
                if (cntsum==0)  lkhdnew = lkhdnew + dlog (ff(i,1)+ff(i,2)*fitnew(i))
                else lkhdnew = lkhdnew + dlog (ff(i,1)+ff(i,2)*fitnew(i)) * dfloat (cnt(i))
            }
            if (cntsum==0)  lkhdnew = lkhdnew / dfloat (nobs)
            else  lkhdnew = lkhdnew / dfloat (cntsum)
            call  dsymv ('u', nr, 1.d0, q, nr, dcnew(2), 1, 0.d0, wk, 1)
            lkhdnew = ddot (nr, dcnew(2), 1, wk, 1) / 2.d0 - lkhdnew
        }
        #   Reset iteration with uniform starting value
        if (flag==1) {
            call  dset (nsr, 0.d0, dc, 1)
            call  dcopy (nqd, qdwt, 1, wt, 1)
            wtsum = dasum (nqd, wt, 1)
            call  dset (nobs, 1.d0/wtsum, fit, 1)
            lkhd = 0.d0
            for (i=1;i<=nobs;i=i+1) {
                if (cntsum==0)  lkhd = lkhd + dlog (ff(i,1)+ff(i,2)*fit(i))
                else  lkhd = lkhd + dlog (ff(i,1)+ff(i,2)*fit(i)) * dfloat (cnt(i))
            }
            if (cntsum==0)  lkhd = lkhd / dfloat (nobs)
            else  lkhd = lkhd / dfloat (cntsum)
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
        call  dset (nobs, 1.d0/wtsum, fit, 1)
        lkhd = 0.d0
        for (i=1;i<=nobs;i=i+1) {
            if (cntsum==0)  lkhd = lkhd + dlog (ff(i,1)+ff(i,2)*fit(i))
            else  lkhd = lkhd + dlog (ff(i,1)+ff(i,2)*fit(i)) * dfloat (cnt(i))
        }
        if (cntsum==0)  lkhd = lkhd / dfloat (nobs)
        else  lkhd = lkhd / dfloat (cntsum)
        iter = 0
        flag = 2
    }
    else {
        info = 2
        break
    }
}

return
end
