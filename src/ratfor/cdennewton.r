
#::::::::::::::::
#   cdennewton
#::::::::::::::::

subroutine  cdennewton (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qdrs, nqd, nx,
                        xxwt, qdwt, prec, maxiter, mchpr, jpvt, wk, info)

integer  nxis, nxi, nobs, nqd, nx, maxiter, jpvt(*), info
double precision  cd(*), q(nxi,*), rs(nxis,*), cntsum, cnt(*), qdrs(nqd,nxis,*),
                  xxwt(*), qdwt(*), prec, mchpr, wk(*)

integer  iwt, iwtsum, imrs, ifit, imu, imuwk, iv, ivwk, icdnew, iwtnew, iwtnewsum,
         ifitnew, iwk

iwt = 1
iwtsum = iwt + nqd*nx
imrs = iwtsum + nx
ifit = imrs + nxis
imu = ifit + nobs
imuwk = imu + nxis
iv = imuwk + nxis
ivwk = iv + nxis*nxis
icdnew = ivwk + nxis*nxis
iwtnew = icdnew + nxis
iwtnewsum = iwtnew + nqd*nx
ifitnew = iwtnewsum + nx
iwk = ifitnew + nobs

call  cdennewton1 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qdrs, nqd, nx, xxwt, qdwt,
                   prec, maxiter, mchpr, wk(iwt), wk(iwtsum), wk(imrs), wk(ifit),
                   wk(imu), wk(imuwk), wk(iv), wk(ivwk), jpvt, wk(icdnew),
                   wk(iwtnew), wk(iwtnewsum), wk(ifitnew), wk(iwk), info)

return
end



#:::::::::::::::::
#   cdennewton1
#:::::::::::::::::

subroutine  cdennewton1 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, qdrs, nqd, nx,
                         xxwt, qdwt, prec, maxiter, mchpr, wt, wtsum, mrs, fit,
                         mu, muwk, v, vwk, jpvt, cdnew, wtnew, wtnewsum, fitnew,
                         wk, info)

integer  nxis, nxi, nobs, nqd, nx, maxiter, jpvt(*), info
double precision  cd(*), q(nxi,*), rs(nxis,*), cntsum, cnt(*), qdrs(nqd,nxis,*),
                  xxwt(*), qdwt(*), prec, mchpr, wt(nqd,*), wtsum(*), mrs(*),
                  fit(*), mu(*), muwk(*), v(nxis,*), vwk(nxis,*), cdnew(*),
                  wtnew(nqd,*), wtnewsum(*), fitnew(*), wk(*)

integer  i, j, k, kk, iter, flag, rkv, idamax, infowk
double precision  norm, tmp, ddot, fitmean, lkhd, mumax, lkhdnew, disc, disc0, trc

#   Calculate constants
info = 0
for (i=1;i<=nxis;i=i+1) {
    mrs(i) = 0.d0
    if (!(cntsum>0.d0)) {
        for (j=1;j<=nobs;j=j+1)  mrs(i) = mrs(i) + rs(i,j)
        mrs(i) = mrs(i) / dble (nobs)
    }
    else {
        for (j=1;j<=nobs;j=j+1)  mrs(i) = mrs(i) + rs(i,j) * cnt(j)
        mrs(i) = mrs(i) / cntsum
    }
}
if (!(cntsum>0.d0))  trc = 1.d0 / dble (nobs)
else  trc = 1.d0 / cntsum
#   Initialization
norm = 0.d0
for (kk=1;kk<=nx;kk=kk+1) {
    wtsum(kk) = 0.d0
    for (i=1;i<=nqd;i=i+1) {
        wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1)) * qdwt(i)
        wtsum(kk) = wtsum(kk) + wt(i,kk)
    }
    norm = norm + xxwt(kk) * dlog (wtsum(kk))
}
fitmean = 0.d0
for (i=1;i<=nobs;i=i+1) {
    tmp = ddot (nxis, rs(1,i), 1, cd, 1)
    fit(i) = dexp (tmp)
    if (cntsum>0.d0)  tmp = tmp * cnt(i)
    fitmean = fitmean + tmp
}
call  dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, wk, 1)
lkhd = ddot (nxi, cd, 1, wk, 1) / 2.d0 - fitmean * trc + norm
iter = 0
flag = 0
#   Iteration
repeat {
    iter = iter + 1
    # Calculate hessian and gradient
    call  dset (nxis, 0.d0, mu, 1)
    call  dset (nxis*nxis, 0.d0, v, 1)
    for (kk=1;kk<=nx;kk=kk+1) {
        for (i=1;i<=nxis;i=i+1)
            muwk(i) = - ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1) / wtsum(kk)
        for (i=1;i<=nxis;i=i+1) {
            for (j=i;j<=nxis;j=j+1) {
                vwk(i,j) = 0.d0
                for (k=1;k<=nqd;k=k+1)
                    vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
                vwk(i,j) = vwk(i,j) / wtsum(kk) - muwk(i) * muwk(j)
            }
        }
        call  daxpy (nxis, xxwt(kk), muwk, 1, mu, 1)
        call  daxpy (nxis*nxis, xxwt(kk), vwk, 1, v, 1)
    }
    for (i=1;i<=nxi;i=i+1) {
        for (j=i;j<=nxi;j=j+1)  v(i,j) = v(i,j) + q(i,j)
    }
    call  daxpy (nxis, 1.d0, mrs, 1, mu, 1)
    call  dsymv ('u', nxi, -1.d0, q, nxi, cd, 1, 1.d0, mu, 1)
    mumax = dabs(mu(idamax(nxis, mu, 1)))
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
        norm = 0.d0
        for (kk=1;kk<=nx;kk=kk+1) {
            wtnewsum(kk) = 0.d0
            for (i=1;i<=nqd;i=i+1) {
                wtnew(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cdnew, 1)) * qdwt(i)
                wtnewsum(kk) = wtnewsum(kk) + wtnew(i,kk)
            }
            norm = norm + xxwt(kk) * dlog (wtnewsum(kk))
        }
        if ((flag==0)|(flag==2)) {
            fitmean = 0.d0
            for (i=1;i<=nobs;i=i+1) {
                tmp = ddot (nxis, rs(1,i), 1, cdnew, 1)
                if (tmp>3.d2) {
                    flag = flag + 1
                    break
                }
                fitnew(i) = dexp (tmp)
                if (cntsum>0.d0)  tmp = tmp * cnt(i)
                fitmean = fitmean + tmp
            }
            call  dsymv ('u', nxi, 1.d0, q, nxi, cdnew, 1, 0.d0, wk, 1)
            lkhdnew = ddot (nxi, cdnew, 1, wk, 1) / 2.d0 - fitmean * trc + norm
        }
        #   Reset iteration with uniform starting value
        if (flag==1) {
            call  dset (nxis, 0.d0, cd, 1)
            for (kk=1;kk<=nx;kk=kk+1)
                call  dcopy (nqd, qdwt, 1, wt(1,kk), 1)
            call  dset (nx, 1.d0, wtsum, 1)
            call  dset (nobs, 1.d0, fit, 1)
            lkhd = 0.d0
            iter = 0
            break
        }
        if (flag==3)  break
        if (lkhdnew-lkhd<1.d1*(1.d0+dabs(lkhd))*mchpr)  break
        call  dscal (nxis, .5d0, mu, 1)
        if (dabs(mu(idamax(nxis, mu, 1))/mumax)<1.d1*mchpr)  break
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
    for (kk=1;kk<=nx;kk=kk+1) {
        for (i=1;i<=nqd;i=i+1)
            disc = dmax1 (disc, dabs(wt(i,kk)-wtnew(i,kk))/(1.d0+dabs(wt(i,kk))))
    }
    for (i=1;i<=nobs;i=i+1)
        disc = dmax1 (disc, dabs(fit(i)-fitnew(i))/(1.d0+dabs(fit(i))))
    disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
    disc0 = dmax1 ((mumax/(1.d0+lkhd))**2, dabs(lkhd-lkhdnew)/(1+dabs(lkhd)))
    #   Set to new values
    call  dcopy (nxis, cdnew, 1, cd, 1)
    call  dcopy (nqd*nx, wtnew, 1, wt, 1)
    call  dcopy (nx, wtnewsum, 1, wtsum, 1)
    call  dcopy (nobs, fitnew, 1, fit, 1)
    lkhd = lkhdnew
    #   Check convergence
    if (disc0<prec)  break
    if (disc<prec)  break
    if (iter<maxiter)  next
    if (flag==0) {
        #   Reset iteration with uniform starting value
        call  dset (nxis, 0.d0, cd, 1)
        for (kk=1;kk<=nx;kk=kk+1)
            call  dcopy (nqd, qdwt, 1, wt(1,kk), 1)
        call  dset (nx, 1.d0, wtsum, 1)
        call  dset (nobs, 1.d0, fit, 1)
        lkhd = 0.d0
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
    if (cntsum>0.d0)  call  dscal (nxis, dsqrt(cnt(i)), rs(1,i), 1)
    call  dtrsl (v, nxis, nxis, rs(1,i), 11, infowk)
}
trc = ddot (nobs*nxis, rs, 1, rs, 1)
if (!(cntsum>0.d0)) {
    trc = trc / dble(nobs) / (dble(nobs)-1.d0)
    lkhd = 0.d0
    for (i=1;i<=nobs;i=i+1)  lkhd = lkhd + dlog (fit(i))
    lkhd = lkhd / dble (nobs)
}
else {
    trc = trc / cntsum / (cntsum-1.d0)
    lkhd = 0.d0
    for (i=1;i<=nobs;i=i+1)  lkhd = lkhd + cnt(i) * dlog (fit(i))
    lkhd = lkhd / cntsum
}
for (kk=1;kk<=nx;kk=kk+1)  lkhd = lkhd - xxwt(kk) * dlog (wtsum(kk))
wtsum(1) = lkhd
wtsum(2) = trc

return
end


#::::::::::::
#   cdenrkl
#::::::::::::

subroutine  cdenrkl (cd, nxis, qdrs, nqd, nx, xxwt, qdwt, wt0, mchpr,
                     wt, wtnew, mu, muwk, v, vwk, jpvt, cdnew,
                     prec, maxiter, info)

integer  nxis, nqd, nx, jpvt(*), maxiter, info
double precision  cd(*), qdrs(nqd,nxis,*), xxwt(*), qdwt(*), wt0(nqd,*),
                  mchpr, wt(nqd,*), wtnew(nqd,*), mu(*), muwk(*), v(nxis,*),
                  vwk(nxis,*), cdnew(*), prec

integer  i, j, k, kk, iter, flag, idamax, infowk
double precision  ddot, dasum, rkl, tmp, mumax, rklnew, disc, disc0


#   Initialization
for (kk=1;kk<=nx;kk=kk+1) {
    for (i=1;i<=nqd;i=i+1) {
        wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1)) * qdwt(i)
    }
    call  dscal (nqd, 1.d0/dasum(nqd,wt(1,kk),1), wt(1,kk), 1)
}
rkl = 0.d0
for (kk=1;kk<=nx;kk=kk+1) {
    tmp = 0.d0
    for (i=1;i<=nqd;i=i+1)  tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
    rkl = rkl + xxwt(kk) * tmp
}
iter = 0
flag = 0
#   Iteration
repeat {
    iter = iter + 1
    # Calculate hessian and gradient
    call  dset (nxis, 0.d0, mu, 1)
    call  dset (nxis*nxis, 0.d0, v, 1)
    for (kk=1;kk<=nx;kk=kk+1) {
        for (i=1;i<=nxis;i=i+1)
            muwk(i) = ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1)
        for (i=1;i<=nxis;i=i+1) {
            for (j=i;j<=nxis;j=j+1) {
                vwk(i,j) = 0.d0
                for (k=1;k<=nqd;k=k+1)
                    vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
                vwk(i,j) = vwk(i,j) - muwk(i) * muwk(j)
            }
            muwk(i) = ddot (nqd, wt0(1,kk), 1, qdrs(1,i,kk), 1) - muwk(i)
        }
        call  daxpy (nxis, xxwt(kk), muwk, 1, mu, 1)
        call  daxpy (nxis*nxis, xxwt(kk), vwk, 1, v, 1)
    }
    mumax = dabs(mu(idamax(nxis, mu, 1)))
    #   Cholesky factorization
    for (i=1;i<=nxis;i=i+1)  jpvt(i) = 0
    call  dmcdc (v, nxis, nxis, cdnew, jpvt, infowk)
    #   Update coefficients
    repeat {
        call  dcopy (nxis, mu, 1, cdnew, 1)
        call  dprmut (cdnew, nxis, jpvt, 0)
        call  dtrsl (v, nxis, nxis, cdnew, 11, infowk)
        call  dtrsl (v, nxis, nxis, cdnew, 01, infowk)
        call  dprmut (cdnew, nxis, jpvt, 1)
        call  daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
        for (kk=1;kk<=nx;kk=kk+1) {
            for (i=1;i<=nqd;i=i+1) {
                wtnew(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cdnew, 1)) * qdwt(i)
            }
            call  dscal (nqd, 1.d0/dasum(nqd,wtnew(1,kk),1), wtnew(1,kk), 1)
        }
        if ((flag==0)|(flag==2)) {
            rklnew = 0.d0
            for (kk=1;kk<=nx;kk=kk+1) {
                tmp = 0.d0
                for (i=1;i<=nqd;i=i+1)  tmp = tmp + dlog(wt0(i,kk)/wtnew(i,kk)) * wt0(i,kk)
                rklnew = rklnew + xxwt(kk) * tmp
            }
        }
        #   Reset iteration with uniform starting value
        if (flag==1) {
            call  dset (nxis, 0.d0, cd, 1)
            for (kk=1;kk<=nx;kk=kk+1) {
                call  dcopy (nqd, qdwt, 1, wt(1,kk), 1)
                call  dscal (nqd, 1.d0/dasum(nqd,wt(1,kk),1), wt(1,kk), 1)
            }
            rkl = 0.d0
            for (kk=1;kk<=nx;kk=kk+1) {
                tmp = 0.d0
                for (i=1;i<=nqd;i=i+1)  tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
                rkl = rkl + xxwt(kk) * tmp
            }
            iter = 0
            break
        }
        if (flag==3)  break
        if (rklnew-rkl<1.d1*(1.d0+dabs(rkl))*mchpr)  break
        call  dscal (nxis, .5d0, mu, 1)
        if (dabs(mu(idamax(nxis, mu, 1))/mumax)<1.d1*mchpr)  break
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
    for (kk=1;kk<=nx;kk=kk+1) {
        for (i=1;i<=nqd;i=i+1)
            disc = dmax1 (disc, dabs(wt(i,kk)-wtnew(i,kk))/(1.d0+dabs(wt(i,kk))))
    }
    disc = dmax1 (disc, (mumax/(1.d0+dabs(rkl)))**2)
    disc0 = dmax1 ((mumax/(1.d0+rkl))**2, dabs(rkl-rklnew)/(1+dabs(rkl)))
    #   Set to new values
    call  dcopy (nxis, cdnew, 1, cd, 1)
    call  dcopy (nqd*nx, wtnew, 1, wt, 1)
    rkl = rklnew
    #   Check convergence
    if (disc0<prec)  break
    if (disc<prec)  break
    if (iter<maxiter)  next
    if (flag==0) {
        #   Reset iteration with uniform starting value
        call  dset (nxis, 0.d0, cd, 1)
        for (kk=1;kk<=nx;kk=kk+1) {
            call  dcopy (nqd, qdwt, 1, wt(1,kk), 1)
            call  dscal (nqd, 1.d0/dasum(nqd,wt(1,kk),1), wt(1,kk), 1)
        }
        rkl = 0.d0
        for (kk=1;kk<=nx;kk=kk+1) {
            tmp = 0.d0
            for (i=1;i<=nqd;i=i+1)  tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
            rkl = rkl + xxwt(kk) * tmp
        }
        iter = 0
        flag = 2
    }
    else {
        info = 2
        break
    }
}
#   calculate rkl
rkl = 0.d0
for (kk=1;kk<=nx;kk=kk+1) {
    tmp = 0.d0
    for (i=1;i<=nqd;i=i+1)  tmp = tmp + dlog(wt0(i,kk)/wt(i,kk)) * wt0(i,kk)
    rkl = rkl + xxwt(kk) * tmp
}
wt(1,1) = rkl

return
end
