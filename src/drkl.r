
#::::::::::
#   drkl
#::::::::::

subroutine  drkl (cd, nxis, qdrs, nqd, qdwt, wt0, mchpr, wt, eta, mu, v,
                  jpvt, wk, cdnew, wtnew, prec, maxiter, info)

integer  nxis, nqd, jpvt(*), maxiter, info
double precision  cd(*), qdrs(nqd,*), qdwt(*), wt0(*), mchpr, wt(*), eta(*),
                  mu(*), v(nxis,*), wk(*), cdnew(*), wtnew(*), prec

integer  i, j, k, iter, flag, idamax, infowk
double precision  wtsum, ddot, rkl, mumax, wtsumnew, rklnew, disc, disc0


#   Initialization
info = 0
wtsum = 0.d0
rkl = 0.d0
for (i=1;i<=nqd;i=i+1) {
    wt0(i) = qdwt(i) * wt0(i)
    eta(i) = ddot (nxis, qdrs(i,1), nqd, cd, 1)
    wt(i) = qdwt(i) * dexp (eta(i))
    wtsum = wtsum + wt(i)
    rkl = rkl + dlog(wt0(i)/wt(i)) * wt0(i)
}
rkl = rkl + dlog(wtsum)
iter = 0
flag = 0
#   Iteration
repeat {
    iter = iter + 1
    # Calculate hessian and gradient
    for (i=1;i<=nxis;i=i+1)
        mu(i) = ddot (nqd, wt, 1, qdrs(1,i), 1) / wtsum
    for (i=1;i<=nxis;i=i+1) {
        for (j=i;j<=nxis;j=j+1) {
            v(i,j) = 0.d0
            for (k=1;k<=nqd;k=k+1)
                v(i,j) = v(i,j) + wt(k) * qdrs(k,i) * qdrs(k,j)
            v(i,j) = v(i,j) / wtsum - mu(i) * mu(j)
        }
    }
    for (i=1;i<=nxis;i=i+1) {
        mu(i) = ddot (nqd, wt0, 1, qdrs(1,i), 1) - mu(i)
    }
    mumax = dabs(mu(idamax(nxis, mu, 1)))
    #   Cholesky factorization
    for (i=1;i<=nxis;i=i+1)  jpvt(i) = 0
    call  dmcdc (v, nxis, nxis, wk, jpvt, infowk)
    #   Update coefficients
    repeat {
        call  dcopy (nxis, mu, 1, cdnew, 1)
        call  dprmut (cdnew, nxis, jpvt, 0)
        call  dtrsl (v, nxis, nxis, cdnew, 11, infowk)
        call  dtrsl (v, nxis, nxis, cdnew, 01, infowk)
        call  dprmut (cdnew, nxis, jpvt, 1)
        call  daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
        wtsumnew = 0.d0
        rklnew = 0.d0
        for (i=1;i<=nqd;i=i+1) {
            eta(i) = ddot (nxis, qdrs(i,1), nqd, cdnew, 1)
            if (eta(i)>3.d2) {
                flag = flag + 1
                break
            }
            wtnew(i) = qdwt(i) * dexp (eta(i))
            wtsumnew = wtsumnew + wtnew(i)
            rklnew = rklnew + dlog(wt0(i)/wtnew(i)) * wt0(i)
        }
        rklnew = rklnew + dlog (wtsumnew)
        if (flag==1) {
            #   Reset iteration with uniform starting value
            call  dset (nxis, 0.d0, cd, 1)
            wtsum = 0.d0
            rkl = 0.d0
            for (i=1;i<=nqd;i=i+1) {
                wt(i) = qdwt(i)
                wtsum = wtsum + wt(i)
                rkl = rkl + dlog(wt0(i)/wt(i)) * wt0(i)
            }
            rkl = rkl + dlog (wtsum)
            iter = 0
            break
        }
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
    for (i=1;i<=nqd;i=i+1)
        disc = dmax1 (disc, dabs(wt(i)-wtnew(i))/(1.d0+dabs(wt(i))))
    disc = dmax1 (disc, (mumax/(1.d0+rkl))**2)
    disc0 = dmax1 ((mumax/(1.d0+rkl))**2, dabs(rkl-rklnew)/(1+dabs(rkl)))
    #   Set to new values
    call  dcopy (nxis, cdnew, 1, cd, 1)
    call  dcopy (nqd, wtnew, 1, wt, 1)
    wtsum = wtsumnew
    rkl = rklnew
    #   Check convergence
    if (disc0<prec)  break
    if (disc<prec)  break
    if (iter<maxiter)  next
    if (flag==0) {
        #   Reset iteration with uniform starting value
        call  dset (nxis, 0.d0, cd, 1)
        wtsum = 0.d0
        rkl = 0.d0
        for (i=1;i<=nqd;i=i+1) {
            wt(i) = qdwt(i)
            wtsum = wtsum + wt(i)
            rkl = rkl + dlog(wt0(i)/wt(i)) * wt0(i)
        }
        rkl = rkl + dlog (wtsum)
        iter = 0
        flag = 2
    }
    else {
        info = 2
        break
    }
}
#   return projection density on the mesh
wtsum = 0.d0
for (i=1;i<=nqd;i=i+1) {
    wt0(i) = dexp (ddot (nxis, qdrs(i,1), nqd, cd, 1))
    wtsum = wtsum + qdwt(i) * wt0(i)
}
call  dscal (nqd, 1.d0/wtsum, wt0, 1)

return
end
