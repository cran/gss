
#:::::::::::::
#   hzdaux1
#:::::::::::::

subroutine  hzdaux1 (cd, nxis, q, nxi, qdrs, nqd, qdwt, nx, mchpr, wt, v, vwk, jpvt)

integer  nxis, nxi, nqd, nx, jpvt(*)
double precision  cd(*), q(nxi,*), qdrs(nqd,nxis,*), qdwt(nqd,*), mchpr, wt(nqd,*),
                  v(nxis,*), vwk(nxis,*)

integer  i, j, k, kk, rkv
double precision  ddot

#   Initialization
for (kk=1;kk<=nx;kk=kk+1) {
    for (i=1;i<=nqd;i=i+1)
        wt(i,kk) = qdwt(i,kk) * dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1))
}
#   H matrix
call  dset (nxis*nxis, 0.d0, v, 1)
for (kk=1;kk<=nx;kk=kk+1) {
    for (i=1;i<=nxis;i=i+1) {
        for (j=i;j<=nxis;j=j+1) {
            vwk(i,j) = 0.d0
            for (k=1;k<=nqd;k=k+1)
                vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
        }
    }
    call  daxpy (nxis*nxis, 1.d0, vwk, 1, v, 1)
}
for (i=1;i<=nxi;i=i+1) {
    for (j=i;j<=nxi;j=j+1)  v(i,j) = v(i,j) + q(i,j)
}
#   Cholesky factorization
for (i=1;i<=nxis;i=i+1)  jpvt(i) = 0
call  dchdc (v, nxis, nxis, vwk, jpvt, 1, rkv)
while (v(rkv,rkv)<v(1,1)*dsqrt(mchpr))  rkv = rkv - 1
for (i=rkv+1;i<=nxis;i=i+1) {
    v(i,i) = v(1,1)
    call  dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
}

return
end


#:::::::::::::
#   hzdaux2
#:::::::::::::

subroutine  hzdaux2 (v, nxis, jpvt, r, nr, se)

double precision  v(nxis,*), r(nxis,*), se(*)
integer  nxis, jpvt(*), nr

double precision  ddot
integer  i, infowk

for (i=1;i<=nr;i=i+1) {
    call  dprmut (r(1,i), nxis, jpvt, 0)
    call  dtrsl (v, nxis, nxis, r(1,i), 11, infowk)
    se(i) = dsqrt (ddot (nxis, r(1,i), 1, r(1,i), 1))
}

return
end
