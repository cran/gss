subroutine  dmudr0 (vmu,_
                    s, lds, nobs, nnull, q, ldqr, ldqc, nq, y,_     # inputs
                    tol, init, prec, maxite,_                       # tune para
                    theta, nlaht, score, varht, c, d,_              # outputs
                    iwk, wk, info)

integer           vmu
integer           lds, nobs, nnull, ldqr, ldqc, nq, init, maxite,_
                  info, iwk(*)

double precision  s(lds,*), q(ldqr,ldqc,*), y(*), tol, prec,_
                  theta(*), nlaht, score, varht, c(*), d(*),_
                  wk(*)

character         vmu1
integer  n, n0
integer  iqraux, itraux, itwk, iqwk, iywk, ithewk, ihes, igra, ihwk1, ihwk2,_
         igwk1, igwk2, ikwk, iwork1, iwork2, ijpvt, ipvtwk

if ( vmu == 1 )  vmu1 = 'v'
if ( vmu == 2 )  vmu1 = 'm'
if ( vmu == 3 )  vmu1 = 'u'

n = nobs
n0 = nnull

iqraux = 1
itraux = iqraux + n0
itwk = itraux + (n-n0-2)
iqwk = itwk + 2 * (n-n0)
iywk = iqwk + n * n
ithewk = iywk + n
ihes = ithewk + nq
igra = ihes + nq * nq
ihwk1 = igra + nq
ihwk2 = ihwk1 + nq * nq
igwk1 = ihwk2 + nq * nq
igwk2 = igwk1 + nq
ikwk = igwk2 + nq
iwork1 = ikwk + (n-n0) * (n-n0) * nq
iwork2 = iwork1 + n
ijpvt = 1
ipvtwk = ijpvt + n0

call  dmudr1 (vmu1,_
             s, lds, nobs, nnull, q, ldqr, ldqc, nq, y,_     # inputs
             tol, init, prec, maxite,_                       # tune para
             theta, nlaht, score, varht, c, d,_              # outputs
             wk(iqraux), iwk(ijpvt), wk(itwk), wk(itraux), wk(iqwk),_
             wk(iywk), wk(ithewk), wk(ihes), wk(igra), wk(ihwk1),_
             wk(ihwk2), wk(igwk1), wk(igwk2), iwk(ipvtwk), wk(ikwk),_
             wk(iwork1), wk(iwork2),_
             info)

return
end
