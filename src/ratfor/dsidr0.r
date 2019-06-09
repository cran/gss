subroutine  dsidr0 (vmu,_
                    s, lds, nobs, nnull, y, q, ldq,_       # data
                    tol, job, limnla,_                     # job requests
                    nlaht, score, varht, c, d,_            # output
                    qraux, jpvt, wk,_                      # work arrays
                    info)                                  # error message

integer           vmu
integer           lds, nobs, nnull, ldq, job, jpvt(*), info
double precision  s(lds,*), y(*), q(ldq,*), tol, limnla(2), nlaht, score(*),_
                  varht, c(*), d(*), qraux(*), wk(*)

character*1       vmu1

if ( vmu == 1 )  vmu1 = 'v'
if ( vmu == 2 )  vmu1 = 'm'
if ( vmu == 3 )  vmu1 = 'u'

info = 0

#   check dimension
if ( nnull < 1 | nnull >= nobs | nobs > lds | nobs > ldq ) {
    info = -1
    return
}

#   main process

call  dstup (s, lds, nobs, nnull, qraux, jpvt, y, q, ldq, nobs, 1, info,_
             wk)
if ( info != 0 )  return

call  dcore (vmu1, q, ldq, nobs, nnull, tol, y, job, limnla, nlaht, score,_
             varht, info, wk, wk(2*nobs+1))
if ( info != 0 )  return

call  dcoef (s, lds, nobs, nnull, qraux, jpvt, y, q, ldq, nlaht, c, d,_
             info, wk)

return
end
