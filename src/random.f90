function random(idum)

    use get_kinds

    implicit none

    integer (kind=int) :: idum
    real (kind=double), parameter :: mbig=4.d+06, mseed=1618033.d00
    real (kind=double), parameter :: mz=0.d00, fac=1.d00/mbig
    integer (kind=int) :: i, iff, ii, inext, inextp, k
    real (kind=double) :: random, mj, mk, ma(55)
    save iff, inext, inextp, ma
    data iff/0/

    if (idum <= 0 .or. iff == 0) then
        iff = 1
        mj = abs(mseed - abs(idum))
        mj = mod(mj, mbig)
        ma(55) = mj
        mk = 1
        do i = 1, 54
            ii = mod(21 * i, 55)
            ma(ii) = mk
            mk = mj - mk
            if (mk < mz) mk = mk + mbig
            mj = ma(ii)
        end do
        do k = 1, 4
            do i = 1, 55
                ma(i) = ma(i) - ma(1 + mod(i + 30, 55))
                if (ma(i) < mz) ma(i) = ma(i) + mbig
            end do
        end do
        inext = 0
        inextp = 31
        idum = 1
    end if

    inext = inext + 1
    if (inext == 56) inext = 1
    inextp = inextp + 1
    if (inextp == 56) inextp = 1
    mj = ma(inext) - ma(inextp)
    if (mj < mz) mj = mj + mbig
    ma(inext) = mj
    random = mj * fac
    return
    
end function random