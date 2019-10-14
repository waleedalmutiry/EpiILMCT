!######################################################################
!# MODULE: logliksinr
!# AUTHORS: 
!#     Waleed Almutiry <walmutir@uoguelph.ca>,  
!#     Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca> and
!#     Rob Deardon <robert.deardon@ucalgary.ca> 
!# 
!# DESCRIPTION:
!#
!#     To calculate the log-likelihood for the SINR continuous-time ILMs:
!#
!#     This program is free software; you can redistribute it and/or
!#     modify it under the terms of the GNU General Public License, 
!#     version 3,  as published by the Free Software Foundation.
!# 
!#     This program is distributed in the hope that it will be useful, 
!#     but without any warranty; without even the implied warranty of
!#     merchantability or fitness for a particular purpose.  See the GNU
!#     General Public License,  version 3,  for more details.
!# 
!#     A copy of the GNU General Public License,  version 3,  is available
!#     at http://www.r-project.org/Licenses/GPL-3
!# 
!# Part of the R/EpiILMCT package
!# Contains: 
!#           loglikcontilmsinr ............ subroutine
!#           gammapdf ..................... function
!#           gammalog ..................... function
!#
!######################################################################

module logliksinr
    use ISO_C_BINDING
    implicit none
    public :: loglikcontilmsinr_f

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 				Log LIKELIHOOD subroutine				 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine loglikcontilmsinr_f(n, ninfected, num, nsuspar, ntranspar, cc, d3, epidat,  &
    & suscov, transcov, suspar, powersus, transpar, powertrans, kernelpar, spark, gamma, deltain1, deltain2, deltanr1,  &
    & deltanr2, likk) bind(C,  name="loglikcontilmsinr_f_")

    external infinity_value

    integer (C_INT) :: j, i, r, m
    integer (C_INT), intent(in) :: n, ninfected, num, nsuspar, ntranspar    !integers
    real (C_DOUBLE), intent(in) :: deltain1, deltain2, deltanr1, deltanr2 ! parameters of the infectious period distribution
    real (C_DOUBLE), intent(in) :: spark, gamma                         ! spark & notification effect parameters
    real (C_DOUBLE), intent(in), dimension(2) :: kernelpar              ! parameter of the kernel function
    real (C_DOUBLE), intent(in), dimension(n, n) :: d3, cc                 ! network & distance matrices
    real (C_DOUBLE), intent(in), dimension(n, 6) :: epidat                ! epidemic data
    real (C_DOUBLE), intent(in), dimension(n, nsuspar) :: suscov          ! susceptibility covariates
    real (C_DOUBLE), intent(in), dimension(n, ntranspar) :: transcov      ! transmissibility covariates
    real (C_DOUBLE), intent(in), dimension(nsuspar) :: suspar, powersus   ! susceptibility parameters
    real (C_DOUBLE), intent(in), dimension(ntranspar) :: transpar, powertrans ! transmissibility parameters
    real (C_DOUBLE), dimension(nsuspar) :: suscov1                      ! suseptible covar. for one individual
    real (C_DOUBLE), dimension(ntranspar) :: transcov1                  ! transmissibilty covar. for one individual
    real (C_DOUBLE), intent(out) :: likk                                ! OUTPUT
    real (C_DOUBLE) :: rate, tt, ss, qa1, qa2, likk1, likk2, Inf
    real (C_DOUBLE), dimension(n) :: gh, sss
    real (C_DOUBLE), dimension(ninfected) :: rt, df1, df2, rrate, gammain, gammanr

    call infinity_value(Inf)

    SELECT CASE (num)
    CASE (1)
! network-based continuous ILMs

! getting t_obs (maximum removal time):
        tt = maxval(epidat(1:ninfected, 2))

!calculate the exponent part of the spark term:
        do r = 1,  n
            sss(r) = spark * (min(tt, epidat(r, 6)) - epidat(1, 6))
        end do
        ss = sum(sss)

! calculate the terms of the exponent part of the double summation:
        df1 = 0.0_c_double
        do i =1, ninfected
            gh = 0.0_c_double
            do j = 1, n
                if (j .le. ninfected) then
                    if (j .ne. i) then
                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do
                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do
                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)
                        rate = qa1 * qa2 * (cc(int(epidat(i, 1)), int(epidat(j, 1))))
                        gh(j) = ( minval( (/ tt, epidat(i, 4), epidat(j, 6) /) ) - min(epidat(i, 6), epidat(j, 6)) ) * (rate)
                    else
                        rate  = 0.0_c_double
                        gh(j) = 0.0_c_double

                    end if
                else
                    do m = 1,  nsuspar
                        suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                    end do
                    do m = 1,  ntranspar
                        transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                    end do
                    qa1 = dot_product(suspar, suscov1)
                    qa2 = dot_product(transpar, transcov1)
                    rate = qa1 * qa2 * (cc(int(epidat(i, 1)), int(epidat(j, 1))))
                    gh(j)   = ( epidat(i, 4) - epidat(i, 6)) * (rate)
                end if
            end do
            df1(i) = sum(gh)
        end do

            df2 = 0.0_c_double
            do i =1, ninfected
                gh = 0.0_c_double
                do j = 1, n
                    if (j .le. ninfected) then
                        if (j .ne. i) then
                            do m = 1,  nsuspar
                                suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                            end do
                            do m = 1,  ntranspar
                                transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                            end do
                            qa1 = dot_product(suspar, suscov1)
                            qa2 = dot_product(transpar, transcov1)
                            rate = gamma * (qa1 * qa2 * (cc(int(epidat(i, 1)), int(epidat(j, 1)))))
                            gh(j) = ( ( minval( (/ tt, epidat(i, 2), epidat(j, 6) /) ) - min(epidat(i, 6), epidat(j, 6)) ) - &
                                & ( minval( (/ tt, epidat(i, 4), epidat(j, 6) /) ) - min(epidat(i, 6), epidat(j, 6)) ) ) * (rate)
                        else
                            rate    = 0.0_c_double
                            gh(j)   = 0.0_c_double
                        end if
                    else
                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do
                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do
                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)
                        rate = gamma * (qa1 * qa2 * (cc(int(epidat(i, 1)), int(epidat(j, 1)))))
                        gh(j) = ( (epidat(i, 2) - epidat(i, 6)) - (epidat(i, 4) - epidat(i, 6)) ) * (rate)
                    end if
                end do
                df2(i) = sum(gh)
            end do
        likk1 = (-(sum(df1)+sum(df2)+ss))

! calculate the first part of the likelihood:

        rt(1) = 1.0_c_double
        do j = 2,  ninfected
            rrate = 0.0_c_double
            do i = 1,  ninfected
                if (i .ne. j) then
                    if ((epidat(j, 6) .gt. epidat(i, 6)) .and.  (epidat(j, 6) .le. epidat(i, 4)) ) then
                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do
                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do
                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)
                        rrate(i) = (qa1 * qa2 * (cc(int(epidat(i, 1)), int(epidat(j, 1)))))
                    else if ((epidat(j, 6) .gt. epidat(i, 4)) .and.  (epidat(j, 6) .le. epidat(i, 2)) ) then
                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do
                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do
                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)
                        rrate(i) = gamma * (qa1 * qa2 * (cc(int(epidat(i, 1)), int(epidat(j, 1)))))
                    else
                        rrate(i)   = 0.0_c_double
                    end if
                else
                    rrate(i)   = 0.0_c_double
                end if
            end do
            rt(j) = sum(rrate) + spark
        end do

        likk2 = sum(log(rt))


! TOTAL LOG-LIKELIHOOD:

        do i = 1,  ninfected
            gammain(i) = gammapdf ( deltain1,  deltain2,  epidat(i, 5) )
        end do
        do i = 1,  ninfected
            gammanr(i) = gammapdf ( deltanr1,  deltanr2,  epidat(i, 3) )
        end do

        likk = likk1 + likk2 + sum(gammain) + sum(gammanr)


!##################################################################
!##################################################################

    CASE (2)
        ! distance-based continuous ILMs (powerlaw kernel)

        ! getting t_obs (maximum removal time):
        tt = maxval(epidat(1:ninfected, 2))

        !calculate the exponent part of the spark term:
        do r = 1,  n
            sss(r) = spark * (min(tt, epidat(r, 6)) - epidat(1, 6))
        end do
        ss = sum(sss)

        ! calculate the terms of the exponent part of the double summation:
        df1 = 0.0_c_double

        do i =1, ninfected
            gh = 0.0_c_double
            do j = 1, n
                if (j .le. ninfected) then
                    if (j .ne. i) then
                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do
                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do
                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)
                        rate = qa1 * qa2 * (d3(int(epidat(i, 1)), int(epidat(j, 1)))**(-kernelpar(1)))
                        gh(j)   = ( minval( (/ tt, epidat(i, 4), epidat(j, 6) /) ) - min(epidat(i, 6), epidat(j, 6)) ) * (rate)
                    else
                        rate    = 0.0_c_double
                        gh(j)   = 0.0_c_double
                    end if
                else
                    do m = 1,  nsuspar
                        suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                    end do
                    do m = 1,  ntranspar
                        transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                    end do
                    qa1 = dot_product(suspar, suscov1)
                    qa2 = dot_product(transpar, transcov1)
                    rate = qa1 * qa2 * (d3(int(epidat(i, 1)), int(epidat(j, 1)))**(-kernelpar(1)))
                    gh(j)   = ( epidat(i, 4) - epidat(i, 6)) * (rate)
                end if
            end do
            df1(i) = sum(gh)
        end do

        df2 = 0.0_c_double
        do i =1, ninfected
            gh = 0.0_c_double
            do j = 1, n
                if (j .le. ninfected) then
                    if (j .ne. i) then
                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do
                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do
                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)
                        rate = gamma * (qa1 * qa2 * (d3(int(epidat(i, 1)), int(epidat(j, 1)))**(-kernelpar(1))))
                        gh(j) = ( ( minval( (/ tt, epidat(i, 2), epidat(j, 6) /) ) - min(epidat(i, 6), epidat(j, 6)) ) - &
                        & ( minval( (/ tt, epidat(i, 4), epidat(j, 6) /) ) - min(epidat(i, 6), epidat(j, 6)) ) ) * (rate)
                    else
                        rate    = 0.0_c_double
                        gh(j)   = 0.0_c_double
                    end if
                else
                    do m = 1,  nsuspar
                        suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                    end do
                    do m = 1,  ntranspar
                        transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                    end do
                    qa1 = dot_product(suspar, suscov1)
                    qa2 = dot_product(transpar, transcov1)
                    rate = gamma * (qa1 * qa2 * (d3(int(epidat(i, 1)), int(epidat(j, 1)))**(-kernelpar(1))))
                    gh(j)   = ( (epidat(i, 2) - epidat(i, 6)) - (epidat(i, 4) - epidat(i, 6)) ) * (rate)
                end if
            end do
            df2(i) = sum(gh)
        end do

        likk1 = (-(sum(df1)+sum(df2)+ss))

! calculate the first part of the likelihood:

        rt(1) = 1.0_c_double
        do j = 2,  ninfected
            rrate = 0.0_c_double
            do i = 1,  ninfected
                if (i .ne. j) then
                    if ((epidat(j, 6) .gt. epidat(i, 6)) .and.  (epidat(j, 6) .le. epidat(i, 4)) ) then
                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do
                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do
                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)
                        rrate(i) = (qa1 * qa2 * (d3(int(epidat(i, 1)), int(epidat(j, 1)))**(-kernelpar(1))))
                    else if ((epidat(j, 6) .gt. epidat(i, 4)) .and.  (epidat(j, 6) .le. epidat(i, 2)) ) then
                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do
                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do
                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)
                        rrate(i) = gamma * (qa1 * qa2 * (d3(int(epidat(i, 1)), int(epidat(j, 1)))**(-kernelpar(1))))
                    else
                        rrate(i)   = 0.0_c_double
                    end if
                else
                    rrate(i)   = 0.0_c_double
                end if
            end do
            rt(j) = sum(rrate) + spark
        end do

        likk2 = sum(log(rt))

! TOTAL LOG-LIKELIHOOD:

        do i = 1,  ninfected
            gammain(i) = gammapdf ( deltain1,  deltain2,  epidat(i, 5) )
        end do
        do i = 1,  ninfected
            gammanr(i) = gammapdf ( deltanr1,  deltanr2,  epidat(i, 3) )
        end do

        likk = likk1 + likk2 + sum(gammain) + sum(gammanr)

    !##################################################################
    !##################################################################

    CASE (3)
! distance-based continuous ILMs (Cauchy kernel)

! getting t_obs (maximum removal time):
        tt = maxval(epidat(1:ninfected, 2))

!calculate the exponent part of the spark term:
        do r = 1,  n
            sss(r) = spark * (min(tt, epidat(r, 6)) - epidat(1, 6))
        end do

        ss = sum(sss)

! calculate the terms of the exponent part of the double summation:
        df1 = 0.0_c_double

        do i =1, ninfected

            gh = 0.0_c_double

            do j = 1, n

                if (j .le. ninfected) then

                    if (j .ne. i) then
                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do
                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do

                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)
                        rate = qa1 * qa2 * &
                        &(kernelpar(1)/((d3(int(epidat(i, 1)), int(epidat(j, 1)))**(2.0_c_double))+ &
                        & (kernelpar(1)**(2.0_c_double))))

                        gh(j) = ( minval( (/ tt, epidat(i, 4), epidat(j, 6) /) ) - &
                            & min(epidat(i, 6), epidat(j, 6)) ) * (rate)
                    else
                        rate    = 0.0_c_double
                        gh(j)   = 0.0_c_double
                    end if

                else
                    do m = 1,  nsuspar
                        suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                    end do
                    do m = 1,  ntranspar
                        transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                    end do
                    qa1 = dot_product(suspar, suscov1)
                    qa2 = dot_product(transpar, transcov1)
                    rate = qa1 * qa2 * &
                    &(kernelpar(1)/((d3(int(epidat(i, 1)), int(epidat(j, 1)))**(2.0_c_double))+ &
                    & (kernelpar(1)**(2.0_c_double))))
                    gh(j) = ( epidat(i, 4) - epidat(i, 6)) * (rate)
                end if
            end do
            df1(i) = sum(gh)
        end do

        df2 = 0.0_c_double

        do i =1, ninfected
            gh = 0.0_c_double
            do j = 1, n
                if (j .le. ninfected) then
                    if (j .ne. i) then
                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do
                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do
                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)
                        rate = gamma * (qa1 * qa2 * &
                            &(kernelpar(1)/((d3(int(epidat(i, 1)), int(epidat(j, 1)))**(2.0_c_double))+ &
                            & (kernelpar(1)**(2.0_c_double)))))
                        gh(j) = ( ( minval( (/ tt, epidat(i, 2), epidat(j, 6) /) ) - &
                            & min(epidat(i, 6), epidat(j, 6)) ) - &
                            & ( minval( (/ tt, epidat(i, 4), epidat(j, 6) /) ) - &
                            & min(epidat(i, 6), epidat(j, 6)) ) ) * (rate)
                    else
                        rate    = 0.0_c_double
                        gh(j)   = 0.0_c_double
                    end if
                else
                    do m = 1,  nsuspar
                        suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                    end do

                    do m = 1,  ntranspar
                        transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                    end do

                    qa1 = dot_product(suspar, suscov1)
                    qa2 = dot_product(transpar, transcov1)
                    rate = gamma * (qa1 * qa2 * &
                        &(kernelpar(1)/((d3(int(epidat(i, 1)), int(epidat(j, 1)))**(2.0_c_double))+ &
                        & (kernelpar(1)**(2.0_c_double)))))
                    gh(j) = ( (epidat(i, 2) - epidat(i, 6)) - (epidat(i, 4) - epidat(i, 6)) ) * (rate)
                end if
            end do
            df2(i) = sum(gh)
        end do

        likk1 = (-(sum(df1)+sum(df2)+ss))


! calculate the first part of the likelihood:

        rt(1) = 1.0_c_double

        do j = 2,  ninfected

            rrate = 0.0_c_double

            do i = 1,  ninfected

                if (i .ne. j) then

                    if ((epidat(j, 6) .gt. epidat(i, 6)) .and.  (epidat(j, 6) .le. epidat(i, 4)) ) then
                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do
                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do
                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)
                        rrate(i) = (qa1 * qa2 * &
                            &(kernelpar(1)/((d3(int(epidat(i, 1)), int(epidat(j, 1)))**(2.0_c_double))+ &
                            & (kernelpar(1)**(2.0_c_double)))))

                    else if ((epidat(j, 6) .gt. epidat(i, 4)) .and.  (epidat(j, 6) .le. epidat(i, 2)) ) then
                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do
                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do
                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)
                        rrate(i) = gamma * (qa1 * qa2 * &
                            &(kernelpar(1)/((d3(int(epidat(i, 1)), int(epidat(j, 1)))**(2.0_c_double))+ &
                            & (kernelpar(1)**(2.0_c_double)))))
                    else
                        rrate(i)   = 0.0_c_double
                    end if
                else
                    rrate(i)   = 0.0_c_double
                end if
            end do
            rt(j) = sum(rrate) + spark
        end do

        likk2 = sum(log(rt))


! TOTAL LOG-LIKELIHOOD:

        do i = 1,  ninfected
            gammain(i) = gammapdf ( deltain1,  deltain2,  epidat(i, 5) )
        end do
        do i = 1,  ninfected
            gammanr(i) = gammapdf ( deltanr1,  deltanr2,  epidat(i, 3) )
        end do

        likk = likk1 + likk2 + sum(gammain) + sum(gammanr)


    !##################################################################
    !##################################################################

    CASE (4)
! distance and network-based continuous ILMs (powerlaw distance kernel)

! getting t_obs (maximum removal time):
        tt = maxval(epidat(1:ninfected, 2))

!calculate the exponent part of the spark term:
        do r = 1,  n
            sss(r) = spark * (min(tt, epidat(r, 6)) - epidat(1, 6))
        end do
        ss = sum(sss)

! calculate the terms of the exponent part of the double summation:
        df1 = 0.0_c_double

        do i =1, ninfected

            gh = 0.0_c_double

            do j = 1, n

                if (j .le. ninfected) then

                    if (j .ne. i) then

                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do
                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do

                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)
                        rate = qa1 * qa2 * &
                            &((d3(int(epidat(i, 1)), int(epidat(j, 1)))**(-kernelpar(1))) + &
                            & (kernelpar(2)*cc(int(epidat(i, 1)), int(epidat(j, 1)))))
                        gh(j)   = ( minval( (/ tt, epidat(i, 4), epidat(j, 6) /) ) - &
                            & min(epidat(i, 6), epidat(j, 6)) ) * (rate)
                    else
                        rate    = 0.0_c_double
                        gh(j)   = 0.0_c_double
                    end if

                else

                    do m = 1,  nsuspar
                        suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                    end do

                    do m = 1,  ntranspar
                        transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                    end do

                    qa1 = dot_product(suspar, suscov1)
                    qa2 = dot_product(transpar, transcov1)
                    rate = qa1 * qa2 * &
                        &((d3(int(epidat(i, 1)), int(epidat(j, 1)))**(-kernelpar(1))) + &
                        & (kernelpar(2)*cc(int(epidat(i, 1)), int(epidat(j, 1)))))

                    gh(j)   = ( epidat(i, 4) - epidat(i, 6)) * (rate)

                end if

            end do

            df1(i) = sum(gh)

        end do

        df2 = 0.0_c_double

        do i =1, ninfected

            gh = 0.0_c_double

            do j = 1, n

                if (j .le. ninfected) then

                    if (j .ne. i) then

                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do

                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do

                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)

                        rate = gamma * (qa1 * qa2 * &
                            &((d3(int(epidat(i, 1)), int(epidat(j, 1)))**(-kernelpar(1))) + &
                            & (kernelpar(2)*cc(int(epidat(i, 1)), int(epidat(j, 1))))))

                        gh(j)   = ( ( minval( (/ tt, epidat(i, 2), epidat(j, 6) /) ) - &
                            & min(epidat(i, 6), epidat(j, 6)) ) - &
                            & ( minval( (/ tt, epidat(i, 4), epidat(j, 6) /) ) - &
                            & min(epidat(i, 6), epidat(j, 6)) ) ) * (rate)
                    else

                        rate    = 0.0_c_double
                        gh(j)   = 0.0_c_double

                    end if

                else

                    do m = 1,  nsuspar
                        suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                    end do

                    do m = 1,  ntranspar
                        transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                    end do

                    qa1 = dot_product(suspar, suscov1)
                    qa2 = dot_product(transpar, transcov1)

                    rate = gamma * (qa1 * qa2 * &
                        &((d3(int(epidat(i, 1)), int(epidat(j, 1)))**(-kernelpar(1))) + &
                        & (kernelpar(2)*cc(int(epidat(i, 1)), int(epidat(j, 1))))))

                    gh(j)   = ( (epidat(i, 2) - epidat(i, 6)) - (epidat(i, 4) - epidat(i, 6)) ) * (rate)

                end if

            end do

            df2(i) = sum(gh)

        end do

        likk1 = (-(sum(df1)+sum(df2)+ss))


! calculate the first part of the likelihood:

        rt(1) = 1.0_c_double

        do j = 2,  ninfected

            rrate = 0.0_c_double

            do i = 1,  ninfected

                if (i .ne. j) then

                    if ((epidat(j, 6) .gt. epidat(i, 6)) .and.  (epidat(j, 6) .le. epidat(i, 4)) ) then

                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do

                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do

                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)

                        rrate(i) = (qa1 * qa2 * &
                            &( (d3(int(epidat(i, 1)), int(epidat(j, 1)))**(-kernelpar(1))) + &
                            & (kernelpar(2)*cc(int(epidat(i, 1)), int(epidat(j, 1))))))

                    else if ((epidat(j, 6) .gt. epidat(i, 4)) .and.  (epidat(j, 6) .le. epidat(i, 2)) ) then

                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do

                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do

                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)

                        rrate(i) = gamma * (qa1 * qa2 * &
                            & ( (d3(int(epidat(i, 1)), int(epidat(j, 1)))**(-kernelpar(1))) + &
                            & (kernelpar(2)*cc(int(epidat(i, 1)), int(epidat(j, 1))))))

                    else

                        rrate(i)   = 0.0_c_double

                    end if

                else

                    rrate(i)   = 0.0_c_double

                end if

            end do

            rt(j) = sum(rrate) + spark

        end do

        likk2 = sum(log(rt))


! TOTAL LOG-LIKELIHOOD:

        do i = 1,  ninfected
            gammain(i) = gammapdf ( deltain1,  deltain2,  epidat(i, 5) )
        end do

        do i = 1,  ninfected
            gammanr(i) = gammapdf ( deltanr1,  deltanr2,  epidat(i, 3) )
        end do

        likk = likk1 + likk2 + sum(gammain) + sum(gammanr)

!##################################################################
!##################################################################

    CASE (5)
! distance and network-based continuous ILMs (Cauchy distance kernel)

! getting t_obs (maximum removal time):
        tt = maxval(epidat(1:ninfected, 2))

!calculate the exponent part of the spark term:
        do r = 1,  n
            sss(r) = spark * (min(tt, epidat(r, 6)) - epidat(1, 6))
        end do
        ss = sum(sss)

! calculate the terms of the exponent part of the double summation:
        df1 = 0.0_c_double

        do i =1, ninfected

            gh = 0.0_c_double

            do j = 1, n

                if (j .le. ninfected) then

                    if (j .ne. i) then

                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do

                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do

                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)

                        rate = qa1 * qa2 * &
                            & ((kernelpar(1)/((d3(int(epidat(i, 1)), int(epidat(j, 1)))**(2.0_c_double))+ &
                            & (kernelpar(1)**(2.0_c_double)))) + &
                            & (kernelpar(2)*cc(int(epidat(i, 1)), int(epidat(j, 1)))))

                        gh(j)   = ( minval( (/ tt, epidat(i, 4), epidat(j, 6) /) ) - &
                            & min(epidat(i, 6), epidat(j, 6)) ) * (rate)

                    else

                        rate    = 0.0_c_double
                        gh(j)   = 0.0_c_double

                    end if

                else

                    do m = 1,  nsuspar
                        suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                    end do

                    do m = 1,  ntranspar
                        transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                    end do

                    qa1 = dot_product(suspar, suscov1)

                    qa2 = dot_product(transpar, transcov1)

                    rate = qa1 * qa2 * &
                        & ((kernelpar(1)/((d3(int(epidat(i, 1)), int(epidat(j, 1)))**(2.0_c_double))+ &
                        & (kernelpar(1)**(2.0_c_double)))) + &
                        & (kernelpar(2)*cc(int(epidat(i, 1)), int(epidat(j, 1)))))

                    gh(j)   = ( epidat(i, 4) - epidat(i, 6)) * (rate)

                end if

            end do

            df1(i) = sum(gh)

        end do

        df2 = 0.0_c_double

        do i =1, ninfected

            gh = 0.0_c_double

            do j = 1, n

                if (j .le. ninfected) then

                    if (j .ne. i) then

                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do

                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do

                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)

                        rate = gamma * (qa1 * qa2 * &
                            & ((kernelpar(1)/((d3(int(epidat(i, 1)), int(epidat(j, 1)))**(2.0_c_double))+ &
                            & (kernelpar(1)**(2.0_c_double)))) + &
                            & (kernelpar(2)*cc(int(epidat(i, 1)), int(epidat(j, 1))))))

                        gh(j)   = ( ( minval( (/ tt, epidat(i, 2), epidat(j, 6) /) ) - &
                            & min(epidat(i, 6), epidat(j, 6)) ) - &
                            & ( minval( (/ tt, epidat(i, 4), epidat(j, 6) /) ) - &
                            & min(epidat(i, 6), epidat(j, 6)) ) ) * (rate)

                    else

                        rate    = 0.0_c_double
                        gh(j)   = 0.0_c_double

                    end if

                else

                    do m = 1,  nsuspar
                        suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                    end do

                    do m = 1,  ntranspar
                        transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                    end do

                    qa1 = dot_product(suspar, suscov1)

                    qa2 = dot_product(transpar, transcov1)

                    rate = gamma * (qa1 * qa2 * &
                        & ((kernelpar(1)/((d3(int(epidat(i, 1)), int(epidat(j, 1)))**(2.0_c_double))+ &
                        & (kernelpar(1)**(2.0_c_double)))) + &
                        & (kernelpar(2)*cc(int(epidat(i, 1)), int(epidat(j, 1))))))

                    gh(j) = ( (epidat(i, 2) - epidat(i, 6)) - (epidat(i, 4) - epidat(i, 6)) ) * (rate)

                end if

            end do

            df2(i) = sum(gh)

        end do

        likk1 = (-(sum(df1)+sum(df2)+ss))


! calculate the first part of the likelihood:

        rt(1) = 1.0_c_double

        do j = 2,  ninfected

            rrate = 0.0_c_double

            do i = 1,  ninfected

                if (i .ne. j) then

                    if ((epidat(j, 6) .gt. epidat(i, 6)) .and.  (epidat(j, 6) .le. epidat(i, 4)) ) then

                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do

                        do m = 1,  ntranspar
                        transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do

                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)

                        rrate(i) = (qa1 * qa2 * &
                            & ((kernelpar(1)/((d3(int(epidat(i, 1)), int(epidat(j, 1)))**(2.0_c_double))+ &
                            & (kernelpar(1)**(2.0_c_double)))) + &
                            & (kernelpar(2)*cc(int(epidat(i, 1)), int(epidat(j, 1))))))

                    else if ((epidat(j, 6) .gt. epidat(i, 4)) .and.  (epidat(j, 6) .le. epidat(i, 2)) ) then

                        do m = 1,  nsuspar
                            suscov1(m) = suscov(int(epidat(j, 1)), m)**powersus(m)
                        end do

                        do m = 1,  ntranspar
                            transcov1(m) = transcov(int(epidat(i, 1)), m)**powertrans(m)
                        end do

                        qa1 = dot_product(suspar, suscov1)
                        qa2 = dot_product(transpar, transcov1)

                        rrate(i) = gamma * (qa1 * qa2 * &
                            & ((kernelpar(1)/((d3(int(epidat(i, 1)), int(epidat(j, 1)))**(2.0_c_double))+ &
                            & (kernelpar(1)**(2.0_c_double)))) + &
                            & (kernelpar(2)*cc(int(epidat(i, 1)), int(epidat(j, 1))))))

                    else

                        rrate(i)   = 0.0_c_double

                    end if

                else

                    rrate(i)   = 0.0_c_double

                end if

            end do

            rt(j) = sum(rrate) + spark

        end do

        likk2 = sum(log(rt))


        ! TOTAL LOG-LIKELIHOOD:

        do i = 1,  ninfected
            gammain(i) = gammapdf ( deltain1,  deltain2,  epidat(i, 5) )
        end do

        do i = 1,  ninfected
            gammanr(i) = gammapdf ( deltanr1,  deltanr2,  epidat(i, 3) )
        end do

        likk = likk1 + likk2 + sum(gammain) + sum(gammanr)

    END SELECT


    end subroutine loglikcontilmsinr_f


!##################################################################
!##################################################################
! Gamma density function
    function gammapdf ( alph,  bet,  rval ) result(value)
    implicit none
    double precision :: alph, bet, valuelog, rval, value

        if ( rval <= 0.0d0 ) then
            value = 0.0d0
        else
            valuelog = gammalog ( alph )
            value = (alph * dlog ( bet )) + (( alph - 1.0d0 ) * dlog ( rval ) ) - &
                  & (bet * rval) - (valuelog)
        end if

    end function gammapdf


    double precision function gammalog ( x )

!*****************************************************************************
!
!  R8_GAMMA_LOG evaluates the logarithm of the gamma function.
!
!   DESCRIPTION:
!
!    This routine calculates the LOG(GAMMA) function for a positive real
!    argument X.  The code is distributed under the GNU LGPL license and
!    was written by by William Cody,  Laura Stoltz (FORTRAN77) and
!    John Burkardt (FORTRAN90). It is released in public domain and
!    modified for the use in this package.
!
!  Parameters:
!    Input,  double precisionX,  the argument of the function.
!    Output,  double precisionR8_GAMMA_LOG,  the value of the function.
!

    implicit none

    double precision,  dimension ( 7 ) :: c = (/ &
    -1.910444077728d-03,  &
    8.4171387781295d-04,  &
    -5.952379913043012d-04,  &
    7.93650793500350248d-04,  &
    -2.777777777777681622553d-03,  &
    8.333333333333333331554247d-02,  &
    5.7083835261d-03 /)
    double precision :: corr
    double precision :: d1 = -5.772156649015328605195174d-01
    double precision :: d2 = 4.227843350984671393993777d-01
    double precision :: d4 = 1.791759469228055000094023d0
    double precision,  parameter :: frtbig = 2.25d+76
    integer :: i
    double precision,  dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888d0,  &
    2.018112620856775083915565d+02,  &
    2.290838373831346393026739d+03,  &
    1.131967205903380828685045d+04,  &
    2.855724635671635335736389d+04,  &
    3.848496228443793359990269d+04,  &
    2.637748787624195437963534d+04,  &
    7.225813979700288197698961d+03 /)
    double precision,  dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064d+00,  &
    5.424138599891070494101986d+02,  &
    1.550693864978364947665077d+04,  &
    1.847932904445632425417223d+05,  &
    1.088204769468828767498470d+06,  &
    3.338152967987029735917223d+06,  &
    5.106661678927352456275255d+06,  &
    3.074109054850539556250927d+06 /)
    double precision,  dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062d+04,  &
    2.426813369486704502836312d+06,  &
    1.214755574045093227939592d+08,  &
    2.663432449630976949898078d+09,  &
    2.940378956634553899906876d+10,  &
    1.702665737765398868392998d+11,  &
    4.926125793377430887588120d+11,  &
    5.606251856223951465078242d+11 /)
    double precision,  dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036d+01,  &
    1.113332393857199323513008d+03,  &
    7.738757056935398733233834d+03,  &
    2.763987074403340708898585d+04,  &
    5.499310206226157329794414d+04,  &
    6.161122180066002127833352d+04,  &
    3.635127591501940507276287d+04,  &
    8.785536302431013170870835d+03 /)
    double precision,  dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942d+02,  &
    7.765049321445005871323047d+03,  &
    1.331903827966074194402448d+05,  &
    1.136705821321969608938755d+06,  &
    5.267964117437946917577538d+06,  &
    1.346701454311101692290052d+07,  &
    1.782736530353274213975932d+07,  &
    9.533095591844353613395747d+06 /)
    double precision,  dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843d+03,  &
    6.393885654300092398984238d+05,  &
    4.135599930241388052042842d+07,  &
    1.120872109616147941376570d+09,  &
    1.488613728678813811542398d+10,  &
    1.016803586272438228077304d+11,  &
    3.417476345507377132798597d+11,  &
    4.463158187419713286462081d+11 /)
    double precision :: res
    double precision,  parameter :: sqrtpi = 0.9189385332046727417803297d0
    double precision :: x
    double precision,  parameter :: xbig = 2.55d+305
    double precision :: xden
    double precision,  parameter :: xinf = 1.79d+308
    double precision :: xm1
    double precision :: xm2
    double precision :: xm4
    double precision :: xnum
    double precision :: y
    double precision :: ysq

        y = x
        if ( 0.0d0 < y .and. y <= xbig ) then
            if ( y <= epsilon ( y ) ) then
                res = - log ( y )
!
!  EPS < X <= 1.5.
!
            else if ( y <= 1.5d0 ) then
                if ( y < 0.6796875d0 ) then
                    corr = -log ( y )
                    xm1 = y
                else
                    corr = 0.0d0
                    xm1 = ( y - 0.5d0 ) - 0.5d0
                end if

                if ( y <= 0.5d0 .or. 0.6796875d0 <= y ) then
                    xden = 1.0d0
                    xnum = 0.0d0
                    do i = 1,  8
                        xnum = xnum * xm1 + p1(i)
                        xden = xden * xm1 + q1(i)
                    end do
                    res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )
                else
                    xm2 = ( y - 0.5d0 ) - 0.5d0
                    xden = 1.0d0
                    xnum = 0.0d0
                    do i = 1,  8
                        xnum = xnum * xm2 + p2(i)
                        xden = xden * xm2 + q2(i)
                    end do
                    res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )
                end if
!
!  1.5 < X <= 4.0.
!
            else if ( y <= 4.0d0 ) then
                xm2 = y - 2.0d0
                xden = 1.0d0
                xnum = 0.0d0
                do i = 1,  8
                    xnum = xnum * xm2 + p2(i)
                    xden = xden * xm2 + q2(i)
                end do
                res = xm2 * ( d2 + xm2 * ( xnum / xden ) )
!
!  4.0 < X <= 12.0.
!
            else if ( y <= 12.0d0 ) then
                xm4 = y - 4.0d0
                xden = -1.0d0
                xnum = 0.0d0
                do i = 1,  8
                    xnum = xnum * xm4 + p4(i)
                    xden = xden * xm4 + q4(i)
                end do
                res = d4 + xm4 * ( xnum / xden )
!
!  Evaluate for 12 <= argument.
!
            else
                res = 0.0d0
                if ( y <= frtbig ) then
                    res = c(7)
                    ysq = y * y
                    do i = 1,  6
                        res = res / ysq + c(i)
                    end do
                end if
                res = res / y
                corr = log ( y )
                res = res + sqrtpi - 0.5d0 * corr
                res = res + y * ( corr - 1.0d0 )
            end if
        else
            res = xinf
        end if

        gammalog = res

    end function gammalog


end module logliksinr
