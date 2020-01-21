!######################################################################
!# MODULE: sinr
!# AUTHORS:
!#     Waleed Almutiry <walmutir@uoguelph.ca>,
!#     Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca> and
!#     Rob Deardon <robert.deardon@ucalgary.ca>
!#
!# DESCRIPTION:
!#
!#     MCMC tools for SINR continuous-time ILMs:
!#
!#     This program is free software; you can redistribute it and/or
!#     modify it under the terms of the GNU General Public License,
!#     version 3, as published by the Free Software Foundation.
!#
!#     This program is distributed in the hope that it will be useful,
!#     but without any warranty; without even the implied warranty of
!#     merchantability or fitness for a particular purpose.  See the GNU
!#     General Public License, version 3, for more details.
!#
!#     A copy of the GNU General Public License, version 3, is available
!#     at http://www.r-project.org/Licenses/GPL-3
!#
!# Part of the R/EpiILMCT package
!# Contains:
!#           mcmcsinr ..................... subroutine
!#           likelihoodsinrsinr2 ........... subroutine
!#           halfnormalden2 ................ function
!#           gammadensity2 ................. function
!#           uniformden2 ................... function
!#           randnormal2 ................... function
!#           randgamma2 .................... function
!#           initrandomseedsinr2 ........... subroutine
!#           ransamsinr2 ................... subroutine
!#
!######################################################################

module sinr2
    use, intrinsic :: iso_c_binding
    implicit none
    public :: mcmcsinr_f

contains

subroutine mcmcsinr_f(n, nsim, ni, num, anum2, nsuspar, ntranspar, net, dis, epidat, blockupdate, &
& priordistsuspar, priordisttranspar, priordistkernelparpar, &
& priordistpowersus, priordistpowertrans, priordistsparkpar, priordistgammapar,&
& suspar, suscov, powersus, transpar, transcov, powertrans, kernelpar, spark, gamma, &
& deltain, deltanr, &
& kernelparproposalvar, sparkproposalvar, gammaproposalvar, &
& susproposalvar, powersusproposalvar, transproposalvar, &
& powertransproposalvar, infperiodproposalin, infperiodproposalnr, &
& priorpar1sus, priorpar2sus, &
& priorpar1powersus, priorpar2powersus, &
& priorpar1trans, priorpar2trans, &
& priorpar1powertrans, priorpar2powertrans, &
& kernelparprior, sparkprior, gammaprior, &
& deltain2prior, deltanr2prior, susparop, powersusparop, transparop, powertransparop, &
& kernelparop, sparkop, gammaop, deltain2op, deltanr2op, &
& epidatmctim, epidatmcrem, loglik) bind(C, name="mcmcsinr_f_")

external infinity_value
external seedin
external seedout
external randomnumber

!Declarations of input variables:
integer (C_INT), intent(in):: n, nsim, ni, nsuspar, ntranspar, num
! epidemic data:
real (C_DOUBLE), intent(in), dimension(n, 6):: epidat
real (C_DOUBLE), dimension(n, 6):: epidat3current, epidat3update, epidat3
! distance and network matrices:
real (C_DOUBLE), intent(in), dimension(n, n):: dis, net
! susceptibility and transmissibility covariates:
real (C_DOUBLE), intent(in), dimension(n, nsuspar):: suscov
real (C_DOUBLE), intent(in), dimension(n, ntranspar):: transcov
! to inndicate which parameter to update:
integer (C_INT), intent(in), dimension(8):: anum2
! prior distributions:
integer (C_INT), intent(in), dimension(ntranspar):: priordisttranspar, priordistpowertrans
integer (C_INT), intent(in), dimension(nsuspar):: priordistsuspar, priordistpowersus
integer (C_INT), intent(in), dimension(2):: priordistkernelparpar
integer (C_INT), intent(in):: priordistsparkpar, priordistgammapar
! for updating infection times/removal times/infectious periods in blocks:
integer (C_INT), intent(in), dimension(2):: blockupdate
integer (C_INT):: obs_inf_time, sizeblock
! The variance of the normal proposal distributions and the parameters of the prior distributions of the model parameters:
real (C_DOUBLE), intent(in), dimension(ntranspar):: transproposalvar, priorpar1trans, priorpar2trans
real (C_DOUBLE), intent(in), dimension(ntranspar):: powertransproposalvar, priorpar1powertrans, priorpar2powertrans
real (C_DOUBLE), intent(in), dimension(nsuspar):: susproposalvar, priorpar1sus, priorpar2sus
real (C_DOUBLE), intent(in), dimension(nsuspar):: powersusproposalvar, priorpar1powersus, priorpar2powersus
real (C_DOUBLE), intent(in), dimension(2):: kernelparproposalvar
real (C_DOUBLE), intent(in):: sparkproposalvar, gammaproposalvar
real (C_DOUBLE), intent(in), dimension(2, 2):: kernelparprior
real (C_DOUBLE), intent(in), dimension(2):: sparkprior, gammaprior
! parameters of the prior distribution of the incubation and delay period rates
! and the parameters of the gamma proposal distributions of the independenc sampler:
real (C_DOUBLE), intent(in), dimension(2):: deltain2prior, deltanr2prior, infperiodproposalin, infperiodproposalnr
! initial values
real (C_DOUBLE), intent(in), dimension(ntranspar):: transpar, powertrans
real (C_DOUBLE), intent(in), dimension(nsuspar):: suspar, powersus
real (C_DOUBLE), intent(in), dimension(2):: kernelpar, deltain, deltanr
real (C_DOUBLE), intent(in):: spark, gamma

! OUTPUT:
real (C_DOUBLE), dimension(nsim, 1), intent(out):: deltain2op, deltanr2op
real (C_DOUBLE), dimension(nsim, n), intent(out):: epidatmctim, epidatmcrem
real (C_DOUBLE), dimension(nsim, 1), intent(out):: loglik
real (C_DOUBLE), dimension(nsim, 2), intent(out):: kernelparop
real (C_DOUBLE), dimension(nsim, 1), intent(out):: sparkop, gammaop
real (C_DOUBLE), dimension(nsim, ntranspar), intent(out) :: transparop, powertransparop
real (C_DOUBLE), dimension(nsim, nsuspar), intent(out) :: susparop, powersusparop

!Declarations of variables used in the subroutine:
real (C_DOUBLE), dimension(2):: infperiodproposalin1, infperiodproposalnr1
real (C_DOUBLE), dimension(nsim, 1):: deltain1op, deltanr1op
real (C_DOUBLE), dimension(nsim, n)::epidatmcperin, epidatmcpernr
real (C_DOUBLE), dimension(nsuspar):: postparsusupdate, postparsuscurrent
real (C_DOUBLE), dimension(ntranspar):: postpartransupdate, postpartranscurrent
real (C_DOUBLE), dimension(2):: postkernelparupdate, postkernelparcurrent
real (C_DOUBLE):: postsparkupdate, postgammaupdate
real (C_DOUBLE):: ratio1, ratio2, ratio3, ratio4, ratio5, ratio6
real (C_DOUBLE):: psi1, psi2, psi3, psi4, psi5, psi6
real (C_DOUBLE):: u1, u2, u3, u4, u5, u6
real (C_DOUBLE):: valueupdate1, valuecurrent1
real (C_DOUBLE):: valueupdate2, valuecurrent2
real (C_DOUBLE):: valueupdate3, valuecurrent3
real (C_DOUBLE):: valueupdate4, valuecurrent4
real (C_DOUBLE):: valueupdate5, valuecurrent5
real (C_DOUBLE):: valueupdate6, valuecurrent6
real (C_DOUBLE):: valueffupdate, valueffcurrent
real (C_DOUBLE):: valuefupdate1, valuefcurrent1
real (C_DOUBLE):: valuefupdate2, valuefcurrent2
real (C_DOUBLE):: valuefupdate3, valuefcurrent3
real (C_DOUBLE):: valuefupdate4, valuefcurrent4
real (C_DOUBLE):: valuefupdate5, valuefcurrent5
real (C_DOUBLE):: valuefupdate6, valuefcurrent6, denupdate, dencurrent
real (C_DOUBLE):: zparsus, zpartrans, zkernelpar, zspark, zgamma
integer (C_INT):: i, j, r, mn, m, nblock
integer (C_INT):: abs32, abs33, abs34, abs35
real (C_DOUBLE):: Inf, loglik66, delaydens, incdens
! for getting different random numbers for each block at each iteration:
integer (C_INT), allocatable, dimension(:, :):: amcmc
integer (C_INT), allocatable, dimension(:):: amcmc1
integer (C_INT), dimension(ni):: it, xxindic


call infinity_value(Inf)
call seedin()

obs_inf_time = blockupdate(1)
sizeblock    = blockupdate(2)

valuefupdate1 = 0.0_c_double
valuefcurrent1 = 0.0_c_double
valueupdate1 = 0.0_c_double
valuecurrent1 = 0.0_c_double

valuefupdate2 = 0.0_c_double
valuefcurrent2 = 0.0_c_double
valueupdate2 = 0.0_c_double
valuecurrent2 = 0.0_c_double

valuefupdate3 = 0.0_c_double
valuefcurrent3 = 0.0_c_double
valueupdate3 = 0.0_c_double
valuecurrent3 = 0.0_c_double

valuefupdate4 = 0.0_c_double
valuefcurrent4 = 0.0_c_double
valueupdate4 = 0.0_c_double
valuecurrent4 = 0.0_c_double

valuefupdate5 = 0.0_c_double
valuefcurrent5 = 0.0_c_double
valueupdate5 = 0.0_c_double
valuecurrent5 = 0.0_c_double

valuefupdate6 = 0.0_c_double
valuefcurrent6 = 0.0_c_double
valueupdate6 = 0.0_c_double
valuecurrent6 = 0.0_c_double

ratio1 = 0.0_c_double
psi1 = 0.0_c_double
ratio2 = 0.0_c_double
psi2 = 0.0_c_double
ratio3 = 0.0_c_double
psi3 = 0.0_c_double
ratio4 = 0.0_c_double
psi4 = 0.0_c_double
ratio5 = 0.0_c_double
psi5 = 0.0_c_double
ratio6 = 0.0_c_double
psi6 = 0.0_c_double

valueffupdate = 0.0_c_double
valueffcurrent = 0.0_c_double

! Initial values for the MCMC updates:

deltain1op(1, 1) = deltain(1)

deltain2op(1, 1) = deltain(2)

deltanr1op(1, 1) = deltanr(1)

deltanr2op(1, 1) = deltanr(2)

epidat3 = epidat

epidatmcperin(1, :) = epidat3(:, 5)
epidatmctim(1, :)   = epidat3(:, 6)
epidatmcpernr(1, :) = epidat3(:, 3)
epidatmcrem(1, :)   = epidat3(:, 2)

kernelparop(1, :)               = kernelpar
sparkop(1, 1)                   = spark
transparop(1, 1:ntranspar)      = transpar
powertransparop(1, 1:ntranspar) = powertrans
susparop(1, 1:nsuspar)          = suspar
powersusparop(1, 1:nsuspar)     = powersus
gammaop(1, 1)                   = gamma


! calculate the loglikelihood with the initial values:

    if (anum2(6) .eq. 1) then

        incdens = 0.0_c_double
        do i = 1, ni
            if (i .gt. obs_inf_time) then
                incdens = incdens + gammadensity2(epidat3(i, 5), deltain1op(1, 1), deltain2op(1, 1))
            else
                incdens = incdens
            end if
        end do

        call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
        & suscov, transcov, susparop(1, :), powersusparop(1, :), &
        & transparop(1, :), powertransparop(1, :), kernelparop(1, :), &
        & sparkop(1, 1), gammaop(1, 1), loglik66)

        loglik(1, 1) = loglik66 + incdens

    else if (anum2(6) .eq. 2) then

        incdens = 0.0_c_double
        delaydens = 0.0_c_double
        do i = 1, ni
            if (i .gt. obs_inf_time) then
                incdens = incdens + gammadensity2(epidat3(i, 5), deltain1op(1, 1), deltain2op(1, 1))
            else
                incdens = incdens
            end if
            delaydens = delaydens + gammadensity2(epidat3(i, 3), deltanr1op(1, 1), deltanr2op(1, 1))
        end do

        call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
        & suscov, transcov, susparop(1, :), powersusparop(1, :), &
        & transparop(1, :), powertransparop(1, :), kernelparop(1, :), &
        & sparkop(1, 1), gammaop(1, 1), loglik66)

        loglik(1, 1) = loglik66 + incdens + delaydens

    else if (anum2(6) .eq. 3) then

        call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
        & suscov, transcov, susparop(1, :), powersusparop(1, :), &
        & transparop(1, :), powertransparop(1, :), kernelparop(1, :), &
        & sparkop(1, 1), gammaop(1, 1), loglik66)

        loglik(1, 1) = loglik66

    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MCMC start
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do j = 1 , (nsim-1)

        call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
        & suscov, transcov, susparop(j, :), powersusparop(j, :), &
        & transparop(j, :), powertransparop(j, :), kernelparop(j, :), &
        & sparkop(j, 1), gammaop(j, 1), valuecurrent1)

    if (anum2(6) .eq. 1) then

    ! if the infection times and incubation periods are assumed unknown:

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: incubation period rate (deltain2) and delay period rate (deltanr2):
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        deltain1op(j+1, 1) = deltain1op(j, 1)
        deltain2op(j+1, 1) = randgamma2(((dble(ni)*deltain1op(j+1, 1))+deltain2prior(1)), &
        & 1.0_c_double/(deltain2prior(2)+sum(epidat3(1:ni, 5))))

        deltanr1op(j+1, 1) = deltanr1op(j, 1)
        deltanr2op(j+1, 1) = deltanr2op(j, 1)

        if (infperiodproposalin(1) .eq. 0.0_c_double) then
            infperiodproposalin1 = (/ deltain1op(j+1, 1), deltain2op(j+1, 1) /)
        else
            infperiodproposalin1 = infperiodproposalin
        end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: infection times/incubation period & removal times/delay periods:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        nblock = floor(real(ni)/real(sizeblock))

        it = (/(r, r = 1 , ni)/)

        xxindic = 0

        epidat3update  = epidat3      ! proposing
        epidat3current = epidat3      ! current

        ! if the sizeblock has a value greater than 1, then update infection times/incubation periods
        ! in blocks with assigning those in each block randomly at each iteration.

        allocate(amcmc(sizeblock, nblock))
        allocate(amcmc1(ni-(sizeblock*nblock)))

        amcmc  = 0
        amcmc1 = 0

        do r = 1 , nblock
            call ransamsinr2(pack(it, (xxindic .ne. 1)), amcmc(:, r), count(xxindic .ne. 1), sizeblock)
            xxindic(amcmc(:, r)) = 1
        end do

        amcmc1 = pack(it, (xxindic .ne. 1))

        mn = (ni-obs_inf_time)-(sizeblock*nblock)

        call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
        & suscov, transcov, susparop(j, :), powersusparop(j, :), &
        & transparop(j, :), powertransparop(j, :), kernelparop(j, :), &
        & sparkop(j, 1), gammaop(j, 1), valuecurrent1)

    do r = 1 , nblock

        do m = 1 , sizeblock
            if (amcmc(m, r) .gt. obs_inf_time) then
                epidat3update(amcmc(m, r), 5)  = randgamma2(infperiodproposalin1(1), &
                & 1.0_c_double/infperiodproposalin1(2))
                epidat3update(amcmc(m, r), 6)  = epidat3update(amcmc(m, r), 4) - &
                & epidat3update(amcmc(m, r), 5)
            else
                epidat3update(amcmc(m, r), 5)  = epidat3current(amcmc(m, r), 5)
                epidat3update(amcmc(m, r), 6)  = epidat3current(amcmc(m, r), 6)
            end if
        end do

        !calculating the log-likelihood function:
        call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3update, &
        & suscov, transcov, susparop(j, :), powersusparop(j, :), &
        & transparop(j, :), powertransparop(j, :), kernelparop(j, :), &
        & sparkop(j, 1), gammaop(j, 1), valueupdate1)

        denupdate      = 0.0_c_double
        dencurrent     = 0.0_c_double
        valuefupdate1  = 0.0_c_double
        valuefcurrent1 = 0.0_c_double

        do m = 1 , sizeblock
            if (amcmc(m, r) .gt. obs_inf_time) then
                denupdate = denupdate + &
                & gammadensity2(epidat3update(amcmc(m, r), 5), deltain1op(j+1, 1), deltain2op(j+1, 1))

                dencurrent = dencurrent + &
                & gammadensity2(epidat3current(amcmc(m, r), 5), deltain1op(j+1, 1), deltain2op(j+1, 1))

                valuefupdate1  = valuefupdate1  + &
                & gammadensity2(epidat3update(amcmc(m, r), 5), infperiodproposalin1(1), infperiodproposalin1(2))

                valuefcurrent1 = valuefcurrent1 + &
                & gammadensity2(epidat3current(amcmc(m, r), 5), infperiodproposalin1(1), infperiodproposalin1(2))
            else
                denupdate = denupdate
                dencurrent = dencurrent
                valuefupdate1  = valuefupdate1
                valuefcurrent1 = valuefcurrent1
            end if
        end do

        valueffupdate  = valueupdate1  + denupdate
        valueffcurrent = valuecurrent1 + dencurrent

        !calculating the log acceptance ratio:
        ratio1    = (valueffupdate+valuefcurrent1)-(valueffcurrent+valuefupdate1)

        psi1      = min(1.0_c_double, dexp(ratio1))

        ! decision: accept/reject the proposed value:
        call randomnumber(u1)

        if ( (u1 .le. psi1) .and. (any(epidat3update(amcmc(:, r), 6) .lt. &
        & epidat3current(obs_inf_time, 6)).eqv. .false.) ) then

            epidat3current(amcmc(:, r), 5) = epidat3update(amcmc(:, r), 5)
            epidat3current(amcmc(:, r), 6) = epidat3update(amcmc(:, r), 6)
                valuecurrent1 = valueupdate1

        else

            epidat3update(amcmc(:, r), 5)  = epidat3current(amcmc(:, r), 5)
            epidat3update(amcmc(:, r), 6)  = epidat3current(amcmc(:, r), 6)
                valuecurrent1 = valuecurrent1

        end if

    end do

    ! in case the last block has less number of infection times than the sizeblock:

    if (mn .ne. 0)then

        do m = 1 , mn
            if (amcmc1(m) .gt. obs_inf_time) then
                epidat3update(amcmc1(m), 5)  = randgamma2(infperiodproposalin1(1), &
                & 1.0_c_double/infperiodproposalin1(2))
                epidat3update(amcmc1(m), 6)  = epidat3update(amcmc1(m), 4) - &
                & epidat3update(amcmc1(m), 5)
            else
                epidat3update(amcmc1(m), 5)  = epidat3current(amcmc1(m), 5)
                epidat3update(amcmc1(m), 6)  = epidat3current(amcmc1(m), 6)
            end if
        end do

        !calculating the log-likelihood function:
        call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3update, &
        & suscov, transcov, susparop(j, :), powersusparop(j, :), &
        & transparop(j, :), powertransparop(j, :), kernelparop(j, :), &
        & sparkop(j, 1), gammaop(j, 1), valueupdate1)

        denupdate      = 0.0_c_double
        dencurrent     = 0.0_c_double
        valuefupdate1  = 0.0_c_double
        valuefcurrent1 = 0.0_c_double

        do m = 1 , mn
            if (amcmc1(m) .gt. obs_inf_time) then
                denupdate  = denupdate  + &
                & gammadensity2(epidat3update(amcmc1(m), 5), deltain1op(j+1, 1), deltain2op(j+1, 1))
                dencurrent = dencurrent + &
                & gammadensity2(epidat3current(amcmc1(m), 5), deltain1op(j+1, 1), deltain2op(j+1, 1))
                valuefupdate1  = valuefupdate1 + &
                & gammadensity2(epidat3update(amcmc1(m), 5), infperiodproposalin1(1), infperiodproposalin1(2))
                valuefcurrent1 = valuefcurrent1+ &
                & gammadensity2(epidat3current(amcmc1(m), 5), infperiodproposalin1(1), infperiodproposalin1(2))
            else
                denupdate  = denupdate
                dencurrent = dencurrent
                valuefupdate1  = valuefupdate1
                valuefcurrent1 = valuefcurrent1
            end if
        end do

        valueffupdate  = valueupdate1  + denupdate
        valueffcurrent = valuecurrent1 + dencurrent

        !calculating the log acceptance ratio:
        ratio1    = (valueffupdate+valuefcurrent1)-(valueffcurrent+valuefupdate1)

        psi1      = min(1.0_c_double, dexp(ratio1))

        ! decision: accept/reject the proposed value:
        call randomnumber(u1)

        if ( (u1 .le. psi1) .and. (any(epidat3update(amcmc1, 6) .lt.&
        & epidat3current(obs_inf_time, 6)).eqv. .false.) ) then
            epidat3current(amcmc1, 5) = epidat3update(amcmc1, 5)
            epidat3current(amcmc1, 6) = epidat3update(amcmc1, 6)
                valuecurrent1 = valueupdate1
        else
            epidat3update(amcmc1, 5)  = epidat3current(amcmc1, 5)
            epidat3update(amcmc1, 6)  = epidat3current(amcmc1, 6)
                valuecurrent1 = valuecurrent1
        end if

    end if

    epidat3 = epidat3current

    deallocate(amcmc)
    deallocate(amcmc1)

    epidatmcperin(j+1, :)   = epidat3(:, 5)
    epidatmctim(j+1, :)     = epidat3(:, 6)
    epidatmcpernr(j+1, :) = epidat3(:, 3)
    epidatmcrem(j+1, :)   = epidat3(:, 2)


    else if (anum2(6) .eq. 2) then

! if the infection times, incubation periods, removal times and delay periods are assumed unknown:

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: incubation period rate (deltain2) & delay period rate (deltanr2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        deltain1op(j+1, 1) = deltain1op(j, 1)

        deltain2op(j+1, 1) = randgamma2(((dble(ni)*deltain1op(j+1, 1))+deltain2prior(1)), &
        & (1.0_c_double/(deltain2prior(2)+sum(epidat3(1:ni, 5)))))

        deltanr1op(j+1, 1) = deltanr1op(j, 1)

        deltanr2op(j+1, 1) = randgamma2(((dble(ni)*deltanr1op(j+1, 1))+deltanr2prior(1)), &
        & (1.0_c_double/(deltanr2prior(2)+sum(epidat3(1:ni, 3)))))

        if (infperiodproposalin(1) .eq. 0.0_c_double) then
            infperiodproposalin1 = (/ deltain1op(j+1, 1), deltain2op(j+1, 1) /)
            infperiodproposalnr1 = (/ deltanr1op(j+1, 1), deltanr2op(j+1, 1) /)
        else
            infperiodproposalin1 = infperiodproposalin
            infperiodproposalnr1 = infperiodproposalnr
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: infection times / incubation period  & removal times / delay period :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        nblock = floor(real(ni)/real(sizeblock))

        it = (/(r, r = 1 , ni)/)

        xxindic = 0

        epidat3update  = epidat3      ! proposing
        epidat3current = epidat3      ! current

! if the sizeblock has a value greater than 1, then update infection times/incubation periods
! and removal times/delay periods in blocks with assigning those in each block randomly
! at each iteration.
! note that: the code is to assign random number for each block that represent event number
! in the epidemic, such as 5 for the fifth infected individual in the epidemic. Thus, the
! infection and removal times of those numbers are also updated together.

        allocate(amcmc(sizeblock, nblock))
        allocate(amcmc1(ni-(sizeblock*nblock)))

        amcmc  = 0
        amcmc1 = 0

        do r = 1 , nblock
            call ransamsinr2(pack(it, (xxindic .ne. 1)), amcmc(:, r), count(xxindic .ne. 1), sizeblock)
            xxindic(amcmc(:, r)) = 1
        end do

        amcmc1 = pack(it, (xxindic .ne. 1))

        mn = ni-(sizeblock*nblock)

            call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
            & suscov, transcov, susparop(j, :), powersusparop(j, :), &
            & transparop(j, :), powertransparop(j, :), kernelparop(j, :), &
            & sparkop(j, 1), gammaop(j, 1), valuecurrent1)

        do r = 1 , nblock

            do m = 1 , sizeblock
                if (amcmc(m, r) .gt. obs_inf_time) then
                    epidat3update(amcmc(m, r), 5)  = randgamma2(infperiodproposalin1(1), &
                    & 1.0_c_double/infperiodproposalin1(2))
                    epidat3update(amcmc(m, r), 6)  = epidat3update(amcmc(m, r), 4) - &
                    & epidat3update(amcmc(m, r), 5)
                    epidat3update(amcmc(m, r), 3)  = randgamma2(infperiodproposalnr1(1), &
                    & 1.0_c_double/infperiodproposalnr1(2))
                    epidat3update(amcmc(m, r), 2)  = epidat3update(amcmc(m, r), 4) + &
                    & epidat3update(amcmc(m, r), 3)
                else
                    epidat3update(amcmc(m, r), 5)  = epidat3current(amcmc(m, r), 5)
                    epidat3update(amcmc(m, r), 6)  = epidat3current(amcmc(m, r), 6)
                    epidat3update(amcmc(m, r), 3)  = randgamma2(infperiodproposalnr1(1), &
                    & 1.0_c_double/infperiodproposalnr1(2))
                    epidat3update(amcmc(m, r), 2)  = epidat3update(amcmc(m, r), 4) + &
                    & epidat3update(amcmc(m, r), 3)
                end if
            end do

            !calculating the log-likelihood function:
            call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3update, &
            & suscov, transcov, susparop(j, :), powersusparop(j, :), &
            & transparop(j, :), powertransparop(j, :), kernelparop(j, :), &
            & sparkop(j, 1), gammaop(j, 1), valueupdate1)

            denupdate      = 0.0_c_double
            dencurrent     = 0.0_c_double
            valuefupdate1  = 0.0_c_double
            valuefcurrent1 = 0.0_c_double

            do m = 1 , sizeblock
                if (amcmc(m, r) .gt. obs_inf_time) then
                    denupdate = denupdate + &
                    & gammadensity2(epidat3update(amcmc(m, r), 5), deltain1op(j+1, 1), deltain2op(j+1, 1)) + &
                    & gammadensity2(epidat3update(amcmc(m, r), 3), deltanr1op(j+1, 1), deltanr2op(j+1, 1))

                    dencurrent = dencurrent + &
                    & gammadensity2(epidat3current(amcmc(m, r), 5), deltain1op(j+1, 1), deltain2op(j+1, 1)) + &
                    & gammadensity2(epidat3current(amcmc(m, r), 3), deltanr1op(j+1, 1), deltanr2op(j+1, 1))

                    valuefupdate1  = valuefupdate1  + &
                    & gammadensity2(epidat3update(amcmc(m, r), 5), infperiodproposalin1(1), infperiodproposalin1(2)) + &
                    & gammadensity2(epidat3update(amcmc(m, r), 3), infperiodproposalnr1(1), infperiodproposalnr1(2))

                    valuefcurrent1 = valuefcurrent1 + &
                    & gammadensity2(epidat3current(amcmc(m, r), 5), infperiodproposalin1(1), infperiodproposalin1(2)) + &
                    & gammadensity2(epidat3current(amcmc(m, r), 3), infperiodproposalnr1(1), infperiodproposalnr1(2))
                else
                    denupdate = denupdate + &
                    & gammadensity2(epidat3update(amcmc(m, r), 3), deltanr1op(j+1, 1), deltanr2op(j+1, 1))

                    dencurrent = dencurrent + &
                    & gammadensity2(epidat3current(amcmc(m, r), 3), deltanr1op(j+1, 1), deltanr2op(j+1, 1))

                    valuefupdate1  = valuefupdate1  + &
                    & gammadensity2(epidat3update(amcmc(m, r), 3), infperiodproposalnr1(1), infperiodproposalnr1(2))

                    valuefcurrent1 = valuefcurrent1 + &
                    & gammadensity2(epidat3current(amcmc(m, r), 3), infperiodproposalnr1(1), infperiodproposalnr1(2))
                end if
            end do

            valueffupdate  = valueupdate1  + denupdate
            valueffcurrent = valuecurrent1 + dencurrent

            !calculating the log acceptance ratio:
            ratio1    = (valueffupdate+valuefcurrent1)-(valueffcurrent+valuefupdate1)

            psi1      = min(1.0_c_double, dexp(ratio1))

            ! decision: accept/reject the proposed value:
            call randomnumber(u1)

            if ( (u1 .le. psi1) .and. (any(epidat3update(amcmc(:, r), 6) .lt. &
            & epidat3current(obs_inf_time, 6)).eqv. .false.) ) then

                epidat3current(amcmc(:, r), 2) = epidat3update(amcmc(:, r), 2)
                epidat3current(amcmc(:, r), 3) = epidat3update(amcmc(:, r), 3)
                epidat3current(amcmc(:, r), 5) = epidat3update(amcmc(:, r), 5)
                epidat3current(amcmc(:, r), 6) = epidat3update(amcmc(:, r), 6)
                valuecurrent1 = valueupdate1

            else

                epidat3update(amcmc(:, r), 2)  = epidat3current(amcmc(:, r), 2)
                epidat3update(amcmc(:, r), 3)  = epidat3current(amcmc(:, r), 3)
                epidat3update(amcmc(:, r), 5)  = epidat3current(amcmc(:, r), 5)
                epidat3update(amcmc(:, r), 6)  = epidat3current(amcmc(:, r), 6)
                valuecurrent1 = valuecurrent1

            end if

        end do

    ! in case the last block has less number of infection times than the sizeblock:

        if (mn .ne. 0)then

            do m = 1 , mn
                if (amcmc1(m) .gt. obs_inf_time) then
                    epidat3update(amcmc1(m), 5)  = randgamma2(infperiodproposalin1(1), &
                    & 1.0_c_double/infperiodproposalin1(2))
                    epidat3update(amcmc1(m), 6)  = epidat3update(amcmc1(m), 4) - &
                    & epidat3update(amcmc1(m), 5)
                    epidat3update(amcmc1(m), 3)  = randgamma2(infperiodproposalnr1(1), &
                    & 1.0_c_double/infperiodproposalnr1(2))
                    epidat3update(amcmc1(m), 2)  = epidat3update(amcmc1(m), 4) + &
                    & epidat3update(amcmc1(m), 3)
                else
                    epidat3update(amcmc1(m), 5)  = epidat3current(amcmc1(m), 5)
                    epidat3update(amcmc1(m), 6)  = epidat3current(amcmc1(m), 6)
                    epidat3update(amcmc1(m), 3)  = randgamma2(infperiodproposalnr1(1), &
                    & 1.0_c_double/infperiodproposalnr1(2))
                    epidat3update(amcmc1(m), 2)  = epidat3update(amcmc1(m), 4) + &
                    & epidat3update(amcmc1(m), 3)
                end if
            end do

            !calculating the log-likelihood function:
            call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3update, &
            & suscov, transcov, susparop(j, :), powersusparop(j, :), &
            & transparop(j, :), powertransparop(j, :), kernelparop(j, :), sparkop(j, 1), gammaop(j, 1), valueupdate1)

            denupdate      = 0.0_c_double
            dencurrent     = 0.0_c_double
            valuefupdate1  = 0.0_c_double
            valuefcurrent1 = 0.0_c_double

            do m = 1 , mn
                if (amcmc1(m) .gt. obs_inf_time) then
                    denupdate  = denupdate  + &
                    & gammadensity2(epidat3update(amcmc1(m), 5), deltain1op(j+1, 1), deltain2op(j+1, 1)) + &
                    & gammadensity2(epidat3update(amcmc1(m), 3), deltanr1op(j+1, 1), deltanr2op(j+1, 1))
                    dencurrent = dencurrent + &
                    & gammadensity2(epidat3current(amcmc1(m), 5), deltain1op(j+1, 1), deltain2op(j+1, 1)) + &
                    & gammadensity2(epidat3current(amcmc1(m), 3), deltanr1op(j+1, 1), deltanr2op(j+1, 1))
                    valuefupdate1  = valuefupdate1 + &
                    & gammadensity2(epidat3update(amcmc1(m), 5), infperiodproposalin1(1), infperiodproposalin1(2)) + &
                    & gammadensity2(epidat3update(amcmc1(m), 3), infperiodproposalnr1(1), infperiodproposalnr1(2))
                    valuefcurrent1 = valuefcurrent1+ &
                    & gammadensity2(epidat3current(amcmc1(m), 5), infperiodproposalin1(1), infperiodproposalin1(2)) + &
                    & gammadensity2(epidat3current(amcmc1(m), 3), infperiodproposalnr1(1), infperiodproposalnr1(2))
                else
                    denupdate  = denupdate  + &
                    & gammadensity2(epidat3update(amcmc1(m), 3), deltanr1op(j+1, 1), deltanr2op(j+1, 1))
                    dencurrent = dencurrent + &
                    & gammadensity2(epidat3current(amcmc1(m), 3), deltanr1op(j+1, 1), deltanr2op(j+1, 1))
                    valuefupdate1  = valuefupdate1 + &
                    & gammadensity2(epidat3update(amcmc1(m), 3), infperiodproposalnr1(1), infperiodproposalnr1(2))
                    valuefcurrent1 = valuefcurrent1+ &
                    & gammadensity2(epidat3current(amcmc1(m), 3), infperiodproposalnr1(1), infperiodproposalnr1(2))
                end if
            end do

            valueffupdate  = valueupdate1  + denupdate
            valueffcurrent = valuecurrent1 + dencurrent

            !calculating the log acceptance ratio:
            ratio1    = (valueffupdate+valuefcurrent1)-(valueffcurrent+valuefupdate1)

            psi1      = min(1.0_c_double, dexp(ratio1))

            ! decision: accept/reject the proposed value:
            call randomnumber(u1)

            if ( (u1 .le. psi1) .and. (any(epidat3update(amcmc1, 6) .lt.&
            & epidat3current(obs_inf_time, 6)).eqv. .false.) ) then
                epidat3current(amcmc1, 2) = epidat3update(amcmc1, 2)
                epidat3current(amcmc1, 3) = epidat3update(amcmc1, 3)
                epidat3current(amcmc1, 5) = epidat3update(amcmc1, 5)
                epidat3current(amcmc1, 6) = epidat3update(amcmc1, 6)
                valuecurrent1 = valueupdate1
            else
                epidat3update(amcmc1, 2)  = epidat3current(amcmc1, 2)
                epidat3update(amcmc1, 3)  = epidat3current(amcmc1, 3)
                epidat3update(amcmc1, 5)  = epidat3current(amcmc1, 5)
                epidat3update(amcmc1, 6)  = epidat3current(amcmc1, 6)
                valuecurrent1 = valuecurrent1
            end if

        end if

        epidat3 = epidat3current

        deallocate(amcmc)
        deallocate(amcmc1)

        epidatmcrem(j+1, :)     = epidat3(:, 2)
        epidatmcpernr(j+1, :)   = epidat3(:, 3)
        epidatmcperin(j+1, :)   = epidat3(:, 5)
        epidatmctim(j+1, :)     = epidat3(:, 6)

    else if (anum2(6) .eq. 3) then
    ! if the infection times/incubation periods , and removal times/delay periods are assumed known:

        deltain1op(j+1, 1) = deltain1op(j, 1)
        deltain2op(j+1, 1) = deltain2op(j, 1)

        deltanr1op(j+1, 1) = deltanr1op(j, 1)
        deltanr2op(j+1, 1) = deltanr2op(j, 1)

        epidatmcpernr(j+1, 1:n) = (/epidat3(1:n, 3)/)
        epidatmcrem(j+1, 1:n) = (/epidat3(1:n, 2)/)
        epidatmcperin(j+1, 1:n) = (/epidat3(1:n, 5)/)
        epidatmctim(j+1, 1:n) = (/epidat3(1:n, 6)/)

    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: susceptibility parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (anum2(1) .eq. 1)then

        postparsusupdate(1:nsuspar)    = susparop(j, 1:nsuspar)
        postparsuscurrent(1:nsuspar)   = susparop(j, 1:nsuspar)

        do r = 1 , nsuspar

            ! proposing candidate and making sure is positive:
            zparsus = 0.0_c_double
            abs32 = 0
            do while(abs32 .eq. 0)
                zparsus = randnormal2(postparsuscurrent(r), susproposalvar(r))
                if (zparsus .gt. 0.0_c_double)then
                    abs32 = 1
                else
                    abs32 = 0
                end if
            end do

            postparsusupdate(r)= zparsus

            !calculating the log-likelihood function:

            call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
            & suscov, transcov, postparsusupdate, powersusparop(j, :), &
            & transparop(j, :), powertransparop(j, :), kernelparop(j, :), &
            & sparkop(j, 1), gammaop(j, 1), valueupdate2)

            !calculating the log-prior:

            if (priordistsuspar(r) .eq. 1)then
                valuefupdate2    = gammadensity2(postparsusupdate(r), &
                & priorpar1sus(r), priorpar2sus(r))
                valuefcurrent2   = gammadensity2(postparsuscurrent(r), &
                & priorpar1sus(r), priorpar2sus(r))
            else if (priordistsuspar(r) .eq. 2)then
                valuefupdate2    = halfnormalden2(postparsusupdate(r), &
                & priorpar1sus(r), priorpar2sus(r))
                valuefcurrent2   = halfnormalden2(postparsuscurrent(r), &
                & priorpar1sus(r), priorpar2sus(r))
            else if (priordistsuspar(r) .eq. 3)then
                valuefupdate2    = uniformden2(postparsusupdate(r), &
                & priorpar1sus(r), priorpar2sus(r))
                valuefcurrent2   = uniformden2(postparsuscurrent(r), &
                & priorpar1sus(r), priorpar2sus(r))
            end if

            !calculating the log acceptance ratio:
            ratio2 = (valueupdate2 + valuefupdate2) - (valuecurrent1 + valuefcurrent2)

            psi2 = min(1.0_c_double, dexp(ratio2))

            ! decision: accept/reject the proposed value:

            call randomnumber(u2)

            if (psi2 .ge. u2) then
                susparop(j+1, r)        = postparsusupdate(r)
                postparsuscurrent(r)   = postparsusupdate(r)
                valuecurrent1 = valueupdate2
            else
                susparop(j+1, r)        = susparop(j, r)
                postparsusupdate(r)    = postparsuscurrent(r)
                valuecurrent1 = valuecurrent1
            end if

        end do   !  end loop r nparsus

    ! if it is not included in the model:
    else
        susparop(j+1, 1:nsuspar)  = susparop(j, 1:nsuspar)
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: susceptibility power parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (anum2(7) .eq. 1)then


        postparsusupdate(1:nsuspar)    = powersusparop(j, 1:nsuspar)
        postparsuscurrent(1:nsuspar)   = powersusparop(j, 1:nsuspar)

        do r = 1 , nsuspar

        ! proposing candidate and making sure is positive:
        zparsus = 0.0_c_double
        abs32 = 0
        do while(abs32 .eq. 0)
            zparsus = randnormal2(postparsuscurrent(r), powersusproposalvar(r))
            if (zparsus .gt. 0.0_c_double)then
                abs32 = 1
            else
                abs32 = 0
            end if
        end do

        postparsusupdate(r)= zparsus

        !calculating the log-likelihood function:

        call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
        & suscov, transcov, susparop(j+1, :), postparsusupdate, &
        & transparop(j, :), powertransparop(j, :), kernelparop(j, :), &
        & sparkop(j, 1), gammaop(j, 1), valueupdate2)

        !calculating the log-prior:

        if (priordistpowersus(r) .eq. 1)then
            valuefupdate2    = gammadensity2(postparsusupdate(r), &
            & priorpar1powersus(r), priorpar2powersus(r))
            valuefcurrent2   = gammadensity2(postparsuscurrent(r), &
            & priorpar1powersus(r), priorpar2powersus(r))
        else if (priordistpowersus(r) .eq. 2)then
            valuefupdate2    = halfnormalden2(postparsusupdate(r), &
            & priorpar1powersus(r), priorpar2powersus(r))
            valuefcurrent2   = halfnormalden2(postparsuscurrent(r), &
            & priorpar1powersus(r), priorpar2powersus(r))
        else if (priordistpowersus(r) .eq. 3)then
            valuefupdate2    = uniformden2(postparsusupdate(r), &
            & priorpar1powersus(r), priorpar2powersus(r))
            valuefcurrent2   = uniformden2(postparsuscurrent(r), &
            & priorpar1powersus(r), priorpar2powersus(r))
        end if

        !calculating the log acceptance ratio:
        ratio2 = (valueupdate2 + valuefupdate2) - (valuecurrent1 + valuefcurrent2)

        psi2 = min(1.0_c_double, dexp(ratio2))

        ! judge accept/reject the proposed value:
        call randomnumber(u2)

        if (psi2 .ge. u2) then
            powersusparop(j+1, r)   = postparsusupdate(r)
            postparsuscurrent(r)   = postparsusupdate(r)
                valuecurrent1 = valueupdate2
        else
            powersusparop(j+1, r)   = powersusparop(j, r)
            postparsusupdate(r)    = postparsuscurrent(r)
                valuecurrent1 = valuecurrent1
        end if

        end do   !  end loop r nparsus

    else
        powersusparop(j+1, 1:nsuspar)  = powersusparop(j, 1:nsuspar)
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: transmissibility parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (anum2(2) .eq. 1)then

        postpartransupdate(1:ntranspar)    = transparop(j, 1:ntranspar)
        postpartranscurrent(1:ntranspar)   = transparop(j, 1:ntranspar)

        do r = 1 , ntranspar

            ! proposing candidate and making sure is positive:
            zpartrans = 0.0_c_double

            abs33 = 0

            do while(abs33 .eq. 0)
                zpartrans = randnormal2(postpartranscurrent(r), transproposalvar(r))
                if (zpartrans .gt. 0.0_c_double)then
                    abs33 = 1
                else
                    abs33 = 0
                end if
            end do

            postpartransupdate(r)= zpartrans

            !calculating the log-likelihood function:

            call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
            & suscov, transcov, susparop(j+1, :), powersusparop(j+1, :), &
            & postpartransupdate, powertransparop(j, :), kernelparop(j, :), &
            & sparkop(j, 1), gammaop(j, 1), valueupdate3)

            !calculating the log-prior:

            if (priordisttranspar(r) .eq. 1)then
                valuefupdate3    = gammadensity2(postpartransupdate(r), &
                & priorpar1trans(r), priorpar2trans(r))
                valuefcurrent3   = gammadensity2(postpartranscurrent(r), &
                & priorpar1trans(r), priorpar2trans(r))
            else if (priordisttranspar(r) .eq. 2)then
                valuefupdate3    = halfnormalden2(postpartransupdate(r), &
                & priorpar1trans(r), priorpar2trans(r))
                valuefcurrent3   = halfnormalden2(postpartranscurrent(r), &
                & priorpar1trans(r), priorpar2trans(r))
            else
                valuefupdate3    = uniformden2(postpartransupdate(r), &
                & priorpar1trans(r), priorpar2trans(r))
                valuefcurrent3   = uniformden2(postpartranscurrent(r), &
                & priorpar1trans(r), priorpar2trans(r))
            end if

            !calculating the log acceptance ratio:
            ratio3 = (valueupdate3 + valuefupdate3) - (valuecurrent1 + valuefcurrent3)

            psi3 = min(1.0_c_double, dexp(ratio3))

            ! judge accept/reject the proposed value:
            call randomnumber(u3)

            if (psi3 .ge. u3) then
                transparop(j+1, r)        = postpartransupdate(r)
                postpartranscurrent(r) = postpartransupdate(r)
                valuecurrent1 = valueupdate3
            else
                transparop(j+1, r)        = transparop(j, r)
                postpartransupdate(r)  = postpartranscurrent(r)
                valuecurrent1 = valuecurrent1
            end if

        end do   !  end loop r nparsus


    else
        transparop(j+1, 1:ntranspar)        = transparop(j, 1:ntranspar)
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: transmissibility power parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (anum2(8) .eq. 1)then

        postpartransupdate(1:ntranspar)    = powertransparop(j, 1:ntranspar)
        postpartranscurrent(1:ntranspar)   = powertransparop(j, 1:ntranspar)

        do r = 1 , ntranspar

            ! proposing candidate and making sure is positive:
            zpartrans = 0.0_c_double

            abs33 = 0

            do while(abs33 .eq. 0)
                zpartrans = randnormal2(postpartranscurrent(r), powertransproposalvar(r))
                if (zpartrans .gt. 0.0_c_double)then
                    abs33 = 1
                else
                    abs33 = 0
                end if
            end do

            postpartransupdate(r)= zpartrans

            !calculating the log-likelihood function:

            call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
            & suscov, transcov, susparop(j+1, :), powersusparop(j+1, :), &
            & transparop(j+1, :), postpartransupdate, kernelparop(j, :), &
            & sparkop(j, 1), gammaop(j, 1), valueupdate3)

            !calculating the log-prior:

            if (priordistpowertrans(r) .eq. 1)then
                valuefupdate3    = gammadensity2(postpartransupdate(r), &
                & priorpar1powertrans(r), priorpar2powertrans(r))
                valuefcurrent3   = gammadensity2(postpartranscurrent(r), &
                & priorpar1powertrans(r), priorpar2powertrans(r))
            else if (priordistpowertrans(r) .eq. 2)then
                valuefupdate3    = halfnormalden2(postpartransupdate(r), &
                & priorpar1powertrans(r), priorpar2powertrans(r))
                valuefcurrent3   = halfnormalden2(postpartranscurrent(r), &
                & priorpar1powertrans(r), priorpar2powertrans(r))
            else if (priordistpowertrans(r) .eq. 3)then
                valuefupdate3    = uniformden2(postpartransupdate(r), &
                & priorpar1powertrans(r), priorpar2powertrans(r))
                valuefcurrent3   = uniformden2(postpartranscurrent(r), &
                & priorpar1powertrans(r), priorpar2powertrans(r))
            end if

            !calculating the log acceptance ratio:
            ratio3 = (valueupdate3 + valuefupdate3) - (valuecurrent1 + valuefcurrent3)

            psi3 = min(1.0_c_double, dexp(ratio3))

            ! judge accept/reject the proposed value:
            call randomnumber(u3)

            if (psi3 .ge. u3) then
                powertransparop(j+1, r) = postpartransupdate(r)
                postpartranscurrent(r) = postpartransupdate(r)
                valuecurrent1 = valueupdate3
            else
                powertransparop(j+1, r) = powertransparop(j, r)
                postpartransupdate(r)  = postpartranscurrent(r)
                valuecurrent1 = valuecurrent1
            end if

        end do   !  end loop r nparsus


    else
        powertransparop(j+1, 1:ntranspar) = powertransparop(j, 1:ntranspar)
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: spark parameter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (anum2(3) .eq. 1)then

        ! proposing candidate and making sure is positive:
        abs34 = 0

        do while(abs34 .eq. 0)
            zspark = randnormal2(sparkop(j, 1), sparkproposalvar)
            if (zspark .gt. 0.0_c_double)then
                abs34 = 1
            else
                abs34 = 0
            end if
        end do        ! end while abs3

        postsparkupdate = zspark

        !calculating the log-likelihood function:

        call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
        & suscov, transcov, susparop(j+1, :), powersusparop(j+1, :), &
        & transparop(j+1, :), powertransparop(j+1, :), kernelparop(j, :), &
        & postsparkupdate, gammaop(j, 1), valueupdate4)

        !calculating the log-prior:

        if (priordistsparkpar .eq. 1)then
            valuefupdate4    = gammadensity2(postsparkupdate, &
            & sparkprior(1), sparkprior(2))
            valuefcurrent4   = gammadensity2(sparkop(j, 1), &
            & sparkprior(1), sparkprior(2))
        else if (priordistsparkpar .eq. 2)then
            valuefupdate4    = halfnormalden2(postsparkupdate, &
            & sparkprior(1), sparkprior(2))
            valuefcurrent4   = halfnormalden2(sparkop(j, 1), &
            & sparkprior(1), sparkprior(2))
        else
            valuefupdate4    = uniformden2(postsparkupdate, &
            & sparkprior(1), sparkprior(2))
            valuefcurrent4   = uniformden2(sparkop(j, 1), &
            & sparkprior(1), sparkprior(2))
        end if

        !calculating the log acceptance ratio:
        ratio4 = (valueupdate4 + valuefupdate4)-(valuecurrent1 + valuefcurrent4)

        psi4 = min(1.0_c_double, dexp(ratio4))

        ! judge accept/reject the proposed value:
        call randomnumber(u4)

        if (psi4 .ge. u4) then
            sparkop(j+1, 1) = postsparkupdate
                valuecurrent1 = valueupdate4
        else
            sparkop(j+1, 1) = sparkop(j, 1)
                valuecurrent1 = valuecurrent1
        end if

    else
        sparkop(j+1, 1) = sparkop(j, 1)
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: kernelpar parameters (spatial and network effect parameters)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (anum2(4) .eq. 1)then
    ! for distance-based continuous-time ILMS:
    ! just the spatial parameter is needed to be updated

    ! proposing candidate and making sure is positive:
        abs35 = 0

        do while(abs35 .eq. 0)
            zkernelpar = randnormal2(kernelparop(j, 1), kernelparproposalvar(1))
            if (zkernelpar .gt. 0.0_c_double)then
                abs35 = 1
            else
                abs35 = 0
            end if
        end do        ! end while abs3

        postkernelparupdate = (/ zkernelpar, kernelparop(j, 2) /)

        ! calculate the log-likelihood function:
        call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
        & suscov, transcov, susparop(j+1, :), powersusparop(j+1, :), &
        & transparop(j+1, :), powertransparop(j+1, :), postkernelparupdate, &
        & sparkop(j+1, 1), gammaop(j, 1), valueupdate5)

        !calculating the log-prior:
        if (priordistkernelparpar(1) .eq. 1)then
            valuefupdate5    = gammadensity2(postkernelparupdate(1), &
            & kernelparprior(1, 1), kernelparprior(1, 2))
            valuefcurrent5   = gammadensity2(kernelparop(j, 1), &
            & kernelparprior(1, 1), kernelparprior(1, 2))
        else if (priordistkernelparpar(1) .eq. 2)then
            valuefupdate5    = halfnormalden2(postkernelparupdate(1), &
            & kernelparprior(1, 1), kernelparprior(1, 2))
            valuefcurrent5   = halfnormalden2(kernelparop(j, 1), &
            & kernelparprior(1, 1), kernelparprior(1, 2))
        else
            valuefupdate5    = uniformden2(postkernelparupdate(1), &
            & kernelparprior(1, 1), kernelparprior(1, 2))
            valuefcurrent5   = uniformden2(kernelparop(j, 1), &
            & kernelparprior(1, 1), kernelparprior(1, 2))
        end if


        !calculating the log acceptance ratio:
        ratio5 = (valueupdate5 + valuefupdate5)-(valuecurrent1 + valuefcurrent5)

        psi5 = min(1.0_c_double, dexp(ratio5))

        ! judge accept/reject the proposed value:
        call randomnumber(u5)

        if (psi5 .ge. u5) then
            kernelparop(j+1, 1) = postkernelparupdate(1)
            kernelparop(j+1, 2) = kernelparop(j, 2)
                valuecurrent1 = valueupdate5
        else
            kernelparop(j+1, 1) = kernelparop(j, 1)
            kernelparop(j+1, 2) = kernelparop(j, 2)
                valuecurrent1 = valuecurrent1
        end if

    ! for combined network- and distance-based continuous-time ILMS:
    else if (anum2(4) .eq. 2)then

    ! updating the spatial parameter:
    ! proposing candidate and making sure is positive:
        abs35 = 0
        do while(abs35 .eq. 0)
            zkernelpar = randnormal2(kernelparop(j, 1), kernelparproposalvar(1))
            if (zkernelpar .gt. 0.0_c_double)then
                abs35 = 1
            else
                abs35 = 0
            end if
        end do

        postkernelparupdate = (/ zkernelpar, kernelparop(j, 2) /)
        postkernelparcurrent = kernelparop(j, :)

        ! calculate the log-likelihood function:
        call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
        & suscov, transcov, susparop(j+1, :), powersusparop(j+1, :), &
        & transparop(j+1, :), powertransparop(j+1, :), postkernelparupdate, &
        & sparkop(j+1, 1), gammaop(j, 1), valueupdate5)

        !calculating the log-prior:
        if (priordistkernelparpar(1) .eq. 1)then
            valuefupdate5    = gammadensity2(postkernelparupdate(1), &
            & kernelparprior(1, 1), kernelparprior(1, 2))
            valuefcurrent5   = gammadensity2(postkernelparcurrent(1), &
            & kernelparprior(1, 1), kernelparprior(1, 2))
        else if (priordistkernelparpar(1) .eq. 2)then
            valuefupdate5    = halfnormalden2(postkernelparupdate(1), &
            & kernelparprior(1, 1), kernelparprior(1, 2))
            valuefcurrent5   = halfnormalden2(postkernelparcurrent(1), &
            & kernelparprior(1, 1), kernelparprior(1, 2))
        else
            valuefupdate5    = uniformden2(postkernelparupdate(1), &
            & kernelparprior(1, 1), kernelparprior(1, 2))
            valuefcurrent5   = uniformden2(postkernelparcurrent(1), &
            & kernelparprior(1, 1), kernelparprior(1, 2))
        end if


        !calculating the log acceptance ratio:
        ratio5 = (valueupdate5 + valuefupdate5)-(valuecurrent1 + valuefcurrent5)

        psi5 = min(1.0_c_double, dexp(ratio5))

        ! judge accept/reject the proposed value:
        call randomnumber(u5)

        if (psi5 .ge. u5) then
            kernelparop(j+1, 1) = postkernelparupdate(1)
                valuecurrent1 = valueupdate5
        else
            kernelparop(j+1, 1) = kernelparop(j, 1)
                valuecurrent1 = valuecurrent1
        end if


        ! updating the network-effect parameter:
        ! proposing candidate and making sure is positive:
        abs35 = 0
        do while(abs35 .eq. 0)
            zkernelpar = randnormal2(kernelparop(j, 2), kernelparproposalvar(2))
            if (zkernelpar .gt. 0.0_c_double)then
                abs35 = 1
            else
                abs35 = 0
            end if
        end do        ! end while abs3

        postkernelparupdate = (/ kernelparop(j+1, 1), zkernelpar /)
        postkernelparcurrent = (/ kernelparop(j+1, 1), kernelparop(j, 2) /)

        ! calculate the log-likelihood function:
        call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
        & suscov, transcov, susparop(j+1, :), powersusparop(j+1, :), &
        & transparop(j+1, :), powertransparop(j+1, :), postkernelparupdate, &
        & sparkop(j+1, 1), gammaop(j, 1), valueupdate5)

        !calculating the log-prior:
        if (priordistkernelparpar(2) .eq. 1)then
            valuefupdate5    = gammadensity2(postkernelparupdate(2), &
            & kernelparprior(2, 1), kernelparprior(2, 2))
            valuefcurrent5   = gammadensity2(postkernelparcurrent(2), &
            & kernelparprior(2, 1), kernelparprior(2, 2))
        else if (priordistkernelparpar(2) .eq. 2)then
            valuefupdate5    = halfnormalden2(postkernelparupdate(2), &
            & kernelparprior(2, 1), kernelparprior(2, 2))
            valuefcurrent5   = halfnormalden2(postkernelparcurrent(2), &
            & kernelparprior(2, 1), kernelparprior(2, 2))
        else
            valuefupdate5    = uniformden2(postkernelparupdate(2), &
            & kernelparprior(2, 1), kernelparprior(2, 2))
            valuefcurrent5   = uniformden2(postkernelparcurrent(2), &
            & kernelparprior(2, 1), kernelparprior(2, 2))
        end if


        !calculating the log acceptance ratio:
        ratio5 = (valueupdate5 + valuefupdate5)-(valuecurrent1 + valuefcurrent5)

        psi5 = min(1.0_c_double, dexp(ratio5))

        ! judge accept/reject the proposed value:

        call randomnumber(u5)

        if (psi5 .ge. u5) then
            kernelparop(j+1, 2) = postkernelparupdate(2)
                valuecurrent1 = valueupdate5
        else
            kernelparop(j+1, 2) = kernelparop(j, 2)
                valuecurrent1 = valuecurrent1
        end if

    ! for network-based continuous-time ILMS:
    else
        kernelparop(j+1, :) = kernelparop(j, :)
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: gamma parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (anum2(5) .eq. 1)then

    ! proposing candidate and making sure is positive:

        abs35 = 0
        do while(abs35 .eq. 0)
            zgamma = randnormal2(gammaop(j, 1), gammaproposalvar)
            if (zgamma .gt. 0.0_c_double)then
                abs35 = 1
            else
                abs35 = 0
            end if
        end do        ! end while abs3

        postgammaupdate = zgamma

        !calculating the log-likelihood function:

        call likelihoodsinrsinr2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
        & suscov, transcov, susparop(j+1, :), powersusparop(j+1, :), &
        & transparop(j+1, :), powertransparop(j+1, :), kernelparop(j+1, :), &
        & sparkop(j+1, 1), postgammaupdate, valueupdate6)

        !calculating the log-prior:

        if (priordistgammapar .eq. 1)then
            valuefupdate6    = gammadensity2(postgammaupdate, &
            & gammaprior(1), gammaprior(2))
            valuefcurrent6   = gammadensity2(gammaop(j, 1), &
            & gammaprior(1), gammaprior(2))
        else if (priordistgammapar .eq. 2)then
            valuefupdate6    = halfnormalden2(postgammaupdate, &
            & gammaprior(1), gammaprior(2))
            valuefcurrent6   = halfnormalden2(gammaop(j, 1), &
            & gammaprior(1), gammaprior(2))
        else if (priordistgammapar .eq. 3)then
            valuefupdate6    = uniformden2(postgammaupdate, &
            & gammaprior(1), gammaprior(2))
            valuefcurrent6   = uniformden2(gammaop(j, 1), &
            & gammaprior(1), gammaprior(2))
        end if


        !calculating the log acceptance ratio:
        ratio6 = (valueupdate6 + valuefupdate6)-(valuecurrent1 + valuefcurrent6)

        psi6 = min(1.0_c_double, dexp(ratio6))

        ! judge accept/reject the proposed value:
        call randomnumber(u6)

        if (psi6 .ge. u6) then
            gammaop(j+1, 1) = postgammaupdate
                valuecurrent1 = valueupdate6
        else
            gammaop(j+1, 1) = gammaop(j, 1)
                valuecurrent1 = valuecurrent1
        end if

    else
        gammaop(j+1, 1) = gammaop(j, 1)
    end if

! to update the log-likelihood based on the updated parameters:
    if (anum2(6) .eq. 1) then

        incdens = 0.0_c_double
        do i = 1, ni
            if (i .gt. obs_inf_time) then
                incdens = incdens + gammadensity2(epidat3(i, 5), deltain1op(j+1, 1), deltain2op(j+1, 1))
            else
                incdens = incdens
            end if
        end do

        loglik(j+1, 1) = valuecurrent1 + incdens

    else if (anum2(6) .eq. 2) then

        incdens = 0.0_c_double
        delaydens = 0.0_c_double
        do i = 1, ni
            if (i .gt. obs_inf_time) then
                incdens = incdens + gammadensity2(epidat3(i, 5), deltain1op(j+1, 1), deltain2op(j+1, 1))
            else
                incdens = incdens
            end if
            delaydens = delaydens + gammadensity2(epidat3(i, 3), deltanr1op(j+1, 1), deltanr2op(j+1, 1))
        end do

        loglik(j+1, 1) = valuecurrent1 + incdens + delaydens

    else if (anum2(6) .eq. 3) then

        loglik(j+1, 1) = valuecurrent1

    end if

!    call iter1(nsim, 10)

end do !end loop j (mcmc)

call seedout()

end subroutine mcmcsinr_f



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                 Log LIKELIHOOD subroutine                 !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine likelihoodsinrsinr2(n, ninfected, num, nsuspar, ntranspar, cc, d3, epidat, &
    & suscov, transcov, suspar, powersus, transpar, powertrans, kernelpar, spark, gamma, likk)

    external infinity_value

    integer:: j, i, r, m
    integer, intent(in)::n, ninfected, num, nsuspar, ntranspar       !integers
    real (C_DOUBLE), intent(in), dimension(2):: kernelpar        ! parameter of the kernel function
    real (C_DOUBLE), intent(in):: spark, gamma                   ! spark & notification effect parameters
    real (C_DOUBLE), intent(in), dimension(n, n)::d3, cc           ! network & distance matrices
    real (C_DOUBLE), intent(in), dimension(n, 6)::epidat          ! epidemic data
    real (C_DOUBLE), dimension(n, 6)::epidat1                     ! copy of sorted epidemic data (regarding to the infection times)
    real (C_DOUBLE), intent(in), dimension(n, nsuspar)::suscov    ! susceptibility covariates
    real (C_DOUBLE), intent(in), dimension(n, ntranspar)::transcov    ! transmissibility covariates
    real (C_DOUBLE), intent(in), dimension(nsuspar)::suspar, powersus ! susceptibility parameters
    real (C_DOUBLE), intent(in), dimension(ntranspar)::transpar, powertrans ! transmissibility parameters
!    real (C_DOUBLE), dimension(nsuspar)::suscov1                ! suseptible covar. for one individual
!    real (C_DOUBLE), dimension(ntranspar)::transcov1            ! transmissibilty covar. for one individual
    real (C_DOUBLE), intent(out)::likk                          ! OUTPUT
    real (C_DOUBLE)::rate, tt, ss, qa1, qa2, likk1, likk2, Inf
    real (C_DOUBLE), dimension(n) :: gh, sss
    real (C_DOUBLE), dimension(ninfected) :: rt, df1, rrate

    call infinity_value(Inf)

    SELECT CASE (num)

    CASE (1)

    ! network-based continuous ILMs

    ! getting t_obs (maximum removal time):

    epidat1 = epidat

    call Sort1(epidat1(:,6), epidat1, n, 6)

    tt = maxval(epidat1(1:ninfected, 2))

    !calculate the exponent part of the spark term:
    do r = 1 , n
        sss(r) = spark * (min(tt, epidat1(r, 6)) - epidat1(1, 6))
    end do

    ss = sum(sss)

    ! calculate the terms of the exponent part of the double summation:
    df1 = 0.0_c_double

    do i =1, ninfected

        qa2 = 0.0_c_double
        do m = 1 , ntranspar
            qa2 = qa2 + (transpar(m) * transcov(int(epidat1(i, 1)), m)**powertrans(m))
        end do

        gh = 0.0_c_double

        do j = i+1, n


            qa1 = 0.0_c_double
            do m = 1 , nsuspar
                qa1 = qa1 + (suspar(m) * suscov(int(epidat1(j, 1)), m)**powersus(m))
            end do

            rate = qa1  * qa2  * (cc(int(epidat1(i, 1)), int(epidat1(j, 1))))

            if (j .le. ninfected) then
                gh(j)   = (( minval( (/ tt, epidat1(i, 4), epidat1(j, 6) /) ) - min(epidat1(i, 6), epidat1(j, 6)) ) * (rate)) + &
                & (( ( minval( (/ tt, epidat1(i, 2), epidat1(j, 6) /) ) - min(epidat1(i, 6), epidat1(j, 6)) ) - &
                & ( minval( (/ tt, epidat1(i, 4), epidat1(j, 6) /) ) - min(epidat1(i, 6), epidat1(j, 6)) ) ) * (gamma * rate))
            else
                gh(j)   = (( epidat1(i, 4) - epidat1(i, 6)) * (rate)) + &
                & (( (epidat1(i, 2) - epidat1(i, 6)) - (epidat1(i, 4) - epidat1(i, 6)) ) * (gamma*rate))
            end if

        end do

        df1(i) = sum(gh)
    end do

    likk1 = (-(sum(df1)+ss))

    ! calculate the first part of the likelihood:

    rt(1) = 1.0_c_double

    do j = 2 , ninfected

            qa1 = 0.0_c_double
            do m = 1 , nsuspar
                qa1 = qa1 + (suspar(m) * suscov(int(epidat1(j, 1)), m)**powersus(m))
            end do

        rrate = 0.0_c_double

        do i = 1 , (j-1)

            qa2 = 0.0_c_double
            do m = 1 , ntranspar
                qa2 = qa2 + (transpar(m) * transcov(int(epidat1(i, 1)), m)**powertrans(m))
            end do

            if ((epidat1(j, 6) .gt. epidat1(i, 6)) .and.  (epidat1(j, 6) .le. epidat1(i, 4)) ) then

                rrate(i) = (qa1 * qa2 * (cc(int(epidat1(i, 1)), int(epidat1(j, 1)))))

            else if ((epidat1(j, 6) .gt. epidat1(i, 4)) .and.  (epidat1(j, 6) .le. epidat1(i, 2)) ) then

                rrate(i) = (gamma * (qa1 * qa2 * (cc(int(epidat1(i, 1)), int(epidat1(j, 1))))))

            else

                rrate(i)   = 0.0_c_double

            end if

        end do


        rt(j) = sum(rrate) + spark

    end do

    likk2 = sum(log(rt))


    ! TOTAL LOG-LIKELIHOOD without the density of the infectious periods:

    likk =  likk1 + likk2

    !##################################################################
    !##################################################################

    CASE (2)

    ! distance-based continuous ILMs (powerlaw kernel)

    ! getting t_obs (maximum removal time):

    epidat1 = epidat

    call Sort1(epidat1(:,6), epidat1, n, 6)

    tt = maxval(epidat1(1:ninfected, 2))

    !calculate the exponent part of the spark term:
    do r = 1 , n
        sss(r) = spark * (min(tt, epidat1(r, 6)) - epidat1(1, 6))
    end do

    ss = sum(sss)

    ! calculate the terms of the exponent part of the double summation:
    df1 = 0.0_c_double

    do i =1, ninfected

        qa2 = 0.0_c_double
        do m = 1 , ntranspar
            qa2 = qa2 + (transpar(m) * transcov(int(epidat1(i, 1)), m)**powertrans(m))
        end do

        gh = 0.0_c_double

        do j = i+1, n


            qa1 = 0.0_c_double
            do m = 1 , nsuspar
                qa1 = qa1 + (suspar(m) * suscov(int(epidat1(j, 1)), m)**powersus(m))
            end do

            rate = qa1  * qa2  * (d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(-kernelpar(1)))

            if (j .le. ninfected) then
                gh(j)   = (( minval( (/ tt, epidat1(i, 4), epidat1(j, 6) /) ) - min(epidat1(i, 6), epidat1(j, 6)) ) * (rate)) + &
                & (( ( minval( (/ tt, epidat1(i, 2), epidat1(j, 6) /) ) - min(epidat1(i, 6), epidat1(j, 6)) ) - &
                & ( minval( (/ tt, epidat1(i, 4), epidat1(j, 6) /) ) - min(epidat1(i, 6), epidat1(j, 6)) ) ) * (gamma * rate))
            else
                gh(j)   = (( epidat1(i, 4) - epidat1(i, 6)) * (rate)) + &
                & (( (epidat1(i, 2) - epidat1(i, 6)) - (epidat1(i, 4) - epidat1(i, 6)) ) * (gamma*rate))
            end if

        end do

        df1(i) = sum(gh)
    end do

    likk1 = (-(sum(df1)+ss))

    ! calculate the first part of the likelihood:

    rt(1) = 1.0_c_double

    do j = 2 , ninfected

            qa1 = 0.0_c_double
            do m = 1 , nsuspar
                qa1 = qa1 + (suspar(m) * suscov(int(epidat1(j, 1)), m)**powersus(m))
            end do

        rrate = 0.0_c_double

        do i = 1 , (j-1)

            qa2 = 0.0_c_double
            do m = 1 , ntranspar
                qa2 = qa2 + (transpar(m) * transcov(int(epidat1(i, 1)), m)**powertrans(m))
            end do

            if ((epidat1(j, 6) .gt. epidat1(i, 6)) .and.  (epidat1(j, 6) .le. epidat1(i, 4)) ) then

                rrate(i) = (qa1 * qa2 * (d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(-kernelpar(1))))

            else if ((epidat1(j, 6) .gt. epidat1(i, 4)) .and.  (epidat1(j, 6) .le. epidat1(i, 2)) ) then

                rrate(i) = (gamma * (qa1 * qa2 * (d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(-kernelpar(1)))))

            else

                rrate(i)   = 0.0_c_double

            end if

        end do


        rt(j) = sum(rrate) + spark

    end do

    likk2 = sum(log(rt))


    ! TOTAL LOG-LIKELIHOOD without the density of the infectious periods:

    likk =  likk1 + likk2

    !##################################################################
    !##################################################################


    CASE (3)

    ! distance-based continuous ILMs (Cauchy kernel)

    ! getting t_obs (maximum removal time):

    epidat1 = epidat

    call Sort1(epidat1(:,6), epidat1, n, 6)

    tt = maxval(epidat1(1:ninfected, 2))

    !calculate the exponent part of the spark term:
    do r = 1 , n
        sss(r) = spark * (min(tt, epidat1(r, 6)) - epidat1(1, 6))
    end do

    ss = sum(sss)

    ! calculate the terms of the exponent part of the double summation:
    df1 = 0.0_c_double

    do i =1, ninfected

        qa2 = 0.0_c_double
        do m = 1 , ntranspar
            qa2 = qa2 + (transpar(m) * transcov(int(epidat1(i, 1)), m)**powertrans(m))
        end do

        gh = 0.0_c_double

        do j = i+1, n


            qa1 = 0.0_c_double
            do m = 1 , nsuspar
                qa1 = qa1 + (suspar(m) * suscov(int(epidat1(j, 1)), m)**powersus(m))
            end do

            rate = qa1  * qa2  * (kernelpar(1)/((d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(2.0_c_double))+ &
                    & (kernelpar(1)**(2.0_c_double))))

            if (j .le. ninfected) then
                gh(j)   = (( minval( (/ tt, epidat1(i, 4), epidat1(j, 6) /) ) - min(epidat1(i, 6), epidat1(j, 6)) ) * (rate)) + &
                & (( ( minval( (/ tt, epidat1(i, 2), epidat1(j, 6) /) ) - min(epidat1(i, 6), epidat1(j, 6)) ) - &
                & ( minval( (/ tt, epidat1(i, 4), epidat1(j, 6) /) ) - min(epidat1(i, 6), epidat1(j, 6)) ) ) * (gamma * rate))
            else
                gh(j)   = (( epidat1(i, 4) - epidat1(i, 6)) * (rate)) + &
                & (( (epidat1(i, 2) - epidat1(i, 6)) - (epidat1(i, 4) - epidat1(i, 6)) ) * (gamma*rate))
            end if

        end do

        df1(i) = sum(gh)
    end do

    likk1 = (-(sum(df1)+ss))

    ! calculate the first part of the likelihood:

    rt(1) = 1.0_c_double

    do j = 2 , ninfected

            qa1 = 0.0_c_double
            do m = 1 , nsuspar
                qa1 = qa1 + (suspar(m) * suscov(int(epidat1(j, 1)), m)**powersus(m))
            end do

        rrate = 0.0_c_double

        do i = 1 , (j-1)

            qa2 = 0.0_c_double
            do m = 1 , ntranspar
                qa2 = qa2 + (transpar(m) * transcov(int(epidat1(i, 1)), m)**powertrans(m))
            end do

            if ((epidat1(j, 6) .gt. epidat1(i, 6)) .and.  (epidat1(j, 6) .le. epidat1(i, 4)) ) then

                rrate(i) = (qa1 * qa2 * (kernelpar(1)/((d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(2.0_c_double))+ &
                    & (kernelpar(1)**(2.0_c_double)))))

            else if ((epidat1(j, 6) .gt. epidat1(i, 4)) .and.  (epidat1(j, 6) .le. epidat1(i, 2)) ) then

                rrate(i) = (gamma * (qa1 * qa2 * (kernelpar(1)/((d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(2.0_c_double))+ &
                    & (kernelpar(1)**(2.0_c_double))))))

            else

                rrate(i)   = 0.0_c_double

            end if

        end do


        rt(j) = sum(rrate) + spark

    end do

    likk2 = sum(log(rt))


    ! TOTAL LOG-LIKELIHOOD without the density of the infectious periods:

    likk =  likk1 + likk2

    !##################################################################
    !##################################################################

    CASE (4)

    ! distance and network-based continuous ILMs (powerlaw distance kernel)

    ! getting t_obs (maximum removal time):

    epidat1 = epidat

    call Sort1(epidat1(:,6), epidat1, n, 6)

    tt = maxval(epidat1(1:ninfected, 2))

    !calculate the exponent part of the spark term:
    do r = 1 , n
        sss(r) = spark * (min(tt, epidat1(r, 6)) - epidat1(1, 6))
    end do

    ss = sum(sss)

    ! calculate the terms of the exponent part of the double summation:
    df1 = 0.0_c_double

    do i =1, ninfected

        qa2 = 0.0_c_double
        do m = 1 , ntranspar
            qa2 = qa2 + (transpar(m) * transcov(int(epidat1(i, 1)), m)**powertrans(m))
        end do

        gh = 0.0_c_double

        do j = i+1, n


            qa1 = 0.0_c_double
            do m = 1 , nsuspar
                qa1 = qa1 + (suspar(m) * suscov(int(epidat1(j, 1)), m)**powersus(m))
            end do

            rate = qa1  * qa2  * ((d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(-kernelpar(1))) + &
                    & (kernelpar(2)*cc(int(epidat1(i, 1)), int(epidat1(j, 1)))))

            if (j .le. ninfected) then
                gh(j)   = (( minval( (/ tt, epidat1(i, 4), epidat1(j, 6) /) ) - min(epidat1(i, 6), epidat1(j, 6)) ) * (rate)) + &
                & (( ( minval( (/ tt, epidat1(i, 2), epidat1(j, 6) /) ) - min(epidat1(i, 6), epidat1(j, 6)) ) - &
                & ( minval( (/ tt, epidat1(i, 4), epidat1(j, 6) /) ) - min(epidat1(i, 6), epidat1(j, 6)) ) ) * (gamma * rate))
            else
                gh(j)   = (( epidat1(i, 4) - epidat1(i, 6)) * (rate)) + &
                & (( (epidat1(i, 2) - epidat1(i, 6)) - (epidat1(i, 4) - epidat1(i, 6)) ) * (gamma*rate))
            end if

        end do

        df1(i) = sum(gh)
    end do

    likk1 = (-(sum(df1)+ss))

    ! calculate the first part of the likelihood:

    rt(1) = 1.0_c_double

    do j = 2 , ninfected

            qa1 = 0.0_c_double
            do m = 1 , nsuspar
                qa1 = qa1 + (suspar(m) * suscov(int(epidat1(j, 1)), m)**powersus(m))
            end do

        rrate = 0.0_c_double

        do i = 1 , (j-1)

            qa2 = 0.0_c_double
            do m = 1 , ntranspar
                qa2 = qa2 + (transpar(m) * transcov(int(epidat1(i, 1)), m)**powertrans(m))
            end do

            if ((epidat1(j, 6) .gt. epidat1(i, 6)) .and.  (epidat1(j, 6) .le. epidat1(i, 4)) ) then

                rrate(i) = (qa1 * qa2 * ((d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(-kernelpar(1))) + &
                    & (kernelpar(2)*cc(int(epidat1(i, 1)), int(epidat1(j, 1))))))

            else if ((epidat1(j, 6) .gt. epidat1(i, 4)) .and.  (epidat1(j, 6) .le. epidat1(i, 2)) ) then

                rrate(i) = (gamma * qa1 * qa2 * (((d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(-kernelpar(1))) + &
                    & (kernelpar(2)*cc(int(epidat1(i, 1)), int(epidat1(j, 1)))))))

            else

                rrate(i)   = 0.0_c_double

            end if

        end do


        rt(j) = sum(rrate) + spark

    end do

    likk2 = sum(log(rt))


    ! TOTAL LOG-LIKELIHOOD without the density of the infectious periods:

    likk =  likk1 + likk2


    !##################################################################
    !##################################################################

    CASE (5)

    ! distance and network-based continuous ILMs (Cauchy distance kernel)

    ! getting t_obs (maximum removal time):

    epidat1 = epidat

    call Sort1(epidat1(:,6), epidat1, n, 6)

    tt = maxval(epidat1(1:ninfected, 2))

    !calculate the exponent part of the spark term:
    do r = 1 , n
        sss(r) = spark * (min(tt, epidat1(r, 6)) - epidat1(1, 6))
    end do

    ss = sum(sss)

    ! calculate the terms of the exponent part of the double summation:
    df1 = 0.0_c_double

    do i =1, ninfected

        qa2 = 0.0_c_double
        do m = 1 , ntranspar
            qa2 = qa2 + (transpar(m) * transcov(int(epidat1(i, 1)), m)**powertrans(m))
        end do

        gh = 0.0_c_double

        do j = i+1, n


            qa1 = 0.0_c_double
            do m = 1 , nsuspar
                qa1 = qa1 + (suspar(m) * suscov(int(epidat1(j, 1)), m)**powersus(m))
            end do

            rate = qa1  * qa2  * ((kernelpar(1)/((d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(2.0_c_double))+ &
                    & (kernelpar(1)**(2.0_c_double)))) + (kernelpar(2)*cc(int(epidat1(i, 1)), int(epidat1(j, 1)))))

            if (j .le. ninfected) then
                gh(j)   = (( minval( (/ tt, epidat1(i, 4), epidat1(j, 6) /) ) - min(epidat1(i, 6), epidat1(j, 6)) ) * (rate)) + &
                & (( ( minval( (/ tt, epidat1(i, 2), epidat1(j, 6) /) ) - min(epidat1(i, 6), epidat1(j, 6)) ) - &
                & ( minval( (/ tt, epidat1(i, 4), epidat1(j, 6) /) ) - min(epidat1(i, 6), epidat1(j, 6)) ) ) * (gamma * rate))
            else
                gh(j)   = (( epidat1(i, 4) - epidat1(i, 6)) * (rate)) + &
                & (( (epidat1(i, 2) - epidat1(i, 6)) - (epidat1(i, 4) - epidat1(i, 6)) ) * (gamma*rate))
            end if

        end do

        df1(i) = sum(gh)
    end do

    likk1 = (-(sum(df1)+ss))

    ! calculate the first part of the likelihood:

    rt(1) = 1.0_c_double

    do j = 2 , ninfected

            qa1 = 0.0_c_double
            do m = 1 , nsuspar
                qa1 = qa1 + (suspar(m) * suscov(int(epidat1(j, 1)), m)**powersus(m))
            end do

        rrate = 0.0_c_double

        do i = 1 , (j-1)

            qa2 = 0.0_c_double
            do m = 1 , ntranspar
                qa2 = qa2 + (transpar(m) * transcov(int(epidat1(i, 1)), m)**powertrans(m))
            end do

            if ((epidat1(j, 6) .gt. epidat1(i, 6)) .and.  (epidat1(j, 6) .le. epidat1(i, 4)) ) then

                rrate(i) = (qa1 * qa2 * ((kernelpar(1)/((d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(2.0_c_double))+ &
                    & (kernelpar(1)**(2.0_c_double)))) + (kernelpar(2)*cc(int(epidat1(i, 1)), int(epidat1(j, 1))))))

            else if ((epidat1(j, 6) .gt. epidat1(i, 4)) .and.  (epidat1(j, 6) .le. epidat1(i, 2)) ) then

                rrate(i) = (gamma * qa1 * qa2 * ((kernelpar(1)/((d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(2.0_c_double))+ &
                    & (kernelpar(1)**(2.0_c_double)))) + (kernelpar(2)*cc(int(epidat1(i, 1)), int(epidat1(j, 1))))))

            else

                rrate(i)   = 0.0_c_double

            end if

        end do


        rt(j) = sum(rrate) + spark

    end do

    likk2 = sum(log(rt))


    ! TOTAL LOG-LIKELIHOOD without the density of the infectious periods:

    likk =  likk1 + likk2

    end SELECT


    end subroutine likelihoodsinrsinr2



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!     Densities functions for different distributions     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !#################### HALF-NORMAL density ######################

    function halfnormalden2(alpha, a, b) RESULT(pdf)
    implicit none
    real (C_DOUBLE), parameter:: pi = 3.141592653589793_c_double
    real (C_DOUBLE) :: alpha, b, a
    real (C_DOUBLE) :: val, pdf

    if (alpha .lt. a)then
        val = 0.0_c_double
    else
        val =  sqrt (2.0_c_double / pi*b) * exp ( - 0.5_c_double * ((alpha)**2 / b) )
    end if

    pdf= log(val)

    end function halfnormalden2 ! returns log value

    !#################### GAMMA density ######################

    function gammadensity2(x, a, b) RESULT(pdf1)
    implicit none
    real (C_DOUBLE) :: x, a, b
    real (C_DOUBLE) :: dn, pdf1

    if (x .le. 0.0_c_double)then
        dn = 0.0_c_double
    else
        dn = (x**(a-1.0_c_double)) * dexp(- (x*b))
    end if
    pdf1= dlog(dn)

    end function gammadensity2 ! returns log value

    !#################### UNifORM density ######################

    function uniformden2(val, a, b) result(pdf)
    implicit none
    real (C_DOUBLE):: val, a, b, pdf1, pdf

    if (val .ge. b) then
        pdf1 = 0.0_c_double
    else if (val .le. a) then
        pdf1 = 0.0_c_double
    else
        pdf1 = (1.0_c_double/(b-a))
    end if

    pdf = dlog(pdf1)

    end function uniformden2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Generating random variables for diffierent distributions !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !####################  NORMAL distribution ######################

    function randnormal2(mean, stdev) RESULT(c)
    implicit none

    real (C_DOUBLE) :: mean, stdev, c, temp(2), r, theta
    real (C_DOUBLE), PARAMETER :: PI=3.141592653589793238462_c_double

    call randomnumber(temp(1))
    call randomnumber(temp(2))
    r=(-2.0_c_double*log(temp(1)))**0.5_c_double
    theta = 2.0_c_double*PI*temp(2)
    c= mean+stdev*r*sin(theta)

    end function randnormal2

    !#################### GAMMA distribution ######################

    RECURSIVE function randgamma2(shape, SCALE) RESULT(ans)
    real (C_DOUBLE):: SHAPE, scale, u, w, d, c, x, xsq, g, ans, v
    !
    ! ## Implementation based on "A Simple Method for Generating Gamma Variables"
    ! ## by George Marsaglia and Wai Wan Tsang.
    ! ## ACM Transactions on Mathematical Software
    ! ## Vol 26, No 3, September 2000, pages 363-372.

    if (shape >= 1.0_c_double) THEN
        d = SHAPE - (1.0_c_double/3.0_c_double)
        c = 1.0_c_double/((9.0_c_double * d)**0.5_c_double)
        do while (.true.)
            x = randnormal2(0.0_c_double, 1.0_c_double)
            v = 1.0 + c*x
            do while (v <= 0.0_c_double)
                x = randnormal2(0.0_c_double, 1.0_c_double)
                v = 1.0_c_double + c*x
            end do

            v = v*v*v
            call randomnumber(u)
            xsq = x*x
            if ((u < 1.0_c_double -.0331_c_double*xsq*xsq) .OR.  &
            (log(u) < 0.5_c_double*xsq + d*(1.0_c_double - v + log(v))) )then
                ans=scale*d*v
                return
            end if
        end do
    else
        g = randgamma2(shape+1.0_c_double, 1.0_c_double)
        call randomnumber(w)
        ans=scale*g*(w**(1.0_c_double/shape))
        return
    end if

    end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!       To generate random samples without replacement       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE ransamsinr2(x, a, n, k)
    implicit none

    integer :: j, m, l
    real (C_DOUBLE) :: u
    integer, intent(in) :: n, k
    integer, intent(in), dimension(n) :: x
    integer, intent(out), dimension(k) :: a

        m=0
        do j = 1, n
            call randomnumber(u)

            l = int(float((n-j+1)) * u) + 1

            if (l .le. (k-m)) then
                m = m + 1
                a(m) = x(j)
            else
                m = m
            end if

            if (m .ge. k) then
                exit
            end if
        end do

    end SUBROUTINE ransamsinr2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! To sort an array into ascending order w.r.t specific column  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine  sort1(x,xx,n,ppp)
    implicit none
    integer, intent(in) :: n, ppp
    real (C_DOUBLE), dimension(n), intent(in) :: x
    real (C_DOUBLE), dimension(n, ppp), intent(inout) :: xx
    integer :: i, j, location
    real (C_DOUBLE),dimension(ppp):: TT
    real (C_DOUBLE) :: minimum

        do i = 1, n-1
            minimum  = x(i)
            location = i
            do j = i+1, n
                if (x(j) < minimum) then
                    minimum  = x(j)
                    location = j
                end if
            end do
            TT  = xx(i, :)
            xx(i, :) = xx(location, :)
            xx(location, :) = TT
        end do
    end subroutine  sort1

end module sinr2
