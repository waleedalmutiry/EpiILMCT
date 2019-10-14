!######################################################################
!# MODULE: sir
!# AUTHORS: 
!#     Waleed Almutiry <walmutir@uoguelph.ca>, 
!#     Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca> and
!#     Rob Deardon <robert.deardon@ucalgary.ca> 
!# 
!# DESCRIPTION:
!#
!#     MCMC tools for SIR continuous-time ILMs:
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
!#           mcmcsir ..................... subroutine
!#           likelihoodSIR2 ............... subroutine
!#           halfnormalden11 ............... function
!#           gammadensity11 ................ function
!#           uniformden11 .................. function
!#           randnormal11 .................. function
!#           randgamma11 ................... function
!#           initrandomseedsir11 ........... subroutine
!#           ransamsir11 ................... subroutine
!#
!######################################################################
module sir11
!    use, intrinsic :: iso_c_binding
    use ISO_C_BINDING
    implicit none
    private
    public :: mcmcsir_f

contains


    subroutine mcmcsir_f(n, nsim, ni, num, anum2, temp, nsuspar, ntranspar, net, dis, epidat, blockupdate, &
    & priordistsuspar, priordisttranspar, priordistkernelparpar, priordistsparkpar, &
    & priordistpowersus, priordistpowertrans, &
    & suspar, suscov, powersus, transpar, transcov, powertrans, kernelpar, spark, delta1, &
    & kernelparproposalvar, sparkproposalvar, susproposalvar, powersusproposalvar, &
    & transproposalvar, powertransproposalvar, infperiodproposal, &
    & priorpar1sus, priorpar2sus, &
    & priorpar1powersus, priorpar2powersus, &
    & priorpar1trans, priorpar2trans, &
    & priorpar1powertrans, priorpar2powertrans, &
    & kernelparprior, sparkprior, delta2prior, susparop, powersusparop, transparop, powertransparop, &
    & kernelparop, sparkop, &
    & delta2op, epidatmctim, epidatmcrem, loglik) bind(C, name="mcmcsir_f_")

    external infinity_value

    !Declarations
    integer(kind = C_INT), intent(in) :: n, nsim, ni, nsuspar, ntranspar, num
    integer(kind = C_INT), intent(in) :: temp
! to inndicate the function for updating each parameter or not:
    integer(kind = C_INT), intent(in), dimension(8) :: anum2
! prior distributions:
    integer(kind = C_INT), intent(in), dimension(ntranspar) :: priordisttranspar, priordistpowertrans
    integer(kind = C_INT), intent(in), dimension(nsuspar) :: priordistsuspar, priordistpowersus
    integer(kind = C_INT), intent(in), dimension(2) :: priordistkernelparpar
    integer(kind = C_INT), intent(in) :: priordistsparkpar
! for updating infection times/infectious periods in blocks:
    integer(kind = C_INT), intent(in), dimension(2) :: blockupdate
    integer(kind = C_INT)  :: obs_inf_time, sizeblock
! epidemic data:
    real(kind = C_DOUBLE), intent(in), dimension(n, 4) :: epidat
    real(kind = C_DOUBLE), dimension(n, 4) :: epidat3current, epidat3update, epidat3
! distance and network matrices:
    real(kind = C_DOUBLE), intent(in), dimension(n, n) :: dis, net
! susceptibility and transmissibility covariates:
    real(kind = C_DOUBLE), intent(in), dimension(n, nsuspar) :: suscov
    real(kind = C_DOUBLE), intent(in), dimension(n, ntranspar) :: transcov
! The variance of the normal proposal distributions and the parameters of the prior distributions of the model parameters:
    real(kind = C_DOUBLE), intent(in), dimension(ntranspar) :: transproposalvar, priorpar1trans, priorpar2trans
    real(kind = C_DOUBLE), intent(in), dimension(ntranspar) :: powertransproposalvar, priorpar1powertrans, priorpar2powertrans
    real(kind = C_DOUBLE), intent(in), dimension(nsuspar) :: susproposalvar, priorpar1sus, priorpar2sus
    real(kind = C_DOUBLE), intent(in), dimension(nsuspar) :: powersusproposalvar, priorpar1powersus, priorpar2powersus
    real(kind = C_DOUBLE), intent(in), dimension(2) :: sparkprior, kernelparproposalvar
    real(kind = C_DOUBLE), intent(in)  :: sparkproposalvar
    real(kind = C_DOUBLE), intent(in), dimension(2, 2) :: kernelparprior
! parameters of the prior distribution of the infectious period rate and gamma proposal distribution of the independenc sampler:
    real(kind = C_DOUBLE), intent(in), dimension(2) :: delta2prior, infperiodproposal
! initial values
    real(kind = C_DOUBLE), intent(in), dimension(ntranspar) :: transpar, powertrans
    real(kind = C_DOUBLE), intent(in), dimension(nsuspar) :: suspar, powersus
    real(kind = C_DOUBLE), intent(in), dimension(2) :: kernelpar
    real(kind = C_DOUBLE), intent(in) :: spark
    real(kind = C_DOUBLE), intent(in), dimension(2) :: delta1
! OUTPUT:
    real(kind = C_DOUBLE), dimension(nsim, 1), intent(out) :: delta2op
    real(kind = C_DOUBLE), dimension(nsim, n), intent(out) :: epidatmctim, epidatmcrem
    real(kind = C_DOUBLE), dimension(nsim, 1), intent(out) :: sparkop, loglik
    real(kind = C_DOUBLE), dimension(nsim, 2), intent(out) :: kernelparop
    real(kind = C_DOUBLE), dimension(nsim, ntranspar), intent(out) :: transparop, powertransparop
    real(kind = C_DOUBLE), dimension(nsim, nsuspar), intent(out) :: susparop, powersusparop

    real(kind = C_DOUBLE), dimension(nsuspar) :: postparsusupdate, postparsuscurrent
    real(kind = C_DOUBLE), dimension(ntranspar) :: postpartransupdate, postpartranscurrent

    real(kind = C_DOUBLE) :: postsparkupdate
    real(kind = C_DOUBLE), dimension(2) :: postkernelparupdate, postkernelparcurrent
    real(kind = C_DOUBLE), dimension(2) :: infperiodproposal1

    real(kind = C_DOUBLE) :: ratio1, ratio2, ratio3, ratio4, ratio5
    real(kind = C_DOUBLE) :: psi1, psi2, psi3, psi4, psi5
    real(kind = C_DOUBLE) :: u1, u2, u3, u4, u5
    real(kind = C_DOUBLE) :: valueupdate1, valuecurrent1
    real(kind = C_DOUBLE) :: valueupdate2, valuecurrent2
    real(kind = C_DOUBLE) :: valueupdate3, valuecurrent3
    real(kind = C_DOUBLE) :: valueupdate4, valuecurrent4
    real(kind = C_DOUBLE) :: valueupdate5, valuecurrent5
    real(kind = C_DOUBLE) :: valueffupdate, valueffcurrent
    real(kind = C_DOUBLE) :: valuefupdate1, valuefcurrent1
    real(kind = C_DOUBLE) :: valuefupdate2, valuefcurrent2
    real(kind = C_DOUBLE) :: valuefupdate3, valuefcurrent3
    real(kind = C_DOUBLE) :: valuefupdate4, valuefcurrent4
    real(kind = C_DOUBLE) :: valuefupdate5, valuefcurrent5, denupdate, dencurrent
    real(kind = C_DOUBLE) :: zparsus, zpartrans, zkernelpar, zspark
    integer(kind = C_INT) :: i, j, r, nblock, m, mn
    integer(kind = C_INT) :: abs32, abs33, abs34, abs35
    real(kind = C_DOUBLE) :: Inf, loglik66, infdens

    integer(kind = C_INT), allocatable, dimension(:, :)  :: amcmc
    integer(kind = C_INT), allocatable, dimension(:)  :: amcmc1
    integer(kind = C_INT), dimension(ni)  :: it, xxindic

    call infinity_value(Inf)

    if(temp .ne. 0) then
       call initrandomseedsir11(temp)
    end if

        obs_inf_time = blockupdate(1)
        sizeblock    = blockupdate(2)

        valuefupdate1 = 0.0_C_DOUBLE
        valuefcurrent1 = 0.0_C_DOUBLE
        valueupdate1 = 0.0_C_DOUBLE
        valuecurrent1 = 0.0_C_DOUBLE

        valuefupdate2 = 0.0_C_DOUBLE
        valuefcurrent2 = 0.0_C_DOUBLE
        valueupdate2 = 0.0_C_DOUBLE
        valuecurrent2 = 0.0_C_DOUBLE

        valuefupdate3 = 0.0_C_DOUBLE
        valuefcurrent3 = 0.0_C_DOUBLE
        valueupdate3 = 0.0_C_DOUBLE
        valuecurrent3 = 0.0_C_DOUBLE

        valuefupdate4 = 0.0_C_DOUBLE
        valuefcurrent4 = 0.0_C_DOUBLE
        valueupdate4 = 0.0_C_DOUBLE
        valuecurrent4 = 0.0_C_DOUBLE

        valuefupdate5 = 0.0_C_DOUBLE
        valuefcurrent5 = 0.0_C_DOUBLE
        valueupdate5 = 0.0_C_DOUBLE
        valuecurrent5 = 0.0_C_DOUBLE

        ratio1 = 0.0_C_DOUBLE
        psi1 = 0.0_C_DOUBLE
        ratio2 = 0.0_C_DOUBLE
        psi2 = 0.0_C_DOUBLE
        ratio3 = 0.0_C_DOUBLE
        psi3 = 0.0_C_DOUBLE
        ratio4 = 0.0_C_DOUBLE
        psi4 = 0.0_C_DOUBLE
        ratio5 = 0.0_C_DOUBLE
        psi5 = 0.0_C_DOUBLE

        valueffupdate = 0.0_C_DOUBLE
        valueffcurrent = 0.0_C_DOUBLE

! Initial values for the MCMC updates:

        delta2op(1, 1) = delta1(2)

        if (anum2(1) .eq. 1) then
        epidat3(:, 1)       = epidat(:, 1)
        epidat3(:, 2)       = epidat(:, 2)
        epidat3(:, 3)       = epidat(:, 3)
        epidat3(:, 4)       = epidat(:, 4)

        else if (anum2(1) .eq. 2) then
        epidat3 = epidat
        end if

        epidatmctim(1, 1:n) = (/epidat3(1:n, 4)/)
        epidatmcrem(1, 1:n) = (/epidat3(1:n, 2)/)

        kernelparop(1, :) = kernelpar
        sparkop(1, 1) = spark
        transparop(1, :) = transpar
        susparop(1, :) = suspar
        powertransparop(1, :) = powertrans
        powersusparop(1, :) = powersus

! calculate the loglikelihood with the initial values:
    if (anum2(1) .eq. 1) then

        infdens = 0.0_C_DOUBLE
        do i = 1, ni
            if (i .gt. obs_inf_time) then
                infdens = infdens + gammadensity11(epidat3(i, 3), delta1(1), delta2op(1, 1))
            else
                infdens = infdens
            end if
        end do

        call likelihoodSIR2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
        & suscov, transcov, susparop(1, :), powersusparop(1, :), &
        & transparop(1, :), powertransparop(1, :), kernelparop(1, :), sparkop(1, 1), loglik66)

        loglik(1, 1) = loglik66 + infdens

    else if (anum2(1) .eq. 2) then

        call likelihoodSIR2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
        & suscov, transcov, susparop(1, :), powersusparop(1, :), &
        & transparop(1, :), powertransparop(1, :), kernelparop(1, :), sparkop(1, 1), loglik66)

        loglik(1, 1) = loglik66

    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MCMC start
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do j = 1, (nsim-1)


            call likelihoodSIR2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
            & suscov, transcov, susparop(j, :), powersusparop(j, :), &
            & transparop(j, :), powertransparop(j, :), kernelparop(j, :), sparkop(j, 1), valuecurrent1)

            if (anum2(1) .eq. 1) then
! if the infection times and infectious periods are assumed unknown:

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: delta2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                delta2op(j+1, 1) = randgamma11(((dble(ni)*delta1(1))+delta2prior(1)), &
                & 1.0_C_DOUBLE/(delta2prior(2)+sum(epidat3(1:ni, 3))))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: infection times and infectious period:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                if (infperiodproposal(1) .eq. 0.0d0) then
                    infperiodproposal1 = (/ delta1(1), delta2op(j+1, 1) /)
                else
                    infperiodproposal1 = infperiodproposal
                end if


                nblock = floor(real(ni-obs_inf_time)/real(sizeblock))

                it = (/(r, r = 1, ni)/)

                xxindic = 0
                xxindic(1:obs_inf_time) = 1

                epidat3update  = epidat3      ! proposing
                epidat3current = epidat3      ! current


                    call likelihoodSIR2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
                        & suscov, transcov, susparop(j, :), powersusparop(j, :), &
                        & transparop(j, :), powertransparop(j, :), kernelparop(j, :), sparkop(j, 1), valuecurrent1)

! if the sizeblock has a value greater than 1, then update infection times/infectious periods
! in blocks with assigning those in each block randomly at each iteration.

                allocate(amcmc(sizeblock, nblock))
                allocate(amcmc1((ni-obs_inf_time)-(sizeblock*nblock)))

                amcmc  = 0
                amcmc1 = 0

                do r = 1, nblock
                    call ransamsir11(pack(it, (xxindic .ne. 1)), amcmc(:, r), count(xxindic .ne. 1), sizeblock)
                    xxindic(amcmc(:, r)) = 1
                end do

                amcmc1 = pack(it, (xxindic .ne. 1))

                mn = (ni-obs_inf_time)-(sizeblock*nblock)


                do r = 1, (nblock)

                    do m = 1, sizeblock
                        epidat3update(amcmc(m, r), 3)  = randgamma11(infperiodproposal1(1), 1.0_C_DOUBLE/infperiodproposal1(2))
                        epidat3update(amcmc(m, r), 4)  = epidat3update(amcmc(m, r), 2) - epidat3update(amcmc(m, r), 3)
                    end do

! calculate the log-likelihood function:

                    call likelihoodSIR2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3update, &
                        & suscov, transcov, susparop(j, :), powersusparop(j, :), &
                        & transparop(j, :), powertransparop(j, :), kernelparop(j, :), sparkop(j, 1), valueupdate1)

                    denupdate      = 0.0_C_DOUBLE
                    dencurrent     = 0.0_C_DOUBLE
                    valuefupdate1  = 0.0_C_DOUBLE
                    valuefcurrent1 = 0.0_C_DOUBLE

                    do m = 1, sizeblock
                        denupdate      = denupdate     + &
                        & gammadensity11(epidat3update(amcmc(m, r), 3), delta1(1), delta2op(j+1, 1))
                        dencurrent     = dencurrent    + &
                        & gammadensity11(epidat3current(amcmc(m, r), 3), delta1(1), delta2op(j+1, 1))
                        valuefupdate1  = valuefupdate1 + &
                        & gammadensity11(epidat3update(amcmc(m, r), 3), infperiodproposal1(1), infperiodproposal1(2))
                        valuefcurrent1 = valuefcurrent1+ &
                        & gammadensity11(epidat3current(amcmc(m, r), 3), infperiodproposal1(1), infperiodproposal1(2))
                    end do

                    valueffupdate  = valueupdate1  + denupdate 
                    valueffcurrent = valuecurrent1 + dencurrent

!calculating the log acceptance ratio:
                    ratio1    = (valueffupdate+valuefcurrent1)-(valueffcurrent+valuefupdate1)

                    psi1      = min(1.0_C_DOUBLE, dexp(ratio1))

! decision: accept/reject the proposed value:
                    call random_number(u1)

                    if ( (u1 .le. psi1) .and. (any(epidat3update(amcmc(:, r), 4) .lt. &
                    & epidat3current(obs_inf_time, 4)).eqv. .false.) ) then
                        epidat3current(amcmc(:, r), 3) = epidat3update(amcmc(:, r), 3)
                        epidat3current(amcmc(:, r), 4) = epidat3update(amcmc(:, r), 4)
                        valuecurrent1 = valueupdate1
                    else
                        epidat3update(amcmc(:, r), 3) = epidat3current(amcmc(:, r), 3)
                        epidat3update(amcmc(:, r), 4) = epidat3current(amcmc(:, r), 4)
                        valuecurrent1 = valuecurrent1
                    end if

                end do

! in case the last block has less number of infection times than the sizeblock:

                if (mn .ne. 0) then

                    do m = 1, mn
                        epidat3update(amcmc1(m), 3)  = randgamma11(infperiodproposal1(1), 1.0_C_DOUBLE/infperiodproposal1(2))
                        epidat3update(amcmc1(m), 4)  = epidat3update(amcmc1(m), 2)- epidat3update(amcmc1(m), 3)
                    end do

                    call likelihoodSIR2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3update, &
                        & suscov, transcov, susparop(j, :), powersusparop(j, :), &
                        & transparop(j, :), powertransparop(j, :), kernelparop(j, :), sparkop(j, 1), valueupdate1)

                    denupdate      = 0.0_C_DOUBLE
                    dencurrent     = 0.0_C_DOUBLE
                    valuefupdate1  = 0.0_C_DOUBLE
                    valuefcurrent1 = 0.0_C_DOUBLE

                    do m = 1, mn
                        denupdate  = denupdate  + &
                        & gammadensity11(epidat3update(amcmc1(m), 3), delta1(1), delta2op(j+1, 1))
                        dencurrent = dencurrent + &
                        & gammadensity11(epidat3current(amcmc1(m), 3), delta1(1), delta2op(j+1, 1))
                        valuefupdate1  = valuefupdate1 + &
                        & gammadensity11(epidat3update(amcmc1(m), 3), infperiodproposal1(1), infperiodproposal1(2))
                        valuefcurrent1 = valuefcurrent1+ &
                        & gammadensity11(epidat3current(amcmc1(m), 3), infperiodproposal1(1), infperiodproposal1(2))
                    end do

                    valueffupdate  = valueupdate1  + denupdate 
                    valueffcurrent = valuecurrent1 + dencurrent

!calculating the log acceptance ratio:
                    ratio1    = (valueffupdate+valuefcurrent1)-(valueffcurrent+valuefupdate1)

                    psi1      = min(1.0_C_DOUBLE, dexp(ratio1))

! decision: accept/reject the proposed value:
                    call random_number(u1)

                    if ( (u1 .le. psi1) .and. (any(epidat3update(amcmc1, 4) .lt. &
                    & epidat3current(obs_inf_time, 4)).eqv. .false.) ) then
                        epidat3current(amcmc1(1:mn), 3) = epidat3update(amcmc1(1:mn), 3)
                        epidat3current(amcmc1(1:mn), 4) = epidat3update(amcmc1(1:mn), 4)
                        valuecurrent1 = valueupdate1
                    else
                        epidat3update(amcmc1(1:mn), 3) = epidat3current(amcmc1(1:mn), 3)
                        epidat3update(amcmc1(1:mn), 4) = epidat3current(amcmc1(1:mn), 4)
                        valuecurrent1 = valuecurrent1
                    end if

                end if

                deallocate(amcmc)
                deallocate(amcmc1)

                epidat3 = epidat3current

                epidatmctim(j+1, 1:n) = (/epidat3(1:n, 4)/)
                epidatmcrem(j+1, 1:n) = (/epidat3(1:n, 2)/)

! if the infection times/infectious period are assumed known:

            else if (anum2(1) .eq. 2) then
                delta2op(j+1, 1) = delta1(2)
                epidat3 = epidat
                epidatmctim(j+1, 1:n) = (/epidat3(:, 4)/)
                epidatmcrem(j+1, 1:n) = (/epidat3(:, 2)/)
                valuecurrent1 = valuecurrent1
            end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: susceptibility parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (anum2(2) .eq. 1) then

                postparsusupdate    = susparop(j, :)
                postparsuscurrent   = susparop(j, :)

                do r = 1, nsuspar

! proposing candidate and making sure is positive:
                    zparsus = 0.0_C_DOUBLE
                    abs32 = 0
                    do while(abs32 .eq. 0)
                        zparsus = randnormal11(postparsuscurrent(r), susproposalvar(r))
                        if (zparsus .gt. 0.0_C_DOUBLE) then
                            abs32 = 1
                        else
                            abs32 = 0
                        end if
                    end do

                    postparsusupdate(r)= zparsus

!calculating the log-likelihood function:

                    call likelihoodSIR2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
                        & suscov, transcov, (/postparsusupdate/), powersusparop(j, :), &
                        & (/transparop(j, :)/), powertransparop(j, :), kernelparop(j, :), sparkop(j, 1), valueupdate2)


!calculating the log-prior:

                    if (priordistsuspar(r) .eq. 1) then
                        valuefupdate2    = gammadensity11(postparsusupdate(r), &
                        & priorpar1sus(r), priorpar2sus(r))
                        valuefcurrent2   = gammadensity11(postparsuscurrent(r), &
                        & priorpar1sus(r), priorpar2sus(r))
                    else if (priordistsuspar(r) .eq. 2) then
                        valuefupdate2    = halfnormalden11(postparsusupdate(r), &
                        & priorpar1sus(r), priorpar2sus(r))
                        valuefcurrent2   = halfnormalden11(postparsuscurrent(r), &
                        & priorpar1sus(r), priorpar2sus(r))
                    else if (priordistsuspar(r) .eq. 3) then
                        valuefupdate2    = uniformden11(postparsusupdate(r), &
                        & priorpar1sus(r), priorpar2sus(r))
                        valuefcurrent2   = uniformden11(postparsuscurrent(r), &
                        & priorpar1sus(r), priorpar2sus(r))
                    end if

!calculating the log acceptance ratio:
                    ratio2 = (valueupdate2 + valuefupdate2) - (valuecurrent1 + valuefcurrent2)

                    psi2 = min(1.0_C_DOUBLE, dexp(ratio2))

! decision: accept/reject the proposed value:

                    call random_number(u2)

                    if (psi2 .ge. u2) then
                        postparsuscurrent(r) = postparsusupdate(r)
                        valuecurrent1 = valueupdate2
                    else
                        postparsusupdate(r)  = postparsuscurrent(r)
                        valuecurrent1 = valuecurrent1
                    end if

                end do   !  end loop r nparsus
                susparop(j+1, :)  = postparsuscurrent

! if it is not included in the model:

            else
                susparop(j+1, 1:nsuspar)  = susparop(j, 1:nsuspar)
                valuecurrent1 = valuecurrent1
            end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: susceptibility power parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (anum2(6) .eq. 1) then

                postparsusupdate    = powersusparop(j, :)
                postparsuscurrent   = powersusparop(j, :)

                do r = 1, nsuspar

! proposing candidate and making sure is positive:
                    zparsus = 0.0_C_DOUBLE
                    abs32 = 0
                    do while(abs32 .eq. 0)
                        zparsus = randnormal11(postparsuscurrent(r), powersusproposalvar(r))
                        if (zparsus .gt. 0.0_C_DOUBLE) then
                            abs32 = 1
                        else
                            abs32 = 0
                        end if
                    end do

                    postparsusupdate(r)= zparsus

!calculating the log-likelihood function:

                    call likelihoodSIR2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
                        & suscov, transcov, susparop(j+1, :), postparsusupdate, &
                        & transparop(j, :), powertransparop(j, :), kernelparop(j, :), sparkop(j, 1), valueupdate2)


!calculating the log-prior:

                    if (priordistpowersus(r) .eq. 1) then
                        valuefupdate2    = gammadensity11(postparsusupdate(r), &
                        & priorpar1powersus(r), priorpar2powersus(r))
                        valuefcurrent2   = gammadensity11(postparsuscurrent(r), &
                        & priorpar1powersus(r), priorpar2powersus(r))
                    else if (priordistpowersus(r) .eq. 2) then
                        valuefupdate2    = halfnormalden11(postparsusupdate(r), &
                        & priorpar1powersus(r), priorpar2powersus(r))
                        valuefcurrent2   = halfnormalden11(postparsuscurrent(r), &
                        & priorpar1powersus(r), priorpar2powersus(r))
                    else if (priordistpowersus(r) .eq. 3) then
                        valuefupdate2    = uniformden11(postparsusupdate(r), &
                        & priorpar1powersus(r), priorpar2powersus(r))
                        valuefcurrent2   = uniformden11(postparsuscurrent(r), &
                        & priorpar1powersus(r), priorpar2powersus(r))
                    end if

!calculating the log acceptance ratio:
                    ratio2 = (valueupdate2 + valuefupdate2) - (valuecurrent1 + valuefcurrent2)

                    psi2 = min(1.0_C_DOUBLE, dexp(ratio2))

! decision: accept/reject the proposed value:

                    call random_number(u2)

                    if (psi2 .ge. u2) then
                        powersusparop(j+1, r) = postparsusupdate(r)
                        postparsuscurrent(r)  = postparsusupdate(r)
                        valuecurrent1 = valueupdate2
                    else
                        powersusparop(j+1, r) = powersusparop(j, r)
                        postparsusupdate(r)  = powersusparop(j, r)
                        valuecurrent1 = valuecurrent1
                    end if

                end do   !  end loop r nparsus
!                powersusparop(j+1, :)  = postparsuscurrent

! if it is not included in the model:
            else
                powersusparop(j+1, 1:nsuspar)  = powersusparop(j, 1:nsuspar)
                valuecurrent1 = valuecurrent1
            end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: transmissibility parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (anum2(3) .eq. 1) then

            postpartransupdate    = transparop(j, :)
            postpartranscurrent   = transparop(j, :)

            do r = 1, ntranspar

! proposing candidate:
                zpartrans = 0.0_C_DOUBLE

                abs33 = 0

                do while(abs33 .eq. 0)
                    zpartrans = randnormal11(postpartranscurrent(r), transproposalvar(r))
                    if (zpartrans .gt. 0.0_C_DOUBLE) then
                        abs33 = 1
                    else
                        abs33 = 0
                    end if
                end do

                postpartransupdate(r)= zpartrans

!calculating the log-likelihoodSIR2:

                call likelihoodSIR2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
                    & suscov, transcov, susparop(j+1, :), powersusparop(j+1, :), &
                    & postpartransupdate, powertransparop(j, :), kernelparop(j, :), sparkop(j, 1), valueupdate3)

!calculating the log-prior:

                if (priordisttranspar(r) .eq. 1) then
                    valuefupdate3    = gammadensity11(postpartransupdate(r), &
                    & priorpar1trans(r), priorpar2trans(r))
                    valuefcurrent3   = gammadensity11(postpartranscurrent(r), &
                    & priorpar1trans(r), priorpar2trans(r))
                else if (priordisttranspar(r) .eq. 2) then
                    valuefupdate3    = halfnormalden11(postpartransupdate(r), &
                    & priorpar1trans(r), priorpar2trans(r))
                    valuefcurrent3   = halfnormalden11(postpartranscurrent(r), &
                    & priorpar1trans(r), priorpar2trans(r))
                else
                    valuefupdate3    = uniformden11(postpartransupdate(r), &
                    & priorpar1trans(r), priorpar2trans(r))
                    valuefcurrent3   = uniformden11(postpartranscurrent(r), &
                    & priorpar1trans(r), priorpar2trans(r))
                end if

!calculating the log acceptance ratio:
                ratio3 = (valueupdate3 + valuefupdate3) - (valuecurrent1 + valuefcurrent3)

                psi3 = min(1.0_C_DOUBLE, dexp(ratio3))

! decision: accept/reject the proposed value:

                call random_number(u3)

                if (psi3 .ge. u3) then
                    postpartranscurrent(r) = postpartransupdate(r)
                    postpartransupdate(r)  = postpartransupdate(r)
                    valuecurrent1 = valueupdate3
                else
                    postpartransupdate(r)  = postpartranscurrent(r)
                    postpartranscurrent(r)  = postpartranscurrent(r)
                    valuecurrent1 = valuecurrent1
                end if

            end do   !  end loop r nparsus

            transparop(j+1, :)        = postpartranscurrent

! if it is not included in the model:
            else
                transparop(j+1, 1:ntranspar)        = transparop(j, 1:ntranspar)
                valuecurrent1 = valuecurrent1
            end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: transmissibility power parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (anum2(7) .eq. 1) then

            postpartransupdate    = powertransparop(j, :)
            postpartranscurrent   = powertransparop(j, :)

            do r = 1, ntranspar

! proposing candidate:
                zpartrans = 0.0_C_DOUBLE

                abs33 = 0

                do while(abs33 .eq. 0)
                    zpartrans = randnormal11(postpartranscurrent(r), powertransproposalvar(r))
                    if (zpartrans .gt. 0.0_C_DOUBLE) then
                        abs33 = 1
                    else
                        abs33 = 0
                    end if
                end do

                postpartransupdate(r)= zpartrans

!calculating the log-likelihoodSIR2:

                call likelihoodSIR2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
                    & suscov, transcov, susparop(j+1, :), powersusparop(j+1, :), &
                    & transparop(j+1, :), postpartransupdate, kernelparop(j, :), sparkop(j, 1), valueupdate3)

!calculating the log-prior:

                if (priordistpowertrans(r) .eq. 1) then
                    valuefupdate3    = gammadensity11(postpartransupdate(r), &
                    & priorpar1powertrans(r), priorpar2powertrans(r))
                    valuefcurrent3   = gammadensity11(postpartranscurrent(r), &
                    & priorpar1powertrans(r), priorpar2powertrans(r))
                else if (priordistpowertrans(r) .eq. 2) then
                    valuefupdate3    = halfnormalden11(postpartransupdate(r), &
                    & priorpar1powertrans(r), priorpar2powertrans(r))
                    valuefcurrent3   = halfnormalden11(postpartranscurrent(r), &
                    & priorpar1powertrans(r), priorpar2powertrans(r))
                else if (priordistpowertrans(r) .eq. 3) then
                    valuefupdate3    = uniformden11(postpartransupdate(r), &
                    & priorpar1powertrans(r), priorpar2powertrans(r))
                    valuefcurrent3   = uniformden11(postpartranscurrent(r), &
                    & priorpar1powertrans(r), priorpar2powertrans(r))
                end if

!calculating the log acceptance ratio:
                ratio3 = (valueupdate3 + valuefupdate3) - (valuecurrent1 + valuefcurrent3)

                psi3 = min(1.0_C_DOUBLE, dexp(ratio3))

! decision: accept/reject the proposed value:

                call random_number(u3)

                if (psi3 .ge. u3) then
                    postpartranscurrent(r) = postpartransupdate(r)
                    postpartransupdate(r)  = postpartransupdate(r)
                    valuecurrent1 = valueupdate3
                else
                    postpartransupdate(r)  = postpartranscurrent(r)
                    postpartranscurrent(r)  = postpartranscurrent(r)
                    valuecurrent1 = valuecurrent1
                end if

            end do   !  end loop r nparsus

            powertransparop(j+1, :)        = postpartranscurrent

! if it is not included in the model:
            else
                powertransparop(j+1, 1:ntranspar)        = powertransparop(j, 1:ntranspar)
                valuecurrent1 = valuecurrent1
            end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: spark parameter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (anum2(4) .eq. 1) then

                abs34 = 0

                do while(abs34 .eq. 0)
                    zspark = randnormal11(sparkop(j, 1), sparkproposalvar)
                    if (zspark .gt. 0.0_C_DOUBLE) then
                        abs34 = 1
                    else
                        abs34 = 0
                    end if
                end do        ! end while abs3

                postsparkupdate = zspark

                call likelihoodSIR2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
                    & suscov, transcov, susparop(j+1, :), powersusparop(j+1, :), &
                    & transparop(j+1, :), powertransparop(j+1, :), kernelparop(j, :), postsparkupdate, valueupdate4)

!calculating the log-prior:
                if (priordistsparkpar .eq. 1) then
                    valuefupdate4    = gammadensity11(postsparkupdate, &
                    & sparkprior(1), sparkprior(2))
                    valuefcurrent4   = gammadensity11(sparkop(j, 1), &
                    & sparkprior(1), sparkprior(2))
                else if (priordistsparkpar .eq. 2) then
                    valuefupdate4    = halfnormalden11(postsparkupdate, &
                    & sparkprior(1), sparkprior(2))
                    valuefcurrent4   = halfnormalden11(sparkop(j, 1), &
                    & sparkprior(1), sparkprior(2))
                else
                    valuefupdate4    = uniformden11(postsparkupdate, &
                    & sparkprior(1), sparkprior(2))
                    valuefcurrent4   = uniformden11(sparkop(j, 1), &
                    & sparkprior(1), sparkprior(2))
                end if

!calculating the log acceptance ratio:
                ratio4 = (valueupdate4 + valuefupdate4)-(valuecurrent1 + valuefcurrent4)

                psi4 = min(1.0_C_DOUBLE, dexp(ratio4))

                ! decision: accept/reject the proposed value:

                call random_number(u4)

                if (psi4 .ge. u4) then
                    sparkop(j+1, 1) = postsparkupdate
                    valuecurrent1 = valueupdate4
                else
                    sparkop(j+1, 1) = sparkop(j, 1)
                    valuecurrent1 = valuecurrent1
                end if

                ! if it is not included in the model:
            else
                sparkop(j+1, 1) = sparkop(j, 1)
                valuecurrent1 = valuecurrent1

            end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE: kernelpar parameters (spatial and network effect parameters)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (anum2(5) .eq. 1) then
! for distance-based continuous-time ILMS:
! just the spatial parameter is needed to be updated

                abs35 = 0
                do while(abs35 .eq. 0)
                    zkernelpar = randnormal11(kernelparop(j, 1), kernelparproposalvar(1))
                    if (zkernelpar .gt. 0.0_C_DOUBLE) then
                        abs35 = 1
                    else
                        abs35 = 0
                    end if
                end do        ! end while abs3

                postkernelparupdate = (/ zkernelpar, kernelparop(j, 2) /)

! calculate the log-likelihood function:
                call likelihoodSIR2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
                    & suscov, transcov, susparop(j+1, :), powersusparop(j+1, :), &
                    & transparop(j+1, :), powertransparop(j+1, :), postkernelparupdate, sparkop(j+1, 1), valueupdate5)

!calculating the log-prior:
                if (priordistkernelparpar(1) .eq. 1) then
                    valuefupdate5    = gammadensity11(postkernelparupdate(1), &
                    & kernelparprior(1, 1), kernelparprior(1, 2))
                    valuefcurrent5   = gammadensity11(kernelparop(j, 1), &
                    & kernelparprior(1, 1), kernelparprior(1, 2))
                else if (priordistkernelparpar(1) .eq. 2) then
                    valuefupdate5    = halfnormalden11(postkernelparupdate(1), &
                    & kernelparprior(1, 1), kernelparprior(1, 2))
                    valuefcurrent5   = halfnormalden11(kernelparop(j, 1), &
                    & kernelparprior(1, 1), kernelparprior(1, 2))
                else
                    valuefupdate5    = uniformden11(postkernelparupdate(1), &
                    & kernelparprior(1, 1), kernelparprior(1, 2))
                    valuefcurrent5   = uniformden11(kernelparop(j, 1), &
                    & kernelparprior(1, 1), kernelparprior(1, 2))
                end if


!calculating the log acceptance ratio:
                ratio5 = (valueupdate5 + valuefupdate5)-(valuecurrent1 + valuefcurrent5)

                psi5 = min(1.0_C_DOUBLE, dexp(ratio5))

! decision: accept/reject the proposed value:

                call random_number(u5)

                if (psi5 .ge. u5) then
                    kernelparop(j+1, :) = postkernelparupdate
                    valuecurrent1 = valueupdate5
                else
                    kernelparop(j+1, :) = kernelparop(j, :)
                    valuecurrent1 = valuecurrent1
                end if



            else if (anum2(5) .eq. 2) then
! for combined network- and distance-based continuous-time ILMS:

! updating the spatial parameter:

                abs35 = 0
                do while(abs35 .eq. 0)
                    zkernelpar = randnormal11(kernelparop(j, 1), kernelparproposalvar(1))
                    if (zkernelpar .gt. 0.0_C_DOUBLE) then
                        abs35 = 1
                    else
                        abs35 = 0
                    end if
                end do        ! end while abs3

                postkernelparupdate = (/ zkernelpar, kernelparop(j, 2) /)
                postkernelparcurrent = kernelparop(j, :)

! calculate the log-likelihood function:
                call likelihoodSIR2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
                    & suscov, transcov, susparop(j+1, :), powersusparop(j+1, :), &
                    & transparop(j+1, :), powertransparop(j+1, :), postkernelparupdate, sparkop(j+1, 1), valueupdate5)

!calculating the log-prior:
                if (priordistkernelparpar(1) .eq. 1) then
                    valuefupdate5    = gammadensity11(postkernelparupdate(1), &
                    & kernelparprior(1, 1), kernelparprior(1, 2))
                    valuefcurrent5   = gammadensity11(kernelparop(j, 1), &
                    & kernelparprior(1, 1), kernelparprior(1, 2))
                else if (priordistkernelparpar(1) .eq. 2) then
                    valuefupdate5    = halfnormalden11(postkernelparupdate(1), &
                    & kernelparprior(1, 1), kernelparprior(1, 2))
                    valuefcurrent5   = halfnormalden11(kernelparop(j, 1), &
                    & kernelparprior(1, 1), kernelparprior(1, 2))
                else
                    valuefupdate5    = uniformden11(postkernelparupdate(1), &
                    & kernelparprior(1, 1), kernelparprior(1, 2))
                    valuefcurrent5   = uniformden11(kernelparop(j, 1), &
                    & kernelparprior(1, 1), kernelparprior(1, 2))
                end if


!calculating the log acceptance ratio:
                ratio5 = (valueupdate5 + valuefupdate5)-(valuecurrent1 + valuefcurrent5)

                psi5 = min(1.0_C_DOUBLE, dexp(ratio5))

! decision: accept/reject the proposed value:

                call random_number(u5)

                if (psi5 .ge. u5) then
                    kernelparop(j+1, 1) = postkernelparupdate(1)
                    valuecurrent1 = valueupdate5
                else
                    kernelparop(j+1, 1) = kernelparop(j, 1)
                    postkernelparupdate(1) = kernelparop(j, 1)
                    valuecurrent1 = valuecurrent1
                end if


! updating the network-effect parameter:
! proposing candidate and making sure is positive:

                abs35 = 0
                do while(abs35 .eq. 0)
                    zkernelpar = randnormal11(kernelparop(j, 2), kernelparproposalvar(2))
                    if (zkernelpar .gt. 0.0_C_DOUBLE) then
                        abs35 = 1
                    else
                        abs35 = 0
                    end if
                end do        ! end while abs3

                postkernelparupdate = (/ kernelparop(j+1, 1), zkernelpar /)
                postkernelparcurrent = (/ kernelparop(j+1, 1), kernelparop(j, 2) /)

! calculate the log-likelihood function:
                call likelihoodSIR2(n, ni, num, nsuspar, ntranspar, net, dis, epidat3, &
                    & suscov, transcov, susparop(j+1, :), powersusparop(j+1, :), &
                    & transparop(j+1, :), powertransparop(j+1, :), postkernelparupdate, sparkop(j+1, 1), valueupdate5)

!calculating the log-prior:
                if (priordistkernelparpar(2) .eq. 1) then
                    valuefupdate5    = gammadensity11(postkernelparupdate(2), &
                    & kernelparprior(2, 1), kernelparprior(2, 2))
                    valuefcurrent5   = gammadensity11(kernelparop(j, 2), &
                    & kernelparprior(2, 1), kernelparprior(2, 2))
                else if (priordistkernelparpar(2) .eq. 2) then
                    valuefupdate5    = halfnormalden11(postkernelparupdate(2), &
                    & kernelparprior(2, 1), kernelparprior(2, 2))
                    valuefcurrent5   = halfnormalden11(kernelparop(j, 2), &
                    & kernelparprior(2, 1), kernelparprior(2, 2))
                else
                    valuefupdate5    = uniformden11(postkernelparupdate(2), &
                    & kernelparprior(2, 1), kernelparprior(2, 2))
                    valuefcurrent5   = uniformden11(kernelparop(j, 2), &
                    & kernelparprior(2, 1), kernelparprior(2, 2))
                end if

!calculating the log acceptance ratio:
                ratio5 = (valueupdate5 + valuefupdate5)-(valuecurrent1 + valuefcurrent5)

                psi5 = min(1.0_C_DOUBLE, dexp(ratio5))

! decision: accept/reject the proposed value:
                call random_number(u5)

                if (psi5 .ge. u5) then
                    kernelparop(j+1, 2) = postkernelparupdate(2)
                    valuecurrent1 = valueupdate5
                else
                    kernelparop(j+1, 2) = kernelparop(j, 2)
                    postkernelparupdate(2) = kernelparop(j, 2)
                    valuecurrent1 = valuecurrent1
                end if


! for network-based continuous-time ILMS:
            else if (anum2(5) .eq. 3) then
                kernelparop(j+1, :) = kernelparop(j, :)
                valuecurrent1 = valuecurrent1
            end if

! to update the log-likelihood based on the updated parameters:
            if (anum2(1) .eq. 1) then

                infdens = 0.0_C_DOUBLE
                do i = 1, ni
                    if (i .gt. obs_inf_time) then
                        infdens = infdens + gammadensity11(epidat3(i, 3), delta1(1), &
                                & delta2op(j+1, 1))
                    else
                        infdens = infdens
                    end if
                end do

                loglik(j+1, 1) = valuecurrent1 + infdens

            else if (anum2(1) .eq. 2) then

                loglik(j+1, 1) = valuecurrent1

            end if

        end do !end loop j (mcmc)

    end subroutine mcmcsir_f






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 				Log likelihoodSIR2 subroutine			 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine likelihoodSIR2(n, ninfected, num, nsuspar, ntranspar, cc, d3, epidat, &
    & suscov, transcov, suspar, powersus, transpar, powertrans, kernelpar, spark, likk)

    external infinity_value

    integer :: j, i, r, m
    integer, intent(in) :: n, ninfected, num, nsuspar, ntranspar           !integers
    double precision, intent(in) :: spark                             ! spark parameter
    double precision, intent(in), dimension(2) :: kernelpar            ! parameter of the kernel function
    double precision, intent(in), dimension(n, n) :: d3, cc             ! distance & network matrices
    double precision, intent(in), dimension(n, 4) :: epidat              ! epidemic data
    double precision, dimension(n, 4) :: epidat1              ! epidemic data
    double precision, intent(in), dimension(n, nsuspar) :: suscov        ! susceptibility covariates
    double precision, intent(in), dimension(n, ntranspar) :: transcov    ! transmissibility covariates
    double precision, intent(in), dimension(nsuspar) :: suspar, powersus ! susceptibility parameters
    double precision, intent(in), dimension(ntranspar) :: transpar, powertrans ! transmissibility parameters
    double precision, intent(out) :: likk                              ! OUTPUT
    double precision :: rate, tt, ss, qa1, qa2, likk1, likk2, Inf
    double precision, dimension(n)  :: gh, sss
    double precision, dimension(ninfected)  :: rt, df, rrate

! defining Infinity
    call infinity_value(Inf)

    SELECT CASE (num)
    CASE (1)

! for contact network-based ILM.

! getting t_obs (maximum removal time):
        
        epidat1 = epidat
        
        call Sort(epidat1(:,4),epidat1,n,4)

        tt = maxval(epidat1(1:ninfected, 2))

!calculate the exponent part of the spark term:
        do r = 1, ninfected
            sss(r) = spark * (min(tt, epidat1(r, 4)) - epidat1(1, 4))
        end do
        do r = (ninfected+1), n
            sss(r) = spark * (tt - epidat1(1, 4))
        end do
        ss = sum(sss)

! calculate the terms of the exponent part of the double summation:
        df = 0.0d0

        do i =1, ninfected
        
            qa2 = 0.0d0
            do m = 1, ntranspar
                qa2 = qa2 + (transpar(m) * (transcov(int(epidat1(i, 1)), m)**powertrans(m)))
            end do
                
            gh = 0.0d0
                 
            do j = (i+1), n
            
                qa1 = 0.0d0
                        
                do m = 1, nsuspar
                    qa1 = qa1 + (suspar(m) * (suscov(int(epidat1(j, 1)), m)**powersus(m)))
                end do

                rate = qa1 * qa2 * (cc(int(epidat1(i, 1)), int(epidat1(j, 1))))

                if (j .le. ninfected) then
                    gh(j) = ( min(epidat1(i, 2), epidat1(j, 4)) - min(epidat1(i, 4), epidat1(j, 4)) ) * (rate)
                else
                    gh(j) = ( epidat1(i, 2) - epidat1(i, 4)) * (rate)
                end if

            end do
            df(i) = sum(gh)
        end do

        likk1 = -(sum(df)+ss)

! calculate the first part of the likelihood:
        rt(1) = 1.0d0

        do j = 2, ninfected

            qa1 = 0.0d0

            do m = 1, nsuspar
                qa1 = qa1 + (suspar(m) * (suscov(int(epidat1(j, 1)), m)**powersus(m)))
            end do

            rrate = 0.0d0

            do i = 1, (j-1)

                    if ((epidat1(j, 4) .gt. epidat1(i, 4)) .and.  (epidat1(j, 4) .le. epidat1(i, 2))) then

!                    if ( (epidat1(j, 4) .le. epidat1(i, 2)) ) then

                        qa2 = 0.0d0
                        
                        do m = 1, ntranspar
                            qa2 = qa2 + (transpar(m) * (transcov(int(epidat1(i, 1)), m)**powertrans(m)))
                        end do

                        rrate(i) = qa1 * qa2 * (cc(int(epidat1(i, 1)), int(epidat1(j, 1))))

                    else

                        rrate(i)   = 0.0d0

                    end if

            end do
            
            rt(j)  = sum(rrate) + spark
            
        end do

        likk2 = sum(log(rt))

! TOTAL LOG-LIKELIHOOD without the density of the infectious periods:
        likk =  likk1 + likk2

!#########################################################
!#########################################################

    CASE (2)
! for distance-based ILM with power-law kernel.

! getting t_obs (maximum removal time):
        
        epidat1 = epidat
        
        call Sort(epidat1(:,4),epidat1,n,4)

        tt = maxval(epidat1(1:ninfected, 2))

!calculate the exponent part of the spark term:
        do r = 1, ninfected
            sss(r) = spark * (min(tt, epidat1(r, 4)) - epidat1(1, 4))
        end do
        do r = (ninfected+1), n
            sss(r) = spark * (tt - epidat1(1, 4))
        end do
        ss = sum(sss)

! calculate the terms of the exponent part of the double summation:
        df = 0.0d0

        do i =1, ninfected
        
            qa2 = 0.0d0
            do m = 1, ntranspar
                qa2 = qa2 + (transpar(m) * (transcov(int(epidat1(i, 1)), m)**powertrans(m)))
            end do
                
            gh = 0.0d0
                 
            do j = (i+1), n
            
                qa1 = 0.0d0
                        
                do m = 1, nsuspar
                    qa1 = qa1 + (suspar(m) * (suscov(int(epidat1(j, 1)), m)**powersus(m)))
                end do

                rate = qa1 * qa2 * (d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(-kernelpar(1)))

                if (j .le. ninfected) then
                    gh(j) = ( min(epidat1(i, 2), epidat1(j, 4)) - min(epidat1(i, 4), epidat1(j, 4)) ) * (rate)
                else
                    gh(j) = ( epidat1(i, 2) - epidat1(i, 4)) * (rate)
                end if

            end do
            df(i) = sum(gh)
        end do

        likk1 = -(sum(df)+ss)

! calculate the first part of the likelihood:
        rt(1) = 1.0d0

        do j = 2, ninfected

            qa1 = 0.0d0

            do m = 1, nsuspar
                qa1 = qa1 + (suspar(m) * (suscov(int(epidat1(j, 1)), m)**powersus(m)))
            end do

            rrate = 0.0d0

            do i = 1, (j-1)

                    if ((epidat1(j, 4) .gt. epidat1(i, 4)) .and.  (epidat1(j, 4) .le. epidat1(i, 2))) then

!                    if ( (epidat1(j, 4) .le. epidat1(i, 2)) ) then

                        qa2 = 0.0d0
                        
                        do m = 1, ntranspar
                            qa2 = qa2 + (transpar(m) * (transcov(int(epidat1(i, 1)), m)**powertrans(m)))
                        end do

                        rrate(i) = qa1 * qa2 * (d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(-kernelpar(1)))

                    else

                        rrate(i)   = 0.0d0

                    end if

            end do
            
            rt(j)  = sum(rrate) + spark
            
        end do

        likk2 = sum(log(rt))

! TOTAL LOG-LIKELIHOOD without the density of the infectious periods:
        likk =  likk1 + likk2

!#########################################################
!#########################################################

    CASE (3)
! for distance-based ILM with Cauchy kernel.

! getting t_obs (maximum removal time):
        
        epidat1 = epidat
        
        call Sort(epidat1(:,4),epidat1,n,4)

        tt = maxval(epidat1(1:ninfected, 2))

!calculate the exponent part of the spark term:
        do r = 1, ninfected
            sss(r) = spark * (min(tt, epidat1(r, 4)) - epidat1(1, 4))
        end do
        do r = (ninfected+1), n
            sss(r) = spark * (tt - epidat1(1, 4))
        end do
        ss = sum(sss)

! calculate the terms of the exponent part of the double summation:
        df = 0.0d0

        do i =1, ninfected
        
            qa2 = 0.0d0
            do m = 1, ntranspar
                qa2 = qa2 + (transpar(m) * (transcov(int(epidat1(i, 1)), m)**powertrans(m)))
            end do
                
            gh = 0.0d0
                 
            do j = (i+1), n
            
                qa1 = 0.0d0
                        
                do m = 1, nsuspar
                    qa1 = qa1 + (suspar(m) * (suscov(int(epidat1(j, 1)), m)**powersus(m)))
                end do

                rate = qa1 * qa2 * (kernelpar(1)/((d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(2.0d0))+ &
                        & (kernelpar(1)**(2.0d0))))
                        
                if (j .le. ninfected) then
                    gh(j) = ( min(epidat1(i, 2), epidat1(j, 4)) - min(epidat1(i, 4), epidat1(j, 4)) ) * (rate)
                else
                    gh(j) = ( epidat1(i, 2) - epidat1(i, 4)) * (rate)
                end if

            end do
            df(i) = sum(gh)
        end do

        likk1 = -(sum(df)+ss)

! calculate the first part of the likelihood:
        rt(1) = 1.0d0

        do j = 2, ninfected

            qa1 = 0.0d0

            do m = 1, nsuspar
                qa1 = qa1 + (suspar(m) * (suscov(int(epidat1(j, 1)), m)**powersus(m)))
            end do

            rrate = 0.0d0

            do i = 1, (j-1)

                    if ((epidat1(j, 4) .gt. epidat1(i, 4)) .and.  (epidat1(j, 4) .le. epidat1(i, 2))) then

                        qa2 = 0.0d0
                        
                        do m = 1, ntranspar
                            qa2 = qa2 + (transpar(m) * (transcov(int(epidat1(i, 1)), m)**powertrans(m)))
                        end do

                        rrate(i) = qa1 * qa2 * (kernelpar(1)/((d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(2.0d0))+ &
                        & (kernelpar(1)**(2.0d0))))
                        
                    else

                        rrate(i)   = 0.0d0

                    end if

            end do
            
            rt(j)  = sum(rrate) + spark
            
        end do

        likk2 = sum(log(rt))

! TOTAL LOG-LIKELIHOOD without the density of the infectious periods:
        likk =  likk1 + likk2

!#########################################################
!#########################################################

    CASE (4)
! for combined network- and distance-based ILM with power-law kernel.

! getting t_obs (maximum removal time):
        
        epidat1 = epidat
        
        call Sort(epidat1(:,4),epidat1,n,4)

        tt = maxval(epidat1(1:ninfected, 2))

!calculate the exponent part of the spark term:
        do r = 1, ninfected
            sss(r) = spark * (min(tt, epidat1(r, 4)) - epidat1(1, 4))
        end do
        do r = (ninfected+1), n
            sss(r) = spark * (tt - epidat1(1, 4))
        end do
        ss = sum(sss)

! calculate the terms of the exponent part of the double summation:
        df = 0.0d0

        do i =1, ninfected
        
            qa2 = 0.0d0
            do m = 1, ntranspar
                qa2 = qa2 + (transpar(m) * (transcov(int(epidat1(i, 1)), m)**powertrans(m)))
            end do
                
            gh = 0.0d0
                 
            do j = (i+1), n
            
                qa1 = 0.0d0
                        
                do m = 1, nsuspar
                    qa1 = qa1 + (suspar(m) * (suscov(int(epidat1(j, 1)), m)**powersus(m)))
                end do

                rate = qa1 * qa2 * ( (d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(-kernelpar(1))) + &
                        &  (kernelpar(2)*cc(int(epidat1(i, 1)), int(epidat1(j, 1)))) )
                        
                if (j .le. ninfected) then
                    gh(j) = ( min(epidat1(i, 2), epidat1(j, 4)) - min(epidat1(i, 4), epidat1(j, 4)) ) * (rate)
                else
                    gh(j) = ( epidat1(i, 2) - epidat1(i, 4)) * (rate)
                end if

            end do
            df(i) = sum(gh)
        end do

        likk1 = -(sum(df)+ss)

! calculate the first part of the likelihood:
        rt(1) = 1.0d0

        do j = 2, ninfected

            qa1 = 0.0d0

            do m = 1, nsuspar
                qa1 = qa1 + (suspar(m) * (suscov(int(epidat1(j, 1)), m)**powersus(m)))
            end do

            rrate = 0.0d0

            do i = 1, (j-1)

                    if ((epidat1(j, 4) .gt. epidat1(i, 4)) .and.  (epidat1(j, 4) .le. epidat1(i, 2))) then

                        qa2 = 0.0d0
                        
                        do m = 1, ntranspar
                            qa2 = qa2 + (transpar(m) * (transcov(int(epidat1(i, 1)), m)**powertrans(m)))
                        end do

                        rrate(i) = qa1 * qa2 * ( (d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(-kernelpar(1))) + &
                        &  (kernelpar(2)*cc(int(epidat1(i, 1)), int(epidat1(j, 1)))) )
                        
                    else

                        rrate(i)   = 0.0d0

                    end if

            end do
            
            rt(j)  = sum(rrate) + spark
            
        end do

        likk2 = sum(log(rt))

! TOTAL LOG-LIKELIHOOD without the density of the infectious periods:
        likk =  likk1 + likk2


!#########################################################
!#########################################################

    CASE (5)
! for combined network- and distance-based ILM with Cauchy kernel.

! getting t_obs (maximum removal time):
        
        epidat1 = epidat
        
        call Sort(epidat1(:,4),epidat1,n,4)

        tt = maxval(epidat1(1:ninfected, 2))

!calculate the exponent part of the spark term:
        do r = 1, ninfected
            sss(r) = spark * (min(tt, epidat1(r, 4)) - epidat1(1, 4))
        end do
        do r = (ninfected+1), n
            sss(r) = spark * (tt - epidat1(1, 4))
        end do
        ss = sum(sss)

! calculate the terms of the exponent part of the double summation:
        df = 0.0d0

        do i =1, ninfected
        
            qa2 = 0.0d0
            do m = 1, ntranspar
                qa2 = qa2 + (transpar(m) * (transcov(int(epidat1(i, 1)), m)**powertrans(m)))
            end do
                
            gh = 0.0d0
                 
            do j = (i+1), n
            
                qa1 = 0.0d0
                        
                do m = 1, nsuspar
                    qa1 = qa1 + (suspar(m) * (suscov(int(epidat1(j, 1)), m)**powersus(m)))
                end do

                rate = qa1 * qa2 * ( (kernelpar(1)/((d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(2.0d0))+ &
                        & kernelpar(1)**(2.0d0))) + (kernelpar(2)*cc(int(epidat1(i, 1)), int(epidat1(j, 1)))) )
                        
                if (j .le. ninfected) then
                    gh(j) = ( min(epidat1(i, 2), epidat1(j, 4)) - min(epidat1(i, 4), epidat1(j, 4)) ) * (rate)
                else
                    gh(j) = ( epidat1(i, 2) - epidat1(i, 4)) * (rate)
                end if

            end do
            df(i) = sum(gh)
        end do

        likk1 = -(sum(df)+ss)

! calculate the first part of the likelihood:
        rt(1) = 1.0d0

        do j = 2, ninfected

            qa1 = 0.0d0

            do m = 1, nsuspar
                qa1 = qa1 + (suspar(m) * (suscov(int(epidat1(j, 1)), m)**powersus(m)))
            end do

            rrate = 0.0d0

            do i = 1, (j-1)

                    if ((epidat1(j, 4) .gt. epidat1(i, 4)) .and.  (epidat1(j, 4) .le. epidat1(i, 2))) then

                        qa2 = 0.0d0
                        
                        do m = 1, ntranspar
                            qa2 = qa2 + (transpar(m) * (transcov(int(epidat1(i, 1)), m)**powertrans(m)))
                        end do

                        rrate(i) = qa1 * qa2 * ( (kernelpar(1)/((d3(int(epidat1(i, 1)), int(epidat1(j, 1)))**(2.0d0))+ &
                        & kernelpar(1)**(2.0d0))) + (kernelpar(2)*cc(int(epidat1(i, 1)), int(epidat1(j, 1)))) )
                        
                    else

                        rrate(i)   = 0.0d0

                    end if

            end do
            
            rt(j)  = sum(rrate) + spark
            
        end do

        likk2 = sum(log(rt))

! TOTAL LOG-LIKELIHOOD without the density of the infectious periods:
        likk =  likk1 + likk2

    END SELECT

    end subroutine likelihoodSIR2


!#########################################################
!#########################################################



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     Densities functions for diffierent distributions     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#################### HALF-NORMAL density ######################

    FUNCTION halfnormalden11(alpha, a, b) RESULT(pdf)
    implicit none
    double precision, parameter :: pi = 3.141592653589793D+00
    double precision  :: alpha, b, a
    double precision  :: val, pdf

        if (alpha .lt. a) then
            val = 0.0d0
        else
            val =  sqrt (2 / pi*b) * exp ( - 0.5D+00 * ((alpha)**2 / b) )
        end if

        pdf= log(val)

    END FUNCTION halfnormalden11 ! returns log value

!#################### GAMMA density ######################

    FUNCTION gammadensity11(x, a, b) RESULT(pdf1)
    implicit none
    double precision  :: x, a, b
    double precision  :: dn, pdf1

        if (x .le. 0.0d0) then
            dn = 0.0d0
        else
            dn = (x**(a-1.0d0)) * dexp(- (x*b))
        end if
        pdf1= dlog(dn)

    END FUNCTION gammadensity11! returns log value

!#################### UNIFORM density ######################

    function uniformden11(val, a, b) result(pdf)
    implicit none
    double precision :: val, a, b, pdf1, pdf

        if (val .ge. b) then
            pdf1 = 0.0d0
        else if (val .le. a) then
            pdf1 = 0.0d0
        else
            pdf1 = (1.0d0/(b-a))
        end if

        pdf = dlog(pdf1)

    end function uniformden11

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Generating random variables for diffierent distributions !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!####################  NORMAL distribution ######################

    FUNCTION randnormal11(mean, stdev) RESULT(c)

    implicit none

    double precision  :: mean, stdev, c, temp(2), r, theta
    double precision, PARAMETER  :: PI=3.141592653589793238462d0

        CALL RANDOM_NUMBER(temp)
        r=(-2.0d0*log(temp(1)))**0.5d0
        theta = 2.0d0*PI*temp(2)
        c= mean+stdev*r*sin(theta)

    END FUNCTION randnormal11

!#################### GAMMA distribution ######################

    RECURSIVE FUNCTION randgamma11(shape, SCALE) RESULT(ans)
    double precision SHAPE, scale, u, w, d, c, x, xsq, g, ans, v
!
! ## Implementation based on "A Simple Method for Generating Gamma Variables"
! ## by George Marsaglia and Wai Wan Tsang.
! ## ACM Transactions on Mathematical Software
! ## Vol 26, No 3, September 2000, pages 363-372.
!
        IF (shape >= 1.0d0) THEN
            d = SHAPE - (1.0d0/3.0d0)
            c = 1.0d0/((9.0d0 * d)**0.5)
            DO while (.true.)
                x = randnormal11(0.0d0, 1.0d0)
                v = 1.0 + c*x
                DO while (v <= 0.0d0)
                    x = randnormal11(0.0d0, 1.0d0)
                    v = 1.0d0 + c*x
                END DO

                v = v*v*v
                CALL RANDOM_NUMBER(u)
                xsq = x*x
                IF ((u < 1.0d0 -.0331d0*xsq*xsq) .OR.  &
                (log(u) < 0.5d0*xsq + d*(1.0d0 - v + log(v))) ) then
                    ans=scale*d*v
                    RETURN
                END IF

            END DO
        ELSE
            g = randgamma11(shape+1.0d0, 1.0d0)
            CALL RANDOM_NUMBER(w)
            ans=scale*g*(w**(1.0d0/shape))
            RETURN
        END IF

    END FUNCTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! To reset the seed for the random number generation method  !! 
!!                     random_number()                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine initrandomseedsir11(temp)
    implicit none
    integer :: n!, clock
    integer, intent(in):: temp
    integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))
        seed = temp
        call random_seed(PUT = seed)
        deallocate(seed)

    end subroutine initrandomseedsir11


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!       To generate random samples without replacement       !! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE ransamsir11(x, a, n, k)
    implicit none
    integer :: j, m, l
    double precision :: u
    integer, intent(in) :: n, k
    integer, intent(in), dimension(n) :: x
    integer, intent(out), dimension(k) :: a

        m=0
        do j = 1, n
            call random_number(u)

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

    end SUBROUTINE ransamsir11


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! To sort an array into ascending order w.r.t specific column  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine  sort(x,xx,n,ppp)
    implicit none
    integer, intent(in) :: n, ppp
    double precision, dimension(n), intent(in) :: x
    double precision, dimension(n, ppp), intent(inout) :: xx
    integer :: i, j, location
    double precision,dimension(ppp):: TT
    double precision :: minimum

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
    end subroutine  sort

end module sir11
