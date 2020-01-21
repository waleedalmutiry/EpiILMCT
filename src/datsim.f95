!######################################################################
!# MODULE: datsim
!# AUTHORS:
!#         Waleed Almutiry <walmutir@uoguelph.ca>,
!#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca> and
!#         Rob Deardon <robert.deardon@ucalgary.ca>
!#
!# DESCRIPTION:
!#
!#     To simulate epidemic from the SIR continuous-time ILMs:
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
!#           datasimulation ............... subroutine
!#           rate ......................... subroutine
!#           randnormal2 .................. function
!#           randgamma2 ................... function
!#           initrandomseed2 .............. subroutine
!#
!######################################################################

module datsim
    use, intrinsic :: iso_c_binding
    implicit none
    public :: datasimulation_f

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 			EPIDEMIC SIMULATION subroutine		 	 	 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine datasimulation_f(n, anum, num, observednum, observedepi, tmax, suspar, nsuspar, powersus, &
    & transpar, ntranspar, powertrans, kernelpar, spark, delta1, delta2, &
    & suscov, transcov, cc, d3, epidat)  bind(C,  name="datasimulation_f_")

    external infinity_value
    external seedin
    external seedout
    external randomnumber

    integer (C_INT), intent(in) :: n, nsuspar, ntranspar, num, anum, observednum ! integers
    real (C_DOUBLE), intent(in), dimension(n, nsuspar) :: suscov             ! susceptibility covariates
    real (C_DOUBLE), intent(in), dimension(n, ntranspar) :: transcov         ! transmissibility covariates
    real (C_DOUBLE), intent(in), dimension(n, n) :: cc, d3                     ! network & distance matrices
    real (C_DOUBLE), intent(in), dimension(nsuspar) :: suspar, powersus       ! susceptibility parameters
    real (C_DOUBLE), intent(in), dimension(ntranspar) :: transpar, powertrans ! transmissibility parameters
    real (C_DOUBLE), intent(in) :: spark, tmax                               ! spark & maximum infection time
    real (C_DOUBLE), intent(in) :: delta1, delta2           ! Parameters of the infectious period distribution
    real (C_DOUBLE), intent(in), dimension(2) :: kernelpar  ! parameter of the kernel function
    real (C_DOUBLE), intent(in), dimension(observednum, 4) :: observedepi     ! observed epidemic to start
    real (C_DOUBLE), dimension(observednum, 4) :: observedepi1     ! observed epidemic to start
    real (C_DOUBLE), dimension(n, 4), intent(out) :: epidat                   ! OUTPUT
    integer (C_INT) :: nnn1                                            ! # of infected by the end of epidemic
    integer (C_INT), dimension(n, 2) :: xx                                   ! Auxiliary variable
    real (C_DOUBLE), dimension(1, 2) :: ts                                   ! Output from the rate subroutine
    real (C_DOUBLE) :: t0                                  ! current infection time during the simulation
    real (C_DOUBLE) :: Inf, u                                                 !defining Infinity
    integer (C_INT) :: ctr, i, j, sdg, mg
    integer (C_INT), allocatable, dimension(:) :: mmg


! defining Infinity
       call infinity_value(Inf)
       call seedin()

       SELECT CASE (anum)

        CASE (1)

! case (1): for contact network or distanse based with spark term.  (SIR)
!			also for distanse based with/without spark term.

! defining auxiliary variable (xx) for the status of each individual:
! 0 : susceptible
! 1 : infectious
! 2 : removed

          xx          = 0
          xx(:, 1)     = (/(j, j=1, n)/)

          epidat      = 0.0_c_double
          observedepi1 = observedepi
! initial observed epidemic:
        if (observednum .eq. 1) then
            if (observedepi1(1,1) .eq. 0) then
                call randomnumber(u)
                observedepi1(1,1) = int(u*n) + 1
            end if
        end if

        do j = 1,  observednum
            if (observedepi1(j, 2) .eq. 0.0_c_double) then
               epidat(j, 1)  = observedepi1(j, 1)
               epidat(j, 4)  = observedepi1(j, 4)
               epidat(j, 3)  = randgamma2(delta1, 1.0_c_double/delta2)
               epidat(j, 2)  = epidat(j, 4) + epidat(j, 3)
            else
               epidat(j, :)  = observedepi1(j, :)
            end if
          end do

! the current infection time to start from:
          t0 = epidat(observednum, 4)
          xx(int(epidat(1:observednum, 1)), 2) = 1

          mg  = size(pack(int(epidat(1:(observednum-1), 1)), epidat(1:(observednum-1), 2) .lt. &
                    & epidat(observednum, 4) ))
          if (mg .gt. 0) then
            allocate(mmg(mg))
            mmg = pack(int(epidat(1:(observednum-1), 1)), epidat(1:(observednum-1), 2) .lt. &
                        & epidat(observednum, 4) )
            xx(mmg, 2) = 2
            deallocate(mmg)
          end if

! starting simulating the epidemic
          ctr = observednum

            do while( (ctr .le. n) )
                ctr = ctr + 1
! to provide the next infected individual with minmum waiting time to infection:
                call rate(n, num, suspar, nsuspar, powersus, transpar, ntranspar, powertrans,  &
                    & kernelpar, spark, xx, suscov, transcov, cc, d3, ts)

! to judge stopping the epidemic or keep generating:
                if ( (ts(1, 2) .ne. Inf) .and. (ts(1, 1) .ne. 0.0_c_double) ) then
                    ts = ts
                else
                    where(xx(:, 2) .eq. 1) xx(:, 2) = 2
                    exit
                end if

!making sure there is still infectious individual that can transmit the disease

                sdg = 0
                do i = 1,  (ctr-1)
                    if ( (epidat(i, 2) .gt. (ts(1, 2)+t0)) ) then
                        sdg = sdg +1
                    else
                        sdg = sdg
                    end if
                end do

! assigning infection time,  infectious period and removal time for the newly infected:

                if (sdg .eq. 0 ) then
                    where(xx(:, 2) .eq. 1) xx(:, 2) = 2
                    exit
                else
                    epidat(ctr, 3) = randgamma2(delta1, 1.0_c_double/delta2)
                    epidat(ctr, 4) = ts(1, 2) + t0
                    epidat(ctr, 2) = epidat(ctr, 3) + epidat(ctr, 4)
                    epidat(ctr, 1) = ts(1, 1)
                    t0 = epidat(ctr, 4)
                    xx(int(epidat(ctr, 1)), 2) = 1
                end if

                if ( (epidat(ctr, 4) .gt. tmax) ) then
                    epidat(ctr, 3) = 0.0_c_double
                    epidat(ctr, 4) = Inf
                    epidat(ctr, 2) = Inf
                    exit
                end if

! update the auxiliary variable of the status of individuals:

                mg  = size(pack(int(epidat(1:ctr-1, 1)), epidat(1:ctr-1, 2) .lt. epidat(ctr, 4) ))
                if (mg .gt. 0) then
                    allocate(mmg(mg))
                        mmg = pack(int(epidat(1:ctr-1, 1)), epidat(1:ctr-1, 2) .lt. epidat(ctr, 4) )
                        xx(mmg, 2) = 2
                    deallocate(mmg)
                end if

            end do !end of while loop.

! assigning infinity values for those uninfected by the end of the epidemic

            nnn1 = count(epidat(:, 2)  .ne. 0.0_c_double)
            do i = (nnn1+1),  n
                do j = 1, n
                    if (all(int(epidat(1:(i-1), 1)) .ne. j)) then
                        epidat(i, 1) = dble(j)
                    end if
                end do
                epidat(i, 2) = Inf
                epidat(i, 4) = Inf
            end do


        CASE (2)

! case (2): for contact without spark term.  (SIR)

! defining auxiliary variable (xx) for the status of each individual:
! 0 : susceptible
! 1 : infectious
! 2 : removed

            xx          = 0
            xx(:, 1)     = (/(j, j=1, n)/)
            epidat      = 0.0_c_double
            observedepi1 = observedepi

! initial observed epidemic:
            if (observednum .eq. 1) then
                if (observedepi(1,1) .eq. 0) then
                    call randomnumber(u)
                    observedepi1(1,1) = int(u*n) + 1
                end if
            end if

            do j = 1,  observednum
                if (observedepi1(j, 2) .eq. 0.0_c_double) then
                    epidat(j, 1)  = observedepi1(j, 1)
                    epidat(j, 4)  = observedepi1(j, 4)
                    epidat(j, 3)  = randgamma2(delta1, 1.0_c_double/delta2)
                    epidat(j, 2)  = epidat(observednum, 4) + epidat(observednum, 3)
                else
                    epidat(j, :)  = observedepi1(j, :)
                end if
            end do

! the current infection time to start from:
            t0 = epidat(observednum, 4)
            xx(int(epidat(1:observednum, 1)), 2) = 1

            mg  = size(pack(int(epidat(1:observednum-1, 1)), epidat(1:observednum-1, 2) .lt. &
                        & epidat(observednum, 4) ))
            if (mg .gt. 0) then
                allocate(mmg(mg))
                    mmg = pack(int(epidat(1:observednum-1, 1)), epidat(1:observednum-1, 2) .lt. &
                        & epidat(observednum, 4) )
                    xx(mmg, 2) = 2
                deallocate(mmg)
            end if

! starting simulating the epidemic

            ctr = observednum

            do while( (ctr .le. n) )

                ctr = ctr + 1

! to provide the next infected individual with minmum waiting time to infection:

                call rate(n, num, suspar, nsuspar, powersus, transpar, ntranspar, powertrans, &
                    &kernelpar, spark, xx, suscov, transcov, cc, d3, ts)

! to judge stopping the epidemic or keep generating:

                if ( (ts(1, 2) .ne. Inf) .and. (ts(1, 1) .ne. 0.0_c_double) ) then
                    ts = ts
                else
                    where(xx(:, 2) .eq. 1) xx(:, 2) = 2
                    exit
                end if

!making sure there is still infectious individual(s) that can transmit the disease

                sdg = 0
                do i = 1,  (ctr-1)
                    if ( (epidat(i, 2) .gt. (ts(1, 2)+t0)) .and. (cc(int(epidat(i, 1)), int(ts(1, 1))) .gt. 0.0_c_double) ) then
                        sdg = sdg +1
                    else
                        sdg = sdg
                    end if
                end do

! assigning infection time,  infectious period and removal time for the newly infected:

                if ( (sdg .eq. 0 ) ) then
                    where(xx(:, 2) .eq. 1) xx(:, 2) = 2
                    exit
                else
                    epidat(ctr, 3) = randgamma2(delta1, 1.0_c_double/delta2)
                    epidat(ctr, 4) = ts(1, 2) + t0
                    epidat(ctr, 2) = epidat(ctr, 3) + epidat(ctr, 4)
                    epidat(ctr, 1) = ts(1, 1)
                    t0 = epidat(ctr, 4)
                    xx(int(epidat(ctr, 1)), 2) = 1
                end if

                if ( (epidat(ctr, 4) .gt. tmax) ) then
                    epidat(ctr, 3) = 0.0_c_double
                    epidat(ctr, 4) = Inf
                    epidat(ctr, 2) = Inf
                    exit
                end if

! update the auxiliary variable of the status of individuals:

                mg  = size(pack(int(epidat(1:ctr-1, 1)), epidat(1:ctr-1, 2) .lt. epidat(ctr, 4) ))
                if (mg .gt. 0) then
                    allocate(mmg(mg))
                        mmg = pack(int(epidat(1:ctr-1, 1)), epidat(1:ctr-1, 2) .lt. epidat(ctr, 4) )
                        xx(mmg, 2) = 2
                    deallocate(mmg)
                end if

            end do !end of while loop.

! assigning infinity values for those uninfected by the end of the epidemic

            nnn1 = count(epidat(:, 2)  .ne. 0.0_c_double)
            do i = (nnn1+1),  n
                do j = 1, n
                    if (all(int(epidat(1:(i-1), 1)) .ne. j)) then
                        epidat(i, 1) = dble(j)
                    end if
                end do
                epidat(i, 2) = Inf
                epidat(i, 4) = Inf
            end do

        END SELECT
        call seedout()
    end subroutine datasimulation_f


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 			INFECTIVITY RATE subroutine 		 	 	 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine rate(n, num, suspar, nsuspar, powersus, transpar, ntranspar, powertrans, kernelpar, spark, xx, &
    & suscov, transcov, cc, d3, mms)

    external infinity_value

    integer :: i, mg, mg1, j, m

    integer, intent(in) :: n, nsuspar, ntranspar, num                    !integers
    real (C_DOUBLE), intent(in), dimension(n, nsuspar) :: suscov      ! susceptibility covariates
    real (C_DOUBLE), intent(in), dimension(n, ntranspar) :: transcov  ! transmissibility covariates
    real (C_DOUBLE), intent(in), dimension(n, n) :: d3, cc              ! network and distance matrices
    real (C_DOUBLE), intent(in), dimension(nsuspar) :: suspar, powersus! susceptibility parameters
    real (C_DOUBLE), intent(in), dimension(ntranspar) :: transpar, powertrans ! transmissibility parameters
    real (C_DOUBLE), intent(in) :: spark                             ! spark parameter
    real (C_DOUBLE), intent(in), dimension(2) :: kernelpar            ! parameters of the kernel function
    integer, intent(in), dimension(n, 2) :: xx                 ! auxiliary variable see above!
    real (C_DOUBLE), dimension(nsuspar) :: suscov1              ! suseptible covar. for one individual
    real (C_DOUBLE), dimension(ntranspar) :: transcov1          ! transmissibilty covar. for one individual
    real (C_DOUBLE) :: Inf
    integer, allocatable, dimension(:) :: mmg, mmg1
    real (C_DOUBLE), allocatable, dimension(:, :) :: rr
    real (C_DOUBLE), allocatable, dimension(:) :: hg
    real (C_DOUBLE), dimension(1, 2) :: mms

    call infinity_value(Inf)

    SELECT CASE (num)

    CASE (1)
! Calculating the infectivity rate of contact network-based ILM with spark term.

! defining infectious individuals:
        mg  = size( pack(xx(:, 1), xx(:, 2).eq. 1.0_c_double ) )
        allocate(mmg(mg))
        mmg = pack(xx(:, 1), xx(:, 2).eq. 1.0_c_double )

! defining susceptible individuals:

        mg1 = size( pack(xx(:, 1), xx(:, 2).eq. 0.0_c_double ) )
        allocate(mmg1(mg1))
        mmg1 = pack(xx(:, 1), xx(:, 2).eq. 0.0_c_double )


! declaring a variable with size equal to the number of susceptible individuals
! the first column is the id number of the susceptible individuals

        allocate(rr(mg1, 2))
        rr(:, 1) = mmg1

! start calculating the infectivity rate for each susceptible individuals:

        do i = 1, mg1

            allocate(hg(mg))

            do j = 1,  mg

                do m = 1,  ntranspar
                    transcov1(m) = transcov(mmg(j), m)**powertrans(m)
                end do

                hg(j) = dot_product(transpar, transcov1)*cc(mmg(j), mmg1(i))

            end do

            do m = 1,  nsuspar
                suscov1(m) = suscov(mmg1(i), m)**powersus(m)
            end do

            rr(i, 2) = (dot_product(suspar, suscov1)*sum(hg))+ spark

            deallocate(hg)

        end do

! assigning waiting time to infection for each susceptible individual:

        do i = 1, mg1
            if (rr(i, 2) .eq. 0.0_c_double) then
                rr(i, 2) = Inf
            else
                rr(i, 2) = randgamma2(1.0_c_double, 1.0_c_double/rr(i, 2))
            end if
        end do

! choose the one with minmum waiting time as the newly infected individual:
        if (all(rr(:, 2) .eq. Inf).eqv. .true.) then
            mms(1, 1)  = 0.0_c_double
            mms(1, 2)  = Inf
        else
            mms(1, 1:2)  = rr(int(minloc(rr(:, 2), 1)), :)
        end if

        deallocate(rr)
        deallocate(mmg1)
        deallocate(mmg)

    CASE (2)
! Calculating the infectivity rate of distance-based ILM with spark term and distance kernel = powerlaw.

! defining infectious individuals:
        mg  = size( pack(xx(:, 1), xx(:, 2).eq. 1.0_c_double ) )
        allocate(mmg(mg))
        mmg = pack(xx(:, 1), xx(:, 2).eq. 1.0_c_double )

! defining susceptible individuals:
        mg1 = size( pack(xx(:, 1), xx(:, 2).eq. 0.0_c_double ) )
        allocate(mmg1(mg1))
        mmg1 = pack(xx(:, 1), xx(:, 2).eq. 0.0_c_double )


! declaring a variable with size equal to the number of susceptible individuals
! the first column is the id number of the susceptible individuals
        allocate(rr(mg1, 2))
        rr(:, 1) = mmg1

! start calculating the infectivity rate for each susceptible individuals:
        do i = 1, mg1
            allocate(hg(mg))
            do j = 1,  mg

                do m = 1,  ntranspar
                transcov1(m) = transcov(mmg(j), m)**powertrans(m)
                end do

                hg(j) = dot_product(transpar, transcov1)*&
                & (d3(mmg(j), mmg1(i))**(-kernelpar(1)))
            end do

            do m = 1,  nsuspar
                suscov1(m) = suscov(mmg1(i), m)**powersus(m)
            end do

            rr(i, 2) = (dot_product(suspar, suscov1)*sum(hg))+spark
            deallocate(hg)
        end do

! assigning waiting time to infection for each susceptible individual:
        do i = 1, mg1
            if (rr(i, 2) .eq. 0.0_c_double) then
                rr(i, 2) = Inf
            else
                rr(i, 2) = randgamma2(1.0_c_double, 1.0_c_double/rr(i, 2))
            end if
        end do

! choose the one with minmum waiting time as the newly infected individual:
        if (all(rr(:, 2) .eq. Inf).eqv. .true.) then
            mms(1, 1)  = 0.0_c_double
            mms(1, 2)  = Inf
        else
            mms(1, 1:2)  = rr(int(minloc(rr(:, 2), 1)), :)
        end if

        deallocate(rr)
        deallocate(mmg1)
        deallocate(mmg)

    CASE (3)
! Calculating the infectivity rate of distance-based ILM with spark term and distance kernel = Cauchy.

! defining infectious individuals:
        mg  = size( pack(xx(:, 1), xx(:, 2).eq. 1.0_c_double ) )
        allocate(mmg(mg))
        mmg = pack(xx(:, 1), xx(:, 2).eq. 1.0_c_double )

! defining susceptible individuals:
        mg1 = size( pack(xx(:, 1), xx(:, 2).eq. 0.0_c_double ) )
        allocate(mmg1(mg1))
        mmg1 = pack(xx(:, 1), xx(:, 2).eq. 0.0_c_double )


! declaring a variable with size equal to the number of susceptible individuals
! the first column is the id number of the susceptible individuals
        allocate(rr(mg1, 2))
        rr(:, 1) = mmg1

! start calculating the infectivity rate for each susceptible individuals:
        do i = 1, mg1
            allocate(hg(mg))
            do j = 1,  mg
                do m = 1,  ntranspar
                    transcov1(m) = transcov(mmg(j), m)**powertrans(m)
                end do
                hg(j) = dot_product(transpar, transcov1)*&
                & (kernelpar(1)/((d3(mmg(j), mmg1(i))**(2.0_c_double)) +&
                & (kernelpar(1)**(2.0_c_double))))
            end do

            do m = 1,  nsuspar
                suscov1(m) = suscov(mmg1(i), m)**powersus(m)
            end do

            rr(i, 2) = (dot_product(suspar, suscov1)*sum(hg))+spark

            deallocate(hg)

        end do

! assigning waiting time to infection for each susceptible individual:
        do i = 1, mg1
            if (rr(i, 2) .eq. 0.0_c_double) then
                rr(i, 2) = Inf
            else
                rr(i, 2) = randgamma2(1.0_c_double, 1.0_c_double/rr(i, 2))
            end if
        end do

! choose the one with minmum waiting time as the newly infected individual:
        if (all(rr(:, 2) .eq. Inf).eqv. .true.) then
            mms(1, 1)  = 0.0_c_double
            mms(1, 2)  = Inf
        else
            mms(1, 1:2)  = rr(int(minloc(rr(:, 2), 1)), :)
        end if

        deallocate(rr)
        deallocate(mmg1)
        deallocate(mmg)


    CASE (4)
! Calculating the infectivity rate of including both distance and network-based ILM with spark term and distance kernel = powerlaw.

! defining infectious individuals:
        mg  = size( pack(xx(:, 1), xx(:, 2).eq. 1.0_c_double ) )
        allocate(mmg(mg))
        mmg = pack(xx(:, 1), xx(:, 2).eq. 1.0_c_double )

! defining susceptible individuals:
        mg1 = size( pack(xx(:, 1), xx(:, 2).eq. 0.0_c_double ) )
        allocate(mmg1(mg1))
        mmg1 = pack(xx(:, 1), xx(:, 2).eq. 0.0_c_double )

! declaring a variable with size equal to the number of susceptible individuals
! the first column is the id number of the susceptible individuals
        allocate(rr(mg1, 2))
        rr(:, 1) = mmg1

! start calculating the infectivity rate for each susceptible individuals:
        do i = 1, mg1
            allocate(hg(mg))
            do j = 1,  mg
                do m = 1,  ntranspar
                    transcov1(m) = transcov(mmg(j), m)**powertrans(m)
                end do
                hg(j) = dot_product(transpar, transcov1)*&
                & (kernelpar(2)*cc(mmg(j), mmg1(i)) + (d3(mmg(j), mmg1(i))**(-kernelpar(1))))
            end do

            do m = 1,  nsuspar
            suscov1(m) = suscov(mmg1(i), m)**powersus(m)
            end do

            rr(i, 2) = (dot_product(suspar, suscov1)*sum(hg))+spark
            deallocate(hg)
        end do

! assigning waiting time to infection for each susceptible individual:
        do i = 1, mg1
            if (rr(i, 2) .eq. 0.0_c_double) then
                rr(i, 2) = Inf
            else
                rr(i, 2) = randgamma2(1.0_c_double, 1.0_c_double/rr(i, 2))
            end if
        end do

! choose the one with minmum waiting time as the newly infected individual:
        if (all(rr(:, 2) .eq. Inf).eqv. .true.) then
            mms(1, 1)  = 0.0_c_double
            mms(1, 2)  = Inf
        else
            mms(1, 1:2)  = rr(int(minloc(rr(:, 2), 1)), :)
        end if

        deallocate(rr)
        deallocate(mmg1)
        deallocate(mmg)

    CASE (5)
! Calculating the infectivity rate of including both distance and network-based ILM with spark term and distance kernel = Cauchy.

! defining infectious individuals:
        mg  = size( pack(xx(:, 1), xx(:, 2).eq. 1.0_c_double ) )
        allocate(mmg(mg))
        mmg = pack(xx(:, 1), xx(:, 2).eq. 1.0_c_double )

! defining susceptible individuals:
        mg1 = size( pack(xx(:, 1), xx(:, 2).eq. 0.0_c_double ) )
        allocate(mmg1(mg1))
        mmg1 = pack(xx(:, 1), xx(:, 2).eq. 0.0_c_double )


! declaring a variable with size equal to the number of susceptible individuals
! the first column is the id number of the susceptible individuals
        allocate(rr(mg1, 2))
        rr(:, 1) = mmg1

! start calculating the infectivity rate for each susceptible individuals:
        do i = 1, mg1
            allocate(hg(mg))
            do j = 1,  mg
                do m = 1,  ntranspar
                    transcov1(m) = transcov(mmg(j), m)**powertrans(m)
                end do
                hg(j) = dot_product(transpar, transcov1)*&
                & (kernelpar(2)*cc(mmg(j), mmg1(i)) + &
                &(kernelpar(1)/((d3(mmg(j), mmg1(i))**(2.0_c_double)) + &
                &(kernelpar(1)**(2.0_c_double)))))
            end do

            do m = 1,  nsuspar
                suscov1(m) = suscov(mmg1(i), m)**powersus(m)
            end do

            rr(i, 2) = (dot_product(suspar, suscov1)*sum(hg))+spark
            deallocate(hg)
        end do

! assigning waiting time to infection for each susceptible individual:
        do i = 1, mg1
            if (rr(i, 2) .eq. 0.0_c_double) then
                rr(i, 2) = Inf
            else
                rr(i, 2) = randgamma2(1.0_c_double, 1.0_c_double/rr(i, 2))
            end if
        end do

! choose the one with minmum waiting time as the newly infected individual:
        if (all(rr(:, 2) .eq. Inf).eqv. .true.) then
            mms(1, 1)  = 0.0_c_double
            mms(1, 2)  = Inf
        else
            mms(1, 1:2)  = rr(int(minloc(rr(:, 2), 1)), :)
        end if

        deallocate(rr)
        deallocate(mmg1)
        deallocate(mmg)


    END SELECT


    end subroutine rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Generating random variables for different distributions !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!####################  NORMAL distribution ######################

    FUNCTION randnormal2(mean, stdev) RESULT(c)
    implicit none

    real (C_DOUBLE) :: mean, stdev, c, temp(2), r, theta
    real (C_DOUBLE),  PARAMETER :: PI=3.141592653589793238462_c_double

    CALL randomnumber(temp(1))
    CALL randomnumber(temp(2))
        r=(-2.0_c_double*log(temp(1)))**0.5_c_double
        theta = 2.0_c_double*PI*temp(2)
        c= mean+stdev*r*sin(theta)

    END FUNCTION randnormal2

!#################### GAMMA distribution ######################

    RECURSIVE FUNCTION randgamma2(shape,  SCALE) RESULT(ans)
    real (C_DOUBLE) :: SHAPE, scale, u, w, d, c, x, xsq, g, ans, v
! DESCRIPTION: Implementation based on "A Simple Method for Generating Gamma Variables"
! by George Marsaglia and Wai Wan Tsang.
! ACM Transactions on Mathematical Software and released in public domain.
! ## Vol 26,  No 3,  September 2000,  pages 363-372.

        IF (shape >= 1.0_c_double) THEN
            d = SHAPE - (1.0_c_double/3.0_c_double)
            c = 1.0_c_double/((9.0_c_double * d)**0.5_c_double)
            DO while (.true.)
                x = randnormal2(0.0_c_double,  1.0_c_double)
                v = 1.0 + c*x
                DO while (v <= 0.0_c_double)
                    x = randnormal2(0.0_c_double,  1.0_c_double)
                    v = 1.0_c_double + c*x
                END DO

            v = v*v*v
            CALL randomnumber(u)
            xsq = x*x
            IF ((u < 1.0_c_double -.0331_c_double*xsq*xsq) .OR.  &
            (log(u) < 0.5_c_double*xsq + d*(1.0_c_double - v + log(v))) ) then
                ans=scale*d*v
                RETURN
            END IF

            END DO
        ELSE
            g = randgamma2(shape+1.0_c_double,  1.0_c_double)
            CALL randomnumber(w)
            ans=scale*g*(w**(1.0_c_double/shape))
            RETURN
        END IF

    END FUNCTION randgamma2

end module datsim
