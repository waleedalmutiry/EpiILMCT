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
!#           datasimulation ............... subroutine
!#           rate ......................... subroutine
!#           randnormal2 .................. function
!#           randgamma2 ................... function
!#           initrandomseed2 .............. subroutine
!#
!######################################################################

module datsim
    use ISO_C_BINDING
    implicit none
    public:: datasimulation

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 			EPIDEMIC SIMULATION subroutine		 	 	 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine datasimulation(n,anum,num,observednum,observedepi,tmax,suspar,nsuspar,powersus,&
    & transpar,ntranspar,powertrans,kernelpar,spark,delta1,delta2,&
    & suscov,transcov,cc,d3,epidat)  bind(C, name="datasimulation_")

    external infinity_value

    integer (C_INT),intent(in)::n,nsuspar,ntranspar,num,anum,observednum ! integers
    real (C_DOUBLE),intent(in),dimension(n,nsuspar):: suscov             ! susceptibility covariates
    real (C_DOUBLE),intent(in),dimension(n,ntranspar):: transcov         ! transmissibility covariates
    real (C_DOUBLE),intent(in),dimension(n,n)::cc,d3                     ! network & distance matrices
    real (C_DOUBLE),intent(in),dimension(nsuspar)::suspar,powersus       ! susceptibility parameters
    real (C_DOUBLE),intent(in),dimension(ntranspar)::transpar,powertrans ! transmissibility parameters
    real (C_DOUBLE),intent(in)::spark,tmax                               ! spark & maximum infection time
    real (C_DOUBLE),intent(in)::delta1,delta2           ! Parameters of the infectious period distribution
    real (C_DOUBLE),intent(in),dimension(2)::kernelpar  ! parameter of the kernel function
    real (C_DOUBLE),intent(in),dimension(observednum,4)::observedepi     ! observed epidemic to start
    real (C_DOUBLE),dimension(n,4),intent(out)::epidat                   ! OUTPUT
    integer (C_INT)::nnn1                                            ! # of infected by the end of epidemic
    integer (C_INT),dimension(n,2)::xx                                   ! Auxiliary variable
    real (C_DOUBLE),dimension(1,2)::ts                                   ! Output from the rate subroutine
    real (C_DOUBLE)::t0                                  ! current infection time during the simulation
    real (C_DOUBLE)::Inf                                                 !defining Infinity
    integer (C_INT)::ctr,i,j,sdg,mg
    integer (C_INT),allocatable,dimension(:)::mmg

       call initrandomseed2()

! defining Infinity
       call infinity_value(Inf)

       SELECT CASE (anum)

        CASE (1)

! case (1): for contact network or distanse based with spark term.  (SIR)
!			also for distanse based with/without spark term.

! defining auxiliary variable (xx) for the status of each individual:
! 0 : susceptible
! 1 : infectious
! 2 : removed

          xx          = 0
          xx(:,1)     = (/(j,j=1,n)/)

          epidat      = 0.0_c_double

! initial observed epidemic:

          do j = 1 , observednum
            if(observedepi(j,2) .eq. 0.0_c_double)then
               epidat(j,1)  = observedepi(j,1)
               epidat(j,4)  = observedepi(j,4)
               epidat(j,3)  = randgamma2(delta1,1.0_c_double/delta2)
               epidat(j,2)  = epidat(observednum,4) + epidat(observednum,3)
            else
               epidat(j,:)  = observedepi(j,:)
            end if
          end do

! the current infection time to start from:
          t0 = epidat(observednum,4)
          xx(int(epidat(1:observednum,1)),2) = 1

          mg  = size(pack(int(epidat(1:(observednum-1),1)),epidat(1:(observednum-1),2) .lt. &
                    & epidat(observednum,4) ))
          if(mg .gt. 0)then
            allocate(mmg(mg))
            mmg = pack(int(epidat(1:(observednum-1),1)),epidat(1:(observednum-1),2) .lt. &
                        & epidat(observednum,4) )
            xx(mmg,2) = 2
            deallocate(mmg)
          end if

! starting simulating the epidemic
          ctr = observednum

            do while( (ctr .le. n) )
                ctr = ctr + 1
! to provide the next infected individual with minmum waiting time to infection:
                call rate(n,num,suspar,nsuspar,powersus,transpar,ntranspar,powertrans, &
                    & kernelpar,spark,xx,suscov,transcov,cc,d3,ts)

! to judge stopping the epidemic or keep generating:
                if( (ts(1,2) .ne. Inf) .and. (ts(1,1) .ne. 0.0_c_double) ) then
                    ts = ts
                else
                    where(xx(:,2) .eq. 1) xx(:,2) = 2
                    exit
                end if

!making sure there is still infectious individual that can transmit the disease

                sdg = 0
                do i = 1 , (ctr-1)
                    if( (epidat(i,2) .gt. (ts(1,2)+t0)) ) then
                        sdg = sdg +1
                    else
                        sdg = sdg
                    end if
                end do

! assigning infection time, infectious period and removal time for the newly infected:

                if(sdg .eq. 0 ) then
                    where(xx(:,2) .eq. 1) xx(:,2) = 2
                    exit
                else
                    epidat(ctr,3) = randgamma2(delta1,1.0_c_double/delta2)
                    epidat(ctr,4) = ts(1,2) + t0
                    epidat(ctr,2) = epidat(ctr,3) + epidat(ctr,4)
                    epidat(ctr,1) = ts(1,1)
                    t0 = epidat(ctr,4)
                    xx(int(epidat(ctr,1)),2) = 1
                end if

                if( (epidat(ctr,4) .gt. tmax) ) then
                    epidat(ctr,3) = 0.0_c_double
                    epidat(ctr,4) = 0.0_c_double
                    epidat(ctr,2) = 0.0_c_double
                    exit
                end if

! update the auxiliary variable of the status of individuals:

                mg  = size(pack(int(epidat(1:ctr-1,1)),epidat(1:ctr-1,2) .lt. epidat(ctr,4) ))
                if(mg .gt. 0)then
                    allocate(mmg(mg))
                        mmg = pack(int(epidat(1:ctr-1,1)),epidat(1:ctr-1,2) .lt. epidat(ctr,4) )
                        xx(mmg,2) = 2
                    deallocate(mmg)
                end if

            end do !end of while loop.

! assigning infinity values for those uninfected by the end of the epidemic

            nnn1 = count(epidat(:,2)  .ne. 0.0_c_double)
            do i = (nnn1+1) , n
                do j = 1,n
                    if(all(int(epidat(1:(i-1),1)) .ne. j))then
                        epidat(i,1) = dble(j)
                    end if
                end do
                epidat(i,2) = Inf
                epidat(i,4) = Inf
            end do


        CASE (2)

! case (2): for contact without spark term.  (SIR)

! defining auxiliary variable (xx) for the status of each individual:
! 0 : susceptible
! 1 : infectious
! 2 : removed

            xx          = 0
            xx(:,1)     = (/(j,j=1,n)/)
            epidat      = 0.0_c_double

! initial observed epidemic:

            do j = 1 , observednum
                if(observedepi(j,2) .eq. 0.0_c_double)then
                    epidat(j,1)  = observedepi(j,1)
                    epidat(j,4)  = observedepi(j,4) 
                    epidat(j,3)  = randgamma2(delta1,1.0_c_double/delta2)
                    epidat(j,2)  = epidat(observednum,4) + epidat(observednum,3)
                else
                    epidat(j,:)  = observedepi(j,:)
                end if
            end do

! the current infection time to start from:
            t0 = epidat(observednum,4)
            xx(int(epidat(1:observednum,1)),2) = 1

            mg  = size(pack(int(epidat(1:observednum-1,1)),epidat(1:observednum-1,2) .lt. &
                        & epidat(observednum,4) ))
            if(mg .gt. 0)then
                allocate(mmg(mg))
                    mmg = pack(int(epidat(1:observednum-1,1)),epidat(1:observednum-1,2) .lt. &
                        & epidat(observednum,4) )
                    xx(mmg,2) = 2
                deallocate(mmg)
            end if

! starting simulating the epidemic

            ctr = observednum

            do while( (ctr .le. n) )

                ctr = ctr + 1

! to provide the next infected individual with minmum waiting time to infection:

                call rate(n,num,suspar,nsuspar,powersus,transpar,ntranspar,powertrans,&
                    &kernelpar,spark,xx,suscov,transcov,cc,d3,ts)

! to judge stopping the epidemic or keep generating:

                if( (ts(1,2) .ne. Inf) .and. (ts(1,1) .ne. 0.0_c_double) ) then
                    ts = ts
                else
                    where(xx(:,2) .eq. 1) xx(:,2) = 2
                    exit
                end if

!making sure there is still infectious individual(s) that can transmit the disease

                sdg = 0
                do i = 1 , (ctr-1)
                    if( (epidat(i,2) .gt. (ts(1,2)+t0)) .and. (cc(int(epidat(i,1)),int(ts(1,1))) .gt. 0.0_c_double) ) then
                        sdg = sdg +1
                    else
                        sdg = sdg
                    end if
                end do

! assigning infection time, infectious period and removal time for the newly infected:

                if( (sdg .eq. 0 ) ) then
                    where(xx(:,2) .eq. 1) xx(:,2) = 2
                    exit
                else
                    epidat(ctr,3) = randgamma2(delta1,1.0_c_double/delta2)
                    epidat(ctr,4) = ts(1,2) + t0
                    epidat(ctr,2) = epidat(ctr,3) + epidat(ctr,4)
                    epidat(ctr,1) = ts(1,1)
                    t0 = epidat(ctr,4)
                    xx(int(epidat(ctr,1)),2) = 1
                end if

                if( (epidat(ctr,4) .gt. tmax) ) then
                    epidat(ctr,3) = 0.0_c_double
                    epidat(ctr,4) = 0.0_c_double
                    epidat(ctr,2) = 0.0_c_double
                    exit
                end if

! update the auxiliary variable of the status of individuals:

                mg  = size(pack(int(epidat(1:ctr-1,1)),epidat(1:ctr-1,2) .lt. epidat(ctr,4) ))
                if(mg .gt. 0)then
                    allocate(mmg(mg))
                        mmg = pack(int(epidat(1:ctr-1,1)),epidat(1:ctr-1,2) .lt. epidat(ctr,4) )
                        xx(mmg,2) = 2
                    deallocate(mmg)
                end if

            end do !end of while loop.

! assigning infinity values for those uninfected by the end of the epidemic

            nnn1 = count(epidat(:,2)  .ne. 0.0_c_double)
            do i = (nnn1+1) , n
                do j = 1,n
                    if(all(int(epidat(1:(i-1),1)) .ne. j))then
                        epidat(i,1) = dble(j)
                    end if
                end do
                epidat(i,2) = Inf
                epidat(i,4) = Inf
            end do

        END SELECT

    end subroutine datasimulation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 			INFECTIVITY RATE subroutine 		 	 	 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine rate(n,num,suspar,nsuspar,powersus,transpar,ntranspar,powertrans,kernelpar,spark,xx,&
    & suscov,transcov,cc,d3,mms)

    external infinity_value

    integer::i,mg,mg1,j,m

    integer,intent(in)::n,nsuspar,ntranspar,num                    !integers
    double precision,intent(in),dimension(n,nsuspar):: suscov      ! susceptibility covariates
    double precision,intent(in),dimension(n,ntranspar):: transcov  ! transmissibility covariates
    double precision,intent(in),dimension(n,n)::d3,cc              ! network and distance matrices
    double precision,intent(in),dimension(nsuspar)::suspar,powersus! susceptibility parameters
    double precision,intent(in),dimension(ntranspar)::transpar,powertrans ! transmissibility parameters
    double precision,intent(in)::spark                             ! spark parameter
    double precision,intent(in),dimension(2)::kernelpar            ! parameters of the kernel function
    integer,intent(in),dimension(n,2)::xx                 ! auxiliary variable see above!
    double precision,dimension(nsuspar)::suscov1              ! suseptible covar. for one individual
    double precision,dimension(ntranspar)::transcov1          ! transmissibilty covar. for one individual
    double precision::Inf
    integer,allocatable,dimension(:)::mmg,mmg1
    double precision,allocatable,dimension(:,:):: rr
    double precision,allocatable,dimension(:):: hg
    double precision,dimension(1,2)::mms

    call infinity_value(Inf)

    SELECT CASE (num)

    CASE (1)
! Calculating the infectivity rate of contact network-based ILM with spark term.

! defining infectious individuals:
        mg  = size( pack(xx(:,1),xx(:,2).eq. 1.0d0 ) )
        allocate(mmg(mg))
        mmg = pack(xx(:,1),xx(:,2).eq. 1.0d0 )

! defining susceptible individuals:

        mg1 = size( pack(xx(:,1),xx(:,2).eq. 0.0d0 ) )
        allocate(mmg1(mg1))
        mmg1 = pack(xx(:,1),xx(:,2).eq. 0.0d0 )


! declaring a variable with size equal to the number of susceptible individuals
! the first column is the id number of the susceptible individuals

        allocate(rr(mg1,2))
        rr(:,1) = mmg1

! start calculating the infectivity rate for each susceptible individuals:

        do i = 1,mg1

            allocate(hg(mg))

            do j = 1 , mg

                do m = 1 , ntranspar
                    transcov1(m) = transcov(mmg(j),m)**powertrans(m)
                end do

                hg(j) = dot_product(transpar,transcov1)*cc(mmg(j),mmg1(i))

            end do

            do m = 1 , nsuspar
                suscov1(m) = suscov(mmg1(i),m)**powersus(m)
            end do

            rr(i,2) = (dot_product(suspar,suscov1)*sum(hg))+ spark

            deallocate(hg)

        end do

! assigning waiting time to infection for each susceptible individual:

        do i = 1,mg1
            if(rr(i,2) .eq. 0.0d0)then
                rr(i,2) = Inf
            else
                rr(i,2) = randgamma2(1.0d0,1.0d0/rr(i,2))
            end if
        end do

! choose the one with minmum waiting time as the newly infected individual:
        if(all(rr(:,2) .eq. Inf).eqv. .true.) then
            mms(1,1)  = 0.0d0
            mms(1,2)  = Inf
        else
            mms(1,1:2)  = rr(int(minloc(rr(:,2),1)),:)
        end if

        deallocate(rr)
        deallocate(mmg1)
        deallocate(mmg)

    CASE (2)
! Calculating the infectivity rate of distance-based ILM with spark term and distance kernel = powerlaw.

! defining infectious individuals:
        mg  = size( pack(xx(:,1),xx(:,2).eq. 1.0d0 ) )
        allocate(mmg(mg))
        mmg = pack(xx(:,1),xx(:,2).eq. 1.0d0 )

! defining susceptible individuals:
        mg1 = size( pack(xx(:,1),xx(:,2).eq. 0.0d0 ) )
        allocate(mmg1(mg1))
        mmg1 = pack(xx(:,1),xx(:,2).eq. 0.0d0 )


! declaring a variable with size equal to the number of susceptible individuals
! the first column is the id number of the susceptible individuals
        allocate(rr(mg1,2))
        rr(:,1) = mmg1

! start calculating the infectivity rate for each susceptible individuals:
        do i = 1,mg1
            allocate(hg(mg))
            do j = 1 , mg

                do m = 1 , ntranspar
                transcov1(m) = transcov(mmg(j),m)**powertrans(m)
                end do

                hg(j) = dot_product(transpar,transcov1)*&
                & (d3(mmg(j),mmg1(i))**(-kernelpar(1)))
            end do

            do m = 1 , nsuspar
                suscov1(m) = suscov(mmg1(i),m)**powersus(m)
            end do

            rr(i,2) = (dot_product(suspar,suscov1)*sum(hg))+spark
            deallocate(hg)
        end do

! assigning waiting time to infection for each susceptible individual:
        do i = 1,mg1
            if(rr(i,2) .eq. 0.0d0)then
                rr(i,2) = Inf
            else
                rr(i,2) = randgamma2(1.0d0,1.0d0/rr(i,2))
            end if
        end do

! choose the one with minmum waiting time as the newly infected individual:
        if(all(rr(:,2) .eq. Inf).eqv. .true.) then
            mms(1,1)  = 0.0d0
            mms(1,2)  = Inf
        else
            mms(1,1:2)  = rr(int(minloc(rr(:,2),1)),:)
        end if

        deallocate(rr)
        deallocate(mmg1)
        deallocate(mmg)

    CASE (3)
! Calculating the infectivity rate of distance-based ILM with spark term and distance kernel = Cauchy.

! defining infectious individuals:
        mg  = size( pack(xx(:,1),xx(:,2).eq. 1.0d0 ) )
        allocate(mmg(mg))
        mmg = pack(xx(:,1),xx(:,2).eq. 1.0d0 )

! defining susceptible individuals:
        mg1 = size( pack(xx(:,1),xx(:,2).eq. 0.0d0 ) )
        allocate(mmg1(mg1))
        mmg1 = pack(xx(:,1),xx(:,2).eq. 0.0d0 )


! declaring a variable with size equal to the number of susceptible individuals
! the first column is the id number of the susceptible individuals
        allocate(rr(mg1,2))
        rr(:,1) = mmg1

! start calculating the infectivity rate for each susceptible individuals:
        do i = 1,mg1
            allocate(hg(mg))
            do j = 1 , mg
                do m = 1 , ntranspar
                    transcov1(m) = transcov(mmg(j),m)**powertrans(m)
                end do
                hg(j) = dot_product(transpar,transcov1)*&
                & (kernelpar(1)/((d3(mmg(j),mmg1(i))**(2.0d0)) +&
                & (kernelpar(1)**(2.0d0))))
            end do

            do m = 1 , nsuspar
                suscov1(m) = suscov(mmg1(i),m)**powersus(m)
            end do

            rr(i,2) = (dot_product(suspar,suscov1)*sum(hg))+spark

            deallocate(hg)

        end do

! assigning waiting time to infection for each susceptible individual:
        do i = 1,mg1
            if(rr(i,2) .eq. 0.0d0)then
                rr(i,2) = Inf
            else
                rr(i,2) = randgamma2(1.0d0,1.0d0/rr(i,2))
            end if
        end do

! choose the one with minmum waiting time as the newly infected individual:
        if(all(rr(:,2) .eq. Inf).eqv. .true.) then
            mms(1,1)  = 0.0d0
            mms(1,2)  = Inf
        else
            mms(1,1:2)  = rr(int(minloc(rr(:,2),1)),:)
        end if

        deallocate(rr)
        deallocate(mmg1)
        deallocate(mmg)


    CASE (4)
! Calculating the infectivity rate of including both distance and network-based ILM with spark term and 
! distance kernel = powerlaw.

! defining infectious individuals:
        mg  = size( pack(xx(:,1),xx(:,2).eq. 1.0d0 ) )
        allocate(mmg(mg))
        mmg = pack(xx(:,1),xx(:,2).eq. 1.0d0 )

! defining susceptible individuals:
        mg1 = size( pack(xx(:,1),xx(:,2).eq. 0.0d0 ) )
        allocate(mmg1(mg1))
        mmg1 = pack(xx(:,1),xx(:,2).eq. 0.0d0 )

! declaring a variable with size equal to the number of susceptible individuals
! the first column is the id number of the susceptible individuals
        allocate(rr(mg1,2))
        rr(:,1) = mmg1

! start calculating the infectivity rate for each susceptible individuals:
        do i = 1,mg1
            allocate(hg(mg))
            do j = 1 , mg
                do m = 1 , ntranspar
                    transcov1(m) = transcov(mmg(j),m)**powertrans(m)
                end do
                hg(j) = dot_product(transpar,transcov1)*&
                & (kernelpar(2)*cc(mmg(j),mmg1(i)) + (d3(mmg(j),mmg1(i))**(-kernelpar(1))))
            end do

            do m = 1 , nsuspar
            suscov1(m) = suscov(mmg1(i),m)**powersus(m)
            end do

            rr(i,2) = (dot_product(suspar,suscov1)*sum(hg))+spark
            deallocate(hg)
        end do

! assigning waiting time to infection for each susceptible individual:
        do i = 1,mg1
            if(rr(i,2) .eq. 0.0d0)then
                rr(i,2) = Inf
            else
                rr(i,2) = randgamma2(1.0d0,1.0d0/rr(i,2))
            end if
        end do

! choose the one with minmum waiting time as the newly infected individual:
        if(all(rr(:,2) .eq. Inf).eqv. .true.) then
            mms(1,1)  = 0.0d0
            mms(1,2)  = Inf
        else
            mms(1,1:2)  = rr(int(minloc(rr(:,2),1)),:)
        end if

        deallocate(rr)
        deallocate(mmg1)
        deallocate(mmg)

    CASE (5)
! Calculating the infectivity rate of including both distance and network-based ILM with spark term and 
! distance kernel = Cauchy.

! defining infectious individuals:
        mg  = size( pack(xx(:,1),xx(:,2).eq. 1.0d0 ) )
        allocate(mmg(mg))
        mmg = pack(xx(:,1),xx(:,2).eq. 1.0d0 )

! defining susceptible individuals:
        mg1 = size( pack(xx(:,1),xx(:,2).eq. 0.0d0 ) )
        allocate(mmg1(mg1))
        mmg1 = pack(xx(:,1),xx(:,2).eq. 0.0d0 )


! declaring a variable with size equal to the number of susceptible individuals
! the first column is the id number of the susceptible individuals
        allocate(rr(mg1,2))
        rr(:,1) = mmg1

! start calculating the infectivity rate for each susceptible individuals:
        do i = 1,mg1
            allocate(hg(mg))
            do j = 1 , mg
                do m = 1 , ntranspar
                    transcov1(m) = transcov(mmg(j),m)**powertrans(m)
                end do
                hg(j) = dot_product(transpar,transcov1)*&
                & (kernelpar(2)*cc(mmg(j),mmg1(i)) + &
                &(kernelpar(1)/((d3(mmg(j),mmg1(i))**(2.0d0)) + &
                &(kernelpar(1)**(2.0d0)))))
            end do

            do m = 1 , nsuspar
                suscov1(m) = suscov(mmg1(i),m)**powersus(m)
            end do

            rr(i,2) = (dot_product(suspar,suscov1)*sum(hg))+spark
            deallocate(hg)
        end do

! assigning waiting time to infection for each susceptible individual:
        do i = 1,mg1
            if(rr(i,2) .eq. 0.0d0)then
                rr(i,2) = Inf
            else
                rr(i,2) = randgamma2(1.0d0,1.0d0/rr(i,2))
            end if
        end do

! choose the one with minmum waiting time as the newly infected individual:
        if(all(rr(:,2) .eq. Inf).eqv. .true.) then
            mms(1,1)  = 0.0d0
            mms(1,2)  = Inf
        else
            mms(1,1:2)  = rr(int(minloc(rr(:,2),1)),:)
        end if

        deallocate(rr)
        deallocate(mmg1)
        deallocate(mmg)


    END SELECT


    end subroutine rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Generating random variables for diffierent distributions !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!####################  NORMAL distribution ######################

    FUNCTION randnormal2(mean,stdev) RESULT(c)

    implicit none

    double precision :: mean,stdev,c,temp(2),r,theta
    double precision, PARAMETER :: PI=3.141592653589793238462d0

        CALL RANDOM_NUMBER(temp)
        r=(-2.0d0*log(temp(1)))**0.5d0
        theta = 2.0d0*PI*temp(2)
        c= mean+stdev*r*sin(theta)

    END FUNCTION randnormal2

!#################### GAMMA distribution ######################

    RECURSIVE FUNCTION randgamma2(shape, SCALE) RESULT(ans)
    double precision:: SHAPE,scale,u,w,d,c,x,xsq,g,ans,v
! DESCRIPTION: Implementation based on "A Simple Method for Generating Gamma Variables"
! by George Marsaglia and Wai Wan Tsang.
! ACM Transactions on Mathematical Software and released in public domain.
! ## Vol 26, No 3, September 2000, pages 363-372.

        IF (shape >= 1.0d0) THEN
            d = SHAPE - (1.0d0/3.0d0)
            c = 1.0d0/((9.0d0 * d)**0.5)
            DO while (.true.)
                x = randnormal2(0.0d0, 1.0d0)
                v = 1.0 + c*x
                DO while (v <= 0.0d0)
                    x = randnormal2(0.0d0, 1.0d0)
                    v = 1.0d0 + c*x
                END DO

            v = v*v*v
            CALL RANDOM_NUMBER(u)
            xsq = x*x
            IF ((u < 1.0d0 -.0331d0*xsq*xsq) .OR.  &
            (log(u) < 0.5d0*xsq + d*(1.0d0 - v + log(v))) )then
                ans=scale*d*v
                RETURN
            END IF

            END DO
        ELSE
            g = randgamma2(shape+1.0d0, 1.0d0)
            CALL RANDOM_NUMBER(w)
            ans=scale*g*(w**(1.0d0/shape))
            RETURN
        END IF

    END FUNCTION randgamma2


!#################### RANDOM SEED ######################

    subroutine initrandomseed2()
!The seed for the random number generation method random_number() has been reset

    implicit none

    integer :: i
    integer :: n
    integer :: clock
    integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(COUNT=clock)

        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(PUT = seed)

        deallocate(seed)
    end subroutine initrandomseed2


end module datsim
