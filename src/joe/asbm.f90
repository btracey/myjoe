!=================================== ASBM =====================================80
!
! Computes the Reynolds stresses using the Algebraic Structure-based Model. 
! On return:
!   ierr = 0 ok
!   ierr = 1 linear solver for AS failed
!   ierr = 2 AS could not be solved in ktr iterations
!   ierr = 3 linear solver for AR failed
!   ierr = 4 AR could not be solved in ktr iterations
!   ierr = 5 negative r_ratio
!   ierr = 6 trace_aa not within bounds
!   ierr = 7 trace of blocking vector out of bounds
!   ierr = 8 flattening parameter is NAN
!
!==============================================================================80

  subroutine asbm(st, wt, wft, bl, as, ar, rey, dmn, cir, ktrmax, ierr)
    implicit none
    integer, parameter  :: dp = selected_real_kind(15, 307) !double precision

    ! INITIAL DECLARATIONS
    ! routine's inputs and outputs
    real(dp), dimension(3,3), intent(in)    :: st  !(Rate of strain)*timescale
    real(dp), dimension(3,3), intent(in)    :: wt  !(Mean rotation)*timescale
    real(dp), dimension(3,3), intent(in)    :: wft !(Frame rotation)*timescale
    real(dp), dimension(3,3), intent(in)    :: bl  !wall blocking tensor
    real(dp), dimension(3,3), intent(inout) :: as  !eddy axis tensor from st
    real(dp), dimension(3,3), intent(inout) :: ar  !eddy axis tensor from st, wt

    real(dp), dimension(3,3), intent(out)   :: rey  !reynolds stresses
    real(dp), dimension(3,3), intent(out)   :: dmn  !dimensionality
    real(dp), dimension(3,3), intent(out)   :: cir  !circulicity

    integer, intent(in)                     :: ktrmax !max iters for N-R
    integer, intent(inout)                  :: ierr   !error flag

    ! constants
    real(dp), parameter                     :: a0 = 1.4_dp
    real(dp), parameter                     :: a01 = (2.1_dp - a0) / 2.0_dp

    real(dp), parameter                     :: zero = 0.0_dp
    real(dp), parameter                     :: fifth = 0.2_dp
    real(dp), parameter                     :: fourth = 0.25_dp
    real(dp), parameter                     :: third = 1.0_dp / 3.0_dp
    real(dp), parameter                     :: half = 0.5_dp
    real(dp), parameter                     :: twoth = 2.0_dp*third
    real(dp), parameter                     :: one = 1.0_dp
    real(dp), parameter, dimension(3,3)     :: delta =                          &
         reshape((/1., 0., 0., 0., 1., 0., 0., 0., 1./),(/3,3/))
    real(dp), parameter, dimension(3,3,3)   :: eps =                            &
         reshape((/0., 0., 0., 0., 0., -1., 0., 1., 0.,                         &
                   0., 0., 1., 0., 0., 0., -1., 0., 0.,                         &
                   0., -1., 0., 1., 0., 0., 0., 0., 0./),(/3,3,3/))
    real(dp), parameter                     :: a_error = 1.0e-10_dp
    real(dp), parameter                     :: r_small = 1.0e-14_dp
    
    ! variables
    integer                  :: i,j,k,l,m,n   !do loop indices
    integer                  :: ktr           !iteration index for N-R
    integer                  :: id, idx       !indeces for N-R vector & matrix
    integer, dimension(3,3)  :: index =                                         & 
      reshape((/1, 2, 3, 2, 4, 5, 3, 5, 6/),(/3, 3/))

    logical                  :: strain        !true for strained flows
    logical                  :: rotation      !true for flows with mean rotation
    logical                  :: rotation_t    !true for flows with total rotation
    logical                  :: converged     !used for N-R solver
    
    real(dp), dimension(3,3) :: a             !eddy axis tensor
    real(dp), dimension(3,3) :: as_old        !input value for as
    real(dp), dimension(3,3) :: ar_old        !input value for ar
    real(dp), dimension(3,3) :: wtt           !frame + mean rotation
    real(dp), dimension(3,3) :: stst          !st_ik*st_kj
    real(dp), dimension(3,3) :: wtwt          !wt_ik*wt_kj
    real(dp), dimension(3,3) :: wtst          !wt_ik*wt_kj
    real(dp)                 :: trace_stst    !st_ik*st_ki
    real(dp)                 :: trace_wtwt    !-wt_ik*wt_ki
    real(dp)                 :: trace_sta     !st_ik*a_ki
    real(dp)                 :: trace_ststa   !st_kp*st_kq*a_pq 
    real(dp)                 :: trace_wtsta   !wt_qr*st_rp*a_pq
    real(dp)                 :: trace_wttwtt  !-wTt_ij*wTt_ji
    real(dp)                 :: trace_aa      !a_ij*a_ji
    real(dp)                 :: trace_ast     !a_ij*st_ji
    real(dp)                 :: trace_bl      !bl_ii
    real(dp)                 :: mag_ststa     !sqrt(st_kp*st_kq*a_pq)
    real(dp)                 :: norm_st       !norm of st_ik       
    real(dp)                 :: norm_wt       !norm of wt_ik
    
    real(dp), dimension(3)   :: vec_wtt       !vector for frame + mean rotation
    real(dp), dimension(3)   :: vec_wdt       !vector for frame - mean rotation
    real(dp)                 :: dot_vec_wtt   !dot product of vec_wtt
    real(dp)                 :: dot_vec_wdt   !dot product of vec_wdt
    real(dp)                 :: mag_vec_wtt   !magnitude of vec_wtt

    real(dp)                 :: num_as        !numerator for as equation
    real(dp)                 :: den_as        !denominator for as equation
    real(dp), dimension(6)   :: x_nr          !rhs and solution for N-R
    real(dp), dimension(6,6) :: a_nr          !matrix for N-R
    real(dp)                 :: r_ratio       !wt_qr*st_rp*a_pq/st_kn*st_nm*a_mk
    real(dp)                 :: r_abs         !abs(r_ratio)
    real(dp)                 :: r_num         !numerator of r_ratio
    real(dp)                 :: alpha         !2nd coeff for rot matrix H
    real(dp)                 :: beta          !3rd coeff for rot matrix H
    real(dp)                 :: alpha2        !alpha squared
    real(dp)                 :: aroot, broot  !needed for alpha & beta
    real(dp)                 :: term          !needed for beta
    real(dp)                 :: c2            !normalized 2nd coeff for H
    real(dp)                 :: c3            !normalized 3rd coeff for H
    real(dp)                 :: coefa
    real(dp)                 :: coefb
    real(dp)                 :: coefp
    real(dp), dimension(3,3) :: h             !rotation matrix H
    real(dp), dimension(3,3) :: p
    real(dp)                 :: norm_x_nr     !norm of x_nr
    real(dp)                 :: min_sw_ws     !min(strain/rot, rot/strain)

    real(dp)                 :: hat_wt        !-a_ij*wt_ik*wt_kj
    real(dp)                 :: hat_wtt       !-a_ij*wTt_ik*wTt_kj  
    real(dp)                 :: hat_st        !a_ij*st_ik*st_kj
    real(dp)                 :: hat_x         !a_ij*wTt_ik*st_kj
    real(dp)                 :: eta_r         !mean rot over strain parameter
    real(dp)                 :: eta_f         !frame rot over strain parameter
    real(dp)                 :: eta_c1        !used to compute eta_r
    real(dp)                 :: eta_c2        !used to compute eta_f
    real(dp)                 :: oma           !one minus trace_aa
    real(dp)                 :: sqamth        !sqrt(trace_aa - third)

    real(dp)                 :: phi           !jettal scalar
    real(dp)                 :: bet           !correlation scalar
    real(dp)                 :: chi           !flattening scalar

    real(dp)                 :: phis          !jettal scalar for any a
    real(dp)                 :: bets          !correlation scalar for any a
    real(dp)                 :: chis          !flattening scalar for any a

    real(dp)                 :: phi1          !jettal scalar for shear
    real(dp)                 :: bet1          !correlation scalar for shear
    real(dp)                 :: chi1          !flattening scalar for shear

    real(dp)                 :: scl_g         !helical scalar
    real(dp), dimension(3)   :: vec_g         !helical vector

    real(dp)                 :: struc_weight  !smoothing parameter
    real(dp)                 :: xp_aa         !extrapolation along trace_aa

    continue

    ! INITIALIZE AND CALCULATE TENSOR PRODUCTS
    ierr = 0
    
    trace_stst = zero
    trace_wtwt = zero
    
    ! initialize default value for A
    a(1,1) = third
    a(2,1) = zero
    a(3,1) = zero

    a(1,2) = zero
    a(2,2) = third
    a(3,2) = zero
    
    a(1,3) = zero
    a(2,3) = zero
    a(3,3) = third
    
    ! tensor products
    stst = zero
    wtwt = zero
    wtst = zero    
  
    do i = 1,3
      do j = 1,3
        do k = 1,3
          stst(i,j) = stst(i,j) + st(i,k)*st(k,j)
          wtwt(i,j) = wtwt(i,j) + wt(i,k)*wt(k,j)
          wtst(i,j) = wtst(i,j) + wt(i,k)*st(k,j)
        end do
      end do

      trace_stst = trace_stst + stst(i,i)
      trace_wtwt = trace_wtwt - wtwt(i,i)
    end do

    if (trace_stst > zero) then
      strain = .true.
      norm_st = sqrt(trace_stst)
    else
      strain  = .false.
      norm_st = zero
    end if

    if (trace_wtwt > zero) then
      rotation = .true.
      norm_wt = sqrt(trace_wtwt)
    else
      rotation = .false.
      norm_wt  = zero
    end if

    ! COMPUTE AS
    ! check for strain
    purestrain: if (strain) then      
      ! reuse as and store original tensor
      a = as
      as_old = as

      ! loop point for Newton-Rhapson (N-R) iteration
      converged = .false.
      ktr = 0

      NewtonRhapson_as: do while (.not.converged)
        ktr = ktr + 1
        ! perturbations
        if (mod(ktr,20) == 0) then
          do i = 1,3
            do j = 1,3
              if (mod(ktr,40) == 0) then
                a(i,j) = 0.5_dp*a(i,j)
              else if (mod(ktr,200) == 0) then
                a(i,j) = a(i,j) + as_old(i,j)
              else 
                a(i,j) = 2.0_dp*a(i,j)
              end if
            end do
          end do
        end if

        ! invariants
        trace_sta = zero
        trace_ststa = zero
        do i = 1,3
          do j = 1,3
            trace_sta = trace_sta + st(i,j)*a(j,i)
            trace_ststa = trace_ststa + stst(i,j)*a(j,i)
          end do
        end do

        ! denominator
        if (trace_ststa > zero) then
          ! proper value
          mag_ststa = sqrt(a01**2 + trace_ststa)
        else
          ! flase value to deal with numberical problems
          mag_ststa = sqrt(trace_stst)
        end if
        den_as = a0 + 2*mag_ststa
       
        ! compute matrix and rhs for as perturbation
        do i = 1,3
          do j = i,3
            ! terms in the as(i,j) equations
            id = index(i,j)

            num_as = - twoth*trace_sta*delta(i,j)
            do k = 1,3
              num_as = num_as + st(i,k)*a(k,j) + st(j,k)*a(k,i)
            end do

            ! right hand side of equation for perturbation in as(i,j)
            x_nr(id) = -a(i,j) + third*delta(i,j) + num_as/den_as

            ! initialize row in the matrix am
            do k = 1,6
              a_nr(id,k) = zero
            end do
            a_nr(id,id) = one

            ! terms for tensor
            do k = 1,3
              idx = index(k,j)
              a_nr(id,idx) = a_nr(id,idx) - st(i,k)/den_as
              idx = index(k,i)
              a_nr(id,idx) = a_nr(id,idx) - st(j,k)/den_as
            end do

            ! terms from deltas
            if (i == j) then
              do n = 1,3
                do m = 1,3
                  idx = index(n,m)
                  a_nr(id,idx) = a_nr(id,idx) + twoth*st(m,n)/den_as
                end do
              end do
            end if

            ! term from denominator (a1 = 2)
            do n = 1,3
              do m = 1,3
                idx = index(n,m)
                a_nr(id,idx) = a_nr(id,idx) +                                   &
                   num_as/den_as*stst(n,m)/(mag_ststa*den_as)
              end do
            end do
             
            ! row for i,j complete
          end do
        end do
        
        ! solve the system
        call linsolver(6,a_nr,x_nr,6,ierr)
        if (ierr /= 0) then
          ierr = 1
          return
        endif

        ! compute the corrected solution
        do i = 1,3
          do j = 1,3
            id = index(i,j)
            a(i,j) = a(i,j) + x_nr(id)
          end do
        end do

        ! check convergence
        converged = .true.

        do i = 1,6
          if (abs(x_nr(i)) > a_error) converged = .false.        
        end do

        if (ktr == ktrmax) then
          converged = .true.
          ierr = 2
          return
        end if

      end do NewtonRhapson_as

    end if purestrain

    ! update as
    as = a

    ! COMPUTE AR
    ! check for rotation
    purerotation: if (rotation) then
      ! reuse ar and store original tensor
      a = ar
      ar_old = ar
     
      ! loop point for Newton-Rhapson (N-R) iteration
      converged = .false.
      ktr = 0

      NewtonRhapson_ar: do while (.not.converged)
        ktr = ktr + 1
        ! perturbations
        if (mod(ktr,20) == 0) then
          do i = 1,3
            do j = 1,3
              if (mod(ktr,40) == 0) then
                a(i,j) = 0.5_dp*a(i,j)
              else if (mod(ktr,200) == 0) then
                a(i,j) = a(i,j) + ar_old(i,j)
              else 
                a(i,j) = 2.0_dp*a(i,j)
              end if
            end do
          end do
        end if

        ! compute rotation parameters
        trace_wtsta = zero
        trace_ststa = zero
        do i = 1,3
          do j = 1,3
            trace_wtsta = trace_wtsta + a(i,j)*wtst(j,i)
            trace_ststa = trace_ststa + stst(i,j)*a(j,i)
          end do
        end do

        if (ktr == 1) then
          ! relative rotation parameter for first trial
          r_ratio = one
        else
          ! relative rotation parameter for following trials
          r_ratio = trace_wtsta/trace_ststa      
          if (abs(r_ratio) < r_small) exit
        end if
        r_abs = abs(r_ratio)
       
        ! compute alpha and beta
        if (r_abs < one) then
          ! hyperbolic mean flow
          aroot = sqrt(one - r_abs)
          alpha2 = one - aroot
        else
          ! elliptic mean flow
          aroot = sqrt(one - one/r_abs)
          alpha2 = one + aroot
        end if
        if (alpha2 < zero) alpha2 = zero
        alpha = sqrt(alpha2)
        
        term = 4.0_dp - 2.0_dp*alpha2
        if (term < zero) term = zero
        broot = sqrt(term)
        beta = 2.0_dp - broot

        c2 = alpha/norm_wt
        c3 = beta/trace_wtwt
        
        ! compute rotation matrix h
        do i = 1,3
          do j = 1,3
            h(i,j) = delta(i,j) + c2*wt(i,j) + c3*wtwt(i,j)
          end do
        end do

        ! compute initial x_nr
        do i = 1,3
          do j = i,3
            id = index(i,j)
            x_nr(id) = zero

            do k = 1,3
              do l = 1,3
                x_nr(id) = x_nr(id) + h(i,k)*h(j,l)*as(k,l)
              end do
            end do

          end do
        end do
       
        ! solution for first trial
        if (ktr == 1) then
          do i = 1,3
            do j = i,3
              id = index(i,j)
              ! first pass solution
              a(i,j) = x_nr(id)
              a(j,i) = a(i,j)
            end do
          end do
          ! done for the first trial
          cycle
        end if

        ! coefficients in the a_nr matrix 
        coefa = fourth/(aroot*alpha*norm_wt)
        coefb = half/(aroot*broot*trace_wtwt)
        
        if (r_abs < one) then
          coefa = coefa*r_ratio
          coefb = coefb*r_ratio
        else
          coefa = coefa/r_ratio
          coefb = coefb/r_ratio
        end if

        ! p(k,l) for nr calculation
        do k = 1,3
          do l = 1,3
            p(k,l) = wtst(k,l)/trace_wtsta - stst(k,l)/trace_ststa
          end do
        end do    

        ! finish x_nr and compute a_nr
        do i = 1,3
          do j = i,3
            id = index(i,j)
            x_nr(id) = x_nr(id) - a(i,j)
              
            ! initialize row in the matrix 
            do k = 1,6
              a_nr(id,k) = zero
            end do
            a_nr(id,id) = one
              
            ! coefficient of the p(m,n) for this i,j
            coefp = zero
            do k = 1,3
              do l = 1,3
                coefp = coefp + ((coefa*wt(i,k) + coefb*wtwt(i,k))*h(j,l) +     &
                                 (coefa*wt(j,l) + coefb*wtwt(j,l))*h(i,k))*     &
                                as(k,l)
              end do
            end do

            do n = 1,3
              do m = 1,3
                idx = index(n,m)
                a_nr(id,idx) = a_nr(id,idx) - coefp*p(m,n)
              end do
            end do
   
            ! row for i,j complete
          end do
        end do

        ! solve the system
        call linsolver(6,a_nr,x_nr,6,ierr)
        if (ierr /= 0) then
          ierr = 3
          return
        end if

        ! compute the corrected solution
        ! The following attempts to limit the size of the correction.
        ! Min norm for aij should be same as the identity matrix.
        ! Using it here as a basis for max size. No real meaning. Could 
        ! use 0.1 too.

        ! check size
        norm_x_nr = x_nr(1)**2 + x_nr(4)**2 + x_nr(6)**2 +                      &
                 2*(x_nr(2)**2 + x_nr(3)**2 + x_nr(5)**2)
        if (norm_x_nr > third) then
          do i = 1,6
            x_nr(i) = third*x_nr(i)/sqrt(norm_x_nr)
          end do
        end if

        ! check numerator of r_ratio
        r_num = -one
        do while (r_num < zero)
          r_num = zero
          
          do i = 1,3
            do j = 1,3 
              id = index(i,j)
              r_num = r_num + (a(i,j) + x_nr(id))*wtst(j,i)
            end do
          end do
          
          if (r_num < zero) then
            r_num = zero
            ! trying to reduce x depending on the flow. Arbitrary
            min_sw_ws = min( (trace_stst/trace_wtwt),(trace_wtwt/trace_stst) )  &
                        **fifth   
            do i = 1,6
              x_nr(i) = min_sw_ws*x_nr(i)
            end do
          end if

        end do
       
        ! compute the corrected solution
        do i = 1,3
          do j = 1,3
            id = index(i,j)
            a(i,j) = a(i,j) + x_nr(id)
          end do
        end do

        ! check convergence
        converged = .true.

        do i = 1,6
          if (abs(x_nr(i)) > a_error) converged = .false.
        end do

        if (ktr == ktrmax) then
          converged = .true.
          ierr = 4
          ! do not return when this error occurs
          write(*,*) 'asbm error, ierr = ', ierr
        end if

      end do NewtonRhapson_ar
      
      ! check that r_ratio is not negative
      if (r_ratio < -zero) then
        ierr = 5
        return
      end if

    end if purerotation

    ! update ar
    ar = a

    ! COMPUTE STRUTURE SCALARS
    trace_wttwtt = zero
    trace_ststa  = zero
    trace_aa     = zero
    trace_ast    = zero
    
    do i = 1,3
      do j = 1,3
        wtt(i,j)     = wft(i,j) + wt(i,j)
        trace_wttwtt = trace_wttwtt + wtt(i,j)*wtt(i,j)
        trace_ststa  = trace_ststa + a(i,j)*stst(i,j)
        trace_aa     = trace_aa + a(i,j)*a(j,i)
        trace_ast    = trace_ast + a(i,j)*st(j,i)
      end do
    end do

    ! bounds check for trace_aa
    if (trace_aa > one) then
      write(*,*) 'trace_aa = ', trace_aa, ' greater than one'
      trace_aa = one
      ierr = 6
    end if
    
    if (trace_aa < third) then
      write(*,*) 'trace_aa = ', trace_aa, ' less than third'
      trace_aa = third
      ierr = 6
    end if

    ! relative rotation rate and frame rotation rate vector
    vec_wtt(1) = wtt(2,3)
    vec_wtt(2) = wtt(3,1)
    vec_wtt(3) = wtt(1,2)

    vec_wdt(1) = wft(3,2) - wt(3,2)
    vec_wdt(2) = wft(1,3) - wt(1,3)
    vec_wdt(3) = wft(2,1) - wt(2,1)
    
    ! their magninuteds
    dot_vec_wtt = vec_wtt(1)**2.0_dp + vec_wtt(2)**2.0_dp + vec_wtt(3)**2.0_dp
    dot_vec_wdt = vec_wdt(1)**2.0_dp + vec_wdt(2)**2.0_dp + vec_wdt(3)**2.0_dp

    ! check for total rotation
    if (dot_vec_wtt > zero) then
      mag_vec_wtt = sqrt(dot_vec_wtt)
      rotation_t = .true.
    else
      mag_vec_wtt = zero
      rotation_t = .false.
    end if

    ! compute quantities based on total rotation rate
    hat_wt  = zero        
    hat_wtt = zero 
    hat_st  = zero                
    hat_x   = zero  
       
    do i = 1,3
      do j = 1,3
        do k = 1,3
          hat_wt  = hat_wt + a(i,j)*wt(i,k)*wt(j,k)
          hat_wtt = hat_wtt + a(i,j)*wtt(i,k)*wtt(j,k)
          hat_st  = hat_st + a(i,j)*st(i,k)*st(j,k)
          hat_x   = hat_x + a(i,j)*wtt(j,k)*st(k,i)
        end do
      end do
    end do

    ! compute structure parameters
    eta_c1 = hat_wt/(hat_st)
    !eta_c1 = hat_wt/hat_st
    eta_c2 = hat_wtt/hat_st
    
    if (eta_c1 < zero) eta_c1 = zero
    if (eta_c2 < zero) eta_c2 = zero

    eta_r = sqrt(eta_c1)
    eta_f = eta_r - sign(sqrt(eta_c2),hat_x)

    ! compute phis, chis, bets
    if (hat_st < zero) then
      ! without strain
      phis = zero
      chis = zero
      bets = one

      if (hat_wtt > zero) then
        ! with rotation
        phis = third
        chis = zero
        bets = zero
      end if
    else
      ! with strain
      oma = 1 - trace_aa
      sqamth = sqrt(trace_aa - third)

      if (eta_r < one) then
        call int_er_lt_one(eta_r,eta_f,oma,sqamth,phis,bets,chis,phi1,bet1,chi1)
      else 
        call int_er_gt_one(eta_r,eta_f,oma,sqamth,phis,bets,chis,phi1,bet1,chi1)
      end if
      struc_weight = exp(-100.0_dp*abs(eta_r - one)**2.0_dp)
      !phis = phis*(one - struc_weight) + phi1*(struc_weight)
      bets = bets*(one - struc_weight) + bet1*(struc_weight)
      chis = chis*(one - struc_weight) + chi1*(struc_weight)
    end if
   
    ! compute phi, chi, bet
    xp_aa = 1.5_dp*(trace_aa - third)
    xp_aa = 0.35_dp*xp_aa**2.5_dp + 0.65_dp*xp_aa**half

    phi = phis*xp_aa
    chi = chis*xp_aa
    if (eta_r < one) then
      bet = bets
    else
      bet = one - max(one - 0.9_dp*(eta_r - one)**0.31_dp, zero)*               &
                  (1.5_dp*(trace_aa - third))**10.0_dp
      !a(1,2) = a(1,2)*(one - 1.5_dp*(trace_aa - third))**0.1_dp
      !a(2,1) = a(1,2)
    end if
    struc_weight = exp(-100.0_dp*abs(eta_r - one)**2.0_dp)
    bet = bet*(one - struc_weight) + bet1*(struc_weight)

    ! compute helical scalar
    scl_g = 2.0_dp*phi*(one - phi)/(one + chi)
    if (scl_g < zero) scl_g = zero
    scl_g = bet*sqrt(scl_g)

    ! IMPOSE BLOCKAGE MODIFICATIONS
    ! blockage correction to eddy-axis tensor
    call blocking(a,bl,delta,ierr)

    ! blockage correction to phi and gamma
    trace_bl = bl(1,1) + bl(2,2) + bl(3,3)

    phi = one - (one - trace_bl)**2 + phi*(one - trace_bl)**2
    scl_g = (one - trace_bl)*scl_g
    !chi = (one - trace_bl)*chi
    
    if (rotation_t) then
    ! compute helical vector
      do k = 1,3
        vec_g(k) = -scl_g*vec_wtt(k)/mag_vec_wtt
      end do
    else
      do k = 1,3
        vec_g(k) = zero
      end do
    end if

    ! COMPUTE STRUCTURE TENSORS
    call structure(rey, dmn, cir, a, phi, chi, vec_g, vec_wdt,           &
                   dot_vec_wdt, rotation_t, delta, eps, ierr)

    ! Output some useful data
    cir(1,1) = phi
    cir(1,2) = bet
    cir(1,3) = chi

    cir(2,1) = eta_r
    cir(2,2) = eta_f
    cir(2,3) = trace_aa
    
    cir(3,1) = hat_wt
    cir(3,2) = hat_st
    cir(3,3) = trace_ststa
  end subroutine asbm

!================================= LINSOLVER ==================================80
!
! Double precision linear algebraic equation solver.
! Solves a(i,j)*x(j) = b(i) i = 1,...,n
! fround is a round-off test factor assuming at least 15 digit accuracy.
! Uses Gaussian elimination with row normalization and selection.
! On return:
!   Solution ok:
!     ierr = 0
!     x is in b
!     a is destroyed
!   Singular matrix or bad dimensioning:
!     ierr = 1
!     a and b are destroyed
!
!==============================================================================80

  subroutine linsolver(ndim,a,b,n,ierr)
    implicit none
    integer, parameter  :: dp = selected_real_kind(15, 307) !double precision

    ! INITIAL DECLARATIONS
    ! routine's inputs and outputs
    integer, intent(in)                           :: ndim   !number of columns
    integer, intent(in)                           :: n      !number of rows
    integer, intent(inout)                        :: ierr   !error flag
    real(dp), dimension(ndim,ndim), intent(inout) :: a      !matrix
    real(dp), dimension(ndim), intent(inout)      :: b      !vector  
    
    ! constants
    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: one = 1.0_dp
    real(dp), parameter :: fround = 1.0e-15_dp

    ! variables
    integer             :: i, j, k
    integer             :: ix = 0.0_dp 
    integer             :: nm1, m, mm1, im1
    real(dp)            :: c, cx, d
    real(dp)            :: termb, term

    continue

    if (n > ndim) then
      ! dimensioning error: treat as program error in calls
      ierr = 1
      return
    end if

    ! form the lower triagonal system
    if (n > 1) then
      nm1 = n - 1
      kloop: do k = 1,nm1
        ! eliminate the last column in the upper M x M matrix
        m = n - k + 1
        mm1 = m - 1

        ! normalize each row on its largest element
        iloop1: do i = 1,m
          c = zero         
          do j = 1,m
            cx = abs(a(i,j))
            if (cx > c) c = cx
          end do

          ! check for singular matrix
          if (c == zero) then
            ierr = 1
            return
          end if

          c = one/c
          do j = 1,m
            a(i,j) = a(i,j)*c
          end do

          b(i) = b(i)*c
        end do iloop1

        ! find the best row ix to eliminate in column M
        c = zero
        iloop2: do i = 1,m
          cx = abs(a(i,m))
          if (cx > c) then
            c = cx
            ix = i
          end if
        end do iloop2

        ! check for singular matrix
        if (c == zero) then
          ierr = 1
          return
        end if

        if (m /= ix) then
          ! switch rows M and ix
          c = b(m)
          b(m) = b(ix)
          b(ix) = c
          do j = 1,m
            c = a(m,j)
            a(m,j) = a(ix,j)
            a(ix,j) = c
          end do
        end if

        ! eliminate last column using the lowest row in the M x M matrix
        
        ! check for singular matrix
        if (a(m,m) == zero) then
          ierr = 1
          return
        end if

        ! column loop
        iloop3: do i = 1,mm1
          ! check for column entry
          if (a(i,m) /= zero) then
            ! eliminate
            c = a(i,m)/a(m,m)
            d = b(i) - c*b(m)
            if (abs(d) < fround*abs(b(i))) d = zero
            b(i) = d
            do j = 1,mm1
              d = a(i,j) - c*a(m,j)
              if (abs(d) < fround*abs(a(i,j))) d = 0
              a(i,j) = d
            end do
          end if
        end do iloop3

      end do kloop  
    end if
    
    ! compute the back solution
    
    ! check for singular matrix
    if (a(1,1) == zero) then
      ierr = 1
      return
    end if

    ! calculate x(1)
    b(1) = b(1)/a(1,1)
    if (n > 1) then
      iloop4: do i = 2,n
        ! calculate x(i)
        im1 = i - 1
        c = b(i)
        termb = abs(c)
        do j = 1,im1
          term = a(i,j)*b(j)
          c = c - term
          if (abs(term) > termb) termb = abs(term)
        end do
      
        if (abs(c) < fround*termb) c = zero
        b(i) = c/a(i,i)
      end do iloop4
    end if

    ! normal exit
    ierr = 0
  end subroutine linsolver

!=============================== INT_ER_LT_ONE ================================80
!
! Computes the structure parameters for an arbitrary a_ij by interpolating 
! between the plane strain and pure shear states.
!
!==============================================================================80

  subroutine int_er_lt_one(eta_r,eta_f,oma,sqamth,phis,bets,chis,phi1,bet1,chi1)
    implicit none
    integer, parameter  :: dp = selected_real_kind(15, 307) !double precision

    ! INITIAL DECLARATIONS
    ! routine's inputs and outputs
    real(dp), intent(in)    :: eta_r
    real(dp), intent(in)    :: eta_f
    real(dp), intent(in)    :: oma
    real(dp), intent(in)    :: sqamth
    real(dp), intent(inout) :: phis
    real(dp), intent(inout) :: bets
    real(dp), intent(inout) :: chis
    real(dp), intent(inout) :: phi1
    real(dp), intent(inout) :: bet1
    real(dp), intent(inout) :: chi1

    ! constants
    real(dp), parameter     :: zero = 0.0_dp
    real(dp), parameter     :: one = 1.0_dp
    
    ! variables
    real(dp)                :: eta_f1, eta_f0
    real(dp)                :: phi0
    real(dp)                :: bet0
    real(dp)                :: chi0
    real(dp)                :: param 
    
    continue

    param = sqrt(3.0_dp)/4.0_dp
    if (eta_f < (param*(eta_r - one))) then 
      eta_f1 = eta_f - param*(eta_r - one)
      eta_f0 = eta_f - param*(eta_r - one) - param
    else if (eta_f < zero) then
      eta_f1 = zero
      eta_f0 = eta_f/(1 - eta_r)
    else if (eta_f < eta_r) then
      eta_f1 = eta_f/eta_r
      eta_f0 = zero
    else if (eta_f < (eta_r + param*(one - eta_r))) then
      eta_f1 = one
      eta_f0 = 1 - (1 - eta_f)/(1 - eta_r)
    else
      eta_f1 = one + eta_f - eta_r + param*(one - eta_r)
      eta_f0 = param + eta_f - eta_r + param*(one - eta_r)
    end if

    ! compute parameters for shear and plane strain states
    call pure_shear(eta_f1,oma,sqamth,phi1,bet1,chi1)
    call plane_strain(eta_f0,oma,sqamth,phi0,bet0,chi0)
    
    ! interpolate along eta_r direction
    phis = phi0 + (phi1 - phi0)*                                                &
                  (0.82_dp*eta_r**2)/(one - (one-0.82_dp)*eta_r**2)
    bets = bet0 + (bet1 - bet0)*eta_r**2.0_dp
    chis = chi0 + (chi1 - chi0)*eta_r**2.0_dp
  
  end subroutine int_er_lt_one

!=============================== INT_ER_GT_ONE ================================80
!
! Computes the structure parameters for an arbitrary a_ij by extrapolating 
! for higher values of eta_r than that of pure shear.
!
!==============================================================================80

  subroutine int_er_gt_one(eta_r,eta_f,oma,sqamth,phis,bets,chis,phi1,bet1,chi1)
    implicit none 
    integer, parameter  :: dp = selected_real_kind(15, 307) !double precision

    ! INITIAL DECLARATIONS
    ! routine's inputs and outputs
    real(dp), intent(in)    :: eta_r
    real(dp), intent(in)    :: eta_f
    real(dp), intent(in)    :: oma
    real(dp), intent(in)    :: sqamth
    real(dp), intent(inout) :: phis
    real(dp), intent(inout) :: bets
    real(dp), intent(inout) :: chis
    real(dp), intent(inout) :: phi1
    real(dp), intent(inout) :: bet1
    real(dp), intent(inout) :: chi1

    ! constants
    real(dp), parameter     :: third = 1.0_dp / 3.0_dp
    real(dp), parameter     :: half = 1.0_dp / 2.0_dp
    real(dp), parameter     :: one = 1.0_dp

    ! variables
    real(dp)                :: aux
    real(dp)                :: aux_shift

    continue

    aux = one/(one + half*(eta_r - one)/(oma)**2.5_dp)
    aux_shift = one/(one + half*(eta_r - one)/(oma+third)**2.5_dp)

    ! compute parameters for shear state
    call pure_shear(eta_f,oma,sqamth,phi1,bet1,chi1)
    
    ! extrapolate along eta_r direction
    phis = third + (phi1 - third)*aux_shift
    bets = bet1*aux
    chis = chi1*aux
  
  end subroutine int_er_gt_one

!================================ PURE_SHEAR ==================================80
!
! Computes reference structure parameters along the shear line.
!
!==============================================================================80

  subroutine pure_shear(eta_f,oma,sqamth,phi1,bet1,chi1)
    implicit none 
    integer, parameter  :: dp = selected_real_kind(15, 307) !double precision

    ! INITIAL DECLARATIONS
    ! routine's inputs and outputs
    real(dp), intent(in)    :: eta_f
    real(dp), intent(in)    :: oma
    real(dp), intent(in)    :: sqamth
    real(dp), intent(inout) :: phi1
    real(dp), intent(inout) :: bet1
    real(dp), intent(inout) :: chi1

    ! constants
    real(dp), parameter     :: zero = 0.0_dp
    real(dp), parameter     :: fifth = 0.2_dp
    real(dp), parameter     :: one = 1.0_dp

    continue

    if (eta_f < zero) then
      phi1 = (eta_f - one)/(3.0_dp*eta_f - one)
      bet1 = one/(one - eta_f*(1 + sqamth)/oma)
      chi1 = fifth*bet1
    else if (eta_f < one) then
      phi1 = one - eta_f
      chi1 = fifth + (one - fifth)*                                             &
                     (one - (one - eta_f)**2/(one + 3.0_dp*eta_f/oma))
      bet1 = one
    else
      phi1 = (eta_f - one)/(3.0_dp*eta_f - one)
      bet1 = one/(one + 0.8_dp*(eta_f - one)*(eta_f*sqamth)/oma)
      chi1 = one - (one - bet1)*(eta_f - one)/(oma + eta_f - one)
    end if

  end subroutine pure_shear

!============================== PLANE_STRAIN ==================================80
!
! Computes reference structure parameters along the shear line.
!
!==============================================================================80

  subroutine plane_strain(eta_f,oma,sqamth,phi0,bet0,chi0)
    implicit none 
    integer, parameter  :: dp = selected_real_kind(15, 307) !double precision

    ! INITIAL DECLARATIONS
    ! routine's inputs and outputs
    real(dp), intent(in)    :: eta_f
    real(dp), intent(in)    :: oma
    real(dp), intent(in)    :: sqamth
    real(dp), intent(inout) :: phi0
    real(dp), intent(inout) :: bet0
    real(dp), intent(inout) :: chi0

    ! constants
    real(dp), parameter     :: zero = 0.0_dp
    real(dp), parameter     :: one = 1.0_dp
    real(dp), parameter     :: two = 2.0_dp
    real(dp), parameter     :: three = 3.0_dp

    ! variables
    real(dp)                :: var1
    real(dp)                :: var2

    continue
    
    var1 = two*eta_f
    var2 = var1*var1

    if ( var2 > 0.75_dp) then
      bet0 = one/(one + (var1 - sqrt(0.75_dp))*(var1*sqamth)/oma)
      phi0 = (one - bet0)/three
      chi0 = -bet0
    else
      phi0 = 0.145_dp*(var2/0.75_dp - (var2/0.75_dp)**9)
      chi0 = -(0.342_dp*(var2/0.75_dp) + (one - 0.342_dp)*(var2/0.75_dp)**6)
      bet0 = one
    end if
    chi0 = zero

  end subroutine plane_strain

!================================ BLOCKING ====================================80
!
! Modifies the eddy axis tensor to account for wall effects.
! Input:
!        a(i,j)       unblocked (homogeneous) tensor
!        bl(i,j)      blocking tensor
! Output
!        a(i,j)       blocked tensor
!
!==============================================================================80

  subroutine blocking(a,bl,delta,ierr)
    implicit none 
    integer, parameter  :: dp = selected_real_kind(15, 307) !double precision

    ! INITIAL DECLARATIONS
    ! routine's inputs and outputs
    real(dp), dimension(3,3), intent(inout) :: a        !eddy axis tensor
    real(dp), dimension(3,3), intent(in)    :: bl       !blocking tensor
    real(dp), dimension(3,3), intent(in)    :: delta    !kronecker delta
 
    integer, intent(inout)                  :: ierr  !error flag

    ! constants
    real(dp)                 :: small = 1.0e-14_dp
    real(dp)                 :: zero = 0.0_dp
    real(dp)                 :: one = 1.0_dp

    ! variables
    integer                  :: i,j,k,l    !do loop indeces
    real(dp), dimension(3,3) :: ah         !homogeneous eddy axis tensor
    real(dp)                 :: trace_bl   !bl_ii
    real(dp)                 :: trace_ahbl !ah_in*bl_ni
    real(dp), dimension(3,3) :: h          !partial projection operator
    real(dp)                 :: d2         !normalizing factor
    real(dp)                 :: dinv       !one over sqrt(d2)
    real(dp)                 :: sum

    continue

    trace_bl = bl(1,1) + bl(2,2) + bl(3,3)

    ! return if no blockage
    if (trace_bl < small) return

    ! check for bad data
    if (trace_bl > one) then
      ierr = 7
      return
    end if

    ah = a

    ! compute normalizing factor
    sum = zero
    do i = 1,3
      do j = 1,3
        sum = sum + ah(i,j)*bl(i,j)
      end do
    end do
    d2 = one - (2.0_dp - trace_bl)*sum
    if (d2 < zero) return
    dinv = one/sqrt(d2);

    ! compute partial projection operator
    do i = 1,3
      do j = 1,3
        h(i,j) = dinv*(delta(i,j) - bl(i,j))
      end do
    end do

    ! apply blockage correction
    a = zero
    do i = 1,3
      do j = i,3
        do k = 1,3
          do l = 1,3
            a(i,j) = a(i,j) + h(i,k)*h(j,l)*ah(k,l)
          end do
        end do
        a(j,i) = a(i,j)
      end do
    end do

  end subroutine blocking

!================================ STRUCTURE ===================================80
!
! Calculates the normalized reynolds stress, dimensionality, and circulicity
! tensors.
!
!==============================================================================80

  subroutine structure(rey, dmn, cir, a, phi, chi, vec_g, vec_wdt,       &
       dot_vec_wdt, rotation_t, delta, eps, ierr)
    implicit none 
    integer, parameter  :: dp = selected_real_kind(15, 307) !double precision

    ! INITIAL DECLARATIONS
    ! routine's inputs and outputs
    real(dp), dimension(3,3), intent(out)   :: rey     !reynolds stresses
    real(dp), dimension(3,3), intent(out)   :: dmn     !dimensionality
    real(dp), dimension(3,3), intent(out)   :: cir     !circulicity
    real(dp), dimension(3,3), intent(in)    :: a       !eddy axis tensor
    real(dp), dimension(3,3), intent(in)    :: delta   !kronecker delta
    real(dp), dimension(3,3,3), intent(in)  :: eps     !eps_ikj

    real(dp), intent(in)               :: phi     !jettal scalar
    real(dp), intent(in)               :: chi     !flattening scalar
    real(dp), dimension(3), intent(in) :: vec_g   !helical vector
    real(dp), dimension(3), intent(in) :: vec_wdt !vec for frame - mean rotation
    real(dp), intent(in)               :: dot_vec_wdt !dot product of vec_wdt

    logical, intent(in)                :: rotation_t

    integer, intent(inout)             :: ierr

    ! constants
    real(dp), parameter      :: zero = 0.0_dp
    real(dp), parameter      :: small = 1.0e-14_dp
    real(dp), parameter      :: half = 0.5_dp
    real(dp), parameter      :: one = 1.0_dp

    ! variables
    integer                  :: i, j, k, l, m
    real(dp), dimension(3,3) :: fl             !flattening tensor
    real(dp), dimension(3,3) :: afl            !a_in*fl_nj
    real(dp)                 :: trace_afl      !a_in*fl_ni
    real(dp)                 :: coef_1, coef_2, coef_3, coef_4, coef_5
    real(dp)                 :: term
    real(dp)                 :: sum, sumi, sumj, sumk, sumg

    continue
    
    ! compute axisymmetric tensors
    do i = 1,3
      do j = i,3
        rey(i,j) = half*(1 - phi)*(delta(i,j) - a(i,j)) + phi*a(i,j)
        dmn(i,j) = half*(delta(i,j) - a(i,j))

        ! check for stropholysis
        if (rotation_t) then
          sumg = zero
          do k = 1,3
            sum = zero
            do l = 1,3
              sum = sum + eps(k,l,i)*a(l,j) + eps(k,l,j)*a(l,i)
            end do
            sumg = sumg + half*vec_g(k)*sum
          end do
          rey(i,j) = rey(i,j) + sumg
        end if

        ! maintain constitutive relation
        cir(i,j) = delta(i,j) - rey(i,j) -dmn(i,j)
        
        ! symmetrize
        rey(j,i) = rey(i,j)
        dmn(j,i) = dmn(i,j)
        cir(j,i) = cir(i,j)
      end do
    end do

    ! compute flattened tensors
    if (chi /= chi ) then
      ierr = 8
      return
    end if

    if (abs(chi) < small) return
     
    if(abs(dot_vec_wdt) < small) then
      fl = zero
    else
      do i = 1,3
        do j = i,3
          fl(i,j) = vec_wdt(i)*vec_wdt(j)/dot_vec_wdt
          fl(j,i) = fl(i,j)
        end do
      end do
    end if

    afl = zero
    trace_afl = zero
    do i = 1,3
      do j = 1,3
        do k = 1,3
          afl(i,j) = afl(i,j) + a(i,k)*fl(k,j)
        end do
      end do

      trace_afl = trace_afl + afl(i,i)
    end do

    coef_1 = half*(one - chi*(one - trace_afl))
    coef_2 = -half*(one - chi*(one + trace_afl))
    coef_3 = -chi*(one - phi)
    coef_4 = (one - phi)*half*(one + chi*(one - trace_afl))
    coef_5 = (one - phi)*(-half*(one + chi*(one + trace_afl))) + phi

    do i = 1,3
      do j = i,3
        term = fl(i,j) - afl(i,j) - afl(j,i)
        
        ! dimensionality
        dmn(i,j) = coef_1*delta(i,j) + coef_2*a(i,j) + chi*term
        
        ! basic terms in Reynolds stress
        rey(i,j) = coef_4*delta(i,j) + coef_5*a(i,j) + coef_3*term
        
        ! check for stropholysis
        if (rotation_t) then
          sumg = zero
          do k = 1,3
            sumk = zero
            do l = 1,3
              sumi = zero
              sumj = zero
              do m = 1,3
                sumi = sumi + eps(m,l,i)*(coef_1*delta(k,m) + chi*fl(k,m) -     &
                                          chi*afl(k,m))
                sumj = sumj + eps(m,l,j)*(coef_1*delta(k,m) + chi*fl(k,m) -     &
                                          chi*afl(k,m))
              end do
              sumk = sumk + sumi*a(l,j) + sumj*a(l,i)
            end do
            sumg = sumg + vec_g(k)*sumk
          end do
          rey(i,j) = rey(i,j) + sumg
        end if

        ! circulicity
        cir(i,j) = delta(i,j) - rey(i,j) - dmn(i,j)

        ! symmetrize
        rey(j,i) = rey(i,j)
        dmn(j,i) = dmn(i,j)
        cir(j,i) = cir(i,j)
      end do
    end do
    
  end subroutine structure
        
    
