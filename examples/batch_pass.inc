  ! ── Batch pass: iterate to converge ──
  b_rms_prev = 1.0d10
  do b_iter = 1, b_n_phase_iter
    b_total_iter = b_total_iter + 1

    x_base_b = x_iter
    jd_current = jd_start
    b_m_total = 0

    do b_ii = 1, b_n_obs_win
      i = batch_obs_idx(b_ii)
      dt_sec = (obs(i)%jd - jd_current) * DAY_S
      if (abs(dt_sec) > 0.1_dp) then
        call propagate_state(x_base_b, dt_sec, x_base_b)
      end if
      jd_current = obs(i)%jd

      call predict_altaz(x_base_b, obs(i)%jd, lat_deg, lon_deg, obs(i)%body, &
                         b_pred_alt, b_pred_az)

      b_m_total = b_m_total + 1
      batch_residuals(1, b_m_total) = obs(i)%alt_obs - b_pred_alt
      batch_residuals(2, b_m_total) = obs(i)%az_obs  - b_pred_az
      if (batch_residuals(2, b_m_total) > 180.0_dp)  batch_residuals(2, b_m_total) = batch_residuals(2, b_m_total) - 360.0_dp
      if (batch_residuals(2, b_m_total) < -180.0_dp) batch_residuals(2, b_m_total) = batch_residuals(2, b_m_total) + 360.0_dp
    end do

    b_rms_alt = 0.0_dp; b_rms_az = 0.0_dp
    do b_ii = 1, b_m_total
      b_rms_alt = b_rms_alt + batch_residuals(1, b_ii)**2
      b_rms_az  = b_rms_az  + batch_residuals(2, b_ii)**2
    end do
    b_rms_alt = sqrt(b_rms_alt / max(real(b_m_total, dp), 1.0_dp)) * 3600.0_dp
    b_rms_az  = sqrt(b_rms_az  / max(real(b_m_total, dp), 1.0_dp)) * 3600.0_dp

    print '(A,I3,A,F10.1,A,F10.1,A,3F9.4)', &
         '  it', b_total_iter, ' RMS:', b_rms_alt, '"', b_rms_az, &
         '" GM%:', &
         (x_iter(19)/GM_SUN_TRUE - 1.0_dp)*100.0_dp, &
         (x_iter(20)/GM_EARTH_TRUE - 1.0_dp)*100.0_dp, &
         (x_iter(21)/GM_MOON_TRUE - 1.0_dp)*100.0_dp

    if (b_rms_alt < 22.0_dp .and. b_rms_az < 35.0_dp) then
      print '(A)', '  *** Converged ***'
      exit
    end if
    if (b_iter > 4 .and. abs(b_rms_alt + b_rms_az - b_rms_prev) / max(b_rms_prev, 1.0_dp) < 5.0d-5) then
      print '(A)', '  (stalled)'
      exit
    end if
    b_rms_prev = b_rms_alt + b_rms_az

    ! Jacobian — only perturb active parameters
    batch_H = 0.0_dp
    do b_kk = 1, b_n_active
      b_col = b_active_idx(b_kk)
      if (b_col <= 18) then
        if (mod(b_col-1, 6) < 3) then
          b_delta_param = 1.0_dp
        else
          b_delta_param = 1.0d-5
        end if
      else
        b_delta_param = abs(x_iter(b_col)) * 1.0d-7
        if (b_delta_param < 0.1_dp) b_delta_param = 0.1_dp
      end if

      x_pert_b = x_iter
      x_pert_b(b_col) = x_pert_b(b_col) + b_delta_param

      x_base_b = x_pert_b
      jd_current = jd_start
      do b_ii = 1, b_n_obs_win
        i = batch_obs_idx(b_ii)
        dt_sec = (obs(i)%jd - jd_current) * DAY_S
        if (abs(dt_sec) > 0.1_dp) then
          call propagate_state(x_base_b, dt_sec, x_base_b)
        end if
        jd_current = obs(i)%jd

        call predict_altaz(x_base_b, obs(i)%jd, lat_deg, lon_deg, obs(i)%body, &
                           b_pred_alt_p, b_pred_az_p)

        b_row = (b_ii - 1) * NM_BATCH
        batch_H(b_row+1, b_col) = (b_pred_alt_p - (obs(batch_obs_idx(b_ii))%alt_obs - batch_residuals(1,b_ii))) / b_delta_param
        batch_H(b_row+2, b_col) = (b_pred_az_p  - (obs(batch_obs_idx(b_ii))%az_obs  - batch_residuals(2,b_ii))) / b_delta_param
      end do
    end do

    ! Normal equations
    batch_HtWH = 0.0_dp
    batch_HtWr = 0.0_dp
    do b_ii = 1, b_m_total
      b_row = (b_ii - 1) * NM_BATCH
      i = batch_obs_idx(b_ii)
      if (obs(i)%body == 1) then
        b_w_alt = 1.0_dp / (10.0_dp / 3600.0_dp)**2   ! Sun: 10"
        b_w_az  = 1.0_dp / (10.0_dp / 3600.0_dp)**2
      else
        b_w_alt = 1.0_dp / (120.0_dp / 3600.0_dp)**2  ! Moon: 120"
        b_w_az  = 1.0_dp / (120.0_dp / 3600.0_dp)**2
      end if

      do b_jj = 1, NX
        batch_HtWr(b_jj) = batch_HtWr(b_jj) &
             + batch_H(b_row+1, b_jj) * b_w_alt * batch_residuals(1, b_ii) &
             + batch_H(b_row+2, b_jj) * b_w_az  * batch_residuals(2, b_ii)
        do b_kk = 1, NX
          batch_HtWH(b_jj,b_kk) = batch_HtWH(b_jj,b_kk) &
               + batch_H(b_row+1, b_jj) * b_w_alt * batch_H(b_row+1, b_kk) &
               + batch_H(b_row+2, b_jj) * b_w_az  * batch_H(b_row+2, b_kk)
        end do
      end do
    end do

    ! LM damping on active params
    do b_kk = 1, b_n_active
      b_jj = b_active_idx(b_kk)
      batch_HtWH(b_jj,b_jj) = batch_HtWH(b_jj,b_jj) * (1.0_dp + 50.0_dp / real(b_iter, dp))
    end do
    ! Freeze inactive parameters
    do b_jj = 1, NX
      if (.not. b_is_active(b_jj)) then
        batch_HtWH(b_jj,b_jj) = 1.0_dp
        batch_HtWr(b_jj) = 0.0_dp
      end if
    end do

    ! GM prior
    do b_jj = 19, 21
      if (b_sigma_prior(b_jj) > 0.0_dp .and. b_is_active(b_jj)) then
        batch_HtWH(b_jj,b_jj) = batch_HtWH(b_jj,b_jj) + 1.0_dp / b_sigma_prior(b_jj)**2
        batch_HtWr(b_jj) = batch_HtWr(b_jj) + (x_prior(b_jj) - x_iter(b_jj)) / b_sigma_prior(b_jj)**2
      end if
    end do

    call solve_symmetric(batch_HtWH, batch_HtWr, NX, batch_dx)

    ! Step limits
    do b_jj = 1, NX
      if (b_jj <= 18) then
        if (mod(b_jj-1,6) < 3) then
          batch_dx(b_jj) = sign(min(abs(batch_dx(b_jj)), b_pos_lim), batch_dx(b_jj))
        else
          batch_dx(b_jj) = sign(min(abs(batch_dx(b_jj)), b_vel_lim), batch_dx(b_jj))
        end if
      else
        batch_dx(b_jj) = sign(min(abs(batch_dx(b_jj)), b_gm_lim * abs(x_iter(b_jj))), batch_dx(b_jj))
      end if
    end do

    x_iter = x_iter + batch_dx
    if (x_iter(19) < 1.0d6)  x_iter(19) = 1.0d6
    if (x_iter(20) < 1.0d2)  x_iter(20) = 1.0d2
    if (x_iter(21) < 1.0d0)  x_iter(21) = 1.0d0
  end do
