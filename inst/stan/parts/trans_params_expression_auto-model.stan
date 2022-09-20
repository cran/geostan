  if (car||sar) {
    log_lambda_mu = rep_vector(intercept, n);
    if (has_re) {
      for (i in 1:n) {
	log_lambda_mu[i] += alpha_re[id[i]];
      }
    }  
    if (dwx) {
	for (i in 1:dwx) log_lambda_mu += csr_matrix_times_vector(n, n, W_w, W_v, W_u, x_all[,wx_idx[i]]) * gamma[i];
      }
    if (dx_all) log_lambda_mu += x_all * beta;
    if (is_auto_gaussian) {
      fitted = offset + log_lambda_mu;
    } else {
      fitted = offset + log_lambda;
    }
  }
