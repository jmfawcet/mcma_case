data {
  int<lower=1> N;
  array[N] int<lower=0> y;
  array[N] int<lower=1> n;
  int<lower=1> M;
  array[N] int<lower=1, upper=M> measure;
  vector[M] prior_mu_se;
  vector[M] prior_mu_sp;
  vector<lower=0>[M] prior_sd_se;
  vector<lower=0>[M] prior_sd_sp;
  real prior_mu_pi;
  real<lower=0> prior_sd_pi;
  real prior_mu_z_rho;
  real<lower=0> prior_sd_z_rho;
}

parameters {
  real mu_pi;
  real<lower=0> tau_pi;
  vector[N] z_pi;
  vector[M] beta_se;
  vector[M] beta_sp;
  real<lower=0> tau_se;
  real<lower=0> tau_sp;
  real mu_z_rho;
  real<lower=0> tau_z_rho;
  vector[N] z_rho_raw;
  matrix[N, 2] z_se_sp;
}

transformed parameters {
  vector[N] logit_pi = mu_pi + tau_pi * z_pi;
  vector[N] z_rho_study = mu_z_rho + tau_z_rho * z_rho_raw;
  vector[N] rho_study = tanh(z_rho_study);
  vector[N] inner_se;
  vector[N] inner_sp;
  vector[N] p_apparent;

  for (i in 1:N) {
    matrix[2, 2] L_corr_i;
    L_corr_i[1, 1] = 1.0;
    L_corr_i[1, 2] = 0.0;
    L_corr_i[2, 1] = rho_study[i];
    L_corr_i[2, 2] = sqrt(1.0 - square(rho_study[i]));

    matrix[2, 2] L_Sigma_i;
    L_Sigma_i[1, 1] = tau_se * L_corr_i[1, 1];
    L_Sigma_i[1, 2] = tau_se * L_corr_i[1, 2];
    L_Sigma_i[2, 1] = tau_sp * L_corr_i[2, 1];
    L_Sigma_i[2, 2] = tau_sp * L_corr_i[2, 2];

    vector[2] u = L_Sigma_i * z_se_sp[i]';
    inner_se[i] = beta_se[measure[i]] + u[1];
    inner_sp[i] = beta_sp[measure[i]] + u[2];

    real pi_i = inv_logit(logit_pi[i]);
    real se_i = inv_logit(inner_se[i]);
    real sp_i = inv_logit(inner_sp[i]);
    p_apparent[i] = pi_i * se_i + (1 - pi_i) * (1 - sp_i);
  }
}

model {
  mu_pi ~ normal(prior_mu_pi, prior_sd_pi);
  tau_pi ~ normal(0, 0.5);
  z_pi ~ std_normal();
  beta_se ~ normal(prior_mu_se, prior_sd_se);
  beta_sp ~ normal(prior_mu_sp, prior_sd_sp);
  tau_se ~ normal(0, 0.5);
  tau_sp ~ normal(0, 0.5);
  mu_z_rho ~ normal(prior_mu_z_rho, prior_sd_z_rho);
  tau_z_rho ~ normal(0, 0.5);
  z_rho_raw ~ std_normal();
  to_vector(z_se_sp) ~ std_normal();
  y ~ binomial(n, p_apparent);
}

generated quantities {
  real pi_pop = inv_logit(mu_pi);
  real rho_pop = tanh(mu_z_rho);
  vector[M] Se_measure;
  vector[M] Sp_measure;
  for (m in 1:M) {
    Se_measure[m] = inv_logit(beta_se[m]);
    Sp_measure[m] = inv_logit(beta_sp[m]);
  }
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = binomial_lpmf(y[i] | n[i], p_apparent[i]);
  }
}
