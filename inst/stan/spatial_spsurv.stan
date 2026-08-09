/Users/sergiobejar/.rvm/scripts/rvm:29: operation not permitted: ps
functions {
  real event_lpdf(real t, real eta, real shape, int baseline) {
    if (baseline == 1) return log(shape) + eta + (shape - 1) * log(t) - exp(eta) * pow(t, shape);
    return log(shape) + eta + (shape - 1) * log(t) - 2 * log1p(exp(eta) * pow(t, shape));
  }
  real log_survival(real t, real eta, real shape, int baseline) {
    if (baseline == 1) return -exp(eta) * pow(t, shape);
    return -log1p(exp(eta) * pow(t, shape));
  }
}
data {
  int<lower=1> N; int<lower=1> M; int<lower=1> K; int<lower=1> Q;
  vector<lower=0>[N] t0; vector<lower=0>[N] time;
  array[N] int<lower=0,upper=1> event; matrix[N,K] X; matrix[M,Q] Z_spell;
  array[M] int<lower=1,upper=N> spell_start; array[M] int<lower=1,upper=N> spell_end;
  array[M] int<lower=0,upper=1> spell_event;
  int<lower=1> J; array[M] int<lower=1,upper=J> spell_group_id;
  int<lower=0> N_edges; array[N_edges] int<lower=1,upper=J> node1;
  array[N_edges] int<lower=1,upper=J> node2;
  int<lower=1> N_components; array[J] int<lower=1,upper=N_components> component_id;
  array[N_components] int<lower=1> component_size; vector<lower=0>[J] icar_scale;
  int<lower=1,upper=2> baseline;
}
parameters {
  vector[K] beta; vector[Q] gamma; real log_shape;
  vector[J] theta_v; vector[J] theta_w; vector[J] psi_v; vector[J] psi_w;
  real log_sigma_v; real log_sigma_w; real logit_phi_v; real logit_phi_w;
}
transformed parameters {
  vector[J] v; vector[J] w;
  real<lower=0> shape = exp(log_shape);
  real<lower=0> sigma_v = exp(log_sigma_v); real<lower=0> sigma_w = exp(log_sigma_w);
  real<lower=0,upper=1> phi_v = inv_logit(logit_phi_v);
  real<lower=0,upper=1> phi_w = inv_logit(logit_phi_w);
  for (j in 1:J) {
    real sv = component_size[component_id[j]] == 1 ? 0 : psi_v[j] / sqrt(icar_scale[j]);
    real sw = component_size[component_id[j]] == 1 ? 0 : psi_w[j] / sqrt(icar_scale[j]);
    v[j] = sigma_v * (sqrt(1 - phi_v) * theta_v[j] + sqrt(phi_v) * sv);
    w[j] = sigma_w * (sqrt(1 - phi_w) * theta_w[j] + sqrt(phi_w) * sw);
  }
}
model {
  beta ~ normal(0, 2); gamma ~ normal(0, 1.5); log_shape ~ normal(0, 0.5);
  theta_v ~ std_normal(); theta_w ~ std_normal();
  log_sigma_v ~ normal(-1, 0.7); log_sigma_w ~ normal(-1, 0.7);
  logit_phi_v ~ normal(-1, 1); logit_phi_w ~ normal(-1, 1);
  psi_v ~ normal(0, 1); psi_w ~ normal(0, 1);
  for (e in 1:N_edges) {
    target += -0.5 * square(psi_v[node1[e]] - psi_v[node2[e]]);
    target += -0.5 * square(psi_w[node1[e]] - psi_w[node2[e]]);
  }
  for (c in 1:N_components) if (component_size[c] > 1) {
    real sum_v = 0; real sum_w = 0;
    for (j in 1:J) if (component_id[j] == c) { sum_v += psi_v[j]; sum_w += psi_w[j]; }
    sum_v ~ normal(0, 0.001 * component_size[c]);
    sum_w ~ normal(0, 0.001 * component_size[c]);
  }
  for (m in 1:M) {
    int j = spell_group_id[m]; real pi = inv_logit(Z_spell[m] * gamma + v[j]); real ll = 0;
    for (n in spell_start[m]:spell_end[m]) {
      real eta = X[n] * beta + w[j];
      if (event[n] == 1) ll += event_lpdf(time[n] | eta, shape, baseline) - log_survival(t0[n], eta, shape, baseline);
      else ll += log_survival(time[n], eta, shape, baseline) - log_survival(t0[n], eta, shape, baseline);
    }
    if (spell_event[m] == 1) target += log(pi) + ll;
    else target += log_sum_exp(log1m(pi), log(pi) + ll);
  }
}
generated quantities {
  vector[M] log_lik;
  for (m in 1:M) {
    int j = spell_group_id[m]; real pi = inv_logit(Z_spell[m] * gamma + v[j]); real ll = 0;
    for (n in spell_start[m]:spell_end[m]) {
      real eta = X[n] * beta + w[j];
      if (event[n] == 1) ll += event_lpdf(time[n] | eta, shape, baseline) - log_survival(t0[n], eta, shape, baseline);
      else ll += log_survival(time[n], eta, shape, baseline) - log_survival(t0[n], eta, shape, baseline);
    }
    log_lik[m] = spell_event[m] == 1 ? log(pi) + ll : log_sum_exp(log1m(pi), log(pi) + ll);
  }
}
