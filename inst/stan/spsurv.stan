functions {
  real event_lpdf(real t, real eta, real shape, int baseline) {
    if (baseline == 1)
      return log(shape) + eta + (shape - 1) * log(t) - exp(eta) * pow(t, shape);
    return log(shape) + eta + (shape - 1) * log(t) - 2 * log1p(exp(eta) * pow(t, shape));
  }
  real log_survival(real t, real eta, real shape, int baseline) {
    if (baseline == 1) return -exp(eta) * pow(t, shape);
    return -log1p(exp(eta) * pow(t, shape));
  }
}
data {
  int<lower=1> N; int<lower=1> K; int<lower=1> Q;
  vector<lower=0>[N] t0; vector<lower=0>[N] time;
  array[N] int<lower=0,upper=1> event;
  matrix[N,K] X; matrix[N,Q] Z;
  int<lower=1> J; array[N] int<lower=1,upper=J> group_id;
  int<lower=0> N_edges; array[N_edges] int<lower=1,upper=J> node1;
  array[N_edges] int<lower=1,upper=J> node2;
  int<lower=1,upper=2> baseline; int<lower=0,upper=2> structure_type;
}
parameters {
  vector[K] beta; vector[Q] gamma; real<lower=0> shape;
  vector[J] v_raw; vector[J] w_raw;
  real<lower=0> sigma_v; real<lower=0> sigma_w;
}
transformed parameters {
  vector[J] v; vector[J] w;
  if (structure_type == 0) {
    v = rep_vector(0, J); w = rep_vector(0, J);
  } else {
    v = sigma_v * (v_raw - mean(v_raw));
    w = sigma_w * (w_raw - mean(w_raw));
  }
}
model {
  beta ~ normal(0, 2.5); gamma ~ normal(0, 2.5);
  shape ~ lognormal(0, 0.5);
  sigma_v ~ normal(0, 1); sigma_w ~ normal(0, 1);
  if (structure_type <= 1) {
    v_raw ~ std_normal(); w_raw ~ std_normal();
  }
  if (structure_type == 2) {
    for (e in 1:N_edges) {
      target += -0.5 * square(v_raw[node1[e]] - v_raw[node2[e]]);
      target += -0.5 * square(w_raw[node1[e]] - w_raw[node2[e]]);
    }
    sum(v_raw) ~ normal(0, 0.001 * J);
    sum(w_raw) ~ normal(0, 0.001 * J);
  }
  for (n in 1:N) {
    real pi = inv_logit(Z[n] * gamma + v[group_id[n]]);
    real eta = X[n] * beta + w[group_id[n]];
    real log_cond_surv = log_survival(time[n], eta, shape, baseline)
                         - log_survival(t0[n], eta, shape, baseline);
    if (event[n] == 1)
      target += log(pi) + event_lpdf(time[n] | eta, shape, baseline)
                - log_survival(t0[n], eta, shape, baseline);
    else
      target += log_sum_exp(log1m(pi), log(pi) + log_cond_surv);
  }
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    real pi = inv_logit(Z[n] * gamma + v[group_id[n]]);
    real eta = X[n] * beta + w[group_id[n]];
    real log_cond_surv = log_survival(time[n], eta, shape, baseline)
                         - log_survival(t0[n], eta, shape, baseline);
    if (event[n] == 1)
      log_lik[n] = log(pi) + event_lpdf(time[n] | eta, shape, baseline)
                   - log_survival(t0[n], eta, shape, baseline);
    else
      log_lik[n] = log_sum_exp(log1m(pi), log(pi) + log_cond_surv);
  }
}
