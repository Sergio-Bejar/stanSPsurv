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
  vector<lower=0>[N] time; array[N] int<lower=0,upper=1> event;
  matrix[N,K] X; matrix[N,Q] Z;
  int<lower=1> J; array[N] int<lower=1,upper=J> group_id;
  int<lower=0> N_edges; array[N_edges] int<lower=1,upper=J> node1;
  array[N_edges] int<lower=1,upper=J> node2;
  int<lower=1,upper=2> baseline; int<lower=0,upper=2> structure_type;
}
parameters {
  vector[K] beta; vector[Q] gamma; real<lower=0> shape;
  vector[J] u_raw; real<lower=0> sigma_u;
}
transformed parameters {
  vector[J] u;
  if (structure_type == 0) u = rep_vector(0, J);
  else u = sigma_u * (u_raw - mean(u_raw));
}
model {
  beta ~ normal(0, 2.5); gamma ~ normal(0, 2.5);
  shape ~ lognormal(0, 0.5); sigma_u ~ normal(0, 1);
  if (structure_type == 1) u_raw ~ std_normal();
  if (structure_type == 2) {
    for (e in 1:N_edges) target += -0.5 * square(u_raw[node1[e]] - u_raw[node2[e]]);
    sum(u_raw) ~ normal(0, 0.001 * J);
  }
  for (n in 1:N) {
    real pi = inv_logit(Z[n] * gamma);
    real eta = X[n] * beta + u[group_id[n]];
    if (event[n] == 1)
      target += log(pi) + event_lpdf(time[n] | eta, shape, baseline);
    else
      target += log_sum_exp(log1m(pi), log(pi) + log_survival(time[n], eta, shape, baseline));
  }
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    real pi = inv_logit(Z[n] * gamma);
    real eta = X[n] * beta + u[group_id[n]];
    if (event[n] == 1)
      log_lik[n] = log(pi) + event_lpdf(time[n] | eta, shape, baseline);
    else
      log_lik[n] = log_sum_exp(log1m(pi), log(pi) + log_survival(time[n], eta, shape, baseline));
  }
}

