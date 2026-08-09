/Users/sergiobejar/.rvm/scripts/rvm:29: operation not permitted: ps
functions {
  real event_lpdf(real t, real eta, real shape, int baseline) {
    if (baseline == 1)
      return log(shape) + eta + (shape - 1) * log(t) - exp(eta) * pow(t, shape);
    return log(shape) + eta + (shape - 1) * log(t)
           - 2 * log1p(exp(eta) * pow(t, shape));
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
  array[M] int<lower=1,upper=N> spell_start;
  array[M] int<lower=1,upper=N> spell_end;
  array[M] int<lower=0,upper=1> spell_event;
  int<lower=1> J; array[M] int<lower=1,upper=J> spell_group_id;
  int<lower=1,upper=2> baseline;
}
parameters {
  vector[K] beta; vector[Q] gamma; real log_shape;
  vector[J] z_v; vector[J] z_w; real log_sigma_v; real log_sigma_w;
}
transformed parameters {
  real<lower=0> shape = exp(log_shape);
  real<lower=0> sigma_v = exp(log_sigma_v);
  real<lower=0> sigma_w = exp(log_sigma_w);
  vector[J] v = sigma_v * (z_v - mean(z_v));
  vector[J] w = sigma_w * (z_w - mean(z_w));
}
model {
  beta ~ normal(0, 2); gamma ~ normal(0, 1.5); log_shape ~ normal(0, 0.5);
  z_v ~ std_normal(); z_w ~ std_normal();
  log_sigma_v ~ normal(-1, 0.7); log_sigma_w ~ normal(-1, 0.7);
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
