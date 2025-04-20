data {
  int<lower=0> N;
  vector[N] age;
  vector[N] weight;
  vector[N] waist;
  vector[N] bp_1d;
  vector[N] bp_1s;
  vector[N] y;
}

parameters {
  real beta_age;
  real beta_weight;
  real beta_waist;
  real beta_bp_1d;
  real beta_bp_1s;
  real<lower=0> sigma;
}

model {
  vector[N] mu;

  mu = beta_age * age +
       beta_weight * weight +
       beta_waist * waist +
       beta_bp_1d * bp_1d +
       beta_bp_1s * bp_1s;

  beta_age ~ normal(0, 5);
  beta_weight ~ normal(0, 5);
  beta_waist ~ normal(0, 5);
  beta_bp_1d ~ normal(0, 5);
  beta_bp_1s ~ normal(0, 5);
  sigma ~ cauchy(0, 2.5);

  y ~ normal(mu, sigma);
}
