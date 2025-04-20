data {
  int<lower=1> N;
  vector[N] age;
  vector[N] weight;
  vector[N] waist;
  vector[N] bp_1d;
  vector[N] bp_1s;
  vector[N] chol;
  vector[N] stab_glu;
  vector[N] hdl;
  vector[N] height;
  vector[N] y;
}

parameters {
  real beta_age;
  real beta_weight;
  real beta_waist;
  real beta_bp_1d;
  real beta_bp_1s;
  real beta_chol;
  real beta_stab_glu;
  real beta_hdl;
  real beta_height;
  real<lower=0> sigma;
}

model {
  beta_age ~ normal(0, 5);
  beta_weight ~ normal(0, 5);
  beta_waist ~ normal(0, 5);
  beta_bp_1d ~ normal(0, 5);
  beta_bp_1s ~ normal(0, 5);
  beta_chol ~ normal(0, 5);
  beta_stab_glu ~ normal(0, 5);
  beta_hdl ~ normal(0, 5);
  beta_height ~ normal(0, 5);
  sigma ~ cauchy(0, 2.5);

  y ~ normal(
    beta_age * age +
    beta_weight * weight +
    beta_waist * waist +
    beta_bp_1d * bp_1d +
    beta_bp_1s * bp_1s +
    beta_chol * chol +
    beta_stab_glu * stab_glu +
    beta_hdl * hdl +
    beta_height * height,
    sigma
  );
}
