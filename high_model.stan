data {
  int<lower=1> N;                   // 样本数
  int<lower=1> G;                   // 分组数
  int<lower=1, upper=G> group[N];   // 每个样本的组别

  vector[N] weight;
  vector[N] waist;

  vector[N] age;
  vector[N] bp_1d;
  vector[N] bp_1s;
  vector[N] chol;
  vector[N] stab_glu;
  vector[N] hdl;
  vector[N] height;
  vector[N] gender;

  vector[N] y;
}

parameters {
  // group-specific parameters
  vector[G] beta_weight;
  vector[G] beta_waist;
  vector[G] alpha;  // intercept for each group

  // shared coefficients
  real beta_age;
  real beta_bp_1d;
  real beta_bp_1s;
  real beta_chol;
  real beta_stab_glu;
  real beta_hdl;
  real beta_height;
  real beta_gender;

  real<lower=0> sigma;
}

model {
  // priors for group-specific
  beta_weight ~ normal(0, 2);
  beta_waist ~ normal(0, 2);
  alpha ~ normal(0, 5);

  // priors for shared
  beta_age ~ normal(0, 2);
  beta_bp_1d ~ normal(0, 2);
  beta_bp_1s ~ normal(0, 2);
  beta_chol ~ normal(0, 2);
  beta_stab_glu ~ normal(0, 2);
  beta_hdl ~ normal(0, 2);
  beta_height ~ normal(0, 2);
  beta_gender ~ normal(0, 2);
  sigma ~ cauchy(0, 2.5);

  // likelihood
  for (n in 1:N) {
    y[n] ~ normal(
      alpha[group[n]] +
      beta_weight[group[n]] * weight[n] +
      beta_waist[group[n]] * waist[n] +
      beta_age * age[n] +
      beta_bp_1d * bp_1d[n] +
      beta_bp_1s * bp_1s[n] +
      beta_chol * chol[n] +
      beta_stab_glu * stab_glu[n] +
      beta_hdl * hdl[n] +
      beta_height * height[n] +
      beta_gender * gender[n],
      sigma
    );
  }
}
