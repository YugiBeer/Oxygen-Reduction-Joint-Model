// Conservative-tracer-only Bayesian river/LSR water-fraction model.
// Same delta18O and chloride likelihood, parameter bounds, and priors
// as the tracer component of stan/04_river_DO_prior_sensitivity.stan.
// No DO, MRT, N2, nitrate or categorical source labels are used.
data {
  int<lower=1> N;
  array[N] int<lower=0,upper=1> has_O18;
  array[N] int<lower=0,upper=1> has_Cl;
  vector[N] O18;
  vector[N] Cl;
  vector<lower=0>[N] Cl_sd;
}
parameters {
  real<lower=-9,upper=-6.5> O18_R;
  real<lower=-7.8,upper=-5.4> O18_LSR;
  real<lower=0,upper=20> Cl_R;
  real<lower=20,upper=100> Cl_LSR;
  vector<lower=0,upper=1>[N] f_LSR;
  real<lower=0> sigma_O18;
  real<lower=0> tau_LSR_O18;
  real<lower=0> sigma_Cl;
}
model {
  O18_R ~ normal(-7.79,0.24);
  O18_LSR ~ normal(-6.58,0.40);
  Cl_R ~ normal(4.9,3);
  Cl_LSR ~ normal(60,15);
  f_LSR ~ beta(1,1);
  sigma_O18 ~ normal(0,0.30);
  tau_LSR_O18 ~ normal(0,0.80);
  sigma_Cl ~ normal(0,10);

  for (i in 1:N) {
    real f = f_LSR[i];
    real muO18 = (1-f)*O18_R + f*O18_LSR;
    real muCl = (1-f)*Cl_R + f*Cl_LSR;
    real sdO18 = sqrt(square(0.1)+square(sigma_O18)+
                      square(f*tau_LSR_O18));
    real sdCl = sqrt(square(Cl_sd[i])+square(sigma_Cl));
    if (has_O18[i]==1)
      O18[i] ~ normal(muO18,sdO18);
    if (has_Cl[i]==1)
      Cl[i] ~ normal(muCl,sdCl);
  }
}
