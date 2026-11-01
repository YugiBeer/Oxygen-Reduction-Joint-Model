// 18_four_tracer_joint_mixing.stan
// Identical two-source Bayesian mixing baseline except for inclusion
// of sodium. 3-tracer control: use_na=0 (O18, Cl, excess air).
// 4-tracer model: use_na=1 (+Na, correlated residuals with Cl).
//
// No groundwater DO, MRT, redox chemistry or source labels in model.
// Delta18O LSR temporal/spatial variation is marginalized as a
// heteroscedastic, zero-mean source random effect f_LSR*tau_LSR_O18.
// This does not recover a time series of rainfall isotope composition.
//
// Sodium and chloride observed together have a BIVARIATE NORMAL
// likelihood with correlated conditional errors. If either ion is
// masked, use its correct one-dimensional marginal. Held-out ions are
// predicted CONDITIONALLY on the other measured ion.
//
// Study-derived endmember samples also feature in likelihood; the
// prior centres are empirical Bayes, NOT external independent data.
// Cl and Na are never duplicated as a ratio and separate ions.
//
data {
  int<lower=1> N;
  array[N] int<lower=0,upper=1> train_O18;
  array[N] int<lower=0,upper=1> train_Cl;
  array[N] int<lower=0,upper=1> train_Na;
  array[N] int<lower=0,upper=1> train_air;
  vector[N] O18;
  vector[N] Cl;
  vector[N] Na;
  vector[N] air;
  vector<lower=0>[N] Cl_sd;
  vector<lower=0>[N] Na_sd;
  vector<lower=0>[N] air_sd;
  int<lower=0,upper=1> use_na;
  real<lower=0> endpoint_prior_scale;
}
parameters {
  real<lower=-9,upper=-6.5> O18_R;
  real<lower=-7.8,upper=-5.4> O18_LSR;
  real<lower=0,upper=20> Cl_R;
  real<lower=20,upper=100> Cl_LSR;
  real<lower=0,upper=25> Na_R;
  real<lower=15,upper=100> Na_LSR;
  real<lower=0,upper=1.5> air_R;
  real<lower=1.5,upper=9> air_LSR;
  vector<lower=0,upper=1>[N] f_LSR;
  real<lower=0> sigma_O18;
  real<lower=0> tau_LSR_O18;
  real<lower=0> sigma_Cl;
  real<lower=0> sigma_Na;
  real<lower=-0.95,upper=0.95> rho_Na_Cl;
  real<lower=0> sigma_air;
}
model {
  // Empirical endmember priors, with explicitly adjustable widths.
  O18_R ~ normal(-7.79,0.24*endpoint_prior_scale);
  O18_LSR ~ normal(-6.58,0.40*endpoint_prior_scale);
  Cl_R ~ normal(4.9,3*endpoint_prior_scale);
  Cl_LSR ~ normal(60,15*endpoint_prior_scale);
  Na_R ~ normal(6.6,3*endpoint_prior_scale);
  Na_LSR ~ normal(44.3333333,15*endpoint_prior_scale);
  air_R ~ normal(0.1221,0.35*endpoint_prior_scale);
  air_LSR ~ normal(4.0434,1.4*endpoint_prior_scale);

  f_LSR ~ beta(1,1);
  sigma_O18 ~ normal(0,0.30);
  tau_LSR_O18 ~ normal(0,0.80);
  sigma_Cl ~ normal(0,10);
  sigma_Na ~ normal(0,5);
  rho_Na_Cl ~ normal(0,0.40);
  sigma_air ~ normal(0,1.2);

  for (i in 1:N) {
    real f = f_LSR[i];
    real muO18 = (1-f)*O18_R + f*O18_LSR;
    real muCl = (1-f)*Cl_R + f*Cl_LSR;
    real muNa = (1-f)*Na_R + f*Na_LSR;
    real muAir = (1-f)*air_R + f*air_LSR;
    real sdO18 = sqrt(square(0.1)+square(sigma_O18)+
                       square(f*tau_LSR_O18));
    real sdCl = sqrt(square(Cl_sd[i])+square(sigma_Cl));
    real sdNa = sqrt(square(Na_sd[i])+square(sigma_Na));
    real sdAir = sqrt(square(air_sd[i])+square(sigma_air));

    if (train_O18[i] == 1)
      O18[i] ~ normal(muO18,sdO18);
    if (train_air[i] == 1)
      air[i] ~ normal(muAir,sdAir);

    if (use_na == 1) {
      if (train_Cl[i] == 1 && train_Na[i] == 1) {
        // Exactly a bivariate normal joint density; the
        // conditionally specified form avoids matrix allocations.
        Cl[i] ~ normal(muCl,sdCl);
        Na[i] ~ normal(muNa +
                       rho_Na_Cl*sdNa/sdCl*(Cl[i]-muCl),
                       sdNa*sqrt(1-square(rho_Na_Cl)));
      } else if (train_Cl[i] == 1) {
        Cl[i] ~ normal(muCl,sdCl);
      } else if (train_Na[i] == 1) {
        Na[i] ~ normal(muNa,sdNa);
      }
    } else {
      if (train_Cl[i] == 1)
        Cl[i] ~ normal(muCl,sdCl);
    }
  }
}
generated quantities {
  vector[N] predicted_O18;
  vector[N] predicted_Cl;
  vector[N] predicted_Na;
  vector[N] predicted_air;
  vector[N] replicated_O18;
  vector[N] replicated_Cl;
  vector[N] replicated_Na;
  vector[N] replicated_air;
  real tau_LSR_O18_median = tau_LSR_O18;

  for (i in 1:N) {
    real f = f_LSR[i];
    real muO18 = (1-f)*O18_R + f*O18_LSR;
    real muCl = (1-f)*Cl_R + f*Cl_LSR;
    real muNa = (1-f)*Na_R + f*Na_LSR;
    real muAir = (1-f)*air_R + f*air_LSR;
    real sdO18 = sqrt(square(0.1)+square(sigma_O18)+
                       square(f*tau_LSR_O18));
    real sdCl = sqrt(square(Cl_sd[i])+square(sigma_Cl));
    real sdNa = sqrt(square(Na_sd[i])+square(sigma_Na));
    real sdAir = sqrt(square(air_sd[i])+square(sigma_air));
    real predictive_sdCl = sdCl;
    real predictive_sdNa = sdNa;

    predicted_O18[i] = muO18;
    predicted_Cl[i] = muCl;
    predicted_Na[i] = muNa;
    predicted_air[i] = muAir;

    // Conditional distribution of an unobserved ion, given the
    // observed other ion. Evaluated at each posterior draw.
    if (use_na == 1 && train_Cl[i] == 0 &&
        train_Na[i] == 1) {
      predicted_Cl[i] = muCl +
        rho_Na_Cl * sdCl/sdNa * (Na[i]-muNa);
      predictive_sdCl = sdCl*sqrt(1-square(rho_Na_Cl));
    }
    if (use_na == 1 && train_Na[i] == 0 &&
        train_Cl[i] == 1) {
      predicted_Na[i] = muNa +
        rho_Na_Cl * sdNa/sdCl * (Cl[i]-muCl);
      predictive_sdNa = sdNa*sqrt(1-square(rho_Na_Cl));
    }

    replicated_O18[i] = normal_rng(predicted_O18[i],sdO18);
    replicated_Cl[i] = normal_rng(predicted_Cl[i],
                                   predictive_sdCl);
    replicated_Na[i] = normal_rng(predicted_Na[i],
                                   predictive_sdNa);
    replicated_air[i] = normal_rng(predicted_air[i],sdAir);
  }
}
