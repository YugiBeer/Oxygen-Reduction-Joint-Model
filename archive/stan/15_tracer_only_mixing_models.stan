// 15_tracer_only_mixing_models.stan
//
// Three Bayesian tracer-only models; deliberately NO groundwater DO,
// MRT, nitrate or reactive excess N2 anywhere in likelihood/prior.
//
// A (variant 1): two endmembers, common tracer composition.
// B (variant 2): A + uncertain LSR isotope composition per recharge
//                location, integrated analytically as additional
//                f_LSR^2 * tau_O18^2 variance.
// C (variant 3): B + positive chloride enrichment not attributed to
//                LSR source mixing. Enrichment ~ HalfNormal(0,tau_Cl),
//                integrated analytically via skew-normal likelihood.
//
// All models share the SAME base endmember, fraction and error priors.
// Pure-endmember prior centres are inferred from this study, so this is
// conditional empirical Bayes rather than external validation.
//
// CV masks one measured tracer per site. Its observation is set to 0
// in the data passed to the model, and its likelihood is omitted.
// It is predicted from the remaining TWO conservative tracers.
//
// With Cl enrichment, additional chloride may reflect land use etc;
// NO process attribution is assumed from concentration data alone.
//
data {
  int<lower=1> N;
  array[N] int<lower=0,upper=1> train_O18;
  array[N] int<lower=0,upper=1> train_Cl;
  array[N] int<lower=0,upper=1> train_air;
  vector[N] O18;
  vector[N] Cl;
  vector[N] air;
  vector<lower=0>[N] Cl_sd;
  vector<lower=0>[N] air_sd;
  int<lower=1,upper=3> variant;
  real<lower=0> endmember_prior_scale;
}
parameters {
  real<lower=-9,upper=-6.5> O18_R;
  real<lower=-7.8,upper=-5.4> O18_LSR;
  real<lower=0,upper=20> Cl_R;
  real<lower=20,upper=100> Cl_LSR;
  real<lower=0,upper=1.5> air_R;
  real<lower=1.5,upper=9> air_LSR;
  vector<lower=0,upper=1>[N] f_LSR;

  // Study-level observation/process discrepancies, not analytical
  // measurement precision.
  real<lower=0> extra_O18_sd;
  real<lower=0> extra_Cl_sd;
  real<lower=0> extra_air_sd;

  // These nuisance priors are common across variants. The coefficients
  // enter the likelihood only when the relevant model is selected.
  real<lower=0> tau_O18_LSR;
  real<lower=0> tau_Cl_enrichment;
}
model {
  O18_R ~ normal(-7.79, 0.24 * endmember_prior_scale);
  O18_LSR ~ normal(-6.58, 0.40 * endmember_prior_scale);
  Cl_R ~ normal(4.9, 3 * endmember_prior_scale);
  Cl_LSR ~ normal(60, 15 * endmember_prior_scale);
  air_R ~ normal(0.122, 0.35 * endmember_prior_scale);
  air_LSR ~ normal(4.04, 1.4 * endmember_prior_scale);
  f_LSR ~ beta(1, 1);

  extra_O18_sd ~ normal(0, 0.30);
  extra_Cl_sd ~ normal(0, 10);
  extra_air_sd ~ normal(0, 1.2);
  tau_O18_LSR ~ normal(0, 0.8);   // permil
  tau_Cl_enrichment ~ normal(0, 12);  // mg/L

  for (i in 1:N) {
    real f = f_LSR[i];
    real mu_O18 = (1-f)*O18_R + f*O18_LSR;
    real mu_Cl = (1-f)*Cl_R + f*Cl_LSR;
    real mu_air = (1-f)*air_R + f*air_LSR;
    real base_o18 = sqrt(square(0.1) + square(extra_O18_sd));
    real base_cl = sqrt(square(Cl_sd[i]) + square(extra_Cl_sd));
    real base_air = sqrt(square(air_sd[i]) + square(extra_air_sd));

    if (train_O18[i] == 1) {
      real effective_o18 = base_o18;
      if (variant >= 2)
        effective_o18 = sqrt(
          square(base_o18) + square(f * tau_O18_LSR)
        );
      O18[i] ~ normal(mu_O18, effective_o18);
    }

    if (train_Cl[i] == 1) {
      if (variant == 3) {
        // sum of independent Normal and positive HalfNormal variates
        real omega = sqrt(
          square(base_cl) + square(tau_Cl_enrichment)
        );
        real shape = tau_Cl_enrichment / base_cl;
        Cl[i] ~ skew_normal(mu_Cl, omega, shape);
      } else
        Cl[i] ~ normal(mu_Cl, base_cl);
    }

    if (train_air[i] == 1)
      air[i] ~ normal(mu_air, base_air);
  }
}
generated quantities {
  vector[N] predicted_O18;
  vector[N] predicted_Cl;
  vector[N] predicted_air;
  vector[N] replicated_O18;
  vector[N] replicated_Cl;
  vector[N] replicated_air;

  for (i in 1:N) {
    real f = f_LSR[i];
    real sd_o18 = sqrt(square(0.1) + square(extra_O18_sd));
    real sd_cl = sqrt(square(Cl_sd[i]) + square(extra_Cl_sd));
    real sd_air = sqrt(square(air_sd[i]) + square(extra_air_sd));

    if (variant >= 2)
      sd_o18 = sqrt(square(sd_o18) +
                    square(f*tau_O18_LSR));

    predicted_O18[i] = (1-f)*O18_R + f*O18_LSR;
    predicted_Cl[i] = (1-f)*Cl_R + f*Cl_LSR;
    predicted_air[i] = (1-f)*air_R + f*air_LSR;

    // Mean of a HalfNormal(0,tau) positive chloride contribution.
    if (variant == 3)
      predicted_Cl[i] += tau_Cl_enrichment * sqrt(2/pi());

    replicated_O18[i] =
      normal_rng(predicted_O18[i], sd_o18);
    replicated_Cl[i] =
      normal_rng((1-f)*Cl_R + f*Cl_LSR, sd_cl);
    if (variant == 3)
      replicated_Cl[i] += fabs(normal_rng(0, tau_Cl_enrichment));
    replicated_air[i] =
      normal_rng(predicted_air[i], sd_air);
  }
}
