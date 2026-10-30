// 12_joint_bayesian_mixing.stan
// Joint mixing + DO + uncertain MRT, with optional reactive excess N2.
//
// Mixing tracers: Cl, delta18O, EXCESS AIR (not excess N2).
// Excess N2 is a separate REACTIVE outcome; its relationship to oxygen
// consumption is an exploratory statistical association, not a
// stoichiometric denitrification model.
//
// The endmember priors use empirical values from the SAME study. This
// is empirical Bayes and can double-use source samples; do not call
// them independent endmember evidence. Test wider priors.
//
// k is an age-associated depletion coefficient, not independently
// identified as a saturated-aquifer reaction rate.

data {
  int<lower=1> N;
  array[N] int<lower=0,upper=1> has_o18;
  array[N] int<lower=0,upper=1> has_cl;
  array[N] int<lower=0,upper=1> has_air;
  array[N] int<lower=0,upper=1> has_n2;
  vector[N] o18;
  vector[N] cl;
  vector[N] air;
  vector[N] do_obs;
  vector[N] n2_obs;
  vector<lower=0>[N] cl_sd;
  vector<lower=0>[N] air_sd;
  vector<lower=0>[N] n2_sd;
  vector<lower=0>[N] mrt_reported;
  real<lower=0> log_mrt_SD;
  int<lower=0,upper=1> use_n2;
  int<lower=0,upper=1> use_gamma;
  int<lower=0,upper=1> vary_lsr_do;
  real<lower=0> gamma_prior_median;
  real<lower=0> tracer_prior_scale;
}
parameters {
  real<lower=-9,upper=-6.5> o18_R;
  real<lower=-7.8,upper=-5.4> o18_LSR;
  real<lower=0,upper=20> cl_R;
  real<lower=20,upper=100> cl_LSR;
  real<lower=0,upper=1.5> air_R;
  real<lower=1.5,upper=9> air_LSR;
  real<lower=0,upper=12> river_DO;
  real<lower=0,upper=9.5> central_LSR_entry_DO;
  real<lower=0> tau_LSR;
  vector[N] z_LSR;
  vector<lower=0,upper=1>[N] f_LSR;
  vector<lower=0>[N] true_MRT;
  real<lower=0,upper=2> k;
  real<lower=0> sigma_DO;
  real<lower=0> sigma_o18;
  real<lower=0> sigma_cl;
  real<lower=0> sigma_air;
  real<lower=0.25,upper=500> alpha;
  real<lower=0> n2_R;
  real<lower=0> n2_LSR;
  real<lower=0> beta_n2;
  real<lower=0> sigma_n2;
}
transformed parameters {
  vector[N] LSR_entry_DO;
  vector[N] mean_DO;
  vector[N] mean_n2;
  vector[N] apparent_age_DO_loss;
  for (i in 1:N) {
    real att;
    real initial_DO;
    if (vary_lsr_do == 1) {
      LSR_entry_DO[i] = 9.5 * inv_logit(
        logit(central_LSR_entry_DO / 9.5) + tau_LSR * z_LSR[i]
      );
    } else {
      LSR_entry_DO[i] = central_LSR_entry_DO;
    }
    initial_DO = (1 - f_LSR[i]) * river_DO +
                 f_LSR[i] * LSR_entry_DO[i];
    if (use_gamma == 1) {
      att = exp(-alpha * log1p(k * true_MRT[i] / alpha));
    } else {
      att = exp(-k * true_MRT[i]);
    }
    mean_DO[i] = initial_DO * att;
    apparent_age_DO_loss[i] = initial_DO - mean_DO[i];
    // Non-stoichiometric association. If N2 reacts through
    // different pathways, beta_n2 is NOT a causal oxygen:N2 yield.
    mean_n2[i] = (1-f_LSR[i]) * n2_R +
                 f_LSR[i] * n2_LSR +
                 beta_n2 * apparent_age_DO_loss[i];
  }
}
model {
  // Tracer endmember priors, loosely centred on study-derived values.
  o18_R ~ normal(-7.79, .24 * tracer_prior_scale);
  o18_LSR ~ normal(-6.58, .40 * tracer_prior_scale);
  cl_R ~ normal(4.9, 3 * tracer_prior_scale);
  cl_LSR ~ normal(60, 15 * tracer_prior_scale);
  air_R ~ normal(.122, .35 * tracer_prior_scale);
  air_LSR ~ normal(4.04, 1.4 * tracer_prior_scale);

  river_DO ~ normal(7.5, .8);
  // central_LSR_entry_DO uniform on 0 to 9.5, by its bounds;
  // no concentrated prior on its initial oxygen loss.
  tau_LSR ~ normal(0, 1.5);
  z_LSR ~ std_normal();
  f_LSR ~ beta(1, 1);

  k ~ normal(0, .15);
  sigma_DO ~ normal(0, 1.5);
  sigma_o18 ~ normal(0, .30);
  sigma_cl ~ normal(0, 10);
  sigma_air ~ normal(0, 1.2);

  // Log-normal prior on gamma shape alpha (including Jacobian).
  target += normal_lpdf(log(alpha) | log(gamma_prior_median), 1.3)
            - log(alpha);

  n2_R ~ normal(.8, .55);
  n2_LSR ~ normal(1.0, .8);
  beta_n2 ~ normal(0, .25);
  sigma_n2 ~ normal(0, 1);

  for (i in 1:N) {
    // Conditional prior for true MRT around each reported MRT.
    // The assumed log SD is sensitivity-tested, not observed.
    if (mrt_reported[i] > 0) {
      true_MRT[i] ~ lognormal(log(mrt_reported[i]), log_mrt_SD);
    } else {
      true_MRT[i] ~ normal(0, .5);
    }

    if (has_o18[i] == 1) {
      o18[i] ~ student_t(
        4,
        (1-f_LSR[i]) * o18_R + f_LSR[i] * o18_LSR,
        sqrt(square(.1) + square(sigma_o18))
      );
    }
    if (has_cl[i] == 1) {
      cl[i] ~ student_t(
        4,
        (1-f_LSR[i]) * cl_R + f_LSR[i] * cl_LSR,
        sqrt(square(cl_sd[i]) + square(sigma_cl))
      );
    }
    if (has_air[i] == 1) {
      air[i] ~ student_t(
        4,
        (1-f_LSR[i]) * air_R + f_LSR[i] * air_LSR,
        sqrt(square(air_sd[i]) + square(sigma_air))
      );
    }
    do_obs[i] ~ normal(mean_DO[i], sqrt(square(.2) + square(sigma_DO)));

    if (use_n2 == 1 && has_n2[i] == 1) {
      n2_obs[i] ~ normal(
        mean_n2[i],
        sqrt(square(n2_sd[i]) + square(sigma_n2))
      );
    }
  }
}
generated quantities {
  real central_LSR_DO_loss = 9.5 - central_LSR_entry_DO;
  real sample_mean_LSR_DO_loss = 9.5 - mean(LSR_entry_DO);
  vector[N] DO_residual = do_obs - mean_DO;
}
