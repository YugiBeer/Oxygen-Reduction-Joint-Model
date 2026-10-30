// 19_joint_three_vs_four_tracer_DO.stan
//
// Fully joint two-source Bayesian mixing + dissolved oxygen depletion.
// Identical 3/4 tracer control apart from optional Na observations.
// use_na=0 : delta18O, Cl, excess air; use_na=1 : plus sodium.
// Chloride and sodium have one JOINT bivariate normal likelihood,
// never treat Na/Cl as independent from chloride.
// Time/space variation in LSR isotope endmember is marginalized as
// extra (fLSR*tau_LSR_O18)^2 observation variance, as in Stan 18.
//
// The source DO boundary estimates are not observations of soil
// infiltration or of the saturated aquifer. k is age-associated only.
// Excess N2 and all other REDOX measurements are intentionally omitted.
// All well-specific DO is withheld by has_DO=0 during validation.
//
// SOURCE WARNING: endpoint priors are empirical from these same wells.
// Pre-infiltration DO reference=9.5 mg/L is ASSUMED, not measured.
data {
  int<lower=1> N;
  array[N] int<lower=0,upper=1> has_O18;
  array[N] int<lower=0,upper=1> has_Cl;
  array[N] int<lower=0,upper=1> has_Na;
  array[N] int<lower=0,upper=1> has_air;
  array[N] int<lower=0,upper=1> has_DO;
  vector[N] O18;
  vector[N] Cl;
  vector[N] Na;
  vector[N] air;
  vector[N] DO_obs;
  vector<lower=0>[N] Cl_sd;
  vector<lower=0>[N] Na_sd;
  vector<lower=0>[N] air_sd;
  vector<lower=0>[N] mrt_reported;
  int<lower=0,upper=1> use_na;
  real<lower=0> endmember_prior_scale;
  real<lower=0> log_mrt_SD;
  real<lower=0> lsr_do_prior_a;
  real<lower=0> lsr_do_prior_b;
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
  vector<lower=0>[N] true_MRT;
  real<lower=0> sigma_O18;
  real<lower=0> tau_LSR_O18;
  real<lower=0> sigma_Cl;
  real<lower=0> sigma_Na;
  real<lower=-0.95,upper=0.95> rho_Na_Cl;
  real<lower=0> sigma_air;
  real<lower=0,upper=12> river_DO;
  real<lower=0,upper=9.5> LSR_entry_DO;
  real<lower=0,upper=2> k;
  real<lower=0> sigma_DO;
}
transformed parameters {
  vector[N] mean_DO;
  vector[N] initial_DO;
  vector[N] age_associated_DO_loss;
  vector[N] LSR_pathway_DO_loss;
  for (i in 1:N) {
    initial_DO[i] = (1-f_LSR[i])*river_DO +
                     f_LSR[i]*LSR_entry_DO;
    mean_DO[i] = initial_DO[i]*exp(-k*true_MRT[i]);
    age_associated_DO_loss[i] = initial_DO[i]-mean_DO[i];
    LSR_pathway_DO_loss[i] = f_LSR[i]*(9.5-LSR_entry_DO);
  }
}
model {
  // Same tracer priors in both scenarios (from source dataset).
  O18_R ~ normal(-7.79, .24*endmember_prior_scale);
  O18_LSR ~ normal(-6.58, .4*endmember_prior_scale);
  Cl_R ~ normal(4.9, 3*endmember_prior_scale);
  Cl_LSR ~ normal(60, 15*endmember_prior_scale);
  Na_R ~ normal(6.6, 3*endmember_prior_scale);
  Na_LSR ~ normal(44.3333333, 15*endmember_prior_scale);
  air_R ~ normal(.1221,.35*endmember_prior_scale);
  air_LSR ~ normal(4.0434,1.4*endmember_prior_scale);
  f_LSR ~ beta(1,1);
  sigma_O18 ~ normal(0,.30);
  tau_LSR_O18 ~ normal(0,.80);
  sigma_Cl ~ normal(0,10);
  sigma_Na ~ normal(0,5);
  rho_Na_Cl ~ normal(0,.40);
  sigma_air ~ normal(0,1.2);
  river_DO ~ normal(7.5,.8);
  (LSR_entry_DO/9.5) ~ beta(lsr_do_prior_a,lsr_do_prior_b);
  k ~ normal(0,.15);
  sigma_DO ~ normal(0,1.5);

  for (i in 1:N) {
    real f = f_LSR[i];
    real muO18 = (1-f)*O18_R + f*O18_LSR;
    real muCl = (1-f)*Cl_R + f*Cl_LSR;
    real muNa = (1-f)*Na_R + f*Na_LSR;
    real muAir = (1-f)*air_R + f*air_LSR;
    real sdO18 = sqrt(square(.1)+square(sigma_O18)+
                      square(f*tau_LSR_O18));
    real sdCl = sqrt(square(Cl_sd[i])+square(sigma_Cl));
    real sdNa = sqrt(square(Na_sd[i])+square(sigma_Na));
    real sdAir = sqrt(square(air_sd[i])+square(sigma_air));

    if (mrt_reported[i] > 0) {
      // Reported MRT centres MEAN, not median, of latent lognormal.
      true_MRT[i] ~ lognormal(
        log(mrt_reported[i])-.5*square(log_mrt_SD),
        log_mrt_SD);
    } else {
      true_MRT[i] ~ normal(0,.5);
    }
    if (has_O18[i] == 1)
      O18[i] ~ normal(muO18,sdO18);
    if (has_air[i] == 1)
      air[i] ~ normal(muAir,sdAir);

    if (use_na == 1) {
      if (has_Cl[i] == 1 && has_Na[i] == 1) {
        Cl[i] ~ normal(muCl,sdCl);
        Na[i] ~ normal(
          muNa+rho_Na_Cl*sdNa/sdCl*(Cl[i]-muCl),
          sdNa*sqrt(1-square(rho_Na_Cl)));
      } else if (has_Cl[i] == 1) {
        Cl[i] ~ normal(muCl,sdCl);
      } else if (has_Na[i] == 1) {
        Na[i] ~ normal(muNa,sdNa);
      }
    } else {
      if (has_Cl[i] == 1)
        Cl[i] ~ normal(muCl,sdCl);
    }
    if (has_DO[i] == 1)
      DO_obs[i] ~ normal(
        mean_DO[i],sqrt(square(.2)+square(sigma_DO)));
  }
}
generated quantities {
  real LSR_DO_deficit = 9.5-LSR_entry_DO;
  real mean_well_LSR_pathway_loss = mean(LSR_pathway_DO_loss);
  real mean_well_age_associated_loss = mean(age_associated_DO_loss);
  real mean_well_LSR_fraction = mean(f_LSR);
  real sampled_well_fraction_LSR_of_loss =
    mean_well_LSR_pathway_loss /
    (mean_well_LSR_pathway_loss+
     mean_well_age_associated_loss+1e-9);
  // THIS IS A NON-VOLUME-WEIGHTED ARITHMETIC WELL SUMMARY.
  // Not a catchment recharge flux budget or a physically located
  // saturated-aquifer oxygen consumption fraction.
}
