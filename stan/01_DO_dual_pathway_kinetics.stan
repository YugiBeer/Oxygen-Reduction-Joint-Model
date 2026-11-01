// 01_DO_dual_pathway_kinetics.stan
// TWO SOURCE END MEMBERS: river recharge (R) and land-surface recharge (LSR).
// Mixing observations: stable-water isotope (delta18O) and chloride ONLY.
// Sodium, excess air, N2 and NO3 are NOT part of the model.
//
// DO is NOT used as a conservative mixing tracer. Rather, each endmember
// is depleted by the SAME aquifer kinetic law along the assumed aquifer MRT,
// and only then are the two remaining DO concentrations mixed.
//
// Boundary DO at rainfall surface = Csat_rain, fixed from assumed T and P.
// LSR entry DO at water table = Csat_rain * survival_LSR (unknown).
// River entry DO = river_DO_boundary, fixed at 7.5 mg/L by the R runner.
//
// aq_order 0: k0 = k_ref*C_ref [mg/L/year], max(0,Cin-k0*t)
// aq_order 1: k1 = k_ref       [1/year],   Cin*exp(-k1*t)
// Same k_ref prior for both orders; fit orders SEPARATELY.
//
// For positive reported MRT, true_MRT is a non-centred lognormal latent
// with prior arithmetic mean = reported MRT and log-SD = log_mrt_SD.
// Zero-MRT proxy 2801 is fixed at age zero; its DO observation is MASKED.
// No separate vadose residence-time model; the LSR water-table DO deficit
// is identified only conditional on the assumed mixing and aquifer ages.

data {
  int<lower=1> N;
  array[N] int<lower=0,upper=1> has_O18;
  array[N] int<lower=0,upper=1> has_Cl;
  array[N] int<lower=0,upper=1> has_DO;
  vector[N] O18;
  vector[N] Cl;
  vector[N] DO_obs;
  vector<lower=0>[N] Cl_sd;
  vector<lower=0>[N] mrt_reported;
  int<lower=0,upper=N> N_mrt_positive;
  array[N] int<lower=0,upper=N_mrt_positive> mrt_index;
  int<lower=0,upper=1> aq_order;
  real<lower=0> Csat_rain;
  real<lower=0> C_ref;
  real<lower=0> endmember_prior_scale;
  real<lower=0> log_mrt_SD;
  real<lower=0> river_DO_boundary;
  real<lower=0> beta_survival_a;
  real<lower=0> beta_survival_b;
}

parameters {
  real<lower=-9,upper=-6.5> O18_R;
  real<lower=-7.8,upper=-5.4> O18_LSR;
  real<lower=0,upper=20> Cl_R;
  real<lower=20,upper=100> Cl_LSR;
  vector<lower=0,upper=1>[N] f_LSR;
  vector[N_mrt_positive] z_mrt;
  real<lower=0> sigma_O18;
  real<lower=0> tau_LSR_O18;
  real<lower=0> sigma_Cl;
  real<lower=0,upper=1> survival_LSR;
  real<lower=0,upper=2> k_ref;
  real<lower=0> sigma_DO;
}

transformed parameters {
  vector[N] true_MRT;
  real river_DO = river_DO_boundary;
  real LSR_DO_at_water_table = Csat_rain*survival_LSR;
  vector[N] mean_DO;
  vector[N] aquifer_DO_loss;
  vector[N] vadose_LSR_DO_loss;
  vector[N] river_DO_after_aquifer;
  vector[N] LSR_DO_after_aquifer;

  for (i in 1:N) {
    if (mrt_index[i] == 0) {
      true_MRT[i] = 0;
    } else {
      true_MRT[i] = exp(
        log(mrt_reported[i]) - 0.5*square(log_mrt_SD) +
        log_mrt_SD*z_mrt[mrt_index[i]]
      );
    }
  }

  for (i in 1:N) {
    real t = true_MRT[i];
    real river_at_well;
    real lsr_at_well;

    if (aq_order == 0) {
      real k0 = k_ref*C_ref;
      river_at_well = fmax(0,river_DO-k0*t);
      lsr_at_well = fmax(0,LSR_DO_at_water_table-k0*t);
    } else {
      river_at_well = river_DO*exp(-k_ref*t);
      lsr_at_well = LSR_DO_at_water_table*exp(-k_ref*t);
    }

    river_DO_after_aquifer[i] = river_at_well;
    LSR_DO_after_aquifer[i] = lsr_at_well;
    mean_DO[i] =
      (1-f_LSR[i])*river_at_well+f_LSR[i]*lsr_at_well;

    // Nominal contributions per unit of sampled groundwater (not a
    // catchment-wide or discharge-weighted oxygen mass budget).
    vadose_LSR_DO_loss[i] =
      f_LSR[i]*(Csat_rain-LSR_DO_at_water_table);
    aquifer_DO_loss[i] =
      (1-f_LSR[i])*(river_DO-river_at_well)+
      f_LSR[i]*(LSR_DO_at_water_table-lsr_at_well);
  }
}

model {
  // Same active priors as the converged no-Na model 23a.
  O18_R ~ normal(-7.79,0.24*endmember_prior_scale);
  O18_LSR ~ normal(-6.58,0.40*endmember_prior_scale);
  Cl_R ~ normal(4.9,3*endmember_prior_scale);
  Cl_LSR ~ normal(60,15*endmember_prior_scale);
  f_LSR ~ beta(1,1);
  z_mrt ~ std_normal();

  sigma_O18 ~ normal(0,0.30);
  tau_LSR_O18 ~ normal(0,0.80);
  sigma_Cl ~ normal(0,10);
  survival_LSR ~ beta(beta_survival_a,beta_survival_b);
  k_ref ~ normal(0,0.15);
  sigma_DO ~ normal(0,1.5);

  for (i in 1:N) {
    real f = f_LSR[i];
    real muO18 = (1-f)*O18_R+f*O18_LSR;
    real muCl = (1-f)*Cl_R+f*Cl_LSR;
    real sdO18 = sqrt(
      square(0.1)+square(sigma_O18)+square(f*tau_LSR_O18)
    );
    real sdCl = sqrt(square(Cl_sd[i])+square(sigma_Cl));

    if (has_O18[i] == 1)
      O18[i] ~ normal(muO18,sdO18);
    if (has_Cl[i] == 1)
      Cl[i] ~ normal(muCl,sdCl);
    if (has_DO[i] == 1)
      DO_obs[i] ~ student_t(
        4,mean_DO[i],sqrt(square(0.2)+square(sigma_DO))
      );
  }
}

generated quantities {
  real aquifer_rate_actual = (aq_order==0) ? k_ref*C_ref : k_ref;
  real LSR_DO_deficit = Csat_rain-LSR_DO_at_water_table;
  real LSR_zero_order_exposure = LSR_DO_deficit;
  real LSR_first_order_exposure = -log(survival_LSR);
  real mean_vadose_LSR_DO_loss = mean(vadose_LSR_DO_loss);
  real mean_aquifer_DO_loss = mean(aquifer_DO_loss);
  real mean_f_LSR = mean(f_LSR);
  real sampled_well_vadose_loss_fraction =
    mean_vadose_LSR_DO_loss/
    (mean_vadose_LSR_DO_loss+mean_aquifer_DO_loss+1e-10);

  vector[N] replicated_DO;
  for (i in 1:N)
    replicated_DO[i] = student_t_rng(
      4,mean_DO[i],sqrt(square(0.2)+square(sigma_DO))
    );
}
