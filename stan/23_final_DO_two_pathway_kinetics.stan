// 23_final_DO_two_pathway_kinetics.stan
//
// PHYSICAL MODEL
// Surface rainfall DO starts at Csat (assumed atmospheric saturation).
// LSR loses DO between land surface and water table; source DO at WT
// is Csat * survival_LSR. River recharge DO at aquifer entry has a
// separate uncertain boundary Criver. Both undergo further aquifer
// DO loss. The two SOURCE PATHS are decayed SEPARATELY, THEN mixed
// at the sampling well; essential for nonlinear second-order decay.
//
// Aquifer kinetics (alternative model choices, not simultaneous):
// aq_order=1: dC/dt = -k_ref*C, k1_aq=k_ref [1/year].
// aq_order=2: dC/dt = -(k_ref/C_ref)*C^2,
//                  k2_aq=k_ref/C_ref [L/(mg year)].
// Both receive the SAME prior on instantaneous fractional rate
// at C_ref=7.5 mg/L; compare grouped held-out DO prediction.
//
// LSR kinetics (NOT identified independently without vadose time):
// survival=s=C_WT/Csat, s in (0,1).
// first-order exposure Da1=-ln(s) = k1_LSR*t_vadose;
// second-order exposure Da2=(1/s-1)=k2_LSR*Csat*t_vadose.
// Both are algebraically equivalent for a single common LSR DO entry
// concentration; data cannot choose the unsaturated reaction order.
// Conversion to conditional rates is done from posterior s draws in R,
// not by pretending t_vadose was observed.
//
// Source fractions jointly constrained by δ18O, Cl, Na. Correlated
// Na/Cl residuals and variable LSR δ18O endmember uncertainty.
// Excess air is OPTIONAL, default excluded (ablation result).
// DO is the only reactive response in the final primary model.
// Excess N2 / reaction progress not used as unproven DO constraints.
// Well 16758 excluded in the matching R runner.
//
// WARNING: reported groundwater MRT is the SAME effective aquifer
// residence time for both recharge sources; no source-specific ages.
// If MRT includes vadose travel, this partition is not identifiable.
// Csat and river_DO priors are assumptions, not measured surface data.
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
  int<lower=1,upper=2> aq_order;
  real<lower=0> Csat_rain;
  real<lower=0> C_ref;
  real<lower=0> endmember_prior_scale;
  real<lower=0> log_mrt_SD;
  real<lower=0> river_DO_prior_mean;
  real<lower=0> river_DO_prior_sd;
  real<lower=0> beta_survival_a;
  real<lower=0> beta_survival_b;
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
  real<lower=0,upper=1> survival_LSR;
  real<lower=0,upper=2> k_ref;
  real<lower=0> sigma_DO;
}
transformed parameters {
  real LSR_DO_at_water_table = Csat_rain*survival_LSR;
  vector[N] mean_DO;
  vector[N] aquifer_DO_loss;
  vector[N] vadose_LSR_DO_loss;
  vector[N] river_DO_after_aquifer;
  vector[N] LSR_DO_after_aquifer;
  for (i in 1:N) {
    real t = true_MRT[i];
    real river_at_well;
    real lsr_at_well;
    if (aq_order == 1) {
      river_at_well = river_DO*exp(-k_ref*t);
      lsr_at_well = LSR_DO_at_water_table*exp(-k_ref*t);
    } else {
      real k2 = k_ref/C_ref;
      river_at_well = river_DO/(1+k2*river_DO*t);
      lsr_at_well = LSR_DO_at_water_table/
                    (1+k2*LSR_DO_at_water_table*t);
    }
    river_DO_after_aquifer[i] = river_at_well;
    LSR_DO_after_aquifer[i] = lsr_at_well;
    mean_DO[i] = (1-f_LSR[i])*river_at_well+
                  f_LSR[i]*lsr_at_well;
    // Arithmetic per-sampled-water contribution decomposition.
    vadose_LSR_DO_loss[i] =
      f_LSR[i]*(Csat_rain-LSR_DO_at_water_table);
    aquifer_DO_loss[i] =
      (1-f_LSR[i])*(river_DO-river_at_well)+
      f_LSR[i]*(LSR_DO_at_water_table-lsr_at_well);
  }
}
model {
  O18_R ~ normal(-7.79,.24*endmember_prior_scale);
  O18_LSR ~ normal(-6.58,.40*endmember_prior_scale);
  Cl_R ~ normal(4.9,3*endmember_prior_scale);
  Cl_LSR ~ normal(60,15*endmember_prior_scale);
  Na_R ~ normal(6.6,3*endmember_prior_scale);
  Na_LSR ~ normal(44.3333333,15*endmember_prior_scale);
  air_R ~ normal(.1221,.35*endmember_prior_scale);
  air_LSR ~ normal(4.0434,1.4*endmember_prior_scale);
  f_LSR ~ beta(1,1);
  sigma_O18 ~ normal(0,.30);
  tau_LSR_O18 ~ normal(0,.80);
  sigma_Cl ~ normal(0,10);
  sigma_Na ~ normal(0,5);
  rho_Na_Cl ~ normal(0,.40);
  sigma_air ~ normal(0,1.2);

  river_DO ~ normal(river_DO_prior_mean,river_DO_prior_sd);
  survival_LSR ~ beta(beta_survival_a,beta_survival_b);
  k_ref ~ normal(0,.15); // common effective fractional-rate prior
  sigma_DO ~ normal(0,1.5);

  for (i in 1:N) {
    real f = f_LSR[i];
    real muO18 = (1-f)*O18_R+f*O18_LSR;
    real muCl = (1-f)*Cl_R+f*Cl_LSR;
    real muNa = (1-f)*Na_R+f*Na_LSR;
    real muAir = (1-f)*air_R+f*air_LSR;
    real sdO18 = sqrt(square(.1)+square(sigma_O18)+
                       square(f*tau_LSR_O18));
    real sdCl = sqrt(square(Cl_sd[i])+square(sigma_Cl));
    real sdNa = sqrt(square(Na_sd[i])+square(sigma_Na));
    real sdAir = sqrt(square(air_sd[i])+square(sigma_air));

    // Reported MRT is the mean of an assumed lognormal age prior.
    // Zero-MRT special case is a nonnegative near-zero-age prior.
    if (mrt_reported[i]>0)
      true_MRT[i] ~ lognormal(
        log(mrt_reported[i])-.5*square(log_mrt_SD),
        log_mrt_SD
      );
    else true_MRT[i] ~ normal(0,.5);

    if (has_O18[i]==1)
      O18[i] ~ normal(muO18,sdO18);
    if (has_air[i]==1)
      air[i] ~ normal(muAir,sdAir);

    // Exactly ONE bivariate-normal joint density for Na & Cl;
    // its conditional form avoids allocating covariance matrices.
    if (has_Cl[i]==1 && has_Na[i]==1) {
      Cl[i] ~ normal(muCl,sdCl);
      Na[i] ~ normal(
        muNa+rho_Na_Cl*sdNa/sdCl*(Cl[i]-muCl),
        sdNa*sqrt(1-square(rho_Na_Cl))
      );
    } else if (has_Cl[i]==1) {
      Cl[i] ~ normal(muCl,sdCl);
    } else if (has_Na[i]==1) {
      Na[i] ~ normal(muNa,sdNa);
    }

    // Robust DO observation to avoid a single anomalous well
    // dominating the common-rate kinetic inference.
    if (has_DO[i]==1)
      DO_obs[i] ~ student_t(
        4,mean_DO[i],sqrt(square(.2)+square(sigma_DO))
      );
  }
}
generated quantities {
  real k_aquifer_first_order_yr_inv = (aq_order==1) ? k_ref : 0;
  real k_aquifer_second_order_L_mg_yr = (aq_order==2)
                                       ? k_ref/C_ref : 0;
  real LSR_DO_deficit = Csat_rain-LSR_DO_at_water_table;
  real LSR_first_order_exposure = -log(survival_LSR);
  real LSR_second_order_exposure = 1/survival_LSR-1;
  real mean_vadose_LSR_DO_loss = mean(vadose_LSR_DO_loss);
  real mean_aquifer_DO_loss = mean(aquifer_DO_loss);
  real mean_f_LSR = mean(f_LSR);
  real sampled_well_vadose_loss_fraction =
    mean_vadose_LSR_DO_loss/
    (mean_vadose_LSR_DO_loss+mean_aquifer_DO_loss+1e-10);
  vector[N] replicated_DO;
  for (i in 1:N) {
    replicated_DO[i] = student_t_rng(
      4,mean_DO[i],sqrt(square(.2)+square(sigma_DO))
    );
  }
  // DO mass-balance attribution here is to nominal pathways, not
  // physical measurements of reaction location. The means above
  // weight every sampled well equally: NOT a catchment flux budget.
}
