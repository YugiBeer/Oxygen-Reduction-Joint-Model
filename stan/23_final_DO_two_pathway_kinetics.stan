// 23_final_DO_two_pathway_kinetics.stan
//
// PHYSICAL MODEL
// Surface rainfall DO starts at Csat (assumed atmospheric saturation).
// LSR loses DO between land surface and water table; source DO at WT
// is Csat * survival_LSR. River recharge DO at aquifer entry has a
// separate uncertain boundary Criver. Both undergo further aquifer
// DO loss. The two SOURCE PATHS are decayed SEPARATELY, THEN mixed
// at the sampling well. Zero-order removal is clipped at zero.
//
// AQUFER kinetics: competing models, NOT fitted simultaneously:
// aq_order=0: dC/dt = -k0, C(t)=max(0,Cin-k0*t),
//             k0=k_ref*C_ref [mg/L/year].
// aq_order=1: dC/dt = -k1*C, C(t)=Cin*exp(-k1*t),
//             k1=k_ref [1/year].
// The two share the same weakly informative prior on INITIAL
// fractional depletion rate at a fixed reference concentration
// C_ref=7.5 mg/L. Compare held-out DO and posterior calibration.
//
// LSR attenuation (no independent vadose transit times):
// survival=s=C_WT/Csat, s in (0,1).
// zero-order exposure ΔC=Csat-C_WT=k0_LSR*t_vadose;
// first-order exposure -log(s)=k1_LSR*t_vadose.
// Both orders can reproduce the SAME water-table boundary DO.
// Conditional LSR rates in R require an explicitly ASSUMED t_vadose;
// DO data alone do NOT identify unsaturated-zone reaction order.
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
// Csat is assumed. The river DO is a FIXED 7.5-mg/L boundary
// supplied by R, not an inferred parameter or a second measurement
// taken from the (model-only) adjusted site 2801 value.
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
  int<lower=0,upper=1> aq_order;
  real<lower=0> Csat_rain;
  real<lower=0> C_ref;
  real<lower=0> endmember_prior_scale;
  real<lower=0> log_mrt_SD;
  real<lower=0> river_DO_boundary; // fixed 7.5 mg/L in default run
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
  real<lower=0,upper=1> survival_LSR;
  real<lower=0,upper=2> k_ref;
  real<lower=0> sigma_DO;
}
transformed parameters {
  real river_DO = river_DO_boundary; // fixed; for output/reporting
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
    if (aq_order == 0) {
      // Exact zero-order kinetics until depletion, never negative.
      real k0 = k_ref*C_ref; // mg/L/year
      river_at_well = fmax(0,river_DO-k0*t);
      lsr_at_well = fmax(0,LSR_DO_at_water_table-k0*t);
    } else {
      river_at_well = river_DO*exp(-k_ref*t);
      lsr_at_well = LSR_DO_at_water_table*exp(-k_ref*t);
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

  // River DO is a hard physical boundary in the supplied data.
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
  // One active kinetic rate per run. Units depend on aq_order:
  // zero = mg/L/year; first = year^-1.
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
  for (i in 1:N) {
    replicated_DO[i] = student_t_rng(
      4,mean_DO[i],sqrt(square(.2)+square(sigma_DO))
    );
  }
  // DO mass-balance attribution here is to nominal pathways, not
  // physical measurements of reaction location. The means above
  // weight every sampled well equally: NOT a catchment flux budget.
}
