// 22_joint_DO_denitrification_progress.stan
//
// TWO recharge sources, three baseline source tracers: delta18O,
// Na and Cl. Optional excess-air tracer. Na and Cl are observed
// jointly with correlated residuals. Variable rainfall isotope
// composition is marginalized (fLSR*tau_LSR_O18)^2.
//
// Paired outcomes:
// 1. Groundwater DO (mg/L): source DO at aquifer entry times
//    exp(-k*true_MRT), with reported MRT uncertainty.
// 2. Denitrification reaction progress in [0,1], derived as
//    excess_N2/(excess_N2+NO3_N+NO2_N).
//
// IMPORTANT: This is an empirically derived reaction-progress
// index, NOT an independent nitrate-mass-balance observation;
// neither excess N2 nor nitrate-derived fields enter as predictors.
// Exact progress zeros get a Bernoulli hurdle; positive progress is
// modeled robustly on logit scale (Student-t), with approximate
// propagated measurement uncertainty from excess-N2 analysis.
// No exact positive observations may equal 1 (must be handled as
// censored if they occur in a future dataset).
//
// Four alternatives (same mixing/DO model):
// use_progress=0 : DO only; process not informed by progress
// progress_link=0: progress ~ source + log1p(age)
// progress_link=1: + apparent cumulative DO depletion
// progress_link=2: + modelled low-DO exposure time
//
// The links are ASSOCIATIONS, not reaction stoichiometry.
// The DO trajectory is a conceptual piston-flow approximation.
// During paired CV has_DO=has_progress=0 at held-out wells.
//
// Modelled DO "LSR pathway loss" is based on assumed 9.5 mg/L
// recharge DO reference. Well-average losses are NOT flow weighted.
data {
  int<lower=1> N;
  array[N] int<lower=0,upper=1> has_O18;
  array[N] int<lower=0,upper=1> has_Cl;
  array[N] int<lower=0,upper=1> has_Na;
  array[N] int<lower=0,upper=1> has_air;
  array[N] int<lower=0,upper=1> has_DO;
  array[N] int<lower=0,upper=1> has_progress;
  vector[N] O18;
  vector[N] Cl;
  vector[N] Na;
  vector[N] air;
  vector[N] DO_obs;
  vector<lower=0,upper=1>[N] progress_obs;
  vector<lower=0>[N] Cl_sd;
  vector<lower=0>[N] Na_sd;
  vector<lower=0>[N] air_sd;
  vector<lower=0>[N] progress_logit_sd;
  vector<lower=0>[N] mrt_reported;
  int<lower=0,upper=1> use_na;
  int<lower=0,upper=1> use_progress;
  int<lower=0,upper=2> progress_link;
  real<lower=0> endmember_prior_scale;
  real<lower=0> log_mrt_SD;
  real<lower=0> lsr_do_prior_a;
  real<lower=0> lsr_do_prior_b;
  real<lower=0> low_DO_threshold;
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

  // A single small zero hurdle: two exact zeros among ~40 records.
  real<lower=0,upper=1> progress_zero_probability;
  real progress_intercept;
  real progress_source_beta;
  real progress_age_beta;
  real<lower=0> progress_link_beta;
  real<lower=0> sigma_progress_logit;
}
transformed parameters {
  vector[N] initial_DO;
  vector[N] mean_DO;
  vector[N] age_associated_DO_loss;
  vector[N] LSR_pathway_DO_loss;
  vector[N] low_DO_duration_years;
  vector[N] progress_logit_center;
  for (i in 1:N) {
    initial_DO[i] = (1-f_LSR[i])*river_DO+
                     f_LSR[i]*LSR_entry_DO;
    mean_DO[i] = initial_DO[i]*exp(-k*true_MRT[i]);
    age_associated_DO_loss[i] = initial_DO[i]-mean_DO[i];
    LSR_pathway_DO_loss[i] = f_LSR[i]*(9.5-LSR_entry_DO);
    if (k>1e-8) {
      real crossing = fmax(
        0,log(fmax(initial_DO[i],1e-10)/low_DO_threshold)/k
      );
      low_DO_duration_years[i] = fmax(0,true_MRT[i]-crossing);
    } else if (initial_DO[i]<=low_DO_threshold) {
      low_DO_duration_years[i] = true_MRT[i];
    } else {
      low_DO_duration_years[i] = 0;
    }
    progress_logit_center[i] =
      progress_intercept+
      progress_source_beta*(f_LSR[i]-.5)+
      progress_age_beta*log1p(true_MRT[i]);
    if (progress_link==1)
      progress_logit_center[i]+=
        progress_link_beta*age_associated_DO_loss[i]/5;
    if (progress_link==2)
      progress_logit_center[i]+=
        progress_link_beta*low_DO_duration_years[i]/20;
  }
}
model {
  // Identical source and DO priors for all candidate formulations.
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
  river_DO ~ normal(7.5,.8);
  (LSR_entry_DO/9.5) ~ beta(lsr_do_prior_a,lsr_do_prior_b);
  k ~ normal(0,.15);
  sigma_DO ~ normal(0,1.5);

  progress_zero_probability ~ beta(1.5,15);
  progress_intercept ~ normal(-.6,1.1);
  progress_source_beta ~ normal(0,.8);
  progress_age_beta ~ normal(0,.5);
  progress_link_beta ~ normal(0,1);
  sigma_progress_logit ~ normal(0,1.4);

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

    if (mrt_reported[i]>0) {
      // Reported MRT is the MEAN of this lognormal age prior.
      true_MRT[i] ~ lognormal(
        log(mrt_reported[i])-.5*square(log_mrt_SD),
        log_mrt_SD
      );
    } else true_MRT[i] ~ normal(0,.5);

    if (has_O18[i]==1) O18[i] ~ normal(muO18,sdO18);
    if (has_air[i]==1) air[i] ~ normal(muAir,sdAir);
    if (use_na==1) {
      if (has_Cl[i]==1 && has_Na[i]==1) {
        Cl[i] ~ normal(muCl,sdCl);
        Na[i] ~ normal(
          muNa+rho_Na_Cl*sdNa/sdCl*(Cl[i]-muCl),
          sdNa*sqrt(1-square(rho_Na_Cl))
        );
      } else if (has_Cl[i]==1) Cl[i] ~ normal(muCl,sdCl);
      else if (has_Na[i]==1) Na[i] ~ normal(muNa,sdNa);
    } else if (has_Cl[i]==1) Cl[i] ~ normal(muCl,sdCl);

    if (has_DO[i]==1)
      DO_obs[i] ~ normal(
        mean_DO[i],sqrt(square(.2)+square(sigma_DO))
      );

    if (use_progress==1 && has_progress[i]==1) {
      if (progress_obs[i]==0) {
        target += bernoulli_lpmf(1 | progress_zero_probability);
      } else {
        target += bernoulli_lpmf(0 | progress_zero_probability);
        // Transform of the FIXED observed fraction. Its Jacobian
        // does not depend on model parameters and cancels in
        // comparisons of posterior parameters.
        target += student_t_lpdf(
          logit(progress_obs[i]) | 4,
          progress_logit_center[i],
          sqrt(square(sigma_progress_logit)+
               square(progress_logit_sd[i]))
        );
      }
    }
  }
}
generated quantities {
  real LSR_DO_deficit = 9.5-LSR_entry_DO;
  real mean_well_LSR_pathway_loss = mean(LSR_pathway_DO_loss);
  real mean_well_age_associated_loss = mean(age_associated_DO_loss);
  real mean_well_LSR_fraction = mean(f_LSR);
  real sampled_well_fraction_LSR_of_loss =
    mean_well_LSR_pathway_loss/
    (mean_well_LSR_pathway_loss+
     mean_well_age_associated_loss+1e-9);
  vector[N] replicated_DO;
  vector[N] replicated_progress;
  for (i in 1:N) {
    replicated_DO[i] = normal_rng(
      mean_DO[i],sqrt(square(.2)+square(sigma_DO))
    );
    // When use_progress=0, these are prior-only simulations,
    // not a trained prediction and should NOT be scored.
    if (bernoulli_rng(progress_zero_probability)==1)
      replicated_progress[i] = 0;
    else
      replicated_progress[i] = inv_logit(student_t_rng(
        4,progress_logit_center[i],
        sqrt(square(sigma_progress_logit)+
             square(progress_logit_sd[i]))
      ));
  }
}
