//--- Time-varying dengue catalytic model ---//
// assumes constant endemic FOI prior to data
// assumes complete immunity after 2nd infection
// assumes equal transmissability of 4 serotypes

data {

  int nA; // N age groups
  int nT; // N time points
  int max_age; //max age of case data
  array[nT, nA, 4] int cases; // reported case data
  matrix[nT, nA] pop; // population data
  array[2, nA] int ageLims; // lower & upper bounds of age groups
  row_vector[max_age] age; // age as sequence 

}

transformed data {

  array[4, nT] int sum_cases; // observed total cases per year
  
  for (t in 1 : nT)
  for(s in 1 : 4){
    sum_cases[s, t] = sum(cases[t,  : , s]);
    }
}

parameters {

  simplex[4] lam_H; // historic average FOI
  array[nT] simplex[4] lam_t; // time varying FOI
  vector<lower=0, upper=1>[2] report;
  real<lower=0, upper = min(1./report)> chi; // relative reporting rate of children aged 1-15 yrs old
  real<lower=0, upper=1> sigma_H;
  array[nT] real<lower=0, upper=1> sigma_t;

}

transformed parameters {

  row_vector<lower=0, upper=1>[max_age] susc0; // proportion susceptible a time 0
  array[4] row_vector<lower=0, upper=1>[max_age] mono0; // proportion monotypic at time 0
  row_vector<lower=0, upper=1>[max_age] multi0; // proportion multitypic at time 0
  matrix<lower=0, upper=1>[nT, max_age] susc; // proportion susceptible
  array[4] matrix<lower=0, upper=1>[nT, max_age] mono; // proportion monotypic
  matrix<lower=0, upper=1>[nT, max_age] multi; // proportion multitypic
  array[4] matrix<lower=0, upper=1>[nT, max_age] inc1; // incidence of primary infections
  array[4] matrix<lower=0, upper=1>[nT, max_age] inc2; // incidence of secondary infections
  array[nT, 4] vector<lower=0>[nA] Ecases; // expected reported cases per year and age group
  array[4] vector<lower=0>[nT] Ecases_sum; // expected reported cases per year
  array[nT, 4] vector<lower=0>[nA] prob_Ecases; // expected probability of reported cases per year and age group

  //--- immune profiles at beginning of time series (time 0)
  susc0 = exp(-sum(sigma_H*lam_H) * age);

  for(s in 1:4){
  mono0[s, ] = exp(-(sum(sigma_H*lam_H) - sigma_H*lam_H[s]) * age) .* (1 - exp(-sigma_H*lam_H[s] * age));
  }
  
  for(a in 1:max_age){
  multi0[a] = 1 - susc0[a] - sum(mono0[, a]);
  }

  //--- infants (assumes no infections in <1 year olds), we have observed cases also for below 1 year

  // for(t in 1:nT) susc[t,1] = 1;
  // for(t in 1:nT) mono[t,1] = 0;
  // for(t in 1:nT) multi[t,1] = 0;
  //--- immune profiles at time 1 
  
  susc[1, 1] = 1;
  
  multi[1, 1] = 0;
  
  susc[1, 2:max_age] = susc0[1:(max_age-1)] - sum(sigma_t[1]*lam_t[1]) * susc0[1:(max_age-1)];
  
  
  for(s in 1:4){
  mono[s, 1, 1] = 0;

  mono[s, 1, 2:max_age] = mono0[s, 1:(max_age-1)] + sigma_t[1]*lam_t[1][s] * susc0[1:(max_age-1)]
                     - (sum(sigma_t[1]*lam_t[1]) - sigma_t[1]*lam_t[1][s]) * mono0[s, 1:(max_age-1)];

  inc1[s, 1, ] = sigma_t[1]*lam_t[1][s] * susc0;
  
  for(a in 1:max_age){
  inc2[s, 1, a] = sigma_t[1]*lam_t[1][s] * (sum(mono0[, a]) - mono0[s,a]);
  
  }}

  for(a in 2:max_age){
    multi[1, a] = multi0[(a-1)] + sum(inc2[, 1, (a-1)]);
  }
  
  
  
  

  //--- immune profiles at subsequent time steps 

if(nT >= 2){
  for (t in 2 : nT) {
    susc[t, 1] = 1;

    multi[t, 1] = 0;
    
    susc[t, 2:max_age] = susc[t - 1, 1:(max_age-1)]
                       - sum(sigma_t[t]*lam_t[t]) * susc[t - 1, 1:(max_age-1)];
                       
    for(s in 1:4){    
    mono[s, t, 1] = 0;

    mono[s, t, 2:max_age] = mono[s, t - 1, 1:(max_age-1)]
                       + sigma_t[t]*lam_t[t][s] * susc[t - 1, 1:(max_age-1)]
                       - (sum(sigma_t[t]*lam_t[t]) - sigma_t[t]*lam_t[t][s]) * mono[s, t - 1, 1:(max_age-1)];

    inc1[s, t, ] = sigma_t[t]*lam_t[t][s] * susc[t - 1, ];
    
    for(a in 1:max_age){
    inc2[s, t, a] = sigma_t[t]*lam_t[t][s] * (sum(mono[, t - 1, a]) - mono[s, t - 1, a]);
    }}
                       
    for(a in 2:max_age){
    multi[t, a] = multi[t - 1, (a - 1)] + sum(inc2[, t, (a - 1)]);
    }                   

  }
}

  //--- expected reported cases for each year and age group

  for (t in 1 : nT) 

    for (a in 1 : nA) {
      
      for(s in 1 : 4){
        
        if (a <= 3){ 

      Ecases[t, s, a] = chi*(report[1]
                     * mean(inc2[s, t, ageLims[1, a] : ageLims[2, a]])
                        + report[2]
                          * mean(inc1[s, t, ageLims[1, a] : ageLims[2, a]]))
                     * pop[t, a];
    }
    else{
      
       Ecases[t, s, a] = (report[1]
                     * mean(inc2[s, t, ageLims[1, a] : ageLims[2, a]])
                        + report[2]
                          * mean(inc1[s, t, ageLims[1, a] : ageLims[2, a]]))
                     * pop[t, a];
    }
      }
    }

  // Apply your conditional logic to ensure that the log likelihood never works with 0 but close (0.0001)

  for (t in 1 : nT) 

    for (a in 1 : nA) {
      
      for(s in 1:4){

      if (Ecases[t, s, a] == 0) {
        Ecases[t, s, a] = 0.0001;

      } else {
        Ecases[t, s, a] = Ecases[t, s, a];
      }

    }
    }

  //---  expected reported cases total per year

  for (t in 1 : nT) {
    for(s in 1:4){
    Ecases_sum[s, t] = sum(Ecases[t, s, : ]);
  }
  }

  //--- Expected probability of reported cases per year and age group

  for (t in 1 : nT) {
    for(s in 1:4){
    prob_Ecases[t, s, : ] = Ecases[t, s, : ] / sum(Ecases[t, s, : ]);
  }
  }

}

model {

  //--- priors
  row_vector[4] alpha = [1, 1, 1, 1];
  for (t in 1:nT) {
    lam_t[t] ~ dirichlet(alpha);
  }
  lam_H ~ dirichlet(alpha);
  report ~ normal(0, 1); 

  sigma_H ~ uniform(0,1);
  sigma_t ~ uniform(0,1);
  chi ~ normal(1, 1); 
  //--- likelihood 
  // poisson likelihood for total per year

  for(s in 1 : 4) for (t in 1 : nT) target += poisson_lpmf(sum_cases[s, t] | Ecases_sum[s, t]);

  // multinomial likelihood for age and year

  for(s in 1 : 4) for (t in 1 : nT) target += multinomial_lpmf(cases[t,  : , s] | prob_Ecases[t, s, : ]);

}

generated quantities {
 
  array[4] vector[nT] log_lik;
  array[4] vector[nT] log_lik1;
  array[4] vector[nT] log_lik2;
 
  for(s in 1 : 4) for (t in 1 : nT)     log_lik1[s, t] = poisson_lpmf(sum_cases[s, t] | Ecases_sum[s, t]);
 
  for(s in 1 : 4) for (t in 1 : nT)     log_lik2[s, t] = multinomial_lpmf(cases[t,  : , s] | prob_Ecases[t, s, : ]);
 
  for(s in 1 : 4)  for (t in 1 : nT)     log_lik[s, t] = log_lik1[s, t] + log_lik2[s, t];
 
}


