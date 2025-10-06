// IMAI MODEL 

data {
  
  int nA; // N age groups
  int nT; // N time points
  array[nT, nA] int cases; // reported case data
  matrix[nT, nA] pop; // population data
  array[2] vector[nA] ageLims; // lower & upper bounds of age groups - NOW NEED upper TO BE 1 higher
  // eg if age group is 0-4 yr olds need lower = 0 and upper = 5 ie "[0,5)" not "[0,4]"
  
}

transformed data {
  
  int sum_cases; // observed total cases per year
  array[nA] int cases_age;
  
  for (j in 1 : nA){
    cases_age[j] = sum(cases[, j]);
  } 
  
  sum_cases = sum(cases_age);
}


parameters {
  
  real<lower=0, upper=0.25> lam; // FOI
  real<lower=0, upper=1> rho; // reporting rate of 2nd infections
  real<lower=0, upper=1> gamma1; // relative reporting rate of 1st infections
  //real<lower=0, upper=1> gamma3; // relative reporting rate of 3rd and 4th infections
  
}



transformed parameters {
  
  vector<lower=0, upper=1>[nA] inc1; // incidence of primary infections
  vector<lower=0, upper=1>[nA] inc2; // incidence of secondary infections
  vector<lower=0, upper=1>[nA] inc3; // incidence of tertiary infections
  vector<lower=0, upper=1>[nA] inc4; // incidence of quaternary infections
  vector<lower=0, upper=1>[nA] D; // average annual incidence per person
  matrix[nT, nA] Ecases_py; // expected reported cases per year perage group 
  vector<lower=0>[nA] Ecases_tot_pa; // expected reported cases per age group
  real Ecases_tot; // expected reported cases in total
  vector<lower=0>[nA] Ecases_pa_prob; // expected probability of reported cases (per age group)
  real gamma3;
  
  gamma3 = 0;
  
    
  for(j in 1:nA){
  
  inc1[j] = exp(-4*lam*ageLims[1, j]) - exp(-4*lam*ageLims[2, j]);
  
  inc2[j] = 4*(exp(-3*lam*ageLims[1, j]) - exp(-3*lam*ageLims[2, j])) - 
    3*(exp(-4*lam*ageLims[1, j]) - exp(-4*lam*ageLims[2, j]));
  
  inc3[j] = 6*(exp(-2*lam*ageLims[1, j]) - exp(-2*lam*ageLims[2, j])) + 
    8*(exp(-3*lam*ageLims[2, j]) - exp(-3*lam*ageLims[1, j])) + 
    3*(exp(-4*lam*ageLims[1, j]) - exp(-4*lam*ageLims[2, j]));
  
  inc4[j] = 4*(exp(-lam*ageLims[1, j]) - exp(-lam*ageLims[2, j]) + 
                 exp(-3*lam*ageLims[1, j]) - exp(-3*lam*ageLims[2, j])) + 
    6*(exp(-2*lam*ageLims[2, j]) - exp(-2*lam*ageLims[1, j])) + 
    exp(-4*lam*ageLims[2, j]) - exp(-4*lam*ageLims[1, j]);
  
  }
  
  D = rho*(inc2 + gamma1*inc1 + gamma3*(inc3 + inc4)) ./ 
    (ageLims[2, ] - ageLims[1, ]);
  
  for(j in 1:nA){
  
  Ecases_py[, j] = D[j]*pop[, j];
  
  Ecases_tot_pa[j] = sum(Ecases_py[,j]);
  
  if (Ecases_tot_pa[j] == 0){
    Ecases_tot_pa[j] = 0.0001;
  } else{
    Ecases_tot_pa[j] = Ecases_tot_pa[j];
  }
  
  }
  
  Ecases_tot = sum(Ecases_tot_pa);
  
  if (Ecases_tot == 0){
    Ecases_tot = 0.0001;
  } else{
    Ecases_tot = Ecases_tot;
  }
  
  Ecases_pa_prob = Ecases_tot_pa / Ecases_tot;
  
  for (j in 1:nA) {
    
    if (Ecases_pa_prob[j] == 0){
      Ecases_pa_prob[j] = 0.0001;
    } else{
      Ecases_pa_prob[j] = Ecases_pa_prob[j];
    }
    
  }
  
}


model {
  
  //--- priors
  
  lam ~ normal(0, 1); 
  rho ~ normal(0, 1);
  gamma1 ~ normal(0, 1);
  //gamma3 ~ normal(0, 1);
  
  //--- likelihood on combined cases across all years!
    
  target += poisson_lpmf(sum_cases | Ecases_tot); //for total number of cases
  
  target += multinomial_lpmf(cases_age | Ecases_pa_prob); //for cases by age-groups
  
}



generated quantities {
  
  real log_lik;
  
  log_lik = poisson_lpmf(sum_cases | Ecases_tot) + multinomial_lpmf(cases_age | Ecases_pa_prob);

}
