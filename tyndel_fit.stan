data{
    array[125842] int time;
     vector[125842] dist;
     vector[125842] social_index_standardized;
     vector[125842] log_time;
    array[125842] int group_index;
    array[125842] int dyad_index;
    array[125842] int j_index;
    array[125842] int i_index;
}
parameters{
     real a_bar;
     real<lower=0> scale;
     real bT;
     real bS;
     real btxs;
     matrix[4,59] z_I;
     vector<lower=0>[4] sigma_I;
     cholesky_factor_corr[4] L_Rho_I;
     matrix[4,145] z_D;
     vector<lower=0>[4] sigma_D;
     cholesky_factor_corr[4] L_Rho_D;
     matrix[4,10] z_G;
     vector<lower=0>[4] sigma_G;
     cholesky_factor_corr[4] L_Rho_G;
}
transformed parameters{
     matrix[59,4] I;
     matrix[145,4] D;
     matrix[10,4] G;
    G = (diag_pre_multiply(sigma_G, L_Rho_G) * z_G)';
    D = (diag_pre_multiply(sigma_D, L_Rho_D) * z_D)';
    I = (diag_pre_multiply(sigma_I, L_Rho_I) * z_I)';
}
model{
     vector[125842] mu;
    L_Rho_G ~ lkj_corr_cholesky( 3 );
    sigma_G ~ exponential( 1 );
    to_vector( z_G ) ~ normal( 0 , 0.5 );
    L_Rho_D ~ lkj_corr_cholesky( 3 );
    sigma_D ~ exponential( 1 );
    to_vector( z_D ) ~ normal( 0 , 0.5 );
    L_Rho_I ~ lkj_corr_cholesky( 3 );
    sigma_I ~ exponential( 1 );
    to_vector( z_I ) ~ normal( 0 , 0.5 );
    btxs ~ normal( 0 , 1 );
    bS ~ normal( 0 , 1 );
    bT ~ normal( 0 , 1 );
    scale ~ exponential( 0.01 );
    a_bar ~ normal( 3.2 , 2 );
    for ( i in 1:125842 ) {
        mu[i] = a_bar + I[i_index[i], 1] + I[j_index[i], 1] + D[dyad_index[i], 1] + G[group_index[i], 1] + (bT + I[i_index[i], 2] + I[j_index[i], 2] + D[dyad_index[i], 2] + G[group_index[i], 2]) * log_time[i] + (bS + I[i_index[i], 3] + I[j_index[i], 3] + D[dyad_index[i], 3] + G[group_index[i], 3]) * social_index_standardized[i] + (btxs + I[i_index[i], 4] + I[j_index[i], 4] + D[dyad_index[i], 4] + G[group_index[i], 4]) * log_time[i] * social_index_standardized[i];
        mu[i] = exp(mu[i]);
    }
    dist ~ gamma( mu/scale , 1/scale );
}
generated quantities{
     matrix[4,4] Rho_I;
     matrix[4,4] Rho_D;
     matrix[4,4] Rho_G;
    Rho_G = multiply_lower_tri_self_transpose(L_Rho_G);
    Rho_D = multiply_lower_tri_self_transpose(L_Rho_D);
    Rho_I = multiply_lower_tri_self_transpose(L_Rho_I);
}
