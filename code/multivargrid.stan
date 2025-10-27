data{
    array[2582] int n_nuts;
    vector[2582] log_follow_length;
    array[2582] int grid_id_follow;
    array[2582] int id;

    vector[15128] height_m;
    array[15128] int nutcracking;
    array[15128] int grid_id_follow_h;
    array[15128] int id_h;

}
parameters{
    matrix[6,17] z_id;
    matrix[6,40] z_grid;
    real bnp;
    real bnl;
    real alh;
    real aph;
    real aln;
    real apn;
    vector<lower=0>[6] sigma_id;
    vector<lower=0>[6] sigma_grid;
    cholesky_factor_corr[6] L_Rho_id;
    cholesky_factor_corr[6] L_Rho_grid;
    real<lower=0> scale;
}
transformed parameters{
    matrix[17,6] a_id;
    matrix[40,6] a_g110;
    a_g110 = (diag_pre_multiply(sigma_grid, L_Rho_grid) * z_grid)';
    a_id = (diag_pre_multiply(sigma_id, L_Rho_id) * z_id)';
}
model{
    vector[2582] lambda;
    vector[2582] p;
    vector[15128] mu;
    vector[15128] pp;

    scale ~ exponential( 1 );
    L_Rho_grid ~ lkj_corr_cholesky( 4 );
    L_Rho_id ~ lkj_corr_cholesky( 4 );
    sigma_grid ~ exponential( 1 );
    sigma_id ~ exponential( 1 );
    apn ~ normal( 0 , 1 );
    aln ~ normal( 0 , 1 );
    aph ~ normal( 0 , 1 );
    alh ~ normal( 0 , 1 );
    bnl ~ normal( 0 , 1 );
    bnp ~ normal( 0 , 1 );
    to_vector( z_grid ) ~ normal( 0 , 1 );
    to_vector( z_id ) ~ normal( 0 , 1 );
    
    for ( i in 1:2582 ){
        p[i] = apn + a_id[ id[i], 1] + a_g110[ grid_id_follow[i], 1];
        p[i] = inv_logit(p[i]);
        lambda[i] = aln + a_id[ id[i], 2] + a_g110[ grid_id_follow[i], 2] + log_follow_length[i];
        lambda[i] = exp(lambda[i]);
        if ( n_nuts[i]==0 ) target += log_mix( p[i] , 0 , poisson_lpmf(0|lambda[i]) );
        if ( n_nuts[i] > 0 ) target += log1m( p[i] ) + poisson_lpmf(n_nuts[i] | lambda[i] );
        }
    for ( j in 1:15128 ){
        pp[j] = aph + a_id[ id_h[j] , 3] + a_g110[ grid_id_follow_h[j], 3] +
        (bnp + a_id[ id_h[j], 4] + a_g110[grid_id_follow_h[j], 4]) * nutcracking[j];
        pp[j] = inv_logit(pp[j]);
        mu[j] = alh + a_id[ id_h[j], 5] + a_g110[grid_id_follow_h[j], 5] + 
        (bnl + a_id[ id_h[j], 6] + a_g110[grid_id_follow_h[j], 6]) * nutcracking[j];
        mu[j] = exp(mu[j]);
        if ( height_m[j] > 0 ) target += bernoulli_lpmf(0 | pp[j]) + gamma_lpdf(height_m[j] | mu[j]/scale, 1/scale);
        if ( height_m[j] == 0 ) target += bernoulli_lpmf(1 | pp[j]);
        }

}
generated quantities{
    matrix[6,6] Rho_id;
    matrix[6,6] Rho_grid;
    Rho_grid = multiply_lower_tri_self_transpose(L_Rho_grid);
    Rho_id = multiply_lower_tri_self_transpose(L_Rho_id);
}
