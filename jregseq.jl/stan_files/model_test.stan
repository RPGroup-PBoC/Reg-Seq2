functions{
    // Linear energy model
     vector energy(vector theta, matrix seq, int N){
      vector[N] E;
      E = seq * theta;
      return E;
    }
    // Kernel density estimation
    real K_E(real E_i, vector E, int N_S, real h, int[] ct, int n){
        vector[N_S] K;
        for (i in 1:N_S){
            K[i] = ct[i]/(h * n) * exp(normal_lpdf((E_i - E[i])/h| 0, 1));
        }
    return sum(K);
    }
    // Mutual information
    real I(vector p_0, vector p_1, int N_E){
        
        real p_mu0 = sum(p_0);
        real p_mu1 = sum(p_1);
 
        vector[2 * N_E] MI;
        for (i in 1:N_E){
            MI[i] = p_0[i] * log2(p_0[i]/((p_0[i] + p_1[i]) * p_mu0));
            MI[i + N_E] = p_1[i] * log2(p_1[i]/((p_0[i] + p_1[i]) * p_mu1));
        }
    return sum(MI);
    }
    
    vector gauge_energy(matrix centering_matrix, vector E, int L_S){
      vector[4*L_S] E_new;
      // centering
      E_new = centering_matrix * E;
      // Sum of column variances equal to length
      E_new = E_new / sqrt(dot_product(E_new',E_new) / 3) * sqrt(L_S);
    return E_new;
    
    }
    
}



data {
  int L_S;
  int N_S;
  matrix[N_S, 4 * L_S] seqs;
  int ct_0[N_S];
  int ct_1[N_S];
  int n;
  matrix[4 * L_S, 4 * L_S] centering_matrix;
}



generated quantities{
    // Model parameters
    vector[L_S*4] theta;
    for (i in 1:(4*L_S)){
        theta[i] = normal_rng(0, 1);
        
    }
    
    // DNA energies
    vector[N_S] E_arr = energy(theta, seqs, L_S);
    theta = gauge_energy(centering_matrix, theta, L_S);
    
    vector[500] p_E_0;
    vector[500] p_E_1;
    real E_min = min(E_arr);
    real E_max = max(E_arr);
    real E_i;
    real d_E = 1/499. * (E_max-E_min);
    for (i in 0:499){
        E_i = E_min + i * d_E;
        p_E_0[i+1] = K_E(E_i, E_arr, N_S, 0.04, ct_0, n);
        p_E_1[i+1] = K_E(E_i, E_arr, N_S, 0.04, ct_1, n);
    }
    p_E_0 += 1e-9;
    p_E_1 += 1e-9;
    real p_mu0 = sum(p_E_0);
    real p_mu1 = sum(p_E_1);
    real norm = (p_mu0 + p_mu1);
    p_E_0 *= 1/norm;
    p_E_1 *= 1/norm;
    real Info = I(p_E_0, p_E_1, 500);
    real L = n * Info * log(2);
    
}