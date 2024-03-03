function S_var = est_var_online (reference, bandw, M, N)


sigma_2_sq = hyy(reference, bandw);
sigma_4_sq = hxxyy(reference, bandw);

temp = (sigma_4_sq/N) +(N-1)*sigma_2_sq/N;
S_var = temp * 2/M/(M-1);

end