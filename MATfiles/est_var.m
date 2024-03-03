function S_var = est_var(reference, bandw, M, N)

% Use Lemma 1

sigma_2_sq = hyy(reference, bandw);
sigma_4_sq = hxxyy(reference, bandw);

B = 2:1:M;
S_var = ( (sigma_4_sq/N) +(N-1)*sigma_2_sq/N ) *2./B./(B-1);


end
