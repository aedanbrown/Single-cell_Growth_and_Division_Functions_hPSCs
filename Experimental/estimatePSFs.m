function [Gamma_est,R_est,Gamma_est_nT,R_est_nT] = estimatePSFs(mu,nb_pdf,total_nb_pdf,div_pdf,total_div_pdf,nb_cdf,total_nb_cdf,div_cdf,total_div_cdf)

Gamma_est = mu.*div_pdf./mean([total_nb_pdf;total_div_pdf],1);
R_est = mu.*(2.*nb_cdf - mean([total_nb_cdf;total_div_cdf],1) - div_cdf)./mean([total_nb_pdf;total_div_pdf],1);

Gamma_est_nT = mu.*div_pdf;
R_est_nT = mu.*(2.*nb_cdf - mean([total_nb_cdf;total_div_cdf],1) - div_cdf);


end