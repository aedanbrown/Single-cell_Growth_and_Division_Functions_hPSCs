function [Gamma_est,R_est,Gamma_est_nT,R_est_nT,x_c,cutoff,all_total,nb_pdf,div_pdf,total_nb_pdf,total_div_pdf,nb,div,total_nb,total_div] = inverseProblem_experimental_IMR90_lactate_func()

clear

%% Constants/Inputs
nb_i = "IMR90";
div_i = "IMR90";
iso_fName = "IMR90CNTRL_Lactate.csv";
mu = 0.0190;

%Using Rice Rule to determine histogram bin number
rice = @(x) ceil(2.*x.^(1/3));

n_x = 5000; %Number of points to use in x-axis


%% Load Data
%Load data and remove anything less than or equal to zeros

%Newborn
fullData = readtable(nb_i + "OCT4Newborn_Lactate.csv");
nb = fullData.OCT4_H;
nb(nb <= 0) = [];

%Total - newborn
fullData = readtable(nb_i + "TotalOCT4Newborn_Lactate.csv");
total_nb = fullData.OCT4_H;
total_nb(total_nb <= 0) = [];

%Isotype
fullData = readtable(iso_fName);
iso = fullData.OCT4_H;
iso(iso <= 0) = [];

%Dividing
fullData = readtable(div_i + "OCT4Dividing_Lactate.csv");
div = fullData.OCT4_H;
div(div <= 0) = [];

%Total - dividing
fullData = readtable(div_i + "TotalOCT4Dividing_Lactate.csv");
total_div = fullData.OCT4_H;
total_div(total_div <= 0) = [];


%% Plot isotypes

%Calculate cutoff
cutoff = prctile(iso,99);





%% Remove data above the cutoff and calculate the truncated KDE

nb_mask = nb > cutoff;
div_mask = div > cutoff;
total_nb_mask = total_nb > cutoff;
total_div_mask = total_div > cutoff;

%Truncated kernel
nb_est = fitdist(nb(nb_mask)-cutoff,'Kernel','Kernel','epanechnikov','support','positive');
div_est = fitdist(div(div_mask)-cutoff,'Kernel','Kernel','epanechnikov','support','positive'); 
total_nb_est = fitdist(total_nb(total_nb_mask)-cutoff,'Kernel','Kernel','epanechnikov','support','positive');
total_div_est = fitdist(total_div(total_div_mask)-cutoff,'Kernel','Kernel','epanechnikov','support','positive');


x_min = min([min(total_nb) min(total_div)]);
x_max = max([max(total_nb) max(total_div)]);


%% Estimate PDFs and CDFs

x_c = linspace(0,1.1*x_max-0.1.*x_min,n_x); %x-axis for the pdfs and cdfs

nb_pdf = pdf(nb_est,x_c).*(1 - length(nb(~nb_mask))./length(nb));
div_pdf = pdf(div_est,x_c).*(1 - length(div(~div_mask))./length(div));   
total_nb_pdf = pdf(total_nb_est,x_c).*(1 - length(total_nb(~total_nb_mask))./length(total_nb));
total_div_pdf = pdf(total_div_est,x_c).*(1 - length(total_div(~total_div_mask))./length(total_div));

nb_cdf = cdf(nb_est,x_c).*(1 - length(nb(~nb_mask))./length(nb)) + length(nb(~nb_mask))./length(nb);
div_cdf = cdf(div_est,x_c).*(1 - length(div(~div_mask))./length(div)) + length(div(~div_mask))./length(div);
total_nb_cdf = cdf(total_nb_est,x_c).*(1 - length(total_nb(~total_nb_mask))./length(total_nb)) + length(total_nb(~total_nb_mask))./length(total_nb);
total_div_cdf = cdf(total_div_est,x_c).*(1 - length(total_div(~total_div_mask))./length(total_div)) + length(total_div(~total_div_mask))./length(total_div);



all_total = vertcat(total_nb,total_div);



%% Calculate and plot PSFs
%All total
% Gamma_est = mu.*div_pdf; %./mean([total_nb_pdf;total_div_pdf],1);
% R_est = mu.*(2.*nb_cdf - mean([total_nb_cdf;total_div_cdf],1) - div_cdf); %./mean([total_nb_pdf;total_div_pdf],1);

[Gamma_est,R_est,Gamma_est_nT,R_est_nT] = estimatePSFs(mu,nb_pdf,total_nb_pdf,div_pdf,total_div_pdf,nb_cdf,total_nb_cdf,div_cdf,total_div_cdf);

end



