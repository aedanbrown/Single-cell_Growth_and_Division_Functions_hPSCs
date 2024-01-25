clear

%% Constants/Inputs
nb_i = "IMR90";
nb_iso_fName = "IMR90CNTRLNewborn.csv";
div_i = "IMR90";
div_iso_fName = "IMR90CNTRLDividing.csv";
mu = 0.0333;

%Using Rice Rule to determine histogram bin number
rice = @(x) ceil(2.*x.^(1/3));

n_x = 5000; %Number of points to use in x-axis


%% Load Data
%Load data and remove anything less than or equal to zeros

%Newborn
fullData = readtable(nb_i + "OCT4Newborn.csv");
nb = fullData.OCT4_R_PE_H;
nb(nb <= 0) = [];

%Total - newborn
fullData = readtable(nb_i + "TotalOCT4Newborn.csv");
total_nb = fullData.OCT4_R_PE_H;
total_nb(total_nb <= 0) = [];

%Newborn isotype
fullData = readtable(nb_iso_fName);
nb_iso = fullData.OCT4_R_PE_H;
nb_iso(nb_iso <= 0) = [];

%Dividing
fullData = readtable(div_i + "OCT4Dividing.csv");
div = fullData.OCT4_H;
div(div <= 0) = [];

%Total - dividing
fullData = readtable(div_i + "TotalOCT4Dividing.csv");
total_div = fullData.OCT4_H;
total_div(total_div <= 0) = [];

%Dividing isotype
fullData = readtable(div_iso_fName);
div_iso = fullData.OCT4_H;
div_iso(div_iso <= 0) = [];

%% Plot isotypes

%Calculate cutoff
cutoff = calcCutoff(nb_iso,div_iso);

%Plot scaled isotypes
x_iso = linspace(min([min(nb_iso) min(div_iso)]), ...
                 max([max(nb_iso) max(div_iso)]), ...
                 max([rice(length(nb_iso)) rice(length(div_iso))]));
    %Need a new x-axis to account for the data shift

figure()
semilogx(x_iso(2:end),histcounts(nb_iso,x_iso(1,:))./(length(nb_iso).*(x_iso(1,2)-x_iso(1,1))),'LineWidth',2,'Color',[0.93 0.69 0.13]); hold on
semilogx(x_iso(2:end),histcounts(div_iso,x_iso(1,:))./(length(div_iso).*(x_iso(1,2)-x_iso(1,1))),'LineWidth',2,'Color',[0.00 0.45 0.74]);
xline(prctile(div_iso,99)) %Plot cutoffs
xline(prctile(nb_iso,99))
xlim([x_iso(2).*0.9 cutoff.*2])
legend(["Newborn Iso";"Dividing Iso"],'Location','Northwest')


%% Calculate bins and make x-axes for plotting

x_ranges = zeros(4,2);

x_ranges(1,1) = min([0,min(nb)]);
x_ranges(1,2) = max(nb);
x_ranges(2,1) = min([0,min(total_nb)]);
x_ranges(2,2) = max(total_nb);

nbins_nb = rice(length(nb));
nbins_tnb = rice(length(total_nb));

x_ranges(3,1) = min([0,min(div)]);
x_ranges(3,2) = max(div);
x_ranges(4,1) = min([0,min(total_div)]);
x_ranges(4,2) = max(total_div);

nbins_div = rice(length(div));
nbins_tdiv = rice(length(total_div));


x_bins{1} = linspace(x_ranges(1,1),x_ranges(1,2),nbins_nb);
x_bins{2} = linspace(x_ranges(2,1),x_ranges(2,2),nbins_tnb);
x_bins{3} = linspace(x_ranges(3,1),x_ranges(3,2),nbins_div);
x_bins{4} = linspace(x_ranges(4,1),x_ranges(4,2),nbins_tdiv);


%% Plot original data

%Plot original data
figure()
subplot(1,4,1)
semilogx(x_bins{1}(2:end),histcounts(nb,x_bins{1})./(length(nb).*(x_bins{1}(2)-x_bins{1}(1))) ...
    ,'LineWidth',2,'Color',[0.93 0.69 0.13]); hold on
xline(cutoff)
subplot(1,4,2)
semilogx(x_bins{2}(2:end),histcounts(total_nb,x_bins{2})./(length(total_nb).*(x_bins{2}(2)-x_bins{2}(1))) ...
    ,'LineWidth',2,'Color',[0.85 0.33 0.10]); hold on
xline(cutoff)

subplot(1,4,3)
semilogx(x_bins{3}(2:end),histcounts(div,x_bins{3})./(length(div).*(x_bins{3}(2)-x_bins{3}(1))) ...
    ,'LineWidth',2,'Color',[0.00 0.45 0.74]); hold on   
xline(cutoff)
subplot(1,4,4)
semilogx(x_bins{4}(2:end),histcounts(total_div,x_bins{4})./(length(total_div).*(x_bins{4}(2)-x_bins{4}(1))) ...
    ,'LineWidth',2,'Color',[0.85 0.33 0.10]); hold on 
xline(cutoff)

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



%% Plot all estimated pdfs and estimated cdfs

figure()
subplot(2,1,1)
xline(cutoff); hold on;
semilogx(x_c+cutoff,nb_pdf,'Color',[0.93,0.69,0.13],'LineWidth',2)
semilogx(x_c+cutoff,div_pdf,'Color',[0.00 0.45 0.74],'LineWidth',2)
semilogx(x_c+cutoff,total_nb_pdf,'Color',[0.85 0.33 0.10],'LineWidth',2)
semilogx(x_c+cutoff,total_div_pdf,'Color',[0.85 0.33 0.10],'LineWidth',2)
set(gca,'XScale','log')
subplot(2,1,2)
xline(cutoff); hold on;
semilogx(x_c+cutoff,nb_cdf,'Color',[0.93,0.69,0.13],'LineWidth',2)
semilogx(x_c+cutoff,div_cdf,'Color',[0.00 0.45 0.74],'LineWidth',2)
semilogx(x_c+cutoff,total_nb_cdf,'Color',[0.85 0.33 0.10],'LineWidth',2)
semilogx(x_c+cutoff,total_div_cdf,'Color',[0.85 0.33 0.10],'LineWidth',2)
set(gca,'XScale','log')

all_total = vertcat(total_nb,total_div);


%% Compare estimated PDFs to the underlying data
figure()
subplot(3,1,1)
xline(cutoff); hold on
semilogx(x_bins{1}(2:end),histcounts(nb,x_bins{1})./(length(nb).*(x_bins{1}(2)-x_bins{1}(1))),'Color',[0.93 0.69 0.13]);
semilogx(x_c+cutoff,nb_pdf,'Color',[0.93,0.69,0.13],'LineWidth',2);
set(gca,'xscale','log')

subplot(3,1,2)
xline(cutoff); hold on
semilogx(x_bins{1}(2:end),histcounts(div,x_bins{1})./(length(div).*(x_bins{1}(2)-x_bins{1}(1))),'Color',[0.00 0.45 0.74]);
semilogx(x_c+cutoff,div_pdf,'Color',[0.00 0.45 0.74],'Linewidth',2);
set(gca,'xscale','log')


subplot(3,1,3)
xline(cutoff); hold on
semilogx(x_bins{1}(2:end),histcounts(total_nb,x_bins{1})./(length(total_nb).*(x_bins{1}(2)-x_bins{1}(1))),'Color',[0.85 0.33 0.10]);
semilogx(x_bins{1}(2:end),histcounts(total_div,x_bins{1})./(length(total_div).*(x_bins{1}(2)-x_bins{1}(1))),'Color',[0.85 0.33 0.10]);
semilogx(x_c+cutoff,mean([total_nb_pdf;total_div_pdf],1),'Color',[0.85 0.33 0.10],'Linewidth',2);
set(gca,'xscale','log')


%% Plot averages pdfs and cdfs
figure()
plot(x_c+cutoff,nb_pdf,'Color',[0.93,0.69,0.13],'LineWidth',2); hold on
plot(x_c+cutoff,div_pdf,'Color',[0.00 0.45 0.74],'LineWidth',2)
plot(x_c+cutoff,mean([total_nb_pdf;total_div_pdf],1),'Color',[0.85 0.33 0.10],'LineWidth',2)
legend(["Newborn";"Dividing";"Total"])
xline(cutoff,"HandleVisibility","off")
xlim([prctile(all_total,5) prctile(all_total,95)])

figure()
plot(x_c+cutoff,nb_cdf,'Color',[0.93,0.69,0.13],'LineWidth',2); hold on
plot(x_c+cutoff,div_cdf,'Color',[0.00 0.45 0.74],'LineWidth',2)
plot(x_c+cutoff,mean([total_nb_cdf;total_div_cdf],1),'Color',[0.85 0.33 0.10],'LineWidth',2)
legend(["Newborn";"Dividing";"Total"])
xline(cutoff,"HandleVisibility","off")
xlim([prctile(all_total,5) prctile(all_total,95)])

%% Calculate and plot PSFs
%All total
Gamma_est = mu.*div_pdf./mean([total_nb_pdf;total_div_pdf],1);
R_est = mu.*(2.*nb_cdf - mean([total_nb_cdf;total_div_cdf],1) - div_cdf)./mean([total_nb_pdf;total_div_pdf],1);


%Plot PSFs
figure()
subplot(1,2,1)
plot(x_c+cutoff,Gamma_est,'LineWidth',2)
xline(cutoff)
xline(prctile(all_total,10))
xline(prctile(all_total,90))
xlim([prctile(all_total,1) prctile(all_total,99)])
subplot(1,2,2)
plot(x_c+cutoff,R_est,'LineWidth',2)
xline(cutoff)
xline(prctile(all_total,10))
xline(prctile(all_total,90))
xlim([prctile(all_total,1) prctile(all_total,99)])

%% Functions

function [cutoff] = calcCutoff(nb_iso,div_iso)

    p = 99;
    cutoff = max([prctile(nb_iso,p) ...
        prctile(div_iso,p)]);
    
end


