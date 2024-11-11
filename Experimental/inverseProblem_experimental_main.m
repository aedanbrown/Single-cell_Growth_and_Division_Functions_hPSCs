clear

H9_p = [0.02645 1.4803e+03];
IMR90_p = [0.03815 1.4785e+03];

nb_color = [0.93 0.69 0.13];
total_color = [0.85 0.33 0.10];
div_color = [0.00 0.45 0.74];
nonlac_psf_color = [0 0.47 0];
lac_psf_color = [1 0.6 0];

[Gamma_est_H9, R_est_H9, Gamma_est_nT_H9, R_est_nT_H9, x_c_H9, cutoff_H9, all_total_H9, ...
    nb_pdf_H9, div_pdf_H9, total_nb_pdf_H9, total_div_pdf_H9, ...
    nb_H9, div_H9, total_nb_H9, total_div_H9] = inverseProblem_experimental_H9_func();

[Gamma_est_IMR90, R_est_IMR90, Gamma_est_nT_IMR90, R_est_nT_IMR90, x_c_IMR90, cutoff_IMR90, all_total_IMR90, ...
    nb_pdf_IMR90, div_pdf_IMR90, total_nb_pdf_IMR90, total_div_pdf_IMR90, ...
    nb_IMR90, div_IMR90, total_nb_IMR90, total_div_IMR90] = inverseProblem_experimental_IMR90_func();

[Gamma_est_H9_lactate, R_est_H9_lactate, Gamma_est_nT_H9_lactate, R_est_nT_H9_lactate, x_c_H9_lactate, cutoff_H9_lactate, all_total_H9_lactate, ...
    nb_pdf_H9_lactate, div_pdf_H9_lactate, total_nb_pdf_H9_lactate, total_div_pdf_H9_lactate, ...
    nb_H9_lactate, div_H9_lactate, total_nb_H9_lactate, total_div_H9_lactate] = inverseProblem_experimental_H9_lactate_func();

[Gamma_est_IMR90_lactate, R_est_IMR90_lactate, Gamma_est_nT_IMR90_lactate, R_est_nT_IMR90_lactate, x_c_IMR90_lactate, cutoff_IMR90_lactate, all_total_IMR90_lactate, ...
    nb_pdf_IMR90_lactate, div_pdf_IMR90_lactate, total_nb_pdf_IMR90_lactate, total_div_pdf_IMR90_lactate, ...
    nb_IMR90_lactate, div_IMR90_lactate, total_nb_IMR90_lactate, total_div_IMR90_lactate] = inverseProblem_experimental_IMR90_lactate_func();

%% Non-lactate estimated pdfs and histograms
nonlac_dists = figure(); 
nonlac_dists.Units = 'inches';
nonlac_dists.Position = [0,0,6,7];
[x,h_nb,h_div,h_total_nb,h_total_div] = histogram_plot(nb_H9,div_H9,total_nb_H9,total_div_H9);
nb_a = subplot(3,2,1,'Units','inches','FontSize',12);
semilogx(x,h_nb,'Color',nb_color); hold on
xline(cutoff_H9,'Alpha',1,'LineWidth',2)
semilogx(x_c_H9+cutoff_H9,nb_pdf_H9,'LineWidth',2,'Color',nb_color)
set(gca,'xscale','log')
title("Newborn - n_{nb}(x)","FontWeight","normal")
xlim([1e3 2e5])
ylim([0 2e-4])
xticks([1e3,1e4,1e5])
tA = text(gca,0,nb_a.Position(4)*1.27,"A - H9 Cells",'Units','inches','FontWeight','bold','FontSize',16);

subplot(3,2,3,'Units','inches','FontSize',12)
semilogx(x,h_div,'Color',div_color); hold on
xline(cutoff_H9,'Alpha',1,'LineWidth',2)
semilogx(x_c_H9+cutoff_H9,div_pdf_H9,'LineWidth',2,'Color',div_color)
set(gca,'xscale','log')
title("Dividing - n_{d}(x)","FontWeight","normal")
xlim([1e3 2e5])
xticks([1e3,1e4,1e5])
ylim([0 2e-4])

subplot(3,2,5,'Units','inches','FontSize',12)
semilogx(x,h_total_nb,'Color',total_color); hold on
semilogx(x,h_total_div,'Color',total_color);
xline(cutoff_H9,'Alpha',1,'LineWidth',2)
semilogx(x_c_H9+cutoff_H9,mean([total_nb_pdf_H9;total_div_pdf_H9],1),'LineWidth',2,'Color',total_color)
set(gca,'xscale','log')
title("Total - n_{T}(x)","FontWeight","normal")
xlim([1e3 2e5])
xticks([1e3,1e4,1e5])
ylim([0 2e-4])
xlabel("OCT4")


[x,h_nb,h_div,h_total_nb,h_total_div] = histogram_plot(nb_IMR90,div_IMR90,total_nb_IMR90,total_div_IMR90);
nb_b = subplot(3,2,2,'Units','inches','FontSize',12);
semilogx(x,h_nb,'Color',nb_color); hold on
xline(cutoff_IMR90,'Alpha',1,'LineWidth',2)
semilogx(x_c_IMR90+cutoff_IMR90,nb_pdf_IMR90,'LineWidth',2,'Color',nb_color)
set(gca,'xscale','log')
title("Newborn - n_{nb}(x)","FontWeight","normal")
xlim([1e3 2e5])
xticks([1e3,1e4,1e5])
ylim([0 2e-4])
tB = text(gca,0,nb_b.Position(4)*1.25,"B - IMR90 Cells",'Units','inches','FontWeight','bold','FontSize',16);

subplot(3,2,4,'Units','inches','FontSize',12)
semilogx(x,h_div,'Color',div_color); hold on
xline(cutoff_IMR90,'Alpha',1,'LineWidth',2)
semilogx(x_c_IMR90+cutoff_IMR90,div_pdf_IMR90,'LineWidth',2,'Color',div_color)
set(gca,'xscale','log')
title("Dividing - n_{d}(x)","FontWeight","normal")
xlim([1e3 2e5])
xticks([1e3,1e4,1e5])
ylim([0 2e-4])

subplot(3,2,6,'Units','inches','FontSize',12)
semilogx(x,h_total_nb,'Color',total_color); hold on
semilogx(x,h_total_div,'Color',total_color);
xline(cutoff_IMR90,'Alpha',1,'LineWidth',2)
semilogx(x_c_IMR90+cutoff_IMR90,mean([total_nb_pdf_IMR90;total_div_pdf_IMR90],1),'LineWidth',2,'Color',total_color)
set(gca,'xscale','log')
title("Total - n_{T}(x)","FontWeight","normal")
xlim([1e3 2e5])
xticks([1e3,1e4,1e5])
ylim([0 2e-4])
xlabel("OCT4")

%% lactate estimated pdfs and histograms
lac_dists = figure(); 
lac_dists.Units = 'inches';
lac_dists.Position = [0,0,6,7];
[x,h_nb,h_div,h_total_nb,h_total_div] = histogram_plot(nb_H9_lactate,div_H9_lactate,total_nb_H9_lactate,total_div_H9_lactate);
nb_a = subplot(3,2,1,'Units','inches','FontSize',12);
semilogx(x,h_nb,'Color',nb_color); hold on
xline(cutoff_H9_lactate,'Alpha',1,'LineWidth',2)
semilogx(x_c_H9_lactate+cutoff_H9_lactate,nb_pdf_H9_lactate,'LineWidth',2,'Color',nb_color)
set(gca,'xscale','log')
title("Newborn - n_{nb}(x)","FontWeight","normal")
xlim([1e3 2e5])
ylim([0 6e-4])
xticks([1e3,1e4,1e5])
tA = text(gca,0,nb_a.Position(4)*1.27,"A - H9 Cells, Lactate",'Units','inches','FontWeight','bold','FontSize',16);

subplot(3,2,3,'Units','inches','FontSize',12)
semilogx(x,h_div,'Color',div_color); hold on
xline(cutoff_H9_lactate,'Alpha',1,'LineWidth',2)
semilogx(x_c_H9_lactate+cutoff_H9_lactate,div_pdf_H9_lactate,'LineWidth',2,'Color',div_color)
set(gca,'xscale','log')
title("Dividing - n_{d}(x)","FontWeight","normal")
xlim([1e3 2e5])
xticks([1e3,1e4,1e5])
ylim([0 6e-4])

subplot(3,2,5,'Units','inches','FontSize',12)
semilogx(x,h_total_nb,'Color',total_color); hold on
semilogx(x,h_total_div,'Color',total_color);
xline(cutoff_H9_lactate,'Alpha',1,'LineWidth',2)
semilogx(x_c_H9_lactate+cutoff_H9_lactate,mean([total_nb_pdf_H9_lactate;total_div_pdf_H9_lactate],1),'LineWidth',2,'Color',total_color)
set(gca,'xscale','log')
title("Total - n_{T}(x)","FontWeight","normal")
xlim([1e3 2e5])
xticks([1e3,1e4,1e5])
ylim([0 6e-4])
xlabel("OCT4")


[x,h_nb,h_div,h_total_nb,h_total_div] = histogram_plot(nb_IMR90_lactate,div_IMR90_lactate,total_nb_IMR90_lactate,total_div_IMR90_lactate);
nb_b = subplot(3,2,2,'Units','inches','FontSize',12);
semilogx(x,h_nb,'Color',nb_color); hold on
xline(cutoff_IMR90_lactate,'Alpha',1,'LineWidth',2)
semilogx(x_c_IMR90_lactate+cutoff_IMR90_lactate,nb_pdf_IMR90_lactate,'LineWidth',2,'Color',nb_color)
set(gca,'xscale','log')
title("Newborn - n_{nb}(x)","FontWeight","normal")
xlim([1e3 2e5])
xticks([1e3,1e4,1e5])
ylim([0 6e-4])
tB = text(gca,0,nb_b.Position(4)*1.25,"B - IMR90 Cells, Lactate",'Units','inches','FontWeight','bold','FontSize',16);

subplot(3,2,4,'Units','inches','FontSize',12)
semilogx(x,h_div,'Color',div_color); hold on
xline(cutoff_IMR90_lactate,'Alpha',1,'LineWidth',2)
semilogx(x_c_IMR90_lactate+cutoff_IMR90_lactate,div_pdf_IMR90_lactate,'LineWidth',2,'Color',div_color)
set(gca,'xscale','log')
title("Dividing - n_{d}(x)","FontWeight","normal")
xlim([1e3 2e5])
xticks([1e3,1e4,1e5])
ylim([0 6e-4])

subplot(3,2,6,'Units','inches','FontSize',12)
semilogx(x,h_total_nb,'Color',total_color); hold on
semilogx(x,h_total_div,'Color',total_color);
xline(cutoff_IMR90_lactate,'Alpha',1,'LineWidth',2)
semilogx(x_c_IMR90_lactate+cutoff_IMR90_lactate,mean([total_nb_pdf_IMR90_lactate;total_div_pdf_IMR90_lactate],1),'LineWidth',2,'Color',total_color)
set(gca,'xscale','log')
title("Total - n_{T}(x)","FontWeight","normal")
xlim([1e3 2e5])
xticks([1e3,1e4,1e5])
ylim([0 6e-4])
xlabel("OCT4")

%% non-lactate PSFs

nonlac_psfs = figure();
nonlac_psfs.Units = 'inches';
nonlac_psfs.Position = [0,0,6,6];

psf_a = subplot(2,2,1,'Units','inches','FontSize',12);
plot(x_c_H9+cutoff_H9,Gamma_est_H9,'LineWidth',2,'Color',nonlac_psf_color)
xline(cutoff_H9,'Alpha',1,'LineWidth',1)
xline(prctile(all_total_H9,90),'Alpha',1,'LineWidth',1,'LineStyle','--')
ylabel("\Gamma(OCT4) 1/hr")
xlim([0 3.5e4])
ylim([0 0.2])
tA = text(psf_a,0,psf_a.Position(4)*1.1,"A - H9 Cells",'Units','inches','FontWeight','bold','FontSize',16);

subplot(2,2,2,'Units','inches','FontSize',12)
plot(x_c_H9+cutoff_H9,R_est_H9,'LineWidth',2,'Color',nonlac_psf_color)
xline(cutoff_H9,'Alpha',1,'LineWidth',1)
xline(prctile(all_total_H9,90),'Alpha',1,'LineWidth',1,'LineStyle','--')
ylabel("R(OCT4) OCT4/hr")
xlim([0 3.5e4])
ylim([0 1200])

psf_b = subplot(2,2,3,'Units','inches','FontSize',12);
plot(x_c_IMR90+cutoff_IMR90,Gamma_est_IMR90,'LineWidth',2,'Color',nonlac_psf_color)
xline(cutoff_IMR90,'Alpha',1,'LineWidth',1)
xline(prctile(all_total_IMR90,90),'Alpha',1,'LineWidth',1,'LineStyle','--')
xlabel("OCT4")
ylabel("\Gamma(OCT4) 1/hr")
xlim([0 3.5e4])
ylim([0 0.2])
tB = text(psf_b,0,psf_b.Position(4)*1.1,"B - IMR90 Cells",'Units','inches','FontWeight','bold','FontSize',16);

subplot(2,2,4,'Units','inches','FontSize',12)
plot(x_c_IMR90+cutoff_IMR90,R_est_IMR90,'LineWidth',2,'Color',nonlac_psf_color)
xline(cutoff_IMR90,'Alpha',1,'LineWidth',1)
xline(prctile(all_total_IMR90,90),'Alpha',1,'LineWidth',1,'LineStyle','--')
xlabel("OCT4")
ylabel("R(OCT4) OCT4/hr")
xlim([0 3.5e4])
ylim([0 1200])

%% lactate PSFs

lac_psfs = figure();
lac_psfs.Units = 'inches';
lac_psfs.Position = [0,0,6,6];

psf_a = subplot(2,2,1,'Units','inches','FontSize',12);
plot(x_c_H9_lactate+cutoff_H9_lactate,Gamma_est_H9_lactate,'LineWidth',2,'Color',lac_psf_color)
xline(cutoff_H9_lactate,'Alpha',1,'LineWidth',1)
xline(prctile(all_total_H9_lactate,90),'Alpha',1,'LineWidth',1,'LineStyle','--')
ylabel("\Gamma(OCT4) 1/hr")
xlim([0 2e4])
ylim([0 0.1])
tA = text(psf_a,0,psf_a.Position(4)*1.1,"A - H9 Cells, Lactate",'Units','inches','FontWeight','bold','FontSize',16);

subplot(2,2,2,'Units','inches','FontSize',12)
plot(x_c_H9_lactate+cutoff_H9_lactate,R_est_H9_lactate,'LineWidth',2,'Color',lac_psf_color)
xline(cutoff_H9_lactate,'Alpha',1,'LineWidth',1)
xline(prctile(all_total_H9_lactate,90),'Alpha',1,'LineWidth',1,'LineStyle','--')
ylabel("R(OCT4) OCT4/hr")
xlim([0 2e4])
ylim([-25 400])

psf_b = subplot(2,2,3,'Units','inches','FontSize',12);
plot(x_c_IMR90_lactate+cutoff_IMR90_lactate,Gamma_est_IMR90_lactate,'LineWidth',2,'Color',lac_psf_color)
xline(cutoff_IMR90_lactate,'Alpha',1,'LineWidth',1)
xline(prctile(all_total_IMR90_lactate,90),'Alpha',1,'LineWidth',1,'LineStyle','--')
xlabel("OCT4")
ylabel("\Gamma(OCT4) 1/hr")
xlim([0 2e4])
ylim([0 0.1])
tB = text(psf_b,0,psf_b.Position(4)*1.1,"B - IMR90 Cells, Lactate",'Units','inches','FontWeight','bold','FontSize',16);

subplot(2,2,4,'Units','inches','FontSize',12)
plot(x_c_IMR90_lactate+cutoff_IMR90_lactate,R_est_IMR90_lactate,'LineWidth',2,'Color',lac_psf_color)
xline(cutoff_IMR90_lactate,'Alpha',1,'LineWidth',1)
xline(prctile(all_total_IMR90_lactate,90),'Alpha',1,'LineWidth',1,'LineStyle','--')
xlabel("OCT4")
ylabel("R(OCT4) OCT4/hr")
xlim([0 2e4])
ylim([-25 400])


%% scaled PSFs

scal_psfs = figure();
scal_psfs.Units = 'inches';
scal_psfs.Position = [0,0,6,6];

psf_a = subplot(2,2,1,'Units','inches','FontSize',12);
plot(x_c_H9+cutoff_H9,Gamma_est_nT_H9,'LineWidth',2,'Color',nonlac_psf_color); hold on
plot(x_c_H9_lactate+cutoff_H9_lactate,Gamma_est_nT_H9_lactate,'LineWidth',2,'Color',lac_psf_color)
ylabel("\Gamma*n_T 1/hr*OCT4")
xlim([0 5e4])
ylim([0 4e-6])
tA = text(psf_a,0,psf_a.Position(4)*1.15,"A - H9 Cells",'Units','inches','FontWeight','bold','FontSize',16);


subplot(2,2,2,'Units','inches','FontSize',12)
plot(x_c_H9+cutoff_H9,R_est_nT_H9,'LineWidth',2,'Color',nonlac_psf_color); hold on
plot(x_c_H9_lactate+cutoff_H9_lactate,R_est_nT_H9_lactate,'LineWidth',2,'Color',lac_psf_color)
ylabel("R*n_T 1/hr")
xlim([0 5e4])
ylim([0 0.04])
legend(["-Lactate";"+Lactate"])

psf_b = subplot(2,2,3,'Units','inches','FontSize',12);
plot(x_c_IMR90+cutoff_IMR90,Gamma_est_nT_IMR90,'LineWidth',2,'Color',nonlac_psf_color); hold on
plot(x_c_IMR90_lactate+cutoff_IMR90_lactate,Gamma_est_nT_IMR90_lactate,'LineWidth',2,'Color',lac_psf_color)
xlabel("OCT4")
ylabel("\Gamma*n_T 1/hr*OCT4")
xlim([0 5e4])
ylim([0 4e-6])
tB = text(psf_b,0,psf_b.Position(4)*1.15,"B - IMR90 Cells",'Units','inches','FontWeight','bold','FontSize',16);

subplot(2,2,4,'Units','inches','FontSize',12)
plot(x_c_IMR90+cutoff_IMR90,R_est_nT_IMR90,'LineWidth',2,'Color',nonlac_psf_color); hold on
plot(x_c_IMR90_lactate+cutoff_IMR90_lactate,R_est_nT_IMR90_lactate,'LineWidth',2,'Color',lac_psf_color)
xlabel("OCT4")
ylabel("R*n_T 1/hr")
xlim([0 5e4])
ylim([0 0.04])



%% Averages calculations

fprintf("OCT4 H9\n")
average_f(x_c_H9+cutoff_H9,mean([total_nb_pdf_H9;total_div_pdf_H9],1).*(x_c_H9+cutoff_H9),nb_pdf_H9.*(x_c_H9+cutoff_H9),div_pdf_H9.*(x_c_H9+cutoff_H9))

fprintf("OCT4 IMR90\n")
average_f(x_c_IMR90+cutoff_IMR90,mean([total_nb_pdf_IMR90;total_div_pdf_IMR90],1).*(x_c_IMR90+cutoff_IMR90),nb_pdf_IMR90.*(x_c_IMR90+cutoff_IMR90),div_pdf_IMR90.*(x_c_IMR90+cutoff_IMR90))

fprintf("OCT4 H9 Lactate\n")
average_f(x_c_H9_lactate+cutoff_H9_lactate,mean([total_nb_pdf_H9_lactate;total_div_pdf_H9_lactate],1).*(x_c_H9_lactate+cutoff_H9_lactate),nb_pdf_H9_lactate.*(x_c_H9_lactate+cutoff_H9_lactate),div_pdf_H9_lactate.*(x_c_H9_lactate+cutoff_H9_lactate))

fprintf("OCT4 IMR90 Lactate\n")
average_f(x_c_IMR90_lactate+cutoff_IMR90_lactate,mean([total_nb_pdf_IMR90_lactate;total_div_pdf_IMR90_lactate],1).*(x_c_IMR90_lactate+cutoff_IMR90_lactate),nb_pdf_IMR90_lactate.*(x_c_IMR90_lactate+cutoff_IMR90_lactate),div_pdf_IMR90_lactate.*(x_c_IMR90_lactate+cutoff_IMR90_lactate))






fprintf("R H9\n")
H9_mask = x_c_H9 + cutoff_H9 < prctile(all_total_H9,90) & isfinite(R_est_H9);
average_f(x_c_H9(H9_mask)+cutoff_H9,mean([total_nb_pdf_H9(H9_mask);total_div_pdf_H9(H9_mask)],1).*R_est_H9(H9_mask),nb_pdf_H9(H9_mask).*R_est_H9(H9_mask),div_pdf_H9(H9_mask).*R_est_H9(H9_mask))

fprintf("R IMR90\n")
IMR90_mask = x_c_IMR90 + cutoff_IMR90 < prctile(all_total_IMR90,90) & isfinite(R_est_IMR90);
average_f(x_c_IMR90(IMR90_mask)+cutoff_IMR90,mean([total_nb_pdf_IMR90(IMR90_mask);total_div_pdf_IMR90(IMR90_mask)],1).*R_est_IMR90(IMR90_mask),nb_pdf_IMR90(IMR90_mask).*R_est_IMR90(IMR90_mask),div_pdf_IMR90(IMR90_mask).*R_est_IMR90(IMR90_mask))

fprintf("R H9 Lactate\n")
H9_mask = x_c_H9_lactate + cutoff_H9_lactate < prctile(all_total_H9_lactate,90) & isfinite(R_est_H9_lactate);
average_f(x_c_H9_lactate(H9_mask)+cutoff_H9_lactate,mean([total_nb_pdf_H9_lactate(H9_mask);total_div_pdf_H9_lactate(H9_mask)],1).*R_est_H9_lactate(H9_mask),nb_pdf_H9_lactate(H9_mask).*R_est_H9_lactate(H9_mask),div_pdf_H9_lactate(H9_mask).*R_est_H9_lactate(H9_mask))

fprintf("R IMR90 Lactate\n")
IMR90_mask = x_c_IMR90_lactate + cutoff_IMR90_lactate < prctile(all_total_IMR90_lactate,90) & isfinite(R_est_IMR90_lactate);
average_f(x_c_IMR90_lactate(IMR90_mask)+cutoff_IMR90_lactate,mean([total_nb_pdf_IMR90_lactate(IMR90_mask);total_div_pdf_IMR90_lactate(IMR90_mask)],1).*R_est_IMR90_lactate(IMR90_mask),nb_pdf_IMR90_lactate(IMR90_mask).*R_est_IMR90_lactate(IMR90_mask),div_pdf_IMR90_lactate(IMR90_mask).*R_est_IMR90_lactate(IMR90_mask))







fprintf("Gamma H9\n")
H9_mask = x_c_H9 + cutoff_H9 < prctile(all_total_H9,90) & isfinite(Gamma_est_H9);
average_f(x_c_H9(H9_mask)+cutoff_H9,mean([total_nb_pdf_H9(H9_mask);total_div_pdf_H9(H9_mask)],1).*Gamma_est_H9(H9_mask),nb_pdf_H9(H9_mask).*Gamma_est_H9(H9_mask),div_pdf_H9(H9_mask).*Gamma_est_H9(H9_mask))

fprintf("Gamma IMR90\n")
IMR90_mask = x_c_IMR90 + cutoff_IMR90 < prctile(all_total_IMR90,90) & isfinite(Gamma_est_IMR90);
average_f(x_c_IMR90(IMR90_mask)+cutoff_IMR90,mean([total_nb_pdf_IMR90(IMR90_mask);total_div_pdf_IMR90(IMR90_mask)],1).*Gamma_est_IMR90(IMR90_mask),nb_pdf_IMR90(IMR90_mask).*Gamma_est_IMR90(IMR90_mask),div_pdf_IMR90(IMR90_mask).*Gamma_est_IMR90(IMR90_mask))

fprintf("Gamma H9 Lactate\n")
H9_mask = x_c_H9_lactate + cutoff_H9_lactate < prctile(all_total_H9_lactate,90) & isfinite(Gamma_est_H9_lactate);
average_f(x_c_H9_lactate(H9_mask)+cutoff_H9_lactate,mean([total_nb_pdf_H9_lactate(H9_mask);total_div_pdf_H9_lactate(H9_mask)],1).*Gamma_est_H9_lactate(H9_mask),nb_pdf_H9_lactate(H9_mask).*Gamma_est_H9_lactate(H9_mask),div_pdf_H9_lactate(H9_mask).*Gamma_est_H9_lactate(H9_mask))

fprintf("Gamma IMR90 Lactate\n")
IMR90_mask = x_c_IMR90_lactate + cutoff_IMR90_lactate < prctile(all_total_IMR90_lactate,90) & isfinite(Gamma_est_IMR90_lactate);
average_f(x_c_IMR90_lactate(IMR90_mask)+cutoff_IMR90_lactate,mean([total_nb_pdf_IMR90_lactate(IMR90_mask);total_div_pdf_IMR90_lactate(IMR90_mask)],1).*Gamma_est_IMR90_lactate(IMR90_mask),nb_pdf_IMR90_lactate(IMR90_mask).*Gamma_est_IMR90_lactate(IMR90_mask),div_pdf_IMR90_lactate(IMR90_mask).*Gamma_est_IMR90_lactate(IMR90_mask))




%% NANOG version - non-lactate PSFs

nonlac_psfs = figure();
nonlac_psfs.Units = 'inches';
nonlac_psfs.Position = [0,0,6,6];

psf_a = subplot(2,2,1,'Units','inches','FontSize',12);
plot((x_c_H9+cutoff_H9).*H9_p(1) + H9_p(2),Gamma_est_H9,'LineWidth',2,'Color',nonlac_psf_color)
xline(cutoff_H9*H9_p(1) + H9_p(2),'Alpha',1,'LineWidth',1)
xline(prctile(all_total_H9,90)*H9_p(1) + H9_p(2),'Alpha',1,'LineWidth',1,'LineStyle','--')
ylabel("\Gamma(NANOG) 1/hr")
xlim([1450 2600])
ylim([0 0.2])
tA = text(psf_a,0,psf_a.Position(4)*1.1,"A - H9 Cells",'Units','inches','FontWeight','bold','FontSize',16);

subplot(2,2,2,'Units','inches','FontSize',12)
plot((x_c_H9+cutoff_H9).*H9_p(1) + H9_p(2),R_est_H9,'LineWidth',2,'Color',nonlac_psf_color)
xline(cutoff_H9*H9_p(1) + H9_p(2),'Alpha',1,'LineWidth',1)
xline(prctile(all_total_H9,90)*H9_p(1) + H9_p(2),'Alpha',1,'LineWidth',1,'LineStyle','--')
ylabel("R(NANOG) NANOG/hr")
xlim([1450 2600])
ylim([0 1200])

psf_b = subplot(2,2,3,'Units','inches','FontSize',12);
plot((x_c_IMR90+cutoff_IMR90).*IMR90_p(1) + IMR90_p(2),Gamma_est_IMR90,'LineWidth',2,'Color',nonlac_psf_color)
xline(cutoff_IMR90.*IMR90_p(1) + IMR90_p(2),'Alpha',1,'LineWidth',1)
xline(prctile(all_total_IMR90,90).*IMR90_p(1) + IMR90_p(2),'Alpha',1,'LineWidth',1,'LineStyle','--')
xlabel("NANOG")
ylabel("\Gamma(NANOG) 1/hr")
xlim([1450 2600])
ylim([0 0.2])
tB = text(psf_b,0,psf_b.Position(4)*1.1,"B - IMR90 Cells",'Units','inches','FontWeight','bold','FontSize',16);

subplot(2,2,4,'Units','inches','FontSize',12)
plot((x_c_IMR90+cutoff_IMR90).*IMR90_p(1) + IMR90_p(2),R_est_IMR90,'LineWidth',2,'Color',nonlac_psf_color)
xline(cutoff_IMR90.*IMR90_p(1) + IMR90_p(2),'Alpha',1,'LineWidth',1)
xline(prctile(all_total_IMR90,90).*IMR90_p(1) + IMR90_p(2),'Alpha',1,'LineWidth',1,'LineStyle','--')
xlabel("NANOG")
ylabel("R(NANOG) NANOG/hr")
xlim([1450 2600])
ylim([0 1200])



%% Non-lactate estimated pdfs and histograms - NANOG
nonlac_dists = figure(); 
nonlac_dists.Units = 'inches';
nonlac_dists.Position = [0,0,6.8,7];
nb_a = subplot(3,2,1,'Units','inches','FontSize',12);
xline(cutoff_H9.*H9_p(1) + H9_p(2),'Alpha',1,'LineWidth',2); hold on
semilogx((x_c_H9+cutoff_H9).*H9_p(1) + H9_p(2),nb_pdf_H9./H9_p(1),'LineWidth',2,'Color',nb_color)
set(gca,'xscale','log')
title("Newborn - n_{nb}(x)","FontWeight","normal")
xlim([1.4e3 4e3])
xticks([1500,2000,3000,4000])
ylim([0 7e-3])
tA = text(gca,0,nb_a.Position(4)*1.27,"A - H9 Cells",'Units','inches','FontWeight','bold','FontSize',16);

subplot(3,2,3,'Units','inches','FontSize',12)
xline(cutoff_H9.*H9_p(1) + H9_p(2),'Alpha',1,'LineWidth',2); hold on
semilogx((x_c_H9+cutoff_H9).*H9_p(1) + H9_p(2),div_pdf_H9./H9_p(1),'LineWidth',2,'Color',div_color)
set(gca,'xscale','log')
title("Dividing - n_{d}(x)","FontWeight","normal")
xlim([1.4e3 4e3])
xticks([1500,2000,3000,4000])
ylim([0 7e-3])

subplot(3,2,5,'Units','inches','FontSize',12)
xline(cutoff_H9.*H9_p(1) + H9_p(2),'Alpha',1,'LineWidth',2); hold on
semilogx((x_c_H9+cutoff_H9).*H9_p(1) + H9_p(2),mean([total_nb_pdf_H9;total_div_pdf_H9],1)./H9_p(1),'LineWidth',2,'Color',total_color)
set(gca,'xscale','log')
title("Total - n_{T}(x)","FontWeight","normal")
xlim([1.4e3 4e3])
xticks([1500,2000,3000,4000])
ylim([0 7e-3])
xlabel("NANOG")


[x,h_nb,h_div,h_total_nb,h_total_div] = histogram_plot(nb_IMR90,div_IMR90,total_nb_IMR90,total_div_IMR90);
nb_b = subplot(3,2,2,'Units','inches','FontSize',12);
xline(cutoff_IMR90.*IMR90_p(1) + IMR90_p(2),'Alpha',1,'LineWidth',2); hold on
semilogx((x_c_IMR90+cutoff_IMR90).*IMR90_p(1) + IMR90_p(2),nb_pdf_IMR90./IMR90_p(1),'LineWidth',2,'Color',nb_color)
set(gca,'xscale','log')
title("Newborn - n_{nb}(x)","FontWeight","normal")
xlim([1.4e3 4e3])
xticks([1500,2000,3000,4000])
ylim([0 7e-3])
tB = text(gca,0,nb_b.Position(4)*1.25,"B - IMR90 Cells",'Units','inches','FontWeight','bold','FontSize',16);

subplot(3,2,4,'Units','inches','FontSize',12)
xline(cutoff_IMR90.*IMR90_p(1) + IMR90_p(2),'Alpha',1,'LineWidth',2); hold on
semilogx((x_c_IMR90+cutoff_IMR90).*IMR90_p(1) + IMR90_p(2),div_pdf_IMR90./IMR90_p(1),'LineWidth',2,'Color',div_color)
set(gca,'xscale','log')
title("Dividing - n_{d}(x)","FontWeight","normal")
xlim([1.4e3 4e3])
xticks([1500,2000,3000,4000])
ylim([0 7e-3])

subplot(3,2,6,'Units','inches','FontSize',12)
xline(cutoff_IMR90.*IMR90_p(1) + IMR90_p(2),'Alpha',1,'LineWidth',2); hold on
semilogx((x_c_IMR90+cutoff_IMR90).*IMR90_p(1) + IMR90_p(2),mean([total_nb_pdf_IMR90;total_div_pdf_IMR90],1)./IMR90_p(1),'LineWidth',2,'Color',total_color)
set(gca,'xscale','log')
title("Total - n_{T}(x)","FontWeight","normal")
xlim([1.4e3 4e3])
xticks([1500,2000,3000,4000])
ylim([0 7e-3])
xlabel("NANOG")

%%
function [bins,h_nb,h_div,h_total_nb,h_total_div] = histogram_plot(nb,div,total_nb,total_div)

    rice = @(c) ceil(2.*c.^(1/3));

    n = rice(length(nb));

    x_lower = min([0,min(nb)]);
    x_upper = max(nb);
    
    x_bins = linspace(x_lower,x_upper,n);

    h_nb = histcounts(nb,x_bins)./(length(nb).*(x_bins(2)-x_bins(1)));
    h_div = histcounts(div,x_bins)./(length(div).*(x_bins(2)-x_bins(1)));
    h_total_nb = histcounts(total_nb,x_bins)./(length(total_nb).*(x_bins(2)-x_bins(1)));
    h_total_div = histcounts(total_div,x_bins)./(length(total_div).*(x_bins(2)-x_bins(1)));
    bins = x_bins(2:end);


end


function average_f(x,f_total,f_nb,f_div)

fprintf("Average Total: %f\n",trapz(x,f_total))
fprintf("Average Newborn: %f\n",trapz(x,f_nb))
fprintf("Average Dividing: %f\n", trapz(x,f_div))

end
