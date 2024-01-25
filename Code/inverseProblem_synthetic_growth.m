function [] = inverseProblem_synthetic_growth(fName,Gamma,R)

    %fName is the root file name for the files that contain the various
    %distributions, cell concentrations, and time values
    %Gamma is a function handle for the Gamma PSF used
    %R is a function handle for the R PSF used

    Xi = readmatrix(sprintf('Xi_%s.csv',fName));
    Cpt = readmatrix(sprintf('Cpt_%s.csv',fName));
    nd = readmatrix(sprintf('nd_%s.csv',fName));
    nnb = readmatrix(sprintf('nnb_%s.csv',fName));
    tauV = readmatrix(sprintf('tauV_%s.csv',fName));

    ndist = 1000;

    %Calculate rate constants
    figure()
    plot(tauV,log(Cpt),'LineWidth',2)
    legend("ln(Cpt)")
    drawnow

    tMin = input("When does linear region begin? ");

    xline(tMin,'LineWidth',2,"HandleVisibility","off");
    exportgraphics(gcf,sprintf("g_%s_%i_concs.png",fName,tMin))
    savefig(sprintf("f_%s_%i_concs.fig",fName,tMin))

    mask = tauV > tMin;

    pf = polyfit(tauV(mask),log(Cpt(mask)),1);

    figure()
    plot(tauV(mask),log(Cpt(mask)),'LineWidth',2)
    hold on
    plot(tauV(mask),pf(1).*tauV(mask)+pf(2),'LineWidth',2)
    legend(["Data";"Linear Fit"])
    exportgraphics(gcf,sprintf("g_%s_%i_muFit.png",fName,tMin))
    savefig(sprintf("f_%s_%i_muFit.fig",fName,tMin))

    covMat = corrcoef(log(Cpt(mask)),pf(1).*tauV(mask)+pf(2));
    fprintf("mu R^2: %.5f\n",covMat(1,2).^2)

    mu = pf(1);

    %Calculate distributions

    nT = reshape(Xi(end,:),[],1);
    fprintf("Number of total cells: %.1f\n",length(nT));

    nT_pd = fitdist(nT,'Kernel','Kernel','epanechnikov','Support','positive');

    nd = nd(nd > 0);
    fprintf("Number of dividing cells: %.1f\n",length(nd));
    if length(nd) > length(nT)
        nd = nd(end-length(nT)+1:end);
        fprintf("Number of dividing cells used: %.1f\n",length(nd));
    end

    nnb = nnb(nnb > 0);
    fprintf("Number of newborn cells: %.1f\n",length(nnb));
    if length(nnb) > 2*length(nT)
        nnb = nnb(end-2*length(nT)+1:end);
        fprintf("Number of newborn cells used: %.1f\n",length(nnb));
    end

    nd_pd = fitdist(nd,'Kernel','Kernel','epanechnikov','support','positive');
    nnb_pd = fitdist(nnb,'Kernel','Kernel','epanechnikov','support','positive');

    %Calculate pdfs

    x = linspace(0.1*min(Xi,[],'all'),1.1.*max(Xi,[],'all'),ndist);

    nd_pdf = pdf(nd_pd,x);
    nT_pdf = pdf(nT_pd,x);
    nnb_pdf = pdf(nnb_pd,x);

    figure()
    plot(x,nd_pdf,'LineWidth',2); hold on
    plot(x,nT_pdf,'LineWidth',2)
    plot(x,nnb_pdf,'LineWidth',2)
    legend(["Dividing";"Total";"Newborn"])
    exportgraphics(gcf,sprintf("g_%s_%i_dists.png",fName,tMin))
    savefig(sprintf("f_%s_%i_dists.fig",fName,tMin)) 

    %Calculate cdfs
    nd_cdf = cdf(nd_pd,x);
    nT_cdf = cdf(nT_pd,x);
    nnb_cdf = cdf(nnb_pd,x);


    %Calculate PSFs

    Gamma_est = mu.*nd_pdf./nT_pdf;
    R_est = mu.*(2.*nnb_cdf - nT_cdf - nd_cdf)./nT_pdf;

    prctile10 = prctile(nT,10);
    prctile90 = prctile(nT,90);

    x_mask = all([x > prctile10; x < prctile90],1);

    figure()
    subplot(1,2,1)
    G_x = Gamma(x);
    plot(x,Gamma_est,'LineWidth',2); hold on; plot(x,G_x,'LineWidth',2);
    ylim([min(G_x).*0.9 max(G_x).*1.1])
    xline(prctile10)
    xline(prctile90)
    subplot(1,2,2)
    R_x = R(x);
    plot(x,R_est,'LineWidth',2); hold on; plot(x,R_x,'LineWidth',2);
    ylim([min(R_x).*0.9 max(R_x).*1.1])
    legend(["Estimated";"Original"])
    xline(prctile10,'HandleVisibility','off')
    xline(prctile90,'HandleVisibility','off')
    exportgraphics(gcf,sprintf("g_%s_%i_PSFs.png",fName,tMin))
    savefig(sprintf("f_%s_%i_PSFs.fig",fName,tMin))

    fprintf("mu: %.5f\n",mu)

    MSE_Gamma = mean((G_x(x_mask) - Gamma_est(x_mask)).^2);
    MSE_R = mean((R_x(x_mask) - R_est(x_mask)).^2);

    fprintf("# of samples in 10-90: %i\n",length(x(x_mask)))
    fprintf("Gamma MSE: %.5e\n",MSE_Gamma)
    fprintf("R MSE: %.5e\n",MSE_R)

    fprintf("10th percentile: %.3f 90th percentile: %.3f\n",prctile10,prctile90)
    fprintf("Gamma_10: %.5f Gamma_90: %.5f\n",Gamma(prctile10),Gamma(prctile90))      
    fprintf("R_10: %.5f R_90: %.5f\n",R(prctile10),R(prctile90))




end