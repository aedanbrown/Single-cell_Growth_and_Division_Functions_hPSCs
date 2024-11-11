function [] = inverseProblem_diffSim(fName,Gamma,R,D)

    %fName is the root file name for the files that contain the various
    %distributions, cell concentrations, and time values
    %Gamma is a function handle for the Gamma PSF used
    %R is a function handle for the R PSF used
    %D is a function handle for the D PSF used

    Xi = readmatrix(sprintf('Xi_%s.csv',fName));
    Cpt = readmatrix(sprintf('Cpt_%s.csv',fName));
    Cdt = readmatrix(sprintf('Cdt_%s.csv',fName));
    nd = readmatrix(sprintf('nd_%s.csv',fName));
    nnb = readmatrix(sprintf('nnb_%s.csv',fName));
    ndiff = readmatrix(sprintf('ndiff_%s.csv',fName));
    tauV = readmatrix(sprintf('tauV_%s.csv',fName));

    ndist = 1000;

    %Calculate rate constants
    figure()
    plot(tauV,log(Cpt),'LineWidth',2)
    hold on
    plot(tauV,log(Cdt),'LineWidth',2)
    legend(["ln(Cpt)";"ln(Cdt)"])

    drawnow

    tMin = input("When does linear region begin? ");

    xline(tMin,'LineWidth',2,'HandleVisibility','off');
    exportgraphics(gcf,sprintf("g_%s_%i_concs.png",fName,tMin))
    savefig(sprintf("f_%s_%i_concs.fig",fName,tMin))

    mask = tauV > tMin;

    pf = polyfit(tauV(mask),log(Cpt(mask)),1);

    figure()
    plot(tauV(mask),log(Cpt(mask)),'LineWidth',2)
    hold on
    plot(tauV(mask),pf(1).*tauV(mask)+pf(2),'LineWidth',2)
    legend(["Data";"Linear Fit"])
    exportgraphics(gcf,sprintf("g_%s_%i_kFit.png",fName,tMin))
    savefig(sprintf("f_%s_%i_kFit.fig",fName,tMin))

    covMat = corrcoef(log(Cpt(mask)),pf(1).*tauV(mask)+pf(2));
    fprintf("k R^2: %.5f\n",covMat(1,2).^2)

    k = pf(1);
    
    pf = polyfit(Cpt(mask),Cdt(mask),1);
    
    figure()
    plot(Cpt(mask),Cdt(mask),'LineWidth',2)
    hold on
    plot(Cpt(mask),pf(1).*Cpt(mask)+pf(2),'LineWidth',2)
    legend(["Data";"Linear Fit"])
    exportgraphics(gcf,sprintf("g_%s_%i_kdiffFit.png",fName,tMin))
    savefig(sprintf("f_%s_%i_kdiffFit.fig",fName,tMin))
    
    covMat = corrcoef(Cdt(mask),pf(1).*Cpt(mask)+pf(2));
    fprintf("k_diff R^2: %.5f\n",covMat(1,2).^2)

    k_diff = pf(1).*k;

    mu = k + k_diff;

    %Calculate distributions

    nT = reshape(Xi(end,:),[],1);
    fprintf("Number of total cells: %.1f\n",length(nT))

    nT_pd = fitdist(nT,'Kernel','Kernel','epanechnikov','Support','positive');
    
    nd = nd(nd > 0);
    fprintf("Number of dividing cells: %.1f\n",length(nd));
    if length(nd) > length(nT)
        nd = nd(end-length(nT)+1:end);
        fprintf("Number of dividing cells used: %.1f\n",length(nd))
    end

    ndiff = ndiff(ndiff > 0);
    fprintf("Number of differentiating cells: %.1f\n",length(ndiff))
    if length(ndiff) > length(nT)
        ndiff = ndiff(end-length(nT)+1:end);
        fprintf("Number of differentiating cells used: %.1f\n",length(ndiff))
    end

    nnb = nnb(nnb > 0);
    fprintf("Number of newborn cells: %.1f\n",length(nnb))
    if length(nnb) > 2*length(nT)
        nnb = nnb(end-2*length(nT)+1:end);
        fprintf("Number of newborn cells used: %.1f\n",length(nnb))
    end

    nd_pd = fitdist(nd,'Kernel','Kernel','epanechnikov','support','positive');
    ndiff_pd = fitdist(ndiff,'Kernel','Kernel','epanechnikov','support','positive');
    nnb_pd = fitdist(nnb,'Kernel','Kernel','epanechnikov','support','positive');



    %Calculate pdfs

    x = linspace(0.1*min(Xi,[],'all'),1.1.*max(Xi,[],'all'),ndist);

    nd_pdf = pdf(nd_pd,x);
    ndiff_pdf = pdf(ndiff_pd,x);
    nT_pdf = pdf(nT_pd,x);
    nnb_pdf = pdf(nnb_pd,x);

    figure()
    plot(x,nd_pdf,'LineWidth',2); hold on
    plot(x,ndiff_pdf,'LineWidth',2)
    plot(x,nT_pdf,'LineWidth',2)
    plot(x,nnb_pdf,'LineWidth',2)
    legend(["Dividing";"Differentiating";"Total";"Newborn"])
    exportgraphics(gcf,sprintf("g_%s_%i_dists.png",fName,tMin))
    savefig(sprintf("f_%s_%i_dists.fig",fName,tMin))   

    %Calculate cdfs
    nd_cdf = cdf(nd_pd,x);
    ndiff_cdf = cdf(ndiff_pd,x);
    nT_cdf = cdf(nT_pd,x);
    nnb_cdf = cdf(nnb_pd,x);


    %Calculate PSFs

    Gamma_est = mu.*nd_pdf./nT_pdf;
    D_est = k_diff.*ndiff_pdf./nT_pdf;
    R_est = mu.*(2.*nnb_cdf - nT_cdf - nd_cdf)./nT_pdf + k_diff.*(nT_cdf - ndiff_cdf)./nT_pdf;

    prctile10 = prctile(nT,10);
    prctile90 = prctile(nT,90);

    x_mask = all([x > prctile10; x < prctile90],1);

    figure()
    subplot(1,3,1)
    G_x = Gamma(x);
    plot(x,Gamma_est,'LineWidth',2); hold on; plot(x,G_x,'LineWidth',2);
    ylim([min(G_x).*0.9 max(G_x).*1.1])
    xline(prctile10)
    xline(prctile90)
    subplot(1,3,2)
    D_x = D(x);
    plot(x,D_est,'LineWidth',2); hold on; plot(x,D_x,'LineWidth',2);
    ylim([min(D_x)*0.9 max(D_x).*1.1])
    xline(prctile10)
    xline(prctile90)
    subplot(1,3,3)
    R_x = R(x);
    plot(x,R_est,'LineWidth',2); hold on; plot(x,R_x,'LineWidth',2);
    ylim([min(R_x).*0.9 max(R_x).*1.1])
    legend(["Estimated";"Original"])
    xline(prctile10,'HandleVisibility','off')
    xline(prctile90,'HandleVisibility','off')
    exportgraphics(gcf,sprintf("g_%s_%i_PSFs.png",fName,tMin))
    savefig(sprintf("f_%s_%i_PSFs.fig",fName,tMin))


    fprintf("mu: %.5f\n",mu)
    fprintf("k: %.5f\n",k)
    fprintf("kdiff: %.5f\n",k_diff)

    MSE_Gamma = mean((G_x(x_mask) - Gamma_est(x_mask)).^2);
    MSE_D = mean((D_x(x_mask) - D_est(x_mask)).^2);
    MSE_R = mean((R_x(x_mask) - R_est(x_mask)).^2);

    fprintf("# of samples in 10-90: %i\n",length(x(x_mask)))
    fprintf("Gamma MSE: %.5e\n",MSE_Gamma)
    fprintf("D MSE: %.5e\n",MSE_D)
    fprintf("R MSE: %.5e\n",MSE_R)

    fprintf("10th percentile: %.3f 90th percentile: %.3f\n",prctile10,prctile90)
    fprintf("Gamma_10: %.5f Gamma_90: %.5f\n",Gamma(prctile10),Gamma(prctile90))
    fprintf("D_10: %.5f D_90: %.5f\n",D(prctile10),D(prctile90))       
    fprintf("R_10: %.5f R_90: %.5f\n",R(prctile10),R(prctile90))


end