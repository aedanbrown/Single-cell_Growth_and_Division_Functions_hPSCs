function [] = inverseProblem_synthetic_diff_multi(fName,nFiles,Gamma,R,D)

    %fName is the root file name for the files that contain the various
    %distributions, cell concentrations, and time values
    %nFiles is the number of files that will be examined
    %Gamma is a function handle for the Gamma PSF used
    %R is a function handle for the R PSF used
    %D is a function handle for the D PSF used

    Xi(:,:,nFiles) = readmatrix(sprintf('Xi_%s_r%i.csv',fName,nFiles));
    Cpt(:,:,nFiles) = readmatrix(sprintf('Cpt_%s_r%i.csv',fName,nFiles));
    Cdt(:,:,nFiles) = readmatrix(sprintf('Cdt_%s_r%i.csv',fName,nFiles));
    nd(:,:,nFiles) = readmatrix(sprintf('nd_%s_r%i.csv',fName,nFiles));
    nnb(:,:,nFiles) = readmatrix(sprintf('nnb_%s_r%i.csv',fName,nFiles));
    ndiff(:,:,nFiles) = readmatrix(sprintf('ndiff_%s_r%i.csv',fName,nFiles));
    tauV(:,:,nFiles) = readmatrix(sprintf('tauV_%s_r%i.csv',fName,nFiles));


    for i = nFiles-1:-1:1
        Xi(:,:,i) = readmatrix(sprintf('Xi_%s_r%i.csv',fName,i));
        Cpt(:,:,i) = readmatrix(sprintf('Cpt_%s_r%i.csv',fName,i));
        Cdt(:,:,i) = readmatrix(sprintf('Cdt_%s_r%i.csv',fName,i));
        nd(:,:,i) = readmatrix(sprintf('nd_%s_r%i.csv',fName,i));
        nnb(:,:,i) = readmatrix(sprintf('nnb_%s_r%i.csv',fName,i));
        ndiff(:,:,i) = readmatrix(sprintf('ndiff_%s_r%i.csv',fName,i));
        tauV(:,:,i) = readmatrix(sprintf('tauV_%s_r%i.csv',fName,i));
    end


    if any(logical(tauV - mean(tauV,3)))
        disp("Times values are not aligned. Cannot proceed.")
    else
    
        tauV = mean(tauV,3); %Only need 1 column of tauV, and they are all the same

        nXi = 50000; %Minimum number of values from Xi we want to use when calculating the distribution
        ndist = 1000;
    
        %Calculate rate constants
        figure()
        plot(tauV,log(mean(Cpt,3)),'LineWidth',2)
        hold on
        plot(tauV,log(mean(Cdt,3)),'LineWidth',2)
        legend(["ln(Cpt)";"ln(Cdt)"])
        drawnow
    
        tMin = input("When does linear region begin? ");
    
        xline(tMin,'LineWidth',2,'HandleVisibility','off');
        exportgraphics(gcf,sprintf("g_%s_%i_concs_multi.png",fName,tMin))
        savefig(sprintf("f_%s_%i_concs_multi.fig",fName,tMin))

        mask = tauV > tMin;
    
        pf = polyfit(tauV(mask),log(mean(Cpt(mask,:,:),3)),1);
    
        figure()
        plot(tauV(mask),log(mean(Cpt(mask,:,:),3)),'LineWidth',2)
        hold on
        plot(tauV(mask),pf(1).*tauV(mask)+pf(2),'LineWidth',2)
        legend(["Data";"Linear Fit"])
        exportgraphics(gcf,sprintf("g_%s_%i_kFit_multi.png",fName,tMin))
        savefig(sprintf("f_%s_%i_kFit_multi.fig",fName,tMin))

        covMat = corrcoef(log(mean(Cpt(mask,:,:),3)),pf(1).*tauV(mask)+pf(2));
        fprintf("k R^2: %.5f\n",covMat(1,2).^2)

        k = pf(1);
        
        pf = polyfit(mean(Cpt(mask,:,:),3),mean(Cdt(mask,:,:),3),1);
        
        figure()
        plot(mean(Cpt(mask,:,:),3),mean(Cdt(mask,:,:),3),'LineWidth',2)
        hold on
        plot(mean(Cpt(mask,:,:),3),pf(1).*mean(Cpt(mask,:,:),3)+pf(2),'LineWidth',2)
        legend(["Data";"Linear Fit"])
        exportgraphics(gcf,sprintf("g_%s_%i_kdiffFit_multi.png",fName,tMin))
        savefig(sprintf("f_%s_%i_kdiffFit_multi.fig",fName,tMin))

        covMat = corrcoef(mean(Cdt(mask,:,:),3),pf(1).*mean(Cpt(mask,:,:),3)+pf(2));
        fprintf("k_diff R^2: %.5f\n",covMat(1,2).^2)
    
        k_diff = pf(1).*k;
    
        mu = k + k_diff;



        %Perform calculations once for the last entry to allocate the
        %correct size, and then go backwards
    
        %Calculate distributions
        
        nT_temp = reshape(Xi(end,:,nFiles),[],1);
        fprintf("Number of total cells: %.1f\n",length(nT_temp))

        nT_pd(nFiles) = fitdist(nT_temp,'Kernel','Kernel','epanechnikov','Support','positive');

        nd_temp = nd(nd(:,:,nFiles) > 0,nFiles);
        fprintf("Number of dividing cells for data %i: %.1f\n",nFiles,length(nd_temp));
        if length(nd_temp) > length(nT_temp)
            nd_temp = nd_temp(end-length(nT_temp)+1:end);
            fprintf("Number of dividing cells for data %i used: %.1f\n",nFiles,length(nd_temp));
        end        

        ndiff_temp = ndiff(ndiff(:,:,nFiles) > 0,nFiles);
        fprintf("Number of differentiating cells for data %i: %.1f\n",nFiles,length(ndiff_temp))
        if length(ndiff_temp) > length(nT_temp)
            ndiff_temp = ndiff_temp(end-length(nT_temp)+1:end);
            fprintf("Number of differentiating cells for data %i used: %.1f\n",nFiles,length(ndiff_temp))
        end

        nnb_temp = nnb(nnb(:,:,nFiles) > 0,nFiles);
        fprintf("Number of newborn cells for data %i: %.1f\n",nFiles,length(nnb_temp))
        if length(nnb_temp) > 2*length(nT_temp)
            nnb_temp = nnb_temp(end-2*length(nT_temp)+1:end);
            fprintf("Number of newborn cells for data %i used: %.1f\n",nFiles,length(nnb_temp))
        end


        nd_pd(nFiles) = fitdist(nd_temp,'Kernel','Kernel','epanechnikov','support','positive');
        ndiff_pd(nFiles) = fitdist(ndiff_temp,'Kernel','Kernel','epanechnikov','support','positive');
        nnb_pd(nFiles) = fitdist(nnb_temp,'Kernel','Kernel','epanechnikov','support','positive');


        %Calculate pdfs
        x = linspace(0.1*min(Xi,[],'all'),1.1.*max(Xi,[],'all'),ndist);
        nx = length(x);

        nd_pdf(nFiles,1:nx) = pdf(nd_pd(nFiles),x);
        ndiff_pdf(nFiles,1:nx) = pdf(ndiff_pd(nFiles),x);
        nT_pdf(nFiles,1:nx) = pdf(nT_pd(nFiles),x);
        nnb_pdf(nFiles,1:nx) = pdf(nnb_pd(nFiles),x);


        %Calculate cdfs
        nd_cdf(nFiles,1:nx) = cdf(nd_pd(nFiles),x);
        ndiff_cdf(nFiles,1:nx) = cdf(ndiff_pd(nFiles),x);
        nT_cdf(nFiles,1:nx) = cdf(nT_pd(nFiles),x);
        nnb_cdf(nFiles,1:nx) = cdf(nnb_pd(nFiles),x);

        for i = nFiles-1:-1:1

            nT_temp = reshape(Xi(end,:,i),[],1);
            fprintf("Number of total cells: %.1f\n",length(nT_temp))
    
            nT_pd(i) = fitdist(nT_temp,'Kernel','Kernel','epanechnikov','Support','positive');
    
            nd_temp = nd(nd(:,:,i) > 0,i);
            fprintf("Number of dividing cells for data %i: %.1f\n",i,length(nd_temp));
            if length(nd_temp) > length(nT_temp)
                nd_temp = nd_temp(end-length(nT_temp)+1:end);
                fprintf("Number of dividing cells for data %i used: %.1f\n",i,length(nd_temp));
            end        
    
            ndiff_temp = ndiff(ndiff(:,:,i) > 0,i);
            fprintf("Number of differentiating cells for data %i: %.1f\n",i,length(ndiff_temp))
            if length(ndiff_temp) > length(nT_temp)
                ndiff_temp = ndiff_temp(end-length(nT_temp)+1:end);
                fprintf("Number of differentiating cells for data %i used: %.1f\n",i,length(ndiff_temp))
            end
    
            nnb_temp = nnb(nnb(:,:,i) > 0,i);
            fprintf("Number of newborn cells for data %i: %.1f\n",i,length(nnb_temp))
            if length(nnb_temp) > 2*length(nT_temp)
                nnb_temp = nnb_temp(end-2*length(nT_temp)+1:end);
                fprintf("Number of newborn cells for data %i used: %.1f\n",i,length(nnb_temp))
            end
    
    
            nd_pd(i) = fitdist(nd_temp,'Kernel','Kernel','epanechnikov','support','positive');
            ndiff_pd(i) = fitdist(ndiff_temp,'Kernel','Kernel','epanechnikov','support','positive');
            nnb_pd(i) = fitdist(nnb_temp,'Kernel','Kernel','epanechnikov','support','positive');
    
    
            %Calculate pdfs
            x = linspace(0.1*min(Xi,[],'all'),1.1.*max(Xi,[],'all'),ndist);
            nx = length(x);
    
            nd_pdf(i,1:nx) = pdf(nd_pd(i),x);
            ndiff_pdf(i,1:nx) = pdf(ndiff_pd(i),x);
            nT_pdf(i,1:nx) = pdf(nT_pd(i),x);
            nnb_pdf(i,1:nx) = pdf(nnb_pd(i),x);
    
    
            %Calculate cdfs
            nd_cdf(i,1:nx) = cdf(nd_pd(i),x);
            ndiff_cdf(i,1:nx) = cdf(ndiff_pd(i),x);
            nT_cdf(i,1:nx) = cdf(nT_pd(i),x);
            nnb_cdf(i,1:nx) = cdf(nnb_pd(i),x);
 

        end

    
        figure()
        plot(x,mean(nd_pdf,1),'LineWidth',2); hold on
        plot(x,mean(ndiff_pdf,1),'LineWidth',2)
        plot(x,mean(nT_pdf,1),'LineWidth',2)
        plot(x,mean(nnb_pdf,1),'LineWidth',2)
        legend(["Dividing";"Differentiating";"Total";"Newborn"])
        exportgraphics(gcf,sprintf("g_%s_%i_dists_multi.png",fName,tMin))
        savefig(sprintf("f_%s_%i_dists_multi.fig",fName,tMin))       
    
    
        %Calculate PSFs
        Gamma_est = mu.*mean(nd_pdf,1)./mean(nT_pdf,1);
        D_est = k_diff.*mean(ndiff_pdf,1)./mean(nT_pdf,1);
        R_est = mu.*(2.*mean(nnb_cdf,1) - mean(nT_cdf,1) - mean(nd_cdf,1))./mean(nT_pdf,1) ...
            + k_diff.*(mean(nT_cdf,1) - mean(ndiff_cdf,1))./mean(nT_pdf,1);
    

        nT_temp = reshape(Xi(end-ceil(nXi/width(Xi)):end,:,:),[],1);
        prctile10 = prctile(nT_temp,10);
        prctile90 = prctile(nT_temp,90);

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
        exportgraphics(gcf,sprintf("g_%s_%i_PSFs_multi.png",fName,tMin))
        savefig(sprintf("f_%s_%i_PSFs_multi.fig",fName,tMin))

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


end