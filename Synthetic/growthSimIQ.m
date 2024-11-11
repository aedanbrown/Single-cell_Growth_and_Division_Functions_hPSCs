rng('shuffle')

fName = "30000_200_r3";
N_p = 30000;
tauV = 1:1:200;
n_dist = 50000;
pd0 = makedist('Normal','mu',50,'sigma',0.5);
pd0 = truncate(pd0,0,100);


%PSF functions
Gamma = @(x) (x./100).^5;
R = @(x) (50 + x.^2)./(100 + x.^2) - 5e-3.*x;

[X_i,C_p_t,n_d,n_nb,tauV] = IQ_growth(N_p,tauV,pd0,Gamma,R,n_dist);


writematrix(X_i,sprintf('Xi_%s.csv',fName));
writematrix(C_p_t,sprintf('Cpt_%s.csv',fName));
writematrix(n_d,sprintf('nd_%s.csv',fName));
writematrix(n_nb,sprintf('nnb_%s.csv',fName));
writematrix(tauV',sprintf('tauV_%s.csv',fName));



function [X_i,C_p_t,n_d,n_nb,tauV] = IQ_growth(N_p,tauV,pd0,Gamma,R,n_dist)


    %tauV must have intervals larger than T
    C_p = N_p; %Initial pluripotent cell concentration


    %Distributions
    n_d = zeros(n_dist,1);
    n_di = 1;
    n_nb = zeros(2.*n_dist,1);
    n_nbi = 1;

    %Timing related variables
    tau_0 = 0; %Initial time
    tauV = [0 tauV]; %Time values to record state at
    tau_i = 2; %Index to keep track of which times have been recorded
    tauMax = max(tauV);
    n_t = length(tauV);

    X_i = zeros(n_t,N_p); %Array to store states
    C_p_t = zeros(n_t,1);
    C_p_t(1) = C_p;

    %Generate initial cell count distribution
    x_i_0 = pd0.icdf(rand(N_p,1)); %Randomly sample to obtain cell counts
    X_i(1,:) = x_i_0;



    %Initial IQ
    [T_0_cdf,x_i_T,~] = calculate_T_growth(tau_0,x_i_0,R,Gamma);


    if tau_0 + T_0_cdf > tauV(tau_i)
        %Overshot time of interest, so need to go back and reintegrate to 
        %get x_i at correct time
        %Division happens at tau_0 + T_0_cdf, which is after tauV(tau_i),
        %so don't need to worry about division
 
        dz = 10.^(ceil(log10(T_0_cdf)) - 1); %Calculate dz

        z = linspace(tau_0,tauV(tau_i),ceil((tauV(tau_i)-tau_0)/dz))';

        x_i_temp = x_i(z,x_i_0,R);
            %Calculate x_i at tauV(tau_i)
        X_i(tau_i,:) = x_i_temp(end,:); %Store x_i

        C_p_t(tau_i) = C_p;

        tau_i = tau_i + 1;


    end

    x_i_tau = x_i_T;

    mother = cdf_sample_disc(@(j) mother_cdf(j,x_i_tau,Gamma),[1 N_p]);
        %Select mother cell for divions

    x_mother = x_i_tau(mother); %Find mother cell content
    
    if n_di <= n_dist
        n_d(n_di) = x_mother;
        n_di = n_di + 1;
    else
        n_d(randi(n_dist)) = x_mother;
    end

    x_daughter_1 = cdf_sample_cont(@(y) partion_cdf(y,x_mother),x_mother/2);
        %Determine daughter 1 cell content using partion dist
    x_daughter_2 = x_mother - x_daughter_1;
        %Determine daughter 2 cell content using mass balance
    
    if n_nbi < 2*n_dist
        n_nb(n_nbi) = x_daughter_1;
        n_nbi = n_nbi + 1;
        n_nb(n_nbi) = x_daughter_2;
        n_nbi = n_nbi + 1;
    else
        n_nb(randi(2*n_dist)) = x_daughter_1;
        n_nb(randi(2*n_dist)) = x_daughter_2;
    end

    x_i_tau(mother) = x_daughter_1;
        %Replace mother

    kick_out = randi(N_p); %Pick a random cell to replace

    x_i_tau(kick_out) = x_daughter_2;
        %Replace the random cell with the second daughter cell




    C_p = (N_p + 1)./N_p.*C_p;
   
    if tau_0 + T_0_cdf == tauV(tau_i) %If we hit the division time, record current state

        X_i(tau_i,:) = x_i_tau;

        C_p_t(tau_i) = C_p;

        tau_i = tau_i + 1;

    end
    
    
    tau = tau_0 + T_0_cdf; %Update tau

    

    while tau < tauMax

        %Calculate length of IQ
        [T,x_i_T,~] = calculate_T_growth(tau,x_i_tau,R,Gamma); %Calculate IQ length


        if tau + T > tauV(tau_i)
            %Overshot time of interest, so need to go back and reintegrate to 
            %get x_i at correct time
            %Division happens at tau + T, which is after tauV(tau_i),
            %so don't need to worry about division

            dz = 10.^(ceil(log10(T)) - 1); %Calculate dz
            z = linspace(tau,tauV(tau_i),ceil((tauV(tau_i)-tau)/dz))';

            x_i_temp = x_i(z,x_i_tau,R);
                %Calculate x_i at tauV(tau_i)
            X_i(tau_i,:) = x_i_temp(end,:); %Store x_i

            C_p_t(tau_i) = C_p;

            tau_i = tau_i + 1;

    
        end

        x_i_tau = x_i_T; %Update cell contents for time change

        mother = cdf_sample_disc(@(j) mother_cdf(j,x_i_tau,Gamma),[1 N_p]);
            %Select mother cell for divions

        x_mother = x_i_tau(mother); %Find mother cell content

        if n_di <= n_dist
            n_d(n_di) = x_mother;
            n_di = n_di + 1;
        else
            n_d(randi(n_dist)) = x_mother;
        end

        x_daughter_1 = cdf_sample_cont(@(y) partion_cdf(y,x_mother),x_mother/2);
            %Determine daughter 1 cell content using partion dist
        x_daughter_2 = x_mother - x_daughter_1;
            %Determine daughter 2 cell content using mass balance

        if n_nbi < 2*n_dist
            n_nb(n_nbi) = x_daughter_1;
            n_nbi = n_nbi + 1;
            n_nb(n_nbi) = x_daughter_2;
            n_nbi = n_nbi + 1;
        else
            n_nb(randi(2*n_dist)) = x_daughter_1;
            n_nb(randi(2*n_dist)) = x_daughter_2;
        end
    
        x_i_tau(mother) = x_daughter_1;
            %Replace mother
    
        kick_out = randi(N_p); %Pick a random cell to replace
    
        x_i_tau(kick_out) = x_daughter_2;
            %Replace the random cell with the second daughter cell
    
        C_p = (N_p + 1)./N_p.*C_p; %Add another cell to account for division
       
        if tau_i < n_t %on the last time, there will be overshoot, increasing tau_i, so need to check index
            if tau + T == tauV(tau_i) %If we hit the division time, record current state

                X_i(tau_i,:) = x_i_tau;

                C_p_t(tau_i) = C_p;

                tau_i = tau_i + 1;
        
            end
        end

        tau = tau + T; %Update time

    end
    

    
end



function [M] = mother_cdf(j,x_i,Gamma)
    
    M = zeros(length(j),1);
    for i = 1:length(j)
        M(i) = sum(Gamma(x_i(1:j(i))))./sum(Gamma(x_i));
    end

end

function [p] = partion_pdf(y,y_)

    q = 40;
    pdf = @(x,x_) (1./x_).*(gamma(q+q)./(gamma(q).*gamma(q))).*((x./x_).^(q-1)).*((1 - x./x_).^(q-1));
        %Baseline pdf
    p = pdf(y,y_);
        %evaluate pdf

end


function [P] = partion_cdf(y,y_)
    
    pdf_y = @(x) partion_pdf(x,y_);

    if y >= 0 && y <= y_
        P = integral(pdf_y,0,y); %Value of cdf (from definition) as y is the value of interest
    elseif y < 0
        %P = pdf_y(0).*y;
        P = y.^3;
    elseif y > y_
        %P = pdf_y(y).*y + 1;
        P = (y - y_).^3 + 1;
    end

end



function [T_,x_i_T_,u] = calculate_T_growth(tau,x_i_tau,R,Gamma)

    u = rand(1,1);

    Ti = 0;
    [Qi,dQdTi,x_i_T_] = Q_dQdT(Ti,u,tau,x_i_tau,R,Gamma);

    while abs(Qi) > 1e-6


        if dQdTi ~= 0 %Confirm derivatie is nonzero
           
            Tip1 = Ti - Qi./dQdTi;

            [Qip1,dQdTip1,x_i_T_] = Q_dQdT(Tip1,u,tau,x_i_tau,R,Gamma);

            Ti = Tip1;
            Qi = Qip1;
            dQdTi = dQdTip1;

        else
            
            %Adjust xi to move away from the dfdx = 0
            if Ti ~= 0
                Ti = Ti.*0.99;
                [Qi,dQdTi,x_i_T_] = Q_dQdT(Ti,u,tau,x_i_tau,R,Gamma);
            else
                Ti = Ti + 1e-4; % "plus" biases towards positive roots
                [Qi,dQdTi,x_i_T_] = Q_dQdT(Ti,u,tau,x_i_tau,R,Gamma);
            end

        end
       

    end


    T_ = Ti;

 
    
end


function [q,dqdt,x_i_T] = Q_dQdT(T,u,tau,x_i_0,R,Gamma)

    if T >= 0
        dz = 10.^(ceil(log10(T)) - 1);
        x_i_t = x_i(linspace(tau,tau+T,ceil(T/dz))',x_i_0,R);
    

        a = alpha(x_i_t,Gamma);

        q = integrate_trap(a, ...
            linspace(0,T,ceil(T/dz))')./log(1-u) + 1;

    
        x_i_T = x_i_t(end,:);
    
        dqdt = a(end)./log(1-u);

    elseif T < 0

        dqdt = alpha(x_i_0',Gamma)./log(1-u);
            %Alpha expects a row vector
        q = T.*dqdt + 1;
        x_i_T = []; %Maybe should just be x_i_0? Maybe it doesn't matter?

    end

end


function [a] = alpha(x_i,Gamma)
    
    a = sum(Gamma(x_i),2);

end


function [x_i_] = x_i(t,x_i_0,R)

    dydt = @(t,x) R(x);

    x_i_ = euler_explicit(dydt,t,x_i_0);


end

function [s,u] = cdf_sample_cont(cdf,x0)


    %cdf_sample samples from a pdf given the cdf. It does this sample by
    %sampling from a uniform distribution on [0,1], and then solving for
    %s = cdf^-1(u)
    %See https://www.mathworks.com/help/stats/generate-random-numbers-using-the-uniform-distribution-inversion-method.html
    %cdf is a function describing the cumulative distribution function of
    %the probability function of interest


    u = rand(1,1); %Generate a random number from the uniform distribution on [0,1]

    f = @(x) cdf(x) - u;

    s = newtonsMethod(f,x0);

end


function [s] = cdf_sample_disc(cdf,range)

    %cdf_sample samples from a discrete pdf given the cdf. It does this 
    %sample by sampling from a uniform distribution on [0,1], and then 
    %finding the minimum value of x s.t. cdf > u
    %See https://www.cse.psu.edu/~rtc12/CSE586/lectures/cse586samplingPreMCMC.pdf
    %cdf is a function describing the cumulative distribution function of
    %the probability function of interest
    %range is the range of values the cdf can take on

    
    u = rand(1,1); %Generate a random number from the uniform distribution 
                   %on [0,1]
    

   %Implementing binary search
    %See https://cp-algorithms.com/num_methods/binary_search.html#lower-bound-and-upper-bound
    %for algorithm
   ran = min(range):max(range);
   
   l = 0;
   r = length(ran)+1;


   while (r - l > 1)

       m = floor((r + l)/2);

       c = cdf(ran(m)); %Evaluate function here to save on function evaluations

       if u <= c

           r = m;

       elseif u > c

            l = m;

       end

   end

   s = r;

end


function [y] = euler_explicit(dydt,tV,y0)
    
    %dydt returns the derivative of each component and is dydt = dydt(t,y)
    %t is a vector of time values that the solution will be found at
    %y0 is the initial condition: y0 = y(t(1))
    %Everything associated with y is a COLUMN vector

    y = zeros(length(tV),length(y0)); %Set up y for storage
    y_i = y0; %Initial y value
    t_i = tV(1); %Initial time
    y(1,:) = y0; %Store initial solution

    for i = 2:length(tV)
        
        t_ip1 = tV(i); %Cycle over times
        
        y_ip1 = y_i + dydt(t_i,y_i).*(t_ip1 - t_i);
            %Integrate using Euler's method
        
        y(i,:) = y_ip1; %Store result
        
        y_i = y_ip1; %Reset values for next step
        t_i = t_ip1;

    end


end


function [x_] = newtonsMethod(f,x0,varargin)


    dfdx_sup = ~isempty(varargin);
    if dfdx_sup
        dfdx = varargin{1};
    end

    xi = x0;
    fi = f(x0);

    while abs(fi) > 1e-6
        
        if dfdx_sup
            dfdxi = dfdx(xi);
        else
            dfdxi = dfdx_num(f,fi,xi);
        end
        
%         fprintf('%.3e %.3e %.3e\n',xi,fi,dfdxi)

        if dfdxi ~= 0 %Confirm derivatie is nonzero
           
            xip1 = xi - fi./dfdxi;
            fip1 = f(xip1);
    
            xi = xip1;
            fi = fip1;

        else
            
            %Adjust xi to move away from the dfdx = 0
            if xi ~= 0
                xi = xi.*1.01;
                fi = f(xi);
            else
                xi = xi + 1e-4; % "plus" biases towards positive roots
                fi = f(xi);
            end

        end



    end
    

    x_ = xi;    
    
end


function [dfdx_] = dfdx_num(f,f_x,x)

    
    dx = 1e-6;
    dfdx_ = (f(x+dx) - f_x)./dx;

end


function [I] =  integrate_trap(f,x)
    
    %f can either be a function handle that takes in a vector, or a array
    %of data corresponding to f(x)


    if class(f) == "double"
        
        f_v = f;

    elseif class(f) == "function_handle"

        f_v = f(x);
        
    end
    
    I = sum((x(2:end) - x(1:end-1)).*(f_v(2:end)+f_v(1:end-1))./2);

end

