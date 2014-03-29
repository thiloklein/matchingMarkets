function [alphadraws,betadraws,deltadraws,postmean,poststd,tstat,inOctave]=smm(selection,NTU,binary,offsetOut,offsetSel,marketFE,censored,dropOnes,repeatRun,niter,data)

% TK Gibbssampling  
%clear;

% -----------------------------------------------------------------------------
% Preliminaries.
% -----------------------------------------------------------------------------
if exist ('OCTAVE_VERSION', 'builtin') > 0 % checks if we are in octave
    isOctave = true;
    pkg load io;
    pkg load statistics;
else
    isOctave = false;
end

rand('state',111);
randn('state',333);

% -----------------------------------------------------------------------------
% Data.
% -----------------------------------------------------------------------------
% D{t}(G)=1 if group G is observed in market t, =0 otherwise. (columnvector)
% here: D{t}(1:2)=1, D{t}(3:end)=0.
% R{t}(G)=1 if group G in market t has a late payment, =0 otherwise. (columnvector)
% W{t}(G,:) contains all characteristics of group G in market t. 
% here: W{t}(1:2,:) are characteristics of the two observed groups in market t.
% X{t}(G,:) contains all observed characteristics of group G in market t.

%% For Transferable Utility game only:
% P{t} index of other group composed of residual borrowers in two-group markets.
% Example with 6 feasible groups AB CD AC AD BD BC: P{t} = 2 1 5 6 3 4.

%% For fixed effects in selection equation only:
% combs{t} matrix with borrower ids of all feasible groups in market t.

% -----------------------------------------------------------------------------
% Arguments.
% -----------------------------------------------------------------------------
%selection = 1; % logical: if TRUE estimate structural model with selection and outcome equation; if FALSE use outcome equation only.
%NTU       = 1; % logical: if TRUE use non-transferable utility (NTU) matching game; if FALSE use transferable utility (TU) matching game.
%binary    = 1; % logical: if TRUE outcome variable
%offsetOut = 0; % vector of integers indicating the indices of columns in X for which coefficients should be forced to 1. Use 0 for none.
%offsetSel = 0; % vector of integers indicating the indices of columns in W for which coefficients should be forced to 1. Use 0 for none.
%marketFE  = 0; % logical: if TRUE use market-level fixed effects in outcome equation; if FALSE don't.
%censored  = 0; % ? delta is 0:not censored, 1:from below, 2:from above
%dropOnes  = 0; % one-group-markets exluded for estimation
%repeatRun = 0; % repeated run
%niter     = 800000; % iterations

% -----------------------------------------------------------------------------
% Read data.
% -----------------------------------------------------------------------------
load(data);

if repeatRun == false
    D = struct2cell(D);
    R = struct2cell(R);
    W = struct2cell(W);
    X = struct2cell(X);
    P = struct2cell(P);
    combs = struct2cell(combs);

    % -----------------------------------------------------------------------------
    % Market identifiers.
    % -----------------------------------------------------------------------------
    T = 1:length(R); % market identifiers.
    One = [];
    Two = [];
    for t=T
        if length(D{t})==1
            One = [One t]; % one-group markets.
        else
            Two = [Two t]; % two-group markets.
        end
    end

    % -----------------------------------------------------------------------------
    % One-group-markets.
    % -----------------------------------------------------------------------------
    if dropOnes == true
        One = []; % drop one-group markets
        T = Two;
    else
        for t=One
            X{t} = [X{t} 1];
        end
        for t=Two
            X{t} = [X{t} zeros(2,1)];
        end
    end

    % -----------------------------------------------------------------------------
    % Fixed effects.
    % -----------------------------------------------------------------------------
    if marketFE == true
        if binary == true
            count1 = 0;
            for t=Two
                if R{t}(1) ~= R{t}(2) % both groups in market have different outcome.
                    count1 = count1+1;
                end
            end
            count2 = 0;
            for t=Two
                X{t} = [X{t} repmat(0,2,count1)];
                if R{t}(1) ~= R{t}(2) % both groups in market have different outcome.
                    X{t}(:,end-count2) = ones(2,1);
                    count2 = count2+1;
                end
            end
            for t=One
                X{t} = [X{t} repmat(0,1,count1)];
            end
        else % continuous outcome variable.
            count2 = 0;
            for t=Two
                X{t} = [X{t} repmat(0,2,length(Two))];
                X{t}(:,end-count2) = ones(2,1);
                count2 = count2+1;
            end
            for t=One
                X{t} = [X{t} repmat(0,1,length(Two))];
            end
        end
    end

    % -----------------------------------------------------------------------------
    % Offset outcome equation.
    % -----------------------------------------------------------------------------
    if sum(offsetOut) ~= 0
        for t=T % set offset and drop column from X
            offOut{t} = sum(X{t}(:,offsetOut),2);
            X{t}(:,offsetOut) = [];
        end
    else
        for t=T % set offset to zero
            offOut{t} = zeros(size(X{t},1),1);
        end
    end

    % -----------------------------------------------------------------------------
    % Offset selection equation.
    % -----------------------------------------------------------------------------
    if sum(offsetSel) ~= 0
        for t=Two % set offset and drop column from W
            offSel{t} = sum(W{t}(:,offsetSel),2);
            W{t}(:,offsetSel) = [];
        end
    else
        for t=Two % set offset to zero
            offSel{t} = zeros(size(W{t},1),1);
        end
    end

    % -----------------------------------------------------------------------------
    % Outcome variable.
    % -----------------------------------------------------------------------------
    if binary == true
        for t=One
            if R{t}==0
                Yupperbar{t} = -offOut{t};
                Ylowerbar{t} = -inf;
            else
                Ylowerbar{t} = -offOut{t};
                Yupperbar{t} = inf;
            end
        end
        for t=Two
            for G=1:2
                if R{t}(G)==0
                    Yupperbar{t}(G) = -offOut{t}(G);
                    Ylowerbar{t}(G) = -inf;
                else
                    Yupperbar{t}(G) = inf;
                    Ylowerbar{t}(G) = -offOut{t}(G);
                end
            end
        end
    else
        for t=T
            Y{t} = -offOut{t} + R{t};
        end
    end

    % -----------------------------------------------------------------------------
    % Size of feasible groups.
    % -----------------------------------------------------------------------------
    Uneq = [];
    for t=Two
        l{t} = length(D{t});
        p{t} = (l{t}./2)+1; % cut point between "ego" and "other" groups
        if min(min(combs{t})) == 0 % unequal group sizes in market.
            Uneq = [Uneq t];
        end
    end

    % -----------------------------------------------------------------------------
    % Indicator for group size difference.
    % -----------------------------------------------------------------------------
    if NTU == true
        s = size(W{1},2);
        for t=Two
            W{t} = [W{t}, [zeros(l{t},length(Uneq))]];
        end
        count = 0;
        for t=Uneq
            count = count + 1;
            W{t}(:,s+count) = [1, 0, ones(1,p{t}-2), zeros(1,p{t}-2)]';
        end
    end

    % -----------------------------------------------------------------------------
    % Priors.
    % -----------------------------------------------------------------------------
    % latent group outcome and valuation
    for t=T
        if binary == true
            Y{t} = zeros(size(R{t},1),1); % outcome
        end
    end
    for t=Two
        V{t} = zeros(size(D{t},1),1); % valuation
    end

    % alpha.
    alphabar = zeros(size(W{1}(1,:),2),1); 
    N = size(vertcat(W{1:end}),1);
    sigmabaralphainverse = vertcat(W{1:end})'*vertcat(W{1:end}) ./ N;
    alphabaroversigma = sigmabaralphainverse*alphabar;

    % beta.
    betabar = zeros(size(X{1}(1,:),2),1); 
    n = size(vertcat(X{1:end}),1);
    sigmabarbetainverse = vertcat(X{1:end})'*vertcat(X{1:end}) ./ n;
    betabaroversigma = sigmabarbetainverse*betabar;

    % delta.
    deltabar = 0;
    sigmabarsquaredelta = 10;

    % -----------------------------------------------------------------------------
    % Prior modes for parameters.
    % -----------------------------------------------------------------------------
    alpha = alphabar;
    beta = betabar;
    delta = deltabar;
    a = 2;
    b = 1;
    sigmasquarexiinverse = b*(a-1); % mode of G(a,b) pdf
    sigmasquarexi = 1./sigmasquarexiinverse;
end

% -----------------------------------------------------------------------------
% Matrices for parameter draws.
% -----------------------------------------------------------------------------
alphadraws = repmat(NaN,length(alpha),niter);
betadraws = repmat(NaN,length(beta),niter);
deltadraws = repmat(NaN,length(delta),niter);
etadraws = repmat(NaN,2*length(Two),niter);

if isOctave == true
    t0 = clock ();
else
    timerID = tic;  %% Start a clock and return the timer ID
end

for iter=1:niter
    iter

    %% Simulate the latent outcomes conditional on the latent group valuations, data, and parameters.
    if binary == true
        for t=One % one-group-market identifiers.        
            Yhat = X{t}*beta;
            uppercdf = normcdf(Yupperbar{t},Yhat,1);
            lowercdf = normcdf(Ylowerbar{t},Yhat,1);
            Y{t} = norminv(lowercdf+(uppercdf-lowercdf).*rand,Yhat,1); % Inverse c.d.f. sampling.
        end
        if selection == true
            for t=Two % two-group-market identifiers.
                for G=1:2 % size(R{t},1)
                    Yhat = X{t}(G,:)*beta + (V{t}(G) - W{t}(G,:)*alpha).*delta;        
                    uppercdf = normcdf(Yupperbar{t}(G),Yhat,1);
                    lowercdf = normcdf(Ylowerbar{t}(G),Yhat,1);
                    Y{t}(G) = norminv(lowercdf+(uppercdf-lowercdf).*rand,Yhat,1); % Inv c.d.f. sampling.
                end
            end
        else
            for t=Two % two-group-market identifiers.
                for G=1:2 % size(R{t},1)
                    Yhat = X{t}(G,:)*beta;
                    uppercdf = normcdf(Yupperbar{t}(G),Yhat,1);
                    lowercdf = normcdf(Ylowerbar{t}(G),Yhat,1);
                    Y{t}(G) = norminv(lowercdf+(uppercdf-lowercdf).*rand,Yhat,1); % Inv c.d.f. sampling.
                end
            end
        end
    end

    if selection == true

        %% Simulate the latent valuations conditional on the latent group outcomes, data, and parameters.

        if NTU == true % non-transferable utility
            for t=Two
                %% considering upper bounds for all unobserved groups in market t.
                Vupperbar = max(V{t}(1:2));
                Vhat = W{t}(3:end,:)*alpha;
                uppercdf = normcdf(Vupperbar,Vhat,1);
                V{t}(3:end) = norminv(uppercdf.*rand([l{t}-2,1]),Vhat,1); % Inv c.d.f. sampling.

                %% considering lower bounds for observed groups in market t.
                for G =1:2 % groups A and B.
                    Vlowerbar = max( V{t}( 1-D{t} & V{t}>V{t}(3-G) ) );
                    Vhat = W{t}(G,:)*alpha + ((Y{t}(G) - X{t}(G,:)*beta).*delta)./(sigmasquarexi + delta.^2);
                    if size(Vlowerbar,1) ~= 0 % to avoid: Empty matrix: 0-by-1.
                        lowercdf = normcdf(Vlowerbar,Vhat, sqrt(sigmasquarexi./(sigmasquarexi+delta.^2)));
                        V{t}(G) = norminv( lowercdf + (1-lowercdf).*rand, Vhat, sqrt(sigmasquarexi./(sigmasquarexi+delta.^2))); % Inv c.d.f. sampling.
                    else % no bounds in this case.
                        V{t}(G) = norminv(rand, Vhat, sqrt(sigmasquarexi./(sigmasquarexi+delta.^2))); % Inv c.d.f. sampling.
                    end
                end
            end

        else % transferable utility
            for t=Two
                %% considering upper bounds for all unobserved groups in market t.
                eqsum = sum( V{t}(1:2) + offSel{t}(1:2) );

                % "ego" groups
                Vupperbar = eqsum - (V{t}(P{t}(3:p{t})) + offSel{t}(P{t}(3:p{t}))) - offSel{t}(3:p{t});
                Vhat = W{t}(3:p{t},:)*alpha;
                uppercdf = normcdf(Vupperbar,Vhat,1);
                V{t}(3:p{t}) = norminv(uppercdf.*rand(p{t}-2,1),Vhat,1); % Inv c.d.f. sampling.

                % "other" groups
                Vupperbar = eqsum - (V{t}(P{t}(p{t}+1:end)) + offSel{t}(P{t}(p{t}+1:end))) - offSel{t}(p{t}+1:end);
                Vhat = W{t}(p{t}+1:end,:)*alpha;
                uppercdf = normcdf(Vupperbar,Vhat,1);
                V{t}(p{t}+1:end) = norminv(uppercdf.*rand(p{t}-2,1),Vhat,1); % Inv c.d.f. sampling.
    
                %% considering lower bounds for observed groups in market t.
                neqmaxsum = max( V{t}(3:end) + V{t}(P{t}(3:end)) + offSel{t}(3:end) + offSel{t}(P{t}(3:end)) );
                for G=1:2
                    Vlowerbar = neqmaxsum - (V{t}(3-G) + offSel{t}(3-G)) - offSel{t}(G);            
                    Vhat = W{t}(G,:)*alpha + ((Y{t}(G) - X{t}(G,:)*beta).*delta)./(sigmasquarexi + delta.^2);            
                    lowercdf = normcdf(Vlowerbar,Vhat,sqrt(sigmasquarexi./(sigmasquarexi+delta.^2)));
                    V{t}(G) = norminv( lowercdf + (1-lowercdf).*rand, Vhat, sqrt(sigmasquarexi./(sigmasquarexi+delta.^2))); % Inv c.d.f. sampling.
                end        
            end
        end

    end
    
    %% Simulate each group of parameters conditional on latents, data, and all other parameters.
    %% beta.
    sum1 = 0;
    sum2 = 0;
    for t=One % market identifiers.
        sum1 = sum1 + X{t}'*X{t};
        sum2 = sum2 + X{t}'*Y{t};
    end
    if selection == true
        for t=Two % market identifiers.
            sum1 = sum1 + X{t}'*X{t};
            sum2 = sum2 + X{t}'*(Y{t} - delta.*(V{t}(1:2)-W{t}(1:2,:)*alpha) );
        end
    else
        for t=Two % market identifiers.
            sum1 = sum1 + X{t}'*X{t};
            sum2 = sum2 + X{t}'*Y{t};
        end
    end
    sigmahatbeta = inv(sigmabarbetainverse + (1./sigmasquarexi).*sum1);
    sigmahatbeta = tril(sigmahatbeta)' +  tril(sigmahatbeta, -1); % matrix is symmetric.
    betahat = -sigmahatbeta * (-betabaroversigma - (1./sigmasquarexi).*sum2);
    beta = mvnrnd(betahat',sigmahatbeta)';

    if selection == true    
        %% alpha.
        sum1 = 0;
        sum2 = 0;
        for t=Two % two-group-market identifiers.
            sum1 = sum1 + W{t}'*W{t} + (delta^2./sigmasquarexi).*W{t}(1:2,:)'*W{t}(1:2,:);
            sum2 = sum2 - W{t}'*V{t} + (delta./sigmasquarexi).*(W{t}(1:2,:)'*(Y{t} - X{t}*beta - V{t}(1:2).*delta));
        end
        sigmahatalpha = inv(sigmabaralphainverse + sum1);
        sigmahatalpha = tril(sigmahatalpha)' +  tril(sigmahatalpha, -1); % matrix is symmetric.
        alphahat = -sigmahatalpha * (-alphabaroversigma + sum2);
        alpha = mvnrnd(alphahat',sigmahatalpha)';

        %% delta.
        sum1 = 0;
        sum2 = 0;
        for t=Two % two-group-market identifiers.
            sum1 = sum1 + sum( (V{t}(1:2) - W{t}(1:2,:)*alpha).^2 );
            sum2 = sum2 + sum( (Y{t} - X{t}*beta) .* (V{t}(1:2) - W{t}(1:2,:)*alpha) );
        end
        sigmahatsquaredelta = 1./(1./sigmabarsquaredelta + (1./sigmasquarexi).*sum1);
        deltahat = -sigmahatsquaredelta.*(-deltabar./sigmabarsquaredelta - (1./sigmasquarexi).*sum2);
        if censored == 1 % from below, positive covariation of residuals
            lowercdf = normcdf(0,deltahat,sqrt(sigmahatsquaredelta));
            delta = norminv(lowercdf+(1-lowercdf).*rand,deltahat,sqrt(sigmahatsquaredelta)); % Inv c.d.f. sampling.
        elseif censored == 2 % from above, negative covariation of residuals
            uppercdf = normcdf(0,deltahat,sqrt(sigmahatsquaredelta));
            delta = norminv(uppercdf.*rand,deltahat,sqrt(sigmahatsquaredelta)); % Inv c.d.f. sampling.
        else % not censored
            delta = norminv(rand,deltahat,sqrt(sigmahatsquaredelta)); % Inv c.d.f. sampling.
        end
    end

    if binary == false
        %% sigma xi.
        sum1 = 0;
        for t=One % one-group-market identifiers
                sum1 = sum1 + sum( (Y{t} - X{t}*beta).^2 );
        end
        if selection == true
            for t=Two % two-group-market identifiers
                sum1 = sum1 + sum( (Y{t} - X{t}*beta - delta.*(V{t}(1:2) - W{t}(1:2,:)*alpha)).^2 );
            end
        else
            for t=Two % two-group-market identifiers
                sum1 = sum1 + sum( (Y{t} - X{t}*beta).^2 );
            end
        end
        ahat = a + n./2;
        bhat = 1./(1./b + sum1./2);
        sigmasquarexiinverse = gamrnd(ahat,bhat);
        sigmasquarexi = (1./sigmasquarexiinverse);
    end

    %% Save this iteration's draws.
    betadraws(:,iter) = beta;
    if selection == true
        alphadraws(:,iter) = alpha;
        deltadraws(:,iter) = delta;
        eta = [];
        for t=Two
            eta = [eta (V{t}(1:2)' - (W{t}(1:2,:)*alpha)')];
        end
        etadraws(:,iter) = eta;
    end

    if isOctave == true
        if iter == niter
            % The last half of all draws are used in approximating the posterior means and the posterior standard deviations.
            startiter = ceil(iter./2);
            if selection == true
                postmean = mean([alphadraws(:,startiter:iter);betadraws(:,startiter:iter);deltadraws(:,startiter:iter)],2);
                poststd = std([alphadraws(:,startiter:iter);betadraws(:,startiter:iter);deltadraws(:,startiter:iter)],0,2);
                for i=1:length(Two)*2
                    eta(i) = mean(etadraws(i,startiter:iter));
                end
            else
                postmean = mean(betadraws(:,startiter:iter),2);
                poststd = std(betadraws(:,startiter:iter),0,2);
            end
            tstat = postmean./poststd;
    
            elapsed = etime (clock (), t0)./(60.*60)

        %    filestr = sprintf('TK_gibbsiter_%d.mat',iter);
        %    save(filestr);
        end
    else
        if iter == niter | toc(timerID)./(60.*60) > 11.9
            % The last half of all draws are used in approximating the posterior means and the posterior standard deviations.
            startiter = ceil(iter./2);
            if selection == true
                postmean = mean([alphadraws(:,startiter:iter);betadraws(:,startiter:iter);deltadraws(:,startiter:iter)],2);
                poststd = std([alphadraws(:,startiter:iter);betadraws(:,startiter:iter);deltadraws(:,startiter:iter)],0,2);
                for i=1:length(Two)*2
                    eta(i) = mean(etadraws(i,startiter:iter));
                end
            else
                postmean = mean(betadraws(:,startiter:iter),2);
                poststd = std(betadraws(:,startiter:iter),0,2);
            end
            tstat = postmean./poststd;
    
            elapsed = toc(timerID)./(60.*60)

        %    filestr = sprintf('TK_gibbsiter_%d.mat',iter);
        %    save(filestr);
            break;
        end
    end

end

