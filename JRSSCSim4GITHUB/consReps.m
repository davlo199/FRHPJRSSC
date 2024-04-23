function FinalParameters = consReps(nBatch, batchSize, nBeta, Nrand, NSIM, alpha, beta, gamma, mu, M0, MMAX, lambda, theta,Tprime,M,c)
    FinalParameters = NaN(10 * batchSize, 1); %MLE and SE 
    
        % Allow for batching.
        for run = 1:batchSize
        % For use as seed.
        nRun=(nBatch-1)*batchSize + run;
        rng(nRun)
        %Simulates data
        SimData = mainSim(NSIM, alpha, beta, gamma, mu, M0, MMAX, lambda, theta,c);
        %simulate data with c=1
        Events = SimData(1, :);
        MAG = SimData(2, :);


            D = Events;

            %% Minimisation
            
            options=optimoptions('fminunc','MaxFunctionEvaluations',5000);
            MLLik=zeros(Nrand,11);
            for k=1:Nrand
            
            x0=2.*rand(1,5)-1; 
            %Minimisation of MLEfracEXP which is the objective function (negative conditional log-likelihood of SFHP)
            [MLE,Llhood,~,~,~,Hessian]=fminunc(@(x)MLEFracEXP(x,D,MAG,M0,Tprime,M),x0,options);  
            err=sqrt(diag(inv(Hessian))).';
            MLLik(k,:)=[MLE,err,Llhood];
            end

            Likelihood = MLLik(:, end);
            MINlik = min(Likelihood);
            idx = Likelihood == MINlik;
            % Find best mle for given beta
            if sum(idx)>1
            AA = MLLik(idx,:);
            BestMLE = AA(1, :); %Vector of MLE and the likelihood
            else
            BestMLE=MLLik(idx,:);
            end
%Transforms MLE back to normal scale from log and logit
        BestMLE(2)=exp(BestMLE(2));
        BestMLE(3)=exp(BestMLE(3));
        BestMLE(5)=exp(BestMLE(5));
        BestMLE(1)=exp(BestMLE(1))/(1+exp(BestMLE(1)));
        BestMLE(4)=exp(BestMLE(4))/(1+exp(BestMLE(4)));

        FinalParameters(10 * run - 9:10 * run, 1) = BestMLE(1:10);
        %returns maximum likelihood estimates for 5 parameters
        end
        
    end

