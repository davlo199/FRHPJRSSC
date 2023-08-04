%%Estimates the parameters RepNum number of times
function CONS = cons(NSIM, parallel, RepNum, nBeta, Nrand, alpha, beta, gamma, mu, M0, MMAX, lambda, theta,Tprime,M,c)
    ltime = tic;

    % Get output as array rather than combined vector to allow for paralelisation
    % Each colums represent a run.
    numBatches = ceil(RepNum / parallel.batchSize);
    finalParametersMatrix = zeros(5 * parallel.batchSize, numBatches);
    
    % What parralel method to use
    %   0 = None
    %   1 = Parpool
    %   2 = Slurm Jobs
    switch parallel.inner
            % Using parpool
        case 1
            parpool(nesiSafeParpool(RepNum));

            parfor nBatch = 1:numBatches
                % in case where RepNum does not evenly divide by batchSize.
                %                 if batchJJ == numBatches
                %                     thisBatchSize=RepNum-(batchJJ*batchSize)+1;
                %                 end
                finalParametersMatrix(:, nBatch) = consReps(nBatch, parallel.batchSize, nBeta, Nrand, NSIM, alpha, beta, gamma, mu, M0, MMAX, lambda, theta,Tprime,M,c);
            end

            % Using multiple SlurmJobs.
        case 2
            % functionHandle is a partial function of consReps with all static variables bound in. 
            functionHandle = @(nBatch)(consReps(nBatch, parallel.batchSize, nBeta, Nrand, NSIM, alpha, beta, gamma, mu, M0, MMAX, lambda, theta,Tprime,M,c));
            finalParametersMatrix = slurmParallel(functionHandle, 1:numBatches, parallel.slurmParameters);
        otherwise
            % Serial is same as doing all in one batch.
            parallel.batchSize=RepNum;
            finalParametersMatrix = consReps(1, RepNum, nBeta, Nrand, NSIM, alpha, beta, gamma, mu, M0, MMAX, lambda, theta,Tprime,M,c);
    end

    BestParam = [NSIM; reshape(finalParametersMatrix, [], 1)];
    %Returns the finished vector with all the best estimates for each RepNum
    %simulation of size NSIM
    fprintf("Sim %g, %g repeats finished in %.2fs\n", NSIM, RepNum, toc(ltime));
    CONS = BestParam;
    fprintf("check sum of best param: %20.15e\n", sum(CONS))

end
