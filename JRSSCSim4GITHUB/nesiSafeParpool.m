function pc = nesiSafeParpool(nWorkers)
    % nesiSafeParpool.m
    % Just creates a 'cluster' object appropriate for running inside slurm job.
    % Should work on cluster, or local, however doesn't do anything essential. Might delay startup by few seconds on local.
    % Inputs:
    %   nWorkers: Number of 'optimal' workers. e.g. max before some would have nothing to do. Does nothing here except throw a warning if mismatched. Feel free to spoof, i just couldn't be bothered varargining
    % Outputs:
    %   pc : parcluster object, can be used to start parpool.

    % To stop annoying TZ error
    setenv('TZ', 'Pacific/Auckland');

    % Delete any existing pools.

    delete(gcp('nocreate'));
    pc = parcluster('local');
    tmpdir = getenv('TMPDIR');
    slurm_cpus = str2num(getenv('SLURM_CPUS_PER_TASK'));

    if tmpdir
        pc.JobStorageLocation = tmpdir;
    end

    if slurm_cpus
        pc.NumThreads = str2num(getenv('SLURM_CPUS_PER_TASK'));
        pc.NumWorkers = str2num(getenv('SLURM_CPUS_PER_TASK')) / 2;
    end

    if nWorkers ~= pc.NumWorkers
        warning("Optimal number of workers would be '%g' ('%g' being used).", nWorkers, pc.NumWorkers)
    end

end
