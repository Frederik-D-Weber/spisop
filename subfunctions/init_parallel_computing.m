function init_parallel_computing(doParallelComputation,numberOfWorkers)

parallelComputation = doParallelComputation; % either 'yes' or 'no' default 'no'
numParallelWorkers = numberOfWorkers;% number of workers, use at maximum one less than the number of real cpu cores in your computer

if matlabpool('size') > 0 % checking to see if my pool is already open
     matlabpool close
end

if strcmp(parallelComputation,'yes') && (matlabpool('size') == 0) % checking to see if my pool is already open
     matlabpool('local',numParallelWorkers)
end

end