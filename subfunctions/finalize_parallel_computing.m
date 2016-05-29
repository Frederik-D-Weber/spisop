function finalize_parallel_computing()

if matlabpool('size') > 0 % checking to see if my pool is already open
     matlabpool close
end

end
