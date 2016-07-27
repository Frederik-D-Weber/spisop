function [path] = find_true_fullpath()
ctfroot
if isdeployed
    if (ispc)
        [s, r] = system('set PATH');
        path = char(regexpi(r, 'Path=(.*?);', 'tokens', 'once'));
    elseif (isunix)
        [s, r] = system('echo $PATH');
        path = char(regexpi(r, 'Path=(.*?):', 'tokens', 'once'));
    end
else
    [path,dummy_name,temp_ext]  = fileparts(mfilename('fullpath'));
end
end