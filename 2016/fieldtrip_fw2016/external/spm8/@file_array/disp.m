function disp(obj)
% Display a file_array object
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: disp.m 2781 2011-02-03 10:48:53Z roboos $


if numel(struct(obj))>1,
    fprintf('       %s object: ', class(obj));
    sz = size(obj);
    if length(sz)>4,
        fprintf('%d-D\n',length(sz));
    else
        for i=1:(length(sz)-1),
            fprintf('%d-by-',sz(i));
        end;
        fprintf('%d\n',sz(end));
    end;
else
    display(mystruct(obj))
end;
return;
%=======================================================================

%=======================================================================
function t = mystruct(obj)
fn = fieldnames(obj);
for i=1:length(fn)
    t.(fn{i}) = subsref(obj,struct('type','.','subs',fn{i}));
end;
return;
%=======================================================================
