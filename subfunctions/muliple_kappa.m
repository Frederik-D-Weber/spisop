function [k,sek,p,z,ci,kj,sekj,zkj,pkj,consensusModusVoteMatrix,consensusCountMatrix,cross_comparisons,chi2_cross,p_cross] = muliple_kappa(alignmentMatrix,alphalevel,referenceOption)
% alignmentMatrix = [0 1 1 1 2   2 3 3 6 1 NaN;...
%                    0 1 1 1 1   2 2 3 6 1 2  ;...
%                    0 1 1 1 NaN 2 2 3 5 1 NaN]
categories = unique(alignmentMatrix);
categories = categories(~isnan(categories));
nCategories = numel(categories);
%referenceOption = 'first'; %either 'first' or 'consenus'

nObservations = size(alignmentMatrix,2);
nSamples = size(alignmentMatrix,1);


consensusCountMatrix = zeros(nCategories,nObservations);
consensusModusVoteMatrix = zeros(1,nObservations);

for iColumn = 1:nObservations
    for iCat = 1:nCategories
        consensusCountMatrix(iCat,iColumn) = sum(alignmentMatrix(:,iColumn) == categories(iCat));
    end
    consensusModusVoteMatrix(iColumn) = categories(find(consensusCountMatrix(:,iColumn) == max(consensusCountMatrix(:,iColumn)),1,'first'));
end

[k,sek,p,z,ci,kj,sekj,zkj,pkj] = fleiss(consensusCountMatrix',alphalevel);

switch referenceOption
    case 'firstInList'
        reference = alignmentMatrix(1,:);
        startCompIndex = 2;
    case 'consenus'
        reference = consensusModusVoteMatrix;
        startCompIndex = 1;
    otherwise
end

cross_comparisons = {};
chi2_cross = [];
p_cross = [];

for iComp = startCompIndex:nSamples
    test = alignmentMatrix(iComp,:);
    [tbl,chi2_cross_temp,p_cross_temp,labels] = crosstab(test,reference);
    row_lab = labels(:,1);
    row_lab_index = cellfun(@(s) ~isempty(s),row_lab,'UniformOutput',false);
    row_lab_index = [row_lab_index{:}];
    row_lab = row_lab(row_lab_index);
    
    col_lab = labels(:,2);
    col_lab_index = cellfun(@(s) ~isempty(s),col_lab,'UniformOutput',false);
    col_lab_index = [col_lab_index{:}];
    col_lab = col_lab(col_lab_index);
    
    row_lab = cellfun(@(s) strrep(s,'-1','unscored'), row_lab,'UniformOutput',false);
    col_lab = cellfun(@(s) strrep(s,'-1','unscored'), col_lab,'UniformOutput',false);


    cross_comparisons{iComp} = array2table(tbl,'RowNames',cellfun(@(s) ['test_' s], row_lab,'UniformOutput',false),'VariableNames',cellfun(@(s) ['ref_' s], col_lab,'UniformOutput',false));
    chi2_cross(iComp) = chi2_cross_temp;
    p_cross(iComp) = p_cross_temp;
end
end


function [k,sek,p,z,ci,kj,sekj,zkj,pkj] = fleiss(varargin)
% FLEISS: compute the Fleiss'es kappa
% Fleiss'es kappa is a generalisation of Scott's pi statistic, a
% statistical measure of inter-rater reliability. It is also related to
% Cohen's kappa statistic. Whereas Scott's pi and Cohen's kappa work for
% only two raters, Fleiss'es kappa works for any number of raters giving
% categorical ratings (see nominal data), to a fixed number of items. It
% can be interpreted as expressing the extent to which the observed amount
% of agreement among raters exceeds what would be expected if all raters
% made their ratings completely randomly. Agreement can be thought of as
% follows, if a fixed number of people assign numerical ratings to a number
% of items then the kappa will give a measure for how consistent the
% ratings are. The scoring range is between 0 and 1.
%
% Syntax: 	fleiss(X,alpha)
%
%     Inputs:
%           X - square data matrix
%           ALPHA - significance level (default = 0.05)
%     Outputs:
%           - kappa value for the j-th category (kj)
%           - kj standard error
%           - z of kj
%           - p-value
%           - Fleiss'es kappa
%           - kappa standard error
%           - kappa confidence interval
%           - k benchmarks by Landis and Koch
%           - z test
%
%      Example:
%An example of the use of Fleiss'es kappa may be the following: Consider
%fourteen psychiatrists are asked to look at ten patients. Each
%psychiatrist gives one of possibly five diagnoses to each patient. The
%Fleiss' kappa can be computed from this matrix to show
%the degree of agreement between the psychiatrists above the level of
%agreement expected by chance.
%x =
%     0     0     0     0    14
%     0     2     6     4     2
%     0     0     3     5     6
%     0     3     9     2     0
%     2     2     8     1     1
%     7     7     0     0     0
%     3     2     6     3     0
%     2     5     3     2     2
%     6     5     2     1     0
%     0     2     2     3     7
%
% So there are 10 rows (1 for each patient) and 5 columns (1 for each
% diagnosis). Each cell represents the number of raters who
% assigned the i-th subject to the j-th category
% x=[0 0 0 0 14; 0 2 6 4 2; 0 0 3 5 6; 0 3 9 2 0; 2 2 8 1 1 ; 7 7 0 0 0;...
% 3 2 6 3 0; 2 5 3 2 2; 6 5 2 1 0; 0 2 2 3 7];
%
%  Calling on Matlab the function: fleiss(x);
%
%           Answer is:
%
% kj:       0.2013    0.0797    0.1716    0.0304    0.5077
%
% s.e.:     0.0331
%
% z:        6.0719    2.4034    5.1764    0.9165   15.3141
%
% p:        0.0000    0.0162    0.0000    0.3594         0
%
% ------------------------------------------------------------
% Fleiss'es (overall) kappa = 0.2099
% kappa error = 0.0170
% kappa C.I. (95%) = 0.1767 	 0.2432
% Fair agreement
% z = 12.3743 	 p = 0.0000
% Reject null hypotesis: observed agreement is not accidental
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2007) Fleisses kappa: compute the Fleiss'es kappa for multiple raters.
% http://www.mathworks.com/matlabcentral/fileexchange/15426

args=cell(varargin);
nu=numel(args);
if isempty(nu)
    error('Warning: Matrix of data is missed...')
elseif nu>2
    error('Warning: Max two input data are required')
end
default.values = {[],0.05};
default.values(1:nu) = args;
[x alpha] = deal(default.values{:});

if isvector(x)
    error('Warning: x must be a matrix, not a vector.');
end
if ~all(isfinite(x(:))) || ~all(isnumeric(x(:))) || isempty(x)
    error('Warning: X data matrix values must be numeric and finite')
end
if ~isequal(x(:),round(x(:)))
    error('Warning: X data matrix values must be whole numbers')
end
n=size(x,1); %subjects
%chech if the raters are the same for each rows
r=sum(x,2);
if any(r-max(r))
    error('The raters are not the same for each rows')
end
if nu==2 %if necessary check alpha
    if ~isscalar(alpha) || ~isnumeric(alpha) || ~isfinite(alpha) || isempty(alpha)
        error('Warning: it is required a numeric, finite and scalar ALPHA value.');
    end
    if alpha <= 0 || alpha >= 1 %check if alpha is between 0 and 1
        error('Warning: ALPHA must be comprised between 0 and 1.')
    end
end
clear args default nu

m=sum(x(1,:)); %raters
a=n*m;
pj=(sum(x)./(a)); %overall proportion of ratings in category j
b=pj.*(1-pj);
c=a*(m-1);
d=sum(b);
kj=1-(sum((x.*(m-x)))./(c.*b)); %the value of kappa for the j-th category
sekj=realsqrt(2/c); %kj standar error
zkj=kj./sekj;
pkj=(1-0.5*erfc(-abs(zkj)/realsqrt(2)))*2;
k=sum(b.*kj)/d; %Fleiss'es (overall) kappa
sek=realsqrt(2*(d^2-sum(b.*(1-2.*pj))))/sum(b.*realsqrt(c)); %kappa standard error
ci=k+([-1 1].*(abs(0.5*erfc(-alpha/2/realsqrt(2)))*sek)); %k confidence interval
z=k/sek; %normalized kappa
p=(1-0.5*erfc(-abs(z)/realsqrt(2)))*2;

% %display results
% fprintf('kj:   '); disp(kj)
% fprintf('s.e.: '); disp(sekj)
% fprintf('z:    '); disp(zkj)
% fprintf('p:    '); disp(pkj)
% disp(repmat('-',1,60))
% fprintf('Fleiss''es (overall) kappa = %0.4f\n',k)
% fprintf('kappa error = %0.4f\n',sek)
% fprintf('kappa C.I. (%d%%) = %0.4f \t %0.4f\n',(1-alpha)*100,ci)
% if k<0
%     disp('Poor agreement')
% elseif k>=0 && k<=0.2
%     disp('Slight agreement')
% elseif k>0.2 && k<=0.4
%     disp('Fair agreement')
% elseif k>0.4 && k<=0.6
%     disp('Moderate agreement')
% elseif k>0.6 && k<=0.8
%     disp('Substantial agreement')
% elseif k>0.8 && k<=1
%     disp('Perfect agreement')
% end
% fprintf('z = %0.4f \t p = %0.4f\n',z,p)
% if p<0.05
%     disp('Reject null hypotesis: observed agreement is not accidental')
% else
%     disp('Accept null hypotesis: observed agreement is accidental')
% end
end