function [allREMsNb funName remEpochs]= allNightAllREMsDensFromExactVct(hip, epochLength, exactREMsVct, sampling)


funName = 'allNightARD';


hipLen = length(hip);


if hipLen > length(exactREMsVct)/sampling/epochLength
    
   fprintf('Vector is shorter than hipnogram!!!\n');     
   hipLen = floor(length(exactREMsVct)/sampling/epochLength);
   
end


remEpochs = 0;
allREMsNb = 0;

for i = 1:1:hipLen
        
    if hip(i) == 5
             
        placeInExtVctBeg = (i-1)*sampling*epochLength + 1;
        placeInExtVctEnd = i*sampling*epochLength;

        
        currentFragment = exactREMsVct(placeInExtVctBeg:placeInExtVctEnd);
        REMsInCurrFrag = 0;
        previousREM = 0;
        
%        currentFragment'
        
        for jj = 1:1:length(currentFragment)
        
            
            if currentFragment(jj) > 0 

                if ~previousREM
            
                    previousREM = 1;                
                    REMsInCurrFrag = REMsInCurrFrag + 1;               

                end
                
            else
                
                previousREM = 0;
                
            end            
            
        end % for jj = 1:1:length(currentFragment)
        
        allREMsNb = allREMsNb + REMsInCurrFrag;
        remEpochs = remEpochs+1;
        
    end % if hip(i) == 5       
    
end % for i = 1:1:hipLen

if remEpochs


    allREMsNb = allREMsNb / remEpochs;

else
    
    allREMsNb = nan;
    
end

end % function [...] = allNightAllREMsFromExactVct(...)
