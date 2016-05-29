% function cellarray = read_mixed_csv(fileName,delimiter)
%     cellarray = dataset2cell(dataset('File',fileName,'Delimiter',delimiter,'ReadVarNames',0));
%     cellarray = cellarray(2:end,:);
% end

 function lineArray = read_mixed_csv(fileName,delimiter,varargin)
 includecomments = false;
 if length(varargin) == 1
     if strcmp(varargin{1},'includecomments')
         includecomments = true;
     else
     error('varargin unkown or not well defined') 
     end
 end
  fid = fopen(fileName,'r');   %# Open the file
  lineArray = cell(1000,1);     %# Preallocate a cell array (ideally slightly
                               %#   larger than is needed)
  lineIndex = 1;               %# Index of cell to place the next line in
  nextLine = fgetl(fid);       %# Read the first line from the file
  while ~isequal(nextLine,-1) %# Loop while not at the end of the file
      if strcmp(nextLine,'')
          
      elseif nextLine(1) == '#' && ~includecomments
          
      else
          lineArray{lineIndex} = nextLine;  %# Add the line to the cell array
          lineIndex = lineIndex+1;
      end
      %# Increment the line index
      nextLine = fgetl(fid);            %# Read the next line from the file
  end
  fclose(fid);                 %# Close the file
  lineArray = lineArray(1:lineIndex-1);  %# Remove empty cells, if needed
  for iLine = 1:lineIndex-1              %# Loop over lines
     lineData = textscan(lineArray{iLine},'%s',...  %# Read strings
                        'Delimiter',delimiter);
    lineData = lineData{1};              %# Remove cell encapsulation
    if strcmp(lineArray{iLine}(end),delimiter)  %# Account for when the line
      lineData{end+1} = '';                     %#   ends with a delimiter
    end
    lineArray(iLine,1:numel(lineData)) = lineData;  %# Overwrite line data
  end
end