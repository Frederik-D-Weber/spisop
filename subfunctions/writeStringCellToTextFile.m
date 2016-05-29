function writeStringCellToTextFile(fid,format,x)
[rows, cols] = size(x);
for i=1:rows
    if cols > 1
        fprintf(fid,['%' format ','],x{i,1:end-1});
        fprintf(fid,['%' format '\n'],x{i,end});
    else
        fprintf(fid,['%' format '\n'],x{i,1:end-1});
    end
    
end
end