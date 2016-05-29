function writeStringCellToTextFileColumsAsLines(fid,format,x)
[rows, cols] = size(x);
for i=1:cols
	fprintf(fid,['%' format '\n'],x{i});
end
end