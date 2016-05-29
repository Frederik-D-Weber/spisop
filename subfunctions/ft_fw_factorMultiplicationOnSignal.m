function data = ft_fw_factorMultiplicationOnSignal(data,fieldname,factor)
for iParts = 1:length(data.(fieldname))
    data.(fieldname){iParts} = data.(fieldname){iParts} .* factor;
end
end
