function st = sleepStage2str_alt2(st)
if (st == 0)
    st = 'Wake';
elseif (st == 1)
    st = 'S1';
elseif (st == 2)
    st = 'NonREM';
elseif (st == 3)
    st = 'NonREM';
elseif (st == 4)
    st = 'NonREM';
elseif (st == 5)
    st = 'REM';
elseif (st == 8)
    st = 'MT';
else
    st = '???';
end