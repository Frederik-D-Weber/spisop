function st = sleepStage2str_alt(st)
if (st == 0)
    st = 'Wake';
elseif (st == 1)
    st = 'S1';
elseif (st == 2)
    st = 'S2';
elseif (st == 3)
    st = 'SWS';
elseif (st == 4)
    st = 'SWS';
elseif (st == 5)
    st = 'REM';
elseif (st == 8)
    st = 'MT';
else
    st = '???';
end