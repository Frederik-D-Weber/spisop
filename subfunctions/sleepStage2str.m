function st = sleepStage2str(st)
if (st == 0)
    st = 'Wake';
elseif (st == 1)
    st = 'S1';
elseif (st == 2)
    st = 'S2';
elseif (st == 3)
    st = 'S3';
elseif (st == 4)
    st = 'S4';
elseif (st == 5)
    st = 'REM';
elseif (st == 8)
    st = 'MT';
else
    st = '???';
end
