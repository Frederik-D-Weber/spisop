%function saves scored rapid eye movements in the way, that its easy to
%read for a human. Its done in simple text file, first column is epoch
%number, second epoch is scored REMs result, but its filled only when
%scored sleep phase in the hipnogram in this place is REM (5)

function saveREMsForPeople(savePath, remDens, meanREMSpeed, meanREMWidth, hip, hipMA)



FileID = fopen(savePath, 'w');

fprintf(FileID, '%s,%s,%s,%s,%s,%s\n', 'epoch','hypnogram','hypnogram_MA','rems_count','rems_mean_speed_microVolt_per_second','rems_mean_duration_seconds');

if length(remDens) > length(hip)
   
    remDens = remDens(1:length(hip));
    
end

for i=1:1:length(remDens)
   
    fprintf(FileID, '%d,', i);
    fprintf(FileID, '%d,', hip(i));
    fprintf(FileID, '%d,', hipMA(i));
    if(hip(i) == 5)
        fprintf(FileID, '%d,', remDens(i));
        fprintf(FileID, '%f,', meanREMSpeed(i));
        fprintf(FileID, '%f', meanREMWidth(i));
    else
        fprintf(FileID, '%d,', 0);
        fprintf(FileID, '%f,', 0);
        fprintf(FileID, '%f', 0);
    end
    if i ~= length(remDens)
        fprintf(FileID, '\n');
    end
end







end