clear FRIPONstations

for ii=1:max(size(StationsFRIPON))
    FRIPONstations(ii).Station_Code=StationsFRIPON(ii,1);
    FRIPONstations(ii).lat=str2double(regexp(StationsFRIPON(ii,2),'([0-9.]+)','match'));
    FRIPONstations(ii).lon=str2double(regexp(StationsFRIPON(ii,3),'([0-9.]+)','match'));
    FRIPONstations(ii).alt=str2double(regexp(StationsFRIPON(ii,4),'([0-9.]+)','match'));
end

save('FRIPONstations.mat','FRIPONstations')