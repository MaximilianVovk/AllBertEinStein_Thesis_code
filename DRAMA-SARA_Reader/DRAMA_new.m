clc
clear all
close all

%%
fileroot='C:\Users\maxiv\Documents\AllBertEinStein\';
files= dir(fullfile([append(fileroot,'SARA\REENTRY\output')],'*.txt'));
%files=dir(fullfile('C:\Users\maxiv\Documents\DRAMA\TEST folder','*.txt'));

AeroThermalHistory={};
Trajectory={};
AeroThermalHistory_files=0;
Trajectory_files=0;
matrix_holder=[];

for ii=1:length(files)
    filename=append(files(ii).folder,'/',files(ii).name);
    fid = fopen( filename );
    lineRead = fgets(fid);
    
    if ~isempty(regexp(files(ii).name,'\w*_Trajectory.txt','match'))
        name_Holder = regexp(files(ii).name,'sara.([\w\-]+)','tokens');
        AeroThermalHistory_files=AeroThermalHistory_files+1;
        while ischar(lineRead)
                if~isempty(regexp(lineRead,'(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)','tokens'))
                    cell_holder=regexp(lineRead,'(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)','tokens');
                    if string(cell_holder{1}(1))~='#'
                        matrix_holder=[matrix_holder; (str2double(cell_holder{1}(:)))'];
                    end
                end
            lineRead = fgets(fid);
        end
        AeroThermalHistory{AeroThermalHistory_files,1} = matrix_holder;
        AeroThermalHistory{AeroThermalHistory_files,2} = name_Holder{1};
        matrix_holder=[];
    elseif ~isempty(regexp(files(ii).name,'\w*_AeroThermalHistory.txt','match'))
        name_Holder = regexp(files(ii).name,'sara.([\w\-]+)','tokens');
        Trajectory_files=Trajectory_files+1;
        while ischar(lineRead)
            if~isempty(regexp(lineRead,'(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)','tokens'))
                cell_holder=regexp(lineRead,'(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)','tokens');
                if string(cell_holder{1}(1))~='#'
                    matrix_holder=[matrix_holder; (str2double(cell_holder{1}(:)))'];
                end
            end
            lineRead = fgets(fid);
        end
        Trajectory{Trajectory_files,1} = matrix_holder;
        Trajectory{Trajectory_files,2} = name_Holder{1};
        matrix_holder=[];
    end
end

clearvars -except AeroThermalHistory Trajectory