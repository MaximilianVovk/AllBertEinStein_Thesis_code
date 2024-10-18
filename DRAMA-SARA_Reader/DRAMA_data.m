clc
clear all
close all
%% Directory of outputs

fileroot='C:\Users\maxiv\Documents\AllBertEinStein\';
files= dir(fullfile([append(fileroot,'SARA\REENTRY\output')],'*.txt'));

% Date = '10:12:25 Tue 11-Dec-18';
% res = datetime(Date ,'InputFormat','HH:mm:ss eee dd-MMM-yy','Format','dd/MM/yyyy HH:mm a');
% DateNumber = datenum(t);

%% Format for each line of text
startRow = 26;
formatSpec = '%10f%16f%11f%11f%17f%17f%11f%11f%11f%11f%11f%19f%15f%[^\n\r]';

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 25);

% Specify range and delimiter
opts.DataLines = [28, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25"];
opts.SelectedVariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "VarName1", "TrimNonNumeric", true);
opts = setvaropts(opts, "VarName1", "ThousandsSeparator", ",");

%% loop to read
jj=1;
qq=1;

for ii=1:length(files)
    filename=append(files(ii).folder,'/',files(ii).name);

    if ~isempty(regexp(files(ii).name,'\w*_Trajectory.txt','match'))  
        %% Open the Aero Thermal History file.
        fileID = fopen(filename,'r');
        
        % Read columns of data according to the format.
        % This call is based on the structure of the file used to generate this code. If an error occurs for a different file, try regenerating the code from the Import Tool.
        dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

        % Close the text file.
        fclose(fileID);

        % Create output variable
        nameHolder = regexp(files(ii).name,'sara.(\w+).','tokens');
        if nameHolder{1}=="Compound_of" && jj>=2

            matrix1=Trajectory{jj-1,1};
            matrix2=[dataArray{1:end-1}];

            if matrix1(end,1)==matrix2(1,1)
                Trajectory{jj-1,1} = [matrix1(1:end-1,:);matrix2];
            elseif matrix1(1,1)==matrix2(end,1)
                Trajectory{jj-1,1} = [matrix2(1:end-1,:);matrix1];
            else
            Trajectory{jj,2} = nameHolder{1};
            Trajectory{jj,1} = [dataArray{1:end-1}];
                jj=jj+1;
            end
        else
            Trajectory{jj,2} = nameHolder{1};
            Trajectory{jj,1} = [dataArray{1:end-1}];
            jj=jj+1;
        end
    elseif ~isempty(regexp(files(ii).name,'\w*_AeroThermalHistory.txt','match'))
        %% Open the Trajectory file.
        sara = readtable(filename, opts);

        % Convert to output type
        nameHolder = regexp(files(ii).name,'sara.(\w+).','tokens'); 
        if qq>=2
            name1=string(AeroThermalHistory{qq-1,2}); 
            name2=string(nameHolder{1});

            matrix1=AeroThermalHistory{qq-1,1};
            matrix2=table2array(sara);

            if name1==name2 && matrix1(end,1)==matrix2(1,1)
                AeroThermalHistory{qq-1,1} = [matrix1(1:end-1,:);matrix2];
            elseif name1==name2 && matrix1(1,1)==matrix2(end,1)
                AeroThermalHistory{qq-1,1} = [matrix2(1:end-1,:);matrix1];
            else
                AeroThermalHistory{qq,2} = nameHolder{1};
                AeroThermalHistory{qq,1} = table2array(sara);
                qq=qq+1;
            end
        else
            AeroThermalHistory{qq,2} = nameHolder{1};
            AeroThermalHistory{qq,1} = table2array(sara);
            qq=qq+1;
        end
    end
end

clearvars -except AeroThermalHistory Trajectory