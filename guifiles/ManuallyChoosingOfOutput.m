function out = ManuallyChoosingOfOutput(History,p,options,fNames,filename);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

options.dir = [options.dir '_analysis'];

%% Create lists
CriteriaList = [];
CriteriaList_names = fieldnames(History(1).results(1));
ParameterList = [];
resultFNames=fieldnames(History(1).results(1));
counter = 0;

for i = 1 : length(History)
    for j = 1 : length(History(1).criteria)
        counter = counter +1;
        %% Create Criteria lists:
        for k = 1 : length(resultFNames)
            CriteriaList(counter,k) = eval(['History(i).results(j).' resultFNames{k} ';']); 
        end
        %% Create Parameter lists:
        tmp = History(i).Pars{j};
        for k = 1 : length(History(1).Pars{1})
            ParameterList(counter,k) = tmp(k);
        end
    end
end
% close all;
History(1).options = options;
save('CW_GUI_DATA.mat','History','CriteriaList','ParameterList');
app = CWGUI;
out= 1;
end


function Minitors = runAnalysis(History,p,options,fNames,filename,testpars)
%% Run an analysis
testcounter = 1;
options.ncount = 1e7;
for i = 1:length(testpars)
    eval(['p.' fNames.variables{i} ' = ''' num2str(testpars(i)) ''';']);
end

thisPars = p;
ParString='';
for j = 1 : length(fNames.all)
    eval(['ParString = [ ParString '' ' fNames.all{j} '=' eval(['p.' fNames.all{j}]) '''];'] ); 
end


runMcStas(ParString,options,fNames,options.analyzefile,testcounter);

Monitors(testcounter).lambda = iLoad(['/home/molsen/Dropbox/CW_Optimize_Project/PGE/PGE1_analysis' '_' num2str(options.generation) '_' num2str(testcounter) '/Lmon_sample_B.dat']);
Monitors(testcounter).x_divx = iLoad(['/home/molsen/Dropbox/CW_Optimize_Project/PGE/PGE1_analysis' '_' num2str(options.generation) '_' num2str(testcounter) '/acceptance_x_divx.dat']);
Monitors(testcounter).y_divy = iLoad(['/home/molsen/Dropbox/CW_Optimize_Project/PGE/PGE1_analysis' '_' num2str(options.generation) '_' num2str(testcounter) '/acceptance_y_divy.dat']);
Monitors(testcounter).div2d = iLoad(['/home/molsen/Dropbox/CW_Optimize_Project/PGE/PGE1_analysis' '_' num2str(options.generation) '_' num2str(testcounter) '/Div2d_sample.dat']);
Monitors(testcounter).hdiv = iLoad(['/home/molsen/Dropbox/CW_Optimize_Project/PGE/PGE1_analysis' '_' num2str(options.generation) '_' num2str(testcounter) '/Hdiv_sample.dat']);
Monitors(testcounter).vdiv = iLoad(['/home/molsen/Dropbox/CW_Optimize_Project/PGE/PGE1_analysis' '_' num2str(options.generation) '_' num2str(testcounter) '/Vdiv_sample.dat']);
Monitors(testcounter).psd = iLoad(['/home/molsen/Dropbox/CW_Optimize_Project/PGE/PGE1_analysis' '_' num2str(options.generation) '_' num2str(testcounter) '/PSD_sample.dat']);
% Monitors(testcounter).lambda = iLoad(['/home/molsen/Dropbox/CW_Optimize_Project/PGE/PGE1_analysis' '_' num2str(options.generation) '_' num2str(testcounter) '/Lmon_sample_B.dat']);
% Monitors(testcounter).lambda = iLoad(['/home/molsen/Dropbox/CW_Optimize_Project/PGE/PGE1_analysis' '_' num2str(options.generation) '_' num2str(testcounter) '/Lmon_sample_B.dat']);
% Monitors(testcounter).lambda = iLoad(['/home/molsen/Dropbox/CW_Optimize_Project/PGE/PGE1_analysis' '_' num2str(options.generation) '_' num2str(testcounter) '/Lmon_sample_B.dat']);
% Monitors(testcounter).lambda = iLoad(['/home/molsen/Dropbox/CW_Optimize_Project/PGE/PGE1_analysis' '_' num2str(options.generation) '_' num2str(testcounter) '/Lmon_sample_B.dat']);
testcounter = testcounter + 1;

end

