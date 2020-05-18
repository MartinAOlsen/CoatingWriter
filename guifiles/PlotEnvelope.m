function [] = PlotEnvelope(refine)
% This function makes an envelope plot of resulting data. 
% If the refine input is marked, all selected points will be rerun with

% the resolution specified. 

%% Params
slopePower = 2;  % The bias of the algorithm to find closer points. If 1 only slope counts - if inf all increases in I will give a new point

%% Load data and handle errors
try
    load('CW_GUI_DATA.mat')
catch ME
  fprintf(2,'Could not find a CW result file.\n This sould be called CW_GUI_DATA.mat and is created after a manual CW run\n Please address any problems at https://github.com/MartinAOlsen/CoatingWriter\n\n') 
  rethrow(ME)
end


%% Unpack AllPars variable
co = 0;
for i = 1:numel(History)
    for j = 1:length(History(i).AllParameterList)
        co = co +1;
        AllPars(co) = History(i).AllParameterList(j);
    end
end

%% Figure out if this data has been refined before
if exist('RefinedList') == 0
    RefinedList = zeros(length(CriteriaList(:,1)),1);
end

if refine ==  0
   refined_done = 1;
else
    refined_done = 0;
end


%% Try to find a brilliance reference
thispath = cd;
cd .. 
cd brilliance_refference
if exist('brilliance1.mat')>0
    UseBT = 1;
    tmp = iLoad('brilliance1/Div2d_sample_B.dat');
    ref = sum(sum(tmp.data));
%     tmp = iLoad('~/Dropbox/CWArticle_analysis/brilliance_th_05/Div2d_sample_B.dat')
%     ref = sum(sum(tmp.data))

    
  
else
    ref = 1;
end

cd(thispath)

%% For debugging
CriteriaList_fieldnames = fieldnames(History(1).results);

%% Choose points
% Find start point(lowest price)
[m,I] = min(CriteriaList(:,2));
PointsList(1) = I;
minPrice = m;


while true
    co = 0;
    for i = 1:length(CriteriaList(:,1)) 
       if (CriteriaList(i,2) > minPrice) && i ~= PointsList(end)
           co = co +1;
           dP = CriteriaList(i,2) - CriteriaList(PointsList(end),2);
           dI = CriteriaList(i,1) - CriteriaList(PointsList(end),1);
           this_slope(i) = dI / (dP^slopePower);
       else 
           this_slope(i) = -1;
       end
       
    end
   if (max(this_slope)< -0.05) 
      break 
   end
   [m,I] = max(this_slope);
   if CriteriaList(I,2) < max(CriteriaList(:,2)*1)
        PointsList(end+1) = I;
        minPrice = CriteriaList(I,2);
   else 
      break 
   end
end
% determine if all points have been refined before
for i = 1:length(PointsList)
   if RefinedList(PointsList(i)) > 0
       refined_done = 0;
       break
   end
end


%% Plot
close all;
figure(1)
scatter(CriteriaList(:,2),CriteriaList(:,1)/ref,5,[0.75,0.75,0.75],'filled')
hold on
errorbar(CriteriaList(PointsList,2),CriteriaList(PointsList,1)/ref,CriteriaList(PointsList,17)/ref)
axis([0,1.1*CriteriaList(PointsList(end),2),0,1.1*CriteriaList(PointsList(end),1)/ref])
xlabel('Price [k€]')
if ref == 1
    ylabel('Intensity []')
else
    ylabel('Brilliance Transfer')
end
if refined_done == 0 
   title('Will be further refined..') 
else
    
end
print(gcf,'plot_Envelope.png','-dpng','-r300')


pause(0.25)
%% Refine
% Not implemented yet
co = 0;
for i = 1:(length(PointsList))
    if RefinedList(PointsList(i)) == 0
        History(1).options.ncount = refine;
        res = runMcStas(AllPars(PointsList(i)),History(1).options,History(1).options.fNames,History(1).options.optimizefile,PointsList(i));
        CriteriaList(PointsList(i),:) = [res{:,2}];
        co = co +1;
        RefinedList(PointsList(i)) = 1;
    end
end
if co == 0
    refined_done  = 1;
end

if exist('Analysis') == 0
    Analysis = zeros(1,2);
end

save('CW_GUI_DATA.mat','Analysis','CriteriaList','History','ParameterList','RefinedList')

if refined_done == 0 
    PlotEnvelope(refine)
else
    title('Envelope plot of result')
end

