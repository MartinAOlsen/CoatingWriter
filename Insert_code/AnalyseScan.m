close all; clear all;
pause(1)
try
   rmdir('CW_analysis_Result');
end
mkdir('CW_analysis_Result');
%% load Results
load('ScanResults.mat')
instrName={};
for i=1:length(Result)
    tmp=(Result(i).instrumentName);
    instrName{i}=[tmp];
    if Result(i).demandPrice > 0
        scanName{i}=[tmp '_price_' num2str(Result(i).demandPrice)];
    else
        scanName=instrName;
    end
        

end
guides=unique(scanName);
for i=1:length(Result)
    monitor_ALLW=Result(i).Monitor;
    Result(i).data.I=(monitor_ALLW(11).data.data(:));        
    Result(i).data.Error=(monitor_ALLW(11).data.errors(:));% Error on Intensity
    Result(i).data.WaveLength=(monitor_ALLW(11).data.x(:));
end

testWL=Result(1).data.WaveLength;
for k=1:length(testWL)
    if testWL(k)<str2num(Result(1).parameterList.WaveMin)
        index_low=k;
    end
    if testWL(k)<str2num(Result(1).parameterList.WaveMax)
        index_high=k;
    end
end
for i=1:length(Result)
   mkdir(['CW_analysis_Result/' scanName{i}]) 
end


  for i=1:numel(guides)
    counter=0;
    for j=1:length(Result)
        if strcmp(Result(j).instrumentName,(guides{i})); 
            counter=counter+1;
%             Instrument_list(counter)=sum(Result{1,j});
            if isfield(Result(j),'demandPrice')
                PriceTarget_list(counter)=Result(j).demandPrice;
            end
            Price_list(counter)=Result(j).resultPrice;
            Intensity_list(counter)=mean(Result(j).data.I(index_low:index_high));
            IntensityError_list(counter)=mean(Result(j).data.Error(index_low:index_high));
%             WL_list(counter)=cell2mat(Result(6,j));
            
        end
    end
  end
    try
        for i=1:length(Result)
            load('output/brill_ref/brilliance_ref.mat');
            BrilRef_list=(monitor_ALLW_ref(11).data.data(:));
            Result(i).data.BT = Result(i).data.I ./ BrilRef_list;
        end
    end



if isfield(Result,'stepMon')
    for j=1:length(Result)
        % Extract data:
        Steps=[];
        for i=1:length(Result(j).stepMon)
            
            if strfind(Result(j).stepMon(i).Title,'Div2d_sample_B')
                if isfield(Steps,'Div2d')
                    Steps.Div2d(end+1).Data=Result(j).stepMon(i);
                    Steps.Div2d_thermal(end+1).Data=Result(j).stepMon_thermal(i);
                    Steps.Div2d_cold(end+1).Data=Result(j).stepMon_cold(i);
                    Steps.Div2d_epithermal(end+1).Data=Result(j).stepMon_epithermal(i);
                else
                    Steps.Div2d(1).Data=Result(j).stepMon(i);
                    Steps.Div2d_thermal(1).Data=Result(j).stepMon_thermal(i);
                    Steps.Div2d_cold(1).Data=Result(j).stepMon_cold(i);
                    Steps.Div2d_epithermal(1).Data=Result(j).stepMon_epithermal(i);
                end
            elseif strfind(Result(j).stepMon(i).Title,'Lmon_sample_')
                if isfield(Steps,'Lmon')
                    Steps.Lmon(end+1).Data=Result(j).stepMon(i);
                    Steps.Lmon_thermal(end+1).Data=Result(j).stepMon_thermal(i);
                    Steps.Lmon_cold(end+1).Data=Result(j).stepMon_cold(i);
                    Steps.Lmon_epithermal(end+1).Data=Result(j).stepMon_epithermal(i);
                else
                    Steps.Lmon(1).Data=Result(j).stepMon(i);
                    Steps.Lmon_thermal(1).Data=Result(j).stepMon_thermal(i);
                    Steps.Lmon_cold(1).Data=Result(j).stepMon_cold(i);
                    Steps.Lmon_epithermal(1).Data=Result(j).stepMon_epithermal(i);
                end
            end
        end
        %% I figure
        figure(8)
        for i=1:length(Steps.Lmon)
            step_I(i)=Steps.Lmon(i).Data.data.Mean;
            step_I_thermal(i)=Steps.Lmon_thermal(i).Data.data.Mean;
            step_I_epithermal(i)=Steps.Lmon_epithermal(i).Data.data.Mean;
            step_I_cold(i)=Steps.Lmon_cold(i).Data.data.Mean;
            step_dist(i)=Steps.Lmon(i).Data.data.position(3); 
        end
        [sorted,order] = sort(step_dist);
        sorted_I=(step_I(order));
        sorted_I_thermal=(step_I_epithermal(order));
        sorted_I_epithermal=(step_I_epithermal(order));
        sorted_I_cold=(step_I_cold(order));
        for i=2:length(step_I)
            lossPrSection(i-1)=100*(sorted_I(i)-sorted_I(i-1))/sorted_I(i-1);
            lossPrSection_thermal(i-1)=100*(sorted_I_thermal(i)-sorted_I_thermal(i-1))/sorted_I_thermal(i-1);
            lossPrSection_epithermal(i-1)=100*(sorted_I_epithermal(i)-sorted_I_epithermal(i-1))/sorted_I_epithermal(i-1);
            lossPrSection_cold(i-1)=100*(sorted_I_cold(i)-sorted_I_cold(i-1))/sorted_I_cold(i-1);
            SegmentName{i-1}=(tmp(i-1))
        end


        subplot(2,1,1);
        scatter(step_dist,step_I);
        title('Intensity along guide')
        xlabel('Distance [m]')
        ylabel('Intensity [arb. unit]')
        subplot(2,1,2);
        bar([lossPrSection_epithermal;lossPrSection_thermal;lossPrSection_cold]')
        set(gca,'xticklabel',SegmentName)
        title('Loss in intensity pr section')
        xlabel('Section')
        ylabel('Change in intensity [%]')
        legend('0.45 Å','2.0 Å','4.0 Å')
        print(['CW_analysis_Result/' Result(j).instrumentName '/stepwise_result'],'-dpng')
        pause(2)
        
        
        %% BT figure
        if exist('monitor_ALLW_ref')
            figure(9)
            for i=1:length(Steps.Div2d)
                step_BT(i)=100*Steps.Div2d(i).Data.data.Mean/monitor_ALLW_ref(1).data.Mean;
                step_BT_thermal(i)=100*Steps.Div2d_thermal(i).Data.data.Mean/monitor_ALLW_ref(1).data.Mean;
                step_BT_epithermal(i)=100*Steps.Div2d_epithermal(i).Data.data.Mean/monitor_ALLW_ref(1).data.Mean;
                step_BT_cold(i)=100*Steps.Div2d_cold(i).Data.data.Mean/monitor_ALLW_ref(1).data.Mean;
                step_dist(i)=Steps.Div2d(i).Data.data.position(3); 
            end
            [sorted,order] = sort(step_dist);
            sorted_BT=(step_BT(order));
            sorted_BT_thermal=(step_BT_thermal(order));
            sorted_BT_epithermal=(step_BT_epithermal(order));
            sorted_BT_cold=(step_BT_cold(order));
            for i=2:length(step_BT)
                lossPrSection(i-1)=(sorted_BT(i)-sorted_BT(i-1));
                lossPrSection_thermal(i-1)=(sorted_BT_thermal(i)-sorted_BT_thermal(i-1));
                lossPrSection_epithermal(i-1)=(sorted_BT_epithermal(i)-sorted_BT_epithermal(i-1));
                lossPrSection_cold(i-1)=(sorted_BT_cold(i)-sorted_BT_cold(i-1));
                SegmentName{i-1}=(tmp(i-1))
            end
            subplot(2,1,1);
            scatter(step_dist,step_BT);
            title('BT along guide')
            xlabel('Distance [m]')
            ylabel('BT [%]')
            subplot(2,1,2);
            bar([lossPrSection_epithermal;lossPrSection_thermal;lossPrSection_cold]')
            set(gca,'xticklabel',SegmentName)
            title('Loss in BT pr section')
            xlabel('Section')
            ylabel('Change in BT [%]')
            legend('0.45 Å','2.0 Å','4.0 Å')
            print(['CW_analysis_Result/' Result(j).instrumentName '/stepwise_BT'],'-dpng')
            pause(2)
            
            
        end

        close all
    end
end

%% Plot guide
% Using plot_guide made by Sandor 
for i=1:length(Result)
    %try
    Profile=figure;
    set(Profile, 'Position', [0 0 1000 600])
    mstruct = read_coating([scanName{i} '/CoatingWriter_output_'  num2str(Result(i).demandPrice) '.txt']);
    listing = dir(scanName{i});
    for j = 1 : length(listing) 
        if sum(strfind(listing(j).name,'geometry.dat'))
            geoFileName = listing(j).name;
        end
    end
    component = read_guide([scanName{i} '/' geoFileName]);
%     if strcmp(scanName{i} , instrName{i}) == 0
%         component = read_guide([scanName{i} '/' instrName{i} '_' num2str(Result(i).demandPrice) '1_geometry.dat']);
%     else
%         component = read_guide([scanName{i} '/' instrName{i} '1_geometry.dat']);
%     end
    plot_guide(component,mstruct, [0,6],jet,1)
    try
        sgtitle(['Guide ' Result(i).instrumentName ' optimized for price:' num2str(Result(i).demandPrice) '. Resulting price:' num2str(Result(i).resultPrice)])
    end
    print(['CW_analysis_Result/' scanName{i} '/Profile_guide_'],'-dpng')
    
    close all;
    pause(1)
    %end
end

%% BT plots for each guide
close all
for i=1:length(Result)
    try
    index_low;
    index_high;
    
    plot(Result(i).data.WaveLength,100*Result(i).data.BT)
    
    xlabel('Wavelength [Å]');
    xlabel('BT [%]');
    ylim([0,100])
    
    title(['Guide ' Result(i).instrumentName ' optimized for price: ' num2str(Result(i).demandPrice) ' k€. Resulting price: ' num2str(Result(i).resultPrice) ' k€'])
    
    hold on
    line([str2num(Result(1).parameterList.WaveMin),str2num(Result(1).parameterList.WaveMin)],[0,100],'Color','black','LineStyle',':')
    line([str2num(Result(1).parameterList.WaveMax),str2num(Result(1).parameterList.WaveMax)],[0,100],'Color','black','LineStyle',':')
    hold off
    
    print(['CW_analysis_Result/' scanName{i} '/BT_lambda_'],'-dpng')
   
    close all;
    pause(0.5)
    end
end


%% Combined BT plot
if strcmp(Result(i).mode,'pricescan')
    figure('Position',[0,0,1000,1500])
    hold on
    line([str2num(Result(1).parameterList.WaveMin),str2num(Result(1).parameterList.WaveMin)],[0,100],'Color','black','LineStyle',':')
    line([str2num(Result(1).parameterList.WaveMax),str2num(Result(1).parameterList.WaveMax)],[0,100],'Color','black','LineStyle',':')
    legendList{1}='fom lower limit';
    legendList{2}='fom upper limit';
    for i=1:length(Result)
       

        plot(Result(i).data.WaveLength,100*Result(i).data.BT)

        xlabel('Wavelength [Å]');
        ylabel('BT [%]');
        ylim([0,100])

        title('All results')
        legendList(i+2)={[instrName{i} ' ' sprintf('%2.0f',Result(i).resultPrice) ' k€']};
        
        legend(legendList);
        print(['CW_analysis_Result/BT_lambda_Compare'],'-dpng')
        pause(0.5)
        


    end
end
close all;


%% Overview BT plot
for i= 1 : length(Result)
    instrumentNames{i} = Result(i).instrumentName;
end
instrList = unique(instrumentNames);
iconList={'o','^','d','s'};
colorList={'r','g','b'};
figure()
hold on
mkrindex = 1; colorindex = 1;
h=[]; 

for i = 1 : length(instrList)
    xlist = [];
    ylist = [];
    first = 1;
    for j = 1 : length(instrumentNames)
        if strcmp(instrumentNames{j},instrList{i})
            scatter(Result(j).resultPrice,100*mean(Result(j).data.BT(index_low:index_high)),colorList{colorindex},iconList{mkrindex})    
            xlist(end+1)=Result(j).resultPrice;
            ylist(end+1)=100*mean(Result(j).data.BT(index_low:index_high));
            if first == 1; h(end+1) = scatter(0,0,colorList{colorindex},iconList{mkrindex} ,'visible', 'off') ; first = 0; end
        end

    end
    [xlist, sortIndex] = sort(xlist)
    ylist = ylist(sortIndex)
    plot(xlist,ylist,[colorList{colorindex} ':']);
    if colorindex == 3; colorindex = 0; mkrindex = mkrindex +1 ;end
    colorindex = colorindex +1;
end
xlabel('Price [k€]')
ylabel('BT [%]')
legend(h,instrList)
pause(0.5)
print(['CW_analysis_Result/BT_price_overview'],'-dpng')
pause(0.5)


close all

%% Envelope and stuff:
try
    

    
    figure(1);
    %% PRISLISTE OG INTENSITETS LISTE MANGLER!!!
    for i = 1:length(Result)
        Price_list(i) = Result(i).resultPrice;
        Intensity_list(i) = sum(Result(i).data.BT);
    end
    scatter(Price_list,100*(Intensity_list)./mean(BrilRef_list(index_low:index_high)))
    title('Brilliance transfer of optimized guides at different prices');
    ylabel('Brilliance Transfer [%]');
    xlabel('Price');
    legends='';
    %legend(guides)
    hold on

end
try
    print('CW_analysis_Result/PriceScan','-dpng')
    pause(2)
end

if isfield(Result,'critList')
    figure(2);
    hold on
    legendList={}
    for j=1:length(Result)
        for k=1:length(Result(j).critList)
            OptimResult(j,k)=Result(j).critList(k);
        end
        if cell2mat(strfind(fieldnames(Result(j)),'OptimizeMode'))
            if strcmp(Result(j).OptimizeMode,'singleParScan')
                tmp=1.5+[1:length(Result(j).critList)-2]/(length(Result(j).critList)-2);
                index=[1,tmp,3];    
                plot(index,abs(OptimResult(j,:)),'--')
            end
        else
            plot(abs(OptimResult(j,:)),'--')
        end
        ylim([0,max(abs(OptimResult(j,:)))])
        legendList(j)={[Result(j).instrumentName]};
        try
            legendList(j)={[legendList{j} '-' num2str(Result(j).demandPrice)]};
        end
    end
    title('change after each run')
    legend(legendList)
    print('CW_analysis_Result/critEvolution','-dpng')
    pause(2)
    ylabel('Intensity')
    xlabel('Optimization number')
    
    figure(3)
    for i=1:length(OptimResult(1,:))
        barResult(i)=mean(OptimResult(:,i));
        bar(abs(barResult),'DisplayName','OptimResult')
        title('average intensity after each optimation')
        ylabel('Intensity')
        xlabel('Optimization number')
    end
    
end
if isfield(Result,'valueList')
%     colorList='gbrmcky';sRGB_Values = uint32([...
colorList={[115,82,68]
[194,150,130]
[98,122,157]
[87,108,67]
[133,128,177]
[103,189,170]
[214,126,44]
[80,91,166]
[193,90,99]
[94,60,108]
[157,188,64]
[224,163,46]
[56,61,150]
[70,148,73]
[175,54,60]
[231,199,31]
[187,86,149]
[8,133,161]
[243,243,242]
[200,200,200]
[160,160,160]
[122,122,121]
[85,85,85]
[52,52,52]};
	
	for i=1:length(Result)
        figure(4)
        pause(0.5)
        hold on
        color=colorList{i}/255;
        if strcmp(Result(i).mode,'coating')
            s = scatter(Result(i).valueList(:,1),abs(Result(i).valueList(:,2)),15,'filled','MarkerFaceColor',color)
            out=maxLine(Result(i).valueList(:,1),abs(Result(i).valueList(:,2)));
            plot(out(:,1),out(:,2),'Color',color)
            pause(0.1)
        else
            s = scatter(Result(i).valueList(:,1),abs(Result(i).valueList(:,2).*Result(i).valueList(:,3)),15,'filled','MarkerFaceColor',color)
            out=maxLine(Result(i).valueList(:,1),abs(Result(i).valueList(:,2).*Result(i).valueList(:,3)));
            plot(out(:,1),out(:,2),'Color',color)
            pause(0.1)
        end
        s.LineWidth = 0.2;
        s.MarkerEdgeColor = 'k';
		xlabel('Price [k€]')
        ylabel('Intensity')
        print('CW_analysis_Result/scan_I_all','-dpng')
        pause(2)
        figure(5)
        hold on
        plot(out(:,1),out(:,2),'Color',color)
        pause(0.1)
        xlabel('Price [k€]')
        ylabel('Intensity')
        legend(Result(:).instrumentName)
        print('CW_analysis_Result/scan_I_outline','-dpng')
        pause(2)
	end
end
if isfield(Result,'valueList')
	try
	for i=1:length(Result)
        figure(6)
        hold on
        color=colorList{i}/255;
        if strcmp(Result(i).mode,'coating')
            s = scatter(Result(i).valueList(:,1),abs(Result(i).valueList(:,2)./Result(i).valueList(:,1)),15,'filled','MarkerFaceColor',color)
            out=maxLine(Result(i).valueList(:,1),abs(Result(i).valueList(:,2)./Result(i).valueList(:,1)));
            plot(out(:,1),out(:,2),'Color',color)
            pause(0.1)
        else
            s = scatter(Result(i).valueList(:,1),Result(i).valueList(:,3)*abs(Result(i).valueList(:,2))./Result(i).valueList(:,1),15,'filled','MarkerFaceColor',color)
            out=maxLine(Result(i).valueList(:,1),Result(i).valueList(:,3)*abs(Result(i).valueList(:,2))./Result(i).valueList(:,1));
            plot(out(:,1),out(:,2),'Color',color)
            pause(0.1)
        end
        s.LineWidth = 0.2;
        s.MarkerEdgeColor = 'k';
		xlabel('Price [k€]')
        ylabel('Intensity/k€')
        print('CW_analysis_Result/scan_Iprice_all','-dpng')
        pause(2)
        figure(7)
        hold on
        plot(out(:,1),out(:,2),'Color',color)
        pause(0.1)
        xlabel('Price [k€]')
        ylabel('Intensity/k€')
        legend(Result(:).instrumentName)
        print('CW_analysis_Result/scan_Iprice_outline','-dpng')
        pause(2)
    end
    end
end


close all;