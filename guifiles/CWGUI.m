classdef CWGUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                 matlab.ui.Figure
        UIAxes                   matlab.ui.control.UIAxes
        Slider                   matlab.ui.control.Slider
        Parameter1Label          matlab.ui.control.Label
        Parameter1DropDown       matlab.ui.control.DropDown
        Slider_2                 matlab.ui.control.Slider
        Slider_3                 matlab.ui.control.Slider
        Parameter3DropDownLabel  matlab.ui.control.Label
        Parameter3DropDown       matlab.ui.control.DropDown
        ParameterWeightsLabel    matlab.ui.control.Label
        Intensitywillhavetheweight1Label  matlab.ui.control.Label
        Parameter2DropDownLabel  matlab.ui.control.Label
        Parameter2DropDown       matlab.ui.control.DropDown
        UpdatePlotButton         matlab.ui.control.Button
        BestoptionLabel          matlab.ui.control.Label
        TextArea                 matlab.ui.control.TextArea
        QuickAnalyzeButton       matlab.ui.control.Button
        ChoosethisButton         matlab.ui.control.Button
        PlotaxisLabel            matlab.ui.control.Label
        xAxisDropDownLabel       matlab.ui.control.Label
        xAxisDropDown            matlab.ui.control.DropDown
        yAxisDropDown            matlab.ui.control.DropDown
        yAxisDropDownLabel       matlab.ui.control.Label
        Slider_4                 matlab.ui.control.Slider
        Parameter4DropDownLabel  matlab.ui.control.Label
        Parameter4DropDown       matlab.ui.control.DropDown
        EnvelopePlotButton matlab.ui.control.Button
        RefineList  matlab.ui.control.DropDown
    end



    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app,History,CriteriaList,ParameterList)

            % Create UIFigure
            app.UIFigure = uifigure;
            app.UIFigure.Position = [100 100 1154 669];
            app.UIFigure.Name = 'UI Figure';

            

            % Create Slider
            app.Slider = uislider(app.UIFigure);
            app.Slider.Limits = [0 1];
            app.Slider.Position = [780 551 275 3];
            app.Slider.Value = 1;

            % Create Parameter1Label
            app.Parameter1Label = uilabel(app.UIFigure);
            app.Parameter1Label.HorizontalAlignment = 'right';
            app.Parameter1Label.Position = [723 563 72 22];
            app.Parameter1Label.Text = 'Parameter 1';

            % Create Parameter1DropDown
            app.Parameter1DropDown = uidropdown(app.UIFigure);
            app.Parameter1DropDown.Enable = 'off';
            app.Parameter1DropDown.Position = [810 563 100 22];
            app.Parameter1DropDown.Items = fieldnames(History(1).results(1));
            app.Parameter1DropDown.Value = 'intensity';

            % Create Slider_2
            app.Slider_2 = uislider(app.UIFigure);
            app.Slider_2.Limits = [0 1];
            app.Slider_2.Position = [780 454 275 3];
            app.Slider_2.Value = 1;

            % Create Slider_3
            app.Slider_3 = uislider(app.UIFigure);
            app.Slider_3.Limits = [0 1];
            app.Slider_3.Position = [780 353 275 3];

            % Create Parameter3DropDownLabel
            app.Parameter3DropDownLabel = uilabel(app.UIFigure);
            app.Parameter3DropDownLabel.HorizontalAlignment = 'right';
            app.Parameter3DropDownLabel.Position = [723 368 72 22];
            app.Parameter3DropDownLabel.Text = 'Parameter 3';

            % Create Parameter3DropDown
            app.Parameter3DropDown = uidropdown(app.UIFigure);
            app.Parameter3DropDown.Position = [810 368 100 22];
            app.Parameter3DropDown.Items = fieldnames(History(1).results(1));
            app.Parameter3DropDown.Value = 'none';

            % Create ParameterWeightsLabel
            app.ParameterWeightsLabel = uilabel(app.UIFigure);
            app.ParameterWeightsLabel.FontSize = 18;
            app.ParameterWeightsLabel.FontWeight = 'bold';
            app.ParameterWeightsLabel.Position = [726 623 316 22];
            app.ParameterWeightsLabel.Text = 'Parameter Weights';

            % Create Intensitywillhavetheweight1Label
%             app.Intensitywillhavetheweight1Label = uilabel(app.UIFigure);
%             app.Intensitywillhavetheweight1Label.Position = [726 602 168 22];
%             app.Intensitywillhavetheweight1Label.Text = 'Intensity will have the weight 1';

            % Create Parameter2DropDownLabel
            app.Parameter2DropDownLabel = uilabel(app.UIFigure);
            app.Parameter2DropDownLabel.HorizontalAlignment = 'right';
            app.Parameter2DropDownLabel.Position = [723 466 72 22];
            app.Parameter2DropDownLabel.Text = 'Parameter 2';

            % Create Parameter2DropDown
            app.Parameter2DropDown = uidropdown(app.UIFigure);
            app.Parameter2DropDown.Position = [810 466 100 22];
            app.Parameter2DropDown.Items = fieldnames(History(1).results(1));
            app.Parameter2DropDown.Value = 'price';

            % Create UpdatePlotButton
            app.UpdatePlotButton = uibutton(app.UIFigure, 'push');
            app.UpdatePlotButton.Position = [961 186 100 22];
            app.UpdatePlotButton.Text = 'Update Plot';
            app.UpdatePlotButton.ButtonPushedFcn = createCallbackFcn(app,@UpdateFigure);

            % Create BestoptionLabel
            app.BestoptionLabel = uilabel(app.UIFigure);
            app.BestoptionLabel.FontWeight = 'bold';
            app.BestoptionLabel.Position = [55 165 71 22];
            app.BestoptionLabel.Text = 'Best option';

            % Create TextArea
            app.TextArea = uitextarea(app.UIFigure);
            app.TextArea.Position = [55 19 332 147];

            % Create QuickAnalyzeButton
            app.QuickAnalyzeButton = uibutton(app.UIFigure, 'push');
            app.QuickAnalyzeButton.Position = [401 19 105 22];
            app.QuickAnalyzeButton.Text = 'Quick Analyze';
            app.QuickAnalyzeButton.ButtonPushedFcn = createCallbackFcn(app,@RunQuickAnalyze);
            
            
            % Create EnvelopePlotButton
            app.EnvelopePlotButton = uibutton(app.UIFigure, 'push');
            app.EnvelopePlotButton.Position = [681 38 105 22];
            app.EnvelopePlotButton.Text = 'Plot Pareto Front';
            app.EnvelopePlotButton.ButtonPushedFcn = createCallbackFcn(app,@EnvelopeButton);
            
            % Create RefineList
            app.RefineList = uidropdown(app.UIFigure);
            app.RefineList.Position = [681 19 105 22];
            app.RefineList.Items = {'0','1e6','1e7','1e8','1e9'};
            app.RefineList.Value = '0';
            

            % Create ChoosethisButton
            app.ChoosethisButton = uibutton(app.UIFigure, 'push');
            app.ChoosethisButton.Position = [961 19 100 22];
            app.ChoosethisButton.Text = 'Choose this';

            % Create PlotaxisLabel
            app.PlotaxisLabel = uilabel(app.UIFigure);
            app.PlotaxisLabel.FontWeight = 'bold';
            app.PlotaxisLabel.Position = [401 165 55 22];
            app.PlotaxisLabel.Text = 'Plot axis';

            % Create xAxisDropDownLabel
            app.xAxisDropDownLabel = uilabel(app.UIFigure);
            app.xAxisDropDownLabel.HorizontalAlignment = 'right';
            app.xAxisDropDownLabel.Position = [418 127 48 22];
            app.xAxisDropDownLabel.Text = 'x - Axis ';

            % Create xAxisDropDown
            app.yAxisDropDown = uidropdown(app.UIFigure);
            app.yAxisDropDown.Position = [481 127 100 22];
            app.yAxisDropDown.Items = fieldnames(History(1).results(1));
            app.yAxisDropDown.Value = 'price';

            % Create yAxisDropDown
            app.xAxisDropDown = uidropdown(app.UIFigure);
            app.xAxisDropDown.Position = [481 93 100 22];
            app.xAxisDropDown.Items = fieldnames(History(1).results(1));
            app.xAxisDropDown.Value = 'intensity';

            % Create yAxisDropDownLabel
            app.yAxisDropDownLabel = uilabel(app.UIFigure);
            app.yAxisDropDownLabel.HorizontalAlignment = 'right';
            app.yAxisDropDownLabel.Position = [421 93 45 22];
            app.yAxisDropDownLabel.Text = 'y - Axis';

            % Create Slider_4
            app.Slider_4 = uislider(app.UIFigure);
            app.Slider_4.Limits = [0 1];
            app.Slider_4.Position = [780 256 275 3];

            % Create Parameter4DropDownLabel
            app.Parameter4DropDownLabel = uilabel(app.UIFigure);
            app.Parameter4DropDownLabel.HorizontalAlignment = 'right';
            app.Parameter4DropDownLabel.Position = [723 271 72 22];
            app.Parameter4DropDownLabel.Text = {'Parameter 4'; ''};

            % Create Parameter4DropDown
            app.Parameter4DropDown = uidropdown(app.UIFigure);
            app.Parameter4DropDown.Position = [810 271 100 22];
            app.Parameter4DropDown.Items = fieldnames(History(1).results(1));
            app.Parameter4DropDown.Value = 'none';
            
            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Solutions')
            xlabel(app.UIAxes, app.xAxisDropDown.Value)
            ylabel(app.UIAxes, app.yAxisDropDown.Value)
            app.UIAxes.Position = [14 186 713 469];
            
           
        end
    end

    methods (Access = public)

        % Construct app
        function app = CWGUI
            load('CW_GUI_DATA.mat');
            
            for i = 1:length(History)
                for j = 1:length(History(1).results)
                    History(i).results(j).none=0;
                end
            end
            Analysis = cell(length(History(1).results) * length(History),1);
            save('CW_GUI_DATA.mat','History','CriteriaList','ParameterList','Analysis')
            
            % Create and configure components
            createComponents(app,History,CriteriaList,ParameterList)
            
            % Update scatterplot
            %markerID = FindBestCrit(app,CriteriaList,History);
            app = UpdateFigure(app);

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
%                 clear app
            end
              %% Save some data and set values
            ids = 0:0.005:1;
            app.Slider_2.Value = 0.5;
            for i = 1 : length(ids)
                app.Slider.Value = ids(i);
                markerID = FindBestCrit(app,CriteriaList,History);
                binCrit(i,1) = CriteriaList(markerID,1);
                binCrit(i,2) = CriteriaList(markerID,2);
            end
            save('binCrits.mat','binCrit');
            
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
        
        function markerID = FindBestCrit(app,CriteriaList,History)
            
            fNames=fieldnames(History(1).results);
           
            %% Par 1
            if strcmp(app.Parameter1DropDown.Value,'none') == 0
                % This parameter index
                for i = 1:length(fNames)
                    if strcmp(app.Parameter1DropDown.Value,fNames{i})
                        thisParId = i;
                    end
                end
                
                minPar = min(CriteriaList(:,thisParId));
                maxPar = max(CriteriaList(:,thisParId));
                norm = app.Slider.Value*99 / (maxPar-minPar) ;
                for i = 1:size(CriteriaList,1)
                    WeightedArray(i) =  1 + norm * (CriteriaList(i,thisParId)-minPar);
                    
                end
            end
            
             %% Par 2
              if strcmp(app.Parameter2DropDown.Value,'none') == 0
                % This parameter index
                for i = 1:length(fNames)
                    if strcmp(app.Parameter2DropDown.Value,fNames{i})
                        thisParId = i;
                    end
                end
                
                minPar = min(CriteriaList(:,thisParId));
                maxPar = max(CriteriaList(:,thisParId));
                norm = app.Slider_2.Value*99 / (maxPar-minPar) ;
                for i = 1:size(CriteriaList,1)
                    WeightedArray(i) = WeightedArray(i) - norm * (CriteriaList(i,thisParId)-minPar);
                    
                end
              end
            
              
              %% Par 3
              if strcmp(app.Parameter3DropDown.Value,'none') == 0
                % This parameter index
                for i = 1:length(fNames)
                    if strcmp(app.Parameter3DropDown.Value,fNames{i})
                        thisParId = i;
                    end
                end
                
                minPar = min(CriteriaList(:,thisParId));
                maxPar = max(CriteriaList(:,thisParId));
                norm = app.Slider_3.Value*99 / (maxPar-minPar) ;
                for i = 1:size(CriteriaList,1)
                    WeightedArray(i) = WeightedArray(i) - norm * (CriteriaList(i,thisParId)-minPar);
                    
                end
              end
            
              %% Par 4
              if strcmp(app.Parameter4DropDown.Value,'none') == 0
                % This parameter index
                for i = 1:length(fNames)
                    if strcmp(app.Parameter4DropDown.Value,fNames{i})
                        thisParId = i;
                    end
                end
                
                minPar = min(CriteriaList(:,thisParId));
                maxPar = max(CriteriaList(:,thisParId));
                norm = app.Slider_4.Value*99 / (maxPar-minPar) ;
                for i = 1:size(CriteriaList,1)
                    WeightedArray(i) = WeightedArray(i) - norm * (CriteriaList(i,thisParId)-minPar);
                    
                end
            end
            
            
            [maxVal,markerID] = max(WeightedArray);
            
        end
        
        function app = UpdateFigure(app)
            load('CW_GUI_DATA.mat');
            markerID = FindBestCrit(app,CriteriaList,History);
            
            fnames = fieldnames(History(1).results);
            for i = 1:length(fnames)
                if strcmp(app.yAxisDropDown.Value,fnames{i})
                    x_crit_index = i;
                elseif strcmp(app.xAxisDropDown.Value,fnames{i})
                    y_crit_index = i;
                end
            end
            
            
            title(app.UIAxes, 'Solutions')
            xlabel(app.UIAxes, app.yAxisDropDown.Value)
            ylabel(app.UIAxes, app.xAxisDropDown.Value)
            app.UIAxes.Position = [14 186 713 469];
            scatter(app.UIAxes,CriteriaList(:,x_crit_index),CriteriaList(:,y_crit_index),20,[0.75,0.75,0.75],'filled');
            hold(app.UIAxes,'on')
            scatter(app.UIAxes,CriteriaList(markerID,x_crit_index),CriteriaList(markerID,y_crit_index),20,[1,0,0],'filled');
            scatter(app.UIAxes,CriteriaList(markerID,x_crit_index),CriteriaList(markerID,y_crit_index),60,[1,0,0]);
            hold(app.UIAxes,'off')
            
            %% Update text area:
            app.TextArea.Value = sprintf('Best Option (green):\n  Intensity:   %2.5f\n  Price: %2.0fk€\nBackground: \n  (PSD):\t%2.4f\t%2.0f%%\n  (DIV):\t%2.4f\t%2.0f%%\n  (WAV):\t%2.4f\t%2.0f%%\n',CriteriaList(markerID,1),CriteriaList(markerID,2),CriteriaList(markerID,9), 100*(CriteriaList(markerID,9))/CriteriaList(markerID,1)  ,CriteriaList(markerID,11),100*(CriteriaList(markerID,11))/CriteriaList(markerID,1),  CriteriaList(markerID,3),100*(CriteriaList(markerID,3))/CriteriaList(markerID,1));
            save('CW_GUI_DATA.mat','History','CriteriaList','ParameterList','Analysis');
   
        end
        
        function app = EnvelopeButton(app)
            PlotEnvelope(str2num(app.RefineList.Value));
        end
        
        
        function app = RunQuickAnalyze(app)
            load('CW_GUI_DATA.mat');
            markerID = FindBestCrit(app,CriteriaList,History);
            %% Get options and parameters
            GenSize = length(History(1).results);
            nGen = floor((markerID)/GenSize) +1;
            nThisGen = markerID -((nGen-1)*GenSize);
            p = History(nGen).AllParameterList(nThisGen);
            
            % scan all wl
            p.WaveMin='0.1';
            p.WaveMax='8';
            
            options.gravitation=1;
            options.ncount=1*1e7;
            options.mpi=2;
            try
               options.mpi=feature('numcores');
            end
            options.dir=[cd '/CW_waveALL_' num2str(markerID) ''];
            
            
            %% Make (or load) analysis
            if length(Analysis{markerID}) > 0

            else
                try; unix(['rm -rf ' options.dir]); end
                Analysis{markerID} = mcstas([History(1).options.analyzefile],p,options);
                try; copyfile('CoatingWriter_rawData.txt',['CW_waveALL_' num2str(markerID) '']); end
            end

            save('CW_GUI_DATA.mat','History','CriteriaList','ParameterList','Analysis');
            
            %% Plot analysis
            close all
            fig = figure('Position',[150,20,1800,800]);
            Monitors = Analysis{markerID};
                %% WLB
                subplot(2,4,[1 2])
                plot(Monitors(11))
                try
                    hold on
                    title('Wavelength')
                    line([str2num(History(nGen).AllParameterList(nThisGen).WaveMin),str2num(History(nGen).AllParameterList(nThisGen).WaveMin)],[0,1.3*max(Monitors(11).Data.data)],'Color','red')
                    line([str2num(History(nGen).AllParameterList(nThisGen).WaveMax),str2num(History(nGen).AllParameterList(nThisGen).WaveMax)],[0,1.3*max(Monitors(11).Data.data)],'Color','red')
                    axis tight
                    hold off
                end
                
                %% Axcep pos/pos
                subplot(2,4,[3])
                plot(Monitors(4))
                view(0,90)
                try
                    hold on
                    title('pos/pos')
                    MaxI = max(max(Monitors(13).Data.data)*1.1);
                    line([-str2num(p.sizeX)*50,str2num(p.sizeX)*50],[-str2num(p.sizeY)*50,-str2num(p.sizeY)*50],[MaxI,MaxI],'Color','red','LineWidth',2)
                    line([-str2num(p.sizeX)*50,str2num(p.sizeX)*50],[str2num(p.sizeY)*50,str2num(p.sizeY)*50],[MaxI,MaxI],'Color','red','LineWidth',2)
                    line([-str2num(p.sizeX)*50,-str2num(p.sizeX)*50],[str2num(p.sizeY)*50,-str2num(p.sizeY)*50],[MaxI,MaxI],'Color','red','LineWidth',2)
                    line([str2num(p.sizeX)*50,str2num(p.sizeX)*50],[str2num(p.sizeY)*50,-str2num(p.sizeY)*50],[MaxI,MaxI],'Color','red','LineWidth',2)
                    axis tight
                    hold off
                end
                
                %% Axcep div/div
                subplot(2,4,[4])
                plot(Monitors(12))
                view(0,90)
                try
                    hold on
                    title('div/div')
                    MaxI = max(max(Monitors(3).Data.data)*1.1);
                    line([-str2num(p.divreq_x),str2num(p.divreq_x)],[-str2num(p.divreq_x),-str2num(p.divreq_x)],[MaxI,MaxI],'Color','red','LineWidth',2)
                    line([-str2num(p.divreq_x),str2num(p.divreq_x)],[str2num(p.divreq_x),str2num(p.divreq_x)],[MaxI,MaxI],'Color','red','LineWidth',2)
                    line([-str2num(p.divreq_x),-str2num(p.divreq_x)],[str2num(p.divreq_x),-str2num(p.divreq_x)],[MaxI,MaxI],'Color','red','LineWidth',2)
                    line([str2num(p.divreq_x),str2num(p.divreq_x)],[str2num(p.divreq_x),-str2num(p.divreq_x)],[MaxI,MaxI],'Color','red','LineWidth',2)
                    axis tight
                    hold off
                end
                
                %% Axcep divx/posx
                subplot(2,4,[7])
                plot(Monitors(9))
                view(0,90)
                try
                    hold on
                    title('x_div/x_pos')
                    MaxI = max(max(Monitors(9).Data.data)*1.1);
                    line([-str2num(p.sizeX)/2,str2num(p.sizeX)/2],[-str2num(p.divreq_x),-str2num(p.divreq_x)],[MaxI,MaxI],'Color','red','LineWidth',2)
                    line([-str2num(p.sizeX)/2,str2num(p.sizeX)/2],[str2num(p.divreq_x),str2num(p.divreq_x)],[MaxI,MaxI],'Color','red','LineWidth',2)
                    line([-str2num(p.sizeX)/2,-str2num(p.sizeX)/2],[str2num(p.divreq_x),-str2num(p.divreq_x)],[MaxI,MaxI],'Color','red','LineWidth',2)
                    line([str2num(p.sizeX)/2,str2num(p.sizeX)/2],[str2num(p.divreq_x),-str2num(p.divreq_x)],[MaxI,MaxI],'Color','red','LineWidth',2)
                    axis tight
                    hold off
                end
                
                %% Axcep divy/posy
                subplot(2,4,[8])
                plot(Monitors(10))
                view(0,90)
                try
                    hold on
                    title('y_div/y_pos')
                    MaxI = max(max(Monitors(10).Data.data)*1.1);
                    line([-str2num(p.sizeY)/2,str2num(p.sizeY)/2],[-str2num(p.divreq_y),-str2num(p.divreq_y)],[MaxI,MaxI],'Color','red','LineWidth',2)
                    line([-str2num(p.sizeY)/2,str2num(p.sizeY)/2],[str2num(p.divreq_y),str2num(p.divreq_y)],[MaxI,MaxI],'Color','red','LineWidth',2)
                    line([-str2num(p.sizeY)/2,-str2num(p.sizeY)/2],[str2num(p.divreq_y),-str2num(p.divreq_y)],[MaxI,MaxI],'Color','red','LineWidth',2)
                    line([str2num(p.sizeY)/2,str2num(p.sizeY)/2],[str2num(p.divreq_y),-str2num(p.divreq_y)],[MaxI,MaxI],'Color','red','LineWidth',2)
                    axis tight
                    hold off
                end
            
                 %% INFO
                subplot(2,4,[5 6])
                Info = History(nGen).results(nThisGen);
                InfoString{1} ='';
                InfoString{end+1} = ['Guide summary:'];
                InfoString{end+1} = '';
                InfoString{end+1} = sprintf('Price:      %2.0f k€',Info.price);
                InfoString{end+1} = sprintf('Intensity: %2f',Info.intensity);
                InfoString{end+1} = '';
                InfoString{end+1} = 'Background';
                InfoString{end+1} = sprintf('(PSD):  %2.3f      %2.0f%%',CriteriaList(markerID,9), 100*(CriteriaList(markerID,9))/CriteriaList(markerID,1));
                InfoString{end+1} = sprintf('(DIV):   %2.3f     %2.0f%%',CriteriaList(markerID,11), 100*(CriteriaList(markerID,11))/CriteriaList(markerID,1));
                %InfoString{end+1} = sprintf('(LAM):  %2.3f     %2.0f%%',CriteriaList(markerID,3), 100*(CriteriaList(markerID,3))/CriteriaList(markerID,1));
                
                try
                   text( 0,0.5 ,InfoString); axis off
                end
                %sprintf('Best Option (green):\n  Intensity:   %2.5f\n  Price: %2.0fk€\nBackground: \n  (PSD):\t%2.4f\t%2.0f%%\n  (DIV):\t%2.4f\t%2.0f%%\n  (WAV):\t%2.4f\t%2.0f%%\n',CriteriaList(markerID,1),CriteriaList(markerID,2),CriteriaList(markerID,9), 100*(CriteriaList(markerID,9))/CriteriaList(markerID,1)  ,CriteriaList(markerID,11),100*(CriteriaList(markerID,11))/CriteriaList(markerID,1),  CriteriaList(markerID,3),100*(CriteriaList(markerID,3))/CriteriaList(markerID,1));
                fig2 = figure('Position',[150,20,1800,800]);
                CWfileName = (['CW_waveALL_' num2str(markerID) '/CoatingWriter_rawData.txt']);
                
                Fid = fopen(CWfileName);
                i = 1;
                tline = fgetl(Fid);
                A{i} = tline;
                while ischar(tline)
                    i = i+1;
                    tline = fgetl(Fid);
                    A{i} = tline;
                end
                fclose(Fid)
                
                numseg = 1;
                for i = length(A)-1:-1:1
                    B = textscan(A{i}, '%f,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s');
                    if B{1} > numseg
                       numseg =  B{1};
                    end
                end
               
                
                Pars = History(nGen).AllParameterList(nThisGen);
                
                
                curLen = eval(p.closest_element);
                segNr = 1;
                cMap = jet(12);
               
                for seg = numseg:-1:1
                    for i = 1:length(A)-1

                        B = textscan(A{i}, '%f,%1c,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s');
                        if B{1}==seg
                            if length(B{1}) > 0
                                len = B{3};
                                height1 = B{4};
                                height2 = B{5};
                                width1 = B{6};
                                width2 = B{7};

                                subplot(2,1,1)
                                line([curLen curLen+len],[height1/2,height2/2],'color',cMap(round(B{10}*2),:),'LineWidth',2.25)
                                line([curLen curLen+len],[-height1/2,-height2/2],'color',cMap(round(B{11}*2),:),'LineWidth',2.25)

                                subplot(2,1,2)
                                line([curLen curLen+len],[width1/2,width2/2],'color',cMap(round(B{8}*2),:),'LineWidth',2.25)
                                line([curLen curLen+len],[-width1/2,-width2/2],'color',cMap(round(B{9}*2),:),'LineWidth',2.25)

                                curLen = curLen+len;
                            end
                        end
                    end
                end
                
                subplot(2,1,1)
                title('Vertical')
                xlim([0,curLen+eval(p.sample_dist)])
                line([curLen+eval(p.sample_dist),curLen+eval(p.sample_dist)],[-eval(p.sizeY)/2,eval(p.sizeY)/2],'color',[0,0,1],'LineWidth',1.25)
                line([0,0],[-eval(p.mod_y)/2,eval(p.mod_y)/2],'color',[0,0,1],'LineWidth',1.25)
                xlabel('dist [m]')
                ylabel('height [m]')
                
                
                
                subplot(2,1,2)
                title('Horizontal')
                xlim([0,curLen+eval(p.sample_dist)])
                line([curLen+eval(p.sample_dist),curLen+eval(p.sample_dist)],[-eval(p.sizeX)/2,eval(p.sizeX)/2],'color',[0,0,1],'LineWidth',1.25)
                line([0,0],[-eval(p.mod_x)/2,eval(p.mod_x)/2],'color',[0,0,1],'LineWidth',1.25)
                xlabel('dist [m]')
                ylabel('width [m]')
                
                
                text( 0,0.205 ,'Sketch that shows coating distributions. Gaps missing and parabolic guides will fail');
        end
    end
end
