clear all;clc;close all;
cd /media/martin/HDD1/Dropbox/Dropbox/Dropbox/CoatingWriter_2.3.01b

%% 1 Insert the mode to control what kind of optimization will be made

    % coating   - Optimize brilliance transfer with coating distributions
    % value     - Optimize for best performance/price
    % budget    - Optimize for best brilliance transfer with a certain price
    % pricescan - Run several budget runs and examine BT as price increases
Mode='pricescan'



%% 2 Mode specific options

    % Coating
    if strcmp(Mode,'coating')
        CW_input.scanType='coating'
    % Value
    elseif strcmp(Mode,'value')
        CW_input.scanType='value'
        CW_input.Price=0;                           % A constant price to be added in the calculation. Use to avoid making too cheap guides. Including the total instrument price could make sense 
    % Budget
    elseif strcmp(Mode,'budget')
        CW_input.scanType='budget';
        CW_input.Price=2000;                        % Target price
        CW_input.budgetMaxError=50;                 % Maximum unpunished error in price
        CW_input.punishment='potential';            % 'shutter' cuts neutrons off, 'potential' gives a price-punishment depending of the error.   
        CW_input.potentialPower=1.5;                % IF potential is chosen as punishment, the punishment will be ERROR^potentialPower 
        CW_input.substrateMode='none';              % Either power-function dependant on length or area of the segment, or constant ('none' to ignore all substrate in prices)
        CW_input.substrateArealDependancy='linear'; % linear or power dependancy - disregarded if no areal dependancy
        CW_input.coatingMode='area';                % Either power-function dependant on length or area of the segment, or constant
        CW_input.coatingArealDependancy='SNfit1';   % linear or power dependancy  - disregarded if no areal dependancy
    % Price-Scan
    elseif strcmp(Mode,'pricescan')
        CW_input.scanType='pricescan'
        CW_input.Price=[800:100:1800];              % Target price list
        CW_input.budgetMaxError=5;                  % Maximum unpunished error in price
        CW_input.punishment='potential';            % 'shutter' cuts neutrons off, 'potential' gives a price-punishment depending of the error.   
        CW_input.potentialPower=1.5;                % IF potential is chosen as punishment, the punishment will be ERROR^potentialPower 
        CW_input.substrateMode='none';              % Either power-function dependant on length or area of the segment, or constant ('none' to ignore all substrate in prices)
        CW_input.substrateArealDependancy='linear'; % linear or power dependancy - disregarded if no areal dependancy
        CW_input.coatingMode='area';                % Either power-function dependant on length or area of the segment, or constant
        CW_input.coatingArealDependancy='SNfit1';   % linear or power dependancy  - disregarded if no areal dependancy
    end


%% 3 Segment parameters

    % Segment 1 - furthest from source, closest to sample.
    CW_input.Input1.type='E'             % Guide_bot segmet type 
    CW_input.Input1.mode='specific'      % Determines input variables. "specific" is advised in this version of CoatingWriter
    CW_input.Input1.distribution='power' % The coating distribution. 
    CW_input.Input1.mincenter=0          % center of the distribution in relative coordinates (from 0 to 1)
    CW_input.Input1.maxcenter=1
    CW_input.Input1.fixSides='HV'        % Optimize sides together. Options are: 'HV', 'curve', 'all', 'none'. ('HV' optimizes horizontal and vertical independantly)
    CW_input.Input1.deltaM=0.1           % The lowest amount the m-value can change each step
    CW_input.Input1.mStepLength=0.5      % The length of mirrors where m-value must be kept constant
    CW_input.Input1.remove_m_under=1;    % Hard limit to lowest m-value
    CW_input.Input1.remove_m_over=6;     % Hard limit to highest m-value
    CW_input.Input1.minNpoly=3           % N is the power 
    CW_input.Input1.maxNpoly=3
    CW_input.Input1.substrate='BK'
    
        % Segment 1 - furthest from source, closest to sample.
    CW_input.Input2.type='C'             % Guide_bot segmet type 
    CW_input.Input2.mode='specific'      % Determines input variables. "specific" is advised in this version of CoatingWriter
    CW_input.Input2.distribution='constant' % The coating distribution. 
    CW_input.Input2.fixSides='curve'        % Optimize sides together. Options are: 'HV', 'curve', 'all', 'none'. ('HV' optimizes horizontal and vertical independantly)
    CW_input.Input2.deltaM=0.1           % The lowest amount the m-value can change each step
    CW_input.Input2.mStepLength=0.5      % The length of mirrors where m-value must be kept constant
    CW_input.Input2.remove_m_under=1;    % Hard limit to lowest m-value
    CW_input.Input2.remove_m_over=6;     % Hard limit to highest m-value
    CW_input.Input2.substrate='BK'
    

    % Segment 2 
    CW_input.Input3.type='G'
    CW_input.Input3.mode='specific'
    CW_input.Input3.distribution='none'

    % Segment 3
    CW_input.Input4.type='P'
    CW_input.Input4.mode='specific'
    CW_input.Input4.distribution='power'
    CW_input.Input4.mincenter=0
    CW_input.Input4.maxcenter=1
    CW_input.Input4.fixSides='HV';
    CW_input.Input4.deltaM=0.1
    CW_input.Input4.mStepLength=0.5
    CW_input.Input4.remove_m_under=1;
    CW_input.Input4.remove_m_over=6;
    CW_input.Input4.substrate='RH'



%Hvis linjen herunder udkommenteres, printes CW ikke automatisk til
%instrumentfiler.
CW_input.filePath='/media/martin/HDD1/Dropbox/Dropbox/Dropbox/guideBot_newest/test1/PGCE'

[ifit,Declare,Define,Initialize,Coating,CoatingFirst]=CoatingWriter(CW_input);
%clc
%% Test-print af output
if 0
    for i=1:length(ifit)
        fprintf('%s\n',ifit{i});
    end
    for i=1:length(Declare)
        fprintf('%s\n',Declare{i});
    end
    for i=1:length(Define)
        fprintf('%s\n',Define{i});
    end
    for i=1:length(Initialize)
        fprintf('%s\n',Initialize{i});
    end
    for i=1:length(CoatingFirst)
        fprintf('%s\n',CoatingFirst{i});
    end
    for i=1:length(Coating)
        fprintf('%s\n',Coating{i});
    end
end

%alpharight=alpha1horizontal,alphaleft=alpha1horizontal,alphatop=alpha1vertical,alphabottom=alpha1vertical, mright=mValues1horizontal,mleft=mValues1horizontal,mtop=mValues1vertical,mbottom=mValues1vertical

