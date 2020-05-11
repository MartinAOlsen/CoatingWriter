clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ____            _   _           __        __    _ _             %
%  / ___|___   __ _| |_(_)_ __   __ \ \      / / __(_) |_ ___ _ __  %
% | |   / _ \ / _` | __| | '_ \ / _` \ \ /\ / / '__| | __/ _ \ '__| %
% | |__| (_) | (_| | |_| | | | | (_| |\ V  V /| |  | | ||  __/ |    %
%  \____\___/ \__,_|\__|_|_| |_|\__, | \_/\_/ |_|  |_|\__\___|_|    %
%                               |___/                               %
%                                                                   %  
%  4-step setup for coating- and price optimizations in guide_bot.  %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1 General options
    
    %% CW-mode:
    % coating   - Optimize brilliance transfer with coating distributions
    % value     - Optimize for best performance/price
    % budget    - Optimize for best brilliance transfer with a certain price
    % pricescan - Run several budget runs and examine BT as price increases
    % manual	- Search parameters and make the user choose the best set
    % TEST      - Costum optimizer
    
    Mode='value'
    
    
    %% Price Function:
    % SwissNeutronics2017    - Numbers from Swiss Neutronics early 2017
    
    CW_input.priceFun='SwissNeutronics2017'; 
    
    
    %% Instrument length
    % Determine how much memory is allocated to each section of the guide.
    % Set this to be at least the length of the longest element.
    % Overestimating the length is advised. 200 meter is default.
    CW_input.instrumentLength=200;


%% 2 Mode specific options
    % Depending on the mode chosen in section 1, the options can be set in
    % the following section

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
        CW_input.Price=[2500];                      % Target price
        CW_input.budgetMaxError=25;                 % Maximum unpunished error in price
        CW_input.punishment='potential';            % 'shutter' cuts neutrons off, 'potential' gives a price-punishment depending of the error.   
        CW_input.potentialPower=2;                  % If potential is chosen as punishment, the punishment will be ERROR^potentialPower       
        
    % Price-Scan
    elseif strcmp(Mode,'pricescan')
        CW_input.scanType='pricescan'
        CW_input.Price=[800:100:1800];              % Target price list
        CW_input.budgetMaxError=25;                 % Maximum unpunished error in price
        CW_input.punishment='potential';            % 'shutter' cuts neutrons off, 'potential' gives a price-punishment depending of the error.   
        CW_input.potentialPower=2;                  % If potential is chosen as punishment, the punishment will be ERROR^potentialPower 

    % Manual
    elseif strcmp(Mode,'manual')
        CW_input.scanType='manual'
        CW_input.WavelengthBackgroundMin = 0; 	    % [Å] Minimum wavelength that is counted as background 
        CW_input.WavelengthBackgroundMax = 1; 	    % [Å] Maximum wavelength that is counted as background 

    % Custom scan
    elseif strcmp(Mode,'TEST')
        CW_input.scanType='value'
        CW_input.Price=0;  
        CW_input.criteria=['intensity' , 'price'];  % Choose two criteria from: 'intensity' , 'price' , 'value' , 'background' , 'uniformity'
        CW_input.background_wavelenght=[0,1.25];    % Wavelenghts considered noise
        CW_input.background_divergence = 0;         % [BINARY] toggle if too high divergence is considered noise
        CW_input.background_position = 0;           % [BINARY] toggle if neutrons missing the sample is considered noise

    % Manual
    elseif strcmp(Mode,'manual')
        CW_input.scanType='manual' 
 
    end

    
%% 3 Segment parameters
    % Set parameters for the individual segments in this section
    % Most parameters wil have a min and max value. If a constant value is
    % wanted, set both min and max to the same value
    %
    % Segments are numbered from moderator to sample.
