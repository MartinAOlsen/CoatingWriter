function [Options]=DefaultCoatingOptions


%% General default options
Options.default_Budget.price=1500;

%% Global defaults
Default.makeFree=0;
Default.printResult=1;
Default.substrate='BK';

Options.mode.substrateMode='none'; 
Options.mode.substrateArealDependancy='linear';
Options.mode.coatingMode='area';  
Options.mode.coatingArealDependancy='SNfit1';
Options.options.budgetMaxError=10;                 % Maximum unpunished error in price
Options.options.punishment='potential';            % 'shutter' cuts neutrons off, 'potential' gives a price-punishment depending of the error.   
Options.options.potentialPower=1.5;   


Options.options.AnalyzeMode='standard';
Options.options.OptimizeMode='standard';
Options.options.SpeedScanMode='off';


%% E-segment
Options.default_E_exp.segmentType='E';
Options.default_E_exp.mode='specific';
Options.default_E_exp.minM=1;
Options.default_E_exp.maxM=6;
Options.default_E_exp.minAvgM=2;
Options.default_E_exp.maxAvgM=3;
Options.default_E_exp.minCenter=0.5;
Options.default_E_exp.maxCenter=0.5;
Options.default_E_exp.minA=1;
Options.default_E_exp.maxA=1;
Options.default_E_exp.minC=1;
Options.default_E_exp.maxC=5;
Options.default_E_exp.minB=1;
Options.default_E_exp.maxB=30;
Options.default_E_exp.fixSides='all';
Options.default_E_exp.remove_m_under=0;  
Options.default_E_exp.remove_m_over=0; 
Options.default_E_exp.nSegments=100;
Options.default_E_exp.distribution='exp';


Options.default_E_power.segmentType='E';
Options.default_E_power.mode='specific';
Options.default_E_power.minM=1;
Options.default_E_power.maxM=6;
Options.default_E_power.minAvgM=2;
Options.default_E_power.maxAvgM=3;
Options.default_E_power.minCenter=0.5;
Options.default_E_power.maxCenter=0.5;
Options.default_E_power.minA=1;
Options.default_E_power.maxA=6;
Options.default_E_power.minC=0;
Options.default_E_power.maxC=0;
Options.default_E_power.minB=0;
Options.default_E_power.maxB=6;
Options.default_E_power.fixSides='all';
Options.default_E_power.remove_m_under=0;  
Options.default_E_power.remove_m_over=0; 
Options.default_E_power.nSegments=100;
Options.default_E_power.minPower=2;
Options.default_E_power.maxPower=6;
Options.default_E_power.distribution='power';


Options.default_E_linear.segmentType='E';
Options.default_E_linear.mode='specific';
Options.default_E_linear.minM=1;
Options.default_E_linear.maxM=6;
Options.default_E_linear.minAvgM=2;
Options.default_E_linear.maxAvgM=3;
Options.default_E_linear.minCenter=0.5;
Options.default_E_linear.maxCenter=0.5;
Options.default_E_linear.minA=1;
Options.default_E_linear.maxA=6;
Options.default_E_linear.minC=1;
Options.default_E_linear.maxC=5;
Options.default_E_linear.minB=1;
Options.default_E_linear.maxB=30;
Options.default_E_linear.fixSides='all';
Options.default_E_linear.remove_m_under=0;  
Options.default_E_linear.remove_m_over=0; 
Options.default_E_linear.nSegments=100;
Options.default_E_linear.distribution='linear';

Options.default_E_constant.segmentType='E';
Options.default_E_constant.mode='specific';
Options.default_E_constant.minM=1;
Options.default_E_constant.maxM=6;
Options.default_E_constant.minAvgM=1;
Options.default_E_constant.maxAvgM=5;
Options.default_E_constant.minCenter=0.5;
Options.default_E_constant.maxCenter=0.5;
Options.default_E_constant.minA=0;
Options.default_E_constant.maxA=0;
Options.default_E_constant.minC=0;
Options.default_E_constant.maxC=0;
Options.default_E_constant.minB=0;
Options.default_E_constant.maxB=0;
Options.default_E_constant.fixSides='all';
Options.default_E_constant.remove_m_under=0;  
Options.default_E_constant.remove_m_over=0; 
Options.default_E_constant.nSegments=1;
Options.default_E_constant.distribution='constant';



%% P segment
Options.default_P_exp.segmentType='P';
Options.default_P_exp.mode='specific';
Options.default_P_exp.minM=1;
Options.default_P_exp.maxM=6;
Options.default_P_exp.minAvgM=2;
Options.default_P_exp.maxAvgM=3;
Options.default_P_exp.minCenter=0.5;
Options.default_P_exp.maxCenter=0.5;
Options.default_P_exp.minA=1;
Options.default_P_exp.maxA=1;
Options.default_P_exp.minC=1;
Options.default_P_exp.maxC=5;
Options.default_P_exp.minB=1;
Options.default_P_exp.maxB=30;
Options.default_P_exp.fixSides='all';
Options.default_P_exp.remove_m_under=0;  
Options.default_P_exp.remove_m_over=0; 
Options.default_P_exp.nSegments=100;
Options.default_P_exp.distribution='exp';

Options.default_P_power.segmentType='P';
Options.default_P_power.mode='specific';
Options.default_P_power.minM=1;
Options.default_P_power.maxM=6;
Options.default_P_power.minAvgM=2;
Options.default_P_power.maxAvgM=3;
Options.default_P_power.minCenter=0.5;
Options.default_P_power.maxCenter=0.5;
Options.default_P_power.minA=1;
Options.default_P_power.maxA=6;
Options.default_P_power.minC=0;
Options.default_P_power.maxC=0;
Options.default_P_power.minB=0;
Options.default_P_power.maxB=6;
Options.default_P_power.minPower=2;
Options.default_P_power.maxPower=2;
Options.default_P_power.fixSides='all';
Options.default_P_power.remove_m_under=0;  
Options.default_P_power.remove_m_over=0; 
Options.default_P_power.nSegments=100;
Options.default_P_power.distribution='power';

Options.default_P_linear.segmentType='P';
Options.default_P_linear.mode='specific';
Options.default_P_linear.minM=1;
Options.default_P_linear.maxM=6;
Options.default_P_linear.minAvgM=2;
Options.default_P_linear.maxAvgM=3;
Options.default_P_linear.minCenter=0.5;
Options.default_P_linear.maxCenter=0.5;
Options.default_P_linear.minA=1;
Options.default_P_linear.maxA=6;
Options.default_P_linear.minC=1;
Options.default_P_linear.maxC=5;
Options.default_P_linear.minB=1;
Options.default_P_linear.maxB=30;
Options.default_P_linear.fixSides='all';
Options.default_P_linear.remove_m_under=0;  
Options.default_P_linear.remove_m_over=0; 
Options.default_P_linear.nSegments=100;
Options.default_P_linear.distribution='linear';

%% S segment
Options.default_S_exp.segmentType='S';
Options.default_S_exp.mode='specific';
Options.default_S_exp.minM=1;
Options.default_S_exp.maxM=6;
Options.default_S_exp.minAvgM=1;
Options.default_S_exp.maxAvgM=5;
Options.default_S_exp.minCenter=0.5;
Options.default_S_exp.maxCenter=0.5;
Options.default_S_exp.minA=1;
Options.default_S_exp.maxA=5;
Options.default_S_exp.minC=0;
Options.default_S_exp.maxC=0;
Options.default_S_exp.minB=0;
Options.default_S_exp.maxB=0;
Options.default_S_exp.fixSides='all';
Options.default_S_exp.remove_m_under=0;  
Options.default_S_exp.remove_m_over=0; 
Options.default_S_exp.nSegments=1;
Options.default_S_exp.distribution='exp';

Options.default_S_power.segmentType='S';
Options.default_S_power.mode='specific';
Options.default_S_power.minM=1;
Options.default_S_power.maxM=6;
Options.default_S_power.minAvgM=1;
Options.default_S_power.maxAvgM=5;
Options.default_S_power.minCenter=0.5;
Options.default_S_power.maxCenter=0.5;
Options.default_S_power.minA=1;
Options.default_S_power.maxA=5;
Options.default_S_power.minC=0;
Options.default_S_power.maxC=0;
Options.default_S_power.minB=0;
Options.default_S_power.maxB=0;
Options.default_S_power.minPower=1;
Options.default_S_power.maxPower=1;
Options.default_S_power.fixSides='all';
Options.default_S_power.remove_m_under=0;  
Options.default_S_power.remove_m_over=0; 
Options.default_S_power.nSegments=1;
Options.default_S_power.distribution='power';

Options.default_S_linear.segmentType='S';
Options.default_S_linear.mode='specific';
Options.default_S_linear.minM=1;
Options.default_S_linear.maxM=6;
Options.default_S_linear.minAvgM=1;
Options.default_S_linear.maxAvgM=5;
Options.default_S_linear.minCenter=0.5;
Options.default_S_linear.maxCenter=0.5;
Options.default_S_linear.minA=1;
Options.default_S_linear.maxA=5;
Options.default_S_linear.minC=0;
Options.default_S_linear.maxC=0;
Options.default_S_linear.minB=0;
Options.default_S_linear.maxB=0;
Options.default_S_linear.fixSides='all';
Options.default_S_linear.remove_m_under=0;  
Options.default_S_linear.remove_m_over=0; 
Options.default_S_linear.nSegments=1;
Options.default_S_linear.distribution='linear';

Options.default_S_constant.segmentType='S';
Options.default_S_constant.mode='specific';
Options.default_S_constant.minM=1;
Options.default_S_constant.maxM=6;
Options.default_S_constant.minAvgM=1;
Options.default_S_constant.maxAvgM=5;
Options.default_S_constant.minCenter=0.5;
Options.default_S_constant.maxCenter=0.5;
Options.default_S_constant.minA=0;
Options.default_S_constant.maxA=0;
Options.default_S_constant.minC=0;
Options.default_S_constant.maxC=0;
Options.default_S_constant.minB=0;
Options.default_S_constant.maxB=0;
Options.default_S_constant.fixSides='all';
Options.default_S_constant.remove_m_under=0;  
Options.default_S_constant.remove_m_over=0; 
Options.default_S_constant.nSegments=1;
Options.default_S_constant.distribution='constant';

%% G segment
Options.default_G_none.segmentType='G';
Options.default_G_none.price='0';
Options.default_G_none.distribution='none';
Options.default_G_none.fixSides='all';

%%
Options.default_C_constant.segmentType='C';
Options.default_C_constant.mode='specific';
Options.default_C_constant.minM=1;
Options.default_C_constant.maxM=6;
Options.default_C_constant.minAvgM=1;
Options.default_C_constant.maxAvgM=5;
Options.default_C_constant.minCenter=0.5;
Options.default_C_constant.maxCenter=0.5;
Options.default_C_constant.minA=0;
Options.default_C_constant.maxA=0;
Options.default_C_constant.minC=0;
Options.default_C_constant.maxC=0;
Options.default_C_constant.minB=0;
Options.default_C_constant.maxB=0;
Options.default_C_constant.fixSides='all';
Options.default_C_constant.remove_m_under=0;  
Options.default_C_constant.remove_m_over=0; 
Options.default_C_constant.nSegments=1;
Options.default_C_constant.distribution='constant';


optionNames=fieldnames(Options);
defaultNames=(fieldnames(Default));
for i=1:numel(fieldnames(Options))
    for j=1:numel(fieldnames(Default))
        eval(sprintf('Options.%s.%s=Default.%s',char(optionNames(i)),char(defaultNames(j)),char(defaultNames(j))));
    end
end




end
