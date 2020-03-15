function coatingOptions = PriceFunctionValues(coatingOptions)
%% Price Function Values
%%%%%%%%%%%%%%%
% In this function new price estimates can be made.
% To create a new price function simply make a new entry in the if-statements
% If substrate mode is set to "none" the prices of substrate can be included in the areal dependancy
% If the coating prices are known as a function of area and m-vaule, just copy the 'SwissNeutronics2017' entry and make a new .coatingArealDependancy
% This is done in the file coatingPrice.m
%
% Contact Martin Olsen if you need help making a new price function
% martinolsen.mo@gmail.com
%%%%%%%%%%%%%%%


if strcmp(coatingOptions.priceFun,'SwissNeutronics2017')
        coatingOptions.substrateMode='none';              % Either power-function dependant on length or area of the segment, or constant ('none' to ignore all substrate in prices)
        coatingOptions.substrateArealDependancy='linear'; % linear or power dependancy - disregarded if no areal dependancy
        coatingOptions.coatingMode='area';                % Either power-function dependant on length or area of the segment, or constant
        coatingOptions.coatingArealDependancy='SNfit1';   % linear or power dependancy  - disregarded if no areal dependancy


elseif strcmp(coatingOptions.priceFun,'test')
	coatingOptions.substrateMode='none';            
        coatingOptions.substrateArealDependancy='linear';
        coatingOptions.coatingMode='area';                
        coatingOptions.coatingArealDependancy='SNfit1';   

else
    error('Price function not recognized')
end

