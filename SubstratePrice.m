function [Lines_instr,Lines_ifit] = SubstratePrice(i,n,panel,coatingOptions,Lines_ifit,Lines_instr)
%SUBSTRATEPRICE Summary of this function goes here
%   Detailed explanation goes here

substrateMode=coatingOptions.substrateMode;


if strcmp(substrateMode,'none')
% This is used if the substrate price is not wanted in the price
% calculations
    if i==1
        if n==1
            Lines_instr.COATINGfirst{end+1}=sprintf('// No Substrate Prices')
            Lines_instr.DECLARE{end+1}=sprintf('double TotalSubstratePrice;');
        end
    end
end

if strcmp(substrateMode,'area')
% This mode uses the size of the substrate to approximate a price.     
    %% Approximate geometry for substrate price calculation 
    % Must be before price calculation. Therefore a "COATINGfirst" was
    % needed. This should be addressed in a later build.
    if n==1
        if i==1
             Lines_instr.DECLARE{end+1}=sprintf('double guideSubstratePrice;');
             Lines_instr.COATINGfirst{1}=sprintf('guideSubstratePrice=0;'); 
             Lines_instr.DECLARE{end+1}=sprintf('double counter2;');
             Lines_instr.DECLARE{end+1}=sprintf('double TotalSubstratePrice;');
        end
        if strcmp(coatingOptions.substrateArealDependancy,'linear')
            
            
            Lines_instr.DEFINE{end+1}=sprintf('segment%iSubstrateMaterialPrice=14,',i);
            Lines_instr.DEFINE{end+1}=sprintf('segment%iSubstrateBasePrice=0,',i);
            Lines_ifit{end+1}=sprintf('p.segment%iSubstrateBasePrice=''0''; %% A base price of the substrate',i);
            Lines_ifit{end+1}=sprintf('p.segment%iSubstrateMaterialPrice=''14''; %%Price of the substrate material pr. m^2 in kEuro. Base is 14 for boron and 25 for alumilium',i);
            if i==1
                Lines_instr.COATINGfirst{end+1}=sprintf('TotalSubstratePrice=Segment%iSubstrateSize*segment%iSubstrateMaterialPrice;',i,i);
            else
                Lines_instr.COATINGfirst{end+1}=sprintf('TotalSubstratePrice+=Segment%iSubstrateSize*segment%iSubstrateMaterialPrice;',i,i);
            end
        elseif strcmp(coatingOptions.substrateArealDependancy,'power')
             
            if i==1;Lines_instr.DECLARE{end+1}=sprintf('double TotalSubstratePrice;');end
            Lines_instr.DEFINE{end+1}=sprintf('segment%iSubstrateMaterialPrice=2,',i);
            Lines_instr.DEFINE{end+1}=sprintf('segment%iSubstrateMaterialPricePower=2,',i);
            Lines_instr.DEFINE{end+1}=sprintf('segment%iSubstrateBasePrice=0,',i);
            Lines_ifit{end+1}=sprintf('p.segment%iSubstrateBasePrice=''15''; %%Price of the substrate material pr. m^2 in kEuro. Base is 14 for boron and 25 for alumilium',i);
            Lines_ifit{end+1}=sprintf('p.segment%iSubstrateMaterialPrice=''2''; %%Price of the substrate material. The price is P= BasePrice+MaterialPrice*Areal^MaterialPricePower',i);
            Lines_ifit{end+1}=sprintf('p.segment%iSubstrateMaterialPricePower=''2''; %%',i);
            if i==1
                Lines_instr.COATINGfirst{end+1}=sprintf('TotalSubstratePrice=segment%iSubstrateBasePrice+pow(Segment%iSubstrateSize,segment%iSubstrateMaterialPricePowe)*segment%iSubstrateMaterialPrice;',i,i,i,i);
            else
                Lines_instr.COATINGfirst{end+1}=sprintf('TotalSubstratePrice+=segment%iSubstrateBasePrice+pow(Segment%iSubstrateSize,segment%iSubstrateMaterialPricePowe)*segment%iSubstrateMaterialPrice;',i,i,i,i);
            end
        end
    end
end



end

