function [Lines_instr,Lines_ifit] = coatingPrice(i,n,panel,coatingOptions,Lines_ifit,Lines_instr,numberInThisRun,freePanels)
%COATINGPRICE Summary of this function goes here
%   Detailed explanation goes here

       


if strcmp(coatingOptions.coatingMode,'length')
    if strcmp(coatingOptions.coatingArealDependancy,'power')
        if i==1&&n==1;Lines_ifit{end+1}=sprintf('p.beta=''2.6''; %%Price exponent - Prices use: Price=Length*m^Beta.');else;end
        if i==1&&n==1;Lines_ifit{end+1}=sprintf('p.guidePrice=''0''; %%Base price pr meter of guide without coating.');else;end
        if freePanels > 1
            Lines_instr.COATING{end+1}=sprintf('    sumOfMValues += (elementLength%d[segment-1] * (pow(mValues%i%s[segment-1],beta))/%i)*%i ;',i,i,panel,4,numberInThisRun);  %% Sum of m-value-meter
        else
            Lines_instr.COATING{end+1}=sprintf('    sumOfMValues += (elementLength%d[segment-1] * pow(mValues%i%s[segment-1],beta)) ;',i,i,panel);  %% Sum of m-value-meter
        end
    end
    
    if strcmp(coatingOptions.coatingArealDependancy,'linear')
        if i==1&&n==1;Lines_ifit{end+1}=sprintf('p.beta=''2.6''; %%Price exponent - Prices use: Price=Length*m*beta.');else;end
        if i==1&&n==1;Lines_ifit{end+1}=sprintf('p.guidePrice=''0''; %%Base price pr meter of guide without coating.');else;end
        if freePanels > 1
            Lines_instr.COATING{end+1}=sprintf('    sumOfMValues += (elementLength%d[segment-1] * (mValues%i%s[segment-1]*beta)/%i)*%i ;',i,i,panel,4,numberInThisRun);  %% Sum of m-value-meter
        else
            Lines_instr.COATING{end+1}=sprintf('    sumOfMValues += (elementLength%d[segment-1] * mValues%i%s[segment-1]*beta) ;',i,i,panel);  %% Sum of m-value-meter
        end
    end
    
    
    
elseif strcmp(coatingOptions.coatingMode,'area')
    
    if strcmp(coatingOptions.coatingArealDependancy,'power')
        if i==1&&n==1;Lines_ifit{end+1}=sprintf('p.beta=''2.6''; %%Price exponent - Prices use: Price=Area*m^Beta.');else;end
        if i==1&&n==1;Lines_ifit{end+1}=sprintf('p.guidePrice=''0''; %%Base price pr meter of guide without coating.');else;end
        if freePanels > 1
            if strcmp(panel,'horizontal')
                Lines_instr.COATING{end+1}=sprintf('    sumOfMValues += (2*Segment%iHorizontalArea[segment-1] * (pow(mValues%i%s[segment-1],beta))) ;',i,i,panel);  %% Sum of m-value-meter
            elseif strcmp(panel,'vertical')
                Lines_instr.COATING{end+1}=sprintf('    sumOfMValues += (2*Segment%iVerticalArea[segment-1] * (pow(mValues%i%s[segment-1],beta))) ;',i,i,panel);  %% Sum of m-value-meter
            elseif strcmp(panel,'top') || strcmp(panel,'bottom')
                Lines_instr.COATING{end+1}=sprintf('    sumOfMValues += (Segment%iVerticalArea[segment-1] * (pow(mValues%i%s[segment-1],beta))) ;',i,i,panel);  %% Sum of m-value-meter
            elseif strcmp(panel,'rigt') || strcmp(panel,'left')
                Lines_instr.COATING{end+1}=sprintf('    sumOfMValues += (Segment%iHorizontalArea[segment-1] * (pow(mValues%i%s[segment-1],beta))) ;',i,i,panel);  %% Sum of m-value-meter
            end
        else
            Lines_instr.COATING{end+1}=sprintf('    sumOfMValues += ((2*Segment%iHorizontalArea[segment-1]+2*Segment%iVerticalArea[segment-1]) * pow(mValues%i%s[segment-1],beta)) ;',i,i,panel);  %% Sum of m-value-meter
        end
    end
    
    if strcmp(coatingOptions.coatingArealDependancy,'linear')
        if i==1&&n==1;Lines_ifit{end+1}=sprintf('p.beta=''2.6''; %%Price factor - Prices use: Price=Area*m*Beta.');else;end
        if i==1&&n==1;Lines_ifit{end+1}=sprintf('p.guidePrice=''0''; %%Base price pr meter of guide without coating.');else;end
        if freePanels > 1
            if strcmp(panel,'horizontal')
                Lines_instr.COATING{end+1}=sprintf('    sumOfMValues += (2*Segment%iHorizontalArea[segment-1] * (mValues%i%s[segment-1]*beta)) ;',i,i,panel);  %% Sum of m-value-meter
            elseif strcmp(panel,'vertical')
                Lines_instr.COATING{end+1}=sprintf('    sumOfMValues += (2*Segment%iVerticalArea[segment-1] * (mValues%i%s[segment-1]*beta)) ;',i,i,panel);  %% Sum of m-value-meter
            elseif strcmp(panel,'top') || strcmp(panel,'bottom')
                Lines_instr.COATING{end+1}=sprintf('    sumOfMValues += (Segment%iVerticalArea[segment-1] * (mValues%i%s[segment-1]*beta)) ;',i,i,panel);  %% Sum of m-value-meter
            elseif strcmp(panel,'rigt') || strcmp(panel,'left')
                Lines_instr.COATING{end+1}=sprintf('    sumOfMValues += (Segment%iHorizontalArea[segment-1] * (mValues%i%s[segment-1]*beta)) ;',i,i,panel);  %% Sum of m-value-meter
            end
        else
            Lines_instr.COATING{end+1}=sprintf('    sumOfMValues += ((2*Segment%iHorizontalArea[segment-1]+2*Segment%iVerticalArea[segment-1]) * mValues%i%s[segment-1]*beta) ;',i,i,panel);  %% Sum of m-value-meter
        end
    end
    
    
    if strcmp(coatingOptions.coatingArealDependancy,'SNfit1')
        %% This fit is a total coating + substrate price found from numbers from a SwissNeutronics guide estimate. Numbers are ESTIMATES!
        % The fit is on the form (P1 * A) * m-val ^2.5 + (P2 * A + P3 * L)
        % Note that P1 and P3 are constant for all substrate types
        % Each calculation is multiplied by 4 to account for four sides of
        % a neutron guide. If only one side is to be calculated it is
        % multiplied by 0.25.
        
        if n==1;
            Lines_ifit{end+1}=sprintf('p.Substrate%i=''%s''; %% | metal: RH | BorKron: BK | BorFloat: BF | Sodium glass: Na |',i,coatingOptions.substrate);
            Lines_instr.DEFINE{end+1}=sprintf('string Substrate%i="%s",',i,coatingOptions.substrate)
        end
        if strcmp(coatingOptions.segmentType,'S') || strcmp(coatingOptions.segmentType,'C')
            if strcmp(panel,'vertical')
                Lines_instr.COATING{end+1}=sprintf('\tif (strcmp(Substrate%i,"RH") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 * ((2.183*4*Segment%iHorizontalArea)*pow(mValues%i%s,2.5)+(16.520*4*Segment%iHorizontalArea+3.992*length%i))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BK") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 * ((2.183*4*Segment%iHorizontalArea)*pow(mValues%i%s,2.5)+(11.66*4*Segment%iHorizontalArea+3.992*length%i))   ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BF") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 *  ((2.183*4*Segment%iHorizontalArea)*pow(mValues%i%s,2.5)+(2.771*4*Segment%iHorizontalArea+3.992*length%i))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"Na") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 *  ((2.183*4*Segment%iHorizontalArea)*pow(mValues%i%s,2.5)+(1.889*4*Segment%iHorizontalArea+3.992*length%i))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
            elseif strcmp(panel,'horizontal')
                Lines_instr.COATING{end+1}=sprintf('\tif (strcmp(Substrate%i,"RH") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 * ((2.183*4*Segment%iVerticalArea)*pow(mValues%i%s,2.5)+(16.520*4*Segment%iVerticalArea+3.992*length%i)) ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BK") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 * ((2.183*4*Segment%iVerticalArea)*pow(mValues%i%s,2.5)+(11.66*4*Segment%iVerticalArea+3.992*length%i))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BF") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 *  ((2.183*4*Segment%iVerticalArea)*pow(mValues%i%s,2.5)+(2.771*4*Segment%iVerticalArea+3.992*length%i))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"Na") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 *  ((2.183*4*Segment%iVerticalArea)*pow(mValues%i%s,2.5)+(1.889*4*Segment%iVerticalArea+3.992*length%i))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
            elseif strcmp(panel,'left') || strcmp(panel,'right')
                Lines_instr.COATING{end+1}=sprintf('\tif (strcmp(Substrate%i,"RH") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 * ((2.183*4*Segment%iVerticalArea)*pow(mValues%i%s,2.5)+(16.520*4*Segment%iVerticalArea+0.25*3.992*length%i)) ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BK") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 * ((2.183*4*Segment%iVerticalArea)*pow(mValues%i%s,2.5)+(11.66*4*Segment%iVerticalArea+0.25*3.992*length%i))   ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BF") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 *  ((2.183*4*Segment%iVerticalArea)*pow(mValues%i%s,2.5)+(2.771*4*Segment%iVerticalArea+0.25*3.992*length%i))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"Na") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 *  ((2.183*4*Segment%iVerticalArea)*pow(mValues%i%s,2.5)+(1.889*4*Segment%iVerticalArea+0.25*3.992*length%i))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
            elseif strcmp(panel,'top') || strcmp(panel,'bottom')
                Lines_instr.COATING{end+1}=sprintf('\tif (strcmp(Substrate%i,"RH") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 *  ((2.183*4*Segment%iHorizontalArea)*pow(mValues%i%s,2.5)+(16.520*4*Segment%iHorizontalArea+0.25*3.992*length%i))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BK") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 * ((2.183*4*Segment%iHorizontalArea)*pow(mValues%i%s,2.5)+(11.66*4*Segment%iHorizontalArea+0.25*3.992*length%i))   ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BF") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 *  ((2.183*4*Segment%iHorizontalArea)*pow(mValues%i%s,2.5)+(2.771*4*Segment%iHorizontalArea+0.25*3.992*length%i))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"Na") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 *  ((2.183*4*Segment%iHorizontalArea)*pow(mValues%i%s,2.5)+(1.889*4*Segment%iHorizontalArea+0.25*3.992*length%i))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
           elseif strcmp(panel,'inside') || strcmp(panel,'outside')
                Lines_instr.COATING{end+1}=sprintf('\tif (strcmp(Substrate%i,"RH") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 *  ((2.183*4*Segment%iVerticalArea)*pow(mValues%i%s,2.5)+(16.520*4*Segment%iVerticalArea+0.25*3.992*length%i))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BK") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 * ((2.183*4*Segment%iVerticalArea)*pow(mValues%i%s,2.5)+(11.66*4*Segment%iVerticalArea+0.25*3.992*length%i))   ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BF") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 *  ((2.183*4*Segment%iVerticalArea)*pow(mValues%i%s,2.5)+(2.771*4*Segment%iVerticalArea+0.25*3.992*length%i))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"Na") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 *  ((2.183*4*Segment%iVerticalArea)*pow(mValues%i%s,2.5)+(1.889*4*Segment%iVerticalArea+0.25*3.992*length%i))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
            else
                Lines_instr.COATING{end+1}=sprintf('\tif (strcmp(Substrate%i,"RH") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 *  ((2.183*4*Segment%iHorizontalArea)*pow(mValues%i%s,2.5)+(16.520*4*Segment%iHorizontalArea+3.992*length%i))  +4/2 *  ((2.183*4*Segment%iVerticalArea)*pow(mValues%i%s,2.5)+(16.520*Segment%iVerticalArea+3.992*length%i))  ;',i,i,panel,i,i,i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BK") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 *  ((2.183*4*Segment%iHorizontalArea)*pow(mValues%i%s,2.5)+(11.66*4*Segment%iHorizontalArea+3.992*length%i))  +4/2 *  ((2.183*4*Segment%iVerticalArea)*pow(mValues%i%s,2.5)+(11.66*Segment%iVerticalArea+3.992*length%i))  ;',i,i,panel,i,i,i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BF") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 *  ((2.183*4*Segment%iHorizontalArea)*pow(mValues%i%s,2.5)+(2.771*4*Segment%iHorizontalArea+3.992*length%i))  +4/2 *  ((2.183*4*Segment%iVerticalArea)*pow(mValues%i%s,2.5)+(2.771*Segment%iVerticalArea+3.992*length%i))  ;',i,i,panel,i,i,i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"Na") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 *  ((2.183*4*Segment%iHorizontalArea)*pow(mValues%i%s,2.5)+(1.889*4*Segment%iHorizontalArea+3.992*length%i))  +4/2 *  ((2.183*4*Segment%iVerticalArea)*pow(mValues%i%s,2.5)+(1.887*Segment%iVerticalArea+3.992*length%i))  ;',i,i,panel,i,i,i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
            end
        else
            if strcmp(panel,'vertical')
                Lines_instr.COATING{end+1}=sprintf('\tif (strcmp(Substrate%i,"RH") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 * ((2.183*4*Segment%iHorizontalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(16.520*4*Segment%iHorizontalArea[segment-1]+3.992*elementLength%i[segment-1]))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BK") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 * ((2.183*4*Segment%iHorizontalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(11.66*4*Segment%iHorizontalArea[segment-1]+3.992*elementLength%i[segment-1]))   ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BF") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 *  ((2.183*4*Segment%iHorizontalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(2.771*4*Segment%iHorizontalArea[segment-1]+3.992*elementLength%i[segment-1]))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"Na") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 *  ((2.183*4*Segment%iHorizontalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(1.889*4*Segment%iHorizontalArea[segment-1]+3.992*elementLength%i[segment-1]))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
            elseif strcmp(panel,'horizontal')
                Lines_instr.COATING{end+1}=sprintf('\tif (strcmp(Substrate%i,"RH") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 * ((2.183*4*Segment%iVerticalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(16.520*4*Segment%iVerticalArea[segment-1]+3.992*elementLength%i[segment-1])) ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BK") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 * ((2.183*4*Segment%iVerticalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(11.66*4*Segment%iVerticalArea[segment-1]+3.992*elementLength%i[segment-1]))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BF") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 *  ((2.183*4*Segment%iVerticalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(2.771*4*Segment%iVerticalArea[segment-1]+3.992*elementLength%i[segment-1]))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"Na") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 *  ((2.183*4*Segment%iVerticalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(1.889*4*Segment%iVerticalArea[segment-1]+3.992*elementLength%i[segment-1]))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
            elseif strcmp(panel,'left') || strcmp(panel,'right')
                Lines_instr.COATING{end+1}=sprintf('\tif (strcmp(Substrate%i,"RH") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 * ((2.183*4*Segment%iVerticalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(16.520*4*Segment%iVerticalArea[segment-1]+0.25*3.992*elementLength%i[segment-1])) ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BK") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 * ((2.183*4*Segment%iVerticalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(11.66*4*Segment%iVerticalArea[segment-1]+0.25*3.992*elementLength%i[segment-1]))   ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BF") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 *  ((2.183*4*Segment%iVerticalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(2.771*4*Segment%iVerticalArea[segment-1]+0.25*3.992*elementLength%i[segment-1]))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"Na") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 *  ((2.183*4*Segment%iVerticalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(1.889*4*Segment%iVerticalArea[segment-1]+0.25*3.992*elementLength%i[segment-1]))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
            elseif strcmp(panel,'top') || strcmp(panel,'bottom')
                Lines_instr.COATING{end+1}=sprintf('\tif (strcmp(Substrate%i,"RH") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 *  ((2.183*4*Segment%iHorizontalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(16.520*4*Segment%iHorizontalArea[segment-1]+0.25*3.992*elementLength%i[segment-1]))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BK") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 * ((2.183*4*Segment%iHorizontalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(11.66*4*Segment%iHorizontalArea[segment-1]+0.25*3.992*elementLength%i[segment-1]))   ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BF") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 *  ((2.183*4*Segment%iHorizontalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(2.771*4*Segment%iHorizontalArea[segment-1]+0.25*3.992*elementLength%i[segment-1]))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"Na") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.25 *  ((2.183*4*Segment%iHorizontalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(1.889*4*Segment%iHorizontalArea[segment-1]+0.25*3.992*elementLength%i[segment-1]))  ;',i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
            else
                Lines_instr.COATING{end+1}=sprintf('\tif (strcmp(Substrate%i,"RH") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 *  ((2.183*4*Segment%iHorizontalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(16.520*4*Segment%iHorizontalArea[segment-1]+3.992*elementLength%i[segment-1]))  +2 *  ((2.183*4*Segment%iVerticalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(16.520*4*Segment%iVerticalArea[segment-1]+3.992*elementLength%i[segment-1]))  ;',i,i,panel,i,i,i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BK") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 *  ((2.183*4*Segment%iHorizontalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(11.66*4*Segment%iHorizontalArea[segment-1]+3.992*elementLength%i[segment-1]))  +2 *  ((2.183*4*Segment%iVerticalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(11.66*4*Segment%iVerticalArea[segment-1]+3.992*elementLength%i[segment-1]))  ;',i,i,panel,i,i,i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"BF") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 *  ((2.183*4*Segment%iHorizontalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(2.771*4*Segment%iHorizontalArea[segment-1]+3.992*elementLength%i[segment-1]))  +2 *  ((2.183*4*Segment%iVerticalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(2.771*4*Segment%iVerticalArea[segment-1]+3.992*elementLength%i[segment-1]))  ;',i,i,panel,i,i,i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
                Lines_instr.COATING{end+1}=sprintf('\telse if (strcmp(Substrate%i,"Na") == 0)',i);
                Lines_instr.COATING{end+1}=sprintf('\t{');
                Lines_instr.COATING{end+1}=sprintf('\t\tsumOfMValues += 0.5 *  ((2.183*4*Segment%iHorizontalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(1.889*4*Segment%iHorizontalArea[segment-1]+3.992*elementLength%i[segment-1]))  +2 *  ((2.183*4*Segment%iVerticalArea[segment-1])*pow(mValues%i%s[segment-1],2.5)+(1.887*4*Segment%iVerticalArea[segment-1]+3.992*elementLength%i[segment-1]))  ;',i,i,panel,i,i,i,i,panel,i,i);  %% Sum of m-value-meter
                Lines_instr.COATING{end+1}=sprintf('\t}');
            end
        end
     end
end   
end

