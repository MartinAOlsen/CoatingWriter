function [Lines_instr] = Area(i,panel,coatingOptions,Lines_instr)
%AREAL Summary of this function goes here
%   Detailed explanation goes here

% Area of a segment should be of single sides. Meaning that
% true size = 2* Segment%iHorizontalArea + 2* Segment%iVerticalArea 

numSecPrElement=500;
if isfield(coatingOptions,'mStepLength') && isfield(coatingOptions,'instrumentLength')
    numSecPrElement=round(coatingOptions.instrumentLength/coatingOptions.mStepLength)+1;
elseif isfield(coatingOptions,'numSecPrElement')
    numSecPrElement=coatingOptions.numSecPrElement;
end



if i==1; Lines_instr.COATINGfirst{end+1}=sprintf('TotalSubstratePrice=0;');end


if strcmp(coatingOptions.segmentType,'G') == 0 % Don't make geometry approximation of gaps
    Lines_instr.COATINGfirst{end+1}=sprintf('// Create Segment price formula for segment %i %s',i,panel);
    if strcmp(coatingOptions.segmentType,'E') %% Not just an approximation of elliptic guide, but actually calculating the surface of the ellipses
        Lines_instr.DECLARE{end+1}=sprintf('double smallaxis_x%i;',i);
        Lines_instr.DECLARE{end+1}=sprintf('double smallaxis_y%i;',i);
        % majorAxis=sqrt(((length+linin+linout)/2)^2+(smallAxis^2))
        Lines_instr.DECLARE{end+1}=sprintf('double majoraxis_x%i;',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('majoraxis_x%i=sqrt(pow(((length%i+Linx%i+Loutx%i)/2),2)+pow(smallaxis_x%i,2));',i,i,i,i,i);
        Lines_instr.DECLARE{end+1}=sprintf('double majoraxis_y%i;',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('majoraxis_y%i=sqrt(pow(((length%i+Liny%i+Louty%i)/2),2)+pow(smallaxis_y%i,2));',i,i,i,i,i);
        % Length from elliptic start:
        Lines_instr.DECLARE{end+1}=sprintf('double lengthFromStartx%i;',i);
        Lines_instr.DECLARE{end+1}=sprintf('double lengthFromStarty%i;',i);

        Lines_instr.COATINGfirst{end+1}=sprintf('lengthFromStarty%i=(majoraxis_y%i-sqrt(pow(smallaxis_y%i,2)+pow(majoraxis_y%i,2)))+Liny%i;',i,i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('lengthFromStartx%i=(majoraxis_x%i-sqrt(pow(smallaxis_x%i,2)+pow(majoraxis_x%i,2)))+Linx%i;',i,i,i,i,i);

        Lines_instr.DECLARE{end+1}=sprintf('double xpar_x%i;',i);
        Lines_instr.DECLARE{end+1}=sprintf('double xpar_y%i;',i);
%         Lines_instr.DECLARE{end+1}=sprintf('double horizontalSegmentSize [500];');
%         Lines_instr.DECLARE{end+1}=sprintf('double verticalSegmentSize [500];');
        Lines_instr.DECLARE{end+1}=sprintf('double Segment%iHorizontalArea [%i];',i,numSecPrElement);
        Lines_instr.DECLARE{end+1}=sprintf('double Segment%iVerticalArea [%i];',i,numSecPrElement);
        Lines_instr.DECLARE{end+1}=sprintf('double Segment%iSubstrateSize;',i);

	Lines_instr.DECLARE{end+1}=sprintf('double Segment%i_h1[%i];',i,numSecPrElement);
	Lines_instr.DECLARE{end+1}=sprintf('double Segment%i_h2[%i];',i,numSecPrElement);
	Lines_instr.DECLARE{end+1}=sprintf('double Segment%i_w1[%i];',i,numSecPrElement);
	Lines_instr.DECLARE{end+1}=sprintf('double Segment%i_w2[%i];',i,numSecPrElement);

        Lines_instr.COATINGfirst{end+1}=sprintf('double dSegments%i=nSegments%i;',i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('Segment%iSubstrateSize=0;',i);

        Lines_instr.COATINGfirst{end+1}=sprintf('counter2=0;');
        Lines_instr.COATINGfirst{end+1}=sprintf('for (segment =1 ; segment < nSegments%i+1 ; segment++)',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('{');
        Lines_instr.COATINGfirst{end+1}=sprintf('counter2+=1;');
        Lines_instr.COATINGfirst{end+1}=sprintf('\t xpar_x%i=-majoraxis_x%i+lengthFromStartx%i+(counter2-0.5)*(length%i/dSegments%i);',i,i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('\t xpar_y%i=-majoraxis_y%i+lengthFromStarty%i+(counter2-0.5)*(length%i/dSegments%i);',i,i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('\t horizontalSegmentSize[segment-1]=2*(length%i/dSegments%i)*((sqrt(fabs(majoraxis_x%i*majoraxis_x%i-xpar_x%i*xpar_x%i))*smallaxis_x%i)/majoraxis_x%i);',i,i,i,i,i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('\t verticalSegmentSize[segment-1]=2*(length%i/dSegments%i)*((sqrt(fabs(majoraxis_y%i*majoraxis_y%i-xpar_y%i*xpar_y%i))*smallaxis_y%i)/majoraxis_y%i);',i,i,i,i,i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%iHorizontalArea[segment-1]=2*(length%i/dSegments%i)*((sqrt(fabs(majoraxis_x%i*majoraxis_x%i-xpar_x%i*xpar_x%i))*smallaxis_x%i)/majoraxis_x%i);',i,i,i,i,i,i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%iVerticalArea[segment-1]=2*(length%i/dSegments%i)*((sqrt(fabs(majoraxis_y%i*majoraxis_y%i-xpar_y%i*xpar_y%i))*smallaxis_y%i)/majoraxis_y%i);',i,i,i,i,i,i,i,i,i);
% Write heights and widths:
	Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%i_h1[segment-1]=2*((sqrt(fabs(majoraxis_y%i*majoraxis_y%i-(xpar_y%i-(length%i/dSegments%i)*0.5)*(xpar_y%i-(length%i/dSegments%i)*0.5)))*smallaxis_y%i)/majoraxis_y%i);',i,i,i,i,i,i,i,i,i,i,i);
	Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%i_h2[segment-1]=2*((sqrt(fabs(majoraxis_y%i*majoraxis_y%i-(xpar_y%i+(length%i/dSegments%i)*0.5)*(xpar_y%i+(length%i/dSegments%i)*0.5)))*smallaxis_y%i)/majoraxis_y%i);',i,i,i,i,i,i,i,i,i,i,i);
	Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%i_w1[segment-1]=2*((sqrt(fabs(majoraxis_x%i*majoraxis_x%i-(xpar_x%i-(length%i/dSegments%i)*0.5)*(xpar_x%i-(length%i/dSegments%i)*0.5)))*smallaxis_x%i)/majoraxis_x%i);',i,i,i,i,i,i,i,i,i,i,i);
	Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%i_w2[segment-1]=2*((sqrt(fabs(majoraxis_x%i*majoraxis_x%i-(xpar_x%i+(length%i/dSegments%i)*0.5)*(xpar_x%i+(length%i/dSegments%i)*0.5)))*smallaxis_x%i)/majoraxis_x%i);',i,i,i,i,i,i,i,i,i,i,i);




        Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%iSubstrateSize+=2*verticalSegmentSize[segment-1]+2*horizontalSegmentSize[segment-1];',i);

        Lines_instr.COATINGfirst{end+1}=sprintf('}');
        %Debug prints:
        Lines_instr.COATINGfirst{end+1}=sprintf('printf("xpar_x=%%2.2f\\n",xpar_x%i);',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('printf("dSegments=%%2.2f\\n",dSegments%i);',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('printf("length=%%2.2f\\n",length%i);',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('printf("majoraxis_x=%%2.2f\\n",majoraxis_x%i);',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('printf("smallaxis_x=%%2.2f\\n",smallaxis_x%i);',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('printf("sqrt=%%2.2f\\n",(sqrt(majoraxis_x%i*majoraxis_x%i-xpar_x%i*xpar_x%i)*smallaxis_x%i));',i,i,i,i,i);


%         Lines_instr.COATINGfirst{end+1}=sprintf('guideSubstratePrice+=Segment%iSubstrateSize;',i);  
    elseif strcmp(coatingOptions.segmentType,'P') 
        Lines_instr.DECLARE{end+1}=sprintf('double smallaxis_parabolic_x%i;',i);
        Lines_instr.DECLARE{end+1}=sprintf('double smallaxis_parabolic_y%i;',i);
        
        Lines_instr.DECLARE{end+1}=sprintf('double majoraxis_x%i;',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('majoraxis_x%i=sqrt(pow(((length%i+Linx%i+Loutx%i)/2),2)+pow(smallaxis_parabolic_x%i,2));',i,i,i,i,i);
        Lines_instr.DECLARE{end+1}=sprintf('double majoraxis_y%i;',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('majoraxis_y%i=sqrt(pow(((length%i+Liny%i+Louty%i)/2),2)+pow(smallaxis_parabolic_y%i,2));',i,i,i,i,i);
        % Length from elliptic start:
        Lines_instr.DECLARE{end+1}=sprintf('double lengthFromStartx%i;',i);
        Lines_instr.DECLARE{end+1}=sprintf('double lengthFromStarty%i;',i);

        Lines_instr.COATINGfirst{end+1}=sprintf('lengthFromStarty%i=(majoraxis_y%i-sqrt(pow(smallaxis_parabolic_y%i,2)+pow(majoraxis_y%i,2)))+Liny%i;',i,i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('lengthFromStartx%i=(majoraxis_x%i-sqrt(pow(smallaxis_parabolic_x%i,2)+pow(majoraxis_x%i,2)))+Linx%i;',i,i,i,i,i);

        Lines_instr.DECLARE{end+1}=sprintf('double xpar_x%i;',i);
        Lines_instr.DECLARE{end+1}=sprintf('double xpar_y%i;',i);
        Lines_instr.DECLARE{end+1}=sprintf('double Segment%iVerticalArea [%i];',i,numSecPrElement);
        Lines_instr.DECLARE{end+1}=sprintf('double Segment%iVerticalArea [%i];',i,numSecPrElement);
        Lines_instr.DECLARE{end+1}=sprintf('double Segment%iSubstrateSize;',i);
        Lines_instr.DECLARE{end+1}=sprintf('double verticalSubstratePrice;');
        Lines_instr.DECLARE{end+1}=sprintf('double horizontalSubstratePrice;');
        Lines_instr.DECLARE{end+1}=sprintf('double Segment%iHorizontalArea [%i];',i,numSecPrElement);
        Lines_instr.DECLARE{end+1}=sprintf('double Segment%iVerticalArea [%i];',i,numSecPrElement);

	Lines_instr.DECLARE{end+1}=sprintf('double Segment%i_h1[%i];',i,numSecPrElement);
	Lines_instr.DECLARE{end+1}=sprintf('double Segment%i_h2[%i];',i,numSecPrElement);
	Lines_instr.DECLARE{end+1}=sprintf('double Segment%i_w1[%i];',i,numSecPrElement);
	Lines_instr.DECLARE{end+1}=sprintf('double Segment%i_w2[%i];',i,numSecPrElement);

        Lines_instr.COATINGfirst{end+1}=sprintf('double dSegments%i=nSegments%i;',i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('Segment%iSubstrateSize=0;',i);

        
        Lines_instr.COATINGfirst{end+1}=sprintf('for (segment =1 ; segment < nSegments%i+1 ; segment++)',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('{');
        Lines_instr.COATINGfirst{end+1}=sprintf('counter2+=1;');
%         Lines_instr.COATINGfirst{end+1}=sprintf('\t xpar_x%i=-majoraxis_x%i+lengthFromStartx%i+(counter2-0.5)*(length%i/dSegments%i);',i,i,i,i,i);
%         Lines_instr.COATINGfirst{end+1}=sprintf('\t xpar_y%i=-majoraxis_y%i+lengthFromStarty%i+(counter2-0.5)*(length%i/dSegments%i);',i,i,i,i,i);
%         Lines_instr.COATINGfirst{end+1}=sprintf('\t horizontalSegmentSize[segment-1]=2*(length%i/dSegments%i)*((sqrt(fabs(majoraxis_x%i*majoraxis_x%i-xpar_x%i*xpar_x%i))*smallaxis_parabolic_x%i)/majoraxis_x%i);',i,i,i,i,i,i,i,i);
%         Lines_instr.COATINGfirst{end+1}=sprintf('\t verticalSegmentSize[segment-1]=2*(length%i/dSegments%i)*((sqrt(fabs(majoraxis_y%i*majoraxis_y%i-xpar_y%i*xpar_y%i))*smallaxis_parabolic_y%i)/majoraxis_y%i);',i,i,i,i,i,i,i,i);
%         Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%iSubstrateSize+=2*verticalSegmentSize[segment-1]+2*horizontalSegmentSize[segment-1];',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('\t xpar_x%i=-majoraxis_x%i+lengthFromStartx%i+(counter2-0.5)*(length%i/dSegments%i);',i,i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('\t xpar_y%i=-majoraxis_y%i+lengthFromStarty%i+(counter2-0.5)*(length%i/dSegments%i);',i,i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('\t horizontalSegmentSize[segment-1]=2*(length%i/dSegments%i)*((sqrt(fabs(majoraxis_x%i*majoraxis_x%i-xpar_x%i*xpar_x%i))*smallaxis_parabolic_x%i)/majoraxis_x%i);',i,i,i,i,i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('\t verticalSegmentSize[segment-1]=2*(length%i/dSegments%i)*((sqrt(fabs(majoraxis_y%i*majoraxis_y%i-xpar_y%i*xpar_y%i))*smallaxis_parabolic_y%i)/majoraxis_y%i);',i,i,i,i,i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%iHorizontalArea[segment-1]=2*(length%i/dSegments%i)*((sqrt(fabs(majoraxis_x%i*majoraxis_x%i-xpar_x%i*xpar_x%i))*smallaxis_parabolic_x%i)/majoraxis_x%i);',i,i,i,i,i,i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%iVerticalArea[segment-1]=2*(length%i/dSegments%i)*((sqrt(fabs(majoraxis_y%i*majoraxis_y%i-xpar_y%i*xpar_y%i))*smallaxis_parabolic_y%i)/majoraxis_y%i);',i,i,i,i,i,i,i,i,i);

        Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%iHorizontalArea[segment-1]=4*(length%i/dSegments%i)*((sqrt(fabs(majoraxis_x%i*majoraxis_x%i-xpar_x%i*xpar_x%i))*smallaxis_parabolic_x%i)/majoraxis_x%i);',i,i,i,i,i,i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%iVerticalArea[segment-1]=4*(length%i/dSegments%i)*((sqrt(fabs(majoraxis_y%i*majoraxis_y%i-xpar_y%i*xpar_y%i))*smallaxis_parabolic_y%i)/majoraxis_y%i);',i,i,i,i,i,i,i,i,i);

% % Write heights and widths:
% 	Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%i_h1[segment-1]=((sqrt(fabs(majoraxis_y%i*majoraxis_y%i-(xpar_y%i-(length%i/dSegments%i)*0.5)*(xpar_y%i-(length%i/dSegments%i)*0.5)))*smallaxis_parabolic_y%i)/majoraxis_y%i);',i,i,i,i,i,i,i,i,i,i,i);
% 	Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%i_h2[segment-1]=((sqrt(fabs(majoraxis_y%i*majoraxis_y%i-(xpar_y%i+(length%i/dSegments%i)*0.5)*(xpar_y%i+(length%i/dSegments%i)*0.5)))*smallaxis_parabolic_y%i)/majoraxis_y%i);',i,i,i,i,i,i,i,i,i,i,i);
% 	Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%i_w1[segment-1]=((sqrt(fabs(majoraxis_x%i*majoraxis_x%i-(xpar_x%i-(length%i/dSegments%i)*0.5)*(xpar_x%i-(length%i/dSegments%i)*0.5)))*smallaxis_parabolic_x%i)/majoraxis_x%i);',i,i,i,i,i,i,i,i,i,i,i);
% 	Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%i_w2[segment-1]=((sqrt(fabs(majoraxis_x%i*majoraxis_x%i-(xpar_x%i+(length%i/dSegments%i)*0.5)*(xpar_x%i+(length%i/dSegments%i)*0.5)))*smallaxis_parabolic_x%i)/majoraxis_x%i);',i,i,i,i,i,i,i,i,i,i,i);
	
% Write heights and widths:
	Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%i_h1[segment-1]=2*((sqrt(fabs(majoraxis_y%i*majoraxis_y%i-(xpar_y%i-(length%i/dSegments%i)*0.5)*(xpar_y%i-(length%i/dSegments%i)*0.5)))*smallaxis_parabolic_y%i)/majoraxis_y%i);',i,i,i,i,i,i,i,i,i,i,i);
	Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%i_h2[segment-1]=2*((sqrt(fabs(majoraxis_y%i*majoraxis_y%i-(xpar_y%i+(length%i/dSegments%i)*0.5)*(xpar_y%i+(length%i/dSegments%i)*0.5)))*smallaxis_parabolic_y%i)/majoraxis_y%i);',i,i,i,i,i,i,i,i,i,i,i);
	Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%i_w1[segment-1]=2*((sqrt(fabs(majoraxis_x%i*majoraxis_x%i-(xpar_x%i-(length%i/dSegments%i)*0.5)*(xpar_x%i-(length%i/dSegments%i)*0.5)))*smallaxis_parabolic_x%i)/majoraxis_x%i);',i,i,i,i,i,i,i,i,i,i,i);
	Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%i_w2[segment-1]=2*((sqrt(fabs(majoraxis_x%i*majoraxis_x%i-(xpar_x%i+(length%i/dSegments%i)*0.5)*(xpar_x%i+(length%i/dSegments%i)*0.5)))*smallaxis_parabolic_x%i)/majoraxis_x%i);',i,i,i,i,i,i,i,i,i,i,i);


	Lines_instr.COATINGfirst{end+1}=sprintf('printf("w1=%%f, w2=%%1.6f, h1=%%1.6f, h2=%%1.6f\\n",Segment%i_h1[segment-1],Segment%i_h2[segment-1],Segment%i_w1[segment-1],Segment%i_w2[segment-1]);',i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('printf("segmentsize=%%2.2f\\n",horizontalSegmentSize[segment-1]);');

        Lines_instr.COATINGfirst{end+1}=sprintf('}');
        %Debug prints:
        Lines_instr.COATINGfirst{end+1}=sprintf('printf("xpar_x=%%2.2f\\n",xpar_x%i);',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('printf("dSegments=%%2.2f\\n",dSegments%i);',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('printf("length=%%2.2f\\n",length%i);',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('printf("majoraxis_x=%%2.2f\\n",majoraxis_x%i);',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('printf("smallaxis_x=%%2.2f\\n",smallaxis_parabolic_x%i);',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('printf("sqrt=%%2.2f\\n",(sqrt(majoraxis_x%i*majoraxis_x%i-xpar_x%i*xpar_x%i)*smallaxis_parabolic_x%i));',i,i,i,i,i);
	

    elseif strcmp(coatingOptions.segmentType,'S') 
        Lines_instr.COATINGfirst{end+1}=sprintf('Segment%iSubstrateSize=(length%i*(starty%i+endy%i))+(length%i*(startx%i+endx%i));',i,i,i,i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%iHorizontalArea=(length%i*(startx%i+endx%i));',i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%iVerticalArea=(length%i*(starty%i+endy%i));',i,i,i,i);
        Lines_instr.DECLARE{end+1}=sprintf('double Segment%iHorizontalArea;',i);
        Lines_instr.DECLARE{end+1}=sprintf('double Segment%iVerticalArea;',i);

    elseif strcmp(coatingOptions.segmentType,'C') 
       % Use f.eks. curve_small_radius2 and curve_radius2
        
        Lines_instr.INITIALIZE{end+1}=sprintf('// Create Segment price formula for segment %i %s',i,panel);
        Lines_instr.COATINGfirst{end+1}=sprintf('Segment%iSubstrateSize=(length%i*(starty%i+endy%i))+(length%i*(startx%i+endx%i));',i,i,i,i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%iHorizontalArea=(length%i*(startx%i+endx%i));',i,i,i,i);
        Lines_instr.COATINGfirst{end+1}=sprintf('\t Segment%iVerticalArea=(length%i*(starty%i+endy%i));',i,i,i,i);
        Lines_instr.DECLARE{end+1}=sprintf('double Segment%iHorizontalArea;',i);
        Lines_instr.DECLARE{end+1}=sprintf('double Segment%iVerticalArea;',i);

    end

end

