function [loadStepScan] = StepScanLoader(coatingOptions)
%STEPSCANLOADER makes the lines needed to loat the stepscan files in
%intruments.

Fnames=fieldnames(coatingOptions);

loadStepScan{1}='';
loadStepScan{end+1}=sprintf('// Load stepScan array');
% Make counters
for i = 1:length(Fnames)
   if strfind(Fnames{i},'Input')
         try 
            tmp=Fnames{i};
            thisFixes=eval(['coatingOptions.' Fnames{i} '.fixSides']);
            if str2num(tmp(6:end))>0
                thisNumber=str2num(tmp(6:end));
                loadStepScan{end+1}=sprintf('\t\tint SScounter%i = 0;',thisNumber);
            end
            switch thisFixes
                    case 'none'
                        loadStepScan{end+1}=sprintf('\t\tdouble correctionMValue%i%s[500];',thisNumber,'top');
                        loadStepScan{end+1}=sprintf('\t\tdouble correctionMValue%i%s[500];',thisNumber,'left');
                        loadStepScan{end+1}=sprintf('\t\tdouble correctionMValue%i%s[500];',thisNumber,'bottom');
                        loadStepScan{end+1}=sprintf('\t\tdouble correctionMValue%i%s[500];',thisNumber,'right');
                    case 'all'
                        loadStepScan{end+1}=sprintf('\t\tdouble correctionMValue%i[500];',thisNumber);
                    case 'HV'
                        loadStepScan{end+1}=sprintf('\t\tdouble correctionMValue%i%s[500];',thisNumber,'vertical');
                        loadStepScan{end+1}=sprintf('\t\tdouble correctionMValue%i%s[500];',thisNumber,'horizontal');
                    case 'curve'
                        loadStepScan{end+1}=sprintf('\t\tdouble correctionMValue%i%s[500];',thisNumber,'vertical');
                        loadStepScan{end+1}=sprintf('\t\tdouble correctionMValue%i%s[500];',thisNumber,'inside');
                        loadStepScan{end+1}=sprintf('\t\tdouble correctionMValue%i%s[500];',thisNumber,'outside');
                end
         end
   end
end

loadStepScan{end+1}=sprintf('if (stepScan==1)');
loadStepScan{end+1}=sprintf('{');
loadStepScan{end+1}=sprintf('\tsprintf(summaryName,"CoatingWriter_rawData%%s.txt",scanname);');
loadStepScan{end+1}=sprintf('\tFILE *FscanStep = fopen(summaryName, "r");');
loadStepScan{end+1}=sprintf('\tdouble number;');
loadStepScan{end+1}=sprintf('\tchar line [100];');
loadStepScan{end+1}=sprintf('\twhile(fgets(line, sizeof line, FscanStep) != NULL)');
loadStepScan{end+1}=sprintf('\t{');
loadStepScan{end+1}=sprintf('\t\tchar segmentNr[25]; char segmentType[25]; char segmentLength[25]; char segmentstarty2[25]; char segmentendy2[25]; char segmentstartx2[25]; char segmentendx2[25]; char segmentMtop[25]; char segmentMbottom[25]; char segmentMleft[25]; char segmentMright[25]; char segmentSubstrate[25];');
loadStepScan{end+1}=sprintf('\t\tsscanf(line, "%%[^,],%%[^,],%%[^,],%%[^,],%%[^,],%%[^,],%%[^,],%%[^,],%%[^,],%%[^,],%%[^,],%%[^,]", segmentNr, segmentType, segmentLength, segmentstarty2, segmentendy2, segmentstartx2, segmentendx2, segmentMtop, segmentMbottom, segmentMleft, segmentMright, segmentSubstrate);');
loadStepScan{end+1}=sprintf('\t\tprintf("test= %%s\\n",segmentMtop);');
% Load into correction arrays:

for i = 1:length(Fnames)
   if strfind(Fnames{i},'Input')
        try % Skip if something else has "input" but no relevant options
            thisFixes=eval(['coatingOptions.' Fnames{i} '.fixSides']);
            thisSegmentType=eval(['coatingOptions.' Fnames{i} '.type']);
            tmp=Fnames{i};
            thisNumber=str2num(tmp(6:end));
            loadStepScan{end+1}=sprintf('\t\tif (atof(segmentNr)==%i) // Add to %s -segment coating correction array',thisNumber,thisSegmentType);
            loadStepScan{end+1}=sprintf('\t\t{');
            %%loadStepScan{end+1}=sprintf('\t\tSScounter%i;',thisNumber);
            
            switch thisFixes
                case 'none'
                    loadStepScan{end+1}=sprintf('\t\t\tcorrectionMValue%i%s[SScounter%i]=atof(segmentMtop);',thisNumber,'top',thisNumber);
                    loadStepScan{end+1}=sprintf('\t\t\tcorrectionMValue%i%s[SScounter%i]=atof(segmentMleft);',thisNumber,'left',thisNumber);
                    loadStepScan{end+1}=sprintf('\t\t\tcorrectionMValue%i%s[SScounter%i]=atof(segmentMbottom);',thisNumber,'bottom',thisNumber);
                    loadStepScan{end+1}=sprintf('\t\t\tcorrectionMValue%i%s[SScounter%i]=atof(segmentMright);',thisNumber,'right',thisNumber);
                case 'all'
                    loadStepScan{end+1}=sprintf('\t\t\tcorrectionMValue%i[SScounter%i]=atof(segmentMtop;)',thisNumber,thisNumber);
                case 'HV'
                    loadStepScan{end+1}=sprintf('\t\t\tcorrectionMValue%i%s[SScounter%i]=atof(segmentMtop);',thisNumber,'vertical',thisNumber);
                    loadStepScan{end+1}=sprintf('\t\t\tcorrectionMValue%i%s[SScounter%i]=atof(segmentMleft);',thisNumber,'horizontal',thisNumber);
                case 'curve'
                    loadStepScan{end+1}=sprintf('\t\t\tcorrectionMValue%i%s[SScounter%i]=atof(segmentMtop);',thisNumber,'vertical',thisNumber);
                    loadStepScan{end+1}=sprintf('\t\t\tcorrectionMValue%i%s[SScounter%i]=atof(segmentMleft);',thisNumber,'inside',thisNumber);
                    loadStepScan{end+1}=sprintf('\t\t\tcorrectionMValue%i%s[SScounter%i]=atof(segmentMright);',thisNumber,'outside',thisNumber);
            end
            loadStepScan{end+1}=sprintf('\t\t\tSScounter%i += 1;',thisNumber);
            loadStepScan{end+1}=sprintf('\t\t}');
            
            
        end
   end
end

for i = 1:length(loadStepScan)
   fprintf('%s\n',loadStepScan{i}) 
end


loadStepScan{end+1}=sprintf('\t\t');
% loadStepScan{end+1}=sprintf('\t\t');
% loadStepScan{end+1}=sprintf('\t\t');
% loadStepScan{end+1}=sprintf('\t\t');
% loadStepScan{end+1}=sprintf('\t\t');
loadStepScan{end+1}=sprintf('\t}');
loadStepScan{end+1}=sprintf('}');


end

