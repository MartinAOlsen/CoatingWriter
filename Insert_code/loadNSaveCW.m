%% get price from summary
Summary=fileread(['CoatingWriter_output' scanname '.txt'])
i=0;
while i<length(Summary);
    i=i+1;
    if strcmp(Summary(i:i+4),'Price')
        id_in=i;
        in=1;
    end
    if exist('in')
       if strcmp(Summary(i),',') || strcmp(Summary(i),'-')
           id_end=i-1;
           break
       end
    end
end
eval(Summary(id_in:id_end))



%% Write data to result-file
cd ..
s = dir('ScanResults.mat');
if length(s)>0
    load('ScanResults.mat');
    makeNew=0;
else
    Result=struct;
    makeNew=1;
end

Result_this.instrumentName=instrument_name;   % instrument
if isfield(p,'totalPrice')
    Result_this.demandPrice=str2num(p.totalPrice);        % Price
end
Result_this.resultPrice=Price;        % Price
Result_this.parameterList=p; 
Result_this.Monitor=monitor_ALLW;

if makeNew==0
    Result(end+1)=Result_this;
else
    Result=Result_this;
end

save('ScanResults.mat','Result');
if select==1  && strcmp(mode,'pricescan') == 0 %Dont do this in price scan since it will clear variables...
    try
        run('AnalyseScan.m')
    end
end
