function CW_output_excel(filename,varargin)
%% Begin
fprintf('________________________\nConverting CW output to excel file\n________________________\n');
if nargin == 0
   listing=dir;
   counter=0;
   for i=1:length(listing)
       if strfind(listing(i).name,'CoatingWriter_rawData')     
            counter=counter+1;
            filename=listing(i).name;
       end
   end
   fprintf('Found %i output files\n',counter)
end
fprintf('using file: %s\n\n',filename);

%% Load CW file:

F=fopen(filename);
fgetl(F); %discard first line
List=[]
while true
    tmp=fgetl(F);
    if length(tmp)<2
       break 
    end
    tmp=strrep(tmp,',',' ');
    tmp=strrep(tmp,'BK','1');tmp=strrep(tmp,'BF','2');tmp=strrep(tmp,'RH','3');tmp=strrep(tmp,'Na','4');
    tmp=strrep(tmp,'P','1');tmp=strrep(tmp,'G','2');tmp=strrep(tmp,'K','3');tmp=strrep(tmp,'C','4');tmp=strrep(tmp,'S','5');tmp=strrep(tmp,'E','6');
    List(end+1,:)=strread(tmp);
end

%% Sort List
oldList=List;
clear List;
List=[];
for i=max(oldList(:,1)):-1:1
    for j=1:length(oldList)
        if oldList(j,1)== i
            List(end+1,:)=oldList(j,:);
        end
    end
end


%% Component classification:
for i = 1:max(oldList(:,1))
    LocalList=[];
    for k=1:length(List(:,1))
        if List(k,1)==i
            LocalList(end+1,:)=List(k,:);
        end
    end
    
    
    if LocalList(1,2)==1 %% Parabolic
        if i==max(oldList(:,1)) && (LocalList(1,4)+List(1,6))>(LocalList(end,5)+List(end,7))*1.5 && sum(LocalList(:,3))<15
            component{i}='Parabolic Feeder';
        elseif (LocalList(1,4)+List(1,6))>(LocalList(end,5)+List(end,7))*1.5
            component{i}='Expanding Parabolic';
        elseif (LocalList(1,4)+List(1,6))*1.5<(LocalList(end,5)+List(end,7))
            component{i}='Focusing Parabolic';
        else
            component{i}='Parabolic'; 
        end
    
    elseif LocalList(1,2)==2 %% Gap
        component{i}='Gap'; 
        
    elseif LocalList(1,2)==3 %% Kink
        component{i}='Gap (kink)'; 
    elseif LocalList(1,2)==4 %% Curve
        component{i}='Curved section'; 
        
    elseif LocalList(1,2)==5 %% Straight
        if LocalList(1,4)==LocalList(end,5) && List(1,6) == List(end,7)
            component{i}='True Straight'; 
        else
            component{i}='Straight'; 
        end
        
    elseif LocalList(1,2)==6 %% Ellipse
        if i==max(oldList(:,1)) && (LocalList(1,4)+List(1,6))>(LocalList(end,5)+List(end,7))*1.5 && sum(LocalList(:,3))<15
            component{i}='Elliptic Feeder';
        elseif (LocalList(1,4)+List(1,6))>(LocalList(end,5)+List(end,7))*1.5
            component{i}='Expanding Ellipse';
        elseif (LocalList(1,4)+List(1,6))*1.5<(LocalList(end,5)+List(end,7))
            component{i}='Focusing Ellipse';
        else
            component{i}='Ellipse'; 
        end
    else
        component{i}='unknown';
    end
       
    
end



%% Shorten and make table
A={'Component','Length [m]','Substrate','Coating left(outside)','Coating right(inside)','Coating top','Coating bottom','Distance at end [m]','Height-1 [cm]','Height-2','Width-1','Width-2'};
NewLine=1;
 A(end+1,1)={'Moderator'}
for i=1:length(List)
   
    
    switch List(i,12)
        case 1
            substrate='Borkron';
        case 2
            substrate='Borfloat';
        case 3    
            substrate='Metal';
        case 4    
            substrate='Sodium';
    end
    if i>1 && NewLine==0
       if List(i,1) ~= List(i-1,1); NewLine=1; end %If new segment
       if List(i,8) ~= List(i-1,8); NewLine=1; end % If new coating
       if List(i,9) ~= List(i-1,9); NewLine=1; end
       if List(i,10) ~= List(i-1,10); NewLine=1; end
       if List(i,11) ~= List(i-1,11); NewLine=1; end
    end
    
    if NewLine==1
       if i>1
           A(end+1,:)={component(List(i-1,1)),Length,substrate,mLeft,mRight,mTop,mBottom,'0',startY,endY,startX,endX}
       end
        
       NewLine=0;
       segNr=List(i,1); 
       startY=List(i,4)*100;
       startX=List(i,6)*100;
       Length=List(i,3);
       endY=List(i,5)*100;
       endX=List(i,7)*100;
       mTop=List(i,8);
       mBottom=List(i,9);
       mLeft=List(i,10);
       mRight=List(i,11);
    else
       endY=List(i,5)*100;
       endX=List(i,7)*100;
       Length=Length+List(i,3);
    end
    
    
    
    if i== length(List) % Print last
        A(end+1,:)={component(List(i,1)),Length,substrate,mLeft,mRight,mTop,mBottom,'0',startY,endY,startX,endX};
    end
    
end
A(end+1,1)={'Sample'}


%% Export to excel
% make all cells to strings
for i=1:size(A,1)
    for j=1:size(A,2)
        if isfloat(A{i,j}) == 1
            A{i,j}=num2str(A{i,j});
        end
        if iscell(A{i,j}) == 1
            A{i,j}=char(A{i,j});
        end
    end
end
xlswrite([strrep(filename(1:end-4),'rawData','easyToRead') '.xlsx'],char(A),1,'A1')

%% Export to .csv (in case excel is not installed)

csvF=fopen([strrep(filename(1:end-4),'rawData','easyToRead') '.csv'],'w');
fprintf(csvF,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',A{1,:});
for i=3:size(A,1)-1
    fprintf(csvF,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',A{i,:});
end
fclose(csvF);

end

