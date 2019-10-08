function out=maxLine(x,y)
    % needs a step at 1/stepN x:


    C = [x,y];


    

    C_sorted = sortrows(C);
    
    Xt=C_sorted(:,1);
    Yt=C_sorted(:,2);
    
    nSteps=8;%round(length(Xt)/100);
    iterSteps=floor(length(Xt)/nSteps);
    Y(1)=Yt(1);
    X(1)=Xt(1);
    for i=1:nSteps
        [Y(i+1),id]=max(Yt(1+(i-1)*iterSteps:i*iterSteps))
        X(i+1)=Xt((i-1)*iterSteps+id);
    end
    if nSteps*iterSteps<length(Xt)
        [Y(end+1),id]=max(Yt(1+nSteps*iterSteps:end))
        X(end+1)=Xt(nSteps*iterSteps+id);
    end
    if X(end)<Xt(end)
        Y(end+1)=Yt(end);
        X(end+1)=Xt(end);
    end
    
    method=2
    if method==1
        for i=1:15
            idmax=round((i/15)*length(X))-1;
            idmin=round(((i-1)/15)*length(X))+1;
            [M,I] = max(Y(idmin:idmax));
            out(i,:)=[X(I+idmin),M];
        end
    elseif method==2
        prevMax=0;
        List(1)=1;
        List(2)=2;
        List(3)=3;
        for i=4:length(X)
           % Draw line to latest maxes
           canKill=[];
          
           
          for j=2:length(List) 
              ListTmp=List;
              %calculate line between points
              coefficients = polyfit([X(List(j-1)), X(i)], [Y(List(j-1)), Y(i)], 1);
              a = coefficients (1) ;             
              b = coefficients (2);
              
%                for k=i:-1: 1
                  if Y(List(j))<= a*(X(List(j)))+b 
                      canKill(end+1)=(j);
                      j=length(List) ;
                  end
%               end
 
          end
          %if length(canKill)>0
          %   ListTmp(min(canKill):end)=[];  % remove non-optimal points
          %end
          List=ListTmp;
          List(end+1)=i;
        end
    end
    out(1,:)=[X(1),Y(1)]
    for i=1:length(List)
        out(i+1,:)=[X(List(i)),Y(List(i))];
    end
    
end