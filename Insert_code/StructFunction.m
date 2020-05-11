function pReeval=StructFunction(listofpar,Runpars,p,pars_i);
if length(Runpars)<length(listofpar);
    for k=(length(listofpar)-length(Runpars)):length(listofpar)
        Runpars(k)=0;
    end
end

for k=1:length(listofpar)
    if Runpars(k)==0;
        if isa(eval(sprintf('p.%s',char(listofpar(k)))),'double')
            eval(sprintf('pReeval.%s=''%d'';',char(listofpar(k)),eval(sprintf('p.%s(2)',char(listofpar(k))))))
        else
            eval(sprintf('pReeval.%s=''%s'';',char(listofpar(k)),eval(sprintf('p.%s',char(listofpar(k))))))
        end
    else
        eval(sprintf('pReeval.%s=''%d'';',char(listofpar(k)),eval(sprintf('pars_i.%s',char(listofpar(k))))))
    end
end



end