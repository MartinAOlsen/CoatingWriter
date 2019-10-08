function plot_guide(component,mValue, mLim, cmap,toCut)
% plots guide geometry from GuideBot/McStas structured data
%
% PLOT_GUIDE(component,mValue, mLim, cmap,toCut)
%
% Plots guide geometry with m-values onto an existing figure. The function
% uses path() for fast plotting.
%
% Input:
%
% component     Components imported using the guide_reader() function.
% mValue        mValue imported using the read_coatingwriter() function.
% mLim          Limits for the mvalue for the color plot.
% cmap          Colormap to show the mvalue with colors.
% toCut         If true, the ellipses are cut at the connection points,
%               default is true.
%

% use the existing figure
nP = 2001;
col = get(gca,'colororder');
lineSpec  = {'lineStyle' '-','LineWidth',3};
patchSpec = {'lineStyle' '-','LineWidth',3,'EdgeColor','flat'};
% cut extra ellipse parts
if nargin < 5
    toCut = true;
end

% colormap
if nargin < 4
    cmap = jet;
end

% If Los calculation, draw line
LoS = cellfun(@(C)C.los,component,'uniformoutput',false);
if find(ismember(LoS,'end'))

    % combine components around LoS markers
    
    sIdx = find(ismember(LoS,'start'));
    eIdx = find(ismember(LoS,'end'));

    los.hA = component{sIdx}.x1(3:4)';
    los.hB = component{sIdx}.x1(1:2)';
    los.hC = component{eIdx}.x1(3:4)';
    los.hD = component{eIdx}.x1(1:2)';

    los.vA = component{sIdx}.y1(3:4)';
    los.vB = component{sIdx}.y1(1:2)';
    los.vC = component{eIdx}.y1(3:4)';
    los.vD = component{eIdx}.y1(1:2)';

    % LoS start combine previous component with current
    cpre = component{sIdx-1};
    curr = component{sIdx};
    curr.x1 = cpre.x1;
    curr.y1 = cpre.y1;
    component{sIdx} = curr;

    cpre = component{eIdx-1};
    curr = component{eIdx};

    cpre.x2 = curr.x2;
    cpre.y2 = curr.y2;
    component{eIdx-1} = cpre;

    component([sIdx-1 eIdx]) = [];



    % draw los
    subplot(2,1,1)
    line([los.hA(1) los.hB(1)],[los.hA(2) los.hB(2)],'color','r','linewidth',2)
    hold on
    line([los.hC(1) los.hD(1)],[los.hC(2) los.hD(2)],'color','g','linewidth',2)
    subplot(2,1,2)
    line([los.vA(1) los.vB(1)],[los.vA(2) los.vB(2)],'color','r','linewidth',2)
    hold on
    line([los.vC(1) los.vD(1)],[los.vC(2) los.vD(2)],'color','g','linewidth',2)

end
% counter for m value list
midx = 1;

% max/min m-values
if (nargin < 3 || isempty(mLim)) && nargin>1
    mmin = cellfun(@(C)min(reshape(C.dat(:,3:6),1,[])),mValue,'uniformoutput',false);
    mmax = cellfun(@(C)max(reshape(C.dat(:,3:6),1,[])),mValue,'uniformoutput',false);
    mLim = [min([mmin{:}]) max([mmax{:}])];
end

if nargin < 2
    plotm = false;
else
    plotm = true;
end


for ii = 1:numel(component)
    csel = component{ii};
    
    switch csel.type
        case 'E'
            % draw ellipse
            % ellipse coordinate transformation
            h  = csel.xdir';
            h1 = h;
            h2 = [-h(2);h(1)];
            Th = [h1 h2];
            
            v  = csel.ydir';
            v1 = v;
            v2 = [-v(2);v(1)];
            Tv = [v1 v2];
            
            % length of component
            L = csel.xpars(1);
            
            % focus point in global coordinate system
            F1h = csel.xpos' - csel.xpars(3)*h;
            F2h = csel.xpos' + h*L + csel.xpars(4)*h;
            F1v = csel.ypos' - csel.ypars(3)*v;
            F2v = csel.ypos' + v*L + csel.ypars(4)*v;
            
            % center of ellipse
            Oh  = (F1h+F2h)/2;
            Ov  = (F1v+F2v)/2;
            
            % ellipse parameters
            ch = norm(F2h-F1h)/2;
            bh = csel.xpars(2);
            ah = norm([bh ch]);
            cv = norm(F2v-F1v)/2;
            bv = csel.ypars(2);
            av = norm([bv cv]);
            
            % draw ellipse in local coordinate system
            xh  = linspace(-ah,ah,nP);
            yh  = bh/ah*sqrt(ah^2-xh.^2);
            xv  = linspace(-av,av,nP);
            yv  = bv/av*sqrt(av^2-xv.^2);
            
            % end points
            hA = csel.x1(3:4)';
            hB = csel.x1(1:2)';
            hC = csel.x2(3:4)';
            hD = csel.x2(1:2)';
            hP = [hA hB hC hD];
            vA = csel.y1(3:4)';
            vB = csel.y1(1:2)';
            vC = csel.y2(3:4)';
            vD = csel.y2(1:2)';
            vP = [vA vB vC vD];
            
            % convert end points to ellipse coordinate system
            rPloch = inv(Th)*bsxfun(@minus,hP,Oh); %#ok<MINV>
            rPlocv = inv(Tv)*bsxfun(@minus,vP,Ov); %#ok<MINV>
            
            % upper curve
            if toCut
                selu = (xh >= rPloch(1,1)) & (xh<= rPloch(1,3));
                selb = (xh >= rPloch(1,2)) & (xh<= rPloch(1,4));
                sell = (xv >= rPlocv(1,1)) & (xv<= rPlocv(1,3));
                selr = (xv >= rPlocv(1,2)) & (xv<= rPlocv(1,4));
            else
                selu = true(size(xh));
                selb = true(size(xh));
                sell = true(size(xv));
                selr = true(size(xv));
            end
            
            rE1h = [xh(selu); yh(selu)];
            rE2h = [xh(selb);-yh(selb)];
            rEh  = [rE1h nan(2,1) rE2h(:,end:-1:1)];
            rE1v = [xv(sell); yv(sell)];
            rE2v = [xv(selr);-yv(selr)];
            rEv  = [rE1v nan(2,1) rE2v(:,end:-1:1)];
            
            
            %rE = [[x x(end:-1:1)];[y -y(end:-1:1)]];
            % transform to global coordinate system
            rGh = bsxfun(@plus,Th*rEh,Oh);
            rGv = bsxfun(@plus,Tv*rEv,Ov);
            
            % draw ellipse
            if plotm
                % horizontal
                subplot(2,1,1)
                nanidx = find(isnan(rGh(1,:)));
                mValue{midx}.mRight=flip(mValue{midx}.mRight);
                if min(mValue{midx}.mRight)==max(mValue{midx}.mRight)
                    lCData = cval(max(mValue{midx}.mRight),nanidx);
                else
                    lCData = cval(mValue{midx}.mRight,nanidx);
                end
                
                patch(rGh(1,1:nanidx),rGh(2,1:nanidx),lCData,patchSpec{:})
%                 mValue{midx}.mRight=flip(mValue{midx}.mRight);
%                 mValue{midx}.mLeft=flip(mValue{midx}.mLeft);
                if min(mValue{midx}.mLeft)==max(mValue{midx}.mLeft)
                    lCData = cval(max(mValue{midx}.mLeft),nanidx);
                else
                    lCData = cval(mValue{midx}.mLeft,nanidx);
                end
                patch(rGh(1,nanidx:end),rGh(2,nanidx:end),lCData,patchSpec{:})
                
                % vertical
                subplot(2,1,2)
                nanidx = find(isnan(rGv(1,:)));
                mValue{midx}.mTop=flip(mValue{midx}.mTop);
                  if min(mValue{midx}.mTop)==max(mValue{midx}.mTop)
                    lCData = cval(max(mValue{midx}.mTop),nanidx);
                else
                    lCData = cval(mValue{midx}.mTop,nanidx);
                end
                patch(rGv(1,1:nanidx),rGv(2,1:nanidx),lCData,patchSpec{:})
                mValue{midx}.mBottom=flipud(mValue{midx}.mBottom);
                if min(mValue{midx}.mBottom)==max(mValue{midx}.mBottom)
                    lCData = cval(max(mValue{midx}.mBottom),nanidx);
                else
                    lCData = cval(mValue{midx}.mBottom,nanidx);
                end
                patch(rGv(1,nanidx:end),rGv(2,nanidx:end),lCData,patchSpec{:})
                
            else
                % find the separator between the 2 sides
                subplot(2,1,1)
                line(rGh(1,:),rGh(2,:),lineSpec{:},'color',col(ii,:));
                subplot(2,1,2)
                line(rGv(1,:),rGv(2,:),lineSpec{:},'color',col(ii,:));
            end
            
            % draw end points
            subplot(2,1,1)
            % end points
%             plot3(hP(1,:),hP(2,:),hP(2,:)*0+10,'ok','markerfacecolor',col(ii,:))
            % focus points
%             plot3([F1h(1) F2h(1)],[F1h(2) F2h(2)],[10 10],'ko','markerfacecolor',col(ii,:))
            
            subplot(2,1,2)
            % end points
%             plot3(vP(1,:),vP(2,:),vP(2,:)*0+10,'ok','markerfacecolor',col(ii,:))
            % focus points
%             plot3([F1v(1) F2v(1)],[F1v(2) F2v(2)],[10 10],'ko','markerfacecolor',col(ii,:))
            
            % set the line color with colormap
            midx = midx+1;
        case 'S'
            % draw curved section
            % entrance/exit position
            hA = csel.x1(3:4)';
            hB = csel.x1(1:2)';
            hC = csel.x2(3:4)';
            hD = csel.x2(1:2)';
            vA = csel.y1(3:4)';
            vB = csel.y1(1:2)';
            vC = csel.y2(3:4)';
            vD = csel.y2(1:2)';
            
            % draw 
            if plotm
                subplot(2,1,1)
                straightXline1=linspace(csel.x1(1),csel.x2(1),nP);
                straightXline2=linspace(csel.x1(2),csel.x2(2),nP);
                straightXline3=linspace(csel.x1(3),csel.x2(3),nP);
                straightXline4=linspace(csel.x1(4),csel.x2(4),nP);
                straightYline1=linspace(csel.y1(1),csel.y2(1),nP);
                straightYline2=linspace(csel.y1(2),csel.y2(2),nP);
                straightYline3=linspace(csel.y1(3),csel.y2(3),nP);
                straightYline4=linspace(csel.y1(4),csel.y2(4),nP);
                
                
                lCData = cval(mValue{midx}.mRight,nP);
                patch([straightXline3 nan],[straightXline4 nan],[lCData lCData(end)],patchSpec{:})
                lCData = cval(mValue{midx}.mLeft,nP);
                patch([straightXline1 nan],[straightXline2 nan],[lCData lCData(end)],patchSpec{:})
                
                subplot(2,1,2)
                
                lCData = cval(mValue{midx}.mTop,nP);
                patch([straightYline3 nan],[straightYline4 nan],[lCData lCData(end)]',patchSpec{:})
                lCData = cval(mValue{midx}.mBottom,nP);
                patch([straightYline1 nan],[straightYline2 nan],[lCData lCData(end)]',patchSpec{:})
                
            else
                subplot(2,1,1)
                line([rh1(1,:) nan rh2(1,:)],[rh1(2,:) nan rh2(2,:)],lineSpec{:},'color',col(ii,:));
                subplot(2,1,2)
                line([vA(1) vC(1) nan vB(1) vD(1)],[vA(2) vC(2) nan vB(2) vD(2)],lineSpec{:},'color',col(ii,:));
                
            end
            
            midx = midx+1;
        case 'C'            
            % draw straight section
            % entrance/exit position
            hA = csel.x1(3:4)';
            hB = csel.x1(1:2)';
            hC = csel.x2(3:4)';
            hD = csel.x2(1:2)';
            vA = csel.y1(3:4)';
            vB = csel.y1(1:2)';
            vC = csel.y2(3:4)';
            vD = csel.y2(1:2)';
            
            
            if isfield(csel,'x')
                rot='H'
                % circle center
                Oh = csel.x(3:4)';
                % circle radius
                R2 = csel.x(2);
                R1 = 2*csel.x(1)-R2;
            else
                rot='V'
                % circle center
                Oh = csel.y(3:4)';
                % circle radius
                R2 = csel.y(2);
                R1 = 2*csel.y(1)-R2;
            end
            

            
            % draw circle
            if plotm
                if strcmp(rot,'H')
                    subplot(2,1,1)
                     % circle angle vectors
                    phi1v = linspace(atan2v(hA-Oh),atan2v(hC-Oh),nP);
                    phi2v = linspace(atan2v(hB-Oh),atan2v(hD-Oh),nP);

                    rh1 = bsxfun(@plus,R1*[sin(phi1v);cos(phi1v)],Oh);
                    rh2 = bsxfun(@plus,R2*[sin(phi2v);cos(phi2v)],Oh);
                    lCData = cval(mValue{midx}.mRight,size(rh1,2));
                    patch([rh1(1,:) nan],[rh1(2,:) nan],[lCData lCData(end)],patchSpec{:})
                    lCData = cval(mValue{midx}.mLeft,size(rh2,2));
                    patch([rh2(1,:) nan],[rh2(2,:) nan],[lCData lCData(end)],patchSpec{:})

                    subplot(2,1,2)

                    lCData = cval(mValue{midx}.mTop,2);
                    patch([vA(1) vC(1) nan],[vA(2) vC(2) nan],[lCData lCData(end)]',patchSpec{:})
                    lCData = cval(mValue{midx}.mBottom,2);
                    patch([vB(1) vD(1) nan],[vB(2) vD(2) nan],[lCData lCData(end)]',patchSpec{:})
                else
                    phi1v = linspace(atan2v(vA-Oh),atan2v(vC-Oh),nP);
                    phi2v = linspace(atan2v(vB-Oh),atan2v(vD-Oh),nP);
                    rh1 = bsxfun(@plus,R1*[sin(phi1v);cos(phi1v)],Oh);
                    rh2 = bsxfun(@plus,R2*[sin(phi2v);cos(phi2v)],Oh);

                    subplot(2,1,1)
                    
                    lCData = cval(mValue{midx}.mTop,2);
                    patch([hA(1) hC(1) nan],[hA(2) hC(2) nan],[lCData lCData(end)]',patchSpec{:})
                    lCData = cval(mValue{midx}.mBottom,2);
                    patch([hB(1) hD(1) nan],[hB(2) hD(2) nan],[lCData lCData(end)]',patchSpec{:})

                    subplot(2,1,2)

                    lCData = cval(mValue{midx}.mRight,size(rh1,2));
                    patch([rh1(1,:) nan],[rh1(2,:) nan],[lCData lCData(end)]',patchSpec{:})
                    lCData = cval(mValue{midx}.mLeft,size(rh1,2));
                    patch([rh2(1,:) nan],[rh2(2,:) nan],[lCData lCData(end)]',patchSpec{:})
                end                
            else
                subplot(2,1,1)
                line([rh1(1,:) nan rh2(1,:)],[rh1(2,:) nan rh2(2,:)],lineSpec{:},'color',col(ii,:));
                subplot(2,1,2)
                line([vA(1) vC(1) nan vB(1) vD(1)],[vA(2) vC(2) nan vB(2) vD(2)],lineSpec{:},'color',col(ii,:));
                
            end
            
            midx = midx+1;
        case 'G'
            midx = midx + 1;
        case 'Sample'
            
            h  = csel.xd';
            h1 = h;
            h2 = [-h(2);h(1)];
            Th = [h1 h2];
            
            v  = csel.yd';
            v1 = v;
            v2 = [-v(2);v(1)];
            Tv = [v1 v2];
            
            Oh = csel.xp';
            Ov = csel.yp';
            
            r1 = bsxfun(@plus,Th*[csel.r(3)*[1 1];csel.r(1)*[-1 1]/2],Oh);
            r2 = bsxfun(@plus,Tv*[csel.r(3)*[1 1];csel.r(2)*[-1 1]/2],Ov);
            
            
            subplot(2,1,1)
            line(r1(1,:),r1(2,:),'color',col(end,:),'linewidth',4)
            axis tight
            subplot(2,1,2)
            line(r2(1,:),r2(2,:),'color',col(end,:),'linewidth',4)
            axis tight
        case 'Moderator'
            
            subplot(2,1,1)
            line([0 0],csel.r(1)*[-1/2 1/2],'color',col(end,:),'linewidth',4)
            axis tight
            subplot(2,1,2)
            line([0 0],csel.r(2)*[-1/2 1/2],'color',col(end,:),'linewidth',4)
            axis tight
            
            
        case 'P'
            % draw ellipse
            % ellipse coordinate transformation
            h  = csel.xdir';
            h1 = h;
            h2 = [-h(2);h(1)];
            Th = [h1 h2];
            
            v  = csel.ydir';
            v1 = v;
            v2 = [-v(2);v(1)];
            Tv = [v1 v2];
            
            % length of component
            L = csel.xpars(1);
            
            % focus point in global coordinate system
            F1h = csel.xpos' - csel.xpars(3)*h;
            F2h = csel.xpos' + h*L + csel.xpars(4)*h;
            F1v = csel.ypos' - csel.ypars(3)*v;
            F2v = csel.ypos' + v*L + csel.ypars(4)*v;
            
            % center of ellipse
            Oh  = (F1h+F2h)/2;
            Ov  = (F1v+F2v)/2;
            
            % ellipse parameters
            ch = norm(F2h-F1h)/2;
            bh = csel.xpars(2);
            ah = norm([bh ch]);
            cv = norm(F2v-F1v)/2;
            bv = csel.ypars(2);
            av = norm([bv cv]);
            
            % draw ellipse in local coordinate system
            xh  = linspace(-ah,ah,nP);
            yh  = bh/ah*sqrt(ah^2-xh.^2);
            xv  = linspace(-av,av,nP);
            yv  = bv/av*sqrt(av^2-xv.^2);
            
            % end points
            hA = csel.x1(3:4)';
            hB = csel.x1(1:2)';
            hC = csel.x2(3:4)';
            hD = csel.x2(1:2)';
            hP = [hA hB hC hD];
            vA = csel.y1(3:4)';
            vB = csel.y1(1:2)';
            vC = csel.y2(3:4)';
            vD = csel.y2(1:2)';
            vP = [vA vB vC vD];
            
            % convert end points to ellipse coordinate system
            rPloch = inv(Th)*bsxfun(@minus,hP,Oh); %#ok<MINV>
            rPlocv = inv(Tv)*bsxfun(@minus,vP,Ov); %#ok<MINV>
            
            % upper curve
            if toCut
                selu = (xh >= rPloch(1,1)) & (xh<= rPloch(1,3));
                selb = (xh >= rPloch(1,2)) & (xh<= rPloch(1,4));
                sell = (xv >= rPlocv(1,1)) & (xv<= rPlocv(1,3));
                selr = (xv >= rPlocv(1,2)) & (xv<= rPlocv(1,4));
            else
                selu = true(size(xh));
                selb = true(size(xh));
                sell = true(size(xv));
                selr = true(size(xv));
            end
            
            rE1h = [xh(selu); yh(selu)];
            rE2h = [xh(selb);-yh(selb)];
            rEh  = [rE1h nan(2,1) rE2h(:,end:-1:1)];
            rE1v = [xv(sell); yv(sell)];
            rE2v = [xv(selr);-yv(selr)];
            rEv  = [rE1v nan(2,1) rE2v(:,end:-1:1)];
            
            
            %rE = [[x x(end:-1:1)];[y -y(end:-1:1)]];
            % transform to global coordinate system
            rGh = bsxfun(@plus,Th*rEh,Oh);
            rGv = bsxfun(@plus,Tv*rEv,Ov);
            
            % draw ellipse
            if plotm
                % horizontal
                subplot(2,1,1)
                nanidx = find(isnan(rGh(1,:)));
                mValue{midx}.mRight=flip(mValue{midx}.mRight);
                if min(mValue{midx}.mRight)==max(mValue{midx}.mRight)
                    lCData = cval(max(mValue{midx}.mRight),nanidx);
                else
                    lCData = cval(mValue{midx}.mRight,nanidx);
                end
                
                patch(rGh(1,1:nanidx),rGh(2,1:nanidx),lCData,patchSpec{:})
%                 mValue{midx}.mRight=flip(mValue{midx}.mRight);
%                 mValue{midx}.mLeft=flip(mValue{midx}.mLeft);
                if min(mValue{midx}.mLeft)==max(mValue{midx}.mLeft)
                    lCData = cval(max(mValue{midx}.mLeft),nanidx);
                else
                    lCData = cval(mValue{midx}.mLeft,nanidx);
                end
                patch(rGh(1,nanidx:end),rGh(2,nanidx:end),lCData,patchSpec{:})
                
                % vertical
                subplot(2,1,2)
                nanidx = find(isnan(rGv(1,:)));
                mValue{midx}.mTop=flip(mValue{midx}.mTop);
                  if min(mValue{midx}.mTop)==max(mValue{midx}.mTop)
                    lCData = cval(max(mValue{midx}.mTop),nanidx);
                else
                    lCData = cval(mValue{midx}.mTop,nanidx);
                end
                patch(rGv(1,1:nanidx),rGv(2,1:nanidx),lCData,patchSpec{:})
                mValue{midx}.mBottom=flipud(mValue{midx}.mBottom);
                if min(mValue{midx}.mBottom)==max(mValue{midx}.mBottom)
                    lCData = cval(max(mValue{midx}.mBottom),nanidx);
                else
                    lCData = cval(mValue{midx}.mBottom,nanidx);
                end
                patch(rGv(1,nanidx:end),rGv(2,nanidx:end),lCData,patchSpec{:})
                
            else
                % find the separator between the 2 sides
                subplot(2,1,1)
                line(rGh(1,:),rGh(2,:),lineSpec{:},'color',col(ii,:));
                subplot(2,1,2)
                line(rGv(1,:),rGv(2,:),lineSpec{:},'color',col(ii,:));
            end
            
            % draw end points
            subplot(2,1,1)
            % end points
%             plot3(hP(1,:),hP(2,:),hP(2,:)*0+10,'ok','markerfacecolor',col(ii,:))
            % focus points
%             plot3([F1h(1) F2h(1)],[F1h(2) F2h(2)],[10 10],'ko','markerfacecolor',col(ii,:))
            
            subplot(2,1,2)
            % end points
%             plot3(vP(1,:),vP(2,:),vP(2,:)*0+10,'ok','markerfacecolor',col(ii,:))
            % focus points
%             plot3([F1v(1) F2v(1)],[F1v(2) F2v(2)],[10 10],'ko','markerfacecolor',col(ii,:))
            
            % set the line color with colormap
            midx = midx+1;    
        otherwise
            error('guide_plotter:WrongInput','The given component type: ''%s'' is not implemented yet',component{ii}.type)
    end
    
end

% add colorbar and labels
if plotm
    subplot(2,1,1)
    box on
    grid on
    colormap(cmap);
    hColorbar = colorbar;
    hColorbar.Label.String = 'm-values';
    caxis(mLim)
    title('Horizontal plane')
    xlabel('z (m)')
    ylabel('x (m)')
    x0 = xlim;
    xlim([-0.5 x0(2)])
    subplot(2,1,2)
    box on
    grid on
    colormap(cmap);
    hColorbar = colorbar;
    hColorbar.Label.String = 'm-values';
    xlabel('z (m)')
    ylabel('y (m)')
    title('Vertical plane')
    xlim([-0.5 x0(2)])
    caxis(mLim)
end

end

function phi = atan2v(r)
% atan2 of 2 element vector

phi = atan2(r(1),r(2));

end

function xval = cval(x,nPoint)
% create line colordata

% create color for each line piece
if numel(x)==1
    xval = repmat(x,1,nPoint);
else
    list=round(linspace(0.5,numel(x)+1/2-10*eps,nPoint));
    if max(list)>length(x)
       list(end)=length(x); 
    end
    xval = x(:,list);
end

end