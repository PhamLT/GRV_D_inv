function GRV_D_inv
%%% Code by Oksum E. Pham. L.T.,(2021) 
%%% Main GUI builder function
clc;clear all;clear global;delete(findobj('Type','figure'));
f=emptW([0 0 .15 1]);
set(f,'Name','GRV_D_inv')
m1=uicontrol('Parent',f,'string','MODE-INVERSION','units','normalized',...
    'position',[0 0.5 1 .5],'CallBack',{@selectmode,-1});
m2=uicontrol('Parent',f,'string','MODE-FORWARD','units','normalized',...
    'position',[0 0 1 .5],'CallBack',{@selectmode,1});
set([m1 m2],'Fontweight','bold')
crt_FRW_window
crt_INV_window
save('storeF.mat')
end

function f=emptW(pos)
f=figure('MenuBar','none','NumberTitle','off','DockControls','off',...
'Color','w','units','normalized','outerposition',pos,'resize','off');
end

function selectmode(~,~,modS)
switch modS
    case 1
    r=findobj('Type','figure','name','FORWARD MODE');    
    if isempty(r);crt_FRW_window;else uistack(r,'top');end
    case -1
    r=findobj('Type','figure','name','INVERSION MODE');    
    if isempty(r);crt_INV_window;else uistack(r,'top');end   
end
r=findobj('Type','figure','name','GRV_D_inv');
set(r,'outerposition',[0 0 .15 1])
end

function crt_INV_window
f1=emptW([0.15 0 .85 1]);set(f1,'Name','INVERSION MODE')
M=uimenu('Parent',f1,'Label','LOAD DATA');
m1=uimenu('Parent',M,'Label','Import Grid Data (Gravity file)');
p1=uipanel(f1,'units','normalized','Position',[.01 .95 .12 .05],...
'BackgroundColor','w','Title','Density Contrast (gr/cc)','FontWeight','bold',...
'TitlePosition','lefttop','ForegroundColor','k');
p2=uipanel(f1,'units','normalized','Position',[.13 .95 .12 .05],...
'BackgroundColor','w','Title','Mean Depth (km)','FontWeight','bold',...
'TitlePosition','lefttop','ForegroundColor','k');
p3=uipanel(f1,'units','normalized','Position',[.26 .95 .12 .05],...
'BackgroundColor','w','Title','SET Filter WH','FontWeight','bold',...
'TitlePosition','lefttop','ForegroundColor','b');
p4=uipanel(f1,'units','normalized','Position',[.39 .95 .12 .05],...
'BackgroundColor','w','Title','SET Filter SH','FontWeight','bold',...
'TitlePosition','lefttop','ForegroundColor','r');
p5=uipanel(f1,'units','normalized','Position',[.52 .95 .12 .05],...
'BackgroundColor','w','Title','RMS Criterio','FontWeight','bold',...
'TitlePosition','lefttop','ForegroundColor','k');
p6=uipanel(f1,'units','normalized','Position',[.65 .95 .12 .05],...
'BackgroundColor','w','Title','Max. Iteration','FontWeight','bold',...
'TitlePosition','lefttop','ForegroundColor','k');
ed1=uicontrol(p1,'style','edit','units','normalized','Position',[0.1 .1 .8 .8],...
'FontWeight','bold','BackgroundColor','w','string','0.4');
ed2=uicontrol(p2,'style','edit','units','normalized','Position',[0.1 .1 .8 .8],...
'FontWeight','bold','BackgroundColor','w','string','20');
ed3=uicontrol(p3,'style','edit','units','normalized','Position',[0.1 .1 .8 .8],...
'FontWeight','bold','BackgroundColor','w','string','0.02','Tag','edWH');
ed4=uicontrol(p4,'style','edit','units','normalized','Position',[0.1 .1 .8 .8],...
'FontWeight','bold','BackgroundColor','w','string','0.04','Tag','edSH');
set(ed3,'CallBack',{@setlinWHSH,ed3,'b'})
set(ed4,'CallBack',{@setlinWHSH,ed4,'r'})
ed5=uicontrol(p5,'style','edit','units','normalized','Position',[0.1 .1 .8 .8],...
'FontWeight','bold','BackgroundColor','w','string','0.0001');
ed6=uicontrol(p6,'style','edit','units','normalized','Position',[0.1 .1 .8 .8],...
'FontWeight','bold','BackgroundColor','w','string','100');
p7=uipanel(f1,'units','normalized','Position',[0.01 0.01 .489 .94],...
'BackgroundColor','w','Title','INPUTS','FontWeight','bold',...
'TitlePosition','centertop','ForegroundColor','k');
p8=uipanel(f1,'units','normalized','Position',[0.491 0.01 .489 .94],...
'BackgroundColor','w','Title','OUTPUTS','FontWeight','bold',...
'TitlePosition','centertop','ForegroundColor','k');
ax11=axes('Parent',p7,'units','normalized','Position',...
[0.2 0.15 0.6 0.4]);axis off
ax12=axes('Parent',p7,'units','normalized','Position',...
[0.1 0.65 0.85 0.3]);axis off
ax2=axes('Parent',p8,'units','normalized','Position',...
[0.15 0.15 0.8 0.8]);axis off
c2 = uicontextmenu;
uicontrol(p8,'UIContextMenu',c2,'string','Export (right click)',...
'units','normalized','ForeGroundColor','k','FontWeight','bold',...
'Position',[0,.975,.2,.025]);
outpush21=uimenu('Parent',c2,'Label','ALL as Image','Separator','on','enable','off');
outpush22=uimenu('Parent',c2,'Label','ALL as Data','Separator','on','enable','off');
%%%%%%%%%%
c3 = uicontextmenu;
uicontrol(p8,'UIContextMenu',c3,'string',' SELECT OUTPUT (right click)',...
'units','normalized','ForeGroundColor','r','FontWeight','bold',...
'Position',[.21,.975,.78,.025]);
plo1=uimenu('Parent',c3,'Label','RMS','Separator','on','enable','off');
plo2=uimenu('Parent',c3,'Label','Zcalc','Separator','on','enable','off');
plo3=uimenu('Parent',c3,'Label','Gcalc','Separator','on','enable','off');
plo4=uimenu('Parent',c3,'Label','Gobs-Gcalc','Separator','on','enable','off');
outpushall=[outpush21 outpush22 plo1 plo2 plo3 plo4];
ax1=[ax11 ax12];
set(m1,'Callback',{@call_data,-1,M,ax1,ax2,outpushall})
uicontrol(f1,'style','pushbutton','units','normalized','Position',[0.78 .95 .24 .04],...
'String','INVERSION','enable','off','Callback',{@START_func,-1,ax2,ed1,ed2,ed3,ed4,ed5,ed6,outpushall})
end

function crt_FRW_window
f2=emptW([0.15 0 .85 1]);set(f2,'Name','FORWARD MODE')
M=uimenu('Parent',f2,'Label','LOAD DATA');
m1=uimenu('Parent',M,'Label','Import Grid Data (Depth file)');
p1=uipanel(f2,'units','normalized','Position',[.01 .95 .12 .05],...
'BackgroundColor','w','Title','Density Contrast (gr/cc)','FontWeight','bold',...
'TitlePosition','lefttop','ForegroundColor','b');
p2=uipanel(f2,'units','normalized','Position',[.13 .95 .12 .05],...
'BackgroundColor','w','Title','Mean Depth (km)','FontWeight','bold',...
'TitlePosition','lefttop','ForegroundColor','b');
ed1=uicontrol(p1,'style','edit','units','normalized','Position',[0.1 .1 .8 .8],...
'FontWeight','bold','BackgroundColor','w','string','0.4');
ed2=uicontrol(p2,'style','edit','units','normalized','Position',[0.1 .1 .8 .8],...
'FontWeight','bold','BackgroundColor','w','string','20');
p3=uipanel(f2,'units','normalized','Position',[0.01 0.01 .489 .94],...
'BackgroundColor','w','Title','Depth Grid','FontWeight','bold',...
'TitlePosition','centertop','ForegroundColor','k');
p4=uipanel(f2,'units','normalized','Position',[0.491 0.01 .489 .94],...
'BackgroundColor','w','Title','Calculated Gravity Grid','FontWeight','bold',...
'TitlePosition','centertop','ForegroundColor','k');
ax1=axes('Parent',p3,'units','normalized','Position',...
[0.1 0.15 0.8 0.8]);axis off
ax2=axes('Parent',p4,'units','normalized','Position',...
[0.1 0.15 0.8 0.8]);axis off
c1 = uicontextmenu;c2 = uicontextmenu;
uicontrol(p3,'UIContextMenu',c1,'string','Export (right click)',...
'units','normalized','ForeGroundColor','k',...
'Position',[0,.965,.2,.035]);
outpush1=uimenu('Parent',c1,'Label','Image','Separator','on','enable','off');
uicontrol(p4,'UIContextMenu',c2,'string','Export (right click)',...
'units','normalized','ForeGroundColor','k',...
'Position',[0,.965,.2,.035]);
outpush21=uimenu('Parent',c2,'Label','Image','Separator','on','enable','off');
outpush22=uimenu('Parent',c2,'Label','Grid Data','Separator','on','enable','off');
outpushall=[outpush1 outpush21 outpush22];
set(m1,'Callback',{@call_data,1,M,ax1,ax2,outpushall})
uicontrol(f2,'style','pushbutton','units','normalized','Position',[0.25 .95 .24 .04],...
'String','Forward CALC','enable','off','Callback',{@START_func,1,ax2,ed1,ed2,0,0,0,0,outpushall})
end

function img_export_map(x,y,matrix,unitt,tit,sorcf)
tfg=figure('MenuBar','none','NumberTitle','off','Resize','off',...
'Color','w','units','normalized','outerposition',[0 0 1 1],...
'DockControls','off','visible','off');
contourf(x,y,matrix,18);shading flat;axis equal;axis tight
title(tit,'FontSize',16,'FontWeight','normal')
xlabel('X (km)','FontSize',16,'FontWeight','normal');
ylabel('Y (km)','FontSize',16,'FontWeight','normal');
xlim([min(min(x)) max(max(x))])
ylim([min(min(y)) max(max(y))])
set(gca,'FontSize',16,'FontWeight','normal')
h=colorbar('eastoutside');
title(h,unitt,'FontSize',16,'FontWeight','normal');
set(h,'FontSize',16,'FontWeight','normal');
print(tfg,sorcf, '-dpng', '-r300');
delete(tfg);
end


function sav_all_images(~,~,typ,x_inv,y_inv,zcalc_INVf,gcalc_INVf,gdiff_INVf,data,rms_it)
r=findobj('type','uimenu','ForeGroundColor','g');
fil=get(r,'Label');[~,fil]=fileparts(fil);
[filename, pathname] = uiputfile('png','Export',[fil '_outs']);
sorcf=[pathname filename];[~,sorc]=fileparts(sorcf);
xmin=min(x_inv(:));xmax=max(x_inv(:));
ymin=min(y_inv(:));ymax=max(y_inv(:));
switch typ
    case '*.png'
img_export_map(x_inv,y_inv,-zcalc_INVf,'km','',[pathname sorc '_zcalc.png'])
img_export_map(x_inv,y_inv,gcalc_INVf,'mGal','',[pathname sorc '_gcalc.png'])
img_export_map(x_inv,y_inv,gdiff_INVf,'mGal','',[pathname sorc '_gdiff.png'])
img_export_map(x_inv,y_inv,data,'mGal','',[pathname sorc '_gobs.png'])
img_export_vctr (rms_it,[pathname sorc '_rms.png'])
    case '*.grd'
grdout(zcalc_INVf,xmin,xmax,ymin,ymax,[pathname sorc '_zcalc.grd']) 
grdout(gcalc_INVf,xmin,xmax,ymin,ymax,[pathname sorc '_gcalc.grd']) 
grdout(gdiff_INVf,xmin,xmax,ymin,ymax,[pathname sorc '_gdiff.grd']) 
grdout(data,xmin,xmax,ymin,ymax,[pathname sorc '_gobs.grd']) 
rmsvec=[(1:numel(rms_it))' rms_it'];
save([pathname sorc '_rms.dat'],'rmsvec','-ascii');
end
end
function rms_plot(~,~,ax2,rms_it) % plot the rms graphic 
drawnow;set(gcf,'CurrentAxes',ax2);refresh
plot(rms_it,'-ko','MarkerFaceColor','r','MarkerSize',5);
set(gca,'FontSize',10,'FontWeight','bold')
sent1=['RMS(1)=' num2str(rms_it(1))];
sent2=['RMS(end)=' num2str(rms_it(end))];
title ([sent1 '    ' sent2])
xlabel('Iteration number');ylabel('RMS (km)');
xlim([1 numel(rms_it)])
grid on
pbaspect([1 .9 1]);
end

function img_export_vctr (rms_it,sorcf)
%%% export vector data
tfg=figure('MenuBar','none','NumberTitle','off','Resize','off',...
'Color','w','units','normalized','outerposition',[0 0 1 1],...
'DockControls','off','visible','off');
drawnow;plot(rms_it,'-ko','MarkerFaceColor','r','MarkerSize',5);
xlabel('Iteration number','FontSize',16);ylabel('RMS (km)','FontSize',16);
xlim([1 numel(rms_it)])
grid on
pbaspect([1 .8 1]);
sent1=['RMS(1)=' num2str(rms_it(1))];
sent2=['RMS(end)=' num2str(rms_it(end))];
title ([sent1 '    ' sent2])
drawnow;set(gca,'FontSize',16,'FontWeight','normal')
drawnow;print(tfg,sorcf, '-dpng', '-r300');
delete(tfg);
end

function setlinWHSH(~,~,ed,clr)
x=str2double(get(ed,'string'));
set(findobj('type','line','color',clr),'Xdata',[x x]);
end

function sorc_out(~,~,X,Y,D,unitt,tit,typ) %% export forward output data to file
[filename, pathname] = uiputfile(typ,'Set filename');
sorcf=[pathname filename];
if ischar([pathname filename])
switch typ
    case '*.png'
    img_export_map(X,Y,D,unitt,tit,sorcf)
    case '*.grd'
    xmin=min(X(:));xmax=max(X(:));ymin=min(Y(:));ymax=max(Y(:));
    grdout(D,xmin,xmax,ymin,ymax,sorcf)
end
end
end

function rapsplot(wn,pwr,ax12)
%%% plot of raps
drawnow;set(gcf,'CurrentAxes',ax12);
plot(wn,pwr,'k','LineWidth',2);
pbaspect([1 1 1])
grid on
xlabel('k/2*pi');ylabel('Log(P)');
pbaspect([2 1 1]);
end
%%% Import Data Functions
function call_data(~,~,modeP,mess,ax1,ax2,outpushall) 
%interactively import gravity or depth data
clc;drawnow
[filename, pathname] = uigetfile('*.grd', 'Import Golden Software Binary/Text grid (*.grd)');
if ischar([pathname filename])
switch modeP
    case -1
    [data_inv,x_inv,y_inv,nx_inv,ny_inv,~,~,~,~,dx_inv,dy_inv]=gridform([pathname filename]);
    errcode=error_checker(data_inv,dx_inv,dy_inv,nx_inv,ny_inv,-1);
    if errcode>0;
    uu=uicontrol('Parent',gcf,'style','text','units','normalized',...
        'BackGroundColor','k','ForeGroundColor','w',...
        'position',[0 0 1 1],'string','incompatible data format, please check...',...
        'FontSize',24,'Fontweight','bold');
    pause(3);delete(uu);drawnow        
    return;
    end
    drawnow;delete(findobj('type','line','color','b'));delete(findobj('type','line','color','r'))
    [pwr,wn]= raps_data(data_inv,dx_inv);
    save('storeF.mat','data_inv','x_inv','y_inv','dx_inv','dy_inv','-append');
    drawnow;mapviewer(x_inv,y_inv,data_inv,'mGal','',ax1(1));
    drawnow;rapsplot(wn,pwr,ax1(2));axis on
    xlim([0 max(wn)]);ylim([min(pwr) max(pwr)]);YL=get(gca,'Ylim');
    edwh=findobj(gcf,'Tag','edWH');
    edsh=findobj(gcf,'Tag','edSH');
    L1=line([0 0],YL,'Color','b','linewidth',4);
    L2=line([max(wn) max(wn)],YL,'Color','r','linewidth',4);
    set(L1,'ButtonDownFcn',{@draglin,L1,edwh,ax1(2)})
    set(L2,'ButtonDownFcn',{@draglin,L2,edsh,ax1(2)})
    title('Graphical Option: Click and drag Line (Blue-WH or Red-SH) to set Filter freq.',...
        'FontWeight','bold')
    cla(ax2);drawnow;
    set(outpushall,'enable','off')
    set(mess,'label',filename,'ForeGroundColor','g')
    %%%%%%%%%%%% enable func button on
    fb=findobj('style','pushbutton','string','INVERSION');
    set(fb,'enable','on')
    %%%%%%%%%%%
    case 1
[data_frw,x_frw,y_frw,nx_frw,ny_frw,~,~,~,~,dx_frw,dy_frw]=gridform([pathname filename]);
    errcode=error_checker(data_frw,dx_frw,dy_frw,nx_frw,ny_frw,1);
    index=find(data_frw<0);if numel(index)>1;errcode=5;end
    if errcode>0;
    uu=uicontrol('Parent',gcf,'style','text','units','normalized',...
        'BackGroundColor','k','ForeGroundColor','w',...
        'position',[0 0 1 1],'string','incompatible data format, please check...',...
        'FontSize',24,'Fontweight','bold');
    pause(3);delete(uu);drawnow        
    return;
    end
    save('storeF.mat','data_frw','x_frw','y_frw','dx_frw','dy_frw','-append');
    mapviewer(x_frw,y_frw,-data_frw,'km','',ax1)
    cla(ax2);drawnow;
    set(outpushall(1),'enable','on','Callback',{@sorc_out,x_frw,y_frw,-data_frw,'km','Depth Model','*.png'})
    set(outpushall(2:3),'enable','off')
    set(mess,'label',filename)
    %%%%%%%%%%%% enable func button on
    fb=findobj('style','pushbutton','string','Forward CALC');
    set(fb,'enable','on')
    %%%%%%%%%%%
end
end
end


function draglin(~,~,L,ed,ax12)
drawnow;set(gcf,'CurrentAxes',ax12);
clr=get(L,'color'); 
if clr==[0 0 1];tit='WH';end
if clr==[1 0 0];tit='SH';end
set(L,'visible','off');drawnow
set(ed,'BackgroundColor','y');drawnow
title(['mouse click to approve desired position of [ ' tit ' ]'], 'color',clr);drawnow
set(gcf,'WindowButtonMotionFcn',{@mouseMove,1,L,ed});
set(gcf,'WindowButtonDownFcn',{@mouseMove,0,L,ed});
set(gcf,'Pointer','crosshair');drawnow
end

function mouseMove(~,~,s,L,ed)
 C = get(gca,'currentpoint'); XL=xlim(gca);
 x = C(1);
 if x<XL(1);x=0;end
 if x>XL(2);x=XL(2);end
 if s==1;
     drawnow;
     set(ed,'string',num2str(x))
 else
 set(L,'XData',[x x],'visible','on');set(ed,'string',num2str(x));drawnow 
 set(ed,'BackgroundColor','w');drawnow
 set(gcf,'WindowButtonMotionFcn','');
 set(gcf,'WindowButtonDownFcn','');
 title('')
 set(gcf,'Pointer','arrow');drawnow
 end
 end


function START_func(~,~,modeP,ax2,ed1,ed2,ed3,ed4,ed5,ed6,outpushall) %%% START an inverion or forward procedure
try
switch modeP
   case 1
   z0=str2double(get(ed2,'string'));z0=abs(z0);r0=str2double(get(ed1,'string'));
   load('storeF.mat','data_frw','dx_frw','dy_frw');
   data=data_frw;dx=dx_frw;dy=dy_frw;
   [ny0,nx0]=size(data); 
   n=10;
   hwtb = waitbar(.2,'processing');
   %%% run forward procedure and memorize outputs
   gcalc_FRWf=FORWARD_func(data,nx0,ny0,dx,dy,r0,z0,n); 
   load('storeF.mat','x_frw','y_frw')
   waitbar(.6,hwtb);pause(.1);waitbar(1,hwtb);delete(hwtb);
   drawnow;mapviewer(x_frw,y_frw,gcalc_FRWf,'mGal','',ax2);
   set(outpushall(2),'enable','on','Callback',{@sorc_out,x_frw,y_frw,gcalc_FRWf,'mGal','Model Gravity','*.png'})
   set(outpushall(3),'enable','on','Callback',{@sorc_out,x_frw,y_frw,gcalc_FRWf,'mGal','Model Gravity','*.grd'})
   case -1
   r0=str2double(get(ed1,'string'));z0=str2double(get(ed2,'string'));z0=abs(z0);
   WH=str2double(get(ed3,'string'));SH=str2double(get(ed4,'string'));
   criterio=str2double(get(ed5,'string'));mxit=str2double(get(ed6,'string'));
   n=10;
   WHSH=unique([WH SH]);WH=WHSH(1);SH=WHSH(2);
   load('storeF.mat','data_inv','dx_inv','dy_inv');
   data=data_inv;dx=dx_inv;dy=dy_inv;
   [ny0,nx0]=size(data);
    %%% run inverse procedure and memorize outputs
    [zcalc_INVf,rms_stor]=INVERSION_func(data,nx0,ny0,dx,dy,r0,z0,WH,SH,criterio,mxit);
    gcalc_INVf=FORWARD_func(zcalc_INVf,nx0,ny0,dx,dy,r0,z0,n);
    gdiff_INVf=data-gcalc_INVf;
    if numel(rms_stor)>1
    load('storeF.mat','x_inv','y_inv')
    drawnow;mapviewer(x_inv,y_inv,-zcalc_INVf,'km','Calculated Depth',ax2)
    set(outpushall(1),'enable','on','Callback',{@sav_all_images,'*.png',x_inv,y_inv,zcalc_INVf,gcalc_INVf,gdiff_INVf,data,rms_stor})
    set(outpushall(2),'enable','on','Callback',{@sav_all_images,'*.grd',x_inv,y_inv,zcalc_INVf,gcalc_INVf,gdiff_INVf,data,rms_stor})
    set(outpushall(3),'enable','on','Callback',{@rms_plot,ax2,rms_stor})
    set(outpushall(4),'enable','on','Callback',{@plotout_inv,ax2,x_inv,y_inv,-zcalc_INVf,'km','Calculated Depth'})
    set(outpushall(5),'enable','on','Callback',{@plotout_inv,ax2,x_inv,y_inv,gcalc_INVf,'mGal','Calculated Gravity' })
    set(outpushall(6),'enable','on','Callback',{@plotout_inv,ax2,x_inv,y_inv,gdiff_INVf,'mGal','Gravity Diff.'})
    else
    hwtb = waitbar(0,'Iteration Failed...i=1');pause(2);delete(hwtb);    
    end
   end
catch 
    hwtb = waitbar(0,'Process Failed ...');pause(2);delete(hwtb);
end
end

function plotout_inv(~,~,ax2,x,y,matrix,unitt,tit)
drawnow;mapviewer(x,y,matrix,unitt,tit,ax2)
end

function g=FW_PRKR3D(z,z0,r0,nx,ny,nx0,ny0,n,k)
%%% forward gravity calculation
hs=-2*pi*20/3.*r0.*exp(-abs(k).*z0);
tongF=0;
for m=1:n;
     tongF=tongF+((-abs(k)).^(m-1))./(factorial(m)).*fft2(z.^m);
end;
Fg=hs.*tongF;
g0=(ifft2(Fg));
g0(1,1)=0;
g=real(g0);
g=g(ny/2+1:ny/2+ny0,nx/2+1:nx/2+nx0);
end

function gcalc=FORWARD_func(data,nx0,ny0,dx,dy,r0,z0,n)
%%% FORWARD PROCEDURE
z=data-z0; % extract mean depth
[ze,nxe,nye]=paddData(nx0,ny0,z); % enlarge data
nxm=2*nxe; nym=2*nye;
k=getfreqs(nxm,nym,dx,dy); % obtain wavenumbers
%%% run forward gravity calculation
gcalc=FW_PRKR3D(ze,z0,r0,nxe,nye,nx0,ny0,n,k); % calculate gravity
end

function [matrix,nx,ny]=paddData(nx,ny,matrix)
% PADDING DATA
matrix(1,nx+floor(nx/2))=0;
matrix(ny+floor(ny/2),1)=0;
matrix=rot90(rot90(matrix));
matrix(1,nx+2*floor(nx/2))=0;
matrix(ny+2*floor(ny/2),1)=0;
matrix=rot90(rot90(matrix));
if (mod(nx,2)~=0) nx=nx-1; matrix(:,end)=[]; end
if (mod(ny,2)~=0) ny=ny-1; matrix(end,:)=[]; end
end

function [k]=getfreqs(nxm,nym,dx,dy) 
% WAVENUMBERS
dkx= 2.*pi./((nxm-1).*dx);
dky= 2.*pi./((nym-1).*dy);
nyqx= (nxm/2)+1;
nyqy= (nym/2)+1; 
[kx,ky]=meshgrid([(0:nxm/2) (nyqx+1:nxm)-(nxm+1)].*dkx,...
    [((0:nym/2)) (nyqy+1:nym)-(nym+1)].*dky);
k= sqrt(bsxfun(@plus,kx.^2,ky.^2));
k(1,1)=0.00000001;
end


function mapviewer(x,y,matrix,unitt,tit,ax)
% plot of selected map

drawnow;set(gcf,'CurrentAxes',ax);refresh
contourf(x,y,matrix,18);shading flat;rotate3d off;axis equal;axis tight
set(gca,'FontSize',10,'FontWeight','bold')
h=colorbar('eastoutside');
title(h,unitt,'FontWeight','bold');
set(h,'FontSize',10,'FontWeight','bold')
xlabel('X (km)');ylabel('Y (km)');title(tit)
xlim([min(min(x)) max(max(x))])
ylim([min(min(y)) max(max(y))])
box on
end

function [matrix,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=gridform(k);
%%grid data loader
fidc=fopen(k);header= fread(fidc,4,'*char' )';fclose(fidc);
c1=strcmp(header,'DSAA');c2=strcmp(header,'DSRB');
sumc=sum([c1 c2]);
if sumc>0
switch c1
    case 1
[matrix,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd6txt(k); %format surfer6 text
    case 0
[matrix,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd7bin(k); %surfer7 binary       
end
else
matrix=0;x=0;y=0;nx=0;ny=0;xmin=0;xmax=0;ymin=0;ymax=0;dx=0;dy=0;     
end
end

function [T,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd6txt(k)
%%%%%reader, Surfer 6 text grid(*.grd)
surfergrd=fopen(k,'r'); % Open *.grid file
dsaa=fgetl(surfergrd);  % Header
% Get the map dimension [NX: East NY: North];
datasize=str2num(fgetl(surfergrd)); nx=datasize(1); ny=datasize(2);
% Map limits: xmin, xmax, ymin ymax
xcoor=str2num(fgetl(surfergrd)); xmin=xcoor(1); xmax=xcoor(2);
ycoor=str2num(fgetl(surfergrd)); ymin=ycoor(1); ymax=ycoor(2);
% check intervals in x and y direction 
dx=(xmax-xmin)/(nx-1);dx=abs(dx);
dy=(ymax-ymin)/(ny-1);dy=abs(dy);
% data limits
anom=str2num(fgetl(surfergrd)); t0min=anom(1); t0max=anom(2);
% data matrix 
[T,numb] = fscanf(surfergrd, '%f', [nx,ny]);
T=T'; % Traspose matrix
fclose(surfergrd);
% map coordinate matrix
[x,y]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);
end

function [T,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy] = lodgrd7bin(filename)
%reader, Surfer 7 Binary grid
fid= fopen(filename);
fread(fid,4,'*char' )';
fread(fid,1,'uint32');fread(fid,1,'uint32');
fread(fid,4,'*char' )';fread(fid,1,'uint32');
ny= fread(fid,1,'uint32'); nx= fread(fid,1,'uint32');
xmin= fread(fid,1,'double'); ymin= fread(fid,1,'double');
dx= fread(fid,1,'double'); dy= fread(fid,1,'double');
fread(fid,1,'double');fread(fid,1,'double');
fread(fid,1,'double');
parm= fread(fid,1,'double');
fread(fid,4,'*char' )';
nn= fread(fid,1,'uint32');
if ny*nx ~= nn/8 ; error('error') ;end
T= nan(nx,ny);
T(1:end) = fread(fid,numel(T),'double');
T=T';
fclose(fid);
T(T==parm) = nan;
xv = xmin + (0:nx-1)*dx;
yv = ymin + (0:ny-1)*dy;
[x,y]=meshgrid(xv,yv);
xmax=xv(end);
ymax=yv(end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% check and controls
function errcode=error_checker(data,dx,dy,nx,ny,modeP)
%%% error message displayer
errcode=0;
if dx~=dy;errcode=1;STR='Message [error]: equal spaced grid dx=dy is required !';end
if any(isnan(data(:)));errcode=2;STR='Message [error]: Blanked grid not supported !';end
if nx==0;errcode=3;STR='Message [error]: File format not supported !';end
if ny==0;errcode=3;STR='Message [error]: File format not supported !';end
end


function grdout(matrix,xmin,xmax,ymin,ymax,namefile)
%%%%%%%%%%% Function for output of a GRID
%Get grid dimensions
aux=size(matrix);
nx=aux(2);ny=aux(1);
grdfile=fopen(namefile,'w');                % Open file
fprintf(grdfile,'%c','DSAA');               % Header code
fprintf(grdfile,'\n %i %i',[nx ny]);        % Grid size
fprintf(grdfile,'\n %f %f',[xmin xmax]);    % X limits
fprintf(grdfile,'\n %f %f',[ymin ymax]);    % Y limits
fprintf(grdfile,'\n %f %f',[min(min(matrix)) max(max(matrix))]); % Z limits
fprintf(grdfile,'\n');
for jj=1:ny                                 % Write matrix
    for ii=1:nx
        fprintf(grdfile,'%g %c',matrix(jj,ii),' ');
    end
    fprintf(grdfile,'\n');
end
fclose(grdfile);
end

function [pwr,wn]= raps_data(T,dx)
% get radially averaged spectrum vs frequency 
[ny,nx]=size(T); 
nrow=2*floor(ny/2); ncol=2*floor(nx/2); 
maxdim=max([nrow ncol]);
np2=2^nextpow2(maxdim);  
rowdiff=round((np2-nrow)/2); 
coldiff=round((np2-ncol)/2);
T=T(1:nrow,1:ncol);
T=T-mean(T(:)); 
wf=tukeywin(nrow,.05)*tukeywin(ncol,.05)'; %truncation window 5% from edges
dw=T.*wf;
TT=zeros(np2); 
TT(rowdiff+1:rowdiff+nrow,coldiff+1:coldiff+ncol)=dw;  
spec =(abs(fftshift(fft2(TT))));
spec(spec<1.0e-13)=1.0e-13;
spec=log(spec);
% create cartesian coordinate system (unit) and
% shifting the origin to zero frequency.
[xo,yo]=meshgrid((1:np2)-np2/2,(1:np2)-np2/2); 
[~,L]=cart2pol(xo,yo);L=round(L);% calculate distances to the center
halfspec=(np2/2)-1; %considering the half of the spectrum
wn=(1:halfspec)/(np2*dx); % freq axis
for i=1:halfspec
pwr(i)=mean(spec(L==i));%find data of equal distances and get mean  
end
end


function [zcalc,rmsset]=INVERSION_func(data,nx0,ny0,dx,dy,r0,z0,WH,SH,criterio,mxit)
%%% INVERSION PROCEDURE
[g0,nxe,nye]=paddData(nx0,ny0,data);% enlarge data
nxm=2*nxe; nym=2*nye;    
k=getfreqs(nxm,nym,dx,dy); % obtain wavenumbers   
LPF=LP_filt(k,WH,SH); % filter design
%%% run inversion sheme
[zcalc,rmsset]=INV_sheme(g0,r0,z0,k,LPF,criterio,nxe,nye,nx0,ny0,mxit);
end

function LPF=LP_filt(k,WH,SH)
%%% build filter
LPF=k.*0;
[nym,nxm]=size(k);     
k=k./(2*pi);  
for j=1:nym;
   for i=1:nxm;
      if k(j,i)<WH
      LPF(j,i)=1;  
elseif k(j,i)<SH
      LPF(j,i)=0.5.*(1+cos((((2*pi)*k(j,i))-(2*pi*WH))/(2*(SH-WH))));
else
LPF(j,i)=0;
      end
    end;
end;
end

function [z,r]=INV_sheme(g0,r0,z0,k,LPF,criterio,nx,ny,nx0,ny0,mxit)
%%% the iterative inversion sheme
hwtb = waitbar(0,'Iteration Starting, please wait..');
pause(1)
Fh=-fft2(g0)./(2.*pi.*6.67.*r0.*exp(-k.*z0));
Fh=Fh.*LPF;
Fh(1,1)=0;
h=real(ifft2(Fh)); 
h_old=h;
rms=1000;
iter=0;
m=2;
while rms>criterio
    Fh=Fh-((-k).^(m-1)).*fft2(h.^m)./factorial(m);
    Fh=Fh.*LPF;
    Fh(1,1)=0;
    h=real(ifft2(Fh));
    dh=h-h_old;
    dh2=dh.^2;
    rms=sqrt(sum(sum(dh2))./(numel(dh2)));
    iter=iter+1;
    h_old=h;
    r(iter)=rms;
    waitbar(iter/mxit,hwtb,sprintf('%12.0f',iter))
    if iter==mxit;break;end
    m=m+1;
end
z=h_old+z0;
z=z(ny/2+1:ny/2+ny0,nx/2+1:nx/2+nx0);
pause(.5);delete(hwtb);
end



