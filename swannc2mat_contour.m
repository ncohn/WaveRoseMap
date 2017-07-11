%% swannc2mat_contour.m 
%% Code takes compiled netcdf files from SWAN output from Gabriel Garcia's
%% hindcast simulations (on shusin2) and generates wave time series for
%% of wave direction, height, and period for all nodes along specific
%% depth contours
%%
%% Provides file format suitable for swanmat2waveroseshape.m
%% to generate shapefiles of wave roses from time series for plotting in ArcGIS
%%
%% Note that long time spans will take a long time and alot of memory
%% to process - This is ALOT of data
%%
%% N. Cohn, Oregon State University,last updated 6 Nov. 2014
clear all, close all

%INPUTS
targetContour = 400; %depth contour to extract for wave time series
StartDate = datenum(2006,1,1);
EndDate = datenum(2011,12,31, 23, 0, 0);

%%%%%%%%CODE RUN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timenum = StartDate:1/24:EndDate;
timevec = datevec(timenum);
for t = 1:numel(timevec(:,1))
    if timevec(t,2)<10
        monthdummy = '0';
    else 
        monthdummy = '';
    end  
    if timevec(t,3)<10
        daydummy = '0';
    else
        daydummy = '';
    end  
    if timevec(t,4)<10
        hourdummy = '0';
    else
        hourdummy = '';
    end   
    filename = ['/home/shusin2/shared/swan_pnw_hindcast/GabeHindcast/',num2str(timevec(t,1)),'/SWAN', num2str(timevec(t,1)), monthdummy, num2str(timevec(t,2)), daydummy, num2str(timevec(t,3)), hourdummy, num2str(timevec(t,4)),'.nc'];  
    if t == 1
        xc = ncread(filename, 'xc');
        yc = ncread(filename, 'yc');
        depth = ncread(filename, 'depth');
        pwp = ncread(filename, 'pwp');
        swh = ncread(filename, 'swh');
        mwd = ncread(filename, 'mwd');       
        ycont = yc(1,:);
        z=depth;     
        Hs = nan(numel(ycont), numel(timevec(:,1)));
        Tp = nan(numel(ycont), numel(timevec(:,1)));
        D = nan(numel(ycont), numel(timevec(:,1)));
        for jj=1:numel(ycont)       
            xcont(jj)=linterp(z(:,jj),xc(:,jj),targetContour);
        end
    figure, pcolor(xc,yc,depth), shading interp, axis equal, axis tight, hold on, plot(xcont, ycont, 'k-')            
    else
        pwp = ncread(filename, 'pwp');
        swh = ncread(filename, 'swh');
        mwd = ncread(filename, 'mwd'); 
    end           
   for jj=1:numel(ycont)       
     Hs(jj,t)=linterp(z(:,jj),swh(:,jj),targetContour);
     Tp(jj,t)=linterp(z(:,jj),pwp(:,jj),targetContour);
     D(jj,t)=linterp(z(:,jj),mwd(:,jj),targetContour);  
   end
        clear swh pwp mwd       
end
clear t monthdummy jj hourdummy filename daydummy simLength simLengthHrs StartDate