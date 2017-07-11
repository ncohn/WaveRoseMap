%% swanmat2waveroseshape.m
%% Takes time series of wave properties to generate statistical wave roses
%% and exports those roses as shapefiles for plotting in ArcGIS
%%
%% Pay particular attention to "stride" and "scale_factor", you will need
%% to play around with those values to make a useful product depending
%% on the particualr application
%%
%% N. Cohn, Oregon State University, last updated 6 Nov. 2014

clear all, close all
% Input File Parameters                     
input_filename = 'Gabe1000m_2006-2011.mat';  %expected variables: x, y, timenum, Hs, Tp, D
output_filename = 'Gabe1000m_2006-2011_Feb';

% Geographic Data Selection/Processing Options                                                            
stride = 75; %1 for all values, >1 indicates number of cells to skip such that 1:stride:end for grid cells in analysis
num_angles = 16;            %Number of compass direction wanted for the wind rose
plotvariable = 'Hs'; %variable to make wave rose for. typically Hs (wave height) or Tp (wave period)
bins = [0 1 2 3 4 5 10]; %put largest/lowest value sufficient to capture whole range of data, dont put more than 10 bins
scale_factor = 0.01; %Will increase/decrease the size of the wave rose proportionally 
months2use = 2; %put 1:12 to use all months, can choose individual months to only consider data from those periods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%ONCE ALL PARAMETERS ABOVE ARE SET LEFT CLICK AND PRESS%%%%%%%%%%%%%%%%%% 
%%%EVALUATE CURRENT CELL, SCRIPT WILL RUN AND ALL SHAPEFILES%%%%%%%%%%%%%%%
%%%WILL BE PUT IN THE CURRENT FOLDER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(input_filename)
for station_index = 1:stride:numel(ycont)   
    lon = xcont(station_index);
    lat = ycont(station_index);
    for ii = 1:numel(months2use);
        if ii == 1
            iuse = find(timevec(:,2) == months2use(ii));
        else
            iuse = [iuse; find(timevec(:,2) == months2use(ii))];
        end
    end
    tot_count = numel(iuse);
    eval(['tempdata = ', plotvariable, ';']);
    data.variable = tempdata(station_index, iuse);
    data.dir = D(station_index, iuse);

    %Create angle bins for wave direction
    sizebin = 360/num_angles;
    for idx = 1:num_angles
            eval(['angle', num2str(idx),'.MinAngle = sizebin*idx-(3*sizebin/2);'])
            eval(['angle', num2str(idx),'.MaxAngle = sizebin*idx-(sizebin/2);'])
            if idx == 1
                        eval(['angle', num2str(idx),'.MinAngle = sizebin*idx-(3*sizebin/2)+360;'])
            end
    end

    %Create bins for wave height
    number_bins = numel(bins)-1;
    for idx = 1:number_bins
        eval(['bin', num2str(idx),'.MinH = bins(idx);'])
        eval(['bin', num2str(idx),'.MaxH = bins(idx+1);'])
        if idx ==number_bins
                eval(['bin', num2str(idx),'.MaxH = 1000;'])
        end
    end

    %Angle Calculations - Not very efficient but it works.....
    for idx =1:num_angles
         eval(['angle', num2str(idx),'.idx = find(data.dir>= angle', num2str(idx),'.MinAngle & data.dir< angle', num2str(idx),'.MaxAngle);'])
            if idx ==1
                eval(['angle', num2str(idx),'.idx1a = find(data.dir>0 & data.dir< angle', num2str(idx),'.MaxAngle);'])
                eval(['angle', num2str(idx),'.idx1b = find(data.dir> angle', num2str(idx),'.MinAngle & data.dir<= 360);'])
                angle1.idx = horzcat(angle1.idx1a, angle1.idx1b);
            end
    end
    for idx =1:num_angles
                eval(['angle', num2str(idx),'.count = numel(angle', num2str(idx),'.idx);'])
    end
    for idx =1:num_angles
        for idx2 = 1:number_bins
                eval(['angle', num2str(idx),'.height', num2str(idx2),'idx = find(data.variable(angle', num2str(idx),'.idx)>= bin',num2str(idx2),'.MinH & data.variable(angle', num2str(idx),'.idx)< bin',num2str(idx2),'.MaxH);'])
                eval(['angle', num2str(idx),'.freq', num2str(idx2),' = numel(angle', num2str(idx),'.height', num2str(idx2),'idx)*100/tot_count;'])
        end
    end
    for idx = 1:num_angles
        for idx2 = 1:number_bins
            if idx2 == 1
                 eval(['dist = angle', num2str(idx),'.freq1;'])
            end
            if idx2 == 2
                 eval(['dist = angle', num2str(idx),'.freq1 + angle', num2str(idx),'.freq2;'])
                 eval(['dist2 = angle', num2str(idx),'.freq1;'])
            end
            if idx2 == 3
                 eval(['dist = angle', num2str(idx),'.freq1 + angle', num2str(idx),'.freq2 + angle', num2str(idx),'.freq3;'])
                 eval(['dist2 = angle', num2str(idx),'.freq1 + angle', num2str(idx),'.freq2;'])
            end   
            if idx2 == 4
                 eval(['dist = angle', num2str(idx),'.freq1 + angle', num2str(idx),'.freq2 + angle', num2str(idx),'.freq3 + angle', num2str(idx),'.freq4;'])
                 eval(['dist2 = angle', num2str(idx),'.freq1 + angle', num2str(idx),'.freq2 + angle', num2str(idx),'.freq3;'])
            end    
            if idx2 == 5
                 eval(['dist = angle', num2str(idx),'.freq1 + angle', num2str(idx),'.freq2 + angle', num2str(idx),'.freq3 + angle', num2str(idx),'.freq4 +angle', num2str(idx),'.freq5;'])
                 eval(['dist2 = angle', num2str(idx),'.freq1 + angle', num2str(idx),'.freq2 + angle', num2str(idx),'.freq3 + angle', num2str(idx),'.freq4;'])
            end    
            if idx2 == 6
                 eval(['dist = angle', num2str(idx),'.freq1 + angle', num2str(idx),'.freq2 + angle', num2str(idx),'.freq3 + angle', num2str(idx),'.freq4 + angle', num2str(idx),'.freq5 + angle', num2str(idx),'.freq6;'])
                 eval(['dist2 = angle', num2str(idx),'.freq1 + angle', num2str(idx),'.freq2 + angle', num2str(idx),'.freq3 + angle', num2str(idx),'.freq4 + angle', num2str(idx),'.freq5;'])
            end    
            if idx2 == 7
                 eval(['dist = angle', num2str(idx),'.freq1 + angle', num2str(idx),'.freq2 + angle', num2str(idx),'.freq3 + angle', num2str(idx),'.freq4 + angle', num2str(idx),'.freq5 + angle', num2str(idx),'.freq6 + angle', num2str(idx),'.freq7;'])
                 eval(['dist2 = angle', num2str(idx),'.freq1 + angle', num2str(idx),'.freq2 + angle', num2str(idx),'.freq3 + angle', num2str(idx),'.freq4 + angle', num2str(idx),'.freq5 + angle', num2str(idx),'.freq6;'])
            end    
            if idx2 == 8
                 eval(['dist = angle', num2str(idx),'.freq1 + angle', num2str(idx),'.freq2 + angle', num2str(idx),'.freq3 + angle', num2str(idx),'.freq4 + angle', num2str(idx),'.freq5 + angle', num2str(idx),'.freq6 + angle', num2str(idx),'.freq7 + angle', num2str(idx),'.freq8;'])
                 eval(['dist2 = angle', num2str(idx),'.freq1 + angle', num2str(idx),'.freq2 + angle', num2str(idx),'.freq3 + angle', num2str(idx),'.freq4 + angle', num2str(idx),'.freq5 + angle', num2str(idx),'.freq6 + angle', num2str(idx),'.freq7;'])
            end
            if idx2 == 9
                 eval(['dist = angle', num2str(idx),'.freq1 + angle', num2str(idx),'.freq2 + angle', num2str(idx),'.freq3 + angle', num2str(idx),'.freq4 + angle', num2str(idx),'.freq5 + angle', num2str(idx),'.freq6 + angle', num2str(idx),'.freq7 + angle', num2str(idx),'.freq8 + angle', num2str(idx),'.freq9;'])
                 eval(['dist2 = angle', num2str(idx),'.freq1 + angle', num2str(idx),'.freq2 + angle', num2str(idx),'.freq3 + angle', num2str(idx),'.freq4 + angle', num2str(idx),'.freq5 + angle', num2str(idx),'.freq6 + angle', num2str(idx),'.freq7 + angle', num2str(idx),'.freq8;'])
            end    
            if idx2 == 10
                 eval(['dist = angle', num2str(idx),'.freq1 + angle', num2str(idx),'.freq2 + angle', num2str(idx),'.freq3 + angle', num2str(idx),'.freq4 + angle', num2str(idx),'.freq5 + angle', num2str(idx),'.freq6 + angle', num2str(idx),'.freq7 + angle', num2str(idx),'.freq8 + angle', num2str(idx),'.freq9 + angle', num2str(idx),'.freq10;'])
                 eval(['dist2 = angle', num2str(idx),'.freq1 + angle', num2str(idx),'.freq2 + angle', num2str(idx),'.freq3 + angle', num2str(idx),'.freq4 + angle', num2str(idx),'.freq5 + angle', num2str(idx),'.freq6 + angle', num2str(idx),'.freq7 + angle', num2str(idx),'.freq8 + angle', num2str(idx),'.freq9;'])
            end    
            if idx2 == 1
             eval(['angle', num2str(idx),'.p', num2str(idx2),'.x1 = 0;'])
             eval(['angle', num2str(idx),'.p', num2str(idx2),'.y1 = 0;'])
             eval(['angle', num2str(idx),'.p', num2str(idx2),'.x2 = dist*sin(deg2rad(angle', num2str(idx),'.MinAngle))*scale_factor;'])
             eval(['angle', num2str(idx),'.p', num2str(idx2),'.y2 = dist*cos(deg2rad(angle', num2str(idx),'.MinAngle))*scale_factor;'])
             eval(['angle', num2str(idx),'.p', num2str(idx2),'.x3 = dist*sin(deg2rad(angle', num2str(idx),'.MaxAngle))*scale_factor;'])
             eval(['angle', num2str(idx),'.p', num2str(idx2),'.y3 = dist*cos(deg2rad(angle', num2str(idx),'.MaxAngle))*scale_factor;'])
            end
            if idx2 > 1
             eval(['angle', num2str(idx),'.p', num2str(idx2),'.x1 = dist2*sin(deg2rad(angle', num2str(idx),'.MinAngle))*scale_factor;'])
             eval(['angle', num2str(idx),'.p', num2str(idx2),'.y1 = dist2*cos(deg2rad(angle', num2str(idx),'.MinAngle))*scale_factor;'])
             eval(['angle', num2str(idx),'.p', num2str(idx2),'.x4 = dist2*sin(deg2rad(angle', num2str(idx),'.MaxAngle))*scale_factor;'])
             eval(['angle', num2str(idx),'.p', num2str(idx2),'.y4 = dist2*cos(deg2rad(angle', num2str(idx),'.MaxAngle))*scale_factor;'])
             eval(['angle', num2str(idx),'.p', num2str(idx2),'.x2 = dist*sin(deg2rad(angle', num2str(idx),'.MinAngle))*scale_factor;'])
             eval(['angle', num2str(idx),'.p', num2str(idx2),'.y2 = dist*cos(deg2rad(angle', num2str(idx),'.MinAngle))*scale_factor;'])
             eval(['angle', num2str(idx),'.p', num2str(idx2),'.x3 = dist*sin(deg2rad(angle', num2str(idx),'.MaxAngle))*scale_factor;'])
             eval(['angle', num2str(idx),'.p', num2str(idx2),'.y3 = dist*cos(deg2rad(angle', num2str(idx),'.MaxAngle))*scale_factor;'])
            end
        clear dist dist2
        end
    end
    for idx = 1:num_angles
        for idx2 = 1:number_bins
        val = idx+(idx2-1)*(num_angles);
        if idx2 == 1
        wavedata(val).Geometry = 'Polygon';
        eval(['wavedata(val).BoundingBox(1,1) = min([angle', num2str(idx),'.p', num2str(idx2),'.x1+lon angle', num2str(idx),'.p', num2str(idx2),'.x2+lon angle', num2str(idx),'.p', num2str(idx2),'.x3+lon]);'])
        eval(['wavedata(val).BoundingBox(2,1) = max([angle', num2str(idx),'.p', num2str(idx2),'.x1+lon angle', num2str(idx),'.p', num2str(idx2),'.x2+lon angle', num2str(idx),'.p', num2str(idx2),'.x3+lon]);'])
        eval(['wavedata(val).BoundingBox(1,2) = min([angle', num2str(idx),'.p', num2str(idx2),'.y1+lat angle', num2str(idx),'.p', num2str(idx2),'.y2+lat angle', num2str(idx),'.p', num2str(idx2),'.y3+lat]);'])
        eval(['wavedata(val).BoundingBox(2,2) = max([angle', num2str(idx),'.p', num2str(idx2),'.y1+lat angle', num2str(idx),'.p', num2str(idx2),'.y2+lat angle', num2str(idx),'.p', num2str(idx2),'.y3+lat]);'])
        eval(['wavedata(val).X = [angle', num2str(idx),'.p', num2str(idx2),'.x1+lon angle', num2str(idx),'.p', num2str(idx2),'.x2+lon angle', num2str(idx),'.p', num2str(idx2),'.x3+lon];'])
        eval(['wavedata(val).Y = [angle', num2str(idx),'.p', num2str(idx2),'.y1+lat angle', num2str(idx),'.p', num2str(idx2),'.y2+lat angle', num2str(idx),'.p', num2str(idx2),'.y3+lat];'])
        eval(['wavedata(val).MinH = bin', num2str(idx2),'.MinH;'])
        eval(['wavedata(val).MaxH = bin', num2str(idx2),'.MaxH;'])
        eval(['wavedata(val).MinAngle = angle', num2str(idx),'.MinAngle;'])
        eval(['wavedata(val).MaxAngle = angle', num2str(idx),'.MaxAngle;'])
        end
        if idx2 > 1
        wavedata(val).Geometry = 'Polygon';
        eval(['wavedata(val).BoundingBox(1,1) = min([angle', num2str(idx),'.p', num2str(idx2),'.x1+lon angle', num2str(idx),'.p', num2str(idx2),'.x2+lon angle', num2str(idx),'.p', num2str(idx2),'.x3+lon angle', num2str(idx),'.p', num2str(idx2),'.x4+lon]);'])
        eval(['wavedata(val).BoundingBox(2,1) = max([angle', num2str(idx),'.p', num2str(idx2),'.x1+lon angle', num2str(idx),'.p', num2str(idx2),'.x2+lon angle', num2str(idx),'.p', num2str(idx2),'.x3+lon angle', num2str(idx),'.p', num2str(idx2),'.x4+lon]);'])
        eval(['wavedata(val).BoundingBox(1,2) = min([angle', num2str(idx),'.p', num2str(idx2),'.y1+lat angle', num2str(idx),'.p', num2str(idx2),'.y2+lat angle', num2str(idx),'.p', num2str(idx2),'.y3+lat angle', num2str(idx),'.p', num2str(idx2),'.y4+lat]);'])
        eval(['wavedata(val).BoundingBox(2,2) = max([angle', num2str(idx),'.p', num2str(idx2),'.y1+lat angle', num2str(idx),'.p', num2str(idx2),'.y2+lat angle', num2str(idx),'.p', num2str(idx2),'.y3+lat angle', num2str(idx),'.p', num2str(idx2),'.y4+lat]);'])
        eval(['wavedata(val).X = [angle', num2str(idx),'.p', num2str(idx2),'.x1+lon angle', num2str(idx),'.p', num2str(idx2),'.x2+lon angle', num2str(idx),'.p', num2str(idx2),'.x3+lon angle', num2str(idx),'.p', num2str(idx2),'.x4+lon];'])
        eval(['wavedata(val).Y = [angle', num2str(idx),'.p', num2str(idx2),'.y1+lat angle', num2str(idx),'.p', num2str(idx2),'.y2+lat angle', num2str(idx),'.p', num2str(idx2),'.y3+lat angle', num2str(idx),'.p', num2str(idx2),'.y4+lat];'])
        eval(['wavedata(val).MinH = bin', num2str(idx2),'.MinH;'])
        eval(['wavedata(val).MaxH = bin', num2str(idx2),'.MaxH;'])
        eval(['wavedata(val).MinAngle = angle', num2str(idx),'.MinAngle;'])
        eval(['wavedata(val).MaxAngle = angle', num2str(idx),'.MaxAngle;'])
        end
        end
    end
    for idx = 1:numel(wavedata)      
        wavedata(idx).X = double(wavedata(idx).X);
        wavedata(idx).Y = double(wavedata(idx).Y);      
    end
%     output_filename3 = [output_filename,'_Node', num2str(station_index),'.shp'];
%     if IndYearlyWind ==1
%         shapewrite(wavedata, output_filename3);
%     end
    eval(['wavedata', num2str(station_index),' = wavedata;'])   
    clear wavedata    
end
for station_index = 1:stride:numel(ycont)
    if station_index ==1 
            combined_data = wavedata1;
    end
    if station_index>1
        eval(['combined_data = vertcat(combined_data, wavedata',num2str(station_index),');'])
    end
end
output_filename2 = [output_filename,'.shp'];
shapewrite(combined_data, output_filename2);
   
    