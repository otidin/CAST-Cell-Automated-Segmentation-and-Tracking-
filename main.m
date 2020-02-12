
% failed experiment
%% experiment paramaters
cd /Users/onurtidin/Desktop/desktopvirtual/AllMicroscopeExperiments/20180104new

expe = experimentPara();

mainDir = expe.mainDir;
imgDir = expe.imgDir;
N = expe.numberOfFrames;

% Define some stuff, move into the right folder

movie = 1;

binsize = 1;    

outDir = [mainDir '/movie' num2str(movie) '/'];
cd(outDir)

stim=expe.stim
save stim.mat stim


addpath ../
addpath ../code
addpath ../code/Tracking
addpath ../code/regionbased_seg
set(0,'defaultlinelinewidth',2)


%% reload everything if needed (after closing matlab)

f = dir('*.mat');
for i=1:length(f)
load(f(i).name);
end

% reduce image resolution by binning
% note: binning takes the mean of the pixel, to avoid overflow in 16bis

% doDraw = 1;
% doImageBinning;

% run these commands to delete everything, if you want to reset the movie
% and start from scratch, (command-t to uncomment)

% !rm zStackedThreshCorrected/*.png
% !rm zStackedThreshSplit/*.png
% !rm zStackedThresh/*.png
% !rm zStackedYFP/*.png
% !rm links.mat
% !rm zStackedYFP_unbinned/*.png
% !rm segPara.mat
% !rm segParaFrames.mat
% !rm divMatrixFinal.mat
% !rm peakMatrixFinal.mat

% combine stacks, denoise, ...

Nz = expe.numberOfColors; %number of images to combine

deNoise = {'none','BM3D','median','localNorm'}; %denoise algo on each stack
deNoise = deNoise{1};

medianSize = 3;
weightsSegmentation = [1 0 1]; %weights for summing the different channels
compressionQuantile = 0.999;   %signal above this quantile will be cut off, set to 1 to disable
gaussianFilterSize = 20;       %typycal length of the background

temporalBinning = 1;

doDraw = 0;

mkdirIfNotExist('zStackedYFP');

combineImages;

%just display the combined images

clf;
for k=1:N
    
    disp(100*k/N)                 
    a = imread(['zStackedYFP/' num2str(k) '.png']);
    clf
    imagesc(a);
    %caxis([0 250])        
    pause(0.08); drawnow;          
end

% test segmentation on a few frames 
       
guiSeg

for i=1:500
system('open -a MATLAB_R2016b.app')
pause(30)
end

%% Do the segmentation
inputFolder = 'zStackedYFP/';

load segPara.mat
load segParaFrames.mat

doDraw = 0;
segMethod = 1;
%
for frame=1:N
    
    disp(100*frame/N);
    
    para = [];        
    para.minSize  = interp1(segParaFrames,[segPara(:).minSize  ],frame);
    para.maxSize  = interp1(segParaFrames,[segPara(:).maxSize  ],frame);
    para.thresh   = interp1(segParaFrames,[segPara(:).thresh   ],frame);
    para.openSize = interp1(segParaFrames,[segPara(:).openSize ],frame);
    para.ratioThresh = interp1(segParaFrames,[segPara(:).ratioThresh ],frame);
                
    if(segMethod ==1)
        
        [filters openFilter] = generateFilters(para,0);
        segmentImage(frame,para,inputFolder,filters,openFilter,doDraw);
    else
    
        load ../segmentationCoeff.mat
        seg = segmentLDA(inputFolder,segPara,segParaFrames,frame,coeff,N1,N2);
        clf
        imagesc(seg)
        drawnow    
    end
end


for i=1:500
system('open -a MATLAB_R2016b.app')
pause(30)
end

%%
% do the measures

threshFolder = 'zStackedThresh/';
doMeasures(N,threshFolder,expe);


%
% try to split some merged cells
saveFolder = 'zStackedThreshSplit/';
mkdirIfNotExist(saveFolder)

doDraw = 0;

threshold = 0.5; %low value -> split everything

splitMergeCells;

threshFolder = saveFolder;
doMeasures(N,threshFolder,expe);

% look at difference between splited and original

doDraw = 1;
if doDraw 
    for k=1:N
        
         a = imread(['zStackedThresh/' num2str(k) '.png']);
         b = imread(['zStackedThreshSplit/' num2str(k) '.png']);

         a = double(a);
         b = double(b);

         imagesc(a+b);
         pause(0.03)
    end
end

for i=1:500
system('open -a MATLAB_R2016b.app')
pause(30)
end

%%
% correct segmentation by hand

linksGui

for i=1:500
system('open -a MATLAB_R2016b.app')
pause(30)
end
%% redo the measures with corrected images

threshFolder = 'zStackedThreshCorrected/';
nObj = doMeasures(N,threshFolder,expe);

clf;
plot(nObj)
ylabel('number of objects')

for i=1:500
system('open -a MATLAB_R2016b.app')
pause(30)
end

%% tweak tracking parameters if needed
% main parameters are frame_displacement and split_cost

edit('get_struct.m');

%% load all measures & do the final Tracking

NToTrack = N;

doLinksOnly = 0;

clc
Me = loadMeasures(N);
[tracks, signal, traj, ind, divisions,divPerframe,trajX,trajY] = doTracking(NToTrack, Me,doLinksOnly);

clf; imagesc( signal(:,:,1) ); colormap jet

for i=1:500
system('open -a MATLAB_R2016b.app')
pause(30)
end

%% plot Tracking 

pauseTime = 0.01;
longTracesOnly = 0;

plotTracking;

for i=1:500
system('open -a MATLAB_R2016b.app')
pause(30)
end
%% revert temporal binning: go back to full framerate (cannot undo)

doDraw = 1;
revertTemporalBinning2


%% select good traces, based on length
delete lengthThresh.mat
clf
lengthOfTrace = zeros(size(ind,1),1);
lengthOfGaps  = zeros(size(ind,1),1);

indAnnotation = zeros(size(ind));
% delete lengthThresh.mat
if( exist('lengthThresh.mat','file') )
    load lengthThresh.mat;
else
    lengthThresh = 0.40; %note: to change lengthThresh value you first need to delete the file if it exists: !rm lengthThresh.mat
end

% for i=1:size(ind,1)
   for i=1:120
       
    %look for continous traces
    indAnnotation(i,:) = markTrace(ind(i,:));

    lengthOfGaps(i) = sum( indAnnotation(i,:) == -3);
    lengthOfTrace(i) = sum( indAnnotation(i,:) == 1);
    
   end

longTraces = find( (lengthOfGaps < 10) .* (lengthOfTrace/N>lengthThresh) );
% longTraces = longTraces(1:10)
% clf
% figure
A = signal(longTraces,:);
[tmp,ia,ic] = unique(A,'rows');

longTraces = sort( longTraces(ia) );

imagesc(signal(longTraces,:,1))
% colormap default
% longTraces=good;
% length(longTraces)
colormap jet
% longTraces=longTraces(1:50);
save lengthThresh.mat lengthThresh
save longTraces.mat longTraces
% longTraces = [3 10 11 23 24 26];


%%
clf
% datato2=refinedMean(longTraces(:),:,2)./refinedMeancyto(longTraces(:),:,2);
% datato2=signal(longTraces(:),:,1);
datato2=signal(longTraces(:),:,2);

[data1,ave1]=calcave(datato2, expe.dt);

tstart=1;
tend=N;
annotvec=stim;

doAnnot=1;


hold all
plot(expe.dt*[tstart:tend], data1);
plot(expe.dt*[tstart:tend], ave1,'k','LineWidth',5);
Annot
%% build area, sum of signal, peak and div matrices
close all
minTimeBetweenPeaks = 20;
peakMethod ='diff';
doDraw =0;
makePeakAndDivMatrices
annotvec=stim;
appendvec={};
save stim.mat stim
%%
%%
% % guiTraces
% clear longTraces
% load longTraces
% load goodTraces
% find(goodTraces==1)
% longTraces=longTraces(find(goodTraces==1));
system('open ./figures')

% %% refine area and signal around each cell, and do images for guiTraces
k=1;
% longTraces=longTraces(2:6)
doDrawBkg = 0;  %display there area where the background is measured 
doDraw = 0;     %display area refinement result

bgkSize = 20;    %size around the cell where the background is not quantified
superSampling = 1; %increase the resolution of the image (must be an integer)

NIteration = 20; % Number of iteration of the area refinement algorithm, increase when using temporal binning
dilateSizeAfterRefine = 0.1; %if >=1 dilate a bit the area after the refinement
dilateSizeAfterRefineext = [3 5 8]; %if >=1 dilate a bit the area after the refinement
dilateSizeAfterRefineint = [3];

s = 60; %size of the window around the cells

mkdirIfNotExist('snapShots')

a = imread(['zStackedYFP/' num2str(1) '.png']);

N1=size(a,1);
N2=size(a,2);






%%














%% refine area and signal around each cell, and do images for guiTraces

doDrawBkg = 0;  %display there area where the background is measured 
doDraw = 0;     %display area refinement result

bgkSize = 15;    %size around the cell where the background is not quantified
superSampling = 1; %increase the resolution of the image (must be an integer)

NIteration 
N= 20; % Number of iteration of the area refinement algorithm, increase when using temporal binning
dilateSizeAfterRefine = 0.1; %if >=1 dilate a bit the area after the refinement
dilateSizeAfterRefineext = 3; %if >=1 dilate a bit the area after the refinement

s = 60; %size of the window around the cells

mkdirIfNotExist('snapShots')

a = imread(['zStackedYFP/' num2str(1) '.png']);

N1=size(a,1);
N2=size(a,2);

touchBorder = zeros(size(areaMatrix)); %is equale to one if the cell is touching the border of the image
refinedArea = zeros(size(areaMatrix));

touchBorder = zeros(size(areaMatrix)); %is equale to one if the cell is touching the border of the image
refinedAreaext = zeros(size(areaMatrix));

refinedMean = zeros([size(areaMatrix) size(signal,3)]); %mean of the pixel values 
refinedMeanext = zeros([size(areaMatrix) size(signal,3)]); %mean of the pixel values 
refinedSum = zeros([size(areaMatrix) size(signal,3)]);  %sum of the pixel values 
refinedSumext = zeros([size(areaMatrix) size(signal,3)]);  %sum of the pixel values 
refinedStd  = zeros([size(areaMatrix) size(signal,3)]); %standard deviation of the pixel values 
refinedStdext  = zeros([size(areaMatrix) size(signal,3)]); %standard deviation of the pixel values 

bkg  = zeros([size(areaMatrix) size(signal,3)]);

refineAreaAndSignal

refinedMeanextb=refinedMeanext;
refinedSumextb=refinedSumext;
refinedMeanb=refinedMean;
refinedSumb=refinedSum;

save refinedMeanext.mat refinedMeanext
save refinedSumext.mat refinedSumext
save refinedAreaext.mat refinedAreaext

save refinedMeanextb.mat refinedMeanextb
save refinedSumextb.mat refinedSumextb
save refinedMeanb.mat refinedMeanb
save refinedSumb.mat refinedSumb

%% make small images around each cell for the trace tool 
% use if you don't want to do the refine area thing above

mkdirIfNotExist('snapShots')

doDraw = 0;
useFullSizeImages = 0;
inputFolder = 'zStackedYFP/';
threshFoler = {'zStackedThreshCorrected','zStackedThreshCorrectedRefined'};
threshFoler = threshFoler{2};
  
NToTrack = N;

colorIndex = 2;
s = 40;

traces = 1:length(longTraces);

makeImagesForTraceTool

makeImagesForTraceTool2


%% delete peakMatrixFinal and divMatrixFinal (reset guiTraces)

%!rm divMatrixFinal.mat
%!rm peakMatrixFinal.mat

%% Traces tool, correct divs and peaks

guiTraces

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAKE SOME PLOTS %%%%%%%%%%%%%%%%%%%%%%%

%% which traces are good ?
clear longTraces
load longTraces
load goodTraces
find(goodTraces==1)
longTraces=longTraces(find(goodTraces==1));

%% conditional


% discard=[684 601 481 86];
% outind=find(ismember(longTraces,discard));
% outind=[96 ]
.
mn=[];
for i=1:length(longTraces)
sel = ind(longTraces(i),:) > 0;
mn(i)=mean(refinedMean(longTraces(i),sel,1));
end
length(find(mn<500))
% longTraces(find(mn<500)) = [];
%%
mn=[];
for i=1:length(longTraces)
sel = ind(longTraces(i),:) > 0;
mn(i)=mean(refinedMean(longTraces(i),sel,2));
end
length(find(mn<250))
find(mn<250)
% longTraces(find(mn<400)) = [];
longTraces(find(mn<250))
% longTraces
% size(longTraces)

%%

 nlucfluc=[31 193 224];
% longTraces=longTraces(setdiff(1:length(longTraces),nlucfluc));
 [~,discard]=ismember(nlucfluc,longTraces);
longTraces(discard) = [];
%  longTraces=longTraces(discard);
longTraces
size(longTraces)

%%
%additional discard
adddis=[283];
[~,adddiscard]=ismember(adddis,longTraces);
longTraces(adddiscard) = [];

longTraces
size(longTraces)
save longTracesFinal.mat longTraces 
%% define signals
refinedMeantemp=[];
refinedMeanexttemp=[];
refinedSumtemp=[];
refinedSumexttemp=[];
refinedMeancyto=[];
refinedSumcyto=[];
refinedMeancytotemp=[];
minvalue=50;

refinedMeantemp=refinedMeanb(:,:,:)-bkg(:,:,:);
refinedMeantemp(:,:,2) =refinedMeantemp(:,:,2)-(min(min(refinedMeantemp(longTraces(:),:,2)'))-minvalue);
minmean=min(refinedMeantemp(longTraces(:),:,2)')

refinedMeanexttemp=refinedMeanextb(:,:,:)-bkg(:,:,:);
minmeanext=min(refinedMeanexttemp(longTraces(:),:,2)')
 refinedMeanexttemp(:,:,2) =refinedMeanexttemp(:,:,2)-(min(min(refinedMeanexttemp(longTraces(:),:,2)'))-minvalue);
% 
refinedSumtemp(:,:,1)=refinedSumb(:,:,1)-bkg(:,:,1).*refinedArea(:,:);
refinedSumtemp(:,:,2)=refinedSumb(:,:,2)-bkg(:,:,2).*refinedArea(:,:);
refinedSumtemp(:,:,2) =refinedSumtemp(:,:,2)-(min(min(refinedMeantemp(longTraces(:),:,2)'))-minvalue).*refinedArea(:,:);
% 
refinedSumexttemp(:,:,1)=refinedSumextb(:,:,1)-bkg(:,:,1).*refinedAreaext(:,:);

refinedSumexttemp(:,:,2)=refinedSumextb(:,:,2)-bkg(:,:,2).*refinedAreaext(:,:);
refinedSumexttemp(:,:,2) =refinedSumexttemp(:,:,2)-(min(min(refinedMeanexttemp(longTraces(:),:,2)'))-minvalue).*refinedAreaext(:,:);



figure
plot(refinedSumexttemp(longTraces(1),:,2))
figure
plot(refinedSumtemp(longTraces(1),:,2))

% figure
% plot(refinedSumextb(longTraces(1),:,2))
% figure
% plot(refinedSumb(longTraces(1),:,2))


%assign them 

refinedMean=refinedMeantemp;
refinedMeanext=refinedMeanexttemp;
refinedSum=refinedSumtemp;
refinedSumext=refinedSumexttemp;


refinedMeancyto(:,:,1)= (refinedSumext(:,:,1)-refinedSum(:,:,1))./(refinedAreaext(:,:)-refinedArea(:,:));
refinedMeancyto(:,:,2)= (refinedSumext(:,:,2)-refinedSum(:,:,2))./(refinedAreaext(:,:)-refinedArea(:,:));
refinedSumcyto= refinedSumext(:,:,:)-refinedSum(:,:,:);

minmeancyto=min(refinedMeancyto(longTraces(:),:,2)')
refinedMeancytotemp(:,:,2) =refinedMeancyto(:,:,2)-(min(min(refinedMeancyto(longTraces(:),:,2)'))-minvalue);
minmeancyto=min(refinedMeancytotemp(longTraces(:),:,2)')

refinedMeancyto=refinedMeancytotemp;


%% save after cleaning  objects
 
% annotvec=[23 300];
% stim=[23 300];
stim=14;
annotvec=[14];
appendvec={};
save stim.mat stim
save longTracesFinal.mat longTraces
save refinedMean.mat refinedMean
save refinedSum.mat refinedSum
save refinedMeanext.mat refinedMeanext
save refinedSumext.mat refinedSumext
save refinedAreaext.mat refinedAreaext
save refinedMeancyto.mat refinedMeancyto

%%
%calculate normalized signals


%% plot mean signal of trace i with std

i=1

clfh
sel = ind(longTraces(i),:) > 0;
col = {'r','g'};
for k=1:min(2,expe.numberOfColors)
    errorbar(expe.t(    sel),refinedMean(longTraces(i),sel,k),refinedStd(longTraces(i),sel,k),col{k});
end

% text( 20, 20, 'My Nice Title', 'FontSize', 34)
% title( 10,10,'My Nice Title', 'FontSize', 34)


xlabel({'hello';
    'there sdakdnaskdnaskndsa'; 'd askdnasdjnas dasdasdasdsadasds dsadas das das d sad as da d s sd asd as da d as'})



%% cal

%normalize meansignal background subtracted
close all
datatonorm=zeros(length(longTraces),expe.numberOfFrames);
for i=1:length(longTraces)
colorIndex = 1;
sel = ind(longTraces(i),:) > 0;
datatonorm(i,:)=refinedMean(longTraces(i),:,colorIndex)-bkg(longTraces(i),:,colorIndex);
end
[ave, signnormflucMean]= normsignal(datatonorm, expe.dt);
% clear datatonorm

close all
datatonorm=zeros(length(longTraces),expe.numberOfFrames);
for i=1:length(longTraces)
colorIndex = 2;
sel = ind(longTraces(i),:) > 0;
datatonorm(i,:)=refinedMean(longTraces(i),:,colorIndex)-bkg(longTraces(i),:,colorIndex);
end
[ave, signnormNlucMean]= normsignal(datatonorm, expe.dt);
% plotallgraphs(signnorm, expe.dt)

normMean={signnormflucMean,signnormNlucMean}
%%
%normalize sumsignal background subtracted
close all
datatonorm=zeros(length(longTraces),expe.numberOfFrames);
for i=1:length(longTraces)
colorIndex = 1;
sel = ind(longTraces(i),:) > 0;
datatonorm(i,:)=refinedSum(longTraces(i),:,colorIndex)-bkg(longTraces(i),:,colorIndex).*refinedArea(longTraces(i),:);
end
[ave, signnormflucSum]= normsignal(datatonorm, expe.dt);

close all
datatonorm=zeros(length(longTraces),expe.numberOfFrames);
for i=1:length(longTraces)
colorIndex = 2;
sel = ind(longTraces(i),:) > 0;
datatonorm(i,:)=refinedSum(longTraces(i),:,colorIndex)-bkg(longTraces(i),:,colorIndex).*refinedArea(longTraces(i),:);
end
[ave, signnormNlucSum]= normsignal(datatonorm, expe.dt);

normSum={signnormflucSum,signnormNlucSum}
% normSum{1:44,1:200,1}=signnormflucSum;
% normSum{:,2}=signnormNlucSum;
% clear signnormNlucSum signnormflucSum 

%% plot mean signal, trace i, and export
close all
% figure
doExplain=0;
doAnnot=0;
doOpen=1;
totCell=length(longTraces);
totCellstart=1;
 ps=12; % plots/page
 np=length(longTraces); % # plots
pagecount=0;
for i=totCellstart:totCell
    remx=rem(i-1,ps)+1;
    % ha = tight_subplot(2,2,[.01 .03],[.1 .01],[.01 .01])
    colorIndex = 2;
    hold all
    if remx==1
        pagecount=pagecount+1;
        figure;      
end
subplot(4,3,remx)
sel = ind(longTraces(i),:) > 0;
% sel=sel(1:100);
plot(expe.t(sel),refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex),'r');
%  plot(expe.t(sel),refinedSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel),'r');
% plot(expe.t(sel),normMean{colorIndex}(i,sel),'r');
 
xlim([0 10]);
%  ylim([0 5000])
% w=6000; h=4000;
% set(figure(1),'Position',[300 600 w h])
% set(gcf, 'Position', get(0,'Screensize'));
% plot(expe.t(sel),bkg(longTraces(i),sel,colorIndex),'r--');
xlabel('time [h]');
title(expe.colorNames{colorIndex});
Annot

if doExplain
str = ['Nanoluc-luc dox inducible cell line experiment-Nuclear mean Nanoluciferase signal-xy limits fixed', char(10)...
        'During continuous flow-plots showing the stable region(no dox induction-only in the end)'];
    annotation('textbox', [0 0.9 1 0.1], ...
    'String', str, ...
    'FontSize', 70, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
end

setFonts
paperSize(100,100)
mkdirIfNotExist('figures')

fname =['figures/Mean_' expe.colorNames{colorIndex} '_page' n2s(pagecount) '.pdf'];

 if doOpen
        print('-dpdf',fname)
%        system(['open ' fname])
 end

doAppend=0;
appendName
end

% close all

system(['open ' [pwd '/figures']])

%% plot mean signal, trace i, and export
clfh
% figure
doExplain=1;
doAnnot=1;
doOpen=0;
totCell=length(longTraces);
totCellstart=1;
 ps=12; % plots/page
 np=length(longTraces); % # plots
pagecount=0;
for i=totCellstart:totCell
    remx=rem(i-1,ps)+1;
    % ha = tight_subplot(2,2,[.01 .03],[.1 .01],[.01 .01])
    colorIndex = 1;
    hold all
    if remx==1
        pagecount=pagecount+1;
        figure;      
end
subplot(4,3,remx)
sel = ind(longTraces(i),:) > 0;
plot(expe.t(sel),signnorm(longTraces(i),sel,colorIndex),'r');
% w=6000; h=4000;
% set(figure(1),'Position',[300 600 w h])
% set(gcf, 'Position', get(0,'Screensize'));
% plot(expe.t(sel),bkg(longTraces(i),sel,colorIndex),'r--');
xlabel('time [h]');
title(expe.colorNames{colorIndex});
Annot

if doExplain
str = ['NLS555 experiment to check with TGFb Stimulation. Two time']
annotation('textbox', [0 0.9 1 0.1], ...
    'String', str, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
end

setFonts
paperSize(100,100)
mkdirIfNotExist('figures')

fname =['figures/mean_' expe.colorNames{colorIndex} '_page' n2s(pagecount) '.pdf'];

if doOpen
    print('-dpdf',fname)
    % system(['open ' fname])
end


doAppend=1;
appendName
end

close all
%% %% plot all traces in the same plot new version
clfh
close all
doExplain=0;
doAnnot=1;
doOpen=1;


datato=zeros(length(longTraces),expe.numberOfFrames);
datato2=zeros(length(longTraces),expe.numberOfFrames);


datato=refinedMean(longTraces(:),:,2)./refinedMeancyto(longTraces(:),:,2);
datato2=refinedMean(longTraces(:),:,1);

[data1,ave1]=calcave(datato, expe.dt);
[data2,ave2]=calcave(datato2, expe.dt);

tstart=1;
tend=N;
data1=data1(:,tstart:tend);
data2=data2(:,tstart:tend);
ave1=ave1(tstart:tend);
ave2=ave2(tstart:tend);

hold all
plot(expe.dt*[tstart:tend], data1);
plot(expe.dt*[tstart:tend], ave1,'k','LineWidth',5);
setFonts
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean.pdf'];
print('-dpdf',fname)

figure
hold all
plot(expe.dt*[tstart:tend], data2);
plot(expe.dt*[tstart:tend], ave2,'k','LineWidth',7);
setFonts
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean2.pdf'];
print('-dpdf',fname)

figure;hold on;
y=data1;
x=expe.dt*[tstart:tend];
% H(1) = shadedErrorBar(x, y, {@nanmean, @(x) 2*nanstd(x)  }, '-r', 0);
% H(2) = shadedErrorBar(x, y, {@nanmean, @(x) 1*nanstd(x)  }, '-m', 0);
H(3) = shadedErrorBar(x, y, {@nanmean, @(x) 0.5*nanstd(x)}, {'-b', 'LineWidth', 5}, 0);
legend([H(3).mainLine, H.patch], ...
    '\mu','0.5\sigma', ...
    'Location', 'Northeast');
setFonts
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean_shaded.pdf'];
print('-dpdf',fname)

figure;hold on;
y=data2;
x=expe.dt*[tstart:tend];
% H(1) = shadedErrorBar(x, y, {@nanmean, @(x) 2*nanstd(x)  }, '-r', 0);
% H(2) = shadedErrorBar(x, y, {@nanmean, @(x) 1*nanstd(x)  }, '-m', 0);
H(3) = shadedErrorBar(x, y, {@nanmean, @(x) 0.5*nanstd(x)}, {'-b', 'LineWidth', 5}, 0);
legend([H(3).mainLine, H.patch], ...
    '\mu','0.5\sigma', ...
    'Location', 'Northeast');
setFonts
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean2_shaded.pdf'];
print('-dpdf',fname)


figure
 [hAx,hLine1,hLine2] = plotyy(expe.dt*[tstart:tend],ave1,expe.dt*[tstart:tend],ave2);
% plot(expe.t(sel),corRatio(sel).*(refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex)));
% title(['After Correction Nluc and fluc '])
ylabel(hAx(1),'nuc/cell ratio') % left y-axis
ylabel(hAx(2),'nuc/cyto ratio') % right y-axis
set(hAx,{'ycolor'},{'k';'b'})
set(hLine1,'linewidth',7)
% set(get(hAX(2),'color'),'black');
set(hLine2,'linewidth',7)
set(hLine1,'color','black');
set(hLine2,'color','blue');
xlabel('Time (h)')
setFonts
% ylim(hAx(2), [0 2*max(corRatio(sel).*(refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2)))]);
% xlim(hAx(1), [0.5 12]);
%  xlim(hAx(2), [0.5 12]);
doAnnot=1;
Annot
% xlim([1 12])
% ylim([0 1000])
setFonts
paperSize(50,20)
fname =['figures/alltraces_mean_tworatios.pdf'];
print('-dpdf',fname)

% 
figure
x=expe.dt*[tstart:tend];
shadedErrorBaryy(x,nanmean(data1,1),0.5*nanstd(data1,1),'r',x,nanmean(data2,1),0.5*nanstd(data2,1),'b')
setFonts
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean_shadedtwo.pdf'];
print('-dpdf',fname)
%% %% plot all traces in the same plot
clfh
close all
doExplain=0;
doAnnot=1;
doOpen=1;

datato=zeros(length(longTraces),expe.numberOfFrames);

% totCell=length(longTraces);
% totCellstart=1;
for i=1:length(longTraces)
% subplot(nRow,round(totCell/nRow),i)
colorIndex = 1;
hold all
% clfh
sel = ind(longTraces(i),:) > 0;
% sel(1:10)=0;
%  plot(expe.t(:),refinedMean(longTraces(i),:,colorIndex));
% plot(expe.t(sel),0.1*divMatrix(longTraces(i),sel),'r');
datato(i,:)=refinedMean(longTraces(i),:,colorIndex)-bkg(longTraces(i),:,colorIndex);
% datato(i,:)=refinedSum(longTraces(i),:,colorIndex)-bkg(longTraces(i),:,colorIndex).*refinedArea(longTraces(i),:);

% datato(i,:)=normSum{colorIndex}(i,:);

xlabel('time [h]');
if doExplain
str = ['time [h]' char(10) 'Nanoluciferase-SMAD4 Dox inducible cell line experiment. Nuclear mean Nluc signal(SMAD4)'];
xlabel(str);
end
title(expe.colorNames{colorIndex});
end
% plot(expe.t(sel),mean(refinedMean(longTraces(1:9),sel,colorIndex)),'k','LineWidth',10);

% startt=1;
% endt=296;
ave=plotallgraphs(datato, expe.dt)

datato=zeros(length(longTraces),expe.numberOfFrames);

setFonts
% ylim([0 1000])
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean_' expe.colorNames{colorIndex} '.pdf'];
if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end


figure
% totCell=length(longTraces);
% totCellstart=1;
for i=1:length(longTraces)
% subplot(nRow,round(totCell/nRow),i)
colorIndex = 2;
hold all
sel = ind(longTraces(i),:) > 0;
datato(i,:)=refinedMean(longTraces(i),:,colorIndex)-bkg(longTraces(i),:,colorIndex);
% datato(i,:)=refinedSum(longTraces(i),:,colorIndex)-bkg(longTraces(i),:,colorIndex).*refinedArea(longTraces(i),:);

% datato(i,:)=normSum{colorIndex}(i,:);

xlabel('time [h]');
if doExplain
str = ['time [h]' char(10) 'Nanoluciferase-SMAD4 Dox inducible cell line experiment. Nuclear mean Nluc signal(SMAD4)'];
xlabel(str);
end
title(expe.colorNames{colorIndex});
end
% plot(expe.t(sel),mean(refinedMean(longTraces(1:9),sel,colorIndex)),'k','LineWidth',10);

ave2=plotallgraphs(datato, expe.dt)
%  ylim([0 1500])


% plot(expe.t(:),ave./ave2)
Annot

setFonts

paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean_' expe.colorNames{colorIndex} '.pdf'];
if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end

figure
 [hAx,hLine1,hLine2] = plotyy(expe.dt*[1:N],ave,expe.dt*[1:N],ave2);
% plot(expe.t(sel),corRatio(sel).*(refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex)));
% title(['After Correction Nluc and fluc '])
ylabel(hAx(1),'ctgf activity') % left y-axis
ylabel(hAx(2),'Nuclear SMAD4 concentration') % right y-axis
set(hAx,{'ycolor'},{'k';'b'})
set(hLine1,'linewidth',7)
% set(get(hAX(2),'color'),'black');
set(hLine2,'linewidth',7)
set(hLine1,'color','black');
set(hLine2,'color','blue');
xlabel('Time (h)')
% ylim(hAx(1), [0 2*max(ave1)]);
%  ylim(hAx(2), [0 2*max(ave2)]);
setFonts
% ylim(hAx(2), [0 2*max(corRatio(sel).*(refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2)))]);
% xlim(hAx(1), [0.5 12]);
%  xlim(hAx(2), [0.5 12]);
doAnnot=1;
Annot
% xlim([1 12])
setFonts
% ylim([0 1000])
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean.pdf'];
if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end


% plot(expe.t(:),ave./ave2)


doAppend=0;
appendName

system(['open ' [pwd '/figures']])

%% %% plot all traces in the same plot
clfh
close all
doExplain=0;
doAnnot=1;
doOpen=1;

datato=zeros(length(longTraces),expe.numberOfFrames);

% totCell=length(longTraces);
% totCellstart=1;
for i=1:length(longTraces)
% subplot(nRow,round(totCell/nRow),i)
colorIndex = 1;
hold all
% clfh
sel = ind(longTraces(i),:) > 0;
% sel(1:10)=0;
%  plot(expe.t(:),refinedMean(longTraces(i),:,colorIndex));
% plot(expe.t(sel),0.1*divMatrix(longTraces(i),sel),'r');
% datato(i,:)=(refinedMean(longTraces(i),:,2)-bkg(longTraces(i),:,2))./((refinedSumext(longTraces(i),:,2)-refinedSum(longTraces(i),:,2))./(refinedAreaext(longTraces(i),:)-refinedArea(longTraces(i),:))-bkg(longTraces(i),:,2));
% datato(i,:)=(refinedMean(longTraces(i),:,2)-bkg(longTraces(i),:,2))./((refinedSumext(longTraces(i),:,2)-refinedSum(longTraces(i),:,2))./(refinedAreaext(longTraces(i),:)-refinedArea(longTraces(i),:))-bkg(longTraces(i),:,2));
datato(i,:)=refinedMean(longTraces(i),:,colorIndex)-bkg(longTraces(i),:,colorIndex);

% datato(i,:)=refinedSum(longTraces(i),:,colorIndex)-bkg(longTraces(i),:,colorIndex).*refinedArea(longTraces(i),:);


% datato(i,:)=normSum{colorIndex}(i,:);

xlabel('time [h]');
if doExplain
str = ['time [h]' char(10) 'Nanoluciferase-SMAD4 Dox inducible cell line experiment. Nuclear mean Nluc signal(SMAD4)'];
xlabel(str);
end
title(expe.colorNames{colorIndex});
end
% plot(expe.t(sel),mean(refinedMean(longTraces(1:9),sel,colorIndex)),'k','LineWidth',10);

% startt=1;
% endt=296;
ave=plotallgraphs(datato, expe.dt)

datato=zeros(length(longTraces),expe.numberOfFrames);

setFonts
% ylim([0 1000])
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean_' expe.colorNames{colorIndex} '.pdf'];
if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end


figure
% totCell=length(longTraces);
% totCellstart=1;
for i=1:length(longTraces)
% subplot(nRow,round(totCell/nRow),i)
colorIndex = 2;
hold all
sel = ind(longTraces(i),:) > 0;
 datato(i,:)=medfilt1( (refinedMean(longTraces(i),:,2)-bkg(longTraces(i),:,2))./((refinedSumext(longTraces(i),:,2)-refinedSum(longTraces(i),:,2))./(refinedAreaext(longTraces(i),:)-refinedArea(longTraces(i),:))-bkg(longTraces(i),:,2)),10);
%  (refinedMean(longTraces(i),:,2)-bkg(longTraces(i),:,2))./((refinedSumext(longTraces(i),:,2)-refinedSum(longTraces(i),:,2))./(refinedAreaext(longTraces(i),:)-refinedArea(longTraces(i),:))-bkg(longTraces(i),:,2))

% datato(i,:)=refinedMean(longTraces(i),:,colorIndex)-bkg(longTraces(i),:,colorIndex);
% datato(i,:)=refinedSum(longTraces(i),:,colorIndex)-bkg(longTraces(i),:,colorIndex).*refinedArea(longTraces(i),:);

% datato(i,:)=normSum{colorIndex}(i,:);

xlabel('time [h]');
if doExplain
str = ['time [h]' char(10) 'Nanoluciferase-SMAD4 Dox inducible cell line experiment. Nuclear mean Nluc signal(SMAD4)'];
xlabel(str);
end
title(expe.colorNames{colorIndex});
end
% plot(expe.t(sel),mean(refinedMean(longTraces(1:9),sel,colorIndex)),'k','LineWidth',10);

ave2=plotallgraphs(datato, expe.dt)
%  ylim([0 1500])


% plot(expe.t(:),ave./ave2)
Annot

setFonts

paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean_' expe.colorNames{colorIndex} '.pdf'];
if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end

figure
 [hAx,hLine1,hLine2] = plotyy(expe.dt*[1:N],ave,expe.dt*[1:N],ave2);
% plot(expe.t(sel),corRatio(sel).*(refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex)));
% title(['After Correction Nluc and fluc '])
ylabel(hAx(1),'ctgf activity') % left y-axis
ylabel(hAx(2),'Nuclear SMAD4 concentration') % right y-axis
set(hAx,{'ycolor'},{'k';'b'})
set(hLine1,'linewidth',7)
% set(get(hAX(2),'color'),'black');
set(hLine2,'linewidth',7)
set(hLine1,'color','black');
set(hLine2,'color','blue');
xlabel('Time (h)')
% ylim(hAx(1), [0 2*max(ave1)]);
 ylim(hAx(2), [1 1.5]);
setFonts
% ylim(hAx(2), [0 2*max(corRatio(sel).*(refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2)))]);
% xlim(hAx(1), [0.5 12]);
%  xlim(hAx(2), [0.5 12]);
doAnnot=1;
Annot
% xlim([1 12])
setFonts
% ylim([0 1000])
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean_ratio.pdf'];
if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end


% plot(expe.t(:),ave./ave2)


doAppend=0;
appendName

system(['open ' [pwd '/figures']])

%% %% plot all traces in the same plot
clfh
close all
doExplain=0;
doAnnot=1;
doOpen=1;

datato=zeros(length(longTraces),expe.numberOfFrames);
datato2=zeros(length(longTraces),expe.numberOfFrames);

% totCell=length(longTraces);
% totCellstart=1;
for i=1:length(longTraces)
% subplot(nRow,round(totCell/nRow),i)
colorIndex = 1;
hold all
% clfh
sel = ind(longTraces(i),:) > 0;
% sel(1:10)=0;
%  plot(expe.t(:),refinedMean(longTraces(i),:,colorIndex));
% plot(expe.t(sel),0.1*divMatrix(longTraces(i),sel),'r');
%  datato2(i,:)=(refinedMean(longTraces(i),:,2)-bkg(longTraces(i),:,2))./((refinedSumext(longTraces(i),:,2)-refinedSum(longTraces(i),:,2))./(refinedAreaext(longTraces(i),:)-refinedArea(longTraces(i),:))-bkg(longTraces(i),:,2));
%  datato2(i,:)=(refinedMean(longTraces(i),:,2)-bkg(longTraces(i),:,2))./((refinedSumext(longTraces(i),:,2)-refinedSum(longTraces(i),:,2))./(refinedAreaext(longTraces(i),:)-refinedArea(longTraces(i),:))-bkg(longTraces(i),:,2));
datato2(i,:)=(refinedMean(longTraces(i),:,2)-bkg(longTraces(i),:,2))./(refinedMeanext(longTraces(i),:,2)-bkg(longTraces(i),:,2));
%  datato(i,:)=medfilt1 (refinedSumext(longTraces(i),:,2)-refinedSum(longTraces(i),:,2),10);
% datato2(i,:)=refinedSum(longTraces(i),:,2)-bkg(longTraces(i),:,2).*refinedArea(longTraces(i),:);
 datato(i,:)=medfilt1(datato2(i,:),15);
% datato(i,:)=normSum{colorIndex}(i,:);

xlabel('time [h]');
if doExplain
str = ['time [h]' char(10) 'Nanoluciferase-SMAD4 Dox inducible cell line experiment. Nuclear mean Nluc signal(SMAD4)'];
xlabel(str);
end
title(expe.colorNames{colorIndex});
end
% plot(expe.t(sel),mean(refinedMean(longTraces(1:9),sel,colorIndex)),'k','LineWidth',10);

% startt=1;
% endt=296;
ave=plotallgraphs(datato, expe.dt)

datato=zeros(length(longTraces),expe.numberOfFrames);
datato2=zeros(length(longTraces),expe.numberOfFrames);

setFonts
% ylim([0 1000])
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean_' expe.colorNames{colorIndex} '.pdf'];
if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end


figure
% totCell=length(longTraces);
% totCellstart=1;
for i=1:length(longTraces)
% subplot(nRow,round(totCell/nRow),i)
colorIndex = 2;
hold all
sel = ind(longTraces(i),:) > 0;
 datato(i,:)=medfilt1( (refinedMean(longTraces(i),:,2)-bkg(longTraces(i),:,2))./((refinedSumext(longTraces(i),:,2)-refinedSum(longTraces(i),:,2))./(refinedAreaext(longTraces(i),:)-refinedArea(longTraces(i),:))-bkg(longTraces(i),:,2)),15);
%  (refinedMean(longTraces(i),:,2)-bkg(longTraces(i),:,2))./((refinedSumext(longTraces(i),:,2)-refinedSum(longTraces(i),:,2))./(refinedAreaext(longTraces(i),:)-refinedArea(longTraces(i),:))-bkg(longTraces(i),:,2))

% datato(i,:)=refinedMean(longTraces(i),:,colorIndex)-bkg(longTraces(i),:,colorIndex);
% datato(i,:)=refinedSum(longTraces(i),:,colorIndex)-bkg(longTraces(i),:,colorIndex).*refinedArea(longTraces(i),:);

% datato(i,:)=normSum{colorIndex}(i,:);

xlabel('time [h]');
if doExplain
str = ['time [h]' char(10) 'Nanoluciferase-SMAD4 Dox inducible cell line experiment. Nuclear mean Nluc signal(SMAD4)'];
xlabel(str);
end
title(expe.colorNames{colorIndex});
end
% plot(expe.t(sel),mean(refinedMean(longTraces(1:9),sel,colorIndex)),'k','LineWidth',10);

ave2=plotallgraphs(datato, expe.dt)
%  ylim([0 1500])


% plot(expe.t(:),ave./ave2)
Annot

setFonts

paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean_' expe.colorNames{colorIndex} '.pdf'];
if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end

figure
 [hAx,hLine1,hLine2] = plotyy(expe.dt*[1:N],ave,expe.dt*[1:N],ave2);
% plot(expe.t(sel),corRatio(sel).*(refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex)));
% title(['After Correction Nluc and fluc '])
ylabel(hAx(1),'nuc/cell ratio') % left y-axis
ylabel(hAx(2),'nuc/cyto ratio') % right y-axis
set(hAx,{'ycolor'},{'k';'b'})
set(hLine1,'linewidth',7)
% set(get(hAX(2),'color'),'black');
set(hLine2,'linewidth',7)
set(hLine1,'color','black');
set(hLine2,'color','blue');
xlabel('Time (h)')
 xlim(hAx(1), [2 30]);
  xlim(hAx(2), [2 30]);
 ylim(hAx(1), [1 1.5]);
  ylim(hAx(2), [1 1.5]);
setFonts
% ylim(hAx(2), [0 2*max(corRatio(sel).*(refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2)))]);
% xlim(hAx(1), [0.5 12]);
%  xlim(hAx(2), [0.5 12]);
doAnnot=1;
Annot
% xlim([1 12])
setFonts
% ylim([0 1000])
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean_tworatios.pdf'];
if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end


% plot(expe.t(:),ave./ave2)


doAppend=0;
appendName

system(['open ' [pwd '/figures']])
%% plot area and division trace i 

i=2

clfh
sel = ind(longTraces(i),:) > 0;
plot(expe.t(sel),imnorm(refinedArea(longTraces(i),sel,1)),'k');
plot(expe.t(sel),0.5*divMatrix(longTraces(i),sel),'r');

legend('area','divisions')

xlabel('time [h]')
%% plot mean signal, trace i, and export

close all
doExplain=0;
doAnnot=1;
doOpen=1;

totCell=length(longTraces);
totCellstart=1;
 ps=12; % plots/page
 np=length(longTraces); % # plots
pagecount=0;
mn=[];
for i=totCellstart:totCell
remx=rem(i-1,ps)+1;    

colorIndex = 2;
hold all
sel = ind(longTraces(i),:) > 0;
    if remx==1
        pagecount=pagecount+1;
        figure;      
    end
% sel([1:27])=0;
subplot(4,3,remx)
% [hAx,hLine1,hLine2] = plotyy(expe.t(sel), normMean{1}(i,sel),expe.t(sel),normMean{2}(i,sel));
% [hAx,hLine1,hLine2] = plotyy(expe.t(sel),refinedSum(longTraces(i),sel,1)-bkg(longTraces(i),sel,1),expe.t(sel),refinedSum(longTraces(i),sel,2)-bkg(longTraces(i),sel,2));
% [hAx,hLine1,hLine2] = plotyy(expe.t(sel),refinedMean(longTraces(i),sel,1)-bkg(longTraces(i),sel,1),expe.t(sel),refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2));
[hAx,hLine1,hLine2] = plotyy(expe.t(sel),refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,1),expe.t(sel),(refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex))./(refinedMeanext(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex)));

% (refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex))./(refinedMeanext(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex))

% mn(i)=mean(refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2));
ylabel(hAx(1),'fluc channel') % left y-axis
ylabel(hAx(2),'Nluc channel') % right y-axis
ylim(hAx(1), [0 1000]);
ylim(hAx(2), [0 3]);
%  xlim(hAx(1), [0 10]);
%   xlim(hAx(2), [0 10]);
title(n2s(longTraces(i)))

% xlim(hAx(1), expe.dt*[28 N]);
% xlim(hAx(2), expe.dt*[28 N]);
%  xlim(hAx(2), [0 7]);
% plot(expe.t(sel),bkg(longTraces(i),sel,colorIndex),'r--');
% colorIndex = 1;
% plot(expe.t(sel),refinedMean(longTraces(i),sel,colorIndex));
%  str = ['Luciferase-Nanoluciferase Dox inducible cell line experiment' char(10) 'Nluc-fluc traces together' char(10) 'Before dox induction signal levels are low for most of cells. Only a few cells showing the effect of addition of prosubstaret at t=5h.']

xlabel('time [h]');
Annot


% axis([min(expe.t(sel)) max(expe.t(sel)) 0 1.1*max(refinedMean(longTraces(i),sel,colorIndex))])

setFonts
paperSize(100,100)
mkdirIfNotExist('figures')

if doExplain
str = ['Nanoluc-luc dox inducible cell line experiment-Nuclear mean Nanoluciferase-luciferase signal-xy limits fixed', char(10)...
        'During continuous flow-plots showing the stable region(no dox induction-only in the end)'];
    annotation('textbox', [0 0.9 1 0.1], ...
    'String', str, ...
    'FontSize', 15, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
end


% if doExplain
%  str = ['Luciferase-Nanoluciferase Dox inducible cell line experiment' char(10) 'Nluc-fluc traces together' char(10) 'Before dox induction signal levels are low for most of cells. Only a few cells showing the effect of addition of prosubstaret at t=5h.']
% [ax1,h1]=suplabel('str');
% set(h1,'FontSize',25) 
% end

fname =['figures/sum_sum_signal_tracesNlucandratio' expe.colorNames{colorIndex} '_page' n2s(pagecount) '.pdf'];

if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end



doAppend=1;
appendName
end
% close all
system(['open ' [pwd '/figures']])

% figure
% hist(mn,50)

%% plot mean signal, trace i, and export

close all
doExplain=0;
doAnnot=1;
doOpen=1;

totCell=length(longTraces);
totCellstart=1;
 ps=12; % plots/page
 np=length(longTraces); % # plots
pagecount=0;
mn=[];
for i=totCellstart:totCell
remx=rem(i-1,ps)+1;    

colorIndex = 2;
hold all
sel = ind(longTraces(i),:) > 0;
    if remx==1
        pagecount=pagecount+1;
        figure;      
    end
% sel([1:27])=0;
subplot(4,3,remx)
% [hAx,hLine1,hLine2] = plotyy(expe.t(sel), normMean{1}(i,sel),expe.t(sel),normMean{2}(i,sel));
% [hAx,hLine1,hLine2] = plotyy(expe.t(sel),refinedSum(longTraces(i),sel,1)-bkg(longTraces(i),sel,1),expe.t(sel),refinedSum(longTraces(i),sel,2)-bkg(longTraces(i),sel,2));
% [hAx,hLine1,hLine2] = plotyy(expe.t(sel),refinedMean(longTraces(i),sel,1)-bkg(longTraces(i),sel,1),expe.t(sel),refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2));
% [hAx,hLine1,hLine2] = plotyy(expe.t(sel),refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2),expe.t(sel),(refinedMean(longTraces(i),sel,2))./((refinedSumext(longTraces(i),sel,2)-refinedSum(longTraces(i),sel,2))./(refinedAreaext(longTraces(i),sel)-refinedArea(longTraces(i),sel))));
% [hAx,hLine1,hLine2] = plotyy(expe.t(sel),refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2),expe.t(sel),(refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2))./((refinedSumext(longTraces(i),sel,2)-refinedSum(longTraces(i),sel,2))./(refinedAreaext(longTraces(i),sel)-refinedArea(longTraces(i),sel))-bkg(longTraces(i),sel,2)));
[hAx,hLine1,hLine2] = plotyy(expe.t(sel),refinedMean(longTraces(i),sel,1)-bkg(longTraces(i),sel,1),expe.t(sel),(refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2))./((refinedSumext(longTraces(i),sel,2)-refinedSum(longTraces(i),sel,2))./(refinedAreaext(longTraces(i),sel)-refinedArea(longTraces(i),sel))-bkg(longTraces(i),sel,2)));

% [hAx,hLine1,hLine2] = plotyy(expe.t(sel),refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,1),expe.t(sel),refinedSum(longTraces(i),sel,2)./(refinedSumext(longTraces(i),sel,2)-refinedSum(longTraces(i),sel,2)));

% (refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex))./(refinedMeanext(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex))

% mn(i)=mean(refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2));
ylabel(hAx(1),'fluc channel') % left y-axis
ylabel(hAx(2),'Nluc channel') % right y-axis
%  ylim(hAx(1), [0 2000]);
 ylim(hAx(2), [0 4]);
%  xlim(hAx(1), [0 10]);
%   xlim(hAx(2), [0 10]);
title(n2s(longTraces(i)))

% xlim(hAx(1), expe.dt*[28 N]);
% xlim(hAx(2), expe.dt*[28 N]);
%  xlim(hAx(2), [0 7]);
% plot(expe.t(sel),bkg(longTraces(i),sel,colorIndex),'r--');
% colorIndex = 1;
% plot(expe.t(sel),refinedMean(longTraces(i),sel,colorIndex));
%  str = ['Luciferase-Nanoluciferase Dox inducible cell line experiment' char(10) 'Nluc-fluc traces together' char(10) 'Before dox induction signal levels are low for most of cells. Only a few cells showing the effect of addition of prosubstaret at t=5h.']

xlabel('time [h]');
Annot


% axis([min(expe.t(sel)) max(expe.t(sel)) 0 1.1*max(refinedMean(longTraces(i),sel,colorIndex))])

setFonts
paperSize(100,100)
mkdirIfNotExist('figures')

if doExplain
str = ['Nanoluc-luc dox inducible cell line experiment-Nuclear mean Nanoluciferase-luciferase signal-xy limits fixed', char(10)...
        'During continuous flow-plots showing the stable region(no dox induction-only in the end)'];
    annotation('textbox', [0 0.9 1 0.1], ...
    'String', str, ...
    'FontSize', 15, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
end


% if doExplain
%  str = ['Luciferase-Nanoluciferase Dox inducible cell line experiment' char(10) 'Nluc-fluc traces together' char(10) 'Before dox induction signal levels are low for most of cells. Only a few cells showing the effect of addition of prosubstaret at t=5h.']
% [ax1,h1]=suplabel('str');
% set(h1,'FontSize',25) 
% end

fname =['figures/sum_sum_signal_tracescytonucleusratio' expe.colorNames{colorIndex} '_page' n2s(pagecount) '.pdf'];

if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end



doAppend=1;
appendName
end
% close all
system(['open ' [pwd '/figures']])

% figure
% hist(mn,50)

%% plot mean signal, trace i, and export

close all
doExplain=0;
doAnnot=1;
doOpen=1;

totCell=length(longTraces);
totCellstart=1;
 ps=12; % plots/page
 np=length(longTraces); % # plots
pagecount=0;
mn=[];
for i=totCellstart:totCell
remx=rem(i-1,ps)+1;    

colorIndex = 2;
hold all
sel = ind(longTraces(i),:) > 0;
    if remx==1
        pagecount=pagecount+1;
        figure;      
    end
% sel([1:27])=0;
subplot(4,3,remx)
% [hAx,hLine1,hLine2] = plotyy(expe.t(sel), normMean{1}(i,sel),expe.t(sel),normMean{2}(i,sel));
% [hAx,hLine1,hLine2] = plotyy(expe.t(sel),refinedSum(longTraces(i),sel,1)-bkg(longTraces(i),sel,1),expe.t(sel),refinedSum(longTraces(i),sel,2)-bkg(longTraces(i),sel,2));
% [hAx,hLine1,hLine2] = plotyy(expe.t(sel),refinedMean(longTraces(i),sel,1)-bkg(longTraces(i),sel,1),expe.t(sel),refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2));
% [hAx,hLine1,hLine2] = plotyy(expe.t(sel),refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2),expe.t(sel),(refinedMean(longTraces(i),sel,2))./((refinedSumext(longTraces(i),sel,2)-refinedSum(longTraces(i),sel,2))./(refinedAreaext(longTraces(i),sel)-refinedArea(longTraces(i),sel))));
% [hAx,hLine1,hLine2] = plotyy(expe.t(sel),refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2),expe.t(sel),(refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2))./((refinedSumext(longTraces(i),sel,2)-refinedSum(longTraces(i),sel,2))./(refinedAreaext(longTraces(i),sel)-refinedArea(longTraces(i),sel))-bkg(longTraces(i),sel,2)));
[hAx,hLine1,hLine2] = plotyy(expe.t(sel),refinedMean(longTraces(i),sel,1)-bkg(longTraces(i),sel,1),expe.t(sel),medfilt1((refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2))./((refinedSumext(longTraces(i),sel,2)-refinedSum(longTraces(i),sel,2))./(refinedAreaext(longTraces(i),sel)-refinedArea(longTraces(i),sel))-bkg(longTraces(i),sel,2)),20));

% [hAx,hLine1,hLine2] = plotyy(expe.t(sel),refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,1),expe.t(sel),refinedSum(longTraces(i),sel,2)./(refinedSumext(longTraces(i),sel,2)-refinedSum(longTraces(i),sel,2)));

% (refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex))./(refinedMeanext(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex))

% mn(i)=mean(refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2));
ylabel(hAx(1),'fluc channel') % left y-axis
ylabel(hAx(2),'Nluc channel') % right y-axis
%  ylim(hAx(1), [0 2000]);
%  ylim(hAx(2), [0 4]);
 xlim(hAx(1), [0 30]);
xlim(hAx(2), [0 30]);
title(n2s(longTraces(i)))

% xlim(hAx(1), expe.dt*[28 N]);
% xlim(hAx(2), expe.dt*[28 N]);
%  xlim(hAx(2), [0 7]);
% plot(expe.t(sel),bkg(longTraces(i),sel,colorIndex),'r--');
% colorIndex = 1;
% plot(expe.t(sel),refinedMean(longTraces(i),sel,colorIndex));
%  str = ['Luciferase-Nanoluciferase Dox inducible cell line experiment' char(10) 'Nluc-fluc traces together' char(10) 'Before dox induction signal levels are low for most of cells. Only a few cells showing the effect of addition of prosubstaret at t=5h.']

xlabel('time [h]');
Annot


% axis([min(expe.t(sel)) max(expe.t(sel)) 0 1.1*max(refinedMean(longTraces(i),sel,colorIndex))])

setFonts
paperSize(100,100)
mkdirIfNotExist('figures')

if doExplain
str = ['Nanoluc-luc dox inducible cell line experiment-Nuclear mean Nanoluciferase-luciferase signal-xy limits fixed', char(10)...
        'During continuous flow-plots showing the stable region(no dox induction-only in the end)'];
    annotation('textbox', [0 0.9 1 0.1], ...
    'String', str, ...
    'FontSize', 15, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
end


% if doExplain
%  str = ['Luciferase-Nanoluciferase Dox inducible cell line experiment' char(10) 'Nluc-fluc traces together' char(10) 'Before dox induction signal levels are low for most of cells. Only a few cells showing the effect of addition of prosubstaret at t=5h.']
% [ax1,h1]=suplabel('str');
% set(h1,'FontSize',25) 
% end

fname =['figures/sum_sum_signal_tracescytonucleusratio' expe.colorNames{colorIndex} '_page' n2s(pagecount) '.pdf'];

if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end



doAppend=1;
appendName
end
% close all
system(['open ' [pwd '/figures']])

% figure
% hist(mn,50)
%% 
%calculate Nluc/fluc correction
close all
doExplain=0;
doAnnot=0;
doOpen=1;
totCell=length(longTraces);
totCellstart=1;
 ps=16; % plots/page
 np=length(longTraces); % # plots
pagecount=0;
corRatios=zeros(length(longTraces),expe.numberOfFrames);
for i=totCellstart:totCell
%     for i=1:1
sel = ind(longTraces(i),:) > 0;
% corRatios(i,sel)=(refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2))./(refinedMean(longTraces(i),sel,1)-bkg(longTraces(i),sel,1));
corRatios(i,sel)=(refinedSum(longTraces(i),sel,2)-bkg(longTraces(i),sel,2).*refinedArea(longTraces(i),sel))./(refinedSum(longTraces(i),sel,1)-bkg(longTraces(i),sel,1).*refinedArea(longTraces(i),sel));
end

for j=1:size(corRatios,2)
    tot=0;
    count=0;
    for k=1:size(corRatios,1)
        
    if any(corRatios(k,j)) 
     
%     else 
      tot=tot+corRatios(k,j);  
      count=count+1;
    end
    
    end
    corRatio(j)=tot/count;
    
end
corRatio(:)=1;
corRatios(:,:)=1;
plot(corRatio)
%% plot everything about a cell
close all
% figure
doExplain=0;
doAnnot=0;
doOpen=1;
totCell=length(longTraces);
% totCell=30;
totCellstart=1;
d1=4;
d2=4;
 ps=16; % plots/page
 np=length(longTraces); % # plots
pagecount=0;
for i=totCellstart:totCell
%     for i=45:45
sel = ind(longTraces(i),:) > 0; 

colorIndex=1        
subplot(d1,d2,1)
plot(expe.t(sel),refinedSum(longTraces(i),sel,colorIndex));
title(['sum ' expe.colorNames{colorIndex}])
xlabel('time [h]');
% xlim([0 7])
% ylim([0 2*max(refinedSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel))])
hold all
% plot(expe.t(sel),signalSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel));
% plot(expe.t(sel),signalSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*areaMatrix(longTraces(i),sel));
% plot(expe.t(sel),refinedSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel)-0.05*(refinedSum(longTraces(i),sel,2)-bkg(longTraces(i),sel,2).*refinedArea(longTraces(i),sel)));
plot(expe.t(sel),bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel));
plot(expe.t(sel),refinedSum(longTraces(i),sel,colorIndex));
% plot(expe.t(sel),bkg(longTraces(i),sel,colorIndex));
legend('sig','bkg','raw','Location','northeast')
title(n2s(longTraces(i)))


hold off


subplot(d1,d2,2)
plot(expe.t(sel),refinedMean(longTraces(i),sel,colorIndex));
title(['mean ' expe.colorNames{colorIndex}])
xlabel('time [h]');
% xlim([0 7])
% ylim([0 2*max(refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex))])
hold all
% plot(expe.t(sel),bkg(longTraces(i)));
% plot(expe.t(sel),signal(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex));


hold off

% subplot(d1,d2,3)
% plot(expe.t(sel),normSum{colorIndex}(i,sel));
% title(['sum normalized ' expe.colorNames{colorIndex}])
% xlabel('time [h]');
% % xlim([0 7])
% ylim([0 2*max(normSum{colorIndex}(i,sel))])

colorIndex=2;
subplot(d1,d2,3)
 plot(expe.t(sel),bkg(longTraces(i),sel,1),expe.t(sel),bkg(longTraces(i),sel,2));

% plot(expe.t(sel),(refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex))./(refinedMeanext(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex)));
%  plot(expe.t(sel),((refinedSumext(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedAreaext(longTraces(i),sel))-(refinedSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel)))./(refinedSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel)));
%  plot(expe.t(sel),((refinedSumext(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedAreaext(longTraces(i),sel))-(refinedSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel)))./(refinedSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel)));

 title(['RATIO SUM ' expe.colorNames{colorIndex}])
xlabel('time [h]');
% xlim([0 10])
%  ylim([0 3])

% subplot(d1,d2,4)
% plot(expe.t(sel),normMean{colorIndex}(i,sel));
% title(['mean normalized ' expe.colorNames{colorIndex}])
% xlabel('time [h]');
% % xlim([0 7])
% ylim([0 2*max(normMean{colorIndex}(i,sel))])

colorIndex=2;
 subplot(d1,d2,4)
% sel=25:107;
 plot(expe.t(sel),(refinedMean(longTraces(i),sel,colorIndex))./(refinedMeanext(longTraces(i),sel,colorIndex)));

% plot(expe.t(sel),(refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex))./(refinedMeanext(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex)));
%  plot(expe.t(sel),((refinedSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel))-(refinedSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel)))./(refinedSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel)));
title(['RATIO SUM ' expe.colorNames{colorIndex}])
xlabel('time [h]');
xlim([0 10])
%  ylim([0 3])


colorIndex=2        
subplot(d1,d2,5)
plot(expe.t(sel),refinedSum(longTraces(i),sel,colorIndex));
title(['sum ' expe.colorNames{colorIndex}])
xlabel('time [h]');
% xlim([0 7])
% ylim([0 2*max(refinedSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel))])
hold all
plot(expe.t(sel),bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel));
plot(expe.t(sel),refinedSum(longTraces(i),sel,colorIndex));
plot(expe.t(sel),bkg(longTraces(i),sel,colorIndex));
plot(expe.t(sel),signalSum(longTraces(i),sel,colorIndex));
hold off

subplot(d1,d2,6)
plot(expe.t(sel),refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex));
title(['mean ' expe.colorNames{colorIndex}])
xlabel('time [h]');
% xlim([0 7])
% ylim([0 2*max(refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex))+1]);
hold all
% plot(expe.t(sel),ave)
hold off

subplot(d1,d2,7)
plot(expe.t(sel),normSum{colorIndex}(i,sel));
title(['sum normalized ' expe.colorNames{colorIndex}])
xlabel('time [h]');
% xlim([0 7])
% ylim([0 2*max(normSum{colorIndex}(i,sel))])


subplot(d1,d2,8)
plot(expe.t(sel),normMean{colorIndex}(i,sel));
title(['mean normalized' expe.colorNames{colorIndex}])
xlabel('time [h]');
% xlim([0 7])
% ylim([0 2*max(normMean{colorIndex}(i,sel))])

subplot(d1,d2,9)
plot(expe.t(sel),refinedArea(longTraces(i),sel));
hold all
plot(expe.t(sel),refinedAreaext(longTraces(i),sel));
hold all
plot(expe.t(sel),areaMatrix(longTraces(i),sel));
legend('ref','refext','orig')

title(['Area'])
xlabel('time [h]');
hold all
% plot(expe.t(sel),areaMatrix(longTraces(i),sel));
% xlim([0 7])
% ylim([0 2*max(refinedArea(longTraces(i),sel))])
hold off

subplot(d1,d2,10)
plot(expe.t(sel),normSum{2}(i,sel)./normSum{1}(i,sel));
title(['normalized Sum Nluc/fluc Ratio'])
xlabel('time [h]');
hold all
% plot(expe.t(sel),areaMatrix(longTraces(i),sel));
% xlim([0 7])
% ylim([0 2*max(normSum{2}(i,sel)./normSum{1}(i,sel))])
hold off

subplot(d1,d2,11)
plot(expe.t(sel),normMean{2}(i,sel)./normMean{1}(i,sel));
title(['normalized Mean Nluc/fluc Ratio'])
xlabel('time [h]');
hold all
% plot(expe.t(sel),areaMatrix(longTraces(i),sel));
% xlim([0 7])
% ylim([0 2*max(normMean{2}(i,sel)./normMean{1}(i,sel))])
hold off



subplot(d1,d2,12)
plot(expe.t(sel),(refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2))./(refinedMean(longTraces(i),sel,1)-bkg(longTraces(i),sel,1)));
title(['Mean Nluc/fluc Ratio'])
xlabel('time [h]');
hold all
% plot(expe.t(sel),areaMatrix(longTraces(i),sel));
% xlim([0 7])
% ylim([0 2*max((refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2))./(refinedMean(longTraces(i),sel,1)-bkg(longTraces(i),sel,1)))])
hold off


subplot(d1,d2,13)
plot(expe.t(sel),(refinedSum(longTraces(i),sel,2)-bkg(longTraces(i),sel,2).*refinedArea(longTraces(i),sel))./(refinedSum(longTraces(i),sel,1)-bkg(longTraces(i),sel,1).*refinedArea(longTraces(i),sel)));
title(['Sum Nluc/fluc Ratio'])
xlabel('time [h]');
hold all
% plot(expe.t(sel),areaMatrix(longTraces(i),sel));
xlim([0 35])
% ylim([0 2*max((refinedSum(longTraces(i),sel,2)-bkg(longTraces(i),sel,2).*refinedArea(longTraces(i),sel))./(refinedSum(longTraces(i),sel,1)-bkg(longTraces(i),sel,1).*refinedArea(longTraces(i),sel)))])
% ylim([0 2*max((refinedSum(longTraces(i),sel,1)-bkg(longTraces(i),sel,1).*refinedArea(longTraces(i),sel))./(refinedSum(longTraces(i),sel,2)-bkg(longTraces(i),sel,2).*refinedArea(longTraces(i),sel)))])
hold off


% subplot(5,5,14)
% hold all
% plot(expe.t(sel),(refinedSum(longTraces(i),sel,2)-bkg(longTraces(i),sel,2).*refinedArea(longTraces(i),sel))./(refinedSum(longTraces(i),sel,1)-bkg(longTraces(i),sel,1).*refinedArea(longTraces(i),sel)));
% title(['Sum Ratio All' expe.colorNames{colorIndex}])
% xlabel('time [h]');
% % plot(expe.t(sel),areaMatrix(longTraces(i),sel));
% % xlim([0 7])
% ylim([0 2*max((refinedSum(longTraces(i),sel,2)-bkg(longTraces(i),sel,2).*refinedArea(longTraces(i),sel))./(refinedSum(longTraces(i),sel,1)-bkg(longTraces(i),sel,1).*refinedArea(longTraces(i),sel)))])
% hold off


subplot(d1,d2,14)
plot(expe.t(:),corRatio(:))
% plot(expe.t(sel),(refinedSum(longTraces(i),sel,2)-bkg(longTraces(i),sel,2).*refinedArea(longTraces(i),sel))./(refinedSum(longTraces(i),sel,1)-bkg(longTraces(i),sel,1).*refinedArea(longTraces(i),sel)));
title(['Population Nluc/fluc Correction Ratio'])
xlabel('time [h]');
% plot(expe.t(sel),areaMatrix(longTraces(i),sel));
hold all
% xlim([0 7])
%   ylim([0 2*max(corRatio(sel))])
%   ylim([0 2*max((refinedSum(longTraces(i),sel,2)-bkg(longTraces(i),sel,2).*refinedArea(longTraces(i),sel))./(refinedSum(longTraces(i),sel,1)-bkg(longTraces(i),sel,1).*refinedArea(longTraces(i),sel)))])
 hold off

% subplot(4,4,14)
% colorIndex=2;
% plot(expe.t(sel),corRatio(sel).*(refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex)));
% title(['Mean Corrected' expe.colorNames{colorIndex}])
% xlabel('time [h]');
% hold all
% % plot(expe.t(sel),areaMatrix(longTraces(i),sel));
% xlim([0 7])
% ylim([0 2*max(corRatio(sel).*(refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex)))])
% hold off


subplot(d1,d2,15)
colorIndex=2;
[hAx,hLine1,hLine2] = plotyy(expe.t(sel),refinedMean(longTraces(i),sel,1)-bkg(longTraces(i),sel,1),expe.t(sel),corRatio(sel).*(refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2)));
% plot(expe.t(sel),corRatio(sel).*(refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex)));
title(['After Correction Nluc and fluc '])
ylabel(hAx(1),'ctgf activity') % left y-axis
ylabel(hAx(2),'nuclear SMAD4') % right y-axis
% ylim(hAx(1), [0 2*max(refinedMean(longTraces(i),sel,1)-bkg(longTraces(i),sel,1))]);
% ylim(hAx(2), [0 2*max(corRatio(sel).*(refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2)))]);
% xlim(hAx(1), [0 12]);
% xlim(hAx(2), [0 12]);

 
% subplot(4,4,14) 
%  
%  
% corrcoef(refinedMean(longTraces(i),sel,1)-bkg(longTraces(i),sel,1),corRatio(sel).*(refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2)))
%  

% colorIndex=1
% plot(expe.t(sel),normMean{colorIndex}(i,sel));
% 
% colorIndex=2
% plot(expe.t(sel),normMean{colorIndex}(i,sel));


% subplot(4,3,2)
% plot(expe.t(sel),refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex),'r');
%  plot(expe.t(sel),refinedSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel),'r');

 
%  xlim([0 7])
%  ylim([0 15000])
% w=6000; h=4000;
% set(figure(1),'Position',[300 600 w h])
% set(gcf, 'Position', get(0,'Screensize'));
% plot(expe.t(sel),bkg(longTraces(i),sel,colorIndex),'r--');
% title(expe.colorNames{colorIndex});
% Annot



setFonts
paperSize(80,80)
mkdirIfNotExist('figures')
% n2s(longTraces(i))
fname =['figures/Cell' n2s(longTraces(i)) '.pdf'];

 if doOpen
        print('-dpdf',fname)
%        system(['open ' fname])
 end

doAppend=0;
appendName
end

% close all

system(['open ' [pwd '/figures']])

% clear corRatio corRatios
%% plot trajectories 

clfh
for i=1:size(ind,1)
    sel= ind(i,:)>0;
    plot(trajX(i,sel),trajY(i,sel),'.','color',rand(3,1))
end


%% Load corrected peaks and divs

load divMatrixFinal.mat
load peakMatrixFinal.mat

clf
imagesc(peakMatrix(:,:,2)-0.1*divMatrix)

%% make nice time plot with images
doExplain=0;
doAnnot=0;
doOpen=1;

% for i=1:length(longTraces)
for i=1:length(longTraces)

clfh
idx = longTraces(i);
colorIndex = 2; %which images to show
sel = ind(longTraces(i),:) > 0;
% signalToPlot = refinedSum(:,:,colorIndex)-bkg(:,:,colorIndex).*refinedArea(:,:); %which signal to plot
signalToPlot = refinedMean(:,:,colorIndex)-bkg(:,:,colorIndex); %which signal to plot
% signalToPlot = normMean{colorIndex}(:,:);
showSeg = 1;

gaussianFilterSize = 0; %set to zero to disable
doNormalize = 0;

s  = 15;    % size of window around the cell
nR = 8;    % number of rows in the image

NtoPlot = N;
start =1;

useFullSizeImages = 0;
doDraw = 0;

makeImageTimePlot

mkdirIfNotExist('figures');
fname =['figures/cell' n2s(idx) '_' expe.colorNames{colorIndex} '.pdf'];
setFonts


if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end

doAppend=0;
appendName

end

clear doExplain doAnnot signalToPlot

%% correlation
clfh
doExplain=0;
doAnnot=0;
doOpen=1;


datato=zeros(length(longTraces),2*expe.numberOfFrames-1);

for i=1:length(longTraces)
% c = linspace(1,14,length(signal(longTraces(i),40:81,2)));
% scatter(refinedMean(longTraces(i),40:81,2)./refinedMean(longTraces(i),40:81,1),[],c,'filled')
% colorbar
% datato(i,:)=(refinedMean(longTraces(i),:,2)-bkg(longTraces(i),:,2))./(refinedpMean(longTraces(i),:,1)-bkg(longTraces(i),:,1));
 [r,lags] = xcov(refinedMean(longTraces(i),:,2)-bkg(longTraces(i),:,2),refinedMean(longTraces(i),:,1)-bkg(longTraces(i),:,1),'coeff');
plot(lags,r)
datato(i,:)=r;
% pause(1)
% close
% plot(datato(i,:))
hold all
% pause(2)
end

% plot(-149:149,mean(datato),'k','LineWidth',10)


% plotallgraphs(datato, expe.dt)
Annot
% plot(sum/9,'k','LineWidth',12)
paperSize(50,20)
xlim([-50 50])
xlabel('Lag')
ylabel('Cross_correlation')
setFonts

fname =['figures/cross_correlation' '.pdf'];
if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end

doAppend=0;
appendName

clear datato

%% nluc to luciferase ratio
clfh
doExplain=0;
doAnnot=0;
doOpen=1;


datato=zeros(length(longTraces),expe.numberOfFrames);

for i=1:length(longTraces)
%     clfh
% c = linspace(1,14,length(signal(longTraces(i),40:81,2)));
% scatter(refinedMean(longTraces(i),40:81,2)./refinedMean(longTraces(i),40:81,1),[],c,'filled')
% colorbar
datato(i,:)=(refinedMean(longTraces(i),:,2)-bkg(longTraces(i),:,2))./(refinedMean(longTraces(i),:,1)-bkg(longTraces(i),:,1));
%  [r,lags] = xcov(refinedMean(longTraces(i),:,2)-bkg(longTraces(i),:,2),refinedMean(longTraces(i),:,1)-bkg(longTraces(i),:,1),'coeff')
% plot(lags,r)
% pause(1)
% close
% plot(datato(i,:))
% ylim([0 5])
% xlim([0 181])
% pause(1)
% pause(2)
end


plotallgraphs(datato, expe.dt)
Annot
% plot(sum/9,'k','LineWidth',12)
paperSize(50,20)

xlabel('time [h]')
ylabel('Nluc/fluc')
setFonts
% ylim([0 3])

fname =['figures/ratio' '.pdf'];
if doOpen
print('-dpdf',fname)
system(['open ' fname])
end


%% correlation
clfh
doExplain=0;
doAnnot=0;
doOpen=1;


% getPortion(80,138,refinedMean)

startf=1
endf=296

datato=zeros(length(longTraces),2*(endf-startf)+1);



for i=1:length(longTraces)
% for i=1:10
[r,lags] = xcov(refinedMean(longTraces(i),startf:endf,2)-bkg(longTraces(i),startf:endf,2),refinedMean(longTraces(i),startf:endf,1)-bkg(longTraces(i),startf:endf,1),'coeff')
corrcoef(refinedMean(longTraces(i),startf:endf,2)-bkg(longTraces(i),startf:endf,2),refinedMean(longTraces(i),startf:endf,1)-bkg(longTraces(i),startf:endf,1))
% plot(datato(i,:))
% plot(lags,r)
% hold all
datato(i,:)=r;
% print('-dpdf',fname)
% close
% pause(2)
end

% 
plotallgraphs(datato, expe.dt)
Annot
% plot(sum/9,'k','LineWidth',12)
paperSize(50,20)
xlabel('time [h]')
ylabel('Nluc/fluc')
setFonts
% xlim([-50 50])
fname =['figures/ratio' '.pdf'];
if doOpen
print('-dpdf',fname)
system(['open ' fname])
end

%%
doAppend=0;
if doAppend
appendvec=[appendvec; fname];
end

%% temporal correlation coefficient
corr=[];
clfh
close all
doExplain=0;
doAnnot=0;
doOpen=1;
load longTraces
d1=2;
d2=4;
subplot(d1,d2,1)
for t=1:50
%   corr(t)= corr2(refinedSum(longTraces(:),t,1),refinedSum(longTraces(:),t,2)) 
  corr(t)= corr2(refinedSum(longTraces(:),t,1)-bkg(longTraces(:),t,1).*refinedArea(longTraces(:),t),refinedSum(longTraces(:),t,2)-bkg(longTraces(:),t,2).*refinedArea(longTraces(:),t)) ;
%     corr(t)= corr2(refinedMean(longTraces(:),t,1),refinedMean(longTraces(:),t,2))  
%     corr(t)= corr2(refinedMean(longTraces(:),t,1)-bkg(longTraces(:),t,1),refinedMean(longTraces(:),t,2)-bkg(longTraces(:),t,2))  

%   corr(t)= corr2(refinedSum(longTraces(:),t,1),refinedSum(longTraces(:),t,2))   
% corr(t)= corr2(normMean{1}(1:45,t),normMean{2}(1:45,t)) ; 
end  
plot(expe.t(:),corr)
title(['sum bkg subtracted'])
xlabel('tim
e [h]');
% xlim([0 7])
% ylim([0.7 1])
ylabel('corrrelation coefficient')
% ylim([0 2*max(refinedSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel))])


subplot(d1,d2,2)
for t=1:50
  corr(t)= corr2(refinedSum(longTraces(:),t,1),refinedSum(longTraces(:),t,2)) ;
end  
plot(expe.t(:),corr)
title(['sum bkg not subtracted' ])
xlabel('time [h]');
ylabel('corrrelation coefficient')
% ylim([0.9 1])
% xlim([0 7])
% ylim([0 2*max(refinedSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel))])


subplot(d1,d2,3)
for t=1:50
corr(t)= corr2(refinedMean(longTraces(:),t,1),refinedMean(longTraces(:),t,2)) ;
end  
plot(expe.t(:),corr)
title(['mean bkg not subtracted' ])
xlabel('time [h]');
ylabel('corrrelation coefficient')
% ylim([0.7 1])
% xlim([0 7])
% ylim([0 2*max(refinedSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel))])

% subplot(d1,d2,2)
% for t=1:50
%     corr(t)= corr2(refinedMean(longTraces(:),t,1)-bkg(longTraces(:),t,1),refinedMean(longTraces(:),t,2)-bkg(longTraces(:),t,2));
% end  
% plot(expe.t(:),corr)
% title(['mean bkg subtracted' ])
% xlabel('time [h]');
% ylabel('corrrelation coefficient')
% ylim([0.7 1])
% xlim([0 7])
% ylim([0 2*max(refinedSum(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex).*refinedArea(longTraces(i),sel))])

subplot(d1,d2,5)
for t=1:50
corrtemp= corrcoef([normMean{1}(:,t) normMean{2}(:,t)],'rows','complete') ; 
corr(t)=corrtemp(1,2);
clear corrtemp
end  
plot(expe.t(:),corr)
title(['normalized mean' ])
set(gca,'FontSize',6)
% xlabel('time [h]');
% ylabel('corrrelation coefficient')
% ylim([0.9 1])
% xlim([0 7])

subplot(d1,d2,6)
for t=1:50
corrtemp= corrcoef([normSum{1}(:,t) normSum{2}(:,t)],'rows','complete') ; 
corr(t)=corrtemp(1,2);
clear corrtemp
end  
plot(expe.dt*(1:50),corr)
title(['normalized sum' ])
% xlabel('time [h]');
% ylabel('corrrelation coefficient')
% ylim([0.9 1])
% xlim([0 7])


% setFonts
% set(gca,'FontSize', 18);
paperSize(50,200)
fname =['figures/temporal_correlations' '.pdf'];
if doOpen
print('-dpdf',fname)
system(['open ' fname])
end

% clear corr fname 

%% scatter plot correlation
clfh
doExplain=0;
doAnnot=0;
doOpen=1;

% mean 
for ii=10:10:50
for i=1:length(longTraces)
%  c = linspace(1,14,length(signal(longTraces(i),40:81,2)));
%  scatter(refinedMean(longTraces(i),40:81,2)./refinedMean(longTraces(i),40:81,1),[],c,'filled')
if refinedMean(longTraces(i),ii,2)>0
% scatter(refinedMean(longTraces(i),ii,1)-bkg(longTraces(i),ii,1),refinedMean(longTraces(i),ii,2)-bkg(longTraces(i),ii,2),'LineWidth',4)
 scatter(refinedMean(longTraces(i),ii,1)-bkg(longTraces(i),ii,1),refinedMean(longTraces(i),ii,2)-bkg(longTraces(i),ii,2),'LineWidth' ,4)
%  scatter(refinedSum(longTraces(i),ii,1)-bkg(longTraces(i),ii,1).*refinedArea(longTraces(i),ii),refinedSum(longTraces(i),ii,2)-bkg(longTraces(i),ii,2).*refinedArea(longTraces(i),ii),'LineWidth' ,4)
% scatter(normMean{1}(i,ii),normMean{2}(i,ii), 'LineWidth',4)
% corrcoef(normMean{1}(i,ii),normMean{2}(i,ii)) 
%  hold all
 end

end
 hold off
% xlim([0 6000])
% ylim([0 10000])
% plotallgraphs(datato, expe.dt)
Annot
% plot(sum/9,'k','LineWidth',12)
paperSize(50,20)
xlabel('fluc')
ylabel('Nluc')
setFonts

fname =['figures/ratio mean norm' num2str(expe.dt*ii) '.pdf'];
if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end
clfh
end

clear fname
%%
%%datafit
close all
doOpen=1;
xdata=[];
ydata=[];

% %refinedMean background subtracted
% xdata= refinedSum(longTraces(:),10,1)-bkg(longTraces(:),10,1).*refinedArea(longTraces(:),10);
% ydata= refinedSum(longTraces(:),10,2)-bkg(longTraces(:),10,2).*refinedArea(longTraces(:),10);
% indic=find(500000>xdata & xdata>0)
% xdata=xdata(indic)
% ydata=ydata(indic)

% % %corrected

for t=1:5:50
xdatatemp=[];
xdatatemp= refinedMean(longTraces(:),t,1)-bkg(longTraces(:),t,1);
xdata= [xdata;xdatatemp];

ydatatemp=[];
ydatatemp= (refinedMean(longTraces(:),t,2)-bkg(longTraces(:),t,2))/corRatio(t);
% ydatatemp= corRatio(t)*(refinedMean(longTraces(:),t,2)-bkg(longTraces(:),t,2));
ydata= [ydata;ydatatemp];
% indic=find(5000>xdata & xdata>0)
% xdata=xdata(indic)
% ydata=ydata(indic)
end
% 
% xdata=[refinedSum(longTraces(:),10,1)-bkg(longTraces(:),10,1).*refinedArea(longTraces(:),10);refinedSum(longTraces(:),30,1)-bkg(longTraces(:),30,1).*refinedArea(longTraces(:),30);refinedSum(longTraces(:),50,1)-bkg(longTraces(:),50,1).*refinedArea(longTraces(:),50); refinedSum(longTraces(:),70,1)-bkg(longTraces(:),70,1).*refinedArea(longTraces(:),70); refinedSum(longTraces(:),190,1)-bkg(longTraces(:),190,1).*refinedArea(longTraces(:),190)]
% ydata=[refinedSum(longTraces(:),10,2)-bkg(longTraces(:),10,2).*refinedArea(longTraces(:),10);refinedSum(longTraces(:),30,2)-bkg(longTraces(:),30,2).*refinedArea(longTraces(:),30);refinedSum(longTraces(:),50,2)-bkg(longTraces(:),50,2).*refinedArea(longTraces(:),50);refinedSum(longTraces(:),70,2)-bkg(longTraces(:),70,2).*refinedArea(longTraces(:),70);refinedSum(longTraces(:),190,2)-bkg(longTraces(:),190,2).*refinedArea(longTraces(:),190)]
% indic=find(1000000>xdata & xdata>10000)
% xdata=xdata(indic)
% ydata=ydata(indic)
% [hAx,hLine1,hLine2] = plotyy(expe.t(sel),refinedMean(longTraces(i),sel,1)-bkg(longTraces(i),sel,1),expe.t(sel),corRatio(sel).*(refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2)));

% %refinedMean background subtracted
% xdata= refinedMean(longTraces(:),30,1)-bkg(longTraces(:),30,1);
% ydata= refinedMean(longTraces(:),30,2)-bkg(longTraces(:),30,2);
% indic=find(5000>xdata & xdata>0)
% xdata=xdata(indic)
% ydata=ydata(indic)

f = fittype('a*x');

% [curve2,gof2]
[fit1,gof1] = fit(xdata,ydata,f,'StartPoint',[0]);
fdata = feval(fit1,xdata);
I = abs(fdata - ydata) > 1*std(ydata);
outliers = excludedata(xdata,ydata,'indices',I);

[fit2,gof2]= fit(xdata,ydata,f,'StartPoint',[0],'Exclude',outliers);
gof2
% plot(fit1,'r-',xdata,ydata,'k.',outliers,'m*')
hold on
scatter(xdata(~outliers),ydata(~outliers),'k','Linewidth',1)
hold on
plot(fit2,'c--')
hold on
% xlim([0 400000])
% n2s(longTraces(i))
%  k=1:length(xdata); text(xdata,ydata,num2str(k'))
paperSize(50,20)
xlabel('fluc')
ylabel('Nluc')
setFonts
text(1000, 4000, ['R^2 = ' num2str(getfield(gof2, 'rsquare'))],'FontSize',20)
% legend({'Rsquare=' n2s(getfield(gof2, 'rsquare'))})

fname =['figures/time point fit' num2str(expe.dt*ii) '.pdf'];
if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end


% clear xdata ydata outliers fit2 fit1 gof2 f indic fdata


%% noise between frames
clfh
doExplain=0;
doAnnot=0;
doOpen=1;

for ii=1:10:150
for i=1:length(longTraces)
%  c = linspace(1,14,length(signal(longTraces(i),40:81,2)));
%  scatter(refinedMean(longTraces(i),40:81,2)./refinedMean(longTraces(i),40:81,1),[],c,'filled')
 if refinedMean(longTraces(i),ii,2)>10 
scatter(refinedMean(longTraces(i),ii,1)-bkg(longTraces(i),ii,1),refinedMean(longTraces(i),ii+1,1)-bkg(longTraces(i),ii+1,1),'LineWidth',4)
 end
% pause(2)
end
hold off

% xlim([8000 15000])
% ylim([8000 15000])


% plotallgraphs(datato, expe.dt)
Annot
% plot(sum/9,'k','LineWidth',12)
paperSize(50,20)

xlabel('fluc')
ylabel('Nluc')
setFonts

fname =['figures/noise_frames' num2str(expe.dt*ii) '.pdf'];
if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end
clfh


end

%%
clfh
doExplain=0;
doAnnot=0;
doOpen=0;


datato=zeros(length(longTraces),2*expe.numberOfFrames-1);

for i=1:length(longTraces)
% c = linspace(1,14,length(signal(longTraces(i),40:81,2)));
% scatter(refinedMean(longTraces(i),40:81,2)./refinedMean(longTraces(i),40:81,1),[],c,'filled')
% colorbar
% datato(i,:)=(refinedMean(longTraces(i),:,2)-bkg(longTraces(i),:,2))./(refinedMean(longTraces(i),:,1)-bkg(longTraces(i),:,1));
 [r,lags] = xcov(refinedMean(longTraces(i),:,2)-bkg(longTraces(i),:,2),refinedMean(longTraces(i),:,1)-bkg(longTraces(i),:,1),'coeff')
plot(lags,r)
datato(i,:)=r;
% pause(1)
% close
% plot(datato(i,:))
hold all
% pause(2)
end

plot(-149:149,mean(datato),'k','LineWidth',10)


% plotallgraphs(datato, expe.dt)
Annot
% plot(sum/9,'k','LineWidth',12)
paperSize(50,20)
xlim([-50 50])
xlabel('Lag')
ylabel('Cross_correlation')
setFonts

fname =['figures/cross_correlation' '.pdf'];
if doOpen
print('-dpdf',fname)
% system(['open ' fname])
end

doAppend=0;
appendName

%%
corr=[];
clfh
for t=1:N
  corr(t)= corr2(refinedMean(longTraces(:),t,1)-bkg(longTraces(:),t,1),refinedMean(longTraces(:),t,2)-bkg(longTraces(:),t,2))  
end  
plot(expe.dt*(1:N),corr)
ylim([-1 1])
hold all
corr=[];
for t=1:149
  corr(t)= corr2(refinedMean(longTraces(:),t,1)-bkg(longTraces(:),t,1),refinedMean(longTraces(:),t+1,1)-bkg(longTraces(:),t+1,1))  
end  
plot(expe.dt*(1:N),corr)
ylim([-1 1])
hold all
corr=[];
for t=1:N
  corr(t)= corr2(refinedMean(longTraces(:),t,2)-bkg(longTraces(:),t,2),refinedMean(longTraces(:),t+1,2)-bkg(longTraces(:),t+1,2))  
end  
plot(expe.dt*(1:N),corr)
ylim([-1 1])

%%
%create appended pdf


fnamepdf='../movie1/figures/20151020.pdf';
% system(['rm ' fnamepdf])
append_pdfs(fnamepdf, appendvec{:})
system(['open ' fnamepdf])
appendvec={};

%% 
system(['cp -r ' fnamepdf ' /Users/onurtidin/Dropbox/presentations'])
