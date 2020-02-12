function [ ave ] = plotallgraphs( dat,dt,varargin)
%PLOTALLGRAPHS plots all traces of an analysis in one plot and superimposes
%them.
%   It needs two input parameters from the analysis - longTraces and signal


lengthTrace= size(dat,2); 
traces=size(dat,1); %how many "good traces" did we have?
average=zeros(1,traces);


switch nargin
    case 2      

for i = 1:traces %plot every single cell in the same plot
    S1= dat(i,:);
    S1(find(S1==0))= NaN; %delete not defined values from the plot
 plot(dt*[1:lengthTrace], S1(1:lengthTrace),'LineWidth',2); %5*[1:179] applies for 179 frames of 5 minutes exposure each
    hold all;
%      pause(0.5)
end
ave = nanmean(dat,1);
plot(dt*[1:lengthTrace], ave(1:lengthTrace),'k','LineWidth',7);
setFonts
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean.pdf'];
print('-dpdf',fname)
% system(['open ' fname])


% mean_y = nanmean(dat,1);
% std_y = nanstd(dat,1);
% x=dt*[1:lengthTrace];
% figure;hold on;
% 
% H1 = plot(x, mean_y, 'Color', 'k', 'LineWidth', 5);
% H2 = plot(x, [mean_y - 0.25*std_y; mean_y + 0.25*std_y], 'Color', 'r');
% H3 = plot(x, [mean_y - std_y; mean_y + std_y], 'Color', 'm');
% H4 = plot(x, [mean_y - 0.5*std_y  ; mean_y +   0.5*std_y], 'Color', 'b');
% 
% legend([H1, H2(1), H3(1), H4(1)], ...
%     '\mu', '2\sigma', '\sigma', '0.5\sigma', ...
%     'Location', 'Northwest');

figure;hold on;
y=dat;
x=dt*[1:lengthTrace];
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

 case 3
  
   dat2=varargin{1} ;
   
   for i = 1:traces %plot every single cell in the same plot
    S1= dat(i,:);
    S1(find(S1==0))= NaN; %delete not defined values from the plot
 plot(dt*[1:lengthTrace], S1(1:lengthTrace),'LineWidth',2); %5*[1:179] applies for 179 frames of 5 minutes exposure each
    hold all;
%      pause(0.5)
   end
colorIndex=1;
ave = nanmean(dat,1);
plot(dt*[1:lengthTrace], ave(1:lengthTrace),'k','LineWidth',7);
setFonts
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean_' n2s(1) '.pdf'];
print('-dpdf',fname)


figure;hold on;
y=dat;
x=dt*[1:lengthTrace];
% H(1) = shadedErrorBar(x, y, {@nanmean, @(x) 2*nanstd(x)  }, '-r', 0);
% H(2) = shadedErrorBar(x, y, {@nanmean, @(x) 1*nanstd(x)  }, '-m', 0);
H(3) = shadedErrorBar(x, y, {@nanmean, @(x) 0.5*nanstd(x)}, {'-b', 'LineWidth', 5}, 0);
legend([H(3).mainLine, H.patch], ...
    '\mu','0.5\sigma', ...
    'Location', 'Northeast');
setFonts
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean_shaded_' n2s(1) '.pdf'];
print('-dpdf',fname)


figure
for i = 1:traces %plot every single cell in the same plot
    S2= dat2(i,:);
    S2(find(S2==0))= NaN; %delete not defined values from the plot
 plot(dt*[1:lengthTrace], S2(1:lengthTrace),'LineWidth',2); %5*[1:179] applies for 179 frames of 5 minutes exposure each
    hold all;
%      pause(0.5)
   end
colorIndex=2;
ave2 = nanmean(dat2,1);
plot(dt*[1:lengthTrace], ave2(1:lengthTrace),'k','LineWidth',7);
setFonts
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean_' n2s(2) '.pdf'];
print('-dpdf',fname)
  
figure;hold on;
y=[];
y=dat2;
x=dt*[1:lengthTrace];
% H(1) = shadedErrorBar(x, y, {@nanmean, @(x) 2*nanstd(x)  }, '-r', 0);
% H(2) = shadedErrorBar(x, y, {@nanmean, @(x) 1*nanstd(x)  }, '-m', 0);
H(3) = shadedErrorBar(x, y, {@nanmean, @(x) 0.5*nanstd(x)}, {'-b', 'LineWidth', 5}, 0);
legend([H(3).mainLine, H.patch], ...
    '\mu','0.5\sigma', ...
    'Location', 'Northeast');
setFonts
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean_shaded_' n2s(2) '.pdf'];
print('-dpdf',fname)


       
figure
shadedErrorBaryy(x,nanmean(dat,1),0.5*nanstd(dat,1),'r',x,nanmean(dat2,1),0.5*nanstd(dat2,1),'b')
setFonts
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean_shaded.pdf'];
print('-dpdf',fname)


figure
 [hAx,hLine1,hLine2] = plotyy(dt*[1:lengthTrace],ave,dt*[1:lengthTrace],ave2);
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
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean_both.pdf'];
print('-dpdf',fname)
    
end

system(['open ' [pwd '/figures']])
end

