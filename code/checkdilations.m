ratioscell={};
for i=1:length(dilateSizeAfterRefineext)

    for j=1:length(dilateSizeAfterRefineint)
        
refinedMeanext=squeeze(refinedMeanextall(i,:,:,:));
refinedSumext=squeeze(refinedSumextall(i,:,:,:));
refinedAreaext=squeeze(refinedAreaextall(i,:,:));

refinedMean=squeeze(refinedMeanintall(j,:,:,:));
refinedSum=squeeze(refinedSumintall(j,:,:,:));
refinedArea=squeeze(refinedAreaintall(j,:,:));

refinedMeanextb=refinedMeanext;
refinedSumextb=refinedSumext;
refinedMeanb=refinedMean;
refinedSumb=refinedSum;

% define signals(subtract background and minimum adjustment)
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
% hist(refinedSumtemp(1,:,2))
% 
refinedSumexttemp(:,:,1)=refinedSumextb(:,:,1)-bkg(:,:,1).*refinedAreaext(:,:);
refinedSumexttemp(:,:,2)=refinedSumextb(:,:,2)-bkg(:,:,2).*refinedAreaext(:,:);
refinedSumexttemp(:,:,2) =refinedSumexttemp(:,:,2)-(min(min(refinedMeanexttemp(longTraces(:),:,2)'))-minvalue).*refinedAreaext(:,:);

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

annotvec=[18];
appendvec={};
stim=[18];

% i=23;
% figure
% plot(refinedMeancyto(longTraces(i),:,2))
% hold all
% plot(refinedMean(longTraces(i),:,2))
% legend('cyto','nuc')

%%clfh
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

ratioscell{j+(i-1)*2}=ave1;
plot(expe.dt*[tstart:tend], ave1,'k','LineWidth',5);
j+(i-1)*2
setFonts
paperSize(50,20)
mkdirIfNotExist('figures')
fname =['figures/alltraces_mean.pdf'];
print('-dpdf',fname)

% figure
% hold all
% plot(expe.dt*[tstart:tend], data2);
% plot(expe.dt*[tstart:tend], ave2,'k','LineWidth',7);
% setFonts
% paperSize(50,20)
% mkdirIfNotExist('figures')
% fname =['figures/alltraces_mean2.pdf'];
% print('-dpdf',fname)
% 
% figure;hold on;
% y=data1;
% x=expe.dt*[tstart:tend];
% % H(1) = shadedErrorBar(x, y, {@nanmean, @(x) 2*nanstd(x)  }, '-r', 0);
% % H(2) = shadedErrorBar(x, y, {@nanmean, @(x) 1*nanstd(x)  }, '-m', 0);
% H(3) = shadedErrorBar(x, y, {@nanmean, @(x) 0.5*nanstd(x)}, {'-b', 'LineWidth', 5}, 0);
% legend([H(3).mainLine, H.patch], ...
%     '\mu','0.5\sigma', ...
%     'Location', 'Northeast');
% setFonts
% paperSize(50,20)
% mkdirIfNotExist('figures')
% fname =['figures/alltraces_mean_shaded.pdf'];
% print('-dpdf',fname)
% 
% figure;
% y=data2;
% x=expe.dt*[tstart:tend];
% % H(1) = shadedErrorBar(x, y, {@nanmean, @(x) 2*nanstd(x)  }, '-r', 0);
% % H(2) = shadedErrorBar(x, y, {@nanmean, @(x) 1*nanstd(x)  }, '-m', 0);
% H(3) = shadedErrorBar(x, y, {@nanmean, @(x) 0.5*nanstd(x)}, {'-b', 'LineWidth', 5}, 0);
% legend([H(3).mainLine, H.patch], ...
%     '\mu','0.5\sigma', ...
%     'Location', 'Northeast');
% setFonts
% paperSize(50,20)
% mkdirIfNotExist('figures')
% fname =['figures/alltraces_mean2_shaded.pdf'];
% print('-dpdf',fname)
% 
% 
% figure
%  [hAx,hLine1,hLine2] = plotyy(expe.dt*[tstart:tend],ave1,expe.dt*[tstart:tend],ave2);
% % plot(expe.t(sel),corRatio(sel).*(refinedMean(longTraces(i),sel,colorIndex)-bkg(longTraces(i),sel,colorIndex)));
% % title(['After Correction Nluc and fluc '])
% ylabel(hAx(1),'nuc/cell ratio') % left y-axis
% ylabel(hAx(2),'nuc/cyto ratio') % right y-axis
% set(hAx,{'ycolor'},{'k';'b'})
% set(hLine1,'linewidth',7)
% % set(get(hAX(2),'color'),'black');
% set(hLine2,'linewidth',7)
% set(hLine1,'color','black');
% set(hLine2,'color','blue');
% xlabel('Time (h)')
% setFonts
% % ylim(hAx(2), [0 2*max(corRatio(sel).*(refinedMean(longTraces(i),sel,2)-bkg(longTraces(i),sel,2)))]);
% % ylim(hAx(1), [0.9 1.1]);
% %  ylim(hAx(2), [0.9 1.1]);
% doAnnot=1;
% Annot
% % xlim([1 12])
% % ylim([0 1000])
% setFonts
% paperSize(50,20)
% fname =['figures/alltraces_mean_tworatios.pdf'];
% print('-dpdf',fname)
% % 
% 
% figure
% x=expe.dt*[tstart:tend];
% shadedErrorBaryy(x,nanmean(data1,1),0.5*nanstd(data1,1),'r',x,nanmean(data2,1),0.5*nanstd(data2,1),'b')
% setFonts
% paperSize(50,20)
% mkdirIfNotExist('figures')
% fname =['figures/alltraces_mean_shadedtwo.pdf'];
% print('-dpdf',fname)
% 

close all
if j+(i-1)*j~=1
   openfig(fnamesaved)
end
hold all
% subplot(length(dilateSizeAfterRefineext)*length(dilateSizeAfterRefineint),1,j+(i-1)*j)

plot(expe.dt*[tstart:tend],ave1,'LineWidth',5)
fnamesaved =['figures/trial.fig'];
savefig(fnamesaved)
close(gcf)
% hold off
if j+(i-1)*j==length(dilateSizeAfterRefineext)*length(dilateSizeAfterRefineint)
  openfig(fnamesaved) 
  fnamepdf =['figures/differentdilation.pdf'];
    legend('cyto 3 pixel dilation nuc 1 pixel erosion','cyto 3 pixel dilation nuc 3 pixel erosion',...
        'cyto 8 pixel dilation nuc 1 pixel erosion','cyto 8 pixel dilation nuc 3 pixel erosion')
      paperSize(50,20)
  setFonts
print('-dpdf',fnamepdf)
end

close all
if j+(i-1)*j~=1
   openfig(fnamesavedctgf)
end
hold all
% subplot(length(dilateSizeAfterRefineext)*length(dilateSizeAfterRefineint),1,j+(i-1)*j)

plot(expe.dt*[tstart:tend],ave2)
fnamesavedctgf =['figures/trialctgf.fig'];
savefig(fnamesavedctgf)
close(gcf)
% hold off
if j+(i-1)*j==length(dilateSizeAfterRefineext)*length(dilateSizeAfterRefineint)
  openfig(fnamesavedctgf)  
  legend('','', 'nuc 1 pixel erosion','nuc 3 pixel erosion')
  fnamepdf =['figures/differentdilationctgf.pdf'];
  paperSize(50,20)
  setFonts
print('-dpdf',fnamepdf)
end
    end
end

% system('open .')