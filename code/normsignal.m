function [ ave, S1norm] = normsignal( dat ,dt)
%PLOTALLGRAPHS plots all traces of an analysis in one plot and superimposes
%them.
%   It needs two input parameters from the analysis - longTraces and signal

% close all

lengthTrace= size(dat,2); 
traces=size(dat,1); %how many "good traces" did we have?
% total=zeros(1,lengthTrace);
% total2=zeros(1,lengthTrace);
average=zeros(1,traces);

% paperSize(50,20)

for i = 1:traces %plot every single cell in the same plot
    
%     S= signal(longTraces(i),:); 
%     
%     S(find(S==0))= total(find(S==0))/(i-1);
%     if (i==1)
%     S=zeros(1,lengthTrace);
%     end
%     total=total+S; % calculate the average of the signal
    
    S1= dat(i,:);
    S1(find(S1==0))= NaN; %delete not defined values from the plot
%     plot(dt*[1:lengthTrace], S1(1:lengthTrace),'LineWidth',10); %5*[1:179] applies for 179 frames of 5 minutes exposure each
%     total2=total2+S1; % calculate the average of the signal
    
%     hold all;
%      pause(0.5)
end



for j=1:size(dat,2)
    tot=0;
    count=0;
    for k=1:size(dat,1)
        
    if any(dat(k,j)) 
     
%     else 
      tot=tot+dat(k,j);  
      count=count+1;
    end
    
    end
    average(j)=tot/count;
    
end
ave=average;
% 
% hold all
% plot(dt*[1:lengthTrace],average,'k','LineWidth',10)


for i = 1:traces %plot every single cell in the same plot
    
%     S= signal(longTraces(i),:); 
%     
%     S(find(S==0))= total(find(S==0))/(i-1);
%     if (i==1)
%     S=zeros(1,lengthTrace);
%     end
%     total=total+S; % calculate the average of the signal
    
    S1= dat(i,:);
    S1(find(S1==0))= NaN; %delete not defined values from the plot
%     plot(dt*[1:lengthTrace], S1(1:lengthTrace)./average,'LineWidth',10); %5*[1:179] applies for 179 frames of 5 minutes exposure each
    S1norm(i,:)= S1(1:lengthTrace)./average;
%     total2=total2+S1; % calculate the average of the signal
    
%     hold all;
%      pause(0.5)
end




% figure
% y = mean(dat,1);
% e = std(dat);
% errorbar(y,e)
% xlim([13 25])


% S1= signal(longTraces,:);
% plot(20*[1:lengthTrace],mean(S1))

% average=total/i;
% 
% average2=total2/i;
% plot(20*[1:lengthTrace], average2,'LineWidth',7);

% hold off;
%  ylim([9000 14000]);
%   xlim([3550 4000]);
% set(gca,'fontsize',18);
% ylabel('Luminescence [AU]','FontSize',16);
% xlabel('Time(min)','FontSize',16);
%legend('show');
% set(gca,'fontsize',18);
% xlabel('Time(min)','FontSize',20);
% ylabel('Luminescence Intensity(au)');
% xlim([0 1000]);
% saveas(gcf,['Cell2.jpg'])
%  ylim([0 18000])
% pause(0.5)
end

