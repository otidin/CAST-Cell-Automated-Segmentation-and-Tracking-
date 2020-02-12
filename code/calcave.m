
function [ datout, ave ] = calcave( dat,dt)

lengthTrace= size(dat,2); 
traces=size(dat,1); %how many "good traces" did we have?
average=zeros(1,traces);

for i = 1:traces %plot every single cell in the same plot
    S1= dat(i,:);
    S1(find(S1==0))= NaN; %delete not defined values from the plot
%     plot(dt*[1:lengthTrace], S1(1:lengthTrace),'LineWidth',2); %5*[1:179] applies for 179 frames of 5 minutes exposure each
    hold all;
%      pause(0.5)
dat(i,:)=S1;
end

datout=dat;
ave = nanmean(dat,1);
hold off
