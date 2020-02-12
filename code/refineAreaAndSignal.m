%%create new folder for after refinement figures
mkdirIfNotExist('zStackedThreshCorrectedRefined')
 
%%the same trajectory is assigned to old trajectory
refinedTrajX = trajX;
refinedTrajY = trajY;

%delete existing binary images in refined folder
if( exist(['zStackedThreshCorrectedRefined/' num2str(k) '.png'],'file') )
   system('rm  zStackedThreshCorrectedRefined/*.png');
   system('rm  snapShots/*.png');
end


segext={};
segint={};

for n=1:length(longTraces)
%  for n=1:3   
%     n    
    %n = 3
    idx = longTraces(n);        

    clf;
    disp(100*n/length(longTraces))
   
    %window size adjustment
    w=2*s+1;
    nR = 6;
    tot = N;

    d = round(tot/nR);
    out = zeros(d*w,nR*w);
     mrefinedext={};
            mrefinedint={};
        
    %i = round(traj{idx}(:,1));
    %j = round(traj{idx}(:,2));
    
    for k=1:N
        disp(k);
        
        %%reading files from thresh corrected folder, create a new image
        %%with zeros
        
        
        if( ~exist(['zStackedThreshCorrectedRefined/' num2str(k) '.png'],'file') )
            m = imread(['zStackedThreshCorrected/' num2str(k) '.png']);
            imwrite(m<-1,['zStackedThreshCorrectedRefined/' num2str(k) '.png']);
            for extcount=1:length(dilateSizeAfterRefineext)
            imwrite(m<-1,['zStackedThreshCorrectedRefined/' num2str(k) 'ext' num2str(dilateSizeAfterRefineext(extcount)) '.png']);
            end
            for intcount=1:length(dilateSizeAfterRefineint)
            imwrite(m<-1,['zStackedThreshCorrectedRefined/' num2str(k) 'int' num2str(dilateSizeAfterRefineint(intcount)) '.png']);
            end
        end
        
        if( ind(idx,k)~=0 )
            
            a = imread(['zStackedYFP/' num2str(k) '.png']);

            data = {};
            for j=1:expe.numberOfColors            
                data{j} = imread( ['img/' getImageName(expe.colorNames{j},k)] );
            end
                        
            m = imread(['zStackedThreshCorrected/' num2str(k) '.png']);
            mrefined = imread(['zStackedThreshCorrectedRefined/' num2str(k) '.png']);
           
            for extcount=1:length(dilateSizeAfterRefineext)
            mrefinedext{extcount} = imread(['zStackedThreshCorrectedRefined/' num2str(k) 'ext' num2str(dilateSizeAfterRefineext(extcount)) '.png']);
            end
            for intcount=1:length(dilateSizeAfterRefineint)
            mrefinedint{intcount} = imread(['zStackedThreshCorrectedRefined/' num2str(k) 'int' num2str(dilateSizeAfterRefineint(intcount)) '.png']);
            end

            
            name = ['Measures/' num2str(k) '.mat'];
              
            %take the center position of the object to make a window
            pos = Me{k}(ind(idx,k)).Centroid;
            pos = round(pos);

            i = pos(1);
            j = pos(2);

            %[seli selj] = getNeiInd( i,j,s,N1,N2 );
            
            seli = (i-s):(i+s); 
            selj = (j-s):(j+s); 
            
            valid = seli > 0 & selj >0  & seli <= N2 & selj <= N1;

            %%choose only in border part of the window   
            seli = seli(valid);
            selj = selj(valid);
            
            [Y X] = meshgrid(selj-min(selj),seli-min(seli));

            sub_a = a(selj,seli);
            
            sub_data = {};
            for j=1:expe.numberOfColors  
                sub_data{j} = data{j}(selj,seli);
            end
             
            %imdilate to have a region of background
            bgkMask = m(selj,seli);
            bgkMask = imdilate(bgkMask,strel('disk',bgkSize));   
            
            if(doDrawBkg)
               clf 
               imagesc( imnorm(sub_a).*double(~bgkMask) )
               drawnow                
            end

            for j=1:expe.numberOfColors  

                 if(~isempty(sub_a(~bgkMask)))
%                     bkg(idx,k,j) = median( sub_data{j}(~bgkMask));
                    
                    [ma mi] = sort(sub_data{j}(~bgkMask));
%                     lowest5index = mi(1:5);
                    bkgsize=50;
                    if length(ma)>bkgsize
                    bkg(idx,k,j) = median(ma(1:50));
                    else
                    bkg(idx,k,j) = median(ma(:));
                    end
  
                    
                else
                    bkg(idx,k,j) = -1;
                end

            end
            
            px = Me{k}; %erase other objects
            px=px(ind(idx,k)).PixelIdxList; 
            m=double(m);
            m(px)=2;
            m = m==2;
            
            %before refinement object
            m = m(selj,seli);
            mbefore=m;
            %imagesc(m)
            %drawnow
            
            %imdilate to check if object touch the border of the frame
            updateArea = 1;
            mBord = imdilate(m,strel('disk',1));
            if( (sum(mBord(:,1)) + sum(mBord(:,end)) + sum(mBord(1,:)) + sum(mBord(end,:)))>0 )
                touchBorder(idx,k) = 1;                
                updateArea = 0;
            end
            
%             m = imdilate(m,strel('disk',3));   
    
            if( superSampling > 1)
                sub_a = imresize(sub_a,superSampling);
                
                for j=1:expe.numberOfColors  
                    sub_data{j} = imresize(sub_data{j},superSampling);
                end 
                m = imresize(m,superSampling);
            end
          
%             updateArea=0;
            if(updateArea)
                seg = region_seg(sub_a, m, NIteration,0.6,doDraw && mod(k,1)==0); %-- Run segmentation
                segrefined=seg;
            else
                seg = m;   
            end
          
            
            if( sum(seg(:)) == 0 )
                seg = m;
            end
            
            
            %imdilate for final segmented region
            if( dilateSizeAfterRefine >= 1 )
                seg = imdilate(seg,strel('disk',dilateSizeAfterRefine));
            end

                        
             mrefined(selj,seli) = mrefined(selj,seli)  +  imresize(seg,1./superSampling);
%             imwrite(mrefined,['zStackedThreshCorrectedRefined/' num2str(k) '.png']);
%             
            %update traj
            
            refinedTrajX(idx,k) = min(seli(:)) + mean(X(seg));
            refinedTrajY(idx,k) = min(selj(:)) + mean(Y(seg));
            %
            
            %%
             segext={};
             for extcount=1:length(dilateSizeAfterRefineext)
             segext{extcount} = imdilate(seg,strel('disk',dilateSizeAfterRefineext(extcount)));
%              mrefinedext(selj,seli) = mrefinedext(selj,seli)  +  imresize(segext{extcount},1./superSampling);
             mrefinedext{extcount}(selj,seli) = mrefinedext{extcount}(selj,seli)  +  imresize(segext{extcount},1./superSampling); 
            imwrite( mrefinedext{extcount},['zStackedThreshCorrectedRefined/' num2str(k) 'ext' num2str(dilateSizeAfterRefineext(extcount)) '.png']);
             end
          
             segint={};
             for intcount=1:length(dilateSizeAfterRefineint)
             segint{intcount} = imerode(seg,strel('disk',dilateSizeAfterRefineint(intcount)));
             mrefinedint{intcount}(selj,seli) = mrefinedint{intcount}(selj,seli)  +  imresize(segint{intcount},1./superSampling);
             imwrite(mrefinedint{intcount},['zStackedThreshCorrectedRefined/' num2str(k) 'int' num2str(dilateSizeAfterRefineint(intcount)) '.png']);
             end
            %%
%              if( dilateSizeAfterRefineint >= 1  )
%                 segint =imerode(seg,strel('disk',3));
%              end
%             
%             
%              mrefinedint(selj,seli) = mrefinedint(selj,seli)  +  imresize(segint,1./superSampling);
%              imwrite(mrefinedint,['zStackedThreshCorrectedRefined/' num2str(k) 'int.png']);
% 
%             %update traj
%             
%             refinedTrajXint(idx,k) = min(seli(:)) + mean(X(segint));
%             refinedTrajYint(idx,k) = min(selj(:)) + mean(Y(segint));
%             
            
%%            
            
          
            tmp = imclearborder(seg);

            if( sum(tmp(:)) > 0 )
                 seg = tmp;
            end

            %quantify different colors            
            for j=1:expe.numberOfColors  
                
                tmp = double(sub_data{j}(seg==1));            
                %qu  = quantile(tmp,0.99);    ql = quantile(tmp,0.01);
                %tmp = tmp(  (tmp < qu) & (tmp > ql) );
                refinedMean(idx,k,j)    = mean(tmp);
                refinedSum(idx,k,j)     = sum(tmp);
                refinedStd(idx,k,j)     = std(tmp); 
            end
            
            refinedArea(idx,k) = sum(seg(:))/superSampling^2;
            
            
            %make small image while we are at it   
            b1 = bwmorph(seg,'remove');
            sub_a(b1) = min(sub_a(:));
 
            
            %save image as png for GUI
            imwrite(sub_a,['snapShots/' num2str(idx) '_' num2str(k) '.png']);
            
 %% 
             for extcount=1:length(dilateSizeAfterRefineext)
            tmpext = imclearborder(segext{extcount});
            if( sum(tmpext(:)) > 0 )
                 segext{extcount} = tmpext;
            end
            %quantify different colors            
            for j=1:expe.numberOfColors     
                tmpext = double(sub_data{j}(segext{extcount}==1));            
                %qu  = quantile(tmp,0.99);    ql = quantile(tmp,0.01);
                %tmp = tmp(  (tmp < qu) & (tmp > ql) );
                refinedMeanextall(extcount,idx,k,j)    = mean(tmpext);
                refinedSumextall(extcount,idx,k,j)     = sum(tmpext);
                refinedStdextall(extcount,idx,k,j)     = std(tmpext); 
            end
            refinedAreaextall(extcount,idx,k) = sum(segext{extcount}(:))/superSampling^2;
            %make small image while we are at it   
            b1ext = bwmorph(segext{extcount},'remove');
            sub_a(b1ext) = min(sub_a(:));
            %save image as png for GUI
            imwrite(sub_a,['snapShots/' num2str(idx) '_' num2str(k) '_ext' num2str(extcount) '.png']);
             end
          %%    
            for intcount=1:length(dilateSizeAfterRefineint)
            tmpint = imclearborder(segint{intcount});
            if( sum(tmpint(:)) > 0 )
                 segint{intcount} = tmpint;
            end
            %quantify different colors            
            for j=1:expe.numberOfColors    
                tmpint = double(sub_data{j}(segint{intcount}==1));            
                %qu  = quantile(tmp,0.99);    ql = quantile(tmp,0.01);
                %tmp = tmp(  (tmp < qu) & (tmp > ql) );
                refinedMeanintall(intcount,idx,k,j)    = mean(tmpint);
                refinedSumintall(intcount,idx,k,j)     = sum(tmpint);
                refinedStdintall(intcount,idx,k,j)     = std(tmpint); 
            end 
            refinedAreaintall(intcount,idx,k) = sum(segint{intcount}(:))/superSampling^2;
            %make small image while we are at it   
            b1int = bwmorph(segint{intcount},'remove');
            sub_a(b1int) = min(sub_a(:));
            % imagesc(sub_a);
            %save image as png for GUI
            imwrite(sub_a,['snapShots/' num2str(idx) '_' num2str(k) '_int' num2str(intcount) '.png']);
           end
        end
    end
    
    
    clf;
    subplot_tight(2,1,1,0.1); hold on
        plot(areaMatrix(idx,:),'k');
        plot(refinedArea(idx,:),'r');
        plot(10*touchBorder(idx,:),'g');
        title('area')
    subplot_tight(2,1,2,0.1); hold on
        %plot(refinedMean(idx,:)-refinedStd(idx,:),'r--');
        %plot(refinedMean(idx,:)+refinedStd(idx,:),'r--');

        plot(refinedMean(idx,:,1),'r');        
        plot(signal(idx,:,1),'k')
        plot(bkg(idx,:,1),'r--')
        
        if(size(signal,3)>1)
            plot(signal(idx,:,2),'k')
            plot(refinedMean(idx,:,2),'g');
            plot(bkg(idx,:,2),'g--')
        end
        
        title('Signal and background')
                
    drawnow
    
end



save refinedMean.mat refinedMean
save refinedSum.mat refinedSum
save refinedStd.mat refinedStd
save bkg.mat bkg
save refinedArea.mat refinedArea
save touchBorder.mat touchBorder