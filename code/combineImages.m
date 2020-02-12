background2=zeros(N,4);
background2(:,1)=1:N;
background2(:,3)=1:N;

for k=1:N/temporalBinning

    disp(100*k/N*temporalBinning);
          
    images = {};
    sorted= {};
    
    for i=1:Nz
        for j=1:temporalBinning                        
            images{i} = [outDir 'img/' getImageName(expe.colorNames{i},(k-1)*temporalBinning + j)];
        end
        im=[];
        im=imread(images{i});
        sorted{i}=sort(im(:),'ascend');
        background2(k,2*i)=mean(sorted{i}(1:25));
    end
    
%     images{1}=images{}
    clear im sorted
   
    [out,N1,N2] = combineStack(images,Nz,deNoise,medianSize,compressionQuantile,gaussianFilterSize,weightsSegmentation,doDraw,background,k);
    imwrite(out,['zStackedYFP/' num2str(k) '.png']);       
end


if(expe.hasTrans)
    
    for k=1:N/temporalBinning
        mkdirIfNotExist('zStackedTrans');

        images = {};
        for i=1:Nz
            for j=1:temporalBinning                        
                images{i} = [outDir 'img/' getImageName(expe.transName,(k-1)*temporalBinning + j)];
            end
        end

        [out,N1,N2] = combineStack(images,Nz,deNoise,medianSize,compressionQuantile,gaussianFilterSize,weightsSegmentation,doDraw);
        out = imnorm(out);
        
        imwrite(out,['zStackedTrans/' num2str(k) '.png']);
    end
end

% combine at real framerate for later
if temporalBinning > 1
    mkdirIfNotExist('zStackedYFP_unbinned');    
        
    for k=1:N
        disp(100*k/N);
        images = {};
        for i=1:Nz
            images{i} = [outDir 'img/' getImageName(expe.colorNames{i},k)];
        end

        [out,N1,N2] = combineStack(images,Nz,deNoise,medianSize,compressionQuantile,gaussianFilterSize,weightsSegmentation,0);
        imwrite(out,['zStackedYFP_unbinned/' num2str(k) '.png']);
    end
    
end

N = floor(N/temporalBinning);