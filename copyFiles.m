%%
restoredefaultpath
expe = experimentPara();
cd(expe.mainDir)
addpath ./code/
% doDraw=1;
naming = {'biop','bsf'};
naming = naming{1};




%%

switch naming
    
    case 'biop'

%% copy stacked tif

for i = expe.indexOfFirstMovie:(expe.indexOfFirstMovie + expe.numberOfMovies  -1)
   
    mkdirIfNotExist([expe.mainDir '/movie' num2str(i)]);
    mkdirIfNotExist([expe.mainDir '/movie' num2str(i) '/img/']);
    
    for j = 1:length(expe.colorNames)
        
        colorName = expe.colorNames{j};
        
        fname = getOriginalImageName(expe,colorName,i,1,naming);
    
        fname = ['''' fname ''''];            
        system(['cp -v ' fname ' ' expe.mainDir '/movie' num2str(i) '/img/']);
    
    end
    
    if( expe.hasTrans )
    
        colorName = expe.transName;
        
        fname = getOriginalImageName(expe,colorName,i,1,naming);
    
        fname = ['''' fname ''''];            
        system(['cp -v ' fname ' ' expe.mainDir '/movie' num2str(i) '/img/']);

    end
                
end

%% unpack tif files

doDraw = 0 ;
background=zeros(expe.numberOfFrames,length(expe.colorNames));

for i = expe.indexOfFirstMovie:(expe.indexOfFirstMovie + expe.numberOfMovies  -1)

    disp(i);
    
   

    for j = 1:length(expe.colorNames)
        
        colorName = expe.colorNames{j};
    
        fname = getOriginalImageName(expe,colorName,i,1,naming);

        info = imfinfo(fname);
        num_images = numel(info);
        
        
        for k = 1:expe.numberOfFrames

            zNumber = j;
            fNumber = k;

            A = imread(fname, k);

            
            if j==1
            
            colorName2 = expe.colorNames{2};
            fname2 = getOriginalImageName(expe,colorName2,i,1,naming);
            B=imread(fname2,k);
            sorted=sort(B(:),'ascend');
            background(k,2)=mean(sorted(1:30));
%             if doDraw==0
%             figure
%             imagesc(A)
%              figure
%             imagesc(B)
%             figure
%             imagesc(expe.subtractratio*(B-background(k,2)))
%             end
             A=A-expe.subtractratio*(B-background(k,2)); 
            sorted=sort(A(:),'ascend');
            background(k,1)=mean(sorted(1:30));
            
             clear B sorted fname2 colorName2
%             if doDraw==0
%             figure
%             imagesc(A)
%             pause
%             close all
%             end
            
            end
            

            if(doDraw && mod(k,5)==0)
                clf
                imagesc(A)
                drawnow
            end

            imgName = getImageName(colorName,fNumber);

            name = [expe.mainDir '/movie' num2str(i) '/img/' imgName];
            %name
            imwrite(A,name)

        end
    
    end
    
    save background.mat background
    
    if( expe.hasTrans )
    
        colorName = expe.transName;
        
        if i<10
            fname = [expe.mainDir '/movie' num2str(i) '/img/00' num2str(i) ' ' colorName '.tif'];
        elseif i<100
            fname = [expe.mainDir '/movie' num2str(i) '/img/0' num2str(i) ' ' colorName '.tif'];
        else
            fname = [expe.mainDir '/movie' num2str(i) '/img/' num2str(i) ' ' colorName '.tif'];
        end

        info = imfinfo(fname);
        num_images = numel(info);
        for k = 1:expe.numberOfFrames

            zNumber = j;
            fNumber = k;

            A = imread(fname, k);

            if(doDraw && mod(k,5)==0)
                clf
                imagesc(A)
                drawnow
            end

            imgName = getImageName(colorName,fNumber);

            name = [expe.mainDir '/movie' num2str(i) '/img/' imgName];
            %name
            imwrite(A,name)

        end
        
        %system(['rm "' fname '"']);

    end
    
            
    system(['mv ' expe.mainDir '/movie' num2str(i) '/img/*.tif ~/.Trash' ]);
    system(['mv ' expe.mainDir '/movie' num2str(i) '/img/*.TIF ~/.Trash' ]);

end


disp('done');
pause(5)

case 'bsf'
   
%% unstacked tif


doDraw = 1 ;


for i = expe.indexOfFirstMovie:(expe.indexOfFirstMovie + expe.numberOfMovies  -1)
    
    mkdirIfNotExist([expe.mainDir '/movie' num2str(i)]);
    mkdirIfNotExist([expe.mainDir '/movie' num2str(i) '/img/']);
    
    system(['sudo chmod +w ' expe.mainDir '/movie' num2str(i) '/img/'])
    system(['sudo chmod +w ' expe.mainDir '/movie' num2str(i) '/img/*.png'])
    
    disp(i);
    
    for j = 1:length(expe.colorNames)
        
        colorName = expe.colorNames{j};
    
        for k = 1:expe.numberOfFrames

            k
            fname = getOriginalImageName(expe,colorName,i,k,naming);    
            fname = ['''' fname ''''];  
            
            imgName = getImageName(colorName,k);

            name = [expe.mainDir '/movie' num2str(i) '/img/' imgName];
            name = ['''' name ''''];  
                        
            system(['cp -v ' fname ' ' name]);
        end
                    
    end
    
    if( expe.hasTrans )
        
        colorName = expe.transName
    
        for k = 1:expe.numberOfFrames
                        k
            fname = getOriginalImageName(expe,colorName,i,k,naming);    
            fname = ['''' fname ''''];  
            
            imgName = getImageName(colorName,k);

            name = [expe.mainDir '/movie' num2str(i) '/img/' imgName];
            name = ['''' name ''''];  
                        
            system(['cp -v ' fname ' ' name]);
        end
    end
        
    
end
   


end







