function expe = experimentPara()

%% don't change that

expe=struct();

%% change this

expe.numberOfFrames = 280;
expe.dt = 5/60;

expe.stim=6;

expe.subtractratio=0.0;


expe.numberOfMovies = 1;
expe.indexOfFirstMovie = 1;

expe.numberOfColors = 2;
expe.colorNames = {'fluc','Nluc'};%{'YFP','TexasRed'};

expe.hasTrans = 0;
expe.transName = 'Trans'; %use Trans for bsf, see getOriginalImageName.m

expe.imgDir  = '/Users/onurtidin/Desktop/desktopvirtual/AllMicroscopeExperiments/20180104new';
expe.mainDir = '/Users/onurtidin/Desktop/desktopvirtual/AllMicroscopeExperiments/20180104new';

%% don't change that

expe.t = linspace(0,(expe.numberOfFrames-1) * expe.dt,expe.numberOfFrames);
N = expe.numberOfFrames;



%%
