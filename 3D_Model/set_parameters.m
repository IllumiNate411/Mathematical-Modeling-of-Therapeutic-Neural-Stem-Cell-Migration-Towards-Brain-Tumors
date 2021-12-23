global n_seeds Finaltimestep CONVERT2MICRON coh_limit chmtx_limit d_w cancer_size  cancer_center ...
    has_cancer  modelNum modelType alpha4chmtx beta4dist chemo_sensitivity d_g FolderName1 FolderName2 

%%% Standard values
n_seeds =           1000;               %Total number of paths generated
Finaltimestep =     30000;              %Num of steps for each cell
CONVERT2MICRON =    13.5;               %Avg μm per pixel in our figure
coh_limit =         0.4;                %coherency threshold
chmtx_limit =       0.1;                %chemotaxis threshold
d_w =               1.5;                %Reference step size on white matter
cancer_size =       [100/CONVERT2MICRON, 400/CONVERT2MICRON, 100/CONVERT2MICRON]; %radius [x y z]


%%% Injection Center 
modelType =         'Intranasal'; 
% modelType =         'Intracerebral'; 


%%% Cancer Location
cancer_center =     [230, 300, 175];    % right, front putamen 
% cancer_center =     [120, 400, 150];    % right, back putamen 
% cancer_center =     [230, 300, 175;     % right, front putamen 
%                     400, 360, 155];    % left, middle putamen 


%%% Conditional
has_cancer =        1;
modelNum =          2;


%%% Variables
alpha4chmtx =       1;                  %alpha(beta=1)  parameter for WM vs chemo
beta4dist =         1;                  %beta(alpha=1)  parameter for distance
chemo_sensitivity = .5;                 %chmtx_bias sensitivity
d_g =               1;                  %reference step size on grey matter


FolderName1 = '/Figures/CancerNA/';     %Save plots to this folder when there is no cancer (Must pre-exist)
FolderName2 = '/Figures/Cancer/';       %Save plots to this folder when there is cancer (Must pre-exist)
