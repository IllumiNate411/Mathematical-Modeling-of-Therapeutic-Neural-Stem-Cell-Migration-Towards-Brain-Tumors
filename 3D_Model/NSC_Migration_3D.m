clear all;
close all;

disp('Gathering data...')

%% Set Parameters
% Global vars initialized:  n_seeds Finaltimestep CONVERT2MICRON coh_limit chmtx_limit d_w cancer_size cancer_center ...
%                           modelType alpha4chmtx beta4dist chemo_sensitivity d_g
set_parameters()

%% Load Anistropy Files
% Global vars initialized:  eigen_map coh_map
load_3D();

%% Set Initial
% Global vars initialized:  seed_sd ubound lbound inj_center has_cancer
[seed_ind,has_cancer] = set_initial();         

%% Set Cancer
% Global vars initialized:  num_cancers curr_site_num
if has_cancer 
    [concentration, cgradX, cgradY, cgradZ] = set_cancer();
end

%% Simulation
global eigen_map coh_map ubound lbound inj_center num_cancers curr_site_num has_cancer

Tstart = tic; 
disp('Running simulation...')

p = struct('coord',{});

for seed_loop = 1:n_seeds
    
    %%%% Initialize particular reference step length
    dmax = betainv(rand(1), 1, beta4dist);
        
    
    %%%% Spawn seed within injection site
    i = 1;
    seed = seed_ind(ceil(rand*length(seed_ind)));                                                           % Choose a random point index to start the simulation
    [seed(1,1),seed(1,2),seed(1,3)] = ind2sub(size(coh_map),seed);

    if seed(3) > ubound(seed(1),seed(2))
        seed(3) = ubound(seed(1),seed(2));
    elseif seed(3) < lbound(seed(1),seed(2))
        seed(3) = lbound(seed(1),seed(2));
    end
    p(seed_loop).coord(:,i) = seed';                                                                        % Coordinate of this starting point

    
    %%%% First step 
    i = 2; 
    eigen_vect_noise = [rand(1)*2-1, rand(1), rand(1)]'; 
    eigen_vect = eigen_vect_noise/norm(eigen_vect_noise); 
    p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + d_g*eigen_vect;                                   % Find next coordinate at step d along EV direction

    seed = round( p(seed_loop).coord(:,i) );
    if seed(3) > ubound(seed(1),seed(2)) || seed(3) < lbound(seed(1),seed(2))
        p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1); 
    end
    
    
    %%%% Algorithm: while NSC stays in bounds
    for i = 3:Finaltimestep                                                                                 % Run the loop a large number of steps            
        
        %%% Find index
        seedx = round( p(seed_loop).coord(1,i-1) );                                                         % Initialize (seedx,seedy,seedz) as (x,y,z) coords of prev step
        seedy = round( p(seed_loop).coord(2,i-1) );
        seedz = round( p(seed_loop).coord(3,i-1) );
        ind = sub2ind(size(coh_map),seedx,seedy,seedz);
        
        
        %%% Calculate eigen vector, step length, and prefered direction +/-
        if coh_map(ind) > coh_limit                                                                         % Get Eigen vector at the previous coordinate point
     
            eigen_vect = squeeze( eigen_map( seedx, seedy, seedz, : ) ); 
            d = dmax * d_w;
            pref_mvmt = sign(dot( p(seed_loop).coord(:,i-1) - p(seed_loop).coord(:,i-2), eigen_vect) + eps );   % Move either in the +/- direction
        else  

            eigen_vect_noise = (rand(3,1)*2-1); 
            eigen_vect = eigen_vect_noise/norm(eigen_vect_noise); 
            d = dmax * d_g;                                                                                 % Move with less magnitude if no WM present 
            pref_mvmt = 1;
        end 

        dstep = d; 
        
        %%% Calculate chemotaxis vector and gradient strength
        chmtx_vect = [0;0;0];
        chmtx_bias = 0;
        if has_cancer
            chmtx_vect = [cgradX(ind); cgradY(ind); cgradZ(ind)];
            chmtx_bias = betainv(rand(1),alpha4chmtx,1) * chemo_sensitivity;
        end
        
        
        %%% Take step: (curr coord) = (prev coord) + (coherency influence) + (cancer influence)
        p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + pref_mvmt*dstep*eigen_vect + chmtx_bias*chmtx_vect;
        
        
        %%% If cell is out of bounds, put it back to previous coordinate
        seednow = round(p(seed_loop).coord(:,i));
        if ~( all(seednow' - size(coh_map)<=-1) && all( seednow'>1 ) ) || ...                                           %Seed reaches edge of coh_map
            ( seednow(3) > ubound(seednow(1),seednow(2)) || seednow(3) < lbound(seednow(1),seednow(2)) ) || ...         %Seed is outside of the lower or upper bounds (z)
            ( isnan(ubound(seednow(1),seednow(2))) || isnan(lbound(seednow(1),seednow(2))) )                            %Seed is outside of the side bounds of brain (x,y)

            p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1); 
        end

        
    end    
end

toc( Tstart ); 




%% Plot Graphs (use smaller p)
disp('Plotting graphs...')

if exist('p_original','var') && exist('acc','var')
    for n=1:n_seeds
        if p(n).coord(:,2) ~= p_original(n).coord(:,1+acc)
            p_original = p; break
        end
    end
    p = p_original;
end

acc = 100;                                                                                                  %Use a smaller timestep (record every acc number steps)
p_original = p;
p_new = struct('coord',{});
for i=1:n_seeds
    p_new(i).coord = p(i).coord(:,1:acc:Finaltimestep); 
end
p = p_new; clear p_new;
xvals = (1:acc:Finaltimestep);

save(savethis('p','mat'), 'p');                                                                             %Save coordinates of each step for all NSC




%% Plot path of seeds
figure;
plotenvironment(.7);

for i=1: n_seeds
    plot3(p(i).coord(1,:),p(i).coord(2,:),p(i).coord(3,:),'linewidth',3);                                       %Plot path
    plot3(p(i).coord(1, end),p(i).coord(2, end),p(i).coord(3, end),'ok','Markersize', 7,'markerfacecolor','w'); %Plot final coords
end

savethis('trajectoryAll','fig');
savethis('trajectoryAll','jpg');




%% Plot final position 
figure;
plotenvironment(.7);

for i=1: n_seeds
    plot3(p(i).coord(1, end),p(i).coord(2, end),p(i).coord(3, end),'ok','Markersize', 7,'markerfacecolor','w'); %Plot final coords
end

savethis('finalPosition','jpg');




%% Boxplots of NSC's distance from center of injection site
clear distInit;
distInit = zeros(n_seeds,Finaltimestep/acc);
for n = 1:n_seeds
    distInit(n,:) = sqrt( sum ( (p(n).coord - inj_center'*ones(1,size(p(n).coord,2))).^2, 1 ) );
end

num_plots = 10;
data = distInit(:, (1:num_plots)*floor(Finaltimestep/(num_plots*acc)))*CONVERT2MICRON;
labels = {num_plots};
for i=1:num_plots
    labels{i} = num2str( floor(Finaltimestep/num_plots)*i );
end
figure;     boxplot( data, labels );
hold on;    plot( (1:acc:Finaltimestep)*(num_plots/Finaltimestep), median( distInit )*CONVERT2MICRON );

ylim([0 max(max(data))+500]);
xlim([0 num_plots+1]);
xlabel( 'Time Intervals' ); ylabel( 'Distance from injection site (\mu)' );
savethis('distboxplot','jpg');




%% Percent of NSC on WM
numAtWM = zeros(Finaltimestep/acc,1);
for i = 1:n_seeds
    seedx = round(p(i).coord(1,:)); 
    seedy = round(p(i).coord(2,:)); 
    seedz = round(p(i).coord(3,:)); 
    ind = sub2ind(size(coh_map),seedx, seedy, seedz); 
    numAtWM = numAtWM + double( coh_map(ind) > coh_limit )'; 
end
percOnWM = numAtWM / n_seeds * 100;

figure;   plot(xvals,percOnWM,'r');

xlabel('Time Intervals');   ylabel('Percent of NSC on WM');
savethis('percentOnWM','jpg');




%% Determine whether NSC made it within a certain radius to cancer center
if has_cancer
    for k=1: num_cancers
        percAtCancerGraph = zeros(Finaltimestep/acc,1);
        for i = 1:Finaltimestep/acc
            numAtCancer = 0;
            for j = 1:n_seeds
                if all( abs(cancer_center(k,:) - p(j).coord(:,i)') < cancer_size )
                    numAtCancer = numAtCancer + 1;
                end
            end
            percAtCancerGraph(i,1) = numAtCancer;        
        end
        percAtCancerGraph = 100*percAtCancerGraph/n_seeds;

        indAtCancer = []; 
        for i = 1:n_seeds
            if all( abs(cancer_center - p(i).coord(:,end)') < cancer_size )
                indAtCancer = [indAtCancer,i]; 
            end
        end
        
        figure; plot(xvals,percAtCancerGraph,'r');
        xlabel('Time Intervals');   ylabel('Percent of NSC that Reach Cancer Site');
        
        %Must specify curr_site_num for correct naming before calling savethis
        curr_site_num=k;
        savethis('percentAtCancer','jpg');
    end
end




%% Min/Max Path
stepDistance = zeros(n_seeds, Finaltimestep/acc - 1);
for i = 1:n_seeds
    for j = 2:Finaltimestep/acc
        % record distance traveled each step
        stepDistance(i, j - 1) = sqrt((p(i).coord(1,j) - p(i).coord(1,j-1)).^2 + (p(i).coord(2,j) - p(i).coord(2,j-1)).^2 + (p(i).coord(3,j) - p(i).coord(3,j-1)).^2);
    end
end

% compute cell trajectory lengths 
stepDistanceTotal = sum( stepDistance, 2 ); 

% find minimum path among all cells
minPath = min( stepDistanceTotal(:) ); 
minPath_seed = find( stepDistanceTotal == minPath );

% find maximum path among all cells
maxPath = max( stepDistanceTotal(:) ); 
maxPath_seed = find( stepDistanceTotal == maxPath );




%% Plot Min/Max
if has_cancer
    for i=1: num_cancers
        cellind = [];
        for j = 1:n_seeds
            if (all( abs(cancer_center(i,:) - p(j).coord(:,end)') < cancer_size ))
                cellind = [cellind, j];
            end
        end 
        
        if ~isempty(cellind)
            dist2cancer = zeros(n_seeds, 1); 

            % Min Calculations
            % find minimum path among those that are close to cancer 
            minPath_cancer = min( stepDistanceTotal(cellind) ); 
            minPath_cancer_seed = find( stepDistanceTotal == minPath_cancer ); 

            % find cells that are close enough to minimum path to cancer 
            mincloseindex = find( abs( stepDistanceTotal-minPath ) < 100 ); 

            % percentage of cells that are near minimum path to cancer
            minclosepercent = length(mincloseindex) / length(stepDistanceTotal);

            % Max Calculations
            % find maximum path among cells that are close to cancer 
            maxPath_cancer = max( stepDistanceTotal(cellind) ); 
            maxPath_cancer_seed = find( stepDistanceTotal == maxPath_cancer ); 

            % find cells that are close enough to maximum path 
            maxcloseindex = find( abs( stepDistanceTotal-maxPath ) < 100 ); 

            % percentage of cells that are near maximum path 
            maxclosepercent = length(maxcloseindex) / length(stepDistanceTotal);

            figure;
            plotenvironment(1);
            
            %Plot cancer Min/Max Paths:
            plot3(p(minPath_cancer_seed).coord(1,:),p(minPath_cancer_seed).coord(2,:),p(minPath_cancer_seed).coord(3,:),'Color','m','linewidth',3);
            plot3(p(minPath_cancer_seed).coord(1, end),p(minPath_cancer_seed).coord(2, end),p(minPath_cancer_seed).coord(3, end),'ok','Markersize', 7,'markerfacecolor','w');

            plot3(p(maxPath_cancer_seed).coord(1,:),p(maxPath_cancer_seed).coord(2,:),p(maxPath_cancer_seed).coord(3,:),'y','linewidth',3);
            plot3(p(maxPath_cancer_seed).coord(1, end),p(maxPath_cancer_seed).coord(2, end),p(maxPath_cancer_seed).coord(3, end),'ok','Markersize', 7,'markerfacecolor','w');
            
            %Must specify curr_site_num for correct naming before calling savethis
            curr_site_num=i;
            savethis('minmaxpath','jpg');
        end
    end
end 

figure; 
plotenvironment(1);

%Plot overall Min/Max Paths:
plot3(p(minPath_seed).coord(1,:),p(minPath_seed).coord(2,:),p(minPath_seed).coord(3,:),'Color','m','linewidth',3);
plot3(p(minPath_seed).coord(1, end),p(minPath_seed).coord(2, end),p(minPath_seed).coord(3, end),'ok','Markersize', 7,'markerfacecolor','w');

plot3(p(maxPath_seed).coord(1,:),p(maxPath_seed).coord(2,:),p(maxPath_seed).coord(3,:),'Color','y','linewidth',3);
plot3(p(maxPath_seed).coord(1, end),p(maxPath_seed).coord(2, end),p(maxPath_seed).coord(3, end),'ok','Markersize', 7,'markerfacecolor','w');

savethis('minmaxpathOverall','jpg');




%% Entire procedure finished
disp( 'done' );

%% Functions
function load_3D()
    global eigen_map coh_map
    % OPEN AND LOAD EIGENVECTOR AND FRACTIONAL ANISOTROPY FILE
    EV_file = strcat('Anistropy_Files/CTRLP60_avg_vec.img');
    FA_file = strcat('Anistropy_Files/CTRLP60_avg_fa.img');
    fid = fopen(EV_file,'r');
    EV = fread(fid,'float32','ieee-le');
    EV = permute(reshape(EV,3,200,280,128),[2,3,4,1]);
    fclose(fid);

    % fractional anisotropy 
    fid = fopen(FA_file,'r');
    FA = fread(fid,'float32','ieee-le');
    FA = reshape(FA,200,280,128);
    fclose(fid);

    % INTERPOLATE TO THE APPROPRIATE SIZE
    limits = [21,180,21,160,36,105];    scale = 2;
    EV = EV(limits(1):limits(2),limits(3):limits(4),limits(5):limits(6),:);
    EV_x = interp3(EV(:,:,:,1),scale);  EV_y = interp3(EV(:,:,:,2),scale);  EV_z = interp3(EV(:,:,:,3),scale);
    clear EV; 
    EV(:,:,:,1) = EV_x; EV(:,:,:,2) = EV_y; EV(:,:,:,3) = EV_z;    
    clear EV_x EV_y EV_z

    FA = interp3(FA(limits(1):limits(2),limits(3):limits(4),limits(5):limits(6)),scale);
    
    eigen_map = EV;
    coh_map = FA;
end

function [seed_ind, has_cancer] = set_initial() 
    global coh_map modelType seed_sd cancer_center ubound lbound inj_center has_cancer
        
    switch modelType
    case 'Intranasal'
        inj_center =    [319, 100, 66];
    case 'Intracerebral'
        inj_center =    [450, 270, 170];
    end
    
    has_cancer = ~isempty(cancer_center);
    
    [Y,X,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));
    
    seed_sd = 10;
    seed_ROI = sqrt((X-inj_center(1)).^2 + (Y-inj_center(2)).^2 + (Z-inj_center(3)).^2)<seed_sd;   
    seed_ind = find(seed_ROI);
    
    %Create better upper and lower bounds (z-axis) for mouse brain
    for n = 1:size(coh_map,1)
      for m = 1:size(coh_map,2)
        ind = find( squeeze(coh_map(n,m,:)) ~=0 );
        if( isempty(ind) )
            indup(n,m) = NaN;
            inddown(n,m) = NaN;  
        else
            indup(n,m) = max(ind);
            inddown(n,m) = min(ind);
        end
      end
    end
    ubound = indup;
    lbound = inddown;
end

function [concentration, cgradX, cgradY, cgradZ] = set_cancer()
    global coh_map cancer_center cancer_size num_cancers curr_site_num
    [Y,X,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));

    % If multiple cancer sites exist, we assume they are all the same size
    num_cancers = size(cancer_center,1);
    curr_site_num = 0;

    concen_sd = [25 45 25]; 
    concen_param = 2; 

    % If multiple cancer sites exist, take average of concentrations
    concentration=0;
    for i=1:num_cancers
        concentration = concentration + 1./(1+(sqrt(((X - cancer_center(i,1))/concen_sd(1)).^2+((Y - cancer_center(i,2))/concen_sd(1)).^2+((Z - cancer_center(i,3))/concen_sd(1)).^2)).^concen_param);    
    end
    concentration = concentration/num_cancers;

    % valcap = concentration at border of cancer
    valcap = 1./(1+(sqrt(((cancer_size(1))/concen_sd(1)).^2+((0)/concen_sd(1)).^2+((0)/concen_sd(1)).^2)).^concen_param);      
    concentration(concentration > valcap) = valcap;
    concentration( coh_map==0 ) = 0; 

    cgradX = zeros(size(X));
    cgradX(2:end-1,:,:) = concentration(3:end,:,:) - concentration(1:end-2,:,:);

    cgradY = zeros(size(X));
    cgradY(:,2:end-1,:) = concentration(:,3:end,:) - concentration(:,1:end-2,:);

    cgradZ = zeros(size(X));
    cgradZ(:,:,2:end-1) = concentration(:,:,3:end) - concentration(:,:,1:end-2);

    % find largest vector magnitude
    L = (cgradX.^2 + cgradY.^2 + cgradZ.^2).^(.5);
    maxL = max(max(L));

    % normalize gradient
    cgradX = cgradX./maxL;
    cgradY = cgradY./maxL;
    cgradZ = cgradZ./maxL;   
end

function filename = savethis(casename,type)
    global modelType d_w d_g chemo_sensitivity alpha4chmtx beta4dist cancer_center has_cancer num_cancers curr_site_num;    
    
    if has_cancer
        if num_cancers > 1
            cancer_loc = strcat('[',num2str(num_cancers),']');
            %If curr_site_num is not specified before calling savethis,
            %include ALL cancer locations in title/naming
            if ~curr_site_num
                for i=1:num_cancers
                    cancer_loc = strcat('[',num2str(cancer_center(i,1)),',',...
                                            num2str(cancer_center(i,2)),',',...
                                            num2str(cancer_center(i,3)),']',cancer_loc);
                end
            else
                cancer_loc = strcat('[',num2str(cancer_center(curr_site_num,1)),',',...
                                        num2str(cancer_center(curr_site_num,2)),',',...
                                        num2str(cancer_center(curr_site_num,3)),']',cancer_loc);
            end
            curr_site_num=0;
        else
            cancer_loc = strcat('[',num2str(cancer_center(1)),',',...
                                    num2str(cancer_center(2)),',',...
                                    num2str(cancer_center(3)),']');
        end
        
        filename = strcat( '3D_', casename,'_',modelType, '_d',num2str(d_w), '_dg',num2str(d_g), ...
            '_a',num2str(alpha4chmtx), '_b',num2str(beta4dist), '_c',num2str(chemo_sensitivity), '_',cancer_loc, '.',type);
        titletext = strcat( casename,', ',modelType,', dw',num2str(d_w),', dg',num2str(d_g),', a',num2str(alpha4chmtx), ...
                ', b',num2str(beta4dist),', c',num2str(chemo_sensitivity),', ',cancer_loc);
    else
        cancer_loc = '[NA]';
        filename = strcat( '3D_', casename,'_',modelType, '_dw',num2str(d_w), '_dg',num2str(d_g), ...
            '_b',num2str(beta4dist), '_',cancer_loc, '.',type);
        titletext = strcat( casename,', ',modelType,', dw',num2str(d_w),', dg',num2str(d_g),', b',num2str(beta4dist),', ',cancer_loc);
    end
        
    if ~strcmp(type,'mat')
        title( titletext );
        saveas(gcf, filename);
    end
end

function plotenvironment(coherency_value)
    global coh_map ubound lbound seed_sd inj_center cancer_size cancer_center has_cancer num_cancers;
    [Y,X,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));
    
    hold on; surf( ubound', 'linestyle', 'none' , 'FaceColor', [0 0.4470 0.7410] , 'facealpha', 0.3 );      %Plot upper bound of mouse brain
    hold on; surf( lbound', 'linestyle', 'none' , 'FaceColor', [0 0.4470 0.7410] , 'facealpha', 0.3 );      %Plot lower bound of mouse brain

    [x,y,z] = sphere;
    x = x*seed_sd + inj_center(1);
    y = y*seed_sd + inj_center(2);
    z = z*seed_sd + inj_center(3);
    h = surf(x, y, z);                                                                                      %Plot injection site
    set(h,'FaceColor',[1 0 1],'FaceAlpha',0.4,'FaceLighting','gouraud','EdgeColor','none');

    if has_cancer
        for i=1:num_cancers
            [x,y,z] = sphere;
            x = x*cancer_size(1) + cancer_center(i,1);
            y = y*cancer_size(2) + cancer_center(i,2);
            z = z*cancer_size(3) + cancer_center(i,3);
            h = surf(x, y, z);                                                                              %Plot cancer site
            set(h,'FaceColor',[1 .7 0],'FaceAlpha',0.4,'FaceLighting','gouraud','EdgeColor','none');
        end
    end
    
    ii = find( coh_map > coherency_value ); 
    hold on; plot3( X(ii(1:100:end)),Y(ii(1:100:end)),Z(ii(1:100:end)),'b.');                               %Plot WM track
    
    xlabel('x'); ylabel('y'); zlabel('z'); grid on; axis equal; daspect([1 1 1]); camlight;
end
