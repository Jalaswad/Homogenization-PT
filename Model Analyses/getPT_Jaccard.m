% Script to compute Jaccard index from CESM P/Tr climate simulations
clear all;
close all; 

% Jaccard index 
% J(x,y) = m / (a + b + m )
% where m = #  genera in both x,y
% a = # unique genera in cell x 
% b = # unique genera in cell y 

% fname 
load('PTr_Jaccard_input.mat')

% LAT LONG
load LATLON.mat

% Earth's Geometry for distance calculation
S = referenceSphere('earth');

% shelf depth
zshelf = 20;

% extinction threshold
Vcrit = 50;

% sampling fraction for fossil record
fbias = 0.4;

% migrate (unused)
migrate=1; 

% Benthic Masks
% Bottom grid cells
% fname 
load permian_grid.mat

Bot = grid.M3d*NaN;
for zz = 1:60
    for ii = 1:116
        for jj = 1:100
        if grid.M3d(zz,ii,jj) ==1 && grid.M3d(zz+1,ii,jj) ~= 1
        Bot(zz,ii,jj) = 1;
        end
        end
    end
end

% mask for shelf regions (1 = ocean above shelf, 0 = land or deeper)
par.shelf = squeeze(sum(Bot(1:zshelf,:,:),'omitnan'));
par.shelf(par.shelf<1)=nan;

% Extinction index (0 = extinct; 1 = survivor)
T.em= T.delVm<Vcrit; % habitat loss below extinction threshold?

%% Species weight * initial habitat in Permian
for aa = 1:12
    for bb = 1:12
        for cc =1:13
            for dd = 1:11
                P.w(aa,bb,cc,dd) = par.Tmw(aa).*par.Pcw(bb).*par.Eow(cc).*par.Aow(dd).*(squeeze(P.V(aa,bb,cc,dd))>0);
            end
        end
    end
end

%% 3D to 2D habitat  
% permian, triassic w/ and w/o extinction (X)
for aa = 1:12
    for bb = 1:12
        for cc =1:13
            for dd = 1:11
                % presence/absence above shelf depth?    
                T.msk2d(aa,bb,cc,dd,:,:) = squeeze(max(T.msk(aa,bb,cc,dd,1:zshelf,:,:),[],5)>0);
                
                % triassic (survivors only)
                T.msk2dX(aa,bb,cc,dd,:,:) = (squeeze(max(T.msk(aa,bb,cc,dd,1:zshelf,:,:),[],5)>0).*T.em(1,aa,bb,cc,dd)>0);

                % permian (survivors only or survivors + victims)
                P.msk2d(aa,bb,cc,dd,:,:) = (squeeze(max(P.msk(aa,bb,cc,dd,1:zshelf,:,:),[],5)>0));
                
            end
        end
    end
end

%% fossil record random sampling scheme
% open space 
Rmat = zeros(12,12,13,11,116,100);

% Permian
for yy = 1:116
for xx = 1:100
    % random sample 
        sizeR = [12    12    13    11];
        num1 = fbias*sizeR(1)*sizeR(2)*sizeR(3)*sizeR(4);% number of fossilized species
        num1=round(num1); % to ensure integer number of indices
        R = zeros(sizeR);  % set all to zero
        ix = randperm(numel(R)); % randomize the linear indices
        ix = ix(1:num1); % select the first set for fossilization 
        R(ix) = 1; % set the corresponding indices to 1
        R = R>0; % make logical
        Rmat(:,:,:,:,yy,xx)=R;% sampled ecophysiotypes for this grid cell
end
end
Pc.Rmat = Rmat;
clear Rmat

% open space 
Rmat = zeros(12,12,13,11,116,100);
% Triassic
for yy = 1:116
for xx = 1:100
    % random sample 
        sizeR = [12    12    13    11];
        num1 = fbias*sizeR(1)*sizeR(2)*sizeR(3)*sizeR(4); % number of fossilized species
        num1=round(num1); % to ensure interger number of indices
        R = zeros(sizeR);  % set all to zero
        ix = randperm(numel(R)); % randomize the linear indices
        ix = ix(1:num1); % select the first set for fossilization
        R(ix) = 1; % set the corresponding indices to 1
        R = R>0; % make logical
        Rmat(:,:,:,:,yy,xx)=R;
end
end
Tc.Rmat = Rmat;
clear Rmat


%% Jaccard index  
tic
% open space for this cell
Pc.J = zeros(116,100,116,100);
Tc.J=Pc.J;
par.distance = Pc.J;

% this grid cell
for xx = 1:100
    for yy = 1:116

       
        % compare this grid cell to all others
        for jj = 1:100
            for ii = 1:116
                
                % compute similarity between shelf regions
                if ~isnan(par.shelf(yy,xx)) && ~isnan(par.shelf(ii,jj))
                    % Permian
                    % species in this grid cell
                    sp_this = squeeze(P.msk2d(:,:,:,:,yy,xx));

                    % randomly sample species
                    sp_this = sp_this>0 & squeeze(Pc.Rmat(:,:,:,:,yy,xx))>0;
        
                    % species in other grid cell
                    sp_other = squeeze(P.msk2d(:,:,:,:,ii,jj));
                    
                    % randomly sample species
                    sp_other = sp_other>0 & squeeze(Pc.Rmat(:,:,:,:,ii,jj))>0;
                    
                    % unique species in this grid cell
                    a = P.w(sp_this>0 & sp_other < 1);% species weigthts
                    a = sum(a(:));

                    % unique species in other grid cell
                    b = P.w(sp_this<1 & sp_other > 0);% species weigthts
                    b = sum(b(:));

                    % species in both
                    m = P.w(sp_this>0 & sp_other > 0);% species weigthts
                    m = sum(m(:));

                    % Jaccard
                    Pc.J(yy,xx,ii,jj) = m./(a+b+m);

                    % Triassic

                    % species in this grid cell
                    sp_this = squeeze(T.msk2dX(:,:,:,:,yy,xx));

                    % randomly sample species
                    sp_this = sp_this>0 & squeeze(Tc.Rmat(:,:,:,:,yy,xx))>0;
        
                    % species in other grid cells
                    sp_other = squeeze(T.msk2dX(:,:,:,:,ii,jj));

                    % randomly sample species
                    sp_other = sp_other>0 & squeeze(Tc.Rmat(:,:,:,:,ii,jj))>0;
      
                    % unique species in this grid cell
                    a = P.w(sp_this>0 & sp_other < 1);% species weigthts
                    a = sum(a(:));

                    % unique species in other grid cell
                    b = P.w(sp_this<1 & sp_other > 0);% species weigthts
                    b = sum(b(:));

                    % species in both
                    m = P.w(sp_this>0 & sp_other > 0);% species weigthts
                    m = sum(m(:));

                    % Jaccard
                    Tc.J(yy,xx,ii,jj) = m./(a+b+m);

                    % Distance (km); distance(lat1,lon1,lat2,lon2,S)
                    % this lat, lon
                    lat1 = TLAT(yy,xx);
                    lon1 = TLONG(yy,xx);

                    % that lat, lon
                    lat2 = TLAT(ii,jj);
                    lon2 = TLONG(ii,jj);

                    % save distances (km)
                    par.distance(yy,xx,ii,jj)= distance(lat1,lon1,lat2,lon2,S).*1e-3;
                else
                    % non-shelf regions
                    Pc.J(yy,xx,ii,jj) = nan;
                    Tc.J(yy,xx,ii,jj) = nan;
                    par.distance(yy,xx,ii,jj)= nan;
                end

            end
        end
    end
end

% output 
par.fbias=fbias;
par.runtime=toc;
par.Vcrit = Vcrit;
par.zshelf =zshelf;
par.migrate=migrate;

% output
save('PTr_Jaccard_shelf_40prctfossil.mat','Tc','Pc','par','-v7.3')
