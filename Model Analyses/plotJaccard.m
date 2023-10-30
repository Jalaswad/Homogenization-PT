clear; close all;
load('PTr_Jaccard_output_shelf_40prctfossil.mat')
load('permian_grid.mat')

% distance bins (km)
dbin = 0e3:2e3:2e4;

% mask out same grid cell (Jaccard = 1)
for yy = 1:116
    for xx = 1:100

        % permian
        Pc.J(yy,xx,yy,xx)=nan;
        % triassic
        Tc.J(yy,xx,yy,xx)=nan;
        
    end
end

% open space 
Pc.J_d = zeros(116,100,length(dbin))*nan;
Tc.J_d = zeros(116,100,length(dbin))*nan;

% Compute mean Jaccard vs distance relationship for each shelf location
for xx = 1:100
    for yy = 1:116
        
        % Jaccard for this grid cell
        this_cell = squeeze(Pc.J(yy,xx,:,:));
        
        % distances to this cell (km)
        this_dist = squeeze(par.distance(yy,xx,:,:));

        % land mask & grid cell area
        Pc.mask = this_cell*nan;Pc.mask(~isnan(this_cell))=1;
        Pc.V = grid.A*nan;Pc.V(Pc.mask>0)=grid.A(Pc.mask>0);
        Tc.mask = this_cell*nan;Tc.mask(~isnan(this_cell))=1;
        Tc.V = grid.A*nan;Tc.V(Tc.mask>0)=grid.A(Tc.mask>0);

        % bin into distance intervals
        for dd = 1:length(dbin)-1
            
            this_bin = dbin(dd);
            if dd == length(dbin)-1
                this_bin_up = dbin(dd+1)+20;% inclusive
            else
                this_bin_up = dbin(dd+1);
            end
            % indices of grid cells @ distance dd
            idx = find(this_dist >= this_bin & this_dist < this_bin_up);

            % average 
            Pc.J_d(yy,xx,dd) = sum(this_cell(idx).*Pc.V(idx),'omitnan')./sum(Pc.V(idx),'omitnan'); % area-weighted 
            
        end

        % Jaccard for this grid cell
        this_cell = squeeze(Tc.J(yy,xx,:,:));
        
        % distances to this cell (km)
        this_dist = squeeze(par.distance(yy,xx,:,:));

        % bin into distance intervals
        for dd = 1:length(dbin)-1
            
            this_bin = dbin(dd);
            if dd == length(dbin)-1
                this_bin_up = dbin(dd+1)+20;% inclusive
            else
                this_bin_up = dbin(dd+1);
            end
            % indices of grid cells @ distance dd
            idx = find(this_dist >= this_bin & this_dist < this_bin_up);

            % average 
            Tc.J_d(yy,xx,dd) = sum(this_cell(idx).*Tc.V(idx),'omitnan')./sum(Tc.V(idx),'omitnan');% area-weighted 
           
        end

    end
end

% Global (shelf) mean relationship 
for dd = 1:length(dbin)-1

    % land mask & grid cell area
    this_cell =squeeze(Pc.J_d(:,:,dd)); 
    Pc.mask = this_cell*nan;Pc.mask(~isnan(this_cell))=1;
    Pc.V = grid.A*nan;Pc.V(Pc.mask>0)=grid.A(Pc.mask>0);
    this_cell =squeeze(Tc.J_d(:,:,dd)); 
    Tc.mask = this_cell*nan;Tc.mask(~isnan(this_cell))=1;
    Tc.V = grid.A*nan;Tc.V(Tc.mask>0)=grid.A(Tc.mask>0);

       
    % Permian
    this_d = Pc.J_d(:,:,dd);
    Pc.Jd_mn(dd)= sum(this_d(:).*Pc.V(:),'omitnan')./sum(Pc.V(:),'omitnan');%weighted ave
    Pc.Jd_sd(dd)= std(this_d(:),'omitnan'); % std 
    
    % Triassic
    this_d = Tc.J_d(:,:,dd); 
    Tc.Jd_mn(dd)= sum(this_d(:).*Tc.V(:),'omitnan')./sum(Tc.V(:),'omitnan');% weighted
    Tc.Jd_sd(dd)= std(this_d(:),'omitnan'); %std 
    
end


%% Plot similarity vs distance 
figure
dbin = dbin(2:end);
c1=[.1 .8 1];
c2=[0 .2 .5];

% permian 
plot(dbin,Pc.Jd_mn,'color',c1,'linewidth',3)
hold on 
plot(dbin,Pc.Jd_mn,'color',c2,'linewidth',3)
patch_x = [dbin(1:end-2) fliplr(dbin(1:end-2))];
patch_y = [Pc.Jd_mn(1:end-2)-Pc.Jd_sd(1:end-2) fliplr(Pc.Jd_mn(1:end-2)+Pc.Jd_sd(1:end-2))];
p=patch(patch_x,patch_y,'b');set(p,'facecolor',c1,'edgecolor',[0.1 .1 0.8],'facealpha',.1,'edgealpha',0);hold on
plot(dbin,Pc.Jd_mn,'color',c1,'linewidth',3)

% triassic 
patch_x = [dbin(1:end-2) fliplr(dbin(1:end-2))];
patch_y = [Tc.Jd_mn(1:end-2)-Tc.Jd_sd(1:end-2) fliplr(Tc.Jd_mn(1:end-2)+Tc.Jd_sd(1:end-2))];
p=patch(patch_x,patch_y,'b');set(p,'facecolor',c2,'edgecolor',[0.8 .6 0.1],'facealpha',.3,'edgealpha',0);hold on
hold on
plot(dbin,Tc.Jd_mn,'color',c2,'linewidth',3)
xlim([0 1.25e4])
ylim([0 1])
xlabel('Distance (km)')
ylabel('Similarity (Jaccard)')
set(gca,'fontsize',18)
box on
xlim([1.5e3 13e3])

% observations 
load('PT_J_fossil_obs.mat')
scatter(Pobs.x(1:end-3),Tobs.J(1:end-3),Tobs.n(1:end-3)*10,c2,'filled')
hold on
scatter(Pobs.x(1:end-3),Pobs.J(1:end-3),Pobs.n(1:end-3)*10,c1,'filled')
errorbar(Pobs.x(1:end-3),Pobs.J(1:end-3),Pobs.Jsd(1:end-3),'color',c1)
errorbar(Pobs.x(1:end-3),Tobs.J(1:end-3),Tobs.Jsd(1:end-3),'color',c2)
plot(Pobs.x(1:end),Tobs.J(1:end),'color',c2,'linewidth',1)
plot(Pobs.x(1:end),Pobs.J(1:end),'color',c1,'linewidth',1)
legend('Before','After')
