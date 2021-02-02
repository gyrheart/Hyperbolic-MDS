
%% This script reproduces the hyperbolic MSD results in Fig.3 in the manuscript 
% Reproduce Figure 3 in the manuscript.

% Perform MDS for samples taken from a k-means cluster or randomly from
% whole batch.

% Here you can load the data and directly get the plots, you can also uncomment the
% codes and run the results yourself

% To load the E_MTAB_62.mat file, you need to download and transform the
% data first using the code in ../Saved_data/SaveLukkData.m
%%
clearvars
clc
addpath('../Hyperbolic_functions/')
addpath('../Saved_data/')
%% Load data
load E_MTAB_62 Data
%%

[r,c] = size(Data);
Rmax = 2.6;
dist_thres = 1e-4;
Min_ratio = 0;
D_Embed = 5;
options.MaxIter =1e3;
fit_model_power = 'power1';
fontsize = 12;
color = [0.6,0.6,0.6];
markersize = 2;
linewidth = 1.0;
%% Samples in a cluster
clc
r = 22283;
ClusterEmbedPower = zeros(4,2);
hFig = figure(2);
set(hFig,'units','centimeters','position',[0,0,15 12])
load HMDSGeneDataClusterIndex ClusterIndex % generated from k-means clustering
for count_dimension = 1:4
    switch count_dimension
    case 1
    D_Sample = 20;
    case 2
    D_Sample = 100;
    case 3
    D_Sample = 1000;
    case 4
    D_Sample = r;            
    end

%     sample_data = Data(:,ClusterIndex);
%     rand_index = randperm(r);
%     gene_index = rand_index(1:D_Sample);
%     dist1 = pdist(sample_data(gene_index,:)')/sqrt(D_Sample);
%     save(['GeneDataClusterSample',num2str(D_Sample),'DDataDistN100.mat'],'dist1');

    load(['GeneDataClusterSample',num2str(D_Sample),'DDataDistN100.mat'],'dist1');

%     [Y,stress, disparity] = mdscale(squareform(dist1),D_Embed,'start','random','Options',options);
%     dist2 = pdist(Y); 
%     save(['GeneDataClusterSample',num2str(D_Sample),'DEucEmbedN100.mat'],'Y');   

    load(['GeneDataClusterSample',num2str(D_Sample),'DEucEmbedN100.mat'],'Y');
     dist2 = pdist(Y); 

    subaxis(4,4,4*(count_dimension-1)+1,'Spacing',0.0,'Padding',0.005,'Margin', 0.00,'MarginTop',0.0);
    [sort_dist1,b] = sort(dist1);
     sort_dist2 = dist2(b);
     x = sort_dist1';
     x = x - min(x)+dist_thres;

     y = sort_dist2';
    f = fit(x,y,fit_model_power);
    plot(x,y,'k.','color',color,'markersize',markersize)
    hold on
    plot(x,f(x),'b-','linewidth',linewidth)
    hold off
    xrange0 = 0;
    yrange0 = 0;
    xrange = 0.8;
    yrange = 1.1;
    axis([xrange0 xrange0+xrange yrange0 yrange0+yrange]);
    set(gca,'fontsize',15,'xtick',[],'ytick',[])

    power = f.b;
    ClusterEmbedPower(count_dimension,1) = power;
    text(xrange0+0.05*xrange,yrange0+0.75*yrange,['\kappa = ',num2str(round(power,2)-1)],'fontsize',fontsize)
    text(xrange0+0.05*xrange,yrange0+0.9*yrange,[num2str(D_Sample),' Genes'],'fontsize',fontsize)

    % hyperbolic mds 
%     iters = 0;
%     count = 0;
%     max_iter = 0;
%     while iters < options.MaxIter
%     
%        try
%             [Y_hyper,stress, disparity,iters] = fmdscale_hyperbolic(squareform(dist1),D_Embed,Rmax,Min_ratio,'Options',options);
%        catch
%            try
%             [Y_hyper,stress, disparity,iters] = fmdscale_hyperbolic(squareform(dist1),D_Embed,Rmax,Min_ratio,'Options',options);
%            catch
%                warning('errors')
%                 break
%            end
%        end
%         iters
%        if iters > max_iter
%             opt_pos = Y_hyper;
%             max_iter = iters; 
%         end
%         count = count + 1;
%         if count > 3e1
%             Y_hyper = opt_pos;
%             break
%         end
%     end
%     Y = fRescaleAngleGeneral(Y_hyper);
%     dist2 = squareform(fDistNative(Y));
%     save(['GeneDataClusterSample',num2str(D_Sample),'DFinalR2_6HyperEmbedN100.mat'],'Y');    

%     % 
    load(['GeneDataClusterSample',num2str(D_Sample),'DFinalR2_6HyperEmbedN100.mat'],'Y');    
    dist2 = squareform(fDistNative(Y));

    subaxis(4,4,4*(count_dimension-1)+2,'Spacing',0.0,'Padding',0.005,'Margin', 0.00,'MarginTop',0.0);
    [sort_dist1,b] = sort(dist1);
    sort_dist2 = dist2(b);
    x = sort_dist1';  
    x = x - min(x)+dist_thres;
    y = sort_dist2';
    f = fit(x,y,'power1');
    power = f.b;
    plot(x,y,'.','color',color,'markersize',markersize)
    hold on
    plot(x,f(x),'r-','linewidth',linewidth)
    hold off
    xrange0 = 0.0;
    yrange0 = 0.0;
    xrange = 0.8;
    yrange = 7.5;
    axis([xrange0 xrange0+xrange yrange0 yrange0+yrange]);
    set(gca,'fontsize',15,'xtick',[],'ytick',[])
    text(xrange0+0.05*xrange,yrange0+0.75*yrange,['\kappa = ',num2str(round(power,2)-1)],'fontsize',fontsize)
    text(xrange0+0.05*xrange,yrange0+0.9*yrange,[num2str(D_Sample),' Genes'],'fontsize',fontsize)
end

%% Random samples
[r,c] = size(Data);
R = 2.6;
Min_ratio = 0;
D_Embed = 5;
options.MaxIter =1e3;
RandomEmbedPower = zeros(4,2);
load HMDSGeneDataRepeatRandomIndex RandomIndex
for count_dimension = 1:4
    switch count_dimension
    case 1
    D_Sample = 20;
    case 2
    D_Sample = 100;
    case 3
    D_Sample = 1000;
    case 4
    D_Sample = r;            
    end
    
%     sample_data = Data(:,RandomIndex);
%     rand_index = randperm(r);
%     gene_index = rand_index(1:D_Sample);
%     dist1 = pdist(sample_data(gene_index,:)')/sqrt(D_Sample);
%     save(['GeneDataRandomSample',num2str(D_Sample),'DDataDistN100.mat'],'dist1');

    load(['GeneDataRandomSample',num2str(D_Sample),'DDataDistN100.mat'],'dist1');

%     [Y,stress, disparity] = mdscale(squareform(dist1),D_Embed,'start','random','Options',options);
%     dist2 = pdist(Y);
%     save(['GeneDataRandomSample',num2str(D_Sample),'DEucEmbedN100.mat'],'Y');   
%     
    load(['GeneDataRandomSample',num2str(D_Sample),'DEucEmbedN100.mat'],'Y');
    dist2 = pdist(Y); 

    subaxis(4,4,4*(count_dimension-1)+3,'Spacing',0.0,'Padding',0.005,'Margin', 0.00,'MarginTop',0.0);
    [sort_dist1,b] = sort(dist1);
    sort_dist2 = dist2(b);
    x = sort_dist1';
    x = x - min(x)+dist_thres;
    y = sort_dist2';
    f = fit(x,y,fit_model_power);
    plot(x,y,'.','color',color,'markersize',markersize)
    hold on
    plot(x,f(x),'b-','linewidth',linewidth)
    hold off
    xrange0 = 0;
    yrange0 = 0;
    xrange = 2.2;
    yrange = 2.4;
    axis([xrange0 xrange0+xrange yrange0 yrange0+yrange]);
    set(gca,'fontsize',15,'xtick',[],'ytick',[])
    power = f.b;
    RandomEmbedPower(count_dimension,1) = power;
    text(xrange0+0.05*xrange,yrange0+0.75*yrange,['\kappa = ',num2str(round(power,2)-1)],'fontsize',fontsize)
    text(xrange0+0.05*xrange,yrange0+0.9*yrange,[num2str(D_Sample),' Genes'],'fontsize',fontsize)

%     iters = 0;
%     count = 0;
%     max_iter = 0;
%     while iters < options.MaxIter
%     iters
%        try
%             [Y_hyper,stress, disparity,iters] = fmdscale_hyperbolic(squareform(dist1),D_Embed,Rmax,Min_ratio,'Options',options);
%        catch
%            try
%             [Y_hyper,stress, disparity,iters] = fmdscale_hyperbolic(squareform(dist1),D_Embed,Rmax,Min_ratio,'Options',options);
%            catch
%                warning('errors')
%                 break
%            end
%        end
% 
%        if iters > max_iter
%             opt_pos = Y_hyper;
%             max_iter = iters; 
%         end
%         count = count + 1;
%         if count > 3e1
%             Y_hyper = opt_pos;
%             break
%         end
%     end
%     Y = fRescaleAngleGeneral(Y_hyper);
%     dist2 = squareform(fDistNative(Y));
%     save(['GeneDataRandomSample',num2str(D_Sample),'DFinalR2_6HyperEmbedN100.mat'],'Y');    

    load(['GeneDataRandomSample',num2str(D_Sample),'DFinalR2_6HyperEmbedN100.mat'],'Y');    
    dist2 = squareform(fDistNative(Y));

    subaxis(4,4,4*(count_dimension-1)+4,'Spacing',0.0,'Padding',0.005,'Margin', 0.00,'MarginTop',0.00);
    [sort_dist1,b] = sort(dist1);
    sort_dist2 = dist2(b);
    x = sort_dist1';
    x = x - min(x)+dist_thres;
    y = sort_dist2';
    f = fit(x,y,fit_model_power);
    power = f.b;
    plot(x,y,'.','color',color,'markersize',markersize)
    hold on
    plot(x,f(x),'r-','linewidth',linewidth)
    hold off
    xrange0 = 0;
    yrange0 = 0;
    xrange = 2;
    yrange = 8;
    axis([xrange0 xrange0+xrange yrange0 yrange0+yrange]);
    set(gca,'fontsize',15,'xtick',[],'ytick',[])
    text(xrange0+0.05*xrange,yrange0+0.75*yrange,['\kappa = ',num2str(round(power,2)-1)],'fontsize',fontsize)
    text(xrange0+0.05*xrange,yrange0+0.9*yrange,[num2str(D_Sample),' Genes'],'fontsize',fontsize)
end
set(gcf, 'Renderer', 'Painters');
