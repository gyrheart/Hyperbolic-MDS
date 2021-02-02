%% Test hyperbolic MDS on synthetic hyperbolic data with known dimension and radius

% Here you can load the data and get the plots immediately, you can also uncomment the
% codes and run the results yourself


clearvars
clc
addpath('../Hyperbolic_functions/')
addpath('../Saved_data/')

%% Generate synthetic hyperbolic data

d = 5;
N = 50;
R = 3;
MinRatio = 0;
% sample_data = fSamplingHyperbolicUnif(N,d,R,MinRatio); % sampling the data
% save(['HMDSSynDataHyperToHyperN50D5R3.mat'],'sample_data')
%% Do hyperbolic MDS embedding for the data

load(['HMDSSynDataHyperToHyperN50D5R3.mat'],'sample_data')

% d = 5;
% R = 3;
% MinRatio = 0;
% options.MaxIter  = 2e3;
% iter_thres = 1e3;
% Dist = fDistNative(sample_data); % hyperbolic distances
% iters = 0;
% count = 0;
% max_iter = 0;
% while iters < iter_thres
% 
%    try
%         [Y_hyper,stress, disparity,iters] = fmdscale_hyperbolic(Dist,d,R,MinRatio,'UpdateRadius',true,'Options',options);
%    catch
%        try
%         [Y_hyper,stress, disparity,iters] = fmdscale_hyperbolic(Dist,d,R,MinRatio,'UpdateRadius',true,'Options',options);
%        catch
%            warning('errors')
%             break
%        end
%    end
%    disp(iters)
%    if iters > max_iter
%         opt_pos = Y_hyper;
%         max_iter = iters; 
%     end
%     count = count + 1;
%     if count > 1e1
%         Y_hyper = opt_pos;
%         break
%     end
% end
% embed_data = Y_hyper;
% 
% save(['HMDSSynEmbedHyperToHyperN50D5_R3to3.mat'],'embed_data')


%% Compare the distances and radii distribution of embedding points with the data

load(['HMDSSynDataHyperToHyperN50D5R3.mat'],'sample_data')
load(['HMDSSynEmbedHyperToHyperN50D5_R3to3.mat'],'embed_data')


Dist_data = squareform(fDistNative(sample_data));
R_data = sample_data(:,1);
Dist_model = squareform(fDistNative(embed_data));
R_model = embed_data(:,1);
nbins = 15;
fontsize = 12;
hFig = figure(1);
set(hFig,'units','centimeters','position',[0,0,18,8]);
subplot(1,2,1) % plot embedding distances vesus data distances
plot(Dist_data,Dist_model,'.','color', [0.5 0.5 0.5])
f = fit(Dist_data',Dist_model','power1');
hold on
plot(sort(Dist_data),f(sort(Dist_data)),'r-','linewidth',1)
hold off
xrange0 = 1;
yrange0 = 1;
xrange = 5.5;
yrange = 5.5;
axis([xrange0,xrange0+xrange,yrange0,yrange0+yrange])
xlabel('Data dist.')
ylabel('Embedding dist. ')
set(gca,'fontsize',10,'xtick',1:6,'ytick',1:6,'ticklength',[0.04,0.04])
title('Pairwise distances plot')
text(xrange0+0.1*xrange,yrange0+0.9*yrange,'R_{data} = 3, R_{mds} = 3','fontsize',fontsize)

subplot(1,2,2) % plot embedding radii versus data radii
histogram(R_data,nbins)
hold on
histogram(R_model,nbins)
hold off
xlim([0 3.0])
ylim([0 30])
xlabel('r')
ylabel('Counts')
set(gca,'fontsize',10,'xtick',1:6,'ytick',0:10:40,'ticklength',[0.04,0.04])
set(gcf,'renderer','painters')
title('Radii distributions')
legend('data', 'embedding')
