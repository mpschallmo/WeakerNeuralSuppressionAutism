function model = Normalization_Models(which_version)
% Usage: model = Normalization_Models(which_version)
%
% Inputs: which_version - integer, acceptable values are:
%   1 - Weaker normalization model, depicted in Figure 6A, B, & C
%   2 - Larger excitatory spatial filter model, from Figure 6D, E, & F
%   3 - Narrower spatial top-down modulation model, from Figure 6G, H, & I
%   4 - Extra-narrow top-down modulation, from Supp. Figure 3A, B, & C
%   5 - Ultra-narrow top-down modulation, from Supp. Figure 3D, E, & F
%
% Outputs: model - structure, with fields:
%   p = model parameters
%   stim = stimulus parameters
%   resp = model peak response (arbitrary units), used to calculate 
%          duration thresholds
%   thresh = model duration threshold
%   SI = model size indices
%   variable_param = structure, with field = parameter that varies in the
%                    chosen model version
%   which_version = which model version (input, as above)
%   dimension_key = key for dimension size of resp & thresh
%
% Authors: Michael-Paul Schallmo, Scott O. Murray
% University of Washington, Dept. of Psychology
% c. June 2019
%
% This code implements five versions of a normalization model to describe
% spatial suppression during motion direction discrimination.
% Specifically, the model variants are those described in the manuscript 
% entitled "Weaker neural suppression in autism" by Schallmo et al., 
% available at: https://www.biorxiv.org/content/10.1101/645846v1
%
% Uses code by Geoff M. Boynton 
% University of Washington, Dept. of Psychology
% - Specifically: makeNeuralImage.m & convolveImage.m
% - These are included in the distribution, along with this function.
% - These functions will need to be in your path to run this function.
%
% Adapted from code by Reynolds & Heeger (2009), from their paper entitled
% "The Normalization Model of Attention" in Neuron, 61(2), p. 168-185 
% Their original code is available from:
% http://www.snl-r.salk.edu/~reynolds/Normalization_Model_of_Attention/
% http://www.cns.nyu.edu/heegerlab/?page=software

%% check input
if ~exist('which_version','var')
    which_version = input('Which model version? (1-5): ');
end
if ~isnumeric(which_version) || isempty(which_version)
    error('which_version must be an integer 1-5');
elseif which_version ~= 1 && which_version ~= 2 && which_version ~= 3 && ...
        which_version ~= 4 && which_version ~= 5
    error('which_version must be an integer 1-5');
end

%% set up model parameters

p.x = -200:2:200;                   % list of stimulus positions (arbitrary coordinates)
p.theta = -180:2:180;               % list of orientations (degrees)

p.use_e_width = [4 4                % excitatory SF width; x_w_e,
                 5 4                %   5 different model versions x 2 groups (ASD, NT)
                 4 4
                 4 4
                 4 4];
p.e.theta.width = 25;               % excitatory orientation tuning width; theta_w_e

p.s.x.width = 40;                   % suppressive spatial pooling width; x_w_s
p.s.theta.width = 25;               % suppressive orientation pooling width; theta_w_s
p.S_gain = [0.75 1                  % suppressive gain, scaling factor; S_g,
            1 1                     %   5 different model versions x 2 groups (ASD, NT)
            1 1
            1 1
            1 1];          

p.sigma = .0002;                    % semi-saturation constant; sigma
  
stim.x.center = 0;                  % stimulus position
stim.all_xWidths = [1 2 12];        % stimulus spatial widths (multiple modeled), arbitrary units
stim.theta.center = -90;            % stimulus orientation 
stim.theta.width = 3;               % stimulus orientation width
stim.all_contrasts = [.03 .98];     % stimulus contrasts (multiple modeled)

topdown.x.center = 0;               % the center of spatial top-down component
p.MW = [14 14                       % the width of spatial top-down component,
        14 14                       %   5 different model versions x 2 groups (ASD, NT)
        6  14
        2  6
        1  6];
topdown.theta.center = -90;         % the center of "feature-based" top-down component
topdown.theta.width = inf;          % the width of "feature-based" top-down component
topdown.range = [1,4];              % sets the top-down gain (min & max); M_g
topdown.method = 'multiply';        % method for combining x & theta when modeling top-down gain

p.response_region_width = 1;        % e.g., 0 = peak response only, 2 = +/-2 from peak, etc.; r_w

p.criterion = 600;                  % arbitrary criterion to convert response to duration threshold

%% generate model pieces

% model response, pre-allocate
resp = nan(length(stim.all_xWidths), size(p.use_e_width,2), ...
    length(stim.all_contrasts));

for k = 1:length(stim.all_contrasts)
    stim.contrast = stim.all_contrasts(k);
    
    for j = 1:size(p.use_e_width,2) % 2 groups
        
        % set excitatory SF width, varies in which_version = 2
        p.e.x.width = p.use_e_width(which_version, j);
        
        % set width of spatial top-down component, varies in which_version = 3, 4, 5
        topdown.x.width = p.MW(which_version, j);
        
        for i = 1:length(stim.all_xWidths)
            stim.x.width = stim.all_xWidths(i);
            
            %(Stim) Stimulus Image
            Stim = makeNeuralImage(p, stim);
            
            %(E) Convolve stimulus with SF to get 'stimulus drive'
            E = convolveImage(p, Stim, p.e.x.width, p.e.theta.width);
            
            %(M) Top-down Image
            M = makeNeuralImage(p, topdown);
            
            % Top-down Gain
            E = E.*M;
            
            %(S) Convolve excitatory drive by suppressive filters for 'suppressive drive'
            S = convolveImage(p, E, p.s.x.width, p.s.theta.width);
            
            % apply suppressive gain, varies in which_version = 1
            S = S*p.S_gain(which_version, j); 
            
            %(R) Pass total excitation through normalization stage to get the population response
            R = (E)./(S + p.sigma);
            
            % find max response index
            [~, max_idx] = max(R(:));
            [max_a, max_b] = ind2sub(size(R), max_idx);
            
            % find peak response region
            get_R_region = R(max_a-p.response_region_width:max_a+p.response_region_width, ...
                max_b-p.response_region_width:max_b+p.response_region_width);
            
            % record mean response within peak region
            resp(i, j, k) = mean(get_R_region(:));
            
        end
    end
end

% convert from response to duration threshold
thresh = p.criterion/resp;

%% plot thresholds

% size of small, medium, and big grating stimuli from psychophysics
use_stim_sizes = [1 2 12];

% in case we model other sizes too, pick these out
small_x_idx = stim.all_xWidths == use_stim_sizes(1);
med_x_idx = stim.all_xWidths == use_stim_sizes(2);
big_x_idx = stim.all_xWidths == use_stim_sizes(3);

% 2 groups
asd_model_idx = 1;
nt_model_idx = 2;

% 2 stimulus contrasts; 3% and 98% grating stimuli from psychophysics
low_cont_idx = 1;
high_cont_idx = 2;

use_axis = [0.6 20 20 200];

model_legends = {'Weak normalization','Strong normalization';
                 'Large excitatory SF','Small excitatory SF';
                 'Narrow top-down','Broad top-down';
                 ['M width = ' num2str(p.MW(which_version, 1))], ...
                    ['M width = ' num2str(p.MW(which_version, 2))];
                 ['M width = ' num2str(p.MW(which_version, 1))], ...
                    ['M width = ' num2str(p.MW(which_version, 2))]};

% low contrast thresholds
figure
subplot(1,2,2)
hold on
plot(use_stim_sizes, [thresh(small_x_idx, asd_model_idx, low_cont_idx) ...
    thresh(med_x_idx, asd_model_idx, low_cont_idx) ...
    thresh(big_x_idx, asd_model_idx, low_cont_idx)], ...
    'b-s','MarkerSize', 8, 'MarkerFaceColor','b','linewidth',1)
plot(use_stim_sizes, [thresh(small_x_idx, nt_model_idx, low_cont_idx) ...
    thresh(med_x_idx, nt_model_idx, low_cont_idx) ...
    thresh(big_x_idx, nt_model_idx, low_cont_idx)], ...
    'b-s','MarkerSize', 8, 'MarkerFaceColor','w','linewidth',1)
set(gca,'XScale','log','XTick',[1 10],'Xticklabel',{'1','10'},...
    'YScale','log','YTick',[10 25 50 100 200],'fontsize',18,...
    'XColor','k','YColor','k')
axis(use_axis)
ax = gca; ax.XColor = [0 0 0]; ax.YColor = [0 0 0]; ax.FontName = 'Arial';
title('3% Contrast','color','k')
box off

% high contrast thresholds
subplot(1,2,1)
hold on
plot(use_stim_sizes, [thresh(small_x_idx, asd_model_idx, high_cont_idx) ...
    thresh(med_x_idx, asd_model_idx, high_cont_idx) ...
    thresh(big_x_idx, asd_model_idx, high_cont_idx)], ...
    'r-s','MarkerSize', 8, 'MarkerFaceColor','r','linewidth',1)
plot(use_stim_sizes, [thresh(small_x_idx, nt_model_idx, high_cont_idx) ...
    thresh(med_x_idx, nt_model_idx, high_cont_idx) ...
    thresh(big_x_idx, nt_model_idx, high_cont_idx)], ...
    'r-s','MarkerSize', 8, 'MarkerFaceColor','w','linewidth',1)
set(gca,'XScale','log','XTick',[1 10],'Xticklabel',{'1','10'},...
    'YScale','log','YTick',[10 25 50 100 200],'fontsize',18,...
    'XColor','k','YColor','k')
axis(use_axis)
ax = gca; ax.XColor = [0 0 0]; ax.YColor = [0 0 0]; ax.FontName = 'Arial';
box off
set(gcf,'color','w','POS',[10   150   670   400])
xlabel('Stim. size (arb. units)','color','k')
ylabel('Threshold (arb. units)','color','k')
h = legend(model_legends{which_version,1}, model_legends{which_version,2});
set(h,'TextColor','k')
title('98% Contrast','color','k')

%% plot size indices (SIs)

SI.low = [log10(thresh(med_x_idx, asd_model_idx, low_cont_idx)) - ...
    log10(thresh(big_x_idx, asd_model_idx, low_cont_idx)); ...
    log10(thresh(med_x_idx, nt_model_idx, low_cont_idx)) - ...
    log10(thresh(big_x_idx, nt_model_idx, low_cont_idx))];
SI.high = [log10(thresh(med_x_idx, asd_model_idx, high_cont_idx)) - ...
    log10(thresh(big_x_idx, asd_model_idx, high_cont_idx)); ...
    log10(thresh(med_x_idx, nt_model_idx, high_cont_idx)) - ...
    log10(thresh(big_x_idx, nt_model_idx, high_cont_idx))];

offset_x = 0.075;

figure
hold on
plot([0 5],[0 0],'k--','linewidth',1)
plot(2+offset_x,SI.low(asd_model_idx), ...
    'sb','markerfacecolor','b','MarkerSize', 8,'linewidth',1);
plot(2-offset_x,SI.low(nt_model_idx), ...
    'sb','markerfacecolor','w','MarkerSize', 8,'linewidth',1)
plot(1+offset_x,SI.high(asd_model_idx), ...
    'sr','markerfacecolor','r','MarkerSize', 8,'linewidth',1);
plot(1-offset_x,SI.high(nt_model_idx), ...
    'sr','markerfacecolor','w','MarkerSize', 8,'linewidth',1)
axis([0.35 2.65 -0.5 0.25])
set(gca,'XTick',[1 2],'XTickLabel',...
    {'98%','3%'},...
    'YTick',[-1:.25:1],'fontsize',18)
set(gcf,'color','w')
ax = gca; ax.XColor = [0 0 0]; ax.YColor = [0 0 0];
title('Medium vs. big')
ylabel('Model size index','color','k')
xlabel('Contrast','color','k')

box off
set(gcf,'color','w','POS',[680   150   320   400])

%% define output

model.p = p;
model.p.use_e_width = model.p.use_e_width(which_version, :);
model.p.S_gain = model.p.S_gain(which_version, :);
model.p.MW = model.p.MW(which_version, :);
model.stim = stim;
model.resp = resp;
model.thresh = thresh;
model.SI = SI;
if which_version == 1
    model.variable_param.S_gain = p.S_gain(which_version, :);
elseif which_version == 2
    model.variable_param.use_e_width = p.use_e_width(which_version, :);
elseif which_version == 3
    model.variable_param.MW = p.MW(which_version, :);
elseif which_version == 4
    model.variable_param.MW = p.MW(which_version, :);
elseif which_version == 5
    model.variable_param.MW = p.MW(which_version, :);

end
model.which_version = which_version;
model.dimension_key = [num2str(size(resp,1)) ' stimulus sizes (1 = small, ' ...
    '2 = medium, 3 = big) x ' num2str(size(resp,2)) ' groups (1 = ASD model, ' ...
    '2 = NT model) x ' num2str(size(resp,3)) ' stimulus contrasts (1 = low, 2 = high)'];

end