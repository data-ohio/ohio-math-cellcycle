%% =========================
%  Gen point scatter plots
%  =========================

%% Set-up

% Adjust the save folder for the figures
% Adjust marker size and figure scale as needed


%% Settings

% Marker size
markerSize = 15;

% Font / axes settings
fontName      = 'Arial';
fontSizeAxis  = 11;
fontSizeLabel = 12;
fontSizeTitle = 12;
legendSize    = 10;
axisLineWidth = 1.0;

% Figure size in inches
figWidth  = 5.5;
figHeight = 4.2;

% Save folder (UPDATE BEFORE RUNNING)
saveFolder = 'C:\';
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

% Color palette
c1 = [0.0000, 0.4470, 0.7410];   % blue
c2 = [0.8500, 0.3250, 0.0980];   % orange
c3 = [0.4940, 0.1840, 0.5560];   % purple
c4 = [0.4660, 0.6740, 0.1880];   % green
c5 = [0.3010, 0.7450, 0.9330];   % cyan
c6 = [0.6350, 0.0780, 0.1840];   % dark red
c7 = [0.25,   0.25,   0.25  ];   % dark gray

%% Default graphics settings
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesFontName', fontName);
set(groot, 'defaultAxesFontSize', fontSizeAxis);
set(groot, 'defaultAxesLineWidth', axisLineWidth);
set(groot, 'defaultAxesBox', 'on');
set(groot, 'defaultAxesLayer', 'top');
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultLineLineWidth', 1.0);

% Uncomment these if your MATLAB version supports them well
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% Helper function for figure creation
makeFigure = @() figure( ...
    'Units', 'inches', ...
    'Position', [1, 1, figWidth, figHeight], ...
    'PaperUnits', 'inches', ...
    'PaperSize', [figWidth, figHeight], ...
    'PaperPosition', [0, 0, figWidth, figHeight], ...
    'Color', 'w');

%% -------------------------
% Figure 1: Min-Max Density
% -------------------------
fig1 = makeFigure();
ax1 = axes(fig1);
hold(ax1, 'on');

scatter(gen_pt, MaxDensity, markerSize, ...
    'MarkerEdgeColor', c1, ...
    'MarkerFaceColor', c1);

scatter(gen_pt, MinDensity, markerSize, ...
    'MarkerEdgeColor', c2, ...
    'MarkerFaceColor', c2);

title('Min--Max Density', 'FontSize', fontSizeTitle, 'FontWeight', 'normal');
xlabel('$g$', 'FontSize', fontSizeLabel);
ylabel('Density', 'FontSize', fontSizeLabel);
legend({'Max density', 'Min density'}, ...
    'Location', 'best', ...
    'FontSize', legendSize, ...
    'Box', 'off');
grid on;
xlim([min(gen_pt), max(gen_pt)]);
hold(ax1, 'off');

%% -------------------------
% Figure 2: Min-Max Budding Index
% -------------------------
fig2 = makeFigure();
ax2 = axes(fig2);
hold(ax2, 'on');

scatter(gen_pt, MaxBud, markerSize, ...
    'MarkerEdgeColor', c1, ...
    'MarkerFaceColor', c1);

scatter(gen_pt, MinBud, markerSize, ...
    'MarkerEdgeColor', c2, ...
    'MarkerFaceColor', c2);

title('Min--Max Budding Index', 'FontSize', fontSizeTitle, 'FontWeight', 'normal');
xlabel('$g$', 'FontSize', fontSizeLabel);
ylabel('Budding index', 'FontSize', fontSizeLabel);
ylim([0, 1]);
xlim([min(gen_pt), max(gen_pt)]);
legend({'Max budding index', 'Min budding index'}, ...
    'Location', 'best', ...
    'FontSize', legendSize, ...
    'Box', 'off');
grid on;
hold(ax2, 'off');

%% -------------------------
% Figure 3: Cluster SD
% -------------------------
fig3 = makeFigure();
ax3 = axes(fig3);
hold(ax3, 'on');

scatter(gen_pt, ClusterSD, markerSize, ...
    'MarkerEdgeColor', c3, ...
    'MarkerFaceColor', c3);

title(sprintf('Weighted Average SD of Cell Location in %d Clusters', kclus), ...
    'FontSize', fontSizeTitle, 'FontWeight', 'normal');
xlabel('$g$', 'FontSize', fontSizeLabel);
ylabel('Average SD', 'FontSize', fontSizeLabel);
ylim([0, 0.12]);
xlim([min(gen_pt), max(gen_pt)]);
grid on;
hold(ax3, 'off');

%% -------------------------
% Figure 4: Average fraction in each generation
% -------------------------
fig4 = figure( ...
    'Units', 'inches', ...
    'Position', [1, 1, 6.5, 4.6], ...
    'PaperUnits', 'inches', ...
    'PaperSize', [6.5, 4.6], ...
    'PaperPosition', [0, 0, 6.5, 4.6], ...
    'Color', 'w');

ax4 = axes(fig4);
hold(ax4, 'on');

scatter(gen_pt, AvgD,  markerSize, 'MarkerEdgeColor', c7, 'MarkerFaceColor', c7);
scatter(gen_pt, AvgM1, markerSize, 'MarkerEdgeColor', c1, 'MarkerFaceColor', c1);
scatter(gen_pt, AvgM2, markerSize, 'MarkerEdgeColor', c2, 'MarkerFaceColor', c2);
scatter(gen_pt, AvgM3, markerSize, 'MarkerEdgeColor', c3, 'MarkerFaceColor', c3);
scatter(gen_pt, AvgM4, markerSize, 'MarkerEdgeColor', c4, 'MarkerFaceColor', c4);
scatter(gen_pt, AvgM5, markerSize, 'MarkerEdgeColor', c5, 'MarkerFaceColor', c5);
scatter(gen_pt, AvgM6, markerSize, 'MarkerEdgeColor', c6, 'MarkerFaceColor', c6);

title('Average Fraction in Each Generation Over Final 2 Cycles', ...
    'FontSize', fontSizeTitle, 'FontWeight', 'normal');
xlabel('$g$', 'FontSize', fontSizeLabel);
ylabel('Average fraction of cells', 'FontSize', fontSizeLabel);
ylim([0, 1]);
xlim([min(gen_pt), max(gen_pt)]);
legend({'D', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6+'}, ...
    'Location', 'eastoutside', ...
    'FontSize', legendSize, ...
    'Box', 'off');
grid on;
hold(ax4, 'off');

%% -------------------------
% Figure 5: Density ratio
% -------------------------
fig5 = makeFigure();
ax5 = axes(fig5);
hold(ax5, 'on');

scatter(gen_pt, MaxDensity ./ MinDensity, markerSize, ...
    'MarkerEdgeColor', c6, ...
    'MarkerFaceColor', c6);

title('Max/Min Density Ratio', 'FontSize', fontSizeTitle, 'FontWeight', 'normal');
xlabel('$g$', 'FontSize', fontSizeLabel);
ylabel('Density ratio', 'FontSize', fontSizeLabel);
xlim([min(gen_pt), max(gen_pt)]);
grid on;
hold(ax5, 'off');

%% -------------------------
% Export figures
% -------------------------
exportgraphics(fig1, fullfile(saveFolder, 'MinMaxDensity.pdf'), 'ContentType', 'vector');
exportgraphics(fig2, fullfile(saveFolder, 'MinMaxBuddingIndex.pdf'), 'ContentType', 'vector');
exportgraphics(fig3, fullfile(saveFolder, 'ClusterSD.pdf'), 'ContentType', 'vector');
exportgraphics(fig4, fullfile(saveFolder, 'AvgGenerationFractions.pdf'), 'ContentType', 'vector');
exportgraphics(fig5, fullfile(saveFolder, 'DensityRatio.pdf'), 'ContentType', 'vector');

exportgraphics(fig1, fullfile(saveFolder, 'MinMaxDensity.png'), 'Resolution', 600);
exportgraphics(fig2, fullfile(saveFolder, 'MinMaxBuddingIndex.png'), 'Resolution', 600);
exportgraphics(fig3, fullfile(saveFolder, 'ClusterSD.png'), 'Resolution', 600);
exportgraphics(fig4, fullfile(saveFolder, 'AvgGenerationFractions.png'), 'Resolution', 600);
exportgraphics(fig5, fullfile(saveFolder, 'DensityRatio.png'), 'Resolution', 600);