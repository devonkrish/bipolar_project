

%% NORMALIZING BY FREQ BAND AND BIPOLAR DISTANCE
%For High Density Grids Only

bandpowers = load('/home/devkrish/bipolar_project/2023/normalizedbandpowers.mat');

delta = bandpowers.delta;
theta = bandpowers.theta;
alpha = bandpowers.alpha;
beta = bandpowers.beta;
gamma = bandpowers.gamma;
highgamma = bandpowers.highgamma;

colors = {[0, 1.0, 0], [0, 0.5, 1.0], [1.0, 0, 0], [0.75, 0, 0.75], [1.0, 0.65, 0], [0, 1.0, 1.0], [0, 0, 0.5], [1.0, 0, 1.0], [0, 0.75, 0.25]};

aggregate = {delta, theta, alpha, beta, gamma, highgamma};
subtitles = {'Delta (2-4Hz)', 'Theta (4-8Hz)', 'Alpha (8-13Hz)', 'Beta (13-25Hz)', 'Gamma (25-50Hz)','High Gamma (50-1200Hz)' };
ticklabels = {'Referential', '4 mm', '8 mm', '12 mm', '16 mm', '20 mm'};

fig = figure();
set(gcf, 'Position', [100, 100, 1100, 750]);

for i = 1:6
    subplot(3, 2, i)
    for j = 1:9
         plot(aggregate{i}(j,:), 'Color', colors{j}, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', colors{j})
         title(subtitles{i}, 'FontSize', 14);
         ylim([-1.5 1.0]);
         yticks(-1.5:0.5:1.0);
         xticks(1:6);
         xticklabels(ticklabels);
         xlabel('Bipolar distance');
         ylabel('Mean log power');
         hold on
         grid on
         
         ax = gca;
         ax.GridColor = [0 0 0];
         ax.GridAlpha = 1;
         
    end   
end

%saveas(fig, '/Users/devonkrish/Desktop/IED/_BipolarReref/Updated figs/Supp Fig 1.png');
%saveas(fig, '/Users/devonkrish/Desktop/IED/_BipolarReref/Updated figs/Supp Fig 1.fig');







