function hf = plot_topology(Pos_bs, Pos_ms, fig_id, AreaLength)
% plot_topology(Pos_bs, Pos_ms)

if nargin < 4
    AreaLength = 20;
end

AnchorSize = 28;
AnchorSizeInner = 5.5;
AgentSize = 25;

Nbs = size(Pos_bs, 1);
Nms = size(Pos_ms, 1);

if nargin < 3
    hf = figure;
else
    hf = figure(fig_id);
end

% hold off
% plot(Pos_bs(1, 1), Pos_bs(1, 2), 'r^', 'MarkerSize', 8);
% hold on
% plot(Pos_ms(1, 1), Pos_ms(1, 2), 'bo', 'MarkerSize', 8);

tx = 0.5 * AreaLength / 20;
ty = - 0.5 * AreaLength / 20;
text(Pos_ms(1, 1) + tx, Pos_ms(1, 2) + ty, '1');

for i = 1:Nbs
    % plot(Pos_bs(i, 1), Pos_bs(i, 2), 'r^', 'MarkerSize', 8);
    plot(Pos_bs(i, 1), Pos_bs(i, 2), 'ro', 'MarkerFaceColor',[1 0 0], 'MarkerSize', 5/9*AnchorSize);
    hold on
    plot(Pos_bs(i, 1), Pos_bs(i, 2), 'wo', 'MarkerFaceColor',[1 1 1], 'MarkerSize', AnchorSizeInner);
    
end

for i = 1:Nms
    % plot(Pos_ms(i, 1), Pos_ms(i, 2), 'bo', 'MarkerSize', 8);
    plot(Pos_ms(i, 1), Pos_ms(i, 2), 'bo', 'MarkerFaceColor',[0 0 1], 'MarkerSize', 5/9*AgentSize);
    if i == 1
        text(Pos_ms(i, 1) + tx, Pos_ms(i, 2) + ty + 3, num2str(i));
    end
end

hold off
set(gca, 'FontSize', 12);
% legend('Anchor', 'Agent');
% xlim([-AreaLength / 2, AreaLength / 2]);
% ylim([-AreaLength / 2, AreaLength / 2]);

ylim([-90 90]);
ylim([-90 90]);
set(gca, 'XTick', -80:40:80);
set(gca, 'YTick', -80:40:80);
set(gca, 'FontSize', 14);

% Configuration
BoxLineWidth = 2;
AspectRatio = 1.26;
XLableDistance_y = - 0.09;     % Units = 'normalized'
YLableDistance_x = - 0.10;      % Units = 'normalized'
CanvasWidth = 560;              % Units = 'pixel'
CanvasHeight = 430;             % Units = 'pixel'
LineWidth = 2;
MarkerSize = 9;
MyFontSize = 14;

% Box
box on
set(gca, 'LineWidth', BoxLineWidth);

axis square
