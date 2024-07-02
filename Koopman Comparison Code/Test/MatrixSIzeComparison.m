%Matrix Size Comparison Figure
% 第一张折线图，有两条折线
figure;
plot([3, 5, 7, 9], [3, 5, 7, 9], '-o','LineWidth',2);
hold on;
plot([3, 5, 7, 9], [3, 5, 7, 9], '--x','LineWidth',2);
title('1D Case','FontSize',16,'Interpreter','latex');
xlabel('Trunction Order N');
ylabel('Matrix Size');
legend('Koopman Spectral', 'Carleman', 'Location', 'northwest');
grid on;
saveas(gcf, '1DMatSizeComp.png');

% 第二张折线图
figure;
plot([3, 5, 7, 9], [9, 25, 49, 81], '-o','LineWidth',2);
hold on;
plot([3,5,7,9],[14,62,254,1022],'--x','LineWidth',2)
title('2D Case','FontSize',16,'Interpreter','latex');
xlabel('Trunction Order N');
ylabel('Matrix Size');
legend('Koopman Spectral', 'Carleman', 'Location', 'northwest');
grid on;
saveas(gcf, '2DMatSizeComp.png');

% 第三张折线图
figure;
plot([3, 5, 7, 9], [27, 125, 343, 729], '-o','LineWidth',2);
hold on;
plot([3, 5, 7, 9], [39, 363, 3279, 29523], '--x','LineWidth',2);
title('3D Case','FontSize',16,'Interpreter','latex');
xlabel('Trunction Order N');
ylabel('Matrix Size');
legend('Koopman Spectral', 'Carleman', 'Location', 'northwest');
grid on;
saveas(gcf, '3DMatSizeComp.png');
