function Time_Series_Main(Axes, PlotParameters)

SystemParameters = PlotParameters.SystemParameters;
FigureParameters = PlotParameters.FigureParameters;
NumericalParameters = PlotParameters.NumericalParameters;

% --- Задание параметров графика --- %

InitializePlot(Axes, SystemParameters, FigureParameters, NumericalParameters);

% --- Нанесение областей D --- %

AreaOfD(Axes, SystemParameters, NumericalParameters);

% --- Нанесение траекторий --- %
Trajectories(Axes, SystemParameters, NumericalParameters);

end

function InitializePlot(ax, SP, FP, NP)

%figure;
cla(ax);
hold(ax, 'on');
grid(ax, 'on');
FontSize = 34;


% D = [pi/2 - SP.Sigma, pi/2 + SP.Sigma];
% 
% if (SP.Sigma > 0) && (SP.Sigma < pi/2)
%     ax.YTick = [0, D(1), pi / 2, D(2), pi, 2 * pi];
%     ax.YTickLabel = {'0', '$\frac{\pi}{2}-\sigma$', '$\frac{\pi}{2}$', '$\frac{\pi}{2}+\sigma$', '$\pi$', '$2\pi$'};
% else if (SP.Sigma == pi/2)
%     ax.YTick = [D(1), pi / 2, D(2), 2 * pi];
%     ax.YTickLabel = {'0', '$\frac{\pi}{2}-\sigma$', '$\frac{\pi}{2}$', '$\frac{\pi}{2}+\sigma$', '$2\pi$'};
%     else if (SP.Sigma > pi/2)
%             ax.YTick = [0, pi / 2, pi, D(2), 2*pi + D(1), 2 * pi];
%             ax.YTickLabel = {'0', '$\frac{\pi}{2}$', '$\pi$', '$\frac{\pi}{2}+\sigma$', '$\frac{\pi}{2}-\sigma$', '$2\pi$'};
%         end
%     end
% end


ax.YTick = [0, pi / 2, pi, 3*pi/2, 2 * pi];
ax.YTickLabel = {'0', '$\frac{\pi}{2}$', '$\pi$', '$\frac{3\pi}{2}$', '$2\pi$'};

ax.TickLabelInterpreter = 'latex';

ax.FontSize = FontSize;

ax.YLim = [0, 2 * pi];
ax.XLim = [0, NP.Tspan];

AxisLabelsFontSize = 30;

ax.YLabel.String = '\phi_1 \phi_2';
ax.YLabel.FontSize = AxisLabelsFontSize;
ax.YAxisLocation = 'left';
%ax.YLabel.Rotation = 0;

ax.XLabel.String = 'Time';
ax.XLabel.FontSize = AxisLabelsFontSize;

end

function AreaOfD(ax, SP, NP)

Sigma = SP.Sigma;
Tspan = NP.Tspan;

D = [pi/2 - Sigma, pi/2 + Sigma];

if(Sigma <= pi/2)
    fill(ax, [0, Tspan, Tspan, 0], [D(1), D(1), D(2), D(2)], [159/255, 1, 160/255], 'EdgeColor', 'none');
else if(Sigma < pi)
        fill(ax, [0, Tspan, Tspan, 0], [0, 0, D(2), D(2)], [159/255, 1, 160/255], 'EdgeColor', 'none');
        fill(ax, [0, Tspan, Tspan, 0], [2*pi + D(1), 2*pi + D(1), 2*pi, 2*pi], [159/255, 1, 160/255], 'EdgeColor', 'none');
    else
        fill(ax, [0, Tspan, Tspan, 0], [0, 0, 2*pi, 2*pi], [97/255, 1, 98/255], 'EdgeColor', 'none');
    end
end

end

function Trajectories(ax, SP, NP)

PlotStep = 1;
LineWidth = 3;

[T, X] = ode23(@(T, X) AdlerSystem(X, SP), [0, NP.Tspan], NP.InitialCondition, NP.Options);

X = mod(X, 2 * pi);
for i = 1:length(X) - 1
   if (abs(X(i, 1) - X(i + 1, 1)) > pi)
       X(i, 1) = NaN;
   end
   if (abs(X(i, 2) - X(i + 1, 2)) > pi)
       X(i, 2) = NaN; 
   end
end
plot(ax, T(1:PlotStep:end), X(1:PlotStep:end, 1), 'b', 'LineWidth', LineWidth);
plot(ax, T(1:PlotStep:end), X(1:PlotStep:end, 2), 'r', 'LineWidth', LineWidth);
end
