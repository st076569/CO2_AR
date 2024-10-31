classdef PlotSettings
    properties
        Range = [1, 1, 65, 40];
        XRange = 1:1:540;
        Path = {'../../QT/build-MixtureSolverCO2Ar-Desktop_Qt_5_14_2_MinGW_64_bit-Debug/results_2024', ...
                'figures'};
        LineStyle = {'-', '--', '-.', ':'};
        Color = {[0.00, 0.00, 1.00], ...
                 [1.00, 0.00, 0.00], ...
                 [0.27, 0.62, 0.19], ...
                 [0.32, 0.59, 1.00], ...
                 [0.20, 0.20, 0.20], ...
                 [0.00, 0.00, 0.00], ...
                 [0.90, 0.40, 0.12]};
        LineWidth = 1;
        FontName = 'default';
        FontSize = 8;
        Position = [2, 2, 8.5, 7];
        Units = "centimeters";
        LabelX = '$\overline{x}$';
        LabelY = {'Temperature, K', ...
                  'Pressure, kPa', ...
                  '$\widetilde{x}_{_{CO_2}}$, \%', ...
                  'Heat flux, kW/m$^{2}$', ...
                  '$\widetilde{H}$, \%', ...
                  '$\zeta/\eta$', ...
                  'Diffusion velocity, m/s', ...
                  '$q$, kW/m$^{2}$', ...
                  '$\pi$, kPa'};
        TicksX = [10^0, 10^1, 10^2, 10^3];
        DeltaInterp = 0.1;
        ItemTokenSize = [15, 18];
    end
    methods
        function SetFont(obj)
            set(0, 'DefaultAxesFontSize', obj.FontSize, ...
                'DefaultAxesFontName', obj.FontName);
            set(0, 'DefaultTextFontSize', obj.FontSize, ...
                'DefaultTextFontName', obj.FontName);
        end
        function SaveFig(obj, name)
            fig = gcf;
            print(fig, obj.Path{2} + "/png/" + name, '-dpng', '-r600')
            print(fig, obj.Path{2} + "/eps/" + name, '-depsc', '-tiff')
            close(fig)
%             set(fig, 'Units', 'centimeters');
%             set(fig, 'PaperPosition', [0, 0, fig.Position(3:4)], ...
%                 'PaperSize', fig.Position(3:4));
%             print(obj.Path{2} + "/" + name, '-dpdf', '-vector');
        end
        function mp = ReadData(obj, name)
            mp = readmatrix(obj.Path{1} + "/" + name, 'Range', obj.Range, ...
                'FileType', 'text');
            mp(:, 3) = mp(:, 3) - mp(1, 3) + 1;
            mp = interp1(1:obj.Range(3), mp, 1:obj.DeltaInterp:obj.Range(3), ...
                'spline');
        end
        function y = GetReadFileName(obj, x, exp)
            y = strcat(x{1}, '_', x{2}, '_M-', x{3}, '_x-', x{4}, '_nT-', x{5}, ...
                '_bV-', x{6}, '_dV-', x{7}, '_p-', x{8}, '_T-', x{9}, ...
                '_T12-', x{10}, '_T3-', x{11}, exp);
        end
        function y = GetSaveFileName(obj, x, innerIndex)
            for i = 1:1:size(innerIndex, 2)
                x{innerIndex(i)} = 'i';
            end
            y = obj.GetReadFileName(x, '');
        end
        function y = GetTitle(obj, x, innerIndex)
            for i = 1:1:size(innerIndex, 2)
                x{innerIndex(i)} = 'Var';
            end
            y = strcat(x{2}, ", CO2:", x{4}, " %, M", x{3}, ", ", x{8}, ...
                " Pa, ", x{9}, ":", x{10}, ":", x{11}, " K");
        end
        function y = GetCellNumbers(obj, x)
            y = [];
            for i = 1:1:size(x, 2)
                y = cat(2, y, size(x{i}, 2));
            end
        end
        function y = GetFromIndex(obj, varInit, index)
            y = {};
            for i = 1:1:size(index, 2)
                y = cat(2, y, varInit{i}{index(i)});
            end
        end
        function y = GetIndexMatrix(obj, varInit)
            s = obj.GetCellNumbers(varInit);
            y = [];
            for i1 = 1:1:s(1)
            for i2 = 1:1:s(2)
            for i3 = 1:1:s(3)
            for i4 = 1:1:s(4)
            for i5 = 1:1:s(5)
            for i6 = 1:1:s(6)
            for i7 = 1:1:s(7)
            for i8 = 1:1:s(8)
            for i9 = 1:1:s(9)
            for i10 = 1:1:s(10)
            for i11 = 1:1:s(11)
                y = cat(1, y, [i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11]);
            end
            end
            end
            end
            end
            end
            end
            end
            end
            end
            end
        end
        
        function PlotFig(obj, position, x, y, color, line, titText, ...
                ylabelText, legendText, colNum, legPos)
            figure("Units", obj.Units, "Position", position)
            p = semilogx(x, y, 'LineWidth', obj.LineWidth);
            for i = 1:1:size(y, 2)
                p(i).LineStyle = line{i};
                p(i).Color = color{i};
            end
            l = legend(legendText, 'Interpreter', 'latex', 'Location', legPos);
            l.NumColumns = colNum;
            l.ItemTokenSize = obj.ItemTokenSize;
            %title(titText, 'FontWeight', 'normal')
            xlabel(obj.LabelX, 'Interpreter', 'latex')
            ylabel(ylabelText, 'Interpreter', 'latex')
            xticks(obj.TicksX)
            axis padded
            grid minor
            set(gca, 'TickLabelInterpreter', 'latex')
        end
        function PlotFigT(obj, fileNames, legendNames, titleName, saveName, position)
            data = cell(size(fileNames));
            x = [];
            y = [];
            col = [];
            lst = {};
            leg = {};
            for i = 1:1:size(fileNames, 2)
                data{i} = obj.ReadData(fileNames{i});
            end
            for i = 1:1:size(fileNames, 2) 
                x = cat(2, x, data{i}(:, 3));
                y = cat(2, y, data{i}(:, 13));
                col = cat(2, col, i);
                lst = cat(2, lst, '-');
                leg = cat(2, leg, '$T$');
            end
            for i = 1:1:size(fileNames, 2)
                x = cat(2, x, data{i}(:, 3));
                y = cat(2, y, data{i}(:, 14));
                col = cat(2, col, i);
                lst = cat(2, lst, '-.');
                leg = cat(2, leg, '$T_{12}$');
            end
            for i = 1:1:size(fileNames, 2) 
                x = cat(2, x, data{i}(:, 3));
                y = cat(2, y, data{i}(:, 15));
                col = cat(2, col, i);
                lst = cat(2, lst, ':');
                leg = cat(2, leg, strcat('$', 'T_{3},\,', legendNames{i}, '$'));
            end
            obj.PlotFig(position, x, y, obj.Color(col), lst, ...
                titleName, obj.LabelY{1}, leg, 3, 'northoutside');
            obj.SaveFig(strcat(saveName,'-T'));
        end
        function PlotFigP(obj, fileNames, legendNames, titleName, saveName, position)
            data = cell(size(fileNames));
            x = [];
            y = [];
            col = [];
            lst = {};
            leg = {};
            for i = 1:1:size(fileNames, 2)
                data{i} = obj.ReadData(fileNames{i});
            end
            for i = 1:1:size(fileNames, 2) 
                x = cat(2, x, data{i}(:, 3));
                y = cat(2, y, data{i}(:, 11) / 1000);
                col = cat(2, col, i);
                lst = cat(2, lst, '-');
                leg = cat(2, leg, '$p$');
            end
            for i = 1:1:size(fileNames, 2) 
                x = cat(2, x, data{i}(:, 3));
                y = cat(2, y, (data{i}(:, 37) - data{i}(:, 11)) / 1000);
                col = cat(2, col, i);
                lst = cat(2, lst, '-.');
                leg = cat(2, leg, strcat('$', '\pi,\,', legendNames{i}, '$'));
            end
            obj.PlotFig(position, x, y, obj.Color(col), lst, ...
                titleName, obj.LabelY{2}, leg, 2, 'northoutside');
            obj.SaveFig(strcat(saveName,'-P'));
        end
        function PlotFigPxx(obj, fileNames, legendNames, titleName, saveName, position)
            data = cell(size(fileNames));
            x = [];
            y = [];
            col = [];
            lst = {};
            leg = {};
            for i = 1:1:size(fileNames, 2)
                data{i} = obj.ReadData(fileNames{i});
            end
            for i = 1:1:size(fileNames, 2) 
                x = cat(2, x, data{i}(obj.XRange, 3));
                y = cat(2, y, (data{i}(obj.XRange, 37) - data{i}(obj.XRange, 11)) / 1000);
                col = cat(2, col, i);
                lst = cat(2, lst, '-');
                leg = cat(2, leg, strcat('$', legendNames{i}, '$'));
            end
            obj.PlotFig(position, x, y, obj.Color(col), lst, ...
                titleName, obj.LabelY{9}, leg, 1, 'northoutside');
            obj.SaveFig(strcat(saveName,'-Pxx'));
        end
        function PlotFigX(obj, fileNames, legendNames, titleName, saveName, position)
            data = cell(size(fileNames));
            x = [];
            y = [];
            col = [];
            lst = {};
            leg = {};
            for i = 1:1:size(fileNames, 2)
                data{i} = obj.ReadData(fileNames{i});
            end
            for i = 1:1:size(fileNames, 2) 
                x = cat(2, x, data{i}(obj.XRange, 3));
                y = cat(2, y, 100 * (data{i}(obj.XRange, 6) - data{i}(1, 6)));
                col = cat(2, col, i);
                lst = cat(2, lst, '-');
                leg = cat(2, leg, strcat('$', legendNames{i}, '$'));
            end
            obj.PlotFig(position, x, y, obj.Color(col), lst, ...
                titleName, obj.LabelY{3}, leg, 1, 'northoutside');
            obj.SaveFig(strcat(saveName,'-XCO2'));
        end
        function PlotFigQ(obj, fileNames, legendNames, titleName, saveName, position)
            data = cell(size(fileNames));
            x = [];
            y = [];
            col = [];
            lst = {};
            leg = {};
            for i = 1:1:size(fileNames, 2)
                data{i} = obj.ReadData(fileNames{i});
            end
            for i = 1:1:size(fileNames, 2) 
                x = cat(2, x, data{i}(obj.XRange, 3));
                y = cat(2, y, sum(data{i}(obj.XRange, 30:34), 2) / 1000);
                col = cat(2, col, i);
                lst = cat(2, lst, '-');
                leg = cat(2, leg, strcat('$', legendNames{i}, '$'));
            end
            obj.PlotFig(position, x, y, obj.Color(col), lst, ...
                titleName, obj.LabelY{8}, leg, 1, 'northoutside');
            obj.SaveFig(strcat(saveName,'-Q'));
        end
        function PlotFigQDV(obj, fileNames, legendNames, titleName, saveName, position)
            data = cell(size(fileNames));
            x = [];
            y = [];
            col = [];
            lst = {};
            leg = {};
            for i = 1:1:size(fileNames, 2)
                data{i} = obj.ReadData(fileNames{i});
            end
            for i = 1:1:size(fileNames, 2)
                x = cat(2, x, data{i}(obj.XRange, 3));
                y = cat(2, y, (data{i}(obj.XRange, 31) + data{i}(obj.XRange, 32)) / 1000);
                col = cat(2, col, i);
                lst = cat(2, lst, '-');
                leg = cat(2, leg, '$q_{vibr}$');
            end
            for i = 1:1:size(fileNames, 2) 
                x = cat(2, x, data{i}(obj.XRange, 3));
                y = cat(2, y, (data{i}(obj.XRange, 34) + data{i}(obj.XRange, 33)) / 1000);
                col = cat(2, col, i);
                lst = cat(2, lst, '-.');
                leg = cat(2, leg, strcat('$', 'q_{diff},\,', legendNames{i}, '$'));
            end
            obj.PlotFig(position, x, y, obj.Color(col), lst, ...
                titleName, obj.LabelY{4}, leg, 2, 'northoutside');
            obj.SaveFig(strcat(saveName,'-QDV'));
        end
        function PlotFigQD(obj, fileNames, legendNames, titleName, saveName, position)
            data = cell(size(fileNames));
            x = [];
            y = [];
            col = [];
            lst = {};
            leg = {};
            for i = 1:1:size(fileNames, 2)
                data{i} = obj.ReadData(fileNames{i});
            end
            for i = 1:1:size(fileNames, 2)
                x = cat(2, x, data{i}(:, 3));
                y = cat(2, y, data{i}(:, 33) / 1000);
                col = cat(2, col, i);
                lst = cat(2, lst, '-');
                leg = cat(2, leg, '$q_{diff}$');
            end
            for i = 1:1:size(fileNames, 2) 
                x = cat(2, x, data{i}(:, 3));
                y = cat(2, y, data{i}(:, 34) / 1000);
                col = cat(2, col, i);
                lst = cat(2, lst, '-.');
                leg = cat(2, leg, strcat('$', 'q_{T_{diff}},\,', legendNames{i}, '$'));
            end
            obj.PlotFig(position, x, y, obj.Color(col), lst, ...
                titleName, obj.LabelY{4}, leg, 2, 'northoutside');
            obj.SaveFig(strcat(saveName,'-QD'));
        end
        function PlotFigQV(obj, fileNames, legendNames, titleName, saveName, position)
            data = cell(size(fileNames));
            x = [];
            y = [];
            col = [];
            lst = {};
            leg = {};
            for i = 1:1:size(fileNames, 2)
                data{i} = obj.ReadData(fileNames{i});
            end
            for i = 1:1:size(fileNames, 2)
                x = cat(2, x, data{i}(:, 3));
                y = cat(2, y, data{i}(:, 31) / 1000);
                col = cat(2, col, i);
                lst = cat(2, lst, '-');
                leg = cat(2, leg, '$q_{12}$');
            end
            for i = 1:1:size(fileNames, 2) 
                x = cat(2, x, data{i}(:, 3));
                y = cat(2, y, data{i}(:, 32) / 1000);
                col = cat(2, col, i);
                lst = cat(2, lst, '-.');
                leg = cat(2, leg, strcat('$', 'q_{3},\,', legendNames{i}, '$'));
            end
            obj.PlotFig(position, x, y, obj.Color(col), lst, ...
                titleName, obj.LabelY{4}, leg, 2, 'northoutside');
            obj.SaveFig(strcat(saveName,'-QV'));
        end
        function PlotFigQTD(obj, fileNames, legendNames, titleName, saveName, position)
            data = cell(size(fileNames));
            x = [];
            y = [];
            col = [];
            lst = {};
            leg = {};
            for i = 1:1:size(fileNames, 2)
                data{i} = obj.ReadData(fileNames{i});
            end
            for i = 1:1:size(fileNames, 2)
                x = cat(2, x, data{i}(:, 3));
                y = cat(2, y, data{i}(:, 30) / 1000);
                col = cat(2, col, i);
                lst = cat(2, lst, '-');
                leg = cat(2, leg, '$q_{tr-rot}$');
            end
            for i = 1:1:size(fileNames, 2) 
                x = cat(2, x, data{i}(:, 3));
                y = cat(2, y, (data{i}(:, 33) + data{i}(:, 34)) / 1000);
                col = cat(2, col, i);
                lst = cat(2, lst, '-.');
                leg = cat(2, leg, strcat('$', 'q_{diff},\,', legendNames{i}, '$'));
            end
            obj.PlotFig(position, x, y, obj.Color(col), lst, ...
                titleName, obj.LabelY{4}, leg, 2, 'northoutside');
            obj.SaveFig(strcat(saveName,'-QTD'));
        end
        function PlotFigH(obj, fileNames, legendNames, titleName, saveName, position)
            data = cell(size(fileNames));
            x = [];
            y = [];
            col = [];
            lst = {};
            leg = {};
            for i = 1:1:size(fileNames, 2)
                data{i} = obj.ReadData(fileNames{i});
            end
            for i = 1:1:size(fileNames, 2) 
                x = cat(2, x, data{i}(obj.XRange, 3));
                h = data{i}(obj.XRange, 38) + 0.5 * data{i}(obj.XRange, 12) .^ 2 + ...
                    data{i}(obj.XRange, 11) ./ (data{i}(obj.XRange, 9) + data{i}(obj.XRange, 10));
                y = cat(2, y, (h - h(1)) / h(1));
                col = cat(2, col, i);
                lst = cat(2, lst, '-');
                leg = cat(2, leg, strcat('$', legendNames{i}, '$'));
            end
            obj.PlotFig(position, x, y, obj.Color(col), lst, ...
                titleName, obj.LabelY{5}, leg, 1, 'northoutside');
            obj.SaveFig(strcat(saveName,'-H'));
        end
        function PlotFigZE(obj, fileNames, legendNames, titleName, saveName, position)
            data = cell(size(fileNames));
            x = [];
            y = [];
            col = [];
            lst = {};
            leg = {};
            for i = 1:1:size(fileNames, 2)
                data{i} = obj.ReadData(fileNames{i});
            end
            for i = 1:1:size(fileNames, 2) 
                x = cat(2, x, data{i}(:, 3));
                y = cat(2, y, data{i}(:, 36) ./ data{i}(:, 35));
                col = cat(2, col, i);
                lst = cat(2, lst, '-');
                leg = cat(2, leg, strcat('$', legendNames{i}, '$'));
            end
            obj.PlotFig(position, x, y, obj.Color(col), lst, ...
                titleName, obj.LabelY{6}, leg, 1, 'northoutside');
            obj.SaveFig(strcat(saveName,'-ZE'));
        end
        function PlotFigVd(obj, fileNames, legendNames, titleName, saveName, position)
            data = cell(size(fileNames));
            x = [];
            y = [];
            col = [];
            lst = {};
            leg = {};
            for i = 1:1:size(fileNames, 2)
                data{i} = obj.ReadData(fileNames{i});
            end
            for i = 1:1:size(fileNames, 2)
                x = cat(2, x, data{i}(:, 3));
                y = cat(2, y, data{i}(:, 19));
                col = cat(2, col, i);
                lst = cat(2, lst, '-');
                leg = cat(2, leg, '$V_{d,\,CO_2}$');
            end
            for i = 1:1:size(fileNames, 2) 
                x = cat(2, x, data{i}(:, 3));
                y = cat(2, y, data{i}(:, 20));
                col = cat(2, col, i);
                lst = cat(2, lst, '-.');
                leg = cat(2, leg, strcat('$', 'V_{d,\,atomic},\,', legendNames{i}, '$'));
            end
            obj.PlotFig(position, x, y, obj.Color(col), lst, ...
                titleName, obj.LabelY{7}, leg, 2, 'northoutside');
            obj.SaveFig(strcat(saveName,'-dV'));
        end
        
        function PlotFigs(obj, varInit, innerIndex, legendNames, plotFig)
            outerIndex = 1:1:size(varInit, 2);
            outerIndex(innerIndex) = [];
            m = obj.GetIndexMatrix(varInit);
            while size(m, 1) > 0
                gi = logical(prod(m(:, outerIndex) == m(1, outerIndex), 2));
                group = m(gi, :);
                y = {};
                for i = 1:1:size(group, 1)
                    sec = obj.GetFromIndex(varInit, group(i, :));
                    y = cat(2, y, obj.GetReadFileName(sec, '.csd'));
                end
                saveName = obj.GetSaveFileName(sec, innerIndex);
                titleName = obj.GetTitle(sec, innerIndex);
                m = m(gi == 0, :);
                plotFig(y, legendNames, titleName, saveName, obj.Position);
            end
        end
        function PlotFigAll(obj, varInit, innerIndex, legendNames)
            obj.PlotFigs(varInit, innerIndex, legendNames, @obj.PlotFigT);
            obj.PlotFigs(varInit, innerIndex, legendNames, @obj.PlotFigP);
            obj.PlotFigs(varInit, innerIndex, legendNames, @obj.PlotFigPxx);
            obj.PlotFigs(varInit, innerIndex, legendNames, @obj.PlotFigX);
            obj.PlotFigs(varInit, innerIndex, legendNames, @obj.PlotFigQ);
            obj.PlotFigs(varInit, innerIndex, legendNames, @obj.PlotFigQDV);
            obj.PlotFigs(varInit, innerIndex, legendNames, @obj.PlotFigQV);
            obj.PlotFigs(varInit, innerIndex, legendNames, @obj.PlotFigH);
        end
    end
end