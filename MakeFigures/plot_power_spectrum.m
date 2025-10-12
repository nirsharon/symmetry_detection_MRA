% script name: "plot_power_spectrum"
% making a nice stem figures from the power spectrums

inds = [5, 20, 25];
% to_save = 0;

if to_save
    folder_name = 'success_rate';
    if isfolder(folder_name)
        cd(folder_name);
    else
        mkdir(folder_name)
        cd(folder_name)
    end
end

d = length(ps_samples(:,1));

for i=1:length(inds)

    figure;
    h = stem(ps_samples(:,inds(i)), 'filled', 'LineWidth', 1.5);  % 'filled' circles and thicker lines
    hold on;

    % Improve marker appearance
    h.Marker = 'o';
    h.MarkerSize = 6;
    h.MarkerFaceColor = [0.2, 0.6, 0.8];  % teal-like color
    h.Color = [0.1, 0.2, 0.5];            % darker stem color

    % Axis and title formatting
    % xlabel('X Axis');
    ylabel('Power Spectrum Value');
    % title('Nice Stem Plot Example');
    grid on;
    ax = gca;
    set(ax, 'FontSize', 18);
    set(findall(ax, '-property', 'Interpreter'), 'Interpreter', 'latex')
    xlim([0, d]);
    if to_save
        name_it = ['power_spectrum_', num2str(i)] ;
        saveas(gcf, name_it ,'fig');
        saveas(gcf,name_it,'jpg');
        exportgraphics(gcf, [name_it,'.pdf'], 'ContentType', 'vector');
    end
end

if to_save
    save('data_for_power_spectrum');
    cd '../'
end