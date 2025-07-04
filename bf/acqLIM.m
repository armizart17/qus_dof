% baseDir = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\SavedDataQUSPhantom';
baseDir = 'D:\emirandaz\qus\data\SavedDataQUSPhantom';

acqNames = ["","_2","_3"];
deadBand = 0.1e-2;

folders = dir(baseDir);
% folders = folders(3:end);
% folders = folders(10); % only homogeneous
% % folders = folders(5:12);

% 261 phantom
folders = folders(3);


resultsDir = fullfile(baseDir, 'bf');
if ~exist(resultsDir,'dir') mkdir(resultsDir); end

for iFolder = 1:length(folders)
    folderStr = folders(iFolder).name;
    subFolderStr = folderStr + "_F";

    for iAcq = 1:length(acqNames)
        samName = subFolderStr + acqNames(iAcq);
        fileName = fullfile(baseDir, folderStr, subFolderStr, samName);
        presetName = fileName + "_preSet.mat";

        if ~exist(presetName,'file'), continue; end

        if strcmp(extractBetween(samName, 1, 3), "261")
            sos = 1509;
            [rf,x,z,fs] = loadAndBfLinear(fileName, presetName, sos);
        else
            [rf,x,z,fs] = loadAndBfLinear(fileName, presetName);
        end
        
        bMode = db(hilbert(rf));
        bMode = bMode - max(bMode(z>deadBand,:),[],'all');

        figure,
        imagesc(x*1e2, z*1e2, bMode);
        axis image
        colormap gray;
        clim([-60 0]);
        colorbar
        ylabel('Depth [cm]', 'FontSize', 12);
        xlabel('Lateral [cm]', 'FontSize', 12);
        ylim([deadBand*100 5])
        caption = strrep(samName, '_', '-');
        title(caption, 'FontSize', 12);

        saveas(gcf, fullfile(resultsDir,samName+'.png'))
        pause(0.5), close,
        save(fullfile(resultsDir,samName),'rf','x','z','fs')
    end
end