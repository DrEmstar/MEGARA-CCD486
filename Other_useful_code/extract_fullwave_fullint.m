function extract_fullwave_fullint(matFileName, outputDir)
    % Load the MAT file
    data = load(matFileName);

    % Extract the 'finaldata' structure
    finaldata = data.finaldata;

    % Extract fullwave and fullint
    fullwave = finaldata.fullwave;          % Vector (1 x N)
    fullint = finaldata.fullint;            % Matrix (125 x N)

    % Check dimensions
    [numSpectra, numPoints] = size(fullint);
    if length(fullwave) ~= numPoints
        error('fullwave length must match the number of columns in fullint');
    end

    % Ensure output directory exists
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % Save each row of fullint along with fullwave and a column of 1s
    for i = 1:numSpectra
        filename = sprintf('spectrum_%03d.dat', i);
        filepath = fullfile(outputDir, filename);

        % Create a matrix with columns: fullwave, fullint(i,:), ones
        spectrumData = [fullwave(:), fullint(i, :)', ones(numPoints, 1)];

        % Write to file
        dlmwrite(filepath, spectrumData, 'delimiter', '\t', 'precision', '%.8f');
    end

    fprintf('Data extraction complete. Files saved in: %s\n', outputDir);
end