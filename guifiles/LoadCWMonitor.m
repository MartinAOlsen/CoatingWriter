function out = LoadCWMonitor(filename)
    %% Open file and get the lines with usefull information
    F = fopen(filename);
    
    reachedData = 0;
    counter = 1;
    tline = fgetl(F);
    
    while ischar(tline)
        tline = fgetl(F);
        
        if reachedData == 1
            try
                format longg
                Data(counter) = str2double(tline(7:end));
            end
            counter = counter + 1;
        end
        if strfind(tline,'# variables:')
            reachedData = 1;
        end
    end
    
    fclose(F);

    %% Put data into correct structure:
    out.intensity = Data(1);
    out.price = Data(2);
    out.invPrice = 1/out.price;
    out.value = out.intensity / out.price ;
    out.neutrons = Data(3);
    out.lambda_background_intensity = Data(4);
    out.lambda_background_neutrons = Data(5);
    out.worldsize = Data(6);
    out.position_background_intensity = Data(7);
    out.position_background_neutrons = Data(8);
    out.divergence_background_intensity = Data(9);
    out.divergence_background_neutrons = Data(10);
    out.combined_background = out.divergence_background_intensity + out.position_background_intensity;
    
    if isnan(out.price)
        out.price = 10000;
        out.intensity = 0;
    end
    
    %% Fix multiplication problems due to MPI sums 
    
    %% Calculate errors
    try; out.lambda_background_error = out.lambda_background_intensity/sqrt(out.lambda_background_neutrons); end
    try; out.divergence_background_error = out.divergence_background_intensity/sqrt(out.divergence_background_neutrons); end
    try; out.position_background_error = out.position_background_intensity/sqrt(out.position_background_neutrons); end
    try; out.intensity_error = out.intensity/sqrt(out.neutrons); end
    try; out.value_error = (out.intensity/sqrt(out.neutrons))/out.price; end




end