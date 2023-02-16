try
        imageObj = ImageDataObject(app.parameters.anatimage_path);
        imageObj = imageObj.readReco;
        imageObj = imageObj.readMethod;
        imageObj = imageObj.readAcqp;
    catch
        warning('Anatomical Image could no be loaded');
    end
    anat_image = squeeze(imageObj.data);
    if (imageObj.Method.PVM_NSPacks == 1)
        anat_image = reshape(anat_image, [size(anat_image) 1]);
    end

    app.parameters.anat_patient_pos = imageObj.Acqp.ACQ_patient_pos;
   

    anat_image = rot90(flip(anat_image, 1));

    app.data.anat_image = anat_image;

    
    app.parameters.mat_anat_image = size(anat_image);
    app.parameters.anat_image_fov_mm = imageObj.Method.PVM_Fov;

    % multislice:
    if imageObj.Method.PVM_SPackArrNSlices > 1
        app.parameters.anat_image_fov_mm = [app.parameters.anat_image_fov_mm ...
            imageObj.Method.PVM_SPackArrNSlices * imageObj.Method.PVM_SliceThick + ...
            (imageObj.Method.PVM_SPackArrNSlices-1) * imageObj.Method.PVM_SPackArrSliceGap];
        % caculate central slice:
        warning('implement method to calculate the central slice:');
        app.parameters.anat_central_slice = round(imageObj.Method.PVM_SPackArrNSlices/2);

    else
        app.parameters.anat_image_fov_mm = [app.parameters.anat_image_fov_mm 0];
    end

    app.parameters.anat_nslices = imageObj.Method.PVM_SPackArrNSlices;

    app.parameters.anat_image_offset_mm = [imageObj.Method.PVM_Phase1Offset(1) ...
                                           imageObj.Method.PVM_Phase2Offset(1) ...
                                           imageObj.Method.PVM_Phase0Offset(1)];

    app.parameters.anat_read_orient = imageObj.Method.PVM_SPackArrReadOrient;
     
%             % define grid that Map sits on:
    [vy, vx, vz, off] = define_grid(app, ...
                        app.parameters.anat_image_fov_mm, ...
                        app.parameters.mat_anat_image, ...                                
                        app.parameters.anat_image_offset_mm, ...
                        app.parameters.anat_patient_pos, ...
                        app.parameters.anat_read_orient);
    app.parameters.anat_image_offset_mm = off;
    app.parameters.anat_gridx = vx;
    app.parameters.anat_gridy = vy;
    app.parameters.anat_gridz = vz;


%             
%             app.parameters.refpow_gridx = vx;
%             app.parameters.refpow_gridy = vy;
    
    
    % calculate refereference power map:
    try
        calc_refpow_anat(app);
    end
