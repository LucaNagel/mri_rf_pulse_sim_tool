classdef simulate_rf_pulse_b1_tool_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        calc_rf_pulse_toolUIFigure     matlab.ui.Figure
        GridLayout                     matlab.ui.container.GridLayout
        LeftPanel                      matlab.ui.container.Panel
        reset_B1map_Button             matlab.ui.control.Button
        startparpoolButton             matlab.ui.control.Button
        progress_sim_repLabel          matlab.ui.control.Label
        SimulateRepsButton             matlab.ui.control.Button
        SimulaterepeatedExcitationFISPbSSFPLabel  matlab.ui.control.Label
        SimulateRepetitionsTabGroup    matlab.ui.container.TabGroup
        BasicParameterTab              matlab.ui.container.Tab
        first_repetition_time_msEditField  matlab.ui.control.NumericEditField
        first_repetition_time_msLabel  matlab.ui.control.Label
        TipbackCheckBox                matlab.ui.control.CheckBox
        SpoilerGradCheckBox            matlab.ui.control.CheckBox
        PlotresultsButton              matlab.ui.control.Button
        freq_1stpulse_HzEditField      matlab.ui.control.NumericEditField
        freq_1stpulse_EditFieldLabel   matlab.ui.control.Label
        amp_1stpulse_EditField         matlab.ui.control.NumericEditField
        FirstpulseampEditFieldLabel    matlab.ui.control.Label
        progress_simRepLabel_10        matlab.ui.control.Label
        progress_simRepLabel_9         matlab.ui.control.Label
        progress_simRepLabel_8         matlab.ui.control.Label
        progress_simRepLabel_7         matlab.ui.control.Label
        progress_simRepLabel_6         matlab.ui.control.Label
        progress_simRepLabel_5         matlab.ui.control.Label
        progress_simRepLabel_4         matlab.ui.control.Label
        progress_simRepLabel_3         matlab.ui.control.Label
        progress_simRepLabel_2         matlab.ui.control.Label
        progress_simRepLabel           matlab.ui.control.Label
        incphaseLabel                  matlab.ui.control.Label
        NRepsLabel                     matlab.ui.control.Label
        NpointsLabel_3                 matlab.ui.control.Label
        TRmsLabel                      matlab.ui.control.Label
        incphaseEditField              matlab.ui.control.NumericEditField
        repetition_time_npointsEditField  matlab.ui.control.NumericEditField
        nrepetitions_EditField         matlab.ui.control.NumericEditField
        repetition_time_msEditField    matlab.ui.control.NumericEditField
        AdvParameterTab                matlab.ui.container.Tab
        phasetbEditField               matlab.ui.control.NumericEditField
        phasetbEditFieldLabel          matlab.ui.control.Label
        alpha2prepCheckBox             matlab.ui.control.CheckBox
        highresoutCheckBox             matlab.ui.control.CheckBox
        homoB1CheckBox                 matlab.ui.control.CheckBox
        altphaseEditField              matlab.ui.control.NumericEditField
        altphaseEditFieldLabel         matlab.ui.control.Label
        freq_2ndpulse_HzEditField      matlab.ui.control.NumericEditField
        fHzEditFieldLabel              matlab.ui.control.Label
        AlterNTimesEditField           matlab.ui.control.NumericEditField
        AlterNTimesLabel               matlab.ui.control.Label
        metabsCheckBox                 matlab.ui.control.CheckBox
        alpha2ndPulseEditField         matlab.ui.control.NumericEditField
        alphaEditFieldLabel            matlab.ui.control.Label
        TRGmsLabel                     matlab.ui.control.Label
        duration_refocus_grad          matlab.ui.control.NumericEditField
        SimrangeTab                    matlab.ui.container.Tab
        plotresCheckBox                matlab.ui.control.CheckBox
        TrfEditFieldLabel              matlab.ui.control.Label
        TrfRange_1EditField            matlab.ui.control.NumericEditField
        TrfRange_2EditField            matlab.ui.control.NumericEditField
        TrfRange_NEditField            matlab.ui.control.NumericEditField
        simrangeCheckBox               matlab.ui.control.CheckBox
        SimRangeButton                 matlab.ui.control.Button
        T2Label_3                      matlab.ui.control.Label
        FARange_1EditField             matlab.ui.control.NumericEditField
        FARange_2EditField             matlab.ui.control.NumericEditField
        FARange_NEditField             matlab.ui.control.NumericEditField
        TRRange_NEditField             matlab.ui.control.NumericEditField
        TRRange_2EditField             matlab.ui.control.NumericEditField
        TRRange_1EditField             matlab.ui.control.NumericEditField
        T2Label_2                      matlab.ui.control.Label
        T2Range_1EditField             matlab.ui.control.NumericEditField
        T2Label                        matlab.ui.control.Label
        T2Range_NEditField             matlab.ui.control.NumericEditField
        T2Range_2EditField             matlab.ui.control.NumericEditField
        T1Range_NEditField             matlab.ui.control.NumericEditField
        T1Range_2EditField             matlab.ui.control.NumericEditField
        T1Range_1EditField             matlab.ui.control.NumericEditField
        T1EditFieldLabel               matlab.ui.control.Label
        nslices_3                      matlab.ui.control.Label
        nslices_2                      matlab.ui.control.Label
        nslices                        matlab.ui.control.Label
        TabGroup                       matlab.ui.container.TabGroup
        findGradTab                    matlab.ui.container.Tab
        find_GradButton                matlab.ui.control.Button
        findGradientstrengthforB1compensationPanel  matlab.ui.container.Panel
        progress_grad_findLabel        matlab.ui.control.Label
        pointsyLabel                   matlab.ui.control.Label
        NpointsLabel_2                 matlab.ui.control.Label
        G2kHzcmLabel                   matlab.ui.control.Label
        G1kHzcmLabel                   matlab.ui.control.Label
        find_Gradients_ind2EditField   matlab.ui.control.NumericEditField
        find_Gradients_ind1EditField   matlab.ui.control.NumericEditField
        BestGradEditField              matlab.ui.control.NumericEditField
        BestGradEditFieldLabel         matlab.ui.control.Label
        find_Gradients_npointsEditField  matlab.ui.control.NumericEditField
        find_Gradients_highEditField   matlab.ui.control.NumericEditField
        find_Gradients_lowEditField    matlab.ui.control.NumericEditField
        RFpowerTab                     matlab.ui.container.Tab
        WButton                        matlab.ui.control.Button
        overlayrefpowCheckBox          matlab.ui.control.CheckBox
        xmmEditField                   matlab.ui.control.NumericEditField
        xmmEditFieldLabel              matlab.ui.control.Label
        RefpowLabel                    matlab.ui.control.Label
        calcrefpowButton               matlab.ui.control.Button
        ymmEditField                   matlab.ui.control.NumericEditField
        ymmEditFieldLabel              matlab.ui.control.Label
        findRFPowerCheckBox            matlab.ui.control.CheckBox
        RFPowerWEditFieldLabel         matlab.ui.control.Label
        RFPowerLabel                   matlab.ui.control.Label
        RFPowerWEditField              matlab.ui.control.NumericEditField
        usemapCheckBox                 matlab.ui.control.CheckBox
        saveTab                        matlab.ui.container.Tab
        save_fileEditField             matlab.ui.control.EditField
        fileEditFieldLabel             matlab.ui.control.Label
        saveButton                     matlab.ui.control.Button
        repEditField_high              matlab.ui.control.NumericEditField
        repEditField_low               matlab.ui.control.NumericEditField
        repEditFieldLabel              matlab.ui.control.Label
        saveoptions_ListBox            matlab.ui.control.ListBox
        b1_mag_khz_EditField           matlab.ui.control.NumericEditField
        alphaLabel_4                   matlab.ui.control.Label
        SimulationParametersPanel      matlab.ui.container.Panel
        valueLabel                     matlab.ui.control.Label
        freq_hz_resLabel               matlab.ui.control.Label
        yrange_cm_resLabel             matlab.ui.control.Label
        slice_thick                    matlab.ui.control.Label
        pulseduration_ms_resLabel      matlab.ui.control.Label
        resolutionLabel                matlab.ui.control.Label
        freq_hz_npointsEditField       matlab.ui.control.EditField
        freq_hz_highEditField          matlab.ui.control.NumericEditField
        freq_hz_lowEditField           matlab.ui.control.NumericEditField
        freqHzEditFieldLabel           matlab.ui.control.Label
        yrange_cm_npointsEditField     matlab.ui.control.EditField
        EditField_2                    matlab.ui.control.EditField
        pulseduration_ms_npointsEditField  matlab.ui.control.EditField
        pointsLabel                    matlab.ui.control.Label
        yrange_cmEditField             matlab.ui.control.NumericEditField
        ycmLabel                       matlab.ui.control.Label
        gradientstrength_khzcmEditField  matlab.ui.control.NumericEditField
        GradientStrengthkHzcmEditFieldLabel  matlab.ui.control.Label
        SimulateButton                 matlab.ui.control.Button
        pulseduration_msEditField      matlab.ui.control.NumericEditField
        PulseDurationmsEditFieldLabel  matlab.ui.control.Label
        SamplePanel                    matlab.ui.container.Panel
        HPEditField                    matlab.ui.control.NumericEditField
        HPEditFieldLabel               matlab.ui.control.Label
        T2_sEditField                  matlab.ui.control.EditField
        T2sEditFieldLabel              matlab.ui.control.Label
        T1_sEditField                  matlab.ui.control.EditField
        T1sEditFieldLabel              matlab.ui.control.Label
        NucleusButtonGroup             matlab.ui.container.ButtonGroup
        CButton                        matlab.ui.control.RadioButton
        HButton                        matlab.ui.control.RadioButton
        SelectAnatImage_Button         matlab.ui.control.Button
        calcButton                     matlab.ui.control.Button
        LoadAnatImage_Button           matlab.ui.control.Button
        AnatImageEditField             matlab.ui.control.EditField
        AnatImageEditFieldLabel        matlab.ui.control.Label
        alphaLabel_3                   matlab.ui.control.Label
        bwfac_HzsEditField             matlab.ui.control.NumericEditField
        sint_EditField                 matlab.ui.control.NumericEditField
        alphaLabel_2                   matlab.ui.control.Label
        flipangle_degEditField         matlab.ui.control.NumericEditField
        alphaLabel                     matlab.ui.control.Label
        LoadB1Profile_Button           matlab.ui.control.Button
        SelectB1Profile_Button         matlab.ui.control.Button
        LoadRFPulse_Button             matlab.ui.control.Button
        SelectRFPulse_Button           matlab.ui.control.Button
        B1ProfileEditField             matlab.ui.control.EditField
        B1ProfileEditFieldLabel        matlab.ui.control.Label
        RFPulseEditField               matlab.ui.control.EditField
        RFPulseEditFieldLabel          matlab.ui.control.Label
        RightPanel                     matlab.ui.container.Panel
        repSlider                      matlab.ui.control.Slider
        repSliderLabel                 matlab.ui.control.Label
        responseprofileCheckBox        matlab.ui.control.CheckBox
        bssfpsignalCheckBox            matlab.ui.control.CheckBox
        refpow_map_caxis               matlab.ui.control.Slider
        homoB1ListBox                  matlab.ui.control.ListBox
        homoB1ListBoxLabel             matlab.ui.control.Label
        inhomoB1ListBox                matlab.ui.control.ListBox
        inhomoB1ListBoxLabel           matlab.ui.control.Label
        UseFrequenciesCheckBox         matlab.ui.control.CheckBox
        FrequenciesPanel               matlab.ui.container.Panel
        HzLabel_2                      matlab.ui.control.Label
        HzLabel                        matlab.ui.control.Label
        Freq2EditField                 matlab.ui.control.NumericEditField
        Freq2EditFieldLabel            matlab.ui.control.Label
        Freq1EditField                 matlab.ui.control.NumericEditField
        Freq1EditFieldLabel            matlab.ui.control.Label
        MetabolitesPanel               matlab.ui.container.Panel
        PulseFreqEditField             matlab.ui.control.NumericEditField
        HzLabel_3                      matlab.ui.control.Label
        PulseFreqEditFieldLabel        matlab.ui.control.Label
        ppmLabel_3                     matlab.ui.control.Label
        BasicFreqEditField             matlab.ui.control.NumericEditField
        BasicFreqEditFieldLabel        matlab.ui.control.Label
        ppmLabel_2                     matlab.ui.control.Label
        ppmLabel                       matlab.ui.control.Label
        Freq2_ppmEditField             matlab.ui.control.NumericEditField
        Freq2EditField_2Label          matlab.ui.control.Label
        Freq1_ppmEditField             matlab.ui.control.NumericEditField
        Freq1EditField_2Label          matlab.ui.control.Label
        SliceSpinner                   matlab.ui.control.Spinner
        SliceSpinnerLabel              matlab.ui.control.Label
        ShowExcitationProfileCheckBox  matlab.ui.control.CheckBox
        Show2FreqsCheckBox             matlab.ui.control.CheckBox
        freqSlider                     matlab.ui.control.Slider
        freqSliderLabel                matlab.ui.control.Label
        UIAxes_anat_image              matlab.ui.control.UIAxes
        UIAxes_mag_vs_freq             matlab.ui.control.UIAxes
        UIAxes_rfpulse                 matlab.ui.control.UIAxes
        UIAxes_b1profile               matlab.ui.control.UIAxes
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    
    properties (Access = private)
        parameters = struct();% Description
        data = struct();% Description
    end
    
    methods (Access = private)
        function [] = the_updater_func(app)
            app.parameters.khz_in_g = app.parameters.gamma / 1000;
            
            % check wether 2 frequency option can be activated:
            if (app.parameters.freq_hz_npoints == 1)
                app.Show2FreqsCheckBox.Enable = 0;
                app.Show2FreqsCheckBox.Value = 0;
            else
                app.Show2FreqsCheckBox.Enable = 1;
            end
            

        end
        
        function [grad_kHz_cm_range_sim, gradient_G_cm_range_sim] = calc_gradients(app, dt_sim, Trf)
            % grad_Gcm = grad_kHzcm / khz_in_g;
            khz_in_g = app.parameters.gamma / 1000;
            gradient_kHz_cm = app.parameters.gradient_kHz_cm;
            gradient_G_cm = app.parameters.gradient_kHz_cm / khz_in_g;

            app.parameters.gradient_G_cm = gradient_G_cm;  
            % range to simulate over
            % range_cm_sim = app.parameters.khz_in_g; % 4 cm
            
            % number of points:
            grad_on_points = round((Trf / dt_sim));
            grad_range_sim = ones(grad_on_points, 1);
            % plus/minus
            grad_range_sim = [grad_range_sim.' -grad_range_sim.'/2]; % rephasing gradients on
            
            % multiply by gradient amplitude (G/cm)
            gradient_G_cm_range_sim = grad_range_sim * gradient_G_cm;

            % multiply by gradient amplitude (kHz/cm)
            grad_kHz_cm_range_sim = grad_range_sim * gradient_kHz_cm;
            
            
            % save in app.data:
            app.data.gradient_G_cm_range_sim = gradient_G_cm_range_sim;
            app.data.grad_kHz_cm_range_sim = grad_kHz_cm_range_sim ;
        end
        
        function [] = calc_freq_hz_res(app)
            
            app.parameters.freq_hz_res = (app.parameters.freq_hz_high + 1 - ...
                                       app.parameters.freq_hz_low) / ...
                                       app.parameters.freq_hz_npoints;
            app.freq_hz_resLabel.Text = num2str(app.parameters.freq_hz_res);

            app.parameters.freq_hz_range = linspace(app.parameters.freq_hz_low, ...
                                       app.parameters.freq_hz_high, ...
                                       app.parameters.freq_hz_npoints);
            % update slider of result 
            app.freqSlider.Limits = [app.freq_hz_lowEditField.Value 
                                     app.freq_hz_highEditField.Value];

        end
        
        function [] = calc_yrange_cm_res(app)
            app.parameters.yrange_cm_res = app.parameters.yrange_cm / ...
                                           app.parameters.yrange_cm_npoints;
            app.yrange_cm_resLabel.Text = num2str(app.parameters.yrange_cm_res);
        end
        
        function [] = calc_pulse_duration_ms_res(app)
            app.parameters.pulse_duration_ms_res = app.parameters.pulse_duration_ms / ...
                                                   app.parameters.pulse_duration_ms_npoints;
            app.pulseduration_ms_resLabel.Text = num2str(app.parameters.pulse_duration_ms_res);
            
        end
        
        function [] = calc_sint(app)
            try
                pulse = app.data.rfpulse;
                sint = sum(real(pulse))/numel(pulse);
                app.parameters.sint =sint;
                app.sint_EditField.Value = sint;
                bwfac = 1000/sint;
                app.parameters.bwfac_hz_ms = bwfac;
            catch
                warning('Integration factor sint couldnt be calculated!');
            end
            app.bwfac_HzsEditField.Value = app.parameters.bwfac_hz_ms;
            app.sint_EditField.Value = app.parameters.sint;
        end
        
        function [] = plot_rfpulse(app)
            % load pulse:
            pulse = app.data.rfpulse;
            % parameters
            % time
            
            Trf = app.parameters.pulse_duration_ms * 1e-3;    
            alpha = app.flipangle_degEditField.Value; % Degrees.
            
            % gamma = 4257.7478518;		% Hz/G. 1H
            gamma = app.parameters.gamma; % Hz/G 13C

            %            N = 100; % 5 x T1 for steady state without alpha/2
 %           Ttot = 10 * Trf;
            sint = app.parameters.sint;
            bwfac =  1000/sint;
            % khz_in_g = 4.25755; % 4.25755 KHz =^=  1 G
            khz_in_g = gamma / 1000;

            % time resolution simulation
            dt_sim = app.parameters.pulse_duration_ms_res*1e-3;
            % time course simulation
            t_rf_sim = 0:dt_sim:Trf-dt_sim;

            % time resolution loaded pulse
            dt_b1 = Trf / numel(pulse(:,1));
            % time course loaded pulse
            t_b1 = 0:dt_b1:Trf-dt_b1;
            
            % interpolate loaded pulse onto to simulate pulse:
            pulse_sim = interp1(t_b1, pulse.', t_rf_sim);
            plot(app.UIAxes_rfpulse,1000 * t_rf_sim, app.b1_mag_khz_EditField.Value * abs(pulse_sim))
            app.UIAxes_rfpulse.YAxis.Label.String = 'Exc. B1 [kHz]';
        end


        function [] = plot_anat_image(app, slice, show_exc_profile)
            if sum(sum(app.data.anat_image))
                cla(app.UIAxes_anat_image, 'reset');
                anat_image_slice = squeeze(app.data.anat_image(:,:, slice));
    
                %%%%%%% proton image
                % to remove white part in image
%                 app.UIAxes_anat_image.XLim=  [1, size(app.data.anat_image,1)];
%                 app.UIAxes_anat_image.YLim=  [1, size(app.data.anat_image,2)];
                cla(app.UIAxes_anat_image, 'reset');
                
                % plot image
                imagesc(app.UIAxes_anat_image, anat_image_slice);              
                colormap(app.UIAxes_anat_image, 'bone');
    
                % plot coverage of the animal with the slcie profiles:
                if (app.ShowExcitationProfileCheckBox.Value && (abs(sum(sum(sum(app.data.mxy_b1)))) > 0))
                    % B1 simulated range:
                    y_sim_mm = 10 * linspace(-app.parameters.yrange_cm/2, ...
                                        app.parameters.yrange_cm/2, ...
                                        app.parameters.yrange_cm_npoints);
    
                    % interpolate profile onto anatomical range:
                    y_range_anat = linspace(-app.parameters.anat_image_fov_mm(1)/2, ...
                                             app.parameters.anat_image_fov_mm(1)/2, ...
                                             size(app.data.anat_image,2));
    
                    freq_range_hz = linspace(app.parameters.freq_hz_low, ...
                                    app.parameters.freq_hz_high, ...
                                    app.parameters.freq_hz_npoints); % Hz
       
                    % use "Hz" settings
                    if app.UseFrequenciesCheckBox.Value
                    % show 2 profiles:
                        if app.Show2FreqsCheckBox.Value
                            freq1_hz = app.Freq1EditField.Value;
                            freq2_hz = app.Freq2EditField.Value;
                            
                            % find freqs closest to set frequencies
                            [~,freq1_ind] = min(abs(freq_range_hz-freq1_hz));
                            [~,freq2_ind] = min(abs(freq_range_hz-freq2_hz));
        
                            mxy_b1_interp1 = interp1(y_sim_mm, squeeze(abs(app.data.mxy_b1(:, end, freq1_ind))), y_range_anat);
                            mxy_b1_interp2 = interp1(y_sim_mm, squeeze(abs(app.data.mxy_b1(:, end, freq2_ind))), y_range_anat);
                            % make same intensity
                            mxy_b1_interp2 = mxy_b1_interp2  / max(mxy_b1_interp2 ) * max(mxy_b1_interp1 );
        
                            cla(app.UIAxes_anat_image);

                            % create a "plane that will cover the anatomical image:
                            % exication profile 1
                            exc_prof1 = zeros(size(app.data.anat_image,1), size(app.data.anat_image,2));
                            exc_prof1(round(end/2)+1:end, :) = repmat(mxy_b1_interp1, [round(size(app.data.anat_image,1)/ 2) 1]);
        
                            % exication profile 1
                            exc_pro2 = zeros(size(app.data.anat_image,1), size(app.data.anat_image,2));
                            exc_pro2(1:round(end/2), :) = repmat(mxy_b1_interp2 , [round(size(app.data.anat_image,1)/ 2) 1]);
        
                            profile_im_1 = exc_prof1 + exc_pro2;
                            profile_im_1 = rot90(profile_im_1, -1);
        
                            BG_1 = ind2rgb(im2uint8(profile_im_1 / max(max(profile_im_1)) ),jet(256));

                            anat_image = mat2gray(anat_image_slice / max(max(anat_image_slice)));
                            
                            % blend parameters
                            threshold = 0.0;
                            alpha = 0.6;
                            amap = (anat_image>=threshold)*alpha; % blend based on gray value!
                            blended = im2double(anat_image).*amap + im2double(BG_1).*(1-amap) ;%+ im2double(BG_2).*(1-amap);
                            
                            imagesc(app.UIAxes_anat_image, blended);
        
        
                        else
                            freq1_hz = app.freqSlider.Value;
                            [~,freq1_ind] = min(abs(freq_range_hz-freq1_hz));
                            mxy_b1_interp1 = interp1(y_sim_mm, squeeze(abs(app.data.mxy_b1(:, end, freq1_ind))), y_range_anat);
        
                            cla(app.UIAxes_anat_image);
                            % create a "plane that will cover the anatomical image:
                            profile_im_1 = repmat(mxy_b1_interp1, [size(app.data.anat_image,1) 1]);
                            profile_im_1 = rot90(profile_im_1, -1);


                            BG_1 = ind2rgb(im2uint8(profile_im_1 / max(max(profile_im_1)) ),jet(256));
                            anat_image = mat2gray(anat_image_slice / max(max(anat_image_slice)));
                            
                            % blend parameters
                            threshold = 0.0;
                            alpha = 0.6;
                            amap = (anat_image>=threshold)*alpha; % blend based on gray value!
                            blended = im2double(anat_image).*amap + im2double(BG_1).*(1-amap);
                            
                            imagesc(app.UIAxes_anat_image, blended);
                        end
    
                    % use ppm settings
                    else
                        % translate ppm into Hz
                        freq1_ppm = app.Freq1_ppmEditField.Value;
                        freq2_ppm = app.Freq2_ppmEditField.Value;
                        freqbasic_ppm = app.BasicFreqEditField.Value;
    
        
                        % ppm to hz depending on nucleus:
                        ppm2hz = app.parameters.gamma*7/100;
        
                        freq1_hz = (freqbasic_ppm - freq1_ppm) * ppm2hz;
                        freq2_hz = (freqbasic_ppm - freq2_ppm) * ppm2hz;
    
                        % add "pulse frequency offset"
                        freq1_hz = freq1_hz + app.PulseFreqEditField.Value;
                        freq2_hz = freq2_hz + app.PulseFreqEditField.Value;
    
                        if app.Show2FreqsCheckBox.Value
                            % find freqs closest to set frequencies
                            [~,freq1_ind] = min(abs(freq_range_hz-freq1_hz));
                            [~,freq2_ind] = min(abs(freq_range_hz-freq2_hz));
        
                            mxy_b1_interp1 = interp1(y_sim_mm, squeeze(abs(app.data.mxy_b1(:, end, freq1_ind))), y_range_anat);
                            mxy_b1_interp2 = interp1(y_sim_mm, squeeze(abs(app.data.mxy_b1(:, end, freq2_ind))), y_range_anat);
                            % make same intensity
                            mxy_b1_interp2 = mxy_b1_interp2  / max(mxy_b1_interp2 ) * max(mxy_b1_interp1 );
        
                            cla(app.UIAxes_anat_image);
                            % create a "plane that will cover the anatomical image:
                            exc_prof1 = zeros(size(app.data.anat_image,1), size(app.data.anat_image,2));
                            exc_prof1(round(end/2)+1:end, :) = repmat(mxy_b1_interp1, [round(size(app.data.anat_image,1)/ 2) 1]);
        
                            exc_pro2 = zeros(size(app.data.anat_image,1), size(app.data.anat_image,2));
                            exc_pro2(1:round(end/2), :) = repmat(mxy_b1_interp2 , [round(size(app.data.anat_image,1)/ 2) 1]);
        
                            profile_im_1 = exc_prof1 + exc_pro2;
                            profile_im_1 = rot90(profile_im_1, -1);

               
                            BG_1 = ind2rgb(im2uint8(profile_im_1 / max(max(profile_im_1)) ),jet(256));

                            anat_image = mat2gray(anat_image_slice / max(max(anat_image_slice)));
                            
                            % blend parameters
                            threshold = 0.0;
                            alpha = 0.6;
                            amap = (anat_image>=threshold)*alpha; % blend based on gray value!
                            blended = im2double(anat_image).*amap + im2double(BG_1).*(1-amap) ;%+ im2double(BG_2).*(1-amap);
                            
                            imagesc(app.UIAxes_anat_image, blended);
        
        
                        else
                            [~,freq1_ind] = min(abs(freq_range_hz-freq1_hz));
                            mxy_b1_interp1 = interp1(y_sim_mm, squeeze(abs(app.data.mxy_b1(:, end, freq1_ind))), y_range_anat);
        
                            cla(app.UIAxes_anat_image);
                            % create a "plane that will cover the anatomical image:
                            profile_im_1 = repmat(mxy_b1_interp1, [size(app.data.anat_image,1) 1]);
                            BG_1 = ind2rgb(im2uint8(profile_im_1 / max(max(profile_im_1)) ),jet(256));
                            anat_image = mat2gray(anat_image_slice / max(max(anat_image_slice)));
                            
                            % blend parameters
                            threshold = 0.0;
                            alpha = 0.6;
                            amap = (anat_image>=threshold)*alpha; % blend based on gray value!
                            blended = im2double(anat_image).*amap + im2double(BG_1).*(1-amap);
                            
                            imagesc(app.UIAxes_anat_image, rot90(blended, -1));
                        end
                    end
                end
                if app.overlayrefpowCheckBox.Value
                    % get refpowmap
                    refpow_map = app.data.refpow_map_calc.';
                    % rgb
                    BG_1 = ind2rgb(im2uint8(refpow_map / max(max(refpow_map))), jet(256));
                    % rgb
                    anat_image = ind2rgb(im2uint8(anat_image_slice / max(max(anat_image_slice))), bone(256));
                   
                    % blend parameters
                    threshold = 0.0;
                    alpha = 0.6;
                    amap = (anat_image>=threshold)*alpha; % blend based on gray value!
                    blended = im2double(anat_image).*amap + im2double(BG_1).*(1-amap) ;%+ im2double(BG_2).*(1-amap);
                    
                    imagesc(app.UIAxes_anat_image, blended);
                    
                    pbaspect(app.UIAxes_anat_image, [size(app.data.anat_image,1) size(app.data.anat_image,2) 1]);

                end
                    sz = size(app.data.anat_image);
                    res = app.parameters.anat_image_fov_mm ./ sz(1:2);
                    
                    % set ticks + labels:
                    xtick = linspace(1,size(app.data.anat_image,2),5);
                    ytick = linspace(1,size(app.data.anat_image,1),5);

                    xticklabel = linspace(-app.parameters.anat_image_fov_mm(2)/2+res(2)/2, ...
                                           app.parameters.anat_image_fov_mm(2)/2-res(2)/2, 5);
                    yticklabel = linspace(-app.parameters.anat_image_fov_mm(1)/2+res(1)/2, ...
                                           app.parameters.anat_image_fov_mm(1)/2-res(1)/2, 5);

                    yticklabel = yticklabel  - app.parameters.anat_image_off_mm(1);
                    xticklabel = xticklabel  - app.parameters.anat_image_off_mm(2);

                    app.UIAxes_anat_image.XTick = xtick;
                    app.UIAxes_anat_image.YTick = ytick;
                    
                    app.UIAxes_anat_image.XTickLabel = xticklabel;
                    app.UIAxes_anat_image.YTickLabel = yticklabel;
                    app.UIAxes_anat_image.XLabel.String = 'x [mm]';
                    app.UIAxes_anat_image.YLabel.String = 'y [mm]';
                    grid1 = [xtick;xtick];
                    grid2 = repmat([ytick(1);ytick(end)], 1, length(xtick));
                    hold(app.UIAxes_anat_image, 'on');
                    plot(app.UIAxes_anat_image, grid1, grid2, 'w-', 'LineWidth', 0.025);
                    grid1 = [ytick;ytick];
                    grid2 = repmat([xtick(1);xtick(end)], 1, length(ytick));
                    plot(app.UIAxes_anat_image, grid2, grid1, 'w-', 'LineWidth', 0.025);
                    hold(app.UIAxes_anat_image, 'off');


                    % plot position of the refpow:
                    [x_pos, y_pos] = calc_mmpos_in_voxpos(app, ...
                                    size(app.data.anat_image(:,:,1)), ...
                                    app.parameters.anat_image_fov_mm, ...
                                    app.parameters.refpow_pos_x_mm, ...
                                    app.parameters.refpow_pos_y_mm, ...
                                    app.parameters.anat_image_off_mm(1:2));            
    
                    hold(app.UIAxes_anat_image, 'on');
                    % plot result
                    plot(app.UIAxes_anat_image, x_pos, y_pos, 'rx');
                    hold(app.UIAxes_anat_image, 'off');
                    % keep image at correct ratio
                    pbaspect(app.UIAxes_anat_image, [size(app.data.anat_image,2) size(app.data.anat_image,1) 1]);
                    yyaxis(app.UIAxes_anat_image, 'right')
                    set(app.UIAxes_anat_image, 'ytick', ytick, 'yticklabels',  ...
                        flip(linspace(-app.parameters.anat_image_fov_mm(2)/2+res(2)/2, ...
                                  app.parameters.anat_image_fov_mm(2)/2-res(2)/2, 5) - app.parameters.anat_image_off_mm(2)), 'Ycolor', 'k')
            end


            app.SliceSpinner.Value = slice;
        end


        function [] = plot_b1profile(app)
            % load pulse:
            % get B1:
            try
                b1profile_log = app.data.b1TransmitProfileFit_log;
            end
            if (sum(b1profile_log) == 0)
                b1profile_log = ones(64, 1);
            end

            % parameters
            % time
             % B1 measured range:
            x_b1_cm = linspace(0.25, 31.75, 64)/10;
            x_b1_cm = x_b1_cm - 1.6;

            % B1 simulated range:
            x_sim_cm = linspace(-app.parameters.yrange_cm/2, ...
                                app.parameters.yrange_cm/2, ...
                                app.parameters.yrange_cm_npoints);

            % interpolate
            b1Transmit_sim_log = interp1(x_b1_cm, b1profile_log, x_sim_cm, 'spline', 'extrap');
            b1Transmit_sim = exp(b1Transmit_sim_log) / max(exp(b1Transmit_sim_log ));
            cla(app.UIAxes_b1profile, 'reset');
            plot(app.UIAxes_b1profile,x_sim_cm,b1Transmit_sim);
        end


        function [] = plot_mag_vs_dist(app, freq)
            % load pulse:
            mxy_b1 = app.data.mxy_b1;             
            mz_b1 = app.data.mz_b1;  

            mxy = app.data.mxy;             
            mz = app.data.mz;  
            
            % parameters
            % time
             % B1 measured range:
%             x_b1_cm = linspace(0.25, 31.75, 64)/10;
%             x_b1_cm = x_b1_cm - 1.6;

            % B1 simulated range:
            if (app.parameters.yrange_cm_npoints > 1)
                x_sim_cm = linspace(-app.parameters.yrange_cm/2, ...
                                    app.parameters.yrange_cm/2, ...
                                    app.parameters.yrange_cm_npoints);
            else
                x_sim_cm = app.parameters.yrange_cm;
            end
            freq_range_hz = linspace(app.parameters.freq_hz_low, ...
                                app.parameters.freq_hz_high, ...
                                app.parameters.freq_hz_npoints); % Hz
            freq1_hz = app.Freq1EditField.Value;
            [~,freq1_ind] = min(abs(freq_range_hz-freq));
            x_1 = size(mxy_b1,2);
            try
            if app.UseFrequenciesCheckBox.Value
                % check if 2 frequencies should be plotted:
                if app.Show2FreqsCheckBox.Value
                    freq2_hz = app.Freq2EditField.Value;

                    % find freqs closest to set frequencies
                    [~,freq1_ind] = min(abs(freq_range_hz-freq1_hz));
                    [~,freq2_ind] = min(abs(freq_range_hz-freq2_hz));
                    
                    % plot stuff
                    cla(app.UIAxes_mag_vs_freq);
                    
                    % surface B1
                    l =[];
                    if (app.parameters.plot_inhomo_mx)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(real(mxy_b1(:, x_1 , freq1_ind))), 'k.-');
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(real(mxy_b1(:, x_1, freq2_ind))), 'k.--');
                        l = [l; "Mx (surface)"; ""];
                    end
                    if (app.parameters.plot_inhomo_my)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(imag(mxy_b1(:, x_1, freq1_ind))), 'r-');
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(imag(mxy_b1(:, x_1, freq2_ind))), 'r--');
                        l = [l; "My (surface)"; ""];
                    end
                    if (app.parameters.plot_inhomo_mxy)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(abs(mxy_b1(:, x_1, freq1_ind))), 'k-');
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(abs(mxy_b1(:, x_1, freq2_ind))), 'k--');
                        l = [l; "|Mxy| (surface)"; ""];
                    end
                    if (app.parameters.plot_inhomo_mz)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze((mz_b1(:, x_1, freq1_ind))), 'b-');
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze((mz_b1(:, x_1, freq2_ind))), 'b--');
                        l = [l; "Mz (surface)"; ""];
                    end
                    if (app.parameters.plot_inhomo_phase)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(angle(mxy_b1(:, x_1, freq1_ind))), 'g-');
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(angle(mxy_b1(:, x_1, freq2_ind))), 'g--');
                        l = [l; "\phi (surface)"; ""];
                    end


                    % uniform B1
                    if (app.parameters.plot_homo_mx)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(real(mxy(:, x_1, freq1_ind))), 'k.-');
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(real(mxy(:, x_1, freq2_ind))), 'k.--');
                        l = [l; "Mx (original)"; ""];
                    end
                    if (app.parameters.plot_homo_my)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(imag(mxy(:, x_1, freq1_ind))), 'r-');
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(imag(mxy(:, x_1, freq2_ind))), 'r--');
                        l = [l; "My (original)"; ""];
                    end
                    if (app.parameters.plot_homo_mxy)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(abs(mxy(:, x_1, freq1_ind))), 'k-');
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(abs(mxy(:, x_1, freq2_ind))), 'k--');
                        l = [l; "|Mxy| (original)"; ""];
                    end
                    if (app.parameters.plot_homo_mz)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze((mz(:, x_1, freq1_ind))), 'b-');
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze((mz(:, x_1, freq2_ind))), 'b--');
                        l = [l; "Mz (original)"; ""];
                    end
                    if (app.parameters.plot_homo_phase)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(angle(mxy(:, x_1, freq1_ind))), 'g-');
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(angle(mxy(:, x_1, freq2_ind))), 'g--');
                        l = [l; "\phi (original)"; ""];
                    end
                    
                    app.UIAxes_mag_vs_freq.Title.String = ['f_1 = ' num2str(freq1_hz) 'Hz & f_2 = ' num2str(freq2_hz) 'Hz'];

                    ylims = app.UIAxes_anat_image.YLim;
%                     plot(app.UIAxes_mag_vs_freq,x_sim_cm(app.parameters.find_grad_ind1)*ones(1,10), linspace(ylims(1),ylims(2),10), 'r-');
%                     plot(app.UIAxes_mag_vs_freq,x_sim_cm(app.parameters.find_grad_ind2)*ones(1,10), linspace(ylims(1),ylims(2),10), 'r-');
                    legend(app.UIAxes_mag_vs_freq, l);
                else
                    cla(app.UIAxes_mag_vs_freq);
                    hold(app.UIAxes_mag_vs_freq, 'on');
                    % surface B1
                    l = [];
                    if (app.parameters.plot_inhomo_mx)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(real(mxy_b1(:, x_1, freq1_ind))), 'k.-');
                        l = [l; "Mx (surface)"];
                    end
                    if (app.parameters.plot_inhomo_my)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(imag(mxy_b1(:, x_1, freq1_ind))), 'r-');
                        l = [l; "My (surface); "];
                    end
                    if (app.parameters.plot_inhomo_mxy)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(abs(mxy_b1(:, x_1, freq1_ind))), 'k-');
                        l = [l; "|Mxy| (surface)"];
                    end
                    if (app.parameters.plot_inhomo_mz)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze((mz_b1(:, x_1, freq1_ind))), 'b-');
                        l = [l; "Mz (surface)"];
                    end
                    if (app.parameters.plot_inhomo_phase)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(angle(mxy_b1(:, x_1, freq1_ind))), 'g-');
                        l = [l; "\phi (surface)"];
                    end

                    % uniform B1
                    if (app.parameters.plot_homo_mx)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,max(squeeze(real(mxy_b1(:, x_1, freq1_ind))))*...
                                                                  squeeze(real(mxy(:, x_1, freq1_ind)))/...
                                                              max(squeeze(real(mxy(:, x_1, freq1_ind)))), 'k.-');
                        l = [l; 'Mx (uniform)'];
                    end
                    if (app.parameters.plot_homo_my)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,max(squeeze(imag(mxy_b1(:, x_1, freq1_ind))))*...
                                                                  squeeze(imag(mxy(:, x_1, freq1_ind)))/...
                                                              max(squeeze(imag(mxy(:, x_1, freq1_ind)))), 'r-');
                        l = [l; "Mx (uniform)"];
                    end
                    if (app.parameters.plot_homo_mxy)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,max(squeeze(abs(mxy_b1(:, x_1, freq1_ind))))*...
                                                                  squeeze(abs(mxy(:, x_1, freq1_ind)))/...
                                                              max(squeeze(abs(mxy(:, x_1, freq1_ind)))), 'k-');
                        l = [l; "|Mxy| (uniform)"];
                    end
                    if (app.parameters.plot_homo_mz)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze((mz(:, x_1, freq1_ind))), 'b-');
                        l = [l; "Mz (uniform)"];
                    end
                    if (app.parameters.plot_homo_phase)
                        plot(app.UIAxes_mag_vs_freq,x_sim_cm,squeeze(angle(mxy(:, x_1, freq1_ind))), 'g-');
                        l = [l; "\phi (uniform)"];
                    end
                   app.UIAxes_mag_vs_freq.Title.String = ['f = ' num2str(freq_range_hz(freq1_ind)) 'Hz'];
                    


%                     
%                     plot(app.UIAxes_mag_vs_freq,max(squeeze(abs(mxy_b1(:, x_1, freq)))) * ...
%                                                               squeeze(abs(mxy(:, x_1, freq))) / ...
%                                                           max(squeeze(abs(mxy(:, x_1, freq)))), 'b--');
%                     plot(app.UIAxes_mag_vs_freq,x_sim_cm(app.parameters.find_grad_ind1)*ones(1,10), linspace(0,0.0001,10), 'r-');
%                     plot(app.UIAxes_mag_vs_freq,x_sim_cm(app.parameters.find_grad_ind2)*ones(1,10), linspace(0,0.0001,10), 'r-');
%     
%     
                    
    %                 app.UIAxes_mag_vs_freq
                legend(app.UIAxes_mag_vs_freq, l);
                end
            else % use ppm
                % translate ppm into Hz
                freq1_ppm = app.parameters.freq_ppm1;
                freq2_ppm = app.parameters.freq_ppm2;
                freqbasic_ppm = app.parameters.freq_ppm_basic;

                % ppm to hz depending on nucleus:
                ppm2hz = app.parameters.gamma*7/100;

                freq1_hz = (freqbasic_ppm - freq1_ppm) * ppm2hz;
                freq2_hz = (freqbasic_ppm - freq2_ppm) * ppm2hz;

                % add "pulse freqencies":
                freq1_hz = freq1_hz + app.parameters.freq_hz_pulse;
                freq2_hz = freq2_hz + app.parameters.freq_hz_pulse;

                freq_range_hz = linspace(app.parameters.freq_hz_low, ...
                                app.parameters.freq_hz_high, ...
                                app.parameters.freq_hz_npoints); % Hz
                % find freqs closest to set frequencies
                [~,freq1_ind] = min(abs(freq_range_hz-freq1_hz));
                [~,freq2_ind] = min(abs(freq_range_hz-freq2_hz));
                
                % plot stuff
                cla(app.UIAxes_mag_vs_freq);
                plot(app.UIAxes_mag_vs_freq,squeeze(abs(mxy_b1(:, x_1, freq1_ind))), 'k');
                hold(app.UIAxes_mag_vs_freq, 'on');
                plot(app.UIAxes_mag_vs_freq,squeeze(abs(mxy_b1(:, x_1, freq2_ind))), 'r');
                app.UIAxes_mag_vs_freq.Title.String = ['f_1 = ' num2str(freq1_hz) 'Hz & f_2 = ' num2str(freq2_hz) 'Hz'];
                plot(app.UIAxes_mag_vs_freq,x_sim_cm(app.parameters.find_grad_ind1)*ones(1,10), linspace(0,0.0001,10), 'ro');
                plot(app.UIAxes_mag_vs_freq,x_sim_cm(app.parameters.find_grad_ind2)*ones(1,10), linspace(0,0.0001,10), 'ro');
            end
            catch
                warning('did you change parameters (e.g. frequency range) after running a simulation?');
            end

        end

        
        function results = find_gradient_strength(app)
            % set parameters values to start with
%             app.parameters.gradient_kHz_cm = app.gradientstrength_khzcmEditField.Value;
%             app.parameters.gamma = 4257.7478518; % [1H]
%             app.parameters.T1 = str2double(app.T1_sEditField.Value);
%             app.parameters.T2 = str2double(app.T2_sEditField.Value);
%             % y range:
%             app.parameters.yrange_cm = app.yrange_cmEditField.Value;
%             app.parameters.yrange_cm_npoints = str2double(app.yrange_cm_npointsEditField.Value);
%             % pulse duration
%             app.parameters.pulse_duration_ms = app.pulseduration_msEditField.Value;
%             app.parameters.pulse_duration_ms_npoints = str2double(app.pulseduration_ms_npointsEditField.Value);
%             % frequency [hz]
%             app.parameters.freq_hz_high = app.freq_hz_highEditField.Value;
%             app.parameters.freq_hz_low = app.freq_hz_lowEditField.Value;
%             app.parameters.freq_hz_npoints = str2double(app.freq_hz_npointsEditField.Value);
%             app.parameters.freq_ind = 1;
%             
            %% gradient stuff
            
            
            % grad_kHz_cm_list = (0.25:0.0005:0.3);
            grad_kHz_cm_list = linspace(app.parameters.find_grad_low, ...
                                        app.parameters.find_grad_high, ...
                                        app.parameters.find_grad_npoints);
            
            % grad_Gcm = grad_kHzcm / khz_in_g;
            khz_in_g = app.parameters.gamma / 1000;
            grad_G_cm_list = grad_kHz_cm_list / khz_in_g;  
           
            % range to simulate over
            range_cm_sim = app.parameters.yrange_cm; % 4 cm
            Trf = app.parameters.pulse_duration_ms;
            dt_sim = app.parameters.pulse_duration_ms / app.parameters.pulse_duration_ms_npoints;
            
            
            grad_on_points = round((Trf / dt_sim));
            grad_kHz_cm_range_sim = ones(grad_on_points, 1);
            grad_kHz_cm_range_sim = [grad_kHz_cm_range_sim.' -grad_kHz_cm_range_sim.'/2 zeros(1, 1000)]; % rephasing gradients on
            % grad_kHz_cm_range_sim = [grad_kHz_cm_range_sim.' zeros(1, 1501)]; % rephasing gradients off
            grad_kHz_cm_range_sim = grad_kHz_cm_range_sim .';
            
            
            
            %% pulse stuff
            x_b1 = linspace(0.25, 31.75, 64)/10;
            x_sim = linspace(-range_cm_sim/2, range_cm_sim/2, app.parameters.yrange_cm_npoints);
            
            x_b1 = x_b1 - 1.6;
            
            
            %% b1 stuff
            b1profile_log = app.data.b1TransmitProfileFit_log;
            
            b1Transmit_sim_log = interp1(x_b1, b1profile_log, x_sim, 'spline', 'extrap');
            b1Transmit_sim = exp(b1Transmit_sim_log) / max(exp(b1Transmit_sim_log ));
            
                        
            t_rf_sim = 0:dt_sim:Trf-dt_sim;
            
            
            %% pulse stuff
            pulse = app.data.rfpulse;
            alpha = app.parameters.alpha;
            sint = app.parameters.sint;
            gamma = app.parameters.gamma;
            
            
            dt_pulse = Trf / numel(pulse(:,1));
            t_b1 = 0:dt_pulse:Trf-dt_pulse;
            pulse_sim = interp1(t_b1, pulse.', t_rf_sim);
            
            
            khz_in_g = gamma / 1000;
            b1_amp_khz  = app.parameters.b1_amp_khz;
            b1_amp_G = b1_amp_khz / khz_in_g;
            
            % get pulse in Gauss
            pulse_sim = pulse_sim * b1_amp_G;
            pulse_sim = [pulse_sim zeros(1, numel(pulse_sim)+1000)];
            pulse_sim = pulse_sim .';
            
            
            
            
            %% parameter stuff:
            freq = app.parameters.freq_hz_range;
            T1 = app.parameters.T1;
            T2 = app.parameters.T2;
            
            
            %% homogeneous B1
            % mxy = zeros(numel(grad_G_cm_list), numel(x_sim), numel(freq));
            % mz = zeros(numel(grad_G_cm_list), numel(x_sim), numel(freq));
            % gc = 0;
            % for g = grad_G_cm_list
            %     gc = gc+1;
            %     [mxx, myx, mzx] = bloch_13c_no_output(pulse_sim, g* grad_kHz_cm_range_sim, dt_sim, T1, T2, freq, x_sim, 2);
            %     mxy(gc, :, :) = mxx(end, :, :) + 1i* myx(end, :, :);
            %     mz(gc, :, :) = mzx(end, :, :);
            % end
            % 
            % 
            % 
            % figure();
            % for i = 1:11
            %     subplot(3,4,i);
            %     plot(abs(squeeze(mxy(:, :, i))).')
            % end
            
            
            pulse_sim_short = interp1(linspace(0,1,numel(pulse_sim)), pulse_sim, linspace(0,1,1000));
            grad_kHz_cm_range_sim_short = interp1(linspace(0,1,numel(pulse_sim)), grad_kHz_cm_range_sim, linspace(0,1,1000));
            dt_sim_short = numel(pulse_sim) / 1000 * dt_sim;

            app.find_GradButton.Text = 'searching for best Grad!';
            app.find_GradButton.FontAngle = 'italic';
            app.find_GradButton.BackgroundColor = [0.00,0.00,1.00];

            pause(0.000001);
            
            
            
            %% inhomogeneous B1
            mxy_x = zeros(numel(grad_G_cm_list), numel(x_sim), numel(pulse_sim_short), numel(freq));
            mz_x = zeros(numel(grad_G_cm_list), numel(x_sim), numel(pulse_sim_short), numel(freq));
            gc = 0;
            if sum(app.CButton.Value) % 
                for g = grad_G_cm_list
                    gc = gc+1;
                    for i = 1:numel(x_sim)
                        [mxx, myx, mzx] = bloch_13c_no_output(pulse_sim_short * b1Transmit_sim(i), g* grad_kHz_cm_range_sim_short, dt_sim_short/1000, T1, T2, freq, x_sim(i), 2);
                        mxy_x(gc, i, :, :) = mxx + 1i* myx;
                        mz_x(gc, i, :, :) = mzx;
                    end
                    pause(0.00000001);
                    app.progress_grad_findLabel.Text = [num2str(gc) '/' num2str(numel(grad_G_cm_list))];
                end
            else
                for g = grad_G_cm_list
                gc = gc+1;
                for i = 1:numel(x_sim)
                    [mxx, myx, mzx] = bloch_no_output(pulse_sim_short * b1Transmit_sim(i), g* grad_kHz_cm_range_sim_short, dt_sim_short/1000, T1, T2, freq, x_sim(i), 2);
                    mxy_x(gc, i, :, :) = mxx + 1i* myx;
                    mz_x(gc, i, :, :) = mzx;
                end
                pause(0.00000001);
                app.progress_grad_findLabel.Text = [num2str(gc) '/' num2str(numel(grad_G_cm_list))];
                end
            end
            app.find_GradButton.Text = 'find Grad!';
            app.find_GradButton.FontAngle = 'normal';
            app.find_GradButton.BackgroundColor = [0.96,0.96,0.96];

            % 
            % figure();
            % for i = 1:20
            %     subplot(4,5,i);
            %     plot(x_sim, abs(squeeze(mxy_x(:, :, end,i))).')
            %     title(['f = ' num2str(freq(i)) 'Hz']);
            % end
            figure();
            subplot(2,1,1);
            imagesc(squeeze(abs(mxy_x(:, :,end,round(end/2)))));
            hold on;
            plot(app.parameters.find_grad_ind1 * ones(1, numel(grad_G_cm_list)),  (1:numel(grad_G_cm_list)),'r-')
            plot(app.parameters.find_grad_ind2 * ones(1, numel(grad_G_cm_list)),  (1:numel(grad_G_cm_list)),'r-')

            
            
            find_plat = squeeze(mxy_x(:, app.parameters.find_grad_ind1:app.parameters.find_grad_ind2, end, round(end/2)));
            find_plat_diff = zeros(numel(grad_G_cm_list), numel(app.parameters.find_grad_ind1:app.parameters.find_grad_ind2) - 1);
            for i = 1:numel(grad_G_cm_list)
               find_plat_diff(i, :) = diff(squeeze(abs(find_plat(i, :))));
            end
            
            find_plat_diff = sum(find_plat_diff, 2);
            
            [~, find_plat_min] = min(abs(find_plat_diff));
            
           
            
            grad_G_cm_best = grad_G_cm_list(find_plat_min);
            grad_kHz_cm_best = grad_kHz_cm_list(find_plat_min);


            subplot(2,1,2);
            plot(grad_kHz_cm_list, abs(find_plat_diff));
            xlabel('kHz/cm')
            title(['Best Grad: ' num2str(grad_G_cm_best) 'G/cm -- ' num2str(grad_kHz_cm_best) 'kHz/cm']);
        
            app.BestGradEditField.Value = grad_kHz_cm_best;   
        end
        
        function [] = simulate_repetitions(app)
            try
            % resolution simulation:
            dt_sim = app.parameters.pulse_duration_ms_res*1e-3;
            trf = app.parameters.pulse_duration_ms*1e-3;
            trep = app.parameters.repetitions_time_ms * 1e-3;
            tref = app.parameters.ref_grad_dur_ms * 1e-3;
            
            trep_npoints = app.parameters.repetitions_time_ms_npoints;
            N = app.parameters.nrepetitions;
            % calculate how many additional round (preparation pulse,
            % tipback, spoiling) need to be simulated:
            nadd = 0;
            % preparation pulse:
            if (app.parameters.first_pulse_amp ~= 1)
                nadd = nadd + 1;
            end
            % tipback
            if strcmp((app.parameters.tipback_pulse), "on")
                nadd = nadd + 1;
            end
            if app.parameters.spoiler_grad
                nadd = nadd + 1;
            end
            N = N + nadd;
            app.parameters.Nsims = N;

            % if tipbackpulse, then do 1 more repetition:
            if strcmp((app.parameters.tipback_pulse), "on")
                % skip for now
%                 N = N + 1;
            end
                        
            T1= app.parameters.T1;
            T2= app.parameters.T2;
            
            %% GRADIENT STUFF
            % number of points
            dt_sim_rep = trep / trep_npoints;

            % if the duration of the rephasing gradient lobe is defined:
            if (tref > 0)
                % number of pointer the pulse is playing ("excitation gradient"):
                grad_exc_points = round(trf/dt_sim_rep);
                % number of points the refocus gradients are on:
                grad_ref_points = round(tref/dt_sim_rep);
                % number of point nothing is on:
                grad_off_points = round(trep/dt_sim_rep) - grad_exc_points - 2*grad_ref_points;
                % nubmer of points any gradient is on:
                grad_on_points = grad_exc_points + 2 * grad_ref_points;

                % to calculate the strength of the refocusing gradients:
                grad_exc_area = grad_exc_points;
                grad_ref_area = grad_ref_points;
                grad_reph_amp = grad_exc_area / grad_ref_area /2 ;
    
                % init empty array
                grad_ss_sim = zeros(1, trep_npoints);
                % fill the refocusing gradient (first half)
                grad_ss_sim(round(grad_off_points/2)+1:round(grad_off_points/2)+round(grad_ref_points)) = -grad_reph_amp*ones(1, grad_ref_points);
                % fill the excitation gradient
                grad_ss_sim(round(grad_off_points/2)+round(grad_ref_points)+1:...
                            round(grad_off_points/2)+round(grad_ref_points)+grad_exc_points) = ...
                            1*ones(1, grad_exc_points);
                % fill the refocusing gradient (second half)
                grad_ss_sim(round(grad_off_points/2)+grad_ref_points+grad_exc_points +1 : ...
                            round(grad_off_points/2)+grad_ref_points+grad_exc_points +grad_ref_points) = ...
                            -grad_reph_amp*ones(1, grad_ref_points);
            else
                % points of RF pulse
                grad_exc_points = round(trf/dt_sim_rep);
                % refocusing points (1 lobe)
                grad_ref_points = round((trep_npoints - grad_exc_points) /2);
                grad_off_points = grad_ref_points ;
                % ratio:
                amp_reph =  grad_exc_points / 2 / grad_ref_points ;

                grad_ss_sim = zeros(trep_npoints, 1);
                grad_ss_sim(1:grad_ref_points) = -amp_reph;
                grad_ss_sim(grad_ref_points+1:grad_ref_points+grad_exc_points) = 1;
                grad_ss_sim(grad_ref_points+grad_exc_points+1:end) = -amp_reph;
            end
%             grad_ss_sim(1:)


%             if (grad_off_points > 0)

%                 % gradient on:
%                 grad_ss_sim = ones(1, round(grad_on_points));
%                 % add gradients (rephase) + zeros
%                 grad_ss_sim = [-grad_ss_sim/2 grad_ss_sim -grad_ss_sim/2 zeros(1, grad_off_points)].'; % rephasing gradients on
%                 
%                 % test if everything adds up to 0 (for proper rephasing)
%                 grad_sum = sum(grad_ss_sim)
            
%             else % have to change duration of rephasing gradients
%             grad_points_left = round(trep/dt_sim_rep - trf/dt_sim_rep);
%             % gradient on during slice selection:
%             grad_ss_sim = ones(1, round(grad_exc_points));
%             % calc area of gradient on:
%             grad_ss_sim_area = sum(grad_ss_sim);
%             % 1 rephasing gradient shoul be 1/2 slice selective gradient (area wise):
%             grad_reph_sim_area = grad_ss_sim_area ;
%             % calc amplitude of rephase gradient:
%             grad_reph_amp = 0.5*grad_reph_sim_area / round(grad_points_left/2);
%             %
%             grad_reph_sim = grad_reph_amp * ones(1, round(grad_points_left/2)); 



            % test 
            grad_sum = sum(grad_ss_sim);
            if (abs(grad_sum) > 0)
                warning(['gradients not balanced: Sum = ' num2str(grad_sum)]);
                % for now continue anyway
                grad_sim = grad_ss_sim;
            else
                grad_sim = grad_ss_sim;
            end



            %% PULSE STUFF
            pulse = app.data.rfpulse;
            alpha = app.parameters.alpha;
            sint = app.parameters.sint;
            gamma = app.parameters.gamma;
            
            dt_pulse = trf / numel(pulse(:,1));
            
            % pulse timing (loaded pulse size scaled to rf pulse duration)
            t_b1 = 0:dt_pulse:trf-dt_pulse;
            % goal of pulse time points
            t_rf_sim = 0:dt_sim_rep:trf-dt_sim_rep;

            pulse_sim = interp1(t_b1, pulse.', t_rf_sim);

            %
            if (tref > 0)
                pulse_range_sim = zeros(size(grad_sim));
                pulse_range_sim(round(grad_off_points/2)+round(grad_ref_points)+1:...
                round(grad_off_points/2)+round(grad_ref_points)+grad_exc_points) = pulse_sim;
            else
                pulse_range_sim = zeros(size(grad_sim));
                pulse_range_sim(grad_ref_points+1:grad_ref_points+grad_exc_points) = pulse_sim;
            end

            %% REDUCE NUMBER OF POINTS TO SPEED UP SIMULATION:
            pulse_sim_short = interp1(linspace(0,1,numel(pulse_range_sim)), pulse_range_sim, linspace(0,1,trep_npoints));
            grad_sim_short = interp1(linspace(0,1,numel(grad_sim)), grad_sim, linspace(0,1,trep_npoints));
            dt_sim_short = trep / numel(grad_sim_short);


            khz_in_g = gamma / 1000;
            b1_amp_khz  = app.parameters.b1_amp_khz;
            b1_amp_G = b1_amp_khz / khz_in_g;
            

            

            % scale gradients:
            grad_amp_kHz_cm = app.parameters.gradient_kHz_cm;
            grad_amp_G_cm  = grad_amp_kHz_cm / khz_in_g;
            grad_range_sim_short_G_cm = grad_sim_short * grad_amp_G_cm;
            
            

            %% RUN SIMULATION
            if (app.parameters.yrange_cm_npoints > 1)
                x_sim = linspace(-app.parameters.yrange_cm/2, ...
                                    app.parameters.yrange_cm/2, ...
                                    app.parameters.yrange_cm_npoints);
            else
                x_sim = app.parameters.yrange_cm;
            end
            

            freq = linspace(app.parameters.freq_hz_low, ...
                            app.parameters.freq_hz_high, ...
                            app.parameters.freq_hz_npoints); 
            
            % frequency offset of the 1st pulse
            freq_offset_1stpulse = app.parameters.freq_1stfreq_hz;

            % frequenct offset of the 2nd pulse:
            freq_offset_2ndpulse = app.parameters.freq_2ndfreq_hz;
            % how often one alternates between 2 frequencies
            n_alternations = app.parameters.alternate_freq_Ntimes_2nd_freq;

            % check if high res sim. data sould be stored:
            if strcmp(app.parameters.save_all_sim_data, 'on')
%                 ByteSize(zeros(N*n_alternations, numel(x_sim), numel(pulse_sim_short), numel(freq)))
                mx_x_b1 = zeros(N*n_alternations, numel(x_sim), numel(pulse_sim_short), numel(freq));
                my_x_b1 = zeros(N*n_alternations, numel(x_sim), numel(pulse_sim_short), numel(freq));
                mz_x_b1 = app.parameters.sample_mag * ones(N*n_alternations, numel(x_sim), numel(pulse_sim_short), numel(freq));
            else
                mx_x_b1 = zeros(N*n_alternations, numel(x_sim), numel(freq));
                my_x_b1 = zeros(N*n_alternations, numel(x_sim), numel(freq));
                mz_x_b1 = app.parameters.sample_mag * ones(N*n_alternations, numel(x_sim), numel(freq));
            end
                % check if 2nd pulse is played):
            if strcmp(app.parameters.simulate_alternating_freqs, 'on')
                freqoff = repmat([freq_offset_1stpulse freq_offset_2ndpulse], [1 n_alternations]);
                b1_amp_G_2ndfreq = (pi * app.parameters.alpha_2ndfreq / 180 ) * 1/((2*pi*gamma*1e4) * trf* sint) * 1e4;% G
                b1_amp_G = repmat([b1_amp_G b1_amp_G_2ndfreq ], [1 n_alternations]);
            else
                freqoff = freq_offset_1stpulse;
            end
            b1Transmit_sim_log = interp1(linspace(0.25, 31.75, 64)/10 - 1.6, ...
                                app.data.b1TransmitProfileFit_log, ...
                                x_sim, "spline", 'extrap');
            b1Transmit_sim_log(isnan(b1Transmit_sim_log)) = 0;
            b1Transmit_sim = exp(b1Transmit_sim_log) / max(exp(b1Transmit_sim_log));
            

            app.SimulateRepsButton.Text = 'simu. reptd. excts!';
            app.SimulateRepsButton.FontAngle = 'italic';
            app.SimulateRepsButton.BackgroundColor = [0.02,0.67,0.01];

            pause(0.00000001);
            catch
                app.repetition_time_npointsEditField.BackgroundColor = 'r';
                warning('probably wrong amount of points, try even number of point.');
                return;
            end

            % get pulse in Gauss
%             pulse_sim_short = pulse_sim_short * b1_amp_G;

            % alternating phase from 1st to second RF Pulse
            phase_diff_pulses = app.parameters.phase_pulse_deg;
            phase_diff_pulses = repmat([0 phase_diff_pulses], [1 n_alternations]);
            % alternating phase between N Repetitions with one frequency
            phase_2nd_pulse =  app.parameters.alternating_phase;
            phase_2nd_pulse = repmat([0 phase_2nd_pulse], [1 n_alternations*N]);
            % factor of alpha prep pulse
            first_pulse_amp = app.parameters.first_pulse_amp;
            % calc the first pulse flipangle (degree)
            fa_first_pulse = repmat(first_pulse_amp * [alpha app.parameters.alpha_2ndfreq], [1 n_alternations]);
            % calc B1 of first pulse:
            results = calc_b1(app, fa_first_pulse);
            % calc ratio of B1 of first pulse to the following pulses:
            first_pulse_amp = results.b1_amp_g ./ b1_amp_G;

%             if app.parameters.alternating_phase
            % calculate the tipbackpulse:
            tipback_pulse_phase=0;
            if strcmp((app.parameters.tipback_pulse), "on")
                tipback_pulse_amp = 0.5;
                tipback_pulse_phase = app.parameters.phase_tipback_pulse;
                %tipback_pulse_amp = 0.3882; % in case of fa = 4 degree and 
                % effective flipangle = 3.1055 @ f=245 (pyruvate)
            end

            % if first pulse amplitude should be played after each
            % repetition
            if app.parameters.first_pulse_always
                first_pulse_amp_nreps = 1; % // divide first amp by 1 --> keeps first amp for the follwing 
                % first pulses
            else
                first_pulse_amp_nreps = first_pulse_amp(1); % // divide first amp by the first amp factor -->
                % following first pulses have alpha amplitude
            end

            if (numel(x_sim) == 1)
                n_par_workers = 0;
            else
                p = gcp();
                if isempty(p)
                    n_par_workers = 0;
                else
                    n_par_workers = p.NumWorkers;
                end
            end
            % give information
            try
                disp(['Running Simulation; expected size:']);
                ByteSize(mx_x_b1)
                % https://stackoverflow.com/questions/4845561/how-to-know-the-size-of-a-variable-in-matlab
            end
            disp(['First TR / TR: ' num2str(app.parameters.first_repetition_time_ms) 'ms / ' num2str(trep*1e3) 'ms']);
            disp(['TRF: ' num2str(trf*1e3) 'ms']);
            disp(['NR per pulse: ' num2str(N)]);
            disp(['# alternations between pulses: ' num2str(n_alternations)]);
            disp(['Excitation frequencies: ' num2str(freqoff) 'Hz']);
            disp(['Phase between 2 diff freq. pulses: ' num2str(phase_2nd_pulse) '°']);
            disp(['Excitation amplitudes: ' num2str(b1_amp_G) 'G']);
            disp(['First Excitation amplitude: ' num2str(first_pulse_amp*100) '%']);
            disp(['Flipangles: ' num2str(app.parameters.alpha) '° / '  num2str(app.parameters.alpha_2ndfreq) '°']);
            disp(['Gradient Strength: ' num2str(grad_amp_G_cm) 'G/cm']);
            disp(['Magnetization: ' num2str(app.parameters.sample_mag) 'G/cm']);
            disp(['T1 / T2: ' num2str(T1) ' / ' num2str(T2) 's']);
            try
                disp(['Simulated range / res: ' num2str(abs(x_sim(end)-x_sim(1))) 'cm / ' num2str(abs(x_sim(2)-x_sim(1))) 'cm']);
            catch
                disp(['Simulated range / res: ' num2str(abs(x_sim)) 'cm / ' num2str(abs(x_sim)) 'cm']);
            end
            disp(['Simulated freq. range / res: ' num2str(abs(freq(end)-freq(1))) 'Hz / ' num2str(abs(freq(2)-freq(1))) 'Hz']);


            %% transform everything into 13C to avoid double coding for both 1H and 13C 
            


        if (sum(isnan(pulse_sim_short)) > 20)
            warning('more than 20 NaN in pulse_sim_short');
                figure()
                plot(abs(pulse_sim_short));
                title(['pulse_sim_short with ' num2str(sum(isnan(pulse_sim_short))) ' NaNs']);
        end
        pulse_sim_short(isnan(pulse_sim_short)) = 0;
        tic

    if (app.parameters.first_repetition_time_ms == ...
        app.parameters.repetitions_time_ms)
        ir = 0;
       
        % use B1 Field instead of homogeneous B1
    %     else
    %         if sum(app.CButton.Value) % 
                % alternating excitation frequency:
        for i_pulse = 1:n_alternations 
            pulse_amp = b1_amp_G(i_pulse);
            freq_offset = freqoff(i_pulse);
            fpamp = first_pulse_amp(i_pulse);
            % repetitions per excitation frequency:
            for i_rep = 1:N
                ir=ir+1;
                phase_rep = phase_2nd_pulse(i_rep) + phase_diff_pulses(i_pulse);
                % very first pulse:
                if ((i_rep == 1) && (i_pulse == 1))
                    v1 = squeeze(mx_x_b1(1, :, :));
                    v2 = squeeze(my_x_b1(1, :, :));
                    v3 = squeeze(mz_x_b1(1, :, :));
    
                    parfor (k = 1:numel(x_sim), n_par_workers)
    % Explanation:
    %  
    % [mxx, myx, mzx] = //output         
    % bloch_13c_no_output // mex file (compiled .c)
    % (first_pulse_amp * ... // amplitude  of first pulse
    % (exp(2*pi*1i*ones(1,numel(pulse_sim_short )) .* ir*phase /360)) .* ...// phase of pulse (might vary) (was removed)
    % pulse_sim_short * ... // pulse shape
    % b1Transmit_sim(i) .*... // B1 field of coil (depends on x = distance to coil, magnitude = 1)
    % exp(- 1i * 2 * pi * linspace(0, trep, trep_npoints) * freq_offset) ... // frequency offset of pulse
    % * pulse_amp .* // pulse amplitude (kHz)
    % exp(2*pi*1i*ones(1,numel(pulse_sim_short)) * phase_rep/360), ... // phase of pulse
    % grad_range_sim_short_G_cm, ... // Gradients
    % dt_sim_short // time resolution
    % , T1, T2, freq, x_sim(i), 2, ... // should be clear
    % v1(i, :), ... // magnetization mx from previous excitation
    % v2(i, :), ... // magnetization my from previous excitation
    % v3(i, :));    // magnetization mz from previous excitation
                        [mxx, myx, mzx] = bloch_13c_no_output(fpamp .* pulse_sim_short * b1Transmit_sim(k) .*...
                                                    exp(- 1i * 2 * pi * linspace(0, trep, trep_npoints) * freq_offset) * pulse_amp .* exp(2*pi*1i*ones(1,numel(pulse_sim_short)) * phase_rep/360), ...
                                                    grad_range_sim_short_G_cm, ...
                                                    dt_sim_short, T1, T2, freq, x_sim(k), 2, ...
                                                    v1(k, :), ...
                                                    v2(k, :), ...
                                                    v3(k, :));
                                                    
                        a(k, :) = mxx(end,:);
                        b(k, :) = myx(end,:);
                        c(k, :) = mzx(end,:);                                
                    end
                % first pulse after alternating frequency (starts with
                % a different pulse amplitude (first_pulse_amp))
                elseif ((i_rep == 1) && (i_pulse ~= 1))
                    v1 = squeeze(mx_x_b1((i_pulse-1)*N, :, :));
                    v2 = squeeze(my_x_b1((i_pulse-1)*N, :, :));
                    v3 = squeeze(mz_x_b1((i_pulse-1)*N, :, :));
    
                    parfor (k = 1:numel(x_sim), n_par_workers)
                        [mxx, myx, mzx] = bloch_13c_no_output(fpamp / first_pulse_amp_nreps .* pulse_sim_short * b1Transmit_sim(k) .*...
                                                    exp(- 1i * 2 * pi * linspace(0, trep, trep_npoints) * freq_offset) * pulse_amp .* exp(2*pi*1i*ones(1,numel(pulse_sim_short)) * phase_rep/360), ...
                                                    grad_range_sim_short_G_cm, ...
                                                    dt_sim_short, T1, T2, freq, x_sim(k), 2, ...
                                                    v1(k, :), ...
                                                    v2(k, :), ...
                                                    v3(k, :));
                                                    
                        a(k, :) = mxx(end,:);
                        b(k, :) = myx(end,:);
                        c(k, :) = mzx(end,:);                                
                    end

                % 2nd last pulse (if a tipback is requested and spoiler )
                elseif ((i_rep == N-1) && strcmp((app.parameters.tipback_pulse), "on") && app.parameters.spoiler_grad)
                    v1 = squeeze(mx_x_b1(ir-1, :, :));
                    v2 = squeeze(my_x_b1(ir-1, :, :));
                    v3 = squeeze(mz_x_b1(ir-1, :, :));
    
                    parfor (k = 1:numel(x_sim), n_par_workers)
                        [mxx, myx, mzx] = bloch_13c_no_output(tipback_pulse_amp  .* pulse_sim_short * b1Transmit_sim(k) .*...
                                                    exp(- 1i * 2 * pi * linspace(0, trep, trep_npoints) * freq_offset) * pulse_amp .* exp(2*pi*1i*ones(1,numel(pulse_sim_short)) * (tipback_pulse_phase+ phase_rep)/360), ...
                                                    grad_range_sim_short_G_cm, ...
                                                    dt_sim_short, T1, T2, freq, x_sim(k), 2, ...
                                                    v1(k, :), ...
                                                    v2(k, :), ...
                                                    v3(k, :));
                                                    
                        a(k, :) = mxx(end,:);
                        b(k, :) = myx(end,:);
                        c(k, :) = mzx(end,:);                                
                    end

                % last pulse (if a tipback is requested and no spoiler )
                elseif ((i_rep == N) && strcmp((app.parameters.tipback_pulse), "on") && ~app.parameters.spoiler_grad)
                    v1 = squeeze(mx_x_b1(ir-1, :, :));
                    v2 = squeeze(my_x_b1(ir-1, :, :));
                    v3 = squeeze(mz_x_b1(ir-1, :, :));
    
                    parfor (k = 1:numel(x_sim), n_par_workers)
                        [mxx, myx, mzx] = bloch_13c_no_output(tipback_pulse_amp  .* pulse_sim_short * b1Transmit_sim(k) .*...
                                                    exp(- 1i * 2 * pi * linspace(0, trep, trep_npoints) * freq_offset) * pulse_amp .* exp(2*pi*1i*ones(1,numel(pulse_sim_short)) * (tipback_pulse_phase+ phase_rep)/360), ...
                                                    grad_range_sim_short_G_cm, ...
                                                    dt_sim_short, T1, T2, freq, x_sim(k), 2, ...
                                                    v1(k, :), ...
                                                    v2(k, :), ...
                                                    v3(k, :));
                                                    
                        a(k, :) = mxx(end,:);
                        b(k, :) = myx(end,:);
                        c(k, :) = mzx(end,:);                                
                    end
                
                

                % was tipped back and will be spoiled --> dont pulse, just
                % use previously simulated data (mx & will be set to 0)
                elseif ((i_rep == N) && strcmp((app.parameters.tipback_pulse), "on") &&(app.parameters.spoiler_grad))
                    a = squeeze(mx_x_b1(ir-1, :, :));
                    b = squeeze(my_x_b1(ir-1, :, :));
                    c = squeeze(mz_x_b1(ir-1, :, :));

                % all other repetitions
                else
                    v1 = squeeze(mx_x_b1(ir-1, :, :));
                    v2 = squeeze(my_x_b1(ir-1, :, :));
                    v3 = squeeze(mz_x_b1(ir-1, :, :));
                    parfor (k = 1:numel(x_sim), n_par_workers)
                       [mxx, myx, mzx] = bloch_13c_no_output(pulse_sim_short * b1Transmit_sim(k) .*...
                                                   exp(- 1i * 2 * pi * linspace(0, trep, trep_npoints) * freq_offset) * pulse_amp .* exp(2*pi*1i*ones(1,numel(pulse_sim_short)) * phase_rep/360), ...
                                                   grad_range_sim_short_G_cm, ...
                                                   dt_sim_short, T1, T2, freq, x_sim(k), 2, ...
                                                   v1(k, :), ...
                                                   v2(k, :), ...
                                                   v3(k, :));
                        a(k, :) = mxx(end,:);
                        b(k, :) = myx(end,:);
                        c(k, :) = mzx(end,:);
                    end
        
        
                end
                mx_x_b1(ir, :, :) = a;
                my_x_b1(ir, :, :) = b;
                mz_x_b1(ir, :, :) = c;
                
                %                         progress_simRep(app, ir, N);
                app.progress_sim_repLabel.Text = [num2str(ir) '/' num2str(N*n_alternations)];
                pause(0.000000000001);
            end
            % "Spoil mxy magnetization"
            if (app.parameters.spoiler_grad)
                mx_x_b1(ir, :, :) = zeros(size(a));
                my_x_b1(ir, :, :) = zeros(size(a));
            end
        end
    % different repetition time (first vs following)
    else
        first_trep = app.parameters.first_repetition_time_ms;
        first_trep_npoints = first_trep * 1e-3/ dt_sim_rep;
        diff_points = round((trep - first_trep*1e-3) / dt_sim_rep);
    
        %% pulse:
        first_pulse = zeros(1, 2*trep_npoints);
        first_pulse(1:trep_npoints) = pulse_sim_short;
        % give first pulse amplitude here:
        first_pulse = first_pulse * first_pulse_amp(1); % might not be 100% precise
        
        first_pulse_2ndpulse = zeros(1, 2*trep_npoints);
        first_pulse_2ndpulse(trep_npoints+1:end) = pulse_sim_short;
    
        first_pulse = circshift(first_pulse, diff_points);
        first_pulse  = first_pulse + first_pulse_2ndpulse;
        %% gradient:
        first_grad = zeros(1, 2*trep_npoints);
        first_grad(1:trep_npoints) = grad_range_sim_short_G_cm;
        
        first_grad_2ndpulse = zeros(1, 2*trep_npoints);
        first_grad_2ndpulse(trep_npoints+1:end) = grad_range_sim_short_G_cm;
    
        first_grad = circshift(first_grad, diff_points);
        first_grad  = first_grad + first_grad_2ndpulse;
    
    
        %% split:
        second_pulse = first_pulse(trep_npoints+1:end);
        first_pulse = first_pulse(1:trep_npoints);
    
        second_grad = first_grad(trep_npoints+1:end);
        first_grad = first_grad(1:trep_npoints);
    
        
        if (sum(isnan(pulse_sim_short)) > 20)
            warning('more than 20 NaN in pulse_sim_short');
                figure()
                plot(abs(pulse_sim_short));
                title(['pulse_sim_short with ' num2str(sum(isnan(pulse_sim_short))) ' NaNs']);
        end
        pulse_sim_short(isnan(pulse_sim_short)) = 0;
         if (sum(isnan(first_pulse)) > 20)
            warning('more than 20 NaN in first_pulse');
                figure()
                plot(abs(first_pulse));
                title(['first_pulse with ' num2str(sum(isnan(first_pulse))) ' NaNs']);
        end
        first_pulse(isnan(first_pulse)) = 0;
        if (sum(isnan(second_pulse)) > 20)
            warning('more than 20 NaN in second pulse');
            figure()
            plot(abs(second_pulse));
            title(['Second pulse with ' num2str(sum(isnan(second_pulse))) ' NaNs']);
        end
        second_pulse(isnan(second_pulse)) = 0;

        ir = 0;
    %     if sum(app.CButton.Value) % 
                % alternating excitation frequency:
                % needs to be pre defined for parfor:
        for i_pulse = 1:n_alternations 
            pulse_amp = b1_amp_G(i_pulse);
            freq_offset = freqoff(i_pulse);
            % repetitions per excitation frequency:
            for i_rep = 1:N
                phase_rep = phase_2nd_pulse(i_rep) + phase_diff_pulses(i_pulse);
                ir=ir+1;
                % very first pulse, first half:
                if ((i_rep == 1) && (i_pulse == 1))
                    v1 = squeeze(mx_x_b1(1, :, :));
                    v2 = squeeze(my_x_b1(1, :, :));
                    v3 = squeeze(mz_x_b1(1, :, :));
    
                    parfor (k = 1:numel(x_sim), n_par_workers)
                        [mxx, myx, mzx] = bloch_13c_no_output(first_pulse * b1Transmit_sim(k) .*...
                                                    exp(- 1i * 2 * pi * linspace(0, trep, trep_npoints) * freq_offset) * pulse_amp.* exp(2*pi*1i*ones(1,numel(pulse_sim_short)) * phase_rep/360), ...
                                                    first_grad, ...
                                                    dt_sim_short, T1, T2, freq, x_sim(k), 2, ...
                                                    v1(k, :), ...
                                                    v2(k, :), ...
                                                    v3(k, :));
                                                    
                        a(k, :) = mxx(end,:);
                        b(k, :) = myx(end,:);
                        c(k, :) = mzx(end,:);                                
                    end
                % very first pulse, second half
                elseif ((i_rep == 2) && (i_pulse == 1))
                    % to avoid changing the phase of the 2nd pulse aswell:
                    phase_term = zeros(size(ones(1,numel(pulse_sim_short))));
                    % phase of second pulse:
                    phase_rep_p2 = phase_2nd_pulse(i_rep) + phase_diff_pulses(i_pulse);
                    % phase of first pulse:
                    phase_rep_p1 = phase_2nd_pulse(i_rep-1) + phase_diff_pulses(i_pulse);
                    
                    % find first point that is zero (not best way but q&d)
                    ind = find(abs(second_pulse)==0);

                    
    
                    % up to first point where the this pulse = 0 (should
                    % mean the first pulse finished: use the phase from
                    % first_pulse
                    phase_term(1:ind(1)) = phase_rep_p1;
                    % afterwards apply the (potential) phase shift for the
                    % second pulse:
                    phase_term(ind(1)+1:end) = phase_rep_p2;

                    % previous magnetization:
                    v1 = squeeze(mx_x_b1(1, :, :));
                    v2 = squeeze(my_x_b1(1, :, :));
                    v3 = squeeze(mz_x_b1(1, :, :));
    
                    parfor (k = 1:numel(x_sim), n_par_workers)
                        [mxx, myx, mzx] = bloch_13c_no_output(second_pulse * b1Transmit_sim(k) .*...
                                                    exp(- 1i * 2 * pi * linspace(0, trep, trep_npoints) * freq_offset) * pulse_amp.* exp(2*pi*1i*ones(1,numel(pulse_sim_short)) .* phase_term/360), ...
                                                    second_grad, ...
                                                    dt_sim_short, T1, T2, freq, x_sim(k), 2, ...
                                                    v1(k, :), ...
                                                    v2(k, :), ...
                                                    v3(k, :));
                                                    
                        a(k, :) = mxx(end,:);
                        b(k, :) = myx(end,:);
                        c(k, :) = mzx(end,:);                                
                    end
                % first pulse of every repetition (besides the very first), first half:
                elseif ((i_rep == 1) && (i_pulse ~= 1))
                    v1 = squeeze(mx_x_b1((i_pulse-1)*N, :, :));
                    v2 = squeeze(my_x_b1((i_pulse-1)*N, :, :));
                    v3 = squeeze(mz_x_b1((i_pulse-1)*N, :, :));
    
                    parfor (k = 1:numel(x_sim), n_par_workers)
                        [mxx, myx, mzx] = bloch_13c_no_output(first_pulse / first_pulse_amp_nreps * b1Transmit_sim(k) .*...
                                                    exp(- 1i * 2 * pi * linspace(0, trep, trep_npoints) * freq_offset) * pulse_amp.* exp(2*pi*1i*ones(1,numel(pulse_sim_short)) * phase_rep/360), ...
                                                    first_grad, ...
                                                    dt_sim_short, T1, T2, freq, x_sim(k), 2, ...
                                                    v1(k, :), ...
                                                    v2(k, :), ...
                                                    v3(k, :));
                                                    
                        a(k, :) = mxx(end,:);
                        b(k, :) = myx(end,:);
                        c(k, :) = mzx(end,:);                                
                    end
                % first pulse of every repetition (besides the very first), second half:
                 elseif ((i_rep == 2) && (i_pulse ~= 1))
                    v1 = squeeze(mx_x_b1((i_pulse-1)*N+1, :, :));
                    v2 = squeeze(my_x_b1((i_pulse-1)*N+1, :, :));
                    v3 = squeeze(mz_x_b1((i_pulse-1)*N+1, :, :));

                    % to avoid changing the phase of the 2nd pulse aswell:
                    phase_term = zeros(size(ones(1,numel(pulse_sim_short))));
                    % phase of second pulse:
                    phase_rep_p2 = phase_2nd_pulse(i_rep) + phase_diff_pulses(i_pulse);
                    % phase of first pulse:
                    phase_rep_p1 = phase_2nd_pulse(i_rep-1) + phase_diff_pulses(i_pulse);
                    
                    % find first point that is zero (not best way but q&d)
                    ind = find(abs(second_pulse)==0);

                    % up to first point where the this pulse = 0 (should
                    % mean the first pulse finished: use the phase from
                    % first_pulse
                    phase_term(1:ind(1)) = phase_rep_p1;
                    % afterwards apply the (potential) phase shift for the
                    % second pulse:
                    phase_term(ind(1)+1:end) = phase_rep_p2;
    
                    parfor (k = 1:numel(x_sim), n_par_workers)
                        [mxx, myx, mzx] = bloch_13c_no_output(second_pulse / first_pulse_amp_nreps * b1Transmit_sim(k) .*...
                                                    exp(- 1i * 2 * pi * linspace(0, trep, trep_npoints) * freq_offset) * pulse_amp.* exp(2*pi*1i*ones(1,numel(pulse_sim_short)) .* phase_term/360), ...
                                                    second_grad, ...
                                                    dt_sim_short, T1, T2, freq, x_sim(k), 2, ...
                                                    v1(k, :), ...
                                                    v2(k, :), ...
                                                    v3(k, :));
                                                    
                        a(k, :) = mxx(end,:);
                        b(k, :) = myx(end,:);
                        c(k, :) = mzx(end,:);                                
                    end
                % 2nd last pulse (if a tipback is requested and spoiler )
                elseif ((i_rep == N-1) && strcmp((app.parameters.tipback_pulse), "on") && app.parameters.spoiler_grad)
                    v1 = squeeze(mx_x_b1(ir-1, :, :));
                    v2 = squeeze(my_x_b1(ir-1, :, :));
                    v3 = squeeze(mz_x_b1(ir-1, :, :));
    
                    parfor (k = 1:numel(x_sim), n_par_workers)
                        [mxx, myx, mzx] = bloch_13c_no_output(tipback_pulse_amp  .* pulse_sim_short * b1Transmit_sim(k) .*...
                                                    exp(- 1i * 2 * pi * linspace(0, trep, trep_npoints) * freq_offset) * pulse_amp .* exp(2*pi*1i*ones(1,numel(pulse_sim_short)) * (tipback_pulse_phase+ phase_rep)/360), ...
                                                    grad_range_sim_short_G_cm, ...
                                                    dt_sim_short, T1, T2, freq, x_sim(k), 2, ...
                                                    v1(k, :), ...
                                                    v2(k, :), ...
                                                    v3(k, :));
                                                    
                        a(k, :) = mxx(end,:);
                        b(k, :) = myx(end,:);
                        c(k, :) = mzx(end,:);                                
                    end

                % last pulse (if a tipback is requested and no spoiler )
                elseif ((i_rep == N) && strcmp((app.parameters.tipback_pulse), "on") && ~app.parameters.spoiler_grad)
                    v1 = squeeze(mx_x_b1(ir-1, :, :));
                    v2 = squeeze(my_x_b1(ir-1, :, :));
                    v3 = squeeze(mz_x_b1(ir-1, :, :));
    
                    parfor (k = 1:numel(x_sim), n_par_workers)
                        [mxx, myx, mzx] = bloch_13c_no_output(tipback_pulse_amp  .* pulse_sim_short * b1Transmit_sim(k) .*...
                                                    exp(- 1i * 2 * pi * linspace(0, trep, trep_npoints) * freq_offset) * pulse_amp .* exp(2*pi*1i*ones(1,numel(pulse_sim_short)) * (tipback_pulse_phase+ phase_rep)/360), ...
                                                    grad_range_sim_short_G_cm, ...
                                                    dt_sim_short, T1, T2, freq, x_sim(k), 2, ...
                                                    v1(k, :), ...
                                                    v2(k, :), ...
                                                    v3(k, :));
                                                    
                        a(k, :) = mxx(end,:);
                        b(k, :) = myx(end,:);
                        c(k, :) = mzx(end,:);                                
                    end

%                 elseif ((i_rep == N) && strcmp((app.parameters.tipback_pulse), "on") &&(app.parameters.spoiler_grad))
%                     a = squeeze(mx_x_b1(ir-1, :, :));
%                     b = squeeze(my_x_b1(ir-1, :, :));
%                     c = squeeze(mz_x_b1(ir-1, :, :));
% 
%                 % keep magnetization from previos (otherwise last pulse
%                 % will flip the magnetization again. this should be
%                 % improved in the future:
%                 elseif ((i_rep == N) && strcmp((app.parameters.tipback_pulse), "on"))
%                     a = squeeze(mx_x_b1(ir-1, :, :));
%                     b = squeeze(my_x_b1(ir-1, :, :));
%                     c = squeeze(mz_x_b1(ir-1, :, :));
                

                % all other repetitions
                else
                    v1 = squeeze(mx_x_b1(ir-1, :, :));
                    v2 = squeeze(my_x_b1(ir-1, :, :));
                    v3 = squeeze(mz_x_b1(ir-1, :, :));
                    parfor (k = 1:numel(x_sim), n_par_workers)
                       [mxx, myx, mzx] = bloch_13c_no_output(pulse_sim_short * b1Transmit_sim(k) .*...
                                                   exp(- 1i * 2 * pi * linspace(0, trep, trep_npoints) * freq_offset) * pulse_amp.* exp(2*pi*1i*ones(1,numel(pulse_sim_short)) * phase_rep/360), ...
                                                   grad_range_sim_short_G_cm, ...
                                                   dt_sim_short, T1, T2, freq, x_sim(k), 2, ...
                                                   v1(k, :), ...
                                                   v2(k, :), ...
                                                   v3(k, :));
                        a(k, :) = mxx(end,:);
                        b(k, :) = myx(end,:);
                        c(k, :) = mzx(end,:);
                    end
                end
                mx_x_b1(ir, :, :) = a;
                my_x_b1(ir, :, :) = b;
                mz_x_b1(ir, :, :) = c;
                
                %                         progress_simRep(app, ir, N);
                app.progress_sim_repLabel.Text = [num2str(ir) '/' num2str(N*n_alternations)];
                pause(0.000000000001);
            end
            % "Spoil mxy magnetization"
            if (app.parameters.spoiler_grad)
                mx_x_b1(ir, :, :) = zeros(size(a));
                my_x_b1(ir, :, :) = zeros(size(b));
            end
        end
    end
          
            app.SimulateRepsButton.Text = 'Simulate Reps!';
            app.SimulateRepsButton.FontAngle = 'normal';
            app.SimulateRepsButton.BackgroundColor = [0.96,0.96,0.96];
            toc;

            mxy_x_b1 = mx_x_b1 + 1i * my_x_b1;
    
            if strcmp(app.parameters.save_all_sim_data, 'on')
                mxy = permute(mxy_x_b1, [3 1 2 4]);
                mz = permute(mz_x_b1, [3 1 2 4]);
                
                sz = size(mxy);

                mxy = reshape(mxy, [sz(1)*sz(2) sz(3) sz(4)]);
                mz= reshape(mz, [sz(1)*sz(2) sz(3) sz(4)]);
                pause();
            else
                
            end
            try
                app.parameters.x_sim = x_sim;
                app.data.mxy_b1_rep = mxy_x_b1;
                app.data.mz_b1_rep = mz_x_b1;

                app.parameters.grad_range_sim_short_G_cm = grad_range_sim_short_G_cm;
                app.parameters.pulse_sim_short = pulse_sim_short;
                app.parameters.b1Transmit_sim = b1Transmit_sim;
                app.parameters.dt_sim_short = dt_sim_short;
            catch
                warning('problems storing data in the app struct, paused!')
                pause();
            end
                if app.simrangeCheckBox.Value
                    save_data(app);
                    
%                     save()
                else
                end
            plot_bssfp_signal(app);

            disp('');

        end

        function results = load_rfpulse(app)
            results = 0;
            try
                rfpulse = load(app.parameters.rfpulse_file);
            catch
                warning('RF pulse could no be loaded');
                app.RFPulseEditField.Value = '';
                app.parameters.rfpulse_file_name = '';
                app.parameters.rfpulse_file_path = '';
            end
            
            % save loaded RF pulse:
            [~, dims] = size(rfpulse);
            try
                if (dims == 2)
                    rfpulse = rfpulse(:,1) .*exp(1i*pi * rfpulse(:, 2) / 180);
                    rfpulse = rfpulse / max(abs(rfpulse));
                    app.data.rfpulse = rfpulse;
                elseif (dims == 4) % probably from RF Pulse Wizard
                    % first column amplitude [au], third column phase
                    % [degree]
                    rfpulse = rfpulse(:,1) .*exp(1i*pi * rfpulse(:, 3) / 180);
                    rfpulse = rfpulse / max(abs(rfpulse));
                    app.data.rfpulse = rfpulse;
                else
                    try
                        rfpulse = rfpulse / max(abs(rfpulse));
                        app.data.rfpulse = rfpulse.';
                    % maybe it's a struct?
                    catch
                        rfpulse = rfpulse.rfpulse;
                        rfpulse = rfpulse / max(abs(rfpulse));
                        app.data.rfpulse = rfpulse.';
                    end
                end
                results = 1;
            catch
                warning('dont know what to do...');
                results = 0;
            end
            
        end
        
        function [] = calc_slice_thick(app)
            bw_hz = app.parameters.bwfac_hz_ms / app.parameters.pulse_duration_ms;
            slice_thick_cm = bw_hz / (app.parameters.gradient_kHz_cm *1e3);
            app.slice_thick.Text = [num2str(slice_thick_cm) ' cm'];
        end
        
        function [] = find_rfpower(app)
            % calc RF power map:
            disp('');
        end
        
        
        
        function results = calc_b1map(app)
            % calculate the b1 map according to the coordinates of the
            % anatomical image:
            if sum(sum(sum(app.data.anat_image))) > 0
                fov_anat = app.parameters.anat_image_fov_mm;
                fov_refpow_map = app.parameters.refpow_map_fov_mm;
                refpow_map = app.data.refpow_map;
            else
            end
            
            disp('');
            
        end
        
        function [] = calc_refpow_anat(app)
            %% calc grid of refpowmap:
            % size reference power map
            s_mat_refpow_map = size(app.data.refpow_map);
            % resolution
            fov_refpow_map_mm = app.parameters.refpow_map_fov_mm;
            res_refpow_map_mm = fov_refpow_map_mm ./ s_mat_refpow_map;


            % size reference power map
            s_mat_anat = size(app.data.anat_image);
            % resolution
            fov_anat_mm = app.parameters.anat_image_fov_mm;
            res_anat_mm = fov_anat_mm ./ s_mat_anat(1:2);
            % offset:
            off_anat_mm = app.parameters.anat_image_off_mm;


            % centered grid (refpow):
%             [y1,~] = ndgrid(linspace((-fov_refpow_map_mm(1)+res_refpow_map_mm(1))/2, ...
%                                       (fov_refpow_map_mm(1)-res_refpow_map_mm(1))/2, ...
%                                        s_mat_refpow_map(1)));
            [x1,y1] = meshgrid(linspace((-fov_refpow_map_mm(1)+res_refpow_map_mm(1))/2, ...
                                      (fov_refpow_map_mm(1)-res_refpow_map_mm(1))/2, ...
                                       s_mat_refpow_map(1)),...
                                       linspace((-fov_refpow_map_mm(2)+res_refpow_map_mm(2))/2, ...
                                      (fov_refpow_map_mm(2)-res_refpow_map_mm(2))/2, ...
                                       s_mat_refpow_map(2)));

            % centered grid (anatomical)
            [x2,y2] = meshgrid(linspace((-fov_anat_mm(1)+res_anat_mm(1))/2, ...
                                      (fov_anat_mm(1)-res_anat_mm(1))/2, ...
                                       s_mat_anat(1)), ...
                                       linspace((-fov_anat_mm(2)+res_anat_mm(2))/2, ...
                                      (fov_anat_mm(2)-res_anat_mm(2))/2, ...
                                       s_mat_anat(2)));

            % coordinates with offset
            x2 = x2 + off_anat_mm(1);
            y2 = y2 + off_anat_mm(2);

            app.data.refpow_map_calc = griddata(x1, y1, app.data.refpow_map, x2, y2);
        end
            

        function [] = plot_bssfp_vs_freq(app, freq)
            % load pulse:
            mxy = app.data.mxy_b1_rep;             
            mz = app.data.mz_b1_rep;  

            mxy_b1 = app.data.mxy_b1_rep;             
            mz_b1 = app.data.mz_b1_rep;  

%             mxy = app.data.mxy;             
%             mz = app.data.mz;  
            
            % parameters
            % time
             % B1 measured range:
%             x_b1_cm = linspace(0.25, 31.75, 64)/10;
%             x_b1_cm = x_b1_cm - 1.6;

            % B1 simulated range:
            x_sim_cm = linspace(-app.parameters.yrange_cm/2, ...
                                app.parameters.yrange_cm/2, ...
                                app.parameters.yrange_cm_npoints);
            freq_range_hz = linspace(app.parameters.freq_hz_low, ...
                                app.parameters.freq_hz_high, ...
                                app.parameters.freq_hz_npoints); % Hz
            freq1_hz = app.Freq1EditField.Value;
            [~,freq1_ind] = min(abs(freq_range_hz-freq));
%             if app.UseFrequenciesCheckBox.Value
                % check if 2 frequencies should be plotted:
                x_1 = 1;
            %     
            if app.Show2FreqsCheckBox.Value
                freq2_hz = app.Freq2EditField.Value;

                % find freqs closest to set frequencies
                [~,freq1_ind] = min(abs(freq_range_hz-freq1_hz));
                [~,freq2_ind] = min(abs(freq_range_hz-freq2_hz));
                
                % plot stuff
                cla(app.UIAxes_mag_vs_freq);
                
                % surface B1
                l =[];
                if (app.parameters.plot_inhomo_mx)
                    plot(app.UIAxes_mag_vs_freq,squeeze(real(mxy_b1(:, x_1, freq1_ind))), 'k.-');
                    plot(app.UIAxes_mag_vs_freq,squeeze(real(mxy_b1(:, x_1, freq2_ind))), 'k.--');
                    l = [l; "Mx (surface)"; ""];
                end
                if (app.parameters.plot_inhomo_my)
                    plot(app.UIAxes_mag_vs_freq,squeeze(imag(mxy_b1(:, x_1, freq1_ind))), 'r-');
                    plot(app.UIAxes_mag_vs_freq,squeeze(imag(mxy_b1(:, x_1, freq2_ind))), 'r--');
                    l = [l; "My (surface)"; ""];
                end
                if (app.parameters.plot_inhomo_mxy)
                    plot(app.UIAxes_mag_vs_freq,squeeze(abs(mxy_b1(:, x_1, freq1_ind))), 'k-');
                    plot(app.UIAxes_mag_vs_freq,squeeze(abs(mxy_b1(:, x_1, freq2_ind))), 'k--');
                    l = [l; "|Mxy| (surface)"; ""];
                end
                if (app.parameters.plot_inhomo_mz)
                    plot(app.UIAxes_mag_vs_freq,squeeze((mz_b1(:, x_1, freq1_ind))), 'b-');
                    plot(app.UIAxes_mag_vs_freq,squeeze((mz_b1(:, x_1, freq2_ind))), 'b--');
                    l = [l; "Mz (surface)"; ""];
                end
                if (app.parameters.plot_inhomo_phase)
                    plot(app.UIAxes_mag_vs_freq,squeeze(angle(mxy_b1(:, x_1, freq1_ind))), 'g-');
                    plot(app.UIAxes_mag_vs_freq,squeeze(angle(mxy_b1(:, x_1, freq2_ind))), 'g--');
                    l = [l; "\phi (surface)"; ""];
                end


                % uniform B1
                if (app.parameters.plot_homo_mx)
                    plot(app.UIAxes_mag_vs_freq,squeeze(real(mxy(:, x_1, freq1_ind))), 'k.-');
                    plot(app.UIAxes_mag_vs_freq,squeeze(real(mxy(:, x_1, freq2_ind))), 'k.--');
                    l = [l; "Mx (original)"; ""];
                end
                if (app.parameters.plot_homo_my)
                    plot(app.UIAxes_mag_vs_freq,squeeze(imag(mxy(:, x_1, freq1_ind))), 'r-');
                    plot(app.UIAxes_mag_vs_freq,squeeze(imag(mxy(:, x_1, freq2_ind))), 'r--');
                    l = [l; "My (original)"; ""];
                end
                if (app.parameters.plot_homo_mxy)
                    plot(app.UIAxes_mag_vs_freq,squeeze(abs(mxy(:, x_1, freq1_ind))), 'b-');
                    plot(app.UIAxes_mag_vs_freq,squeeze(abs(mxy(:, x_1, freq2_ind))), 'b--');
                    l = [l; "|Mxy| (original)"; ""];
                end
                if (app.parameters.plot_homo_mz)
                    plot(app.UIAxes_mag_vs_freq,squeeze((mz(:, x_1, freq1_ind))), 'b-');
                    plot(app.UIAxes_mag_vs_freq,squeeze((mz(:, x_1, freq2_ind))), 'bp--');
                    l = [l; "Mz (original)"; ""];
                end
                if (app.parameters.plot_homo_phase)
                    plot(app.UIAxes_mag_vs_freq,squeeze(angle(mxy(:, x_1, freq1_ind))), 'g-');
                    plot(app.UIAxes_mag_vs_freq,squeeze(angle(mxy(:, x_1, freq2_ind))), 'bh--');
                    l = [l; "\phi (original)"; ""];
                end
                
                app.UIAxes_mag_vs_freq.Title.String = ['f_1 = ' num2str(freq1_hz) 'Hz & f_2 = ' num2str(freq2_hz) 'Hz'];

                ylims = app.UIAxes_anat_image.YLim;
%                     plot(app.UIAxes_mag_vs_freq,x_sim_cm(app.parameters.find_grad_ind1)*ones(1,10), linspace(ylims(1),ylims(2),10), 'r-');
%                     plot(app.UIAxes_mag_vs_freq,x_sim_cm(app.parameters.find_grad_ind2)*ones(1,10), linspace(ylims(1),ylims(2),10), 'r-');
                legend(app.UIAxes_mag_vs_freq, l);
            else
                cla(app.UIAxes_mag_vs_freq);
                hold(app.UIAxes_mag_vs_freq, 'on');
                % surface B1
                l = [];
                if (app.parameters.plot_inhomo_mx)
                    plot(app.UIAxes_mag_vs_freq,squeeze(real(mxy_b1(:, x_1, freq1_ind))), 'k.-');
                    l = [l; "Mx (surface)"];
                end
                if (app.parameters.plot_inhomo_my)
                    plot(app.UIAxes_mag_vs_freq,squeeze(imag(mxy_b1(:, x_1, freq1_ind))), 'r-');
                    l = [l; "My (surface); "];
                end
                if (app.parameters.plot_inhomo_mxy)
                    plot(app.UIAxes_mag_vs_freq,squeeze(abs(mxy_b1(:, x_1, freq1_ind))), 'k-');
                    l = [l; "|Mxy| (surface)"];
                end
                if (app.parameters.plot_inhomo_mz)
                    plot(app.UIAxes_mag_vs_freq,squeeze((mz_b1(:, x_1, freq1_ind))), 'b-');
                    l = [l; "Mz (surface)"];
                end
                if (app.parameters.plot_inhomo_phase)
                    plot(app.UIAxes_mag_vs_freq,squeeze(angle(mxy_b1(:, x_1, freq1_ind))), 'g-');
                    l = [l; "\phi (surface)"];
                end

                % uniform B1
                if (app.parameters.plot_homo_mx)
                    plot(app.UIAxes_mag_vs_freq,max(squeeze(real(mxy_b1(:, x_1, freq1_ind))))*...
                                                              squeeze(real(mxy(:, x_1, freq1_ind)))/...
                                                          max(squeeze(real(mxy(:, x_1, freq1_ind)))), 'k.-');
                    l = [l; 'Mx (uniform)'];
                end
                if (app.parameters.plot_homo_my)
                    plot(app.UIAxes_mag_vs_freq,max(squeeze(imag(mxy_b1(:, x_1, freq1_ind))))*...
                                                              squeeze(imag(mxy(:, x_1, freq1_ind)))/...
                                                          max(squeeze(imag(mxy(:, x_1, freq1_ind)))), 'r-');
                    l = [l; "Mx (uniform)"];
                end
                if (app.parameters.plot_homo_mxy)
                    plot(app.UIAxes_mag_vs_freq,max(squeeze(abs(mxy_b1(:, x_1, freq1_ind))))*...
                                                              squeeze(abs(mxy(:, x_1, freq1_ind)))/...
                                                          max(squeeze(abs(mxy(:, x_1, freq1_ind)))), 'b-');
                    l = [l; "|Mxy| (uniform)"];
                end
                if (app.parameters.plot_homo_mz)
                    plot(app.UIAxes_mag_vs_freq,squeeze((mz(:, x_1, freq1_ind))), 'b-');
                    l = [l; "Mz (uniform)"];
                end
                if (app.parameters.plot_homo_phase)
                    plot(app.UIAxes_mag_vs_freq,squeeze(angle(mxy(:, x_1, freq1_ind))), 'g-');
                    l = [l; "\phi (uniform)"];
                end
               app.UIAxes_mag_vs_freq.Title.String = ['f = ' num2str(freq_range_hz(freq1_ind)) 'Hz'];
               app.UIAxes_mag_vs_freq.XLabel.String = 'Repetitions';
                


%                     
%                     plot(app.UIAxes_mag_vs_freq,max(squeeze(abs(mxy_b1(:, x_1, freq)))) * ...
%                                                               squeeze(abs(mxy(:, x_1, freq))) / ...
%                                                           max(squeeze(abs(mxy(:, x_1, freq)))), 'b--');
%                     plot(app.UIAxes_mag_vs_freq,x_sim_cm(app.parameters.find_grad_ind1)*ones(1,10), linspace(0,0.0001,10), 'r-');
%                     plot(app.UIAxes_mag_vs_freq,x_sim_cm(app.parameters.find_grad_ind2)*ones(1,10), linspace(0,0.0001,10), 'r-');
%     
%     
                
%                 app.UIAxes_mag_vs_freq
            legend(app.UIAxes_mag_vs_freq, l);
            end
        end

        function [] = plot_bssfp_vs_rep(app, rep)
            cla(app.UIAxes_anat_image);
            mxy = app.data.mxy_b1_rep;             
            mz = app.data.mz_b1_rep;  

            mxy_b1 = app.data.mxy_b1_rep;             
            mz_b1 = app.data.mz_b1_rep;  

%             mxy = app.data.mxy;             
%             mz = app.data.mz;  
            
            % parameters
            % time
             % B1 measured range:
%             x_b1_cm = linspace(0.25, 31.75, 64)/10;
%             x_b1_cm = x_b1_cm - 1.6;

            % B1 simulated range:
            x_sim_cm = linspace(-app.parameters.yrange_cm/2, ...
                                app.parameters.yrange_cm/2, ...
                                app.parameters.yrange_cm_npoints);
            freq_range_hz = linspace(app.parameters.freq_hz_low, ...
                                app.parameters.freq_hz_high, ...
                                app.parameters.freq_hz_npoints); % Hz
            l = [];
            x_1 = 1;
            if app.parameters.plot_homo_mx
                plot(app.UIAxes_anat_image, freq_range_hz,squeeze(real(mxy(rep, x_1, :))), 'k.-');
                l = [l; "Mx "; ""];
            end
            if (app.parameters.plot_homo_my)
                plot(app.UIAxes_anat_image, freq_range_hz,squeeze(imag(mxy(rep, x_1, :))), 'r-');
                l = [l; "My "; ""];
            end
            if (app.parameters.plot_homo_mxy)
                plot(app.UIAxes_anat_image, freq_range_hz,squeeze(abs(mxy(rep, x_1, :))), 'b-');
                l = [l; "|Mxy| "; ""];
            end
            if (app.parameters.plot_homo_mz)
                plot(app.UIAxes_anat_image, freq_range_hz,squeeze((mz(rep, x_1, :))), 'b-');
                l = [l; "Mz "; ""];
            end
            if (app.parameters.plot_homo_phase)
                plot(app.UIAxes_anat_image, freq_range_hz,squeeze(angle(mxy(rep, x_1, :))), 'g-');
                l = [l; "\phi"; ""];
            end
            
            app.UIAxes_anat_image.Title.String = ['rep = ' num2str(rep)];
            app.UIAxes_anat_image.XLabel.String = 'f [Hz]';
            ylims = app.UIAxes_anat_image.YLim;
%                     plot(app.UIAxes_mag_vs_freq,x_sim_cm(app.parameters.find_grad_ind1)*ones(1,10), linspace(ylims(1),ylims(2),10), 'r-');
%                     plot(app.UIAxes_mag_vs_freq,x_sim_cm(app.parameters.find_grad_ind2)*ones(1,10), linspace(ylims(1),ylims(2),10), 'r-');
            legend(app.UIAxes_anat_image, l); 

        end

        function results = plot_b1map(app)
            disp('');
            cla(app.UIAxes_b1profile);
            
            
            % to remove white part in image
            app.UIAxes_b1profile.XLim=  [1, size(app.data.refpow_map,1)];
            app.UIAxes_b1profile.YLim=  [1, size(app.data.refpow_map,2)];
            % plot image
            imagesc(app.UIAxes_b1profile, app.data.refpow_map);
%             imagesc(app.UIAxes_b1profile, rot90(anat_image_slice, -1));
            % set ticks + labels:
            sz = size(app.data.refpow_map);
            res = app.parameters.refpow_map_fov_mm ./ sz(1:2);

            xtick = linspace(1,size(app.data.refpow_map,1),5);
            ytick = linspace(1,size(app.data.refpow_map,2),5);

            xticklabel = linspace(-app.parameters.refpow_map_fov_mm(2)/2 + res(2)/2, ...
                                   app.parameters.refpow_map_fov_mm(2)/2 - res(2)/2, 5);
            yticklabel = linspace(-app.parameters.refpow_map_fov_mm(1)/2 + res(1)/2, ...
                                   app.parameters.refpow_map_fov_mm(1)/2 - res(1)/2, 5);

            yticklabel = yticklabel  - app.parameters.refpow_map_off_mm(1);
%             pause();
            % wrong x/y refpowmap offset are wrong!
            xticklabel = xticklabel  - app.parameters.refpow_map_off_mm(2); %anat_image_off_mm(1);

            app.UIAxes_b1profile.XTick = xtick;
            app.UIAxes_b1profile.YTick = ytick;

            app.UIAxes_b1profile.XTickLabel = xticklabel;
            app.UIAxes_b1profile.YTickLabel = yticklabel;
            app.UIAxes_b1profile.XLabel.String = 'x [mm]';
            app.UIAxes_b1profile.YLabel.String = 'y [mm]';

            % keep image at correct ratio
            pbaspect(app.UIAxes_b1profile, [size(app.data.refpow_map,1) size(app.data.refpow_map,2) 1]);
            
            % plot point where ref pow sould be calulated:
            hold(app.UIAxes_b1profile, 'on');
            
            % calc position of voxel (in mm) in voxels:  matrix_size, fov, x_pos_mm, y_pos_mm, offset_mm)
            [x_pos, y_pos] = calc_mmpos_in_voxpos(app, ...
                                    size(app.data.refpow_map), ...
                                    app.parameters.refpow_map_fov_mm, ...
                                    app.parameters.refpow_pos_x_mm, ...
                                    app.parameters.refpow_pos_y_mm, ...
                                    app.parameters.refpow_map_off_mm);
        
            
            % plot result
            plot(app.UIAxes_b1profile, x_pos, y_pos, 'rx');

            % add colorbar
            colorbar(app.UIAxes_b1profile);
            app.refpow_map_caxis.Visible = 'on';
            app.refpow_map_caxis.Limits = [0 10];

            if (app.refpow_map_caxis.Value < 1e-6)
                app.refpow_map_caxis.Value = 1e-6;
            end
            caxis(app.UIAxes_b1profile, [0 app.refpow_map_caxis.Value]);
            
            % Grid
%             grid(app.UIAxes_b1profile, 'on');
%             app.UIAxes_b1profile.GridColor = 'w';
            grid1 = [xtick;xtick];
            grid2 = repmat([ytick(1);ytick(end)], 1, length(xtick));
            plot(app.UIAxes_b1profile, grid1, grid2, 'w-', 'LineWidth', 0.025);
            plot(app.UIAxes_b1profile, grid2, grid1, 'w-', 'LineWidth', 0.025);

        end

        function [] = plot_bssfp_signal(app)
            figure()
            subplot(2,3,1);
            imagesc(squeeze(real(app.data.mxy_b1_rep(:,:,round(end/2)))));
            title(['real - f = ' num2str(app.parameters.freq_hz_range(round(end/2))) 'Hz']);
            colorbar();
            subplot(2,3,2);
            plot(app.parameters.freq_hz_range, abs(squeeze((app.data.mxy_b1_rep(end,round(end/2),:)))));
            title(['spec - y = ' num2str(app.parameters.x_sim(round(end/2))) 'cm']);
            xlabel('f [Hz]')
            
            subplot(2,3,3);
            imagesc(squeeze(abs(app.data.mxy_b1_rep(:,:,round(end/2)))));
            title(['mag - f = ' num2str(app.parameters.freq_hz_range(round(end/2))) 'Hz']);
            colorbar();
            subplot(2,3,4);
            plot(squeeze(real(app.data.mxy_b1_rep(:,round(end/2),round(end/2)))));
            hold on;
            plot(squeeze(imag(app.data.mxy_b1_rep(:,round(end/2),round(end/2)))));
            plot(squeeze(abs(app.data.mxy_b1_rep(:,round(end/2),round(end/2)))));
            title(['f = ' num2str(app.parameters.freq_hz_range(round(end/2))) 'Hz']);
            xlabel('Repetitions');
            subplot(2,3,5);
            plot((app.parameters.repetitions_time_ms:app.parameters.repetitions_time_ms:app.parameters.repetitions_time_ms*numel(squeeze(app.data.mz_b1_rep(:,round(end/2),round(end/2)))))*1e-3, ...
                squeeze((app.data.mz_b1_rep(:,round(end/2),round(end/2)))));
            title(['f = ' num2str(app.parameters.freq_hz_range(round(end/2))) 'Hz']);
            xlabel('t [s]');
            subplot(2,3,6);
            plot((app.parameters.dt_sim_short:app.parameters.dt_sim_short:app.parameters.dt_sim_short*numel(app.parameters.grad_range_sim_short_G_cm))*1e-3, app.parameters.grad_range_sim_short_G_cm);
            hold on;
            plot((app.parameters.dt_sim_short:app.parameters.dt_sim_short:app.parameters.dt_sim_short*numel(app.parameters.grad_range_sim_short_G_cm))*1e-3, abs(app.parameters.pulse_sim_short));
            title('pulse sequence');
            xlabel('t [s]');
        end

        function [] = calc_refpow(app)
            % calc position of voxel (in mm) in voxels:
            [x_pos, y_pos] = calc_mmpos_in_voxpos(app, ...
                                size(app.data.refpow_map), ...
                                app.parameters.refpow_map_fov_mm, ...
                                app.parameters.refpow_pos_x_mm, ...
                                app.parameters.refpow_pos_y_mm, ...
                                app.parameters.refpow_map_off_mm);

            if app.parameters.calc_refpow_outside_map
                % try interpolation along y:
                    refpow_y_measured = app.data.refpow_map(:, round(x_pos));
%                     refpow_y_interp = exp(interp1(linspace(-app.parameters.refpow_map_fov_mm(1)/2, ...
%                                                             app.parameters.refpow_map_fov_mm(1)/2, 3;)
            else
                refpow = app.data.refpow_map(round(y_pos), round(x_pos));
            end
            app.RefpowLabel.Text = [num2str(refpow) 'W'];
        end
        
        function [x_pos, y_pos] = calc_mmpos_in_voxpos(app, matrix_size, fov, x_pos_mm, y_pos_mm, offset_mm)
            % calc position of voxel (in mm) in voxels:
            % size reference power map
            % resolution
            res_mm = fov ./ matrix_size ;
            
            % calc positions:
            % image center in voxel:
            center_vox = matrix_size/2;

            % position of pointer in voxel (x)
            x_pos_vox = x_pos_mm /res_mm(2);
            % offset in voxel:
            x_off_vox = offset_mm(2)/res_mm(2);


            % position of pointer in voxel (x)
            y_pos_vox = y_pos_mm /res_mm(1);
            % offset in voxel:
            y_off_vox = offset_mm(1)/res_mm(1);


            x_pos = center_vox(2) + x_pos_vox  + x_off_vox ;
            y_pos = center_vox(1) + y_pos_vox  + y_off_vox ;
            

            % calc positions + taking into account the offset of the
            % anatomical image:
             
        end

        
        % save data
        function success = save_data(app)
            try
            % check which save options are ticked:
            save_options = app.saveoptions_ListBox.Value;
            save_opts = zeros(1, numel(app.saveoptions_ListBox.Items));
            save_struct = [];
            if sum(strcmp(save_options, 'time'))
                save_opts(1) = 1;
                time_ind = 1:size(app.data.mxy,2);
            else
                save_opts(1) = 0;
                time_ind = size(app.data.mxy,2)
            end
            if sum(strcmp(save_options, 'freqs'))
                save_opts(2) = 1;
                freq_ind = 1:numel(app.parameters.freq_hz_range);
            else
                save_opts(2) = 0;
                freq_ind = round(numel(app.parameters.freq_hz_range)/2);
            end
            if sum(strcmp(save_options, 'dist'))
                save_opts(3) = 1;
                dist_ind = 1:app.parameters.yrange_cm_npoints;
            else
                save_opts(3) = 0;
                dist_ind = round(app.parameters.yrange_cm_npoints/2);
            end
            if sum(strcmp(save_options, 'reps'))
                save_opts(4) = 1;
            end
            if sum(strcmp(save_options, 'plots'))
                save_opts(5) = 1;
            end

            try
                data.mxy = app.data.mxy(dist_ind, time_ind, freq_ind);
                data.mz = app.data.mz(dist_ind, time_ind, freq_ind);
                data.mxy_b1 = app.data.mxy_b1(dist_ind, time_ind, freq_ind);
                data.mz_b1 = app.data.mz_b1(dist_ind, time_ind, freq_ind);
            catch
                warning('could not save single excitation simulation');
            end

            % save repetead excitation simulation
            if save_opts(4)
                try
                data.mxy_rep = app.data.mxy_b1_rep;
                data.mz_rep= app.data.mz_b1_rep;
                catch
                    warning('could not save repeated excitation simulation!');
                end
            end

            % save repetead excitation simulation
            if ~save_opts(5)
                try
                data.mxy = [];
                data.mz= [];
                data.mxy_b1 = [];
                data.mz_b1= [];
                catch
                    warning('could not NOT save single excitation simulation!');
                end
            end

            parameters = app.parameters;

            if strcmp(app.save_fileEditField.Value, '')
                if app.simrangeCheckBox.Value
                    disp('');
                    pathname = app.parameters.save_multi_sim_path;
                    filename = ['\' num2str(app.parameters.multi_sim_n)];
                else
                    filterspec = {'*.mat'};
                    
                    try
                        pulse_path = app.parameters.rfpulse_file_path;
                        [filename, pathname] = uiputfile(filterspec, 'Choose file', fullfile(pulse_path, '\'));
                    catch
                        [filename, pathname] = uiputfile(filterspec, 'Choose file', fullfile(pulse_path, '\'));
                    end
                end
            

            else
%                 absPath = fileparts(mfilename('fullpath'))
                absPath = pwd();
%                 pathname = [absPath  '\simulations\'];
                pathname = [absPath '\'];
                filename = app.save_fileEditField.Value;
                if numel(filename) > 3
                    if ~strcmp(filename(end-3:end), 'mat')
                        filename = [filename '.mat'];
                    end
                end

            end
            % try saving data:
            try
                save([pathname filename], 'data', 'parameters');
            catch
                try
                    save([pathname filename '.mat'], 'data', 'parameters');
                catch
                    warning('data could not be saved');
                end
            end

                success = 1;
            catch
                success = 0;
        end
        end
        

        


        function  progress_simRep(app, ir, N)
            if round(ir/N, 2) < 0.1
                app.progress_simRepLabel_11.BackgroundColor = 'w';
                app.progress_simRepLabel_12.BackgroundColor = 'w';
                app.progress_simRepLabel_13.BackgroundColor = 'w';
                app.progress_simRepLabel_14.BackgroundColor = 'w';
                app.progress_simRepLabel_15.BackgroundColor = 'w';
                app.progress_simRepLabel_16.BackgroundColor = 'w';
                app.progress_simRepLabel_17.BackgroundColor = 'w';
                app.progress_simRepLabel_18.BackgroundColor = 'w';
                app.progress_simRepLabel_19.BackgroundColor = 'w';
                app.progress_simRepLabel_20.BackgroundColor = 'w';
            elseif (round(ir/N, 2) >= 0.1) && (round(ir/N, 2) < 0.2)
                app.progress_simRepLabel_11.BackgroundColor = 'g';
                app.progress_simRepLabel_12.BackgroundColor = 'w';
                app.progress_simRepLabel_13.BackgroundColor = 'w';
                app.progress_simRepLabel_14.BackgroundColor = 'w';
                app.progress_simRepLabel_15.BackgroundColor = 'w';
                app.progress_simRepLabel_16.BackgroundColor = 'w';
                app.progress_simRepLabel_17.BackgroundColor = 'w';
                app.progress_simRepLabel_18.BackgroundColor = 'w';
                app.progress_simRepLabel_19.BackgroundColor = 'w';
                app.progress_simRepLabel_20.BackgroundColor = 'w';
            elseif (round(ir/N, 2) >= 0.2) && (round(ir/N, 2) < 0.3)
                app.progress_simRepLabel_11.BackgroundColor = 'g';
                app.progress_simRepLabel_12.BackgroundColor = 'g';
                app.progress_simRepLabel_13.BackgroundColor = 'w';
                app.progress_simRepLabel_14.BackgroundColor = 'w';
                app.progress_simRepLabel_15.BackgroundColor = 'w';
                app.progress_simRepLabel_16.BackgroundColor = 'w';
                app.progress_simRepLabel_17.BackgroundColor = 'w';
                app.progress_simRepLabel_18.BackgroundColor = 'w';
                app.progress_simRepLabel_19.BackgroundColor = 'w';
                app.progress_simRepLabel_20.BackgroundColor = 'w';
            elseif (round(ir/N, 2) >= 0.3) && (round(ir/N, 2) < 0.4)
                app.progress_simRepLabel_11.BackgroundColor = 'g';
                app.progress_simRepLabel_12.BackgroundColor = 'g';
                app.progress_simRepLabel_13.BackgroundColor = 'g';
                app.progress_simRepLabel_14.BackgroundColor = 'w';
                app.progress_simRepLabel_15.BackgroundColor = 'w';
                app.progress_simRepLabel_16.BackgroundColor = 'w';
                app.progress_simRepLabel_17.BackgroundColor = 'w';
                app.progress_simRepLabel_18.BackgroundColor = 'w';
                app.progress_simRepLabel_19.BackgroundColor = 'w';
                app.progress_simRepLabel_20.BackgroundColor = 'w';
            elseif (round(ir/N, 2) >= 0.4) && (round(ir/N, 2) < 0.5)
                app.progress_simRepLabel_11.BackgroundColor = 'g';
                app.progress_simRepLabel_12.BackgroundColor = 'g';
                app.progress_simRepLabel_13.BackgroundColor = 'g';
                app.progress_simRepLabel_14.BackgroundColor = 'g';
                app.progress_simRepLabel_15.BackgroundColor = 'w';
                app.progress_simRepLabel_16.BackgroundColor = 'w';
                app.progress_simRepLabel_17.BackgroundColor = 'w';
                app.progress_simRepLabel_18.BackgroundColor = 'w';
                app.progress_simRepLabel_19.BackgroundColor = 'w';
                app.progress_simRepLabel_20.BackgroundColor = 'w';
            elseif (round(ir/N, 2) >= 0.5) && (round(ir/N, 2) < 0.6)
                app.progress_simRepLabel_11.BackgroundColor = 'g';
                app.progress_simRepLabel_12.BackgroundColor = 'g';
                app.progress_simRepLabel_13.BackgroundColor = 'g';
                app.progress_simRepLabel_14.BackgroundColor = 'g';
                app.progress_simRepLabel_15.BackgroundColor = 'g';
                app.progress_simRepLabel_16.BackgroundColor = 'w';
                app.progress_simRepLabel_17.BackgroundColor = 'w';
                app.progress_simRepLabel_18.BackgroundColor = 'w';
                app.progress_simRepLabel_19.BackgroundColor = 'w';
                app.progress_simRepLabel_20.BackgroundColor = 'w';
            elseif (round(ir/N, 2) >= 0.6) && (round(ir/N, 2) < 0.7)
                app.progress_simRepLabel_11.BackgroundColor = 'g';
                app.progress_simRepLabel_12.BackgroundColor = 'g';
                app.progress_simRepLabel_13.BackgroundColor = 'g';
                app.progress_simRepLabel_14.BackgroundColor = 'g';
                app.progress_simRepLabel_15.BackgroundColor = 'g';
                app.progress_simRepLabel_16.BackgroundColor = 'g';
                app.progress_simRepLabel_17.BackgroundColor = 'w';
                app.progress_simRepLabel_18.BackgroundColor = 'w';
                app.progress_simRepLabel_19.BackgroundColor = 'w';
                app.progress_simRepLabel_20.BackgroundColor = 'w';
            elseif (round(ir/N, 2) >= 0.7) && (round(ir/N, 2) < 0.8)
                app.progress_simRepLabel_11.BackgroundColor = 'g';
                app.progress_simRepLabel_12.BackgroundColor = 'g';
                app.progress_simRepLabel_13.BackgroundColor = 'g';
                app.progress_simRepLabel_14.BackgroundColor = 'g';
                app.progress_simRepLabel_15.BackgroundColor = 'g';
                app.progress_simRepLabel_16.BackgroundColor = 'g';
                app.progress_simRepLabel_17.BackgroundColor = 'g';
                app.progress_simRepLabel_18.BackgroundColor = 'w';
                app.progress_simRepLabel_19.BackgroundColor = 'w';
                app.progress_simRepLabel_20.BackgroundColor = 'w';
            elseif (round(ir/N, 2) >= 0.8) && (round(ir/N, 2) < 0.9)
                app.progress_simRepLabel_11.BackgroundColor = 'g';
                app.progress_simRepLabel_12.BackgroundColor = 'g';
                app.progress_simRepLabel_13.BackgroundColor = 'g';
                app.progress_simRepLabel_14.BackgroundColor = 'g';
                app.progress_simRepLabel_15.BackgroundColor = 'g';
                app.progress_simRepLabel_16.BackgroundColor = 'g';
                app.progress_simRepLabel_17.BackgroundColor = 'g';
                app.progress_simRepLabel_18.BackgroundColor = 'g';
                app.progress_simRepLabel_19.BackgroundColor = 'w';
                app.progress_simRepLabel_20.BackgroundColor = 'w';
            elseif (round(ir/N, 2) >= 0.9) && (round(ir/N, 2) < 1)
                app.progress_simRepLabel_11.BackgroundColor = 'g';
                app.progress_simRepLabel_12.BackgroundColor = 'g';
                app.progress_simRepLabel_13.BackgroundColor = 'g';
                app.progress_simRepLabel_14.BackgroundColor = 'g';
                app.progress_simRepLabel_15.BackgroundColor = 'g';
                app.progress_simRepLabel_16.BackgroundColor = 'g';
                app.progress_simRepLabel_17.BackgroundColor = 'g';
                app.progress_simRepLabel_18.BackgroundColor = 'g';
                app.progress_simRepLabel_19.BackgroundColor = 'g';
                app.progress_simRepLabel_20.BackgroundColor = 'w';
            elseif round(ir/N, 2) == 1
                app.progress_simRepLabel_11.BackgroundColor = 'g';
                app.progress_simRepLabel_12.BackgroundColor = 'g';
                app.progress_simRepLabel_13.BackgroundColor = 'g';
                app.progress_simRepLabel_14.BackgroundColor = 'g';
                app.progress_simRepLabel_15.BackgroundColor = 'g';
                app.progress_simRepLabel_16.BackgroundColor = 'g';
                app.progress_simRepLabel_17.BackgroundColor = 'g';
                app.progress_simRepLabel_18.BackgroundColor = 'g';
                app.progress_simRepLabel_19.BackgroundColor = 'g';
                app.progress_simRepLabel_20.BackgroundColor = 'g';
            end
            pause(0.00000001);
        end
        
        function results = calc_fa_b1_app(app, fa_changed, b1_changed)
            if fa_changed
                gamma = app.parameters.gamma;
                khz_in_g = gamma / 1000;
                trf = app.parameters.pulse_duration_ms * 1e-3;    
                alpha = app.flipangle_degEditField.Value; % Degrees.
                sint = app.parameters.sint;
                b1_amp_g = (pi * alpha / 180 ) * 1/((2*pi*gamma*1e4) * trf* sint) * 1e4;% G
                % b1_amp_G = pi /2 * 1/((2*pi*gamma*1e4) * Trf* sint) * 1e4 * 4  % G 
                app.parameters.b1_amp_khz = b1_amp_g * khz_in_g;
                app.b1_mag_khz_EditField.Value = app.parameters.b1_amp_khz;
            else
                gamma = app.parameters.gamma;
                khz_in_g = gamma / 1000;
                trf = app.parameters.pulse_duration_ms * 1e-3;
%                 alpha = app.flipangle_degEditField.Value; % Degrees.
                sint = app.parameters.sint;
                app.parameters.b1_amp_khz;
                % b1_amp_G = pi /2 * 1/((2*pi*gamma*1e4) * Trf* sint) * 1e4 * 4  % G 
                b1_amp_g  = app.parameters.b1_amp_khz / khz_in_g ;
                alpha = 180* ((2*pi*gamma*1e4) * trf* sint)* b1_amp_g / (1e4) / pi;% G
                app.parameters.alpha = alpha;
                app.flipangle_degEditField.Value = alpha;
                pause(0.00001);
            end
        end

        function results = calc_b1(app, fa)
            gamma = app.parameters.gamma;
            khz_in_g = gamma / 1000;
            trf = app.parameters.pulse_duration_ms * 1e-3;    
            alpha = fa; % Degrees.
            sint = app.parameters.sint;
            
            b1_amp_g = (pi * alpha / 180 ) * 1/((2*pi*gamma*1e4) * trf* sint) * 1e4;% G
            
            results.alpha = alpha;
            results.b1_amp_khz = b1_amp_g * khz_in_g;
            results.b1_amp_g = b1_amp_g;
        end

        function results = calc_fa(app, b1_amp_khz)
            gamma = app.parameters.gamma;
            khz_in_g = gamma / 1000;
            trf = app.parameters.pulse_duration_ms * 1e-3;
%                 alpha = app.flipangle_degEditField.Value; % Degrees.
            sint = app.parameters.sint;
            % b1_amp_G = pi /2 * 1/((2*pi*gamma*1e4) * Trf* sint) * 1e4 * 4  % G 
            b1_amp_g  = b1_amp_khz / khz_in_g ;
            alpha = 180* ((2*pi*gamma*1e4) * trf* sint)* b1_amp_g / (1e4) / pi;% G

            results.alpha = alpha;
            results.b1_amp_khz = b1_amp_khz;
            results.b1_amp_g = b1_amp_g;
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % change folder to file directory
            if(~isdeployed)
                cd(fileparts(which(mfilename)))
            end
            % set parameters values to start with
            app.parameters.empty = [];
            app.parameters.gradient_kHz_cm = app.gradientstrength_khzcmEditField.Value;
            app.parameters.gamma = 1070.84; % [1H]
            app.parameters.T1 = str2double(app.T1_sEditField.Value);
            app.parameters.T2 = str2double(app.T2_sEditField.Value);
            app.parameters.sample_mag = app.HPEditField.Value;
            
            % y range:
            app.parameters.yrange_cm = app.yrange_cmEditField.Value;
            app.parameters.yrange_cm_npoints = str2double(app.yrange_cm_npointsEditField.Value);
            
            % pulse duration
            app.parameters.pulse_duration_ms = app.pulseduration_msEditField.Value;
            app.parameters.pulse_duration_ms_npoints = str2double(app.pulseduration_ms_npointsEditField.Value);
            app.parameters.pulse_duration_ms_res = app.parameters.pulse_duration_ms / app.parameters.pulse_duration_ms_npoints;
            app.duration_refocus_grad.Value = app.parameters.pulse_duration_ms_res;
            
            % frequency [hz]
            app.parameters.freq_hz_high = app.freq_hz_highEditField.Value;
            app.parameters.freq_hz_low = app.freq_hz_lowEditField.Value;
            app.parameters.freq_hz_npoints = str2double(app.freq_hz_npointsEditField.Value);
            app.parameters.freq_ind = 1;

            % frequency [ppm]
            app.parameters.freq_ppm1 = app.Freq1_ppmEditField.Value;
            app.parameters.freq_ppm2 = app.Freq2_ppmEditField.Value;
            app.parameters.freq_ppm_basic = app.BasicFreqEditField.Value;
            app.parameters.freq_hz_pulse = app.PulseFreqEditField.Value;

            % flipangle
            app.parameters.alpha = app.flipangle_degEditField.Value;
            app.parameters.b1_amp_khz = app.b1_mag_khz_EditField.Value;

            % integration factor
            app.parameters.sint = app.sint_EditField.Value;
            app.parameters.bwfac_hz_ms = 1000/app.parameters.sint;

            % run simulation for finding gradient strength for B1
            % compensation:
            app.parameters.find_grad_low = app.find_Gradients_lowEditField.Value;
            app.parameters.find_grad_high = app.find_Gradients_highEditField.Value;
            app.parameters.find_grad_npoints  = app.find_Gradients_npointsEditField.Value;

            app.parameters.find_grad_ind1 = app.find_Gradients_ind1EditField.Value;
            app.parameters.find_grad_ind2 = app.find_Gradients_ind2EditField.Value;

            % run simulation of repeated excitation (to see the Steady
            % state)
            app.parameters.repetitions_time_ms = app.repetition_time_msEditField.Value;
            app.parameters.repetitions_time_ms_npoints = app.repetition_time_npointsEditField.Value;
            app.parameters.nrepetitions  = app.nrepetitions_EditField.Value;
            app.parameters.alternating_phase  = app.incphaseEditField.Value;
            app.parameters.first_pulse_amp = app.amp_1stpulse_EditField.Value;
            app.parameters.first_pulse_always = app.alpha2prepCheckBox.Value;
            app.parameters.save_all_sim_data = app.highresoutCheckBox.Value;
            app.parameters.ref_grad_dur_ms = app.repetition_time_msEditField.Value - 2*app.pulseduration_msEditField.Value;
            app.parameters.first_repetition_time_ms  = app.parameters.repetitions_time_ms;
            app.first_repetition_time_msEditField.Value = app.parameters.first_repetition_time_ms;

            % 2nd pulse:
            app.parameters.simulate_alternating_freqs = 'off';
            app.parameters.alpha_2ndfreq = app.alpha2ndPulseEditField.Value;
            app.parameters.freq_2ndfreq_hz = app.freq_2ndpulse_HzEditField.Value;
            app.parameters.alternate_freq_Ntimes_2nd_freq = app.AlterNTimesEditField.Value;
            app.parameters.phase_pulse_deg = app.altphaseEditField.Value;

            % 1st pulse 
            app.parameters.freq_1stfreq_hz = app.freq_1stpulse_HzEditField.Value;

            % Spoiler Gradient:
            app.parameters.spoiler_grad = app.SpoilerGradCheckBox.Value;

            % Tipback pulse (after N repetitions)
            app.parameters.tipback_pulse= app.TipbackCheckBox.Value;
            app.parameters.phase_tipback_pulse = app.phasetbEditField.Value;
            
            % number of rounds to simulated (N phase encodes + tipback,
            % spoiling, perparation pulse)
            app.parameters.Nsims = app.parameters.nrepetitions;


            % use map instead of profile:
            app.parameters.use_map = app.usemapCheckBox.Value;
            % pointer position 
            app.parameters.pointer2d = [1 1];
            % wether rf power should be found:
            app.parameters.find_rfpower = app.findRFPowerCheckBox.Value;
            % init refpowmap FOV
            app.parameters.refpow_map_fov_mm = [0 0];
            app.parameters.refpow_pos_x_mm = app.xmmEditField.Value;
            app.parameters.refpow_pos_y_mm = app.ymmEditField.Value;
            app.parameters.mat_refpow_map = [0 0];

            % filenames:
            app.parameters.rfpulse_file_name = '';
            app.parameters.rfpulse_file_path = '';

            % wether to simulate a range of values:
            app.parameters.save_multi_sim_path = '';
            app.parameters.simulate_range_of_vals = app.simrangeCheckBox.Value;
            app.parameters.multi_sim_n =0;
            app.parameters.fa_range_val = 0;
            app.parameters.tr_range_val = 0;
            app.parameters.t1_range_val = 0;
            app.parameters.t2_range_val = 0;

            % anatomical image parameters:
            app.parameters.anat_image_fov_mm = [0 0];
            app.parameters.anat_image_mat= [0 0];
            app.parameters.anat_image_off_mm = [0 0];

            % plotting parameters:
            app.parameters.plot_bssfp_signal = app.bssfpsignalCheckBox.Value;

            app.parameters.plot_homo_mx = 0;
            app.parameters.plot_homo_my = 0;
            app.parameters.plot_homo_mxy = 1;
            app.parameters.plot_homo_mz = 0;
            app.parameters.plot_homo_phase = 0;
            app.parameters.plot_homo_none = 0;


            app.parameters.plot_inhomo_mx = 0;
            app.parameters.plot_inhomo_my = 0;
            app.parameters.plot_inhomo_mxy = 1;
            app.parameters.plot_inhomo_mz = 0;
            app.parameters.plot_inhomo_phase = 0;
            app.parameters.plot_inhomo_none = 0;

            app.parameters.anat_image_slice = 1;
            app.parameters.anat_nslices = 1;

            % data
            app.data.rfpulse = zeros(3,1);
            app.data.b1TransmitProfileFit_log = ones(64, 1);

            app.data.mxy_b1 = zeros(3,3);
            app.data.mz_b1 = zeros(3,3);

            app.data.mxy_b1 = zeros(3,3);

            app.data.mxy = zeros(3,3);
            app.data.mz = zeros(3,3);

            app.data.anat_image = zeros(3,3);

            % 
            app.data.refpow_map = zeros(3,3);
            app.data.refpow_map_calc = zeros(3,3);

            if app.UseFrequenciesCheckBox.Value
                app.Freq1_ppmEditField.Enable = 'off';
                app.Freq2_ppmEditField.Enable = 'off';
                app.BasicFreqEditField.Enable = 'off';

            else
                app.Freq1EditField.Enable = 'off';
                app.Freq2EditField.Enable = 'off';
            end


            % calc resolutions:
            calc_freq_hz_res(app);
            calc_yrange_cm_res(app);
            calc_pulse_duration_ms_res(app);
            calc_slice_thick(app);
            calc_fa_b1_app(app,1,0);

            set(0, 'DefaultLineLineWidth', 1);
            
            %app.parameters.pulse_duration_ms = app.pulseduration_msEditField.Value;
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.calc_rf_pulse_toolUIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {921, 921};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {445, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end

        % Button pushed function: SelectRFPulse_Button
        function SelectRFPulse_ButtonPushed(app, event)
            try 
               [filename, pathname] = uigetfile('*', 'Select Pulse ', 'D:\LRZ Sync+Share\phd\bSSFP\uniform excitation surface');
            catch
               [filename, pathname] = uigetfile('');
            end
            try
                app.parameters.rfpulse_file = [pathname filename];
                app.RFPulseEditField.Value = filename;
                app.parameters.rfpulse_file_name = filename;
                app.parameters.rfpulse_file_path = pathname;
            end
            figure(app.calc_rf_pulse_toolUIFigure);
        end

        % Button pushed function: LoadRFPulse_Button
        function LoadRFPulse_ButtonPushed(app, event)
            % load rf pulse:
            load_rfpulse(app);
            
            % calculate integration factor
            calc_sint(app);
            % calculate slice thickness:
            calc_slice_thick(app);
            
            % plot the pulse
            figure(app.calc_rf_pulse_toolUIFigure);
            plot_rfpulse(app);
        end

        % Button pushed function: SelectB1Profile_Button
        function SelectB1Profile_ButtonPushed(app, event)
            try 
               [filename, pathname] = uigetfile('*.mat', 'Select B1 Profile ', 'D:\LRZ Sync+Share\phd\bSSFP\uniform excitation surface');
            catch
               [filename, pathname] = uigetfile('');
            end
            try
                app.parameters.b1profile_file = [pathname filename];
                app.B1ProfileEditField.Value = filename;
            end
            figure(app.calc_rf_pulse_toolUIFigure);
        end

        % Button pushed function: LoadB1Profile_Button
        function LoadB1Profile_ButtonPushed(app, event)
            % try to load (assume cryo B1)
            try
                load(app.parameters.b1profile_file)
            catch
                warning("can't load B1 profile");
            end

            try
                if (ndims(squeeze(fit_res.amp)) == 2)
                    b1TransmitProfile = fit_res.amp(:, end/2);
                    % fit along high SNR points
                    b1TransmitProfileFit_log = fit((6:32).', log(abs(b1TransmitProfile(6:32))), 'poly1');
                    % use fit curve to extent to 
                    app.data.b1TransmitProfileFit_log = b1TransmitProfileFit_log(1:64);
                    % save reference power map:
                    app.data.refpow_map = squeeze(fit_res.refpow_map(:, :));
                    %% for now
                    app.data.refpow_map_calc = app.data.refpow_map;
                    %%
                    % save fov:
                    app.parameters.refpow_map_fov_mm = fit_res.header.p.FOV;
                    app.parameters.refpow_map_header = fit_res.header;
                    app.parameters.refpow_map_off_mm = [fit_res.header.method.PVM_Phase1Offset ...
                                                        fit_res.header.method.PVM_Phase0Offset];
                    app.parameters.mat_refpow_map = size(app.data.refpow_map);
                    % in case it was red:
                    app.LoadB1Profile_Button.BackgroundColor = [0.96,0.96,0.96];                    % plot
                    plot_b1profile(app);
                else
                    b1TransmitProfile = fit_res.amp(:, end/2, end/2);
                    % fit along high SNR points
                    b1TransmitProfileFit_log = fit((6:32).', log(abs(b1TransmitProfile(6:32))), 'poly1');
                    % use fit curve to extent to 
                    app.data.b1TransmitProfileFit_log = b1TransmitProfileFit_log(1:64);
                    % save reference power map:
                    app.data.refpow_map = squeeze(fit_res.refpow_map(:, :, end/2));
                    %% for now
                    app.data.refpow_map_calc = app.data.refpow_map;
                    %%
                    app.parameters.refpow_map_fov_mm = fit_res.header.p.FOV;
                    app.parameters.refpow_map_header = fit_res.header;
                    app.parameters.refpow_map_off_mm = [fit_res.header.method.PVM_Phase1Offset ...
                                                        fit_res.header.method.PVM_Phase0Offset];
                    app.parameters.mat_refpow_map = size(app.data.refpow_map);
                    % in case it was red:
                    app.LoadB1Profile_Button.BackgroundColor = [0.96,0.96,0.96];
                    % plot
                    plot_b1profile(app);
                end
            catch
                warning('could not load B1 Profile!')
                app.LoadB1Profile_Button.BackgroundColor = 'r';
            end
            try
                % try to interpolate the reference power map onto the
                % anatomical image:
                calc_refpow_anat(app);
            end

            figure(app.calc_rf_pulse_toolUIFigure);

            

        end

        % Selection changed function: NucleusButtonGroup
        function NucleusButtonGroupSelectionChanged(app, event)
            selectedButton = app.NucleusButtonGroup.SelectedObject;
            if sum(app.CButton.Value) % default, 
                app.parameters.gamma = 1070.84; % Hz/G 13C
            else
                app.parameters.gamma = 4257.7478518;		% Hz/G. 1H
            end
            figure(app.calc_rf_pulse_toolUIFigure);
        end

        % Value changed function: yrange_cmEditField
        function yrange_cmEditFieldValueChanged(app, event)
            value = app.yrange_cmEditField.Value;
            app.parameters.yrange_cm = value;
            calc_yrange_cm_res(app);
        end

        % Button pushed function: SimulateButton
        function SimulateButtonPushed(app, event)
            % get pulse:
            try
                pulse = app.data.rfpulse;
            end

            % get B1:
            try
                b1profile_log = app.data.b1TransmitProfileFit_log;
            end
            if (sum(b1profile_log) == 0)
                b1profile_log = ones(64, 1);
            end

            calc_fa_b1_app(app,1,0);

            % perform Bloch Simulation:
            % 1.st get parameters:
            Trf = app.parameters.pulse_duration_ms * 1e-3;    
            dt = app.parameters.pulse_duration_ms_res * 1e-3;
            alpha = app.flipangle_degEditField.Value; % Degrees.
            
            % gamma = 4257.7478518;		% Hz/G. 1H
            gamma = app.parameters.gamma; % Hz/G 13C
            
            t1 = app.parameters.T1;			% Sec.
            t2 = app.parameters.T2;		% Sec.
            freq = linspace(app.parameters.freq_hz_low, ...
                            app.parameters.freq_hz_high, ...
                            app.parameters.freq_hz_npoints); % Hz
%            N = 100; % 5 x T1 for steady state without alpha/2
 %           Ttot = 10 * Trf;
            sint = app.parameters.sint;
            bwfac =  1000/sint;
            % khz_in_g = 4.25755; % 4.25755 KHz =^=  1 G
            khz_in_g = gamma / 1000;

            % time resolution simulation
            dt_sim = app.parameters.pulse_duration_ms_res*1e-3;
            % time course simulation
            t_rf_sim = 0:dt_sim:Trf-dt_sim;

            % time resolution loaded pulse
            dt_b1 = Trf / numel(pulse(:,1));
            % time course loaded pulse
            t_b1 = 0:dt_b1:Trf-dt_b1;
            
            % interpolate loaded pulse onto to simulate pulse:
            pulse_sim = interp1(t_b1, pulse.', t_rf_sim);
            pulse_sim(isnan(pulse_sim)) = 0;

            % calculate pulse amplitude (if B1 amp in settings = 0):
%             if (app.parameters.b1_amp_khz == 0)
             % b1_amp_G = pi /2 * 1/((2*pi*gamma*1e4) * Trf* sint) * 1e4 * 4  % G 
            
            b1_amp_g  = app.parameters.b1_amp_khz / khz_in_g;
%             else
%                 b1_amp_khz  = app.parameters.b1_amp_khz;
%                 b1_amp_g =  b1_amp_khz / khz_in_g;
%                 alpha = b1_amp_g / pi * 180 *((2*pi*gamma*1e4) * Trf* sint) / 1e4;
%                 app.flipangle_degEditField.Value = alpha;
%             end

            % get pulse in Gauss

            pulse_sim_G = pulse_sim * b1_amp_g;
%             pulse_sim_G = flip(pulse_sim_G);

            % calculate gradients:
            [grad_kHz_cm_range_sim, gradient_G_cm_range_sim] = calc_gradients(app, dt_sim, Trf);
            
            % B1 measured range:
            x_b1_cm = linspace(0.25, 31.75, 64)/10;
            x_b1_cm = x_b1_cm - 1.6;

            % B1 simulated range:
            x_sim_cm = linspace(-app.parameters.yrange_cm/2, ...
                                 app.parameters.yrange_cm/2, ...
                                 app.parameters.yrange_cm_npoints).';

            % interpolate
            b1Transmit_sim_log = interp1(x_b1_cm, b1profile_log, x_sim_cm, 'spline', 'extrap');
            b1Transmit_sim = exp(b1Transmit_sim_log) / max(exp(b1Transmit_sim_log ));

            % add 0s to the pulse:
            pulse_sim_G = [pulse_sim_G zeros(1, numel(pulse_sim_G))];

            % Bloch Simulation:
            mxy_x_b1 = zeros(numel(x_sim_cm), numel(pulse_sim_G), numel(freq));
            mz_x_b1 = zeros(numel(x_sim_cm), numel(pulse_sim_G), numel(freq));
            
            app.SimulateButton.Text = 'Simulating!';
            app.SimulateButton.FontAngle = 'italic';
            app.SimulateButton.BackgroundColor = 'r';

            pause(0.00000001);
            
            % inhomogeneous B1
            if sum(app.CButton.Value) % C13
                for i = 1:numel(x_sim_cm)
                    [mxx, myx, mzx] = bloch_13c_no_output(pulse_sim_G * b1Transmit_sim(i), gradient_G_cm_range_sim, dt_sim, t1, t2, freq, x_sim_cm(i), 2);
                    mxy_x_b1(i, :, :) = mxx + 1i* myx;
                    mz_x_b1(i, :, :) = mzx;
                end
            else
                for i = 1:numel(x_sim_cm)
                    [mxx, myx, mzx] = bloch_no_output(pulse_sim_G * b1Transmit_sim(i), gradient_G_cm_range_sim, dt_sim, t1, t2, freq, x_sim_cm(i), 2);
                    mxy_x_b1(i, :, :) = mxx + 1i* myx;
                    mz_x_b1(i, :, :) = mzx;
                end
            end
            % homogeneous B1
            mxy_x = zeros(numel(x_sim_cm), numel(pulse_sim_G), numel(freq));
            mz_x = zeros(numel(x_sim_cm), numel(pulse_sim_G), numel(freq));

            if sum(app.CButton.Value) % C13
                [mxx, myx, mzx] = bloch_13c_no_output(pulse_sim_G, gradient_G_cm_range_sim, dt_sim, t1, t2, freq, x_sim_cm, 2);
                mxy_x = mxx + 1i* myx;
                mz_x = mzx;
            else
                [mxx, myx, mzx] = bloch_no_output(pulse_sim_G, gradient_G_cm_range_sim, dt_sim, t1, t2, freq, x_sim_cm, 2);
                mxy_x = mxx + 1i* myx;
                mz_x = mzx;
            end

            app.SimulateButton.Text = 'Simulate!';
            app.SimulateButton.FontAngle = 'normal';
            app.SimulateButton.BackgroundColor = [0.96,0.96,0.96];
            
            % save simulations:
            app.data.mxy_b1 = mxy_x_b1;
            app.data.mz_b1 = mz_x_b1;

            % permute to match b1 mxy + b1 mz
            app.data.mxy = permute(mxy_x, [2 1 3]);  
            app.data.mz = permute(mz_x, [2 1 3]);  

            % plot simulated data:
            plot_b1profile(app);
            plot_rfpulse(app);
            try
                plot_mag_vs_dist(app, app.parameters.freq_ind);
            catch
                plot_mag_vs_dist(app, 1);
            end

        end

        % Value changed function: gradientstrength_khzcmEditField
        function gradientstrength_khzcmEditFieldValueChanged(app, event)
            value = app.gradientstrength_khzcmEditField.Value;
            app.parameters.gradient_kHz_cm = value;
            calc_slice_thick(app);
        end

        % Value changed function: T1_sEditField
        function T1_sEditFieldValueChanged(app, event)
            value = str2double(app.T1_sEditField.Value);
            app.parameters.T1 = value;
        end

        % Value changed function: T2_sEditField
        function T2_sEditFieldValueChanged(app, event)
            value = str2double(app.T2_sEditField.Value);
            app.parameters.T2 = value;
        end

        % Value changed function: pulseduration_msEditField
        function pulseduration_msEditFieldValueChanged(app, event)
            value = app.pulseduration_msEditField.Value;
            app.parameters.pulse_duration_ms = value;
            % treat pulse duration as change in flipangle:
            calc_fa_b1_app(app, 1, 0);
            calc_pulse_duration_ms_res(app);
            calc_slice_thick(app);
            plot_rfpulse(app);
            
        end

        % Value changed function: pulseduration_ms_npointsEditField
        function pulseduration_ms_npointsEditFieldValueChanged(app, event)
            value = str2double(app.pulseduration_ms_npointsEditField.Value);
            app.parameters.pulse_duration_ms_npoints = value;
            calc_pulse_duration_ms_res(app);
        end

        % Value changed function: yrange_cm_npointsEditField
        function yrange_cm_npointsEditFieldValueChanged(app, event)
            value = str2double(app.yrange_cm_npointsEditField.Value);
            app.parameters.yrange_cm_npoints = value;
            calc_yrange_cm_res(app);
        end

        % Value changed function: freq_hz_npointsEditField
        function freq_hz_npointsEditFieldValueChanged(app, event)
            value = str2double(app.freq_hz_npointsEditField.Value);
            app.parameters.freq_hz_npoints = value;
            calc_freq_hz_res(app);
            app.parameters.freq_ind = round(value/2);
            the_updater_func(app);
        end

        % Value changed function: freq_hz_highEditField
        function freq_hz_highEditFieldValueChanged(app, event)
            value = app.freq_hz_highEditField.Value;
            app.parameters.freq_hz_high = value;
            calc_freq_hz_res(app);
        end

        % Value changed function: freq_hz_lowEditField
        function freq_hz_lowEditFieldValueChanged(app, event)
            value = app.freq_hz_lowEditField.Value;
            app.parameters.freq_hz_low = value;
            calc_freq_hz_res(app);

        end

        % Button pushed function: calcButton
        function calcButtonPushed(app, event)
            calc_sint(app);
            app.bwfac_HzsEditField.Value = app.parameters.bwfac_hz_ms;
            app.sint_EditField.Value = app.parameters.sint;
%             sint = app.
        end

        % Value changed function: sint_EditField
        function sint_EditFieldValueChanged(app, event)
            value = app.sint_EditField.Value;
            app.parameters.sint = value;            
        end

        % Value changed function: flipangle_degEditField
        function flipangle_degEditFieldValueChanged(app, event)
            value = app.flipangle_degEditField.Value;
            calc_fa_b1_app(app, 1, 0);
            app.parameters.alpha = value;
            plot_rfpulse(app);
            
        end

        % Value changing function: freqSlider
        function freqSliderValueChanging(app, event)
            changingValue = event.Value;
            freq = linspace(app.parameters.freq_hz_low, ...
                            app.parameters.freq_hz_high, ...
                            app.parameters.freq_hz_npoints); % Hz
            [~,ind] = min(abs(freq-changingValue));
            app.freqSlider.Value=freq(ind);
        end

        % Value changed function: freqSlider
        function freqSliderValueChanged(app, event)
            value = app.freqSlider.Value;
            freq_range = linspace(app.parameters.freq_hz_low, ...
                            app.parameters.freq_hz_high, ...
                            app.parameters.freq_hz_npoints); % Hz
            [~,freq_ind] = min(abs(freq_range-value));
            % save index:
            app.parameters.freq_ind = freq_ind;
            app.freqSlider.Value = freq_range(freq_ind);
            freq = freq_range(freq_ind);
            % if option of plotting 2 frequencies is acitvated:
            if app.Show2FreqsCheckBox.Value
%                 freq_low  = value + app.Freq1EditField.Value;
%                 freq_high = value + app.Freq2EditField.Value;
%                 if ((freq_low < freq(1)) && ... % bad
%                     freq_high < freq(end)) % good
%                     % calc difference from limit to set value:
%                     app.Freq1EditField.Value = freq(1);
%                     app.Freq2EditField.Value = freq(1) + abs(freq_high-freq_low);
%                     
%                     % adjust 
%                 else
%                     app.Freq1EditField.Value = freq_low;
%                 end
%                 if (freq_high > freq(end))
%                     app.Freq2EditField.Value = freq(end);
%                 else
%                     app.Freq2EditField.Value = freq_high;
%                 end
            end



            % plot magnetization on this frequency:
            if app.parameters.plot_bssfp_signal
                plot_bssfp_vs_freq(app, freq);
            else
                plot_mag_vs_dist(app, freq );
            end
            try
                plot_anat_image(app, app.parameters.anat_image_slice, app.ShowExcitationProfileCheckBox.Value);
            end
            
        end

        % Value changed function: Freq2EditField
        function Freq2EditFieldValueChanged(app, event)
            value = app.Freq2EditField.Value;
            if app.parameters.plot_bssfp_signal
                plot_bssfp_vs_freq(app, value);
            else
                plot_mag_vs_dist(app, value);
            end
            try
                plot_anat_image(app, app.parameters.anat_image_slice, app.ShowExcitationProfileCheckBox.Value);
            end
        end

        % Value changed function: Freq1EditField
        function Freq1EditFieldValueChanged(app, event)
            value = app.Freq1EditField.Value;
            if app.parameters.plot_bssfp_signal
                plot_bssfp_vs_freq(app, value);
            else
                plot_mag_vs_dist(app, value);
            end
            try
                plot_anat_image(app, app.parameters.anat_image_slice, app.ShowExcitationProfileCheckBox.Value);
            end
        end

        % Button pushed function: SelectAnatImage_Button
        function SelectAnatImage_ButtonPushed(app, event)
            try 
               pathname = uigetdir('D:\LRZ Sync+Share\phd\bSSFP\uniform excitation surface');
            catch
               pathname = uigetdir('');
            end
            try
                app.parameters.anatimage_path = pathname;
                app.AnatImageEditField.Value = pathname;
            end
        end

        % Button pushed function: LoadAnatImage_Button
        function LoadAnatImage_ButtonPushed(app, event)
             try
                imageObj = ImageDataObject(app.parameters.anatimage_path);
                imageObj = imageObj.readReco;
                imageObj = imageObj.readMethod;
                imageObj = imageObj.readAcqp;
%                 anatimage = 
            catch
                warning('Anatomical Image could no be loaded');
            end
            anat_image = squeeze(imageObj.data);
            if (imageObj.Method.PVM_NSPacks == 1)
                anat_image = reshape(anat_image, [size(anat_image) 1]);
            end

            patient_pos = imageObj.Acqp.ACQ_patient_pos;
            if strcmp(patient_pos, 'Head_Supine')
                % this is the orientation as "material (PV6) and matches
                % the B1/Refpowmap
                anat_image = rot90(anat_image, -1);
                app.data.anat_image = anat_image;

                app.parameters.mat_anat_image = size(anat_image);
                app.parameters.anat_image_fov_mm = flip(imageObj.Method.PVM_Fov);
                app.parameters.anat_nslices = imageObj.Method.PVM_SPackArrNSlices;
                app.parameters.anat_image_off_mm = [imageObj.Method.PVM_Phase1Offset(1) ...
                    imageObj.Method.PVM_Phase0Offset(1) ...
                    imageObj.Method.PVM_Phase2Offset(1)];
            elseif strcmp(patient_pos, 'Head_Prone')
                 anat_image = rot90(anat_image, -1);
                 anat_image = flip(anat_image, 1);
                 app.data.anat_image = anat_image;

                app.parameters.mat_anat_image = size(anat_image);
                app.parameters.anat_image_fov_mm = flip(imageObj.Method.PVM_Fov);
                app.parameters.anat_nslices = imageObj.Method.PVM_SPackArrNSlices;
                app.parameters.anat_image_off_mm = [imageObj.Method.PVM_Phase1Offset(1) ...
                    imageObj.Method.PVM_Phase0Offset(1) ...
                    imageObj.Method.PVM_Phase2Offset(1)];
            elseif strcmp(patient_pos, 'Foot_Prone')
                anat_image = rot90(anat_image, -1);
                 anat_image = flip(anat_image, 1);
                 app.data.anat_image = anat_image;

                app.parameters.mat_anat_image = size(anat_image);
                app.parameters.anat_image_fov_mm = flip(imageObj.Method.PVM_Fov);
                app.parameters.anat_nslices = imageObj.Method.PVM_SPackArrNSlices;
                app.parameters.anat_image_off_mm = [imageObj.Method.PVM_Phase1Offset(1) ...
                    imageObj.Method.PVM_Phase0Offset(1) ...
                    imageObj.Method.PVM_Phase2Offset(1)];
            end
            
            
            % calculate refereference power map:
            try
                calc_refpow_anat(app);
            end

            % plot the anatomical image
            plot_anat_image(app, round(app.parameters.anat_nslices/2) , 1);
        end

        % Value changed function: ShowExcitationProfileCheckBox
        function ShowExcitationProfileCheckBoxValueChanged(app, event)
            value = app.ShowExcitationProfileCheckBox.Value;
            try
                plot_anat_image(app, app.parameters.anat_image_slice, value);
            end
        end

        % Value changed function: SliceSpinner
        function SliceSpinnerValueChanged(app, event)
            value = app.SliceSpinner.Value;
            if value > app.parameters.anat_nslices
                app.SliceSpinner.Value = app.parameters.anat_nslices;
                value = app.SliceSpinner.Value;
            end

            app.parameters.anat_image_slice = value;
            try
                plot_anat_image(app, app.parameters.anat_image_slice, app.ShowExcitationProfileCheckBox.Value);
            end
        end

        % Value changed function: Show2FreqsCheckBox
        function Show2FreqsCheckBoxValueChanged(app, event)
            value = app.Show2FreqsCheckBox.Value;
            app.parameters.freq_diff = app.Freq2EditField.Value  - app.Freq1EditField.Value;
            if value
                app.freqSlider.Enable = "off";
            else
                app.freqSlider.Enable = "on";
            end
            freqSliderValueChanged(app);
        end

        % Value changed function: UseFrequenciesCheckBox
        function UseFrequenciesCheckBoxValueChanged(app, event)
            value = app.UseFrequenciesCheckBox.Value;
            if value % use frequencies
                app.Freq1_ppmEditField.Enable = 'off';
                app.Freq2_ppmEditField.Enable = 'off';
                app.BasicFreqEditField.Enable = 'off';

                app.Freq1EditField.Enable = 'on';
                app.Freq2EditField.Enable = 'on';
                plot_mag_vs_dist(app, app.parameters.freq_ind);
                plot_anat_image(app, app.parameters.anat_image_slice, app.ShowExcitationProfileCheckBox.Value);
            else % use ppm
                app.Freq1EditField.Enable = 'off';
                app.Freq2EditField.Enable = 'off';

                app.Freq1_ppmEditField.Enable = 'on';
                app.Freq2_ppmEditField.Enable = 'on';
                app.BasicFreqEditField.Enable = 'on';

                app.parameters.freq_ppm2
                app.parameters.freq_ppm2
                
                
                plot_mag_vs_dist(app, app.parameters.freq_ind);
                plot_anat_image(app, app.parameters.anat_image_slice, app.ShowExcitationProfileCheckBox.Value);
            end
            
        end

        % Value changed function: PulseFreqEditField
        function PulseFreqEditFieldValueChanged(app, event)
            value = app.PulseFreqEditField.Value;
            app.parameters.freq_hz_pulse = value;
            plot_mag_vs_dist(app, app.parameters.freq_hz_pulse);
            plot_anat_image(app, app.parameters.anat_image_slice, app.ShowExcitationProfileCheckBox.Value);
        end

        % Value changed function: BasicFreqEditField
        function BasicFreqEditFieldValueChanged(app, event)
            value = app.BasicFreqEditField.Value;
            app.parameters.freq_ppm_basic = value;
            plot_mag_vs_dist(app, app.parameters.freq_ind);
            plot_anat_image(app, app.parameters.anat_image_slice, app.ShowExcitationProfileCheckBox.Value);
        end

        % Value changed function: Freq2_ppmEditField
        function Freq2_ppmEditFieldValueChanged(app, event)
            value = app.Freq2_ppmEditField.Value;
            app.parameters.freq_ppm2 = value;
            plot_mag_vs_dist(app, app.parameters.freq_ind);
            plot_anat_image(app, app.parameters.anat_image_slice, app.ShowExcitationProfileCheckBox.Value);
        end

        % Value changed function: Freq1_ppmEditField
        function Freq1_ppmEditFieldValueChanged(app, event)
            value = app.Freq1_ppmEditField.Value;
            app.parameters.freq_ppm1 = value;
            plot_mag_vs_dist(app, app.parameters.freq_ind);
            plot_anat_image(app, app.parameters.anat_image_slice, app.ShowExcitationProfileCheckBox.Value);
        end

        % Value changed function: find_Gradients_lowEditField
        function find_Gradients_lowEditFieldValueChanged(app, event)
            app.parameters.find_grad_low = app.find_Gradients_lowEditField.Value;
        end

        % Value changed function: find_Gradients_highEditField
        function find_Gradients_highEditFieldValueChanged(app, event)
            app.parameters.find_grad_high = app.find_Gradients_highEditField.Value;
        end

        % Button pushed function: find_GradButton
        function find_GradButtonPushed(app, event)
            find_gradient_strength(app);
        end

        % Value changed function: find_Gradients_npointsEditField
        function find_Gradients_npointsEditFieldValueChanged(app, event)
            app.parameters.find_grad_npoints  = app.find_Gradients_npointsEditField.Value;
        end

        % Value changed function: repetition_time_msEditField
        function repetition_time_msEditFieldValueChanged(app, event)
            value = app.repetition_time_msEditField.Value;
            if (value <= app.parameters.pulse_duration_ms)
                value = app.parameters.pulse_duration_ms + 0.1;
                app.repetition_time_msEditField.Value = value;
            end
            app.parameters.repetitions_time_ms = value;
        end

        % Callback function
        function repetition_time_npointsEditFieldValueChanged(app, event)
            value = app.repetition_time_npointsEditField.Value;
            app.parameters.repetitions_time_ms_npoints = value;            
        end

        % Value changed function: nrepetitions_EditField
        function nrepetitions_EditFieldValueChanged(app, event)
            value = app.nrepetitions_EditField.Value;
            app.parameters.nrepetitions  = value;
            app.repSlider.Limits(2) = value;
        end

        % Button pushed function: SimulateRepsButton
        function SimulateRepsButtonPushed(app, event)
            app.repetition_time_npointsEditField.BackgroundColor = 'w';
            simulate_repetitions(app);
        end

        % Callback function
        function nrepetitions_EditFieldValueChanged2(app, event)
            value = app.nrepetitions_EditField.Value;
            app.parameters.nrepetitions = value;
            app.repetition_time_npointsEditField.BackgroundColor = 'w';
        end

        % Value changed function: find_Gradients_ind2EditField
        function find_Gradients_ind2EditFieldValueChanged(app, event)
            value = app.find_Gradients_ind2EditField.Value;
            if (value >= app.parameters.yrange_cm_npoints)
                app.parameters.find_grad_ind2 = app.parameters.yrange_cm_npoints;
            else
                app.parameters.find_grad_ind2 = value;
            end
            plot_mag_vs_dist(app, 1);
        end

        % Value changed function: find_Gradients_ind1EditField
        function find_Gradients_ind1EditFieldValueChanged(app, event)
            value = app.find_Gradients_ind1EditField.Value;
            if (value >= app.parameters.find_grad_ind2)
                app.parameters.find_grad_ind1 = app.parameters.find_grad_ind2;
            else
                app.parameters.find_grad_ind1 = value;
            end
            
            plot_mag_vs_dist(app, 1);
        end

        % Value changed function: repetition_time_npointsEditField
        function repetition_time_npointsEditFieldValueChanged2(app, event)
            value = app.repetition_time_npointsEditField.Value;
            app.parameters.repetitions_time_ms_npoints = value;
        end

        % Value changed function: HPEditField
        function HPEditFieldValueChanged(app, event)
            value = app.HPEditField.Value;
            app.parameters.sample_mag = value;            
        end

        % Value changed function: incphaseEditField
        function incphaseEditFieldValueChanged(app, event)
            value = app.incphaseEditField.Value;
            app.parameters.alternating_phase = value;            
        end

        % Value changed function: highresoutCheckBox
        function highresoutCheckBoxValueChanged(app, event)
            value = app.highresoutCheckBox.Value;
            app.parameters.save_all_sim_data = value ;
        end

        % Value changed function: amp_1stpulse_EditField
        function amp_1stpulse_EditFieldValueChanged(app, event)
            value = app.amp_1stpulse_EditField.Value;
            app.parameters.first_pulse_amp = value;
        end

        % Value changed function: inhomoB1ListBox
        function inhomoB1ListBoxValueChanged(app, event)
            value = app.inhomoB1ListBox.Value;
            if sum(strcmp(value, 'Mx'))
               app.parameters.plot_inhomo_mx  = 1;
            else
                app.parameters.plot_inhomo_mx  = 0;
            end
            if sum(strcmp(value, 'My'))
               app.parameters.plot_inhomo_my  = 1;
            else
                app.parameters.plot_inhomo_my  = 0;
            end
            if sum(strcmp(value, '|Mxy|'))
               app.parameters.plot_inhomo_mxy  = 1;
            else
                app.parameters.plot_inhomo_mxy  = 0;
            end
            if sum(strcmp(value, 'Mz'))
               app.parameters.plot_inhomo_mz  = 1;
            else
                app.parameters.plot_inhomo_mz  = 0;
            end
            if sum(strcmp(value, 'phi'))
               app.parameters.plot_inhomo_phase  = 1;
            else
                app.parameters.plot_inhomo_phase  = 0;
            end
            if sum(strcmp(value, 'none'))
               app.parameters.plot_inhomo_none = 1;
            else
                app.parameters.plot_inhomo_none = 0;
            end
            freq = app.freqSlider.Value;
            if app.parameters.plot_bssfp_signal
                plot_bssfp_vs_freq(app, freq);
            else
                plot_mag_vs_dist(app, freq);
            end
        end

        % Value changed function: homoB1ListBox
        function homoB1ListBoxValueChanged(app, event)
            value = app.homoB1ListBox.Value;
            if sum(strcmp(value, 'Mx'))
               app.parameters.plot_homo_mx  = 1;
            else
                app.parameters.plot_homo_mx  = 0;
            end
            if sum(strcmp(value, 'My'))
               app.parameters.plot_homo_my  = 1;
            else
                app.parameters.plot_homo_my  = 0;
            end
            if sum(strcmp(value, '|Mxy|'))
               app.parameters.plot_homo_mxy  = 1;
            else
                app.parameters.plot_homo_mxy  = 0;
            end
            if sum(strcmp(value, 'Mz'))
               app.parameters.plot_homo_mz  = 1;
            else
                app.parameters.plot_homo_mz  = 0;
            end
            if sum(strcmp(value, 'phi'))
               app.parameters.plot_homo_phase  = 1;
            else
                app.parameters.plot_homo_phase  = 0;
            end
            if sum(strcmp(value, 'none'))
               app.parameters.plot_homo_none = 1;
            else
                app.parameters.plot_homo_none = 0;
            end
            freq = app.freqSlider.Value;
            if app.parameters.plot_bssfp_signal
                plot_bssfp_vs_freq(app, freq);
            else
                plot_mag_vs_dist(app, freq);
            end
        end

        % Value changed function: b1_mag_khz_EditField
        function b1_mag_khz_EditFieldValueChanged(app, event)
            value = app.b1_mag_khz_EditField.Value;
            app.parameters.b1_amp_khz = value;
            calc_fa_b1_app(app, 0, 1);
        end

        % Value changed function: usemapCheckBox
        function usemapCheckBoxValueChanged(app, event)
            value = app.usemapCheckBox.Value;
            if (value == 1)
                calc_b1map(app);
                plot_b1map(app);
                app.refpow_map_caxis.Visible = 'on';
            else
                plot_b1profile(app);
                app.refpow_map_caxis.Visible = 'off';
            end
        end

        % Button down function: UIAxes_b1profile
        function UIAxes_b1profileButtonDown(app, event)
            pointer = get(gca, 'currentpoint');
            app.parameters.pointer2d = round(pointer(1, 1:2));
            if app.parameters.find_rfpower
                find_rfpower(app);
            end
        end

        % Value changed function: findRFPowerCheckBox
        function findRFPowerCheckBoxValueChanged(app, event)
            value = app.findRFPowerCheckBox.Value;
            app.parameters.find_rfpower = value;
        end

        % Button pushed function: calcrefpowButton
        function calcrefpowButtonPushed(app, event)
            calc_refpow(app);
        end

        % Value changed function: ymmEditField
        function ymmEditFieldValueChanged(app, event)
            value = app.ymmEditField.Value;
            app.parameters.refpow_pos_y_mm = value;
            calc_refpow(app);
            plot_b1map(app);
            try
                plot_anat_image(app, app.parameters.anat_image_slice, app.ShowExcitationProfileCheckBox.Value);
            end
        end

        % Value changed function: xmmEditField
        function xmmEditFieldValueChanged(app, event)
            value = app.xmmEditField.Value;
            app.parameters.refpow_pos_x_mm = value;
            calc_refpow(app);
            plot_b1map(app);
            try
                plot_anat_image(app, app.parameters.anat_image_slice, app.ShowExcitationProfileCheckBox.Value);
            end
        end

        % Value changed function: overlayrefpowCheckBox
        function overlayrefpowCheckBoxValueChanged(app, event)
            value = app.overlayrefpowCheckBox.Value;
            plot_anat_image(app, app.parameters.anat_image_slice, app.ShowExcitationProfileCheckBox.Value)
        end

        % Value changed function: refpow_map_caxis
        function refpow_map_caxisValueChanged(app, event)
            plot_b1map(app);
        end

        % Button pushed function: WButton
        function WButtonPushed(app, event)
            app.parameters.calc_refpow_outside_map = 1;
            calc_refpow(app);
            app.parameters.calc_refpow_outside_map = 0;
        end

        % Button pushed function: saveButton
        function saveButtonPushed(app, event)
            app.saveButton.BackgroundColor = [205, 0, 0]/256;
            app.saveButton.Text = 'saving';
            pause(0.000001);
            success = save_data(app);
            
            if success
                app.saveButton.BackgroundColor = [0, 179, 0]/256;
                app.saveButton.Text = 'save';
            else
                app.saveButton.Text = 'saving failed';
            end
        end

        % Callback function
        function duration_refocus_gradValueChanged(app, event)
            value = app.duration_refocus_grad.Value;
            % if the refocus time of the gradient is shorter than half the
            % time that's left when you subtract the pulse from the
            % repetition time:
            if (value <= (app.parameters.repetitions_time_ms - app.parameters.pulse_duration_ms) /2)
                app.parameters.ref_grad_dur_ms = value;
            else
                value = (app.parameters.repetitions_time_ms - app.parameters.pulse_duration_ms) /2;
                app.parameters.ref_grad_dur_ms = value;
                app.duration_refocus_grad.Value = value;
            end
        end

        % Value changed function: metabsCheckBox
        function metabsCheckBoxValueChanged(app, event)
            value = app.metabsCheckBox.Value;
            if value
                app.AlterNTimesEditField.Enable = 'on';
                app.AlterNTimesLabel.Enable = 'on';
                app.alpha2ndPulseEditField.Enable = 'on';
                app.freq_2ndpulse_HzEditField.Enable = 'on';
                app.parameters.simulate_alternating_freqs = 'on';
                app.parameters.save_all_sim_data = 'off';
                app.alpha2prepCheckBox.Enable = 'on';
            else
                app.AlterNTimesEditField.Enable = 'off';
                app.AlterNTimesLabel.Enable = 'off';
                app.parameters.alternate_freq_Ntimes_2nd_freq = 1;
                app.AlterNTimesEditField.Value = 1;
    
                app.alpha2ndPulseEditField.Enable = 'off';
                app.freq_2ndpulse_HzEditField.Enable = 'off';
                app.parameters.simulate_alternating_freqs = 'off';
                app.alpha2prepCheckBox.Enable = 'off';
            end

        end

        % Value changed function: first_repetition_time_msEditField
        function first_repetition_time_msEditFieldValueChanged(app, event)
            value = app.first_repetition_time_msEditField.Value;
            if (value > app.parameters.pulse_duration_ms)
                app.parameters.first_repetition_time_ms = value;
            else
                disp(['Min time >= ' num2str(app.parameters.pulse_duration_ms) 'ms']);
                app.first_repetition_time_msEditField.Value = app.parameters.pulse_duration_ms;
                app.parameters.first_repetition_time_ms = app.parameters.pulse_duration_ms;
            end
            
%             app.parameters.simulate_alternating_freqs = 'off';
%             app.parameters.alpha_2ndfreq = app.alpha2ndPulseEditField.Value;
%             app.parameters.freq_2ndfreq_hz = app.freq_2ndpulse_HzEditField.Value;
%             app.parameters.change_freq_2nd_freq = app.AlterNTimesEditField.Value;
        end

        % Value changed function: alpha2ndPulseEditField
        function alpha2ndPulseEditFieldValueChanged(app, event)
            value = app.alpha2ndPulseEditField.Value;
            app.parameters.alpha_2ndfreq = value;           
        end

        % Value changed function: freq_2ndpulse_HzEditField
        function freq_2ndpulse_HzEditFieldValueChanged(app, event)
            value = app.freq_2ndpulse_HzEditField.Value;
            app.parameters.freq_2ndfreq_hz = value;
        end

        % Value changed function: AlterNTimesEditField
        function AlterNTimesEditFieldValueChanged(app, event)
            value = app.AlterNTimesEditField.Value;
            app.parameters.alternate_freq_Ntimes_2nd_freq = value;
        end

        % Value changed function: freq_1stpulse_HzEditField
        function freq_1stpulse_HzEditFieldValueChanged(app, event)
            value = app.freq_1stpulse_HzEditField.Value;
            app.parameters.freq_1stfreq_hz = value;
        end

        % Value changed function: SpoilerGradCheckBox
        function SpoilerGradCheckBoxValueChanged(app, event)
            value = app.SpoilerGradCheckBox.Value;
            
            if value
%                 app.TipbackCheckBox.Value = false;
%                 app.TipbackCheckBox.Enable = "off";
%                 app.parameters.tipback_pulse = "off";

                app.SpoilerGradCheckBox.Value = true;
%                 app.SpoilerGradCheckBox.Enable = "on";
                app.parameters.spoiler_grad = value;
            else
                app.TipbackCheckBox.Enable = "on";
%                 app.parameters.tipback_pulse = "off";

%                 app.SpoilerGradCheckBox.Value = false;
%                 app.SpoilerGradCheckBox.Enable = "off";
                app.parameters.spoiler_grad = value;
            end
        end

        % Button pushed function: PlotresultsButton
        function PlotresultsButtonPushed(app, event)
            plot_bssfp_signal(app);
        end

        % Value changed function: altphaseEditField
        function altphaseEditFieldValueChanged(app, event)
            value = app.altphaseEditField.Value;
            app.parameters.phase_pulse_deg = value;
        end

        % Button pushed function: startparpoolButton
        function startparpoolButtonPushed(app, event)
            app.startparpoolButton.BackgroundColor = 'cyan';
            p = gcp('nocreate'); % If no pool, do not create new one.
            if isempty(p)
                poolsize = 0;
                try
                    parpool();
                    fprintf('\nparpool started.\n');
                    app.startparpoolButton.BackgroundColor = 'green';
                catch
                    fprintf('\nparpool could not be started.\n');
                end
            else
                poolsize = p.NumWorkers;
                fprintf('\nparpool already running with %d workers\n', poolsize);
            end
            
        end

        % Button pushed function: reset_B1map_Button
        function reset_B1map_ButtonPushed(app, event)
            % delete B1 map:
            disp('');
            try
                app.parameters.b1Transmit_sim = ones(size(app.parameters.b1Transmit_sim));
            end
            try
                app.parameters.refpow_map_header = ones(size(app.parameters.refpow_map_header));
            end
            try
                app.data.refpow_map_calc = ones(size(app.data.refpow_map_calc));
            end
            try
                app.data.b1TransmitProfileFit_log = zeros(size(app.data.b1TransmitProfileFit_log));
            end
            app.B1ProfileEditField.Value = '';

            % plot
            plot_b1profile(app);

            
        end

        % Value changed function: TipbackCheckBox
        function TipbackCheckBoxValueChanged(app, event)
            value = app.TipbackCheckBox.Value;
            if value
                app.TipbackCheckBox.Value = value;
                app.TipbackCheckBox.Enable = "on";
                app.parameters.tipback_pulse = "on";

%                 app.SpoilerGradCheckBox.Value = false;
%                 app.SpoilerGradCheckBox.Enable = "off";
%                 app.parameters.spoiler_grad = 0;
            else
%                 app.TipbackCheckBox.Value = "off";
                app.parameters.tipback_pulse = "off";

%                 app.SpoilerGradCheckBox.Enable = "on";
            end
        end

        % Value changed function: alpha2prepCheckBox
        function alpha2prepCheckBoxValueChanged(app, event)
            value = app.alpha2prepCheckBox.Value;
            app.parameters.first_pulse_always = value;
        end

        % Value changed function: bssfpsignalCheckBox
        function bssfpsignalCheckBoxValueChanged(app, event)
            value = app.bssfpsignalCheckBox.Value;
            app.parameters.plot_bssfp_signal = value;
        end

        % Value changed function: repSlider
        function repSliderValueChanged(app, event)
            value = round(app.repSlider.Value);
            plot_bssfp_vs_rep(app, value);
        end

        % Button pushed function: SimRangeButton
        function SimRangeButtonPushed(app, event)
            % Simulate Range of values, might be useful for finding best
            % settings:
            try 
               pathname = uigetdir('D:\Data\Simulations\test multi sims\');
            catch
               pathname = uigetdir('');
            end
            app.parameters.save_multi_sim_path = pathname;
            disp('');
            if (app.FARange_NEditField.Value > 1)
                fa_range = linspace(app.FARange_1EditField.Value, ...
                                    app.FARange_2EditField.Value, ...
                                    app.FARange_NEditField.Value);
            else
                fa_range = app.FARange_1EditField.Value;
            end

            if (app.FARange_NEditField.Value > 1)
            tr_range = linspace(app.TRRange_1EditField.Value, ...
                                app.TRRange_2EditField.Value, ...
                                app.TRRange_NEditField.Value);
            else
                fa_range = app.TRRange_1EditField.Value;
            end

            if (app.TrfRange_NEditField.Value > 1)
            trf_range = linspace(app.TrfRange_1EditField.Value, ...
                                 app.TrfRange_2EditField.Value, ...
                                 app.TrfRange_NEditField.Value);
            else
                fa_range = app.TrfRange_1EditField.Value;
            end

            if (app.T1Range_NEditField.Value > 1)
            t1_range = linspace(app.T1Range_1EditField.Value, ...
                                app.T1Range_2EditField.Value, ...
                                app.T1Range_NEditField.Value);

            else
                fa_range = app.T1Range_1EditField.Value;
            end

            if (app.T2Range_NEditField.Value > 1)
            t2_range = linspace(app.T2Range_1EditField.Value, ...
                                app.T2Range_2EditField.Value, ...
                                app.T2Range_NEditField.Value);
            else
                fa_range = app.T2Range_1EditField.Value;
            end

            ntot = numel(fa_range) * numel(tr_range)* numel(t1_range)* numel(t2_range);
            app.parameters.multi_sim_n = 0;
            n=0;
            for fa = 1:numel(fa_range)
            for tr = 1:numel(tr_range)
            for trf =1:numel(trf_range)
            for t1 = 1:numel(t1_range)
            for t2 = 1:numel(t2_range)
            if ((t1 > t2) && (tr > trf))
                app.parameters.fa_range_val = fa_range(fa);
                app.parameters.tr_range_val = tr_range(tr);
                app.parameters.t1_range_val = t1_range(t1);
                app.parameters.t2_range_val = t2_range(t2);
                % for now overwrite the set parameters:
                app.parameters.T1 = t1_range(t1);
                app.parameters.repetitions_time_ms = tr_range(tr);
                app.parameters.T2 = t2_range(t2);
                app.parameters.alpha = fa_range(fa);
                calc_fa_b1_app(app, 1, 0);
                app.parameters.pulse_duration_ms = trf_range(trf);
                calc_pulse_duration_ms_res(app);
                calc_slice_thick(app);
                calc_fa_b1_app(app, 1, 0);

                n=n+1;
                app.parameters.multi_sim_n = n;
                simulate_repetitions(app);
            end
            end
            end
            end
            end
            end
        end

        % Value changed function: simrangeCheckBox
        function simrangeCheckBoxValueChanged(app, event)
            value = app.simrangeCheckBox.Value;
            app.parameters.simulate_range_of_vals = value;
        end

        % Value changed function: TrfRange_NEditField
        function TrfRange_NEditFieldValueChanged(app, event)
            value = app.TrfRange_NEditField.Value;
            if value == 1
                app.TrfRange_2EditField.Enable = "off";
            else
                app.TrfRange_2EditField.Enable = "on";
            end
            
        end

        % Value changed function: FARange_NEditField
        function FARange_NEditFieldValueChanged(app, event)
            value = app.FARange_NEditField.Value;
            if value == 1
                app.FARange_2EditField.Enable = "off";
            else
                app.FARange_2EditField.Enable = "on";
            end
        end

        % Value changed function: TRRange_NEditField
        function TRRange_NEditFieldValueChanged(app, event)
            value = app.TRRange_NEditField.Value;
            if value == 1
                app.TRRange_2EditField.Enable = "off";
            else
                app.TRRange_2EditField.Enable = "on";
            end
        end

        % Value changed function: T2Range_NEditField
        function T2Range_NEditFieldValueChanged(app, event)
            value = app.T2Range_NEditField.Value;
            if value == 1
                app.T2Range_2EditField.Enable = "off";
            else
                app.T2Range_2EditField.Enable = "on";
            end
        end

        % Value changed function: T1Range_NEditField
        function T1Range_NEditFieldValueChanged(app, event)
            value = app.T1Range_NEditField.Value;
            if value == 1
                app.T1Range_2EditField.Enable = "off";
            else
                app.T1Range_2EditField.Enable = "on";
            end
        end

        % Value changed function: phasetbEditField
        function phasetbEditFieldValueChanged(app, event)
            value = app.phasetbEditField.Value;
            app.parameters.phase_tipback_pulse = value;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create calc_rf_pulse_toolUIFigure and hide until all components are created
            app.calc_rf_pulse_toolUIFigure = uifigure('Visible', 'off');
            app.calc_rf_pulse_toolUIFigure.AutoResizeChildren = 'off';
            app.calc_rf_pulse_toolUIFigure.Position = [100 100 1272 921];
            app.calc_rf_pulse_toolUIFigure.Name = 'RF pulse simulator';
            app.calc_rf_pulse_toolUIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.calc_rf_pulse_toolUIFigure);
            app.GridLayout.ColumnWidth = {445, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;
            app.LeftPanel.Scrollable = 'on';

            % Create RFPulseEditFieldLabel
            app.RFPulseEditFieldLabel = uilabel(app.LeftPanel);
            app.RFPulseEditFieldLabel.HorizontalAlignment = 'right';
            app.RFPulseEditFieldLabel.Position = [8 879 55 22];
            app.RFPulseEditFieldLabel.Text = 'RF Pulse';

            % Create RFPulseEditField
            app.RFPulseEditField = uieditfield(app.LeftPanel, 'text');
            app.RFPulseEditField.Tag = 'rfpulse_Button';
            app.RFPulseEditField.Position = [87 879 207 22];
            app.RFPulseEditField.Value = 'path to RF pulse ...';

            % Create B1ProfileEditFieldLabel
            app.B1ProfileEditFieldLabel = uilabel(app.LeftPanel);
            app.B1ProfileEditFieldLabel.HorizontalAlignment = 'right';
            app.B1ProfileEditFieldLabel.Position = [8 816 58 22];
            app.B1ProfileEditFieldLabel.Text = 'B1 Profile';

            % Create B1ProfileEditField
            app.B1ProfileEditField = uieditfield(app.LeftPanel, 'text');
            app.B1ProfileEditField.Position = [87 817 207 22];

            % Create SelectRFPulse_Button
            app.SelectRFPulse_Button = uibutton(app.LeftPanel, 'push');
            app.SelectRFPulse_Button.ButtonPushedFcn = createCallbackFcn(app, @SelectRFPulse_ButtonPushed, true);
            app.SelectRFPulse_Button.Tag = 'selectRF_Button';
            app.SelectRFPulse_Button.Position = [87 858 52 22];
            app.SelectRFPulse_Button.Text = 'Select';

            % Create LoadRFPulse_Button
            app.LoadRFPulse_Button = uibutton(app.LeftPanel, 'push');
            app.LoadRFPulse_Button.ButtonPushedFcn = createCallbackFcn(app, @LoadRFPulse_ButtonPushed, true);
            app.LoadRFPulse_Button.Tag = 'loadRF_Button';
            app.LoadRFPulse_Button.Position = [138 858 52 22];
            app.LoadRFPulse_Button.Text = 'Load';

            % Create SelectB1Profile_Button
            app.SelectB1Profile_Button = uibutton(app.LeftPanel, 'push');
            app.SelectB1Profile_Button.ButtonPushedFcn = createCallbackFcn(app, @SelectB1Profile_ButtonPushed, true);
            app.SelectB1Profile_Button.Position = [87 796 52 22];
            app.SelectB1Profile_Button.Text = 'Select';

            % Create LoadB1Profile_Button
            app.LoadB1Profile_Button = uibutton(app.LeftPanel, 'push');
            app.LoadB1Profile_Button.ButtonPushedFcn = createCallbackFcn(app, @LoadB1Profile_ButtonPushed, true);
            app.LoadB1Profile_Button.Position = [138 796 52 22];
            app.LoadB1Profile_Button.Text = 'Load';

            % Create alphaLabel
            app.alphaLabel = uilabel(app.LeftPanel);
            app.alphaLabel.HorizontalAlignment = 'right';
            app.alphaLabel.Position = [317 879 50 22];
            app.alphaLabel.Text = 'alpha [°]';

            % Create flipangle_degEditField
            app.flipangle_degEditField = uieditfield(app.LeftPanel, 'numeric');
            app.flipangle_degEditField.Limits = [0 Inf];
            app.flipangle_degEditField.ValueChangedFcn = createCallbackFcn(app, @flipangle_degEditFieldValueChanged, true);
            app.flipangle_degEditField.Tag = 'pulse_duration_Edit';
            app.flipangle_degEditField.Position = [370 879 69 22];
            app.flipangle_degEditField.Value = 10;

            % Create alphaLabel_2
            app.alphaLabel_2 = uilabel(app.LeftPanel);
            app.alphaLabel_2.HorizontalAlignment = 'right';
            app.alphaLabel_2.Position = [294 837 75 22];
            app.alphaLabel_2.Text = 'BWfac [Hz s]';

            % Create sint_EditField
            app.sint_EditField = uieditfield(app.LeftPanel, 'numeric');
            app.sint_EditField.Limits = [0.001 100];
            app.sint_EditField.ValueChangedFcn = createCallbackFcn(app, @sint_EditFieldValueChanged, true);
            app.sint_EditField.Tag = 'pulse_duration_Edit';
            app.sint_EditField.Position = [370 816 69 22];
            app.sint_EditField.Value = 0.2518;

            % Create bwfac_HzsEditField
            app.bwfac_HzsEditField = uieditfield(app.LeftPanel, 'numeric');
            app.bwfac_HzsEditField.Limits = [0.001 Inf];
            app.bwfac_HzsEditField.Tag = 'pulse_duration_Edit';
            app.bwfac_HzsEditField.Position = [370 837 69 22];
            app.bwfac_HzsEditField.Value = 3972;

            % Create alphaLabel_3
            app.alphaLabel_3 = uilabel(app.LeftPanel);
            app.alphaLabel_3.HorizontalAlignment = 'right';
            app.alphaLabel_3.Position = [339 816 26 22];
            app.alphaLabel_3.Text = 'Sint';

            % Create AnatImageEditFieldLabel
            app.AnatImageEditFieldLabel = uilabel(app.LeftPanel);
            app.AnatImageEditFieldLabel.HorizontalAlignment = 'right';
            app.AnatImageEditFieldLabel.Position = [8 755 67 22];
            app.AnatImageEditFieldLabel.Text = 'Anat Image';

            % Create AnatImageEditField
            app.AnatImageEditField = uieditfield(app.LeftPanel, 'text');
            app.AnatImageEditField.Position = [87 756 207 22];

            % Create LoadAnatImage_Button
            app.LoadAnatImage_Button = uibutton(app.LeftPanel, 'push');
            app.LoadAnatImage_Button.ButtonPushedFcn = createCallbackFcn(app, @LoadAnatImage_ButtonPushed, true);
            app.LoadAnatImage_Button.Position = [138 735 52 22];
            app.LoadAnatImage_Button.Text = 'Load';

            % Create calcButton
            app.calcButton = uibutton(app.LeftPanel, 'push');
            app.calcButton.ButtonPushedFcn = createCallbackFcn(app, @calcButtonPushed, true);
            app.calcButton.Position = [374 769 62 22];
            app.calcButton.Text = 'calc !';

            % Create SelectAnatImage_Button
            app.SelectAnatImage_Button = uibutton(app.LeftPanel, 'push');
            app.SelectAnatImage_Button.ButtonPushedFcn = createCallbackFcn(app, @SelectAnatImage_ButtonPushed, true);
            app.SelectAnatImage_Button.Position = [87 735 52 22];
            app.SelectAnatImage_Button.Text = 'Select';

            % Create SamplePanel
            app.SamplePanel = uipanel(app.LeftPanel);
            app.SamplePanel.ForegroundColor = [1 0 0];
            app.SamplePanel.Title = 'Sample';
            app.SamplePanel.Position = [11 361 428 109];

            % Create NucleusButtonGroup
            app.NucleusButtonGroup = uibuttongroup(app.SamplePanel);
            app.NucleusButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @NucleusButtonGroupSelectionChanged, true);
            app.NucleusButtonGroup.Enable = 'off';
            app.NucleusButtonGroup.Title = 'Nucleus';
            app.NucleusButtonGroup.Visible = 'off';
            app.NucleusButtonGroup.Tag = 'nuc_select_Buttongroup';
            app.NucleusButtonGroup.Position = [11 8 97 70];

            % Create HButton
            app.HButton = uiradiobutton(app.NucleusButtonGroup);
            app.HButton.Tag = 'HButton';
            app.HButton.Text = {'1H'; ''};
            app.HButton.Position = [11 24 58 22];

            % Create CButton
            app.CButton = uiradiobutton(app.NucleusButtonGroup);
            app.CButton.Tag = 'HButton';
            app.CButton.Text = '13C';
            app.CButton.Position = [11 2 65 22];
            app.CButton.Value = true;

            % Create T1sEditFieldLabel
            app.T1sEditFieldLabel = uilabel(app.SamplePanel);
            app.T1sEditFieldLabel.HorizontalAlignment = 'right';
            app.T1sEditFieldLabel.Position = [117 43 35 22];
            app.T1sEditFieldLabel.Text = 'T1 [s]';

            % Create T1_sEditField
            app.T1_sEditField = uieditfield(app.SamplePanel, 'text');
            app.T1_sEditField.ValueChangedFcn = createCallbackFcn(app, @T1_sEditFieldValueChanged, true);
            app.T1_sEditField.Position = [168 43 45 22];
            app.T1_sEditField.Value = '30';

            % Create T2sEditFieldLabel
            app.T2sEditFieldLabel = uilabel(app.SamplePanel);
            app.T2sEditFieldLabel.HorizontalAlignment = 'right';
            app.T2sEditFieldLabel.Position = [118 11 35 22];
            app.T2sEditFieldLabel.Text = 'T2 [s]';

            % Create T2_sEditField
            app.T2_sEditField = uieditfield(app.SamplePanel, 'text');
            app.T2_sEditField.ValueChangedFcn = createCallbackFcn(app, @T2_sEditFieldValueChanged, true);
            app.T2_sEditField.Position = [169 11 44 22];
            app.T2_sEditField.Value = '10';

            % Create HPEditFieldLabel
            app.HPEditFieldLabel = uilabel(app.SamplePanel);
            app.HPEditFieldLabel.HorizontalAlignment = 'right';
            app.HPEditFieldLabel.Position = [232 43 25 22];
            app.HPEditFieldLabel.Text = 'HP';

            % Create HPEditField
            app.HPEditField = uieditfield(app.SamplePanel, 'numeric');
            app.HPEditField.ValueChangedFcn = createCallbackFcn(app, @HPEditFieldValueChanged, true);
            app.HPEditField.Tooltip = {'Polarization Level'};
            app.HPEditField.Position = [271 43 55 22];
            app.HPEditField.Value = 1;

            % Create SimulationParametersPanel
            app.SimulationParametersPanel = uipanel(app.LeftPanel);
            app.SimulationParametersPanel.Title = 'Simulation Parameters';
            app.SimulationParametersPanel.Position = [11 485 428 216];

            % Create PulseDurationmsEditFieldLabel
            app.PulseDurationmsEditFieldLabel = uilabel(app.SimulationParametersPanel);
            app.PulseDurationmsEditFieldLabel.HorizontalAlignment = 'right';
            app.PulseDurationmsEditFieldLabel.Position = [3 145 110 22];
            app.PulseDurationmsEditFieldLabel.Text = 'Pulse Duration [ms]';

            % Create pulseduration_msEditField
            app.pulseduration_msEditField = uieditfield(app.SimulationParametersPanel, 'numeric');
            app.pulseduration_msEditField.Limits = [0.001 100];
            app.pulseduration_msEditField.ValueChangedFcn = createCallbackFcn(app, @pulseduration_msEditFieldValueChanged, true);
            app.pulseduration_msEditField.Tag = 'pulse_duration_Edit';
            app.pulseduration_msEditField.Position = [157 145 79 22];
            app.pulseduration_msEditField.Value = 5;

            % Create SimulateButton
            app.SimulateButton = uibutton(app.SimulationParametersPanel, 'push');
            app.SimulateButton.ButtonPushedFcn = createCallbackFcn(app, @SimulateButtonPushed, true);
            app.SimulateButton.Tag = 'simulate_Button';
            app.SimulateButton.Position = [296 10 125 22];
            app.SimulateButton.Text = 'Simulate !';

            % Create GradientStrengthkHzcmEditFieldLabel
            app.GradientStrengthkHzcmEditFieldLabel = uilabel(app.SimulationParametersPanel);
            app.GradientStrengthkHzcmEditFieldLabel.HorizontalAlignment = 'right';
            app.GradientStrengthkHzcmEditFieldLabel.Position = [3 113 150 22];
            app.GradientStrengthkHzcmEditFieldLabel.Text = 'Gradient Strength [kHz/cm]';

            % Create gradientstrength_khzcmEditField
            app.gradientstrength_khzcmEditField = uieditfield(app.SimulationParametersPanel, 'numeric');
            app.gradientstrength_khzcmEditField.Limits = [-10000 10000];
            app.gradientstrength_khzcmEditField.ValueDisplayFormat = '%11.7g';
            app.gradientstrength_khzcmEditField.ValueChangedFcn = createCallbackFcn(app, @gradientstrength_khzcmEditFieldValueChanged, true);
            app.gradientstrength_khzcmEditField.Tag = 'grad_strength_Edit';
            app.gradientstrength_khzcmEditField.Position = [157 112 79 22];

            % Create ycmLabel
            app.ycmLabel = uilabel(app.SimulationParametersPanel);
            app.ycmLabel.HorizontalAlignment = 'right';
            app.ycmLabel.Position = [3 77 37 22];
            app.ycmLabel.Text = 'y [cm]';

            % Create yrange_cmEditField
            app.yrange_cmEditField = uieditfield(app.SimulationParametersPanel, 'numeric');
            app.yrange_cmEditField.ValueChangedFcn = createCallbackFcn(app, @yrange_cmEditFieldValueChanged, true);
            app.yrange_cmEditField.Tag = 'grad_strength_Edit';
            app.yrange_cmEditField.Position = [157 77 79 22];
            app.yrange_cmEditField.Value = 3.2;

            % Create pointsLabel
            app.pointsLabel = uilabel(app.SimulationParametersPanel);
            app.pointsLabel.Position = [270 166 38 22];
            app.pointsLabel.Text = 'points';

            % Create pulseduration_ms_npointsEditField
            app.pulseduration_ms_npointsEditField = uieditfield(app.SimulationParametersPanel, 'text');
            app.pulseduration_ms_npointsEditField.ValueChangedFcn = createCallbackFcn(app, @pulseduration_ms_npointsEditFieldValueChanged, true);
            app.pulseduration_ms_npointsEditField.Tag = 'pulse_duration_res';
            app.pulseduration_ms_npointsEditField.Position = [270 145 53 22];
            app.pulseduration_ms_npointsEditField.Value = '1000';

            % Create EditField_2
            app.EditField_2 = uieditfield(app.SimulationParametersPanel, 'text');
            app.EditField_2.Enable = 'off';
            app.EditField_2.Position = [270 112 53 22];
            app.EditField_2.Value = '1';

            % Create yrange_cm_npointsEditField
            app.yrange_cm_npointsEditField = uieditfield(app.SimulationParametersPanel, 'text');
            app.yrange_cm_npointsEditField.ValueChangedFcn = createCallbackFcn(app, @yrange_cm_npointsEditFieldValueChanged, true);
            app.yrange_cm_npointsEditField.Position = [270 77 53 22];
            app.yrange_cm_npointsEditField.Value = '33';

            % Create freqHzEditFieldLabel
            app.freqHzEditFieldLabel = uilabel(app.SimulationParametersPanel);
            app.freqHzEditFieldLabel.HorizontalAlignment = 'right';
            app.freqHzEditFieldLabel.Position = [3 43 51 22];
            app.freqHzEditFieldLabel.Text = 'freq [Hz]';

            % Create freq_hz_lowEditField
            app.freq_hz_lowEditField = uieditfield(app.SimulationParametersPanel, 'numeric');
            app.freq_hz_lowEditField.ValueChangedFcn = createCallbackFcn(app, @freq_hz_lowEditFieldValueChanged, true);
            app.freq_hz_lowEditField.Tag = 'grad_strength_Edit';
            app.freq_hz_lowEditField.Position = [124 43 62 22];
            app.freq_hz_lowEditField.Value = -1000;

            % Create freq_hz_highEditField
            app.freq_hz_highEditField = uieditfield(app.SimulationParametersPanel, 'numeric');
            app.freq_hz_highEditField.ValueChangedFcn = createCallbackFcn(app, @freq_hz_highEditFieldValueChanged, true);
            app.freq_hz_highEditField.Tag = 'grad_strength_Edit';
            app.freq_hz_highEditField.Position = [190 43 62 22];
            app.freq_hz_highEditField.Value = 1000;

            % Create freq_hz_npointsEditField
            app.freq_hz_npointsEditField = uieditfield(app.SimulationParametersPanel, 'text');
            app.freq_hz_npointsEditField.ValueChangedFcn = createCallbackFcn(app, @freq_hz_npointsEditFieldValueChanged, true);
            app.freq_hz_npointsEditField.Position = [270 43 54 22];
            app.freq_hz_npointsEditField.Value = '11';

            % Create resolutionLabel
            app.resolutionLabel = uilabel(app.SimulationParametersPanel);
            app.resolutionLabel.Position = [343 166 58 22];
            app.resolutionLabel.Text = 'resolution';

            % Create pulseduration_ms_resLabel
            app.pulseduration_ms_resLabel = uilabel(app.SimulationParametersPanel);
            app.pulseduration_ms_resLabel.Position = [343 145 76 22];

            % Create slice_thick
            app.slice_thick = uilabel(app.SimulationParametersPanel);
            app.slice_thick.Enable = 'off';
            app.slice_thick.Tooltip = {'slicethickness'};
            app.slice_thick.Position = [342 112 77 22];
            app.slice_thick.Text = '0';

            % Create yrange_cm_resLabel
            app.yrange_cm_resLabel = uilabel(app.SimulationParametersPanel);
            app.yrange_cm_resLabel.Position = [343 77 76 22];

            % Create freq_hz_resLabel
            app.freq_hz_resLabel = uilabel(app.SimulationParametersPanel);
            app.freq_hz_resLabel.Position = [343 43 82 22];

            % Create valueLabel
            app.valueLabel = uilabel(app.SimulationParametersPanel);
            app.valueLabel.Position = [157 166 34 22];
            app.valueLabel.Text = 'value';

            % Create alphaLabel_4
            app.alphaLabel_4 = uilabel(app.LeftPanel);
            app.alphaLabel_4.HorizontalAlignment = 'right';
            app.alphaLabel_4.Position = [317 858 51 22];
            app.alphaLabel_4.Text = 'B1 [kHz]';

            % Create b1_mag_khz_EditField
            app.b1_mag_khz_EditField = uieditfield(app.LeftPanel, 'numeric');
            app.b1_mag_khz_EditField.ValueChangedFcn = createCallbackFcn(app, @b1_mag_khz_EditFieldValueChanged, true);
            app.b1_mag_khz_EditField.Tag = 'pulse_duration_Edit';
            app.b1_mag_khz_EditField.Position = [370 858 69 22];

            % Create TabGroup
            app.TabGroup = uitabgroup(app.LeftPanel);
            app.TabGroup.Position = [11 1 428 172];

            % Create findGradTab
            app.findGradTab = uitab(app.TabGroup);
            app.findGradTab.Title = 'find Grad';

            % Create findGradientstrengthforB1compensationPanel
            app.findGradientstrengthforB1compensationPanel = uipanel(app.findGradTab);
            app.findGradientstrengthforB1compensationPanel.ForegroundColor = [0 0 1];
            app.findGradientstrengthforB1compensationPanel.Title = 'find Gradient strength for B1 compensation:';
            app.findGradientstrengthforB1compensationPanel.Position = [4 16 424 116];

            % Create find_Gradients_lowEditField
            app.find_Gradients_lowEditField = uieditfield(app.findGradientstrengthforB1compensationPanel, 'numeric');
            app.find_Gradients_lowEditField.ValueChangedFcn = createCallbackFcn(app, @find_Gradients_lowEditFieldValueChanged, true);
            app.find_Gradients_lowEditField.Tag = 'grad_strength_Edit';
            app.find_Gradients_lowEditField.Position = [12 47 51 22];

            % Create find_Gradients_highEditField
            app.find_Gradients_highEditField = uieditfield(app.findGradientstrengthforB1compensationPanel, 'numeric');
            app.find_Gradients_highEditField.ValueChangedFcn = createCallbackFcn(app, @find_Gradients_highEditFieldValueChanged, true);
            app.find_Gradients_highEditField.Tag = 'grad_strength_Edit';
            app.find_Gradients_highEditField.Position = [89 47 51 22];
            app.find_Gradients_highEditField.Value = 1;

            % Create find_Gradients_npointsEditField
            app.find_Gradients_npointsEditField = uieditfield(app.findGradientstrengthforB1compensationPanel, 'numeric');
            app.find_Gradients_npointsEditField.ValueChangedFcn = createCallbackFcn(app, @find_Gradients_npointsEditFieldValueChanged, true);
            app.find_Gradients_npointsEditField.Position = [166 47 49 22];
            app.find_Gradients_npointsEditField.Value = 11;

            % Create BestGradEditFieldLabel
            app.BestGradEditFieldLabel = uilabel(app.findGradientstrengthforB1compensationPanel);
            app.BestGradEditFieldLabel.HorizontalAlignment = 'right';
            app.BestGradEditFieldLabel.Position = [304 66 60 22];
            app.BestGradEditFieldLabel.Text = 'Best Grad';

            % Create BestGradEditField
            app.BestGradEditField = uieditfield(app.findGradientstrengthforB1compensationPanel, 'numeric');
            app.BestGradEditField.ValueDisplayFormat = '%11.7g';
            app.BestGradEditField.Editable = 'off';
            app.BestGradEditField.Position = [310 47 74 22];

            % Create find_Gradients_ind1EditField
            app.find_Gradients_ind1EditField = uieditfield(app.findGradientstrengthforB1compensationPanel, 'numeric');
            app.find_Gradients_ind1EditField.Limits = [1 Inf];
            app.find_Gradients_ind1EditField.ValueChangedFcn = createCallbackFcn(app, @find_Gradients_ind1EditFieldValueChanged, true);
            app.find_Gradients_ind1EditField.Position = [231 47 28 22];
            app.find_Gradients_ind1EditField.Value = 1;

            % Create find_Gradients_ind2EditField
            app.find_Gradients_ind2EditField = uieditfield(app.findGradientstrengthforB1compensationPanel, 'numeric');
            app.find_Gradients_ind2EditField.Limits = [2 Inf];
            app.find_Gradients_ind2EditField.ValueChangedFcn = createCallbackFcn(app, @find_Gradients_ind2EditFieldValueChanged, true);
            app.find_Gradients_ind2EditField.Position = [261 47 28 22];
            app.find_Gradients_ind2EditField.Value = 16;

            % Create G1kHzcmLabel
            app.G1kHzcmLabel = uilabel(app.findGradientstrengthforB1compensationPanel);
            app.G1kHzcmLabel.Position = [12 66 71 22];
            app.G1kHzcmLabel.Text = 'G1 [kHz/cm]';

            % Create G2kHzcmLabel
            app.G2kHzcmLabel = uilabel(app.findGradientstrengthforB1compensationPanel);
            app.G2kHzcmLabel.Position = [89 66 71 22];
            app.G2kHzcmLabel.Text = 'G2 [kHz/cm]';

            % Create NpointsLabel_2
            app.NpointsLabel_2 = uilabel(app.findGradientstrengthforB1compensationPanel);
            app.NpointsLabel_2.Position = [166 66 56 22];
            app.NpointsLabel_2.Text = 'N [points]';

            % Create pointsyLabel
            app.pointsyLabel = uilabel(app.findGradientstrengthforB1compensationPanel);
            app.pointsyLabel.Position = [233 66 54 22];
            app.pointsyLabel.Text = 'points [y]';

            % Create progress_grad_findLabel
            app.progress_grad_findLabel = uilabel(app.findGradientstrengthforB1compensationPanel);
            app.progress_grad_findLabel.Position = [169 10 87 22];
            app.progress_grad_findLabel.Text = '';

            % Create find_GradButton
            app.find_GradButton = uibutton(app.findGradTab, 'push');
            app.find_GradButton.ButtonPushedFcn = createCallbackFcn(app, @find_GradButtonPushed, true);
            app.find_GradButton.Tag = 'simulate_Button';
            app.find_GradButton.Position = [293 30 125 22];
            app.find_GradButton.Text = 'find Grad !';

            % Create RFpowerTab
            app.RFpowerTab = uitab(app.TabGroup);
            app.RFpowerTab.Title = 'RF power';

            % Create usemapCheckBox
            app.usemapCheckBox = uicheckbox(app.RFpowerTab);
            app.usemapCheckBox.ValueChangedFcn = createCallbackFcn(app, @usemapCheckBoxValueChanged, true);
            app.usemapCheckBox.Text = 'use map';
            app.usemapCheckBox.Position = [131 61 68 22];

            % Create RFPowerWEditField
            app.RFPowerWEditField = uieditfield(app.RFpowerTab, 'numeric');
            app.RFPowerWEditField.Position = [102 84 100 22];

            % Create RFPowerLabel
            app.RFPowerLabel = uilabel(app.RFpowerTab);
            app.RFPowerLabel.Position = [29 47 59 22];
            app.RFPowerLabel.Text = 'RF Power';

            % Create RFPowerWEditFieldLabel
            app.RFPowerWEditFieldLabel = uilabel(app.RFpowerTab);
            app.RFPowerWEditFieldLabel.HorizontalAlignment = 'right';
            app.RFPowerWEditFieldLabel.Position = [7 84 80 22];
            app.RFPowerWEditFieldLabel.Text = 'RF Power [W]';

            % Create findRFPowerCheckBox
            app.findRFPowerCheckBox = uicheckbox(app.RFpowerTab);
            app.findRFPowerCheckBox.ValueChangedFcn = createCallbackFcn(app, @findRFPowerCheckBoxValueChanged, true);
            app.findRFPowerCheckBox.Text = 'find RF Power';
            app.findRFPowerCheckBox.Position = [132 35 98 22];

            % Create ymmEditFieldLabel
            app.ymmEditFieldLabel = uilabel(app.RFpowerTab);
            app.ymmEditFieldLabel.HorizontalAlignment = 'right';
            app.ymmEditFieldLabel.Position = [364 110 41 22];
            app.ymmEditFieldLabel.Text = 'y [mm]';

            % Create ymmEditField
            app.ymmEditField = uieditfield(app.RFpowerTab, 'numeric');
            app.ymmEditField.ValueChangedFcn = createCallbackFcn(app, @ymmEditFieldValueChanged, true);
            app.ymmEditField.Position = [368 87 37 22];

            % Create calcrefpowButton
            app.calcrefpowButton = uibutton(app.RFpowerTab, 'push');
            app.calcrefpowButton.ButtonPushedFcn = createCallbackFcn(app, @calcrefpowButtonPushed, true);
            app.calcrefpowButton.Position = [320 61 68 22];
            app.calcrefpowButton.Text = 'calc refpow';

            % Create RefpowLabel
            app.RefpowLabel = uilabel(app.RFpowerTab);
            app.RefpowLabel.HorizontalAlignment = 'center';
            app.RefpowLabel.Position = [320 40 99 22];
            app.RefpowLabel.Text = 'Refpow';

            % Create xmmEditFieldLabel
            app.xmmEditFieldLabel = uilabel(app.RFpowerTab);
            app.xmmEditFieldLabel.HorizontalAlignment = 'right';
            app.xmmEditFieldLabel.Position = [320 110 41 22];
            app.xmmEditFieldLabel.Text = 'x [mm]';

            % Create xmmEditField
            app.xmmEditField = uieditfield(app.RFpowerTab, 'numeric');
            app.xmmEditField.ValueChangedFcn = createCallbackFcn(app, @xmmEditFieldValueChanged, true);
            app.xmmEditField.Position = [324 87 37 22];

            % Create overlayrefpowCheckBox
            app.overlayrefpowCheckBox = uicheckbox(app.RFpowerTab);
            app.overlayrefpowCheckBox.ValueChangedFcn = createCallbackFcn(app, @overlayrefpowCheckBoxValueChanged, true);
            app.overlayrefpowCheckBox.Text = 'overlay refpow';
            app.overlayrefpowCheckBox.Position = [328 19 100 22];

            % Create WButton
            app.WButton = uibutton(app.RFpowerTab, 'push');
            app.WButton.ButtonPushedFcn = createCallbackFcn(app, @WButtonPushed, true);
            app.WButton.Position = [391 61 34 22];
            app.WButton.Text = '0W';

            % Create saveTab
            app.saveTab = uitab(app.TabGroup);
            app.saveTab.Title = 'save';

            % Create saveoptions_ListBox
            app.saveoptions_ListBox = uilistbox(app.saveTab);
            app.saveoptions_ListBox.Items = {'time', 'freqs', 'dist', 'reps', 'single', 'plots'};
            app.saveoptions_ListBox.Multiselect = 'on';
            app.saveoptions_ListBox.Tooltip = {'You can use "Ctrl" to (un-) select options. For now it''s recommended to just select all the options'};
            app.saveoptions_ListBox.Position = [16 35 100 97];
            app.saveoptions_ListBox.Value = {'time', 'freqs', 'dist', 'reps', 'single', 'plots'};

            % Create repEditFieldLabel
            app.repEditFieldLabel = uilabel(app.saveTab);
            app.repEditFieldLabel.HorizontalAlignment = 'right';
            app.repEditFieldLabel.Position = [157 49 25 22];
            app.repEditFieldLabel.Text = 'rep';

            % Create repEditField_low
            app.repEditField_low = uieditfield(app.saveTab, 'numeric');
            app.repEditField_low.Limits = [1 Inf];
            app.repEditField_low.Position = [197 49 29 22];
            app.repEditField_low.Value = 1;

            % Create repEditField_high
            app.repEditField_high = uieditfield(app.saveTab, 'numeric');
            app.repEditField_high.Limits = [1 Inf];
            app.repEditField_high.Position = [231 49 29 22];
            app.repEditField_high.Value = 1;

            % Create saveButton
            app.saveButton = uibutton(app.saveTab, 'push');
            app.saveButton.ButtonPushedFcn = createCallbackFcn(app, @saveButtonPushed, true);
            app.saveButton.Position = [319 5 100 22];
            app.saveButton.Text = 'save';

            % Create fileEditFieldLabel
            app.fileEditFieldLabel = uilabel(app.saveTab);
            app.fileEditFieldLabel.HorizontalAlignment = 'right';
            app.fileEditFieldLabel.Position = [203 89 25 22];
            app.fileEditFieldLabel.Text = 'file';

            % Create save_fileEditField
            app.save_fileEditField = uieditfield(app.saveTab, 'text');
            app.save_fileEditField.Position = [243 89 100 22];

            % Create nslices
            app.nslices = uilabel(app.LeftPanel);
            app.nslices.FontSize = 6;
            app.nslices.Position = [205 735 25 22];
            app.nslices.Text = 'nslices';

            % Create nslices_2
            app.nslices_2 = uilabel(app.LeftPanel);
            app.nslices_2.FontSize = 6;
            app.nslices_2.Position = [231 735 25 22];
            app.nslices_2.Text = {'mat'; ''};

            % Create nslices_3
            app.nslices_3 = uilabel(app.LeftPanel);
            app.nslices_3.FontSize = 6;
            app.nslices_3.Position = [257 735 25 22];
            app.nslices_3.Text = 'offset';

            % Create SimulateRepetitionsTabGroup
            app.SimulateRepetitionsTabGroup = uitabgroup(app.LeftPanel);
            app.SimulateRepetitionsTabGroup.Position = [12 184 418 136];

            % Create BasicParameterTab
            app.BasicParameterTab = uitab(app.SimulateRepetitionsTabGroup);
            app.BasicParameterTab.Title = 'Basic parameters';

            % Create repetition_time_msEditField
            app.repetition_time_msEditField = uieditfield(app.BasicParameterTab, 'numeric');
            app.repetition_time_msEditField.Limits = [0.001 100];
            app.repetition_time_msEditField.ValueChangedFcn = createCallbackFcn(app, @repetition_time_msEditFieldValueChanged, true);
            app.repetition_time_msEditField.Tag = 'grad_strength_Edit';
            app.repetition_time_msEditField.Tooltip = {'Repetition Time'};
            app.repetition_time_msEditField.Position = [10 62 54 22];
            app.repetition_time_msEditField.Value = 10;

            % Create nrepetitions_EditField
            app.nrepetitions_EditField = uieditfield(app.BasicParameterTab, 'numeric');
            app.nrepetitions_EditField.Limits = [1 Inf];
            app.nrepetitions_EditField.ValueChangedFcn = createCallbackFcn(app, @nrepetitions_EditFieldValueChanged, true);
            app.nrepetitions_EditField.Tooltip = {'Number of phase encodes'};
            app.nrepetitions_EditField.Position = [125 62 46 22];
            app.nrepetitions_EditField.Value = 100;

            % Create repetition_time_npointsEditField
            app.repetition_time_npointsEditField = uieditfield(app.BasicParameterTab, 'numeric');
            app.repetition_time_npointsEditField.ValueChangedFcn = createCallbackFcn(app, @repetition_time_npointsEditFieldValueChanged2, true);
            app.repetition_time_npointsEditField.Tooltip = {'Number of points that the go into the Bloch Simulation.'};
            app.repetition_time_npointsEditField.Position = [68 62 52 22];
            app.repetition_time_npointsEditField.Value = 512;

            % Create incphaseEditField
            app.incphaseEditField = uieditfield(app.BasicParameterTab, 'numeric');
            app.incphaseEditField.ValueChangedFcn = createCallbackFcn(app, @incphaseEditFieldValueChanged, true);
            app.incphaseEditField.Position = [178 62 72 22];
            app.incphaseEditField.Value = 180;

            % Create TRmsLabel
            app.TRmsLabel = uilabel(app.BasicParameterTab);
            app.TRmsLabel.Position = [10 82 47 22];
            app.TRmsLabel.Text = 'TR [ms]';

            % Create NpointsLabel_3
            app.NpointsLabel_3 = uilabel(app.BasicParameterTab);
            app.NpointsLabel_3.Position = [68 82 56 22];
            app.NpointsLabel_3.Text = 'N [points]';

            % Create NRepsLabel
            app.NRepsLabel = uilabel(app.BasicParameterTab);
            app.NRepsLabel.Position = [125 81 52 22];
            app.NRepsLabel.Text = 'N [Reps]';

            % Create incphaseLabel
            app.incphaseLabel = uilabel(app.BasicParameterTab);
            app.incphaseLabel.HorizontalAlignment = 'right';
            app.incphaseLabel.Position = [172 81 75 22];
            app.incphaseLabel.Text = 'inc. phase [°]';

            % Create progress_simRepLabel
            app.progress_simRepLabel = uilabel(app.BasicParameterTab);
            app.progress_simRepLabel.BackgroundColor = [0.8 0.8 0.8];
            app.progress_simRepLabel.Position = [272 104 5 5];
            app.progress_simRepLabel.Text = '';

            % Create progress_simRepLabel_2
            app.progress_simRepLabel_2 = uilabel(app.BasicParameterTab);
            app.progress_simRepLabel_2.BackgroundColor = [0.8 0.8 0.8];
            app.progress_simRepLabel_2.Position = [289 104 5 5];
            app.progress_simRepLabel_2.Text = '';

            % Create progress_simRepLabel_3
            app.progress_simRepLabel_3 = uilabel(app.BasicParameterTab);
            app.progress_simRepLabel_3.BackgroundColor = [0.8 0.8 0.8];
            app.progress_simRepLabel_3.Position = [306 104 5 5];
            app.progress_simRepLabel_3.Text = '';

            % Create progress_simRepLabel_4
            app.progress_simRepLabel_4 = uilabel(app.BasicParameterTab);
            app.progress_simRepLabel_4.BackgroundColor = [0.8 0.8 0.8];
            app.progress_simRepLabel_4.Position = [323 104 5 5];
            app.progress_simRepLabel_4.Text = '';

            % Create progress_simRepLabel_5
            app.progress_simRepLabel_5 = uilabel(app.BasicParameterTab);
            app.progress_simRepLabel_5.BackgroundColor = [0.8 0.8 0.8];
            app.progress_simRepLabel_5.Position = [340 104 5 5];
            app.progress_simRepLabel_5.Text = '';

            % Create progress_simRepLabel_6
            app.progress_simRepLabel_6 = uilabel(app.BasicParameterTab);
            app.progress_simRepLabel_6.BackgroundColor = [0.8 0.8 0.8];
            app.progress_simRepLabel_6.Position = [356 104 5 5];
            app.progress_simRepLabel_6.Text = '';

            % Create progress_simRepLabel_7
            app.progress_simRepLabel_7 = uilabel(app.BasicParameterTab);
            app.progress_simRepLabel_7.BackgroundColor = [0.8 0.8 0.8];
            app.progress_simRepLabel_7.Position = [372 104 5 5];
            app.progress_simRepLabel_7.Text = '';

            % Create progress_simRepLabel_8
            app.progress_simRepLabel_8 = uilabel(app.BasicParameterTab);
            app.progress_simRepLabel_8.BackgroundColor = [0.8 0.8 0.8];
            app.progress_simRepLabel_8.Position = [388 104 5 5];
            app.progress_simRepLabel_8.Text = '';

            % Create progress_simRepLabel_9
            app.progress_simRepLabel_9 = uilabel(app.BasicParameterTab);
            app.progress_simRepLabel_9.BackgroundColor = [0.8 0.8 0.8];
            app.progress_simRepLabel_9.Position = [404 104 5 5];
            app.progress_simRepLabel_9.Text = '';

            % Create progress_simRepLabel_10
            app.progress_simRepLabel_10 = uilabel(app.BasicParameterTab);
            app.progress_simRepLabel_10.BackgroundColor = [0.8 0.8 0.8];
            app.progress_simRepLabel_10.Position = [420 104 5 5];
            app.progress_simRepLabel_10.Text = '';

            % Create FirstpulseampEditFieldLabel
            app.FirstpulseampEditFieldLabel = uilabel(app.BasicParameterTab);
            app.FirstpulseampEditFieldLabel.HorizontalAlignment = 'right';
            app.FirstpulseampEditFieldLabel.Position = [249 81 87 22];
            app.FirstpulseampEditFieldLabel.Text = 'First pulse amp';

            % Create amp_1stpulse_EditField
            app.amp_1stpulse_EditField = uieditfield(app.BasicParameterTab, 'numeric');
            app.amp_1stpulse_EditField.ValueChangedFcn = createCallbackFcn(app, @amp_1stpulse_EditFieldValueChanged, true);
            app.amp_1stpulse_EditField.Tooltip = {'Amplitude of the first pulse (0 - 1)'};
            app.amp_1stpulse_EditField.Position = [255 62 71 22];
            app.amp_1stpulse_EditField.Value = 0.5;

            % Create freq_1stpulse_EditFieldLabel
            app.freq_1stpulse_EditFieldLabel = uilabel(app.BasicParameterTab);
            app.freq_1stpulse_EditFieldLabel.Tooltip = {'Frequency 1st pulse'};
            app.freq_1stpulse_EditFieldLabel.Position = [344 83 34 22];
            app.freq_1stpulse_EditFieldLabel.Text = 'f [Hz]';

            % Create freq_1stpulse_HzEditField
            app.freq_1stpulse_HzEditField = uieditfield(app.BasicParameterTab, 'numeric');
            app.freq_1stpulse_HzEditField.ValueChangedFcn = createCallbackFcn(app, @freq_1stpulse_HzEditFieldValueChanged, true);
            app.freq_1stpulse_HzEditField.Tooltip = {'Frequency offset.'};
            app.freq_1stpulse_HzEditField.Position = [342 62 56 22];

            % Create PlotresultsButton
            app.PlotresultsButton = uibutton(app.BasicParameterTab, 'push');
            app.PlotresultsButton.ButtonPushedFcn = createCallbackFcn(app, @PlotresultsButtonPushed, true);
            app.PlotresultsButton.Tooltip = {'plot the simulation results again'};
            app.PlotresultsButton.Position = [323 6 91 22];
            app.PlotresultsButton.Text = 'Plot results';

            % Create SpoilerGradCheckBox
            app.SpoilerGradCheckBox = uicheckbox(app.BasicParameterTab);
            app.SpoilerGradCheckBox.ValueChangedFcn = createCallbackFcn(app, @SpoilerGradCheckBoxValueChanged, true);
            app.SpoilerGradCheckBox.Tooltip = {'"Spoil" Mxy Magnetization after N repetitions'};
            app.SpoilerGradCheckBox.Text = 'Spoiler Grad';
            app.SpoilerGradCheckBox.Position = [86 32 90 22];

            % Create TipbackCheckBox
            app.TipbackCheckBox = uicheckbox(app.BasicParameterTab);
            app.TipbackCheckBox.ValueChangedFcn = createCallbackFcn(app, @TipbackCheckBoxValueChanged, true);
            app.TipbackCheckBox.Tooltip = {'Tips the Magnetization after N repetitions back to the z-axis'};
            app.TipbackCheckBox.Text = 'Tipback';
            app.TipbackCheckBox.Position = [86 11 64 22];

            % Create first_repetition_time_msLabel
            app.first_repetition_time_msLabel = uilabel(app.BasicParameterTab);
            app.first_repetition_time_msLabel.Position = [10 37 67 22];
            app.first_repetition_time_msLabel.Text = '1st TR [ms]';

            % Create first_repetition_time_msEditField
            app.first_repetition_time_msEditField = uieditfield(app.BasicParameterTab, 'numeric');
            app.first_repetition_time_msEditField.Limits = [0.001 100];
            app.first_repetition_time_msEditField.ValueChangedFcn = createCallbackFcn(app, @first_repetition_time_msEditFieldValueChanged, true);
            app.first_repetition_time_msEditField.Tag = 'grad_strength_Edit';
            app.first_repetition_time_msEditField.Tooltip = {'Time between first and second pulse (peak to peak).'};
            app.first_repetition_time_msEditField.Position = [10 17 54 22];
            app.first_repetition_time_msEditField.Value = 10;

            % Create AdvParameterTab
            app.AdvParameterTab = uitab(app.SimulateRepetitionsTabGroup);
            app.AdvParameterTab.Title = 'Adv. parameters';

            % Create duration_refocus_grad
            app.duration_refocus_grad = uieditfield(app.AdvParameterTab, 'numeric');
            app.duration_refocus_grad.Limits = [0 100];
            app.duration_refocus_grad.Tag = 'grad_strength_Edit';
            app.duration_refocus_grad.Tooltip = {'duration of the slice refphasing gradient'};
            app.duration_refocus_grad.Position = [9 58 54 22];

            % Create TRGmsLabel
            app.TRGmsLabel = uilabel(app.AdvParameterTab);
            app.TRGmsLabel.Position = [9 78 57 22];
            app.TRGmsLabel.Text = 'TRG [ms]';

            % Create alphaEditFieldLabel
            app.alphaEditFieldLabel = uilabel(app.AdvParameterTab);
            app.alphaEditFieldLabel.HorizontalAlignment = 'right';
            app.alphaEditFieldLabel.Enable = 'off';
            app.alphaEditFieldLabel.Tooltip = {'Simulate 2 metabolites (change excitation frequency after X repetitions, N times)'};
            app.alphaEditFieldLabel.Position = [221 27 50 22];
            app.alphaEditFieldLabel.Text = 'alpha [°]';

            % Create alpha2ndPulseEditField
            app.alpha2ndPulseEditField = uieditfield(app.AdvParameterTab, 'numeric');
            app.alpha2ndPulseEditField.ValueChangedFcn = createCallbackFcn(app, @alpha2ndPulseEditFieldValueChanged, true);
            app.alpha2ndPulseEditField.Enable = 'off';
            app.alpha2ndPulseEditField.Tooltip = {'Flipangle of  the 2nd pulse'};
            app.alpha2ndPulseEditField.Position = [213 8 57 22];
            app.alpha2ndPulseEditField.Value = 10;

            % Create metabsCheckBox
            app.metabsCheckBox = uicheckbox(app.AdvParameterTab);
            app.metabsCheckBox.ValueChangedFcn = createCallbackFcn(app, @metabsCheckBoxValueChanged, true);
            app.metabsCheckBox.Tooltip = {'Simulate 2 metabolites (change excitation frequency after X repetitions, N times)'};
            app.metabsCheckBox.Text = '2 metabs.';
            app.metabsCheckBox.Position = [217 59 75 22];

            % Create AlterNTimesLabel
            app.AlterNTimesLabel = uilabel(app.AdvParameterTab);
            app.AlterNTimesLabel.HorizontalAlignment = 'right';
            app.AlterNTimesLabel.Enable = 'off';
            app.AlterNTimesLabel.Position = [286 32 57 28];
            app.AlterNTimesLabel.Text = {'Alternate '; 'N times'};

            % Create AlterNTimesEditField
            app.AlterNTimesEditField = uieditfield(app.AdvParameterTab, 'numeric');
            app.AlterNTimesEditField.ValueChangedFcn = createCallbackFcn(app, @AlterNTimesEditFieldValueChanged, true);
            app.AlterNTimesEditField.Enable = 'off';
            app.AlterNTimesEditField.Tooltip = {'When to change between metabolites'};
            app.AlterNTimesEditField.Position = [286 8 57 22];
            app.AlterNTimesEditField.Value = 1;

            % Create fHzEditFieldLabel
            app.fHzEditFieldLabel = uilabel(app.AdvParameterTab);
            app.fHzEditFieldLabel.HorizontalAlignment = 'right';
            app.fHzEditFieldLabel.Enable = 'off';
            app.fHzEditFieldLabel.Tooltip = {'Frequency 2nd pulse'};
            app.fHzEditFieldLabel.Position = [380 29 34 22];
            app.fHzEditFieldLabel.Text = 'f [Hz]';

            % Create freq_2ndpulse_HzEditField
            app.freq_2ndpulse_HzEditField = uieditfield(app.AdvParameterTab, 'numeric');
            app.freq_2ndpulse_HzEditField.ValueChangedFcn = createCallbackFcn(app, @freq_2ndpulse_HzEditFieldValueChanged, true);
            app.freq_2ndpulse_HzEditField.Enable = 'off';
            app.freq_2ndpulse_HzEditField.Tooltip = {'Frequency 2nd pulse'};
            app.freq_2ndpulse_HzEditField.Position = [354 8 58 22];
            app.freq_2ndpulse_HzEditField.Value = -900;

            % Create altphaseEditFieldLabel
            app.altphaseEditFieldLabel = uilabel(app.AdvParameterTab);
            app.altphaseEditFieldLabel.HorizontalAlignment = 'right';
            app.altphaseEditFieldLabel.Position = [9 27 58 22];
            app.altphaseEditFieldLabel.Text = 'alt. phase';

            % Create altphaseEditField
            app.altphaseEditField = uieditfield(app.AdvParameterTab, 'numeric');
            app.altphaseEditField.ValueChangedFcn = createCallbackFcn(app, @altphaseEditFieldValueChanged, true);
            app.altphaseEditField.Tooltip = {'changing pulse phase after each "image". From 0 to set phase and back.'};
            app.altphaseEditField.Position = [9 9 52 22];

            % Create homoB1CheckBox
            app.homoB1CheckBox = uicheckbox(app.AdvParameterTab);
            app.homoB1CheckBox.Enable = 'off';
            app.homoB1CheckBox.Visible = 'off';
            app.homoB1CheckBox.Tooltip = {'Simulate repeated Excitation.'; 'Assume homogeneous B1 field.'};
            app.homoB1CheckBox.Text = 'homo. B1';
            app.homoB1CheckBox.Position = [338 84 74 22];
            app.homoB1CheckBox.Value = true;

            % Create highresoutCheckBox
            app.highresoutCheckBox = uicheckbox(app.AdvParameterTab);
            app.highresoutCheckBox.ValueChangedFcn = createCallbackFcn(app, @highresoutCheckBoxValueChanged, true);
            app.highresoutCheckBox.Tooltip = {'Save every simulated timepoint (N [points] per repetition). This dataset might get quite big really fast (N reps * N points * N freq * N y).'};
            app.highresoutCheckBox.Text = 'high res out';
            app.highresoutCheckBox.Position = [332 65 85 22];

            % Create alpha2prepCheckBox
            app.alpha2prepCheckBox = uicheckbox(app.AdvParameterTab);
            app.alpha2prepCheckBox.ValueChangedFcn = createCallbackFcn(app, @alpha2prepCheckBoxValueChanged, true);
            app.alpha2prepCheckBox.Text = 'alpha/2 prep';
            app.alpha2prepCheckBox.Position = [217 81 89 22];
            app.alpha2prepCheckBox.Value = true;

            % Create phasetbEditFieldLabel
            app.phasetbEditFieldLabel = uilabel(app.AdvParameterTab);
            app.phasetbEditFieldLabel.HorizontalAlignment = 'right';
            app.phasetbEditFieldLabel.Position = [108 77 52 22];
            app.phasetbEditFieldLabel.Text = 'phase tb';

            % Create phasetbEditField
            app.phasetbEditField = uieditfield(app.AdvParameterTab, 'numeric');
            app.phasetbEditField.ValueChangedFcn = createCallbackFcn(app, @phasetbEditFieldValueChanged, true);
            app.phasetbEditField.Tooltip = {'phase of the tipback pulse (in degree, difference to normal phase of pulse)'};
            app.phasetbEditField.Position = [171 77 33 22];

            % Create SimrangeTab
            app.SimrangeTab = uitab(app.SimulateRepetitionsTabGroup);
            app.SimrangeTab.Title = 'Sim. range';

            % Create T1EditFieldLabel
            app.T1EditFieldLabel = uilabel(app.SimrangeTab);
            app.T1EditFieldLabel.HorizontalAlignment = 'right';
            app.T1EditFieldLabel.Position = [9 84 25 22];
            app.T1EditFieldLabel.Text = 'T1';

            % Create T1Range_1EditField
            app.T1Range_1EditField = uieditfield(app.SimrangeTab, 'numeric');
            app.T1Range_1EditField.Position = [49 84 49 22];
            app.T1Range_1EditField.Value = 0.2;

            % Create T1Range_2EditField
            app.T1Range_2EditField = uieditfield(app.SimrangeTab, 'numeric');
            app.T1Range_2EditField.Position = [103 84 49 22];
            app.T1Range_2EditField.Value = 30;

            % Create T1Range_NEditField
            app.T1Range_NEditField = uieditfield(app.SimrangeTab, 'numeric');
            app.T1Range_NEditField.ValueChangedFcn = createCallbackFcn(app, @T1Range_NEditFieldValueChanged, true);
            app.T1Range_NEditField.Position = [159 84 49 22];
            app.T1Range_NEditField.Value = 2;

            % Create T2Range_2EditField
            app.T2Range_2EditField = uieditfield(app.SimrangeTab, 'numeric');
            app.T2Range_2EditField.Position = [103 56 49 22];
            app.T2Range_2EditField.Value = 25;

            % Create T2Range_NEditField
            app.T2Range_NEditField = uieditfield(app.SimrangeTab, 'numeric');
            app.T2Range_NEditField.ValueChangedFcn = createCallbackFcn(app, @T2Range_NEditFieldValueChanged, true);
            app.T2Range_NEditField.Position = [159 56 49 22];
            app.T2Range_NEditField.Value = 2;

            % Create T2Label
            app.T2Label = uilabel(app.SimrangeTab);
            app.T2Label.HorizontalAlignment = 'right';
            app.T2Label.Position = [9 57 25 22];
            app.T2Label.Text = 'T2';

            % Create T2Range_1EditField
            app.T2Range_1EditField = uieditfield(app.SimrangeTab, 'numeric');
            app.T2Range_1EditField.Position = [49 57 49 22];
            app.T2Range_1EditField.Value = 0.1;

            % Create T2Label_2
            app.T2Label_2 = uilabel(app.SimrangeTab);
            app.T2Label_2.HorizontalAlignment = 'right';
            app.T2Label_2.Position = [8 30 25 22];
            app.T2Label_2.Text = 'TR';

            % Create TRRange_1EditField
            app.TRRange_1EditField = uieditfield(app.SimrangeTab, 'numeric');
            app.TRRange_1EditField.Position = [48 30 49 22];
            app.TRRange_1EditField.Value = 2;

            % Create TRRange_2EditField
            app.TRRange_2EditField = uieditfield(app.SimrangeTab, 'numeric');
            app.TRRange_2EditField.Position = [102 29 49 22];
            app.TRRange_2EditField.Value = 10;

            % Create TRRange_NEditField
            app.TRRange_NEditField = uieditfield(app.SimrangeTab, 'numeric');
            app.TRRange_NEditField.ValueChangedFcn = createCallbackFcn(app, @TRRange_NEditFieldValueChanged, true);
            app.TRRange_NEditField.Position = [158 29 49 22];
            app.TRRange_NEditField.Value = 2;

            % Create FARange_NEditField
            app.FARange_NEditField = uieditfield(app.SimrangeTab, 'numeric');
            app.FARange_NEditField.ValueChangedFcn = createCallbackFcn(app, @FARange_NEditFieldValueChanged, true);
            app.FARange_NEditField.Position = [159 1 49 22];
            app.FARange_NEditField.Value = 2;

            % Create FARange_2EditField
            app.FARange_2EditField = uieditfield(app.SimrangeTab, 'numeric');
            app.FARange_2EditField.Position = [103 1 49 22];
            app.FARange_2EditField.Value = 90;

            % Create FARange_1EditField
            app.FARange_1EditField = uieditfield(app.SimrangeTab, 'numeric');
            app.FARange_1EditField.Position = [49 2 49 22];
            app.FARange_1EditField.Value = 2;

            % Create T2Label_3
            app.T2Label_3 = uilabel(app.SimrangeTab);
            app.T2Label_3.HorizontalAlignment = 'right';
            app.T2Label_3.Position = [9 2 25 22];
            app.T2Label_3.Text = 'FA';

            % Create SimRangeButton
            app.SimRangeButton = uibutton(app.SimrangeTab, 'push');
            app.SimRangeButton.ButtonPushedFcn = createCallbackFcn(app, @SimRangeButtonPushed, true);
            app.SimRangeButton.Position = [331 1 73 22];
            app.SimRangeButton.Text = 'Sim Range!';

            % Create simrangeCheckBox
            app.simrangeCheckBox = uicheckbox(app.SimrangeTab);
            app.simrangeCheckBox.ValueChangedFcn = createCallbackFcn(app, @simrangeCheckBoxValueChanged, true);
            app.simrangeCheckBox.Text = 'sim. range';
            app.simrangeCheckBox.Position = [244 2 78 22];

            % Create TrfRange_NEditField
            app.TrfRange_NEditField = uieditfield(app.SimrangeTab, 'numeric');
            app.TrfRange_NEditField.ValueChangedFcn = createCallbackFcn(app, @TrfRange_NEditFieldValueChanged, true);
            app.TrfRange_NEditField.Position = [366 84 49 22];
            app.TrfRange_NEditField.Value = 2;

            % Create TrfRange_2EditField
            app.TrfRange_2EditField = uieditfield(app.SimrangeTab, 'numeric');
            app.TrfRange_2EditField.Position = [310 84 49 22];
            app.TrfRange_2EditField.Value = 30;

            % Create TrfRange_1EditField
            app.TrfRange_1EditField = uieditfield(app.SimrangeTab, 'numeric');
            app.TrfRange_1EditField.Position = [256 84 49 22];
            app.TrfRange_1EditField.Value = 0.1;

            % Create TrfEditFieldLabel
            app.TrfEditFieldLabel = uilabel(app.SimrangeTab);
            app.TrfEditFieldLabel.HorizontalAlignment = 'right';
            app.TrfEditFieldLabel.Position = [216 84 25 22];
            app.TrfEditFieldLabel.Text = 'Trf';

            % Create plotresCheckBox
            app.plotresCheckBox = uicheckbox(app.SimrangeTab);
            app.plotresCheckBox.Text = 'plot res.';
            app.plotresCheckBox.Position = [244 20 64 22];

            % Create SimulaterepeatedExcitationFISPbSSFPLabel
            app.SimulaterepeatedExcitationFISPbSSFPLabel = uilabel(app.LeftPanel);
            app.SimulaterepeatedExcitationFISPbSSFPLabel.FontColor = [0 0.7098 0];
            app.SimulaterepeatedExcitationFISPbSSFPLabel.Position = [16 318 240 22];
            app.SimulaterepeatedExcitationFISPbSSFPLabel.Text = 'Simulate repeated Excitation (FISP/bSSFP)';

            % Create SimulateRepsButton
            app.SimulateRepsButton = uibutton(app.LeftPanel, 'push');
            app.SimulateRepsButton.ButtonPushedFcn = createCallbackFcn(app, @SimulateRepsButtonPushed, true);
            app.SimulateRepsButton.Tag = 'simulate_Button';
            app.SimulateRepsButton.Position = [330 322 92 19];
            app.SimulateRepsButton.Text = 'Simulate Reps!';

            % Create progress_sim_repLabel
            app.progress_sim_repLabel = uilabel(app.LeftPanel);
            app.progress_sim_repLabel.Position = [260 321 63 22];
            app.progress_sim_repLabel.Text = '';

            % Create startparpoolButton
            app.startparpoolButton = uibutton(app.LeftPanel, 'push');
            app.startparpoolButton.ButtonPushedFcn = createCallbackFcn(app, @startparpoolButtonPushed, true);
            app.startparpoolButton.Tooltip = {'Starting the parallel pool (CPU cores will run in parallel). Parallel pool will be started in the repeated excitations anyways.'};
            app.startparpoolButton.Position = [337 714 100 22];
            app.startparpoolButton.Text = 'start par. pool';

            % Create reset_B1map_Button
            app.reset_B1map_Button = uibutton(app.LeftPanel, 'push');
            app.reset_B1map_Button.ButtonPushedFcn = createCallbackFcn(app, @reset_B1map_ButtonPushed, true);
            app.reset_B1map_Button.FontColor = [1 0 0];
            app.reset_B1map_Button.Tooltip = {'deletes the B1 profile'};
            app.reset_B1map_Button.Position = [188 796 66 22];
            app.reset_B1map_Button.Text = 'reset';

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;
            app.RightPanel.Scrollable = 'on';

            % Create UIAxes_b1profile
            app.UIAxes_b1profile = uiaxes(app.RightPanel);
            title(app.UIAxes_b1profile, 'B1 profile')
            xlabel(app.UIAxes_b1profile, 'y [cm]')
            ylabel(app.UIAxes_b1profile, 'Y')
            zlabel(app.UIAxes_b1profile, 'Z')
            app.UIAxes_b1profile.XGrid = 'on';
            app.UIAxes_b1profile.YGrid = 'on';
            app.UIAxes_b1profile.ButtonDownFcn = createCallbackFcn(app, @UIAxes_b1profileButtonDown, true);
            app.UIAxes_b1profile.Position = [454 719 288 182];

            % Create UIAxes_rfpulse
            app.UIAxes_rfpulse = uiaxes(app.RightPanel);
            title(app.UIAxes_rfpulse, {'RF Pulse'; ''})
            xlabel(app.UIAxes_rfpulse, {'time [ms]'; ''})
            ylabel(app.UIAxes_rfpulse, 'Y')
            zlabel(app.UIAxes_rfpulse, 'Z')
            app.UIAxes_rfpulse.XGrid = 'on';
            app.UIAxes_rfpulse.YGrid = 'on';
            app.UIAxes_rfpulse.Position = [23 719 284 182];

            % Create UIAxes_mag_vs_freq
            app.UIAxes_mag_vs_freq = uiaxes(app.RightPanel);
            title(app.UIAxes_mag_vs_freq, {'Magnetization vs. distance'; ''})
            xlabel(app.UIAxes_mag_vs_freq, {'y [cm]'; ''})
            ylabel(app.UIAxes_mag_vs_freq, 'Y')
            zlabel(app.UIAxes_mag_vs_freq, 'Z')
            app.UIAxes_mag_vs_freq.XGrid = 'on';
            app.UIAxes_mag_vs_freq.YGrid = 'on';
            app.UIAxes_mag_vs_freq.Position = [23 534 284 182];

            % Create UIAxes_anat_image
            app.UIAxes_anat_image = uiaxes(app.RightPanel);
            title(app.UIAxes_anat_image, 'Anat Image')
            xlabel(app.UIAxes_anat_image, 'X')
            ylabel(app.UIAxes_anat_image, 'Y')
            zlabel(app.UIAxes_anat_image, 'Z')
            app.UIAxes_anat_image.XGrid = 'on';
            app.UIAxes_anat_image.YGrid = 'on';
            app.UIAxes_anat_image.Position = [454 506 288 210];

            % Create freqSliderLabel
            app.freqSliderLabel = uilabel(app.RightPanel);
            app.freqSliderLabel.HorizontalAlignment = 'right';
            app.freqSliderLabel.Position = [15 240 26 22];
            app.freqSliderLabel.Text = 'freq';

            % Create freqSlider
            app.freqSlider = uislider(app.RightPanel);
            app.freqSlider.ValueChangedFcn = createCallbackFcn(app, @freqSliderValueChanged, true);
            app.freqSlider.ValueChangingFcn = createCallbackFcn(app, @freqSliderValueChanging, true);
            app.freqSlider.Position = [58 249 194 3];

            % Create Show2FreqsCheckBox
            app.Show2FreqsCheckBox = uicheckbox(app.RightPanel);
            app.Show2FreqsCheckBox.ValueChangedFcn = createCallbackFcn(app, @Show2FreqsCheckBoxValueChanged, true);
            app.Show2FreqsCheckBox.Text = 'Show 2 Freqs';
            app.Show2FreqsCheckBox.Position = [149 381 96 22];

            % Create ShowExcitationProfileCheckBox
            app.ShowExcitationProfileCheckBox = uicheckbox(app.RightPanel);
            app.ShowExcitationProfileCheckBox.ValueChangedFcn = createCallbackFcn(app, @ShowExcitationProfileCheckBoxValueChanged, true);
            app.ShowExcitationProfileCheckBox.Text = ' Show Excitation Profile';
            app.ShowExcitationProfileCheckBox.Position = [339 506 148 22];

            % Create SliceSpinnerLabel
            app.SliceSpinnerLabel = uilabel(app.RightPanel);
            app.SliceSpinnerLabel.HorizontalAlignment = 'right';
            app.SliceSpinnerLabel.Position = [342 474 35 22];
            app.SliceSpinnerLabel.Text = 'Slice ';

            % Create SliceSpinner
            app.SliceSpinner = uispinner(app.RightPanel);
            app.SliceSpinner.Limits = [1 Inf];
            app.SliceSpinner.ValueChangedFcn = createCallbackFcn(app, @SliceSpinnerValueChanged, true);
            app.SliceSpinner.Position = [410 474 66 22];
            app.SliceSpinner.Value = 1;

            % Create MetabolitesPanel
            app.MetabolitesPanel = uipanel(app.RightPanel);
            app.MetabolitesPanel.Title = 'Metabolites';
            app.MetabolitesPanel.Position = [306 213 208 160];

            % Create Freq1EditField_2Label
            app.Freq1EditField_2Label = uilabel(app.MetabolitesPanel);
            app.Freq1EditField_2Label.HorizontalAlignment = 'right';
            app.Freq1EditField_2Label.Position = [8 108 40 22];
            app.Freq1EditField_2Label.Text = 'Freq 1';

            % Create Freq1_ppmEditField
            app.Freq1_ppmEditField = uieditfield(app.MetabolitesPanel, 'numeric');
            app.Freq1_ppmEditField.ValueChangedFcn = createCallbackFcn(app, @Freq1_ppmEditFieldValueChanged, true);
            app.Freq1_ppmEditField.Position = [76 108 94 22];
            app.Freq1_ppmEditField.Value = 175;

            % Create Freq2EditField_2Label
            app.Freq2EditField_2Label = uilabel(app.MetabolitesPanel);
            app.Freq2EditField_2Label.HorizontalAlignment = 'right';
            app.Freq2EditField_2Label.Position = [8 76 40 22];
            app.Freq2EditField_2Label.Text = 'Freq 2';

            % Create Freq2_ppmEditField
            app.Freq2_ppmEditField = uieditfield(app.MetabolitesPanel, 'numeric');
            app.Freq2_ppmEditField.ValueChangedFcn = createCallbackFcn(app, @Freq2_ppmEditFieldValueChanged, true);
            app.Freq2_ppmEditField.Position = [76 76 95 22];
            app.Freq2_ppmEditField.Value = 183.5;

            % Create ppmLabel
            app.ppmLabel = uilabel(app.MetabolitesPanel);
            app.ppmLabel.Position = [173 108 29 22];
            app.ppmLabel.Text = 'ppm';

            % Create ppmLabel_2
            app.ppmLabel_2 = uilabel(app.MetabolitesPanel);
            app.ppmLabel_2.Position = [174 76 28 22];
            app.ppmLabel_2.Text = 'ppm';

            % Create BasicFreqEditFieldLabel
            app.BasicFreqEditFieldLabel = uilabel(app.MetabolitesPanel);
            app.BasicFreqEditFieldLabel.HorizontalAlignment = 'right';
            app.BasicFreqEditFieldLabel.Position = [4 45 63 22];
            app.BasicFreqEditFieldLabel.Text = 'Basic Freq';

            % Create BasicFreqEditField
            app.BasicFreqEditField = uieditfield(app.MetabolitesPanel, 'numeric');
            app.BasicFreqEditField.ValueChangedFcn = createCallbackFcn(app, @BasicFreqEditFieldValueChanged, true);
            app.BasicFreqEditField.Position = [76 45 94 22];
            app.BasicFreqEditField.Value = 183.5;

            % Create ppmLabel_3
            app.ppmLabel_3 = uilabel(app.MetabolitesPanel);
            app.ppmLabel_3.Position = [174 45 29 22];
            app.ppmLabel_3.Text = 'ppm';

            % Create PulseFreqEditFieldLabel
            app.PulseFreqEditFieldLabel = uilabel(app.MetabolitesPanel);
            app.PulseFreqEditFieldLabel.HorizontalAlignment = 'right';
            app.PulseFreqEditFieldLabel.Position = [4 12 64 22];
            app.PulseFreqEditFieldLabel.Text = 'Pulse Freq';

            % Create HzLabel_3
            app.HzLabel_3 = uilabel(app.MetabolitesPanel);
            app.HzLabel_3.Position = [175 12 25 22];
            app.HzLabel_3.Text = 'Hz';

            % Create PulseFreqEditField
            app.PulseFreqEditField = uieditfield(app.MetabolitesPanel, 'numeric');
            app.PulseFreqEditField.ValueChangedFcn = createCallbackFcn(app, @PulseFreqEditFieldValueChanged, true);
            app.PulseFreqEditField.Position = [76 12 94 22];

            % Create FrequenciesPanel
            app.FrequenciesPanel = uipanel(app.RightPanel);
            app.FrequenciesPanel.Title = 'Frequencies';
            app.FrequenciesPanel.Position = [22 279 226 94];

            % Create Freq1EditFieldLabel
            app.Freq1EditFieldLabel = uilabel(app.FrequenciesPanel);
            app.Freq1EditFieldLabel.HorizontalAlignment = 'right';
            app.Freq1EditFieldLabel.Position = [5 42 52 22];
            app.Freq1EditFieldLabel.Text = 'Freq 1';

            % Create Freq1EditField
            app.Freq1EditField = uieditfield(app.FrequenciesPanel, 'numeric');
            app.Freq1EditField.ValueChangedFcn = createCallbackFcn(app, @Freq1EditFieldValueChanged, true);
            app.Freq1EditField.Position = [66 40 112 22];

            % Create Freq2EditFieldLabel
            app.Freq2EditFieldLabel = uilabel(app.FrequenciesPanel);
            app.Freq2EditFieldLabel.HorizontalAlignment = 'right';
            app.Freq2EditFieldLabel.Position = [-10 10 67 22];
            app.Freq2EditFieldLabel.Text = 'Freq 2';

            % Create Freq2EditField
            app.Freq2EditField = uieditfield(app.FrequenciesPanel, 'numeric');
            app.Freq2EditField.ValueChangedFcn = createCallbackFcn(app, @Freq2EditFieldValueChanged, true);
            app.Freq2EditField.Position = [66 10 112 22];

            % Create HzLabel
            app.HzLabel = uilabel(app.FrequenciesPanel);
            app.HzLabel.Position = [185 42 25 22];
            app.HzLabel.Text = 'Hz';

            % Create HzLabel_2
            app.HzLabel_2 = uilabel(app.FrequenciesPanel);
            app.HzLabel_2.Position = [186 10 28 22];
            app.HzLabel_2.Text = 'Hz';

            % Create UseFrequenciesCheckBox
            app.UseFrequenciesCheckBox = uicheckbox(app.RightPanel);
            app.UseFrequenciesCheckBox.ValueChangedFcn = createCallbackFcn(app, @UseFrequenciesCheckBoxValueChanged, true);
            app.UseFrequenciesCheckBox.Text = 'Use Frequencies';
            app.UseFrequenciesCheckBox.Position = [26 381 113 22];
            app.UseFrequenciesCheckBox.Value = true;

            % Create inhomoB1ListBoxLabel
            app.inhomoB1ListBoxLabel = uilabel(app.RightPanel);
            app.inhomoB1ListBoxLabel.HorizontalAlignment = 'right';
            app.inhomoB1ListBoxLabel.Position = [66 523 66 22];
            app.inhomoB1ListBoxLabel.Text = 'inhomo. B1';

            % Create inhomoB1ListBox
            app.inhomoB1ListBox = uilistbox(app.RightPanel);
            app.inhomoB1ListBox.Items = {'Mx', 'My', '|Mxy|', 'Mz', 'phi', 'none'};
            app.inhomoB1ListBox.Multiselect = 'on';
            app.inhomoB1ListBox.ValueChangedFcn = createCallbackFcn(app, @inhomoB1ListBoxValueChanged, true);
            app.inhomoB1ListBox.Position = [62 414 91 110];
            app.inhomoB1ListBox.Value = {'none'};

            % Create homoB1ListBoxLabel
            app.homoB1ListBoxLabel = uilabel(app.RightPanel);
            app.homoB1ListBoxLabel.HorizontalAlignment = 'right';
            app.homoB1ListBoxLabel.FontColor = [0 0 1];
            app.homoB1ListBoxLabel.Position = [163 525 57 22];
            app.homoB1ListBoxLabel.Text = 'homo. B1';

            % Create homoB1ListBox
            app.homoB1ListBox = uilistbox(app.RightPanel);
            app.homoB1ListBox.Items = {'Mx', 'My', '|Mxy|', 'Mz', 'phi', 'none'};
            app.homoB1ListBox.Multiselect = 'on';
            app.homoB1ListBox.ValueChangedFcn = createCallbackFcn(app, @homoB1ListBoxValueChanged, true);
            app.homoB1ListBox.FontColor = [0 0 1];
            app.homoB1ListBox.Position = [159 414 100 110];
            app.homoB1ListBox.Value = {'|Mxy|'};

            % Create refpow_map_caxis
            app.refpow_map_caxis = uislider(app.RightPanel);
            app.refpow_map_caxis.Limits = [1e-06 100];
            app.refpow_map_caxis.Orientation = 'vertical';
            app.refpow_map_caxis.ValueChangedFcn = createCallbackFcn(app, @refpow_map_caxisValueChanged, true);
            app.refpow_map_caxis.Visible = 'off';
            app.refpow_map_caxis.FontSize = 10;
            app.refpow_map_caxis.Position = [772 736 3 148];
            app.refpow_map_caxis.Value = 1e-06;

            % Create bssfpsignalCheckBox
            app.bssfpsignalCheckBox = uicheckbox(app.RightPanel);
            app.bssfpsignalCheckBox.ValueChangedFcn = createCallbackFcn(app, @bssfpsignalCheckBoxValueChanged, true);
            app.bssfpsignalCheckBox.Tooltip = {'if activated, bssfp signal is plotted'};
            app.bssfpsignalCheckBox.Text = 'bssfp signal';
            app.bssfpsignalCheckBox.Position = [27 174 86 22];

            % Create responseprofileCheckBox
            app.responseprofileCheckBox = uicheckbox(app.RightPanel);
            app.responseprofileCheckBox.Text = 'response profile';
            app.responseprofileCheckBox.Position = [132 174 107 22];

            % Create repSliderLabel
            app.repSliderLabel = uilabel(app.RightPanel);
            app.repSliderLabel.HorizontalAlignment = 'right';
            app.repSliderLabel.Position = [21 140 25 22];
            app.repSliderLabel.Text = 'rep';

            % Create repSlider
            app.repSlider = uislider(app.RightPanel);
            app.repSlider.Limits = [1 100];
            app.repSlider.ValueChangedFcn = createCallbackFcn(app, @repSliderValueChanged, true);
            app.repSlider.Position = [63 149 194 3];
            app.repSlider.Value = 1;

            % Show the figure after all components are created
            app.calc_rf_pulse_toolUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = simulate_rf_pulse_b1_tool_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.calc_rf_pulse_toolUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.calc_rf_pulse_toolUIFigure)
        end
    end
end