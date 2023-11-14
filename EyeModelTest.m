%
% A computational method to distinguish between axial and refractive causes of refractive error.
% LNT, JR, MJ
%


    %% Figure 1 - Axial
    % 
    clc; 
    clear variables, close all
    
    
    fh = figure(1);
    set(fh, 'Units', 'normalized');
    set(fh, 'Position', [0.1 0.1 0.8 0.8]);
    set(fh, 'Name', 'Axial');
    
    % increase ACD makes eye hyperopic
    subplot(2,2,1);
        loopSpan = -2:0.1:2;
        
        refractiveErrors = zeros(length(loopSpan),1);
        dioptricLengths  = zeros(length(loopSpan),1);
        equivalentPowers = zeros(length(loopSpan),1);
        
        eye = EyeModel;
        anteriorChamberDepth = eye.getAnteriorChamberDepth();

        for i = 1:length(loopSpan)
                
            eye = eye.setAnteriorChamberDepth(anteriorChamberDepth + loopSpan(i));

            refractiveErrors(i) = eye.getRefractiveError();
            dioptricLengths(i)  = eye.getDioptricLength();
            equivalentPowers(i) = eye.getEquivalentPower();
            
        end
        
        
        yyaxis left;
        plot(anteriorChamberDepth + loopSpan, refractiveErrors);
        ylabel('Refractive Error [D]');
        xlabel('Anterior Chamber Depth [mm]');
        hold on
        yyaxis right;
        ratio = refractiveErrors./equivalentPowers;
        plot(anteriorChamberDepth + loopSpan, ratio); 
        plot([anteriorChamberDepth anteriorChamberDepth], [min(refractiveErrors) max(refractiveErrors)], 'k'); hold off;
        ylabel('Dioptric Length / Equiv. Power'); ylim([min(ratio) max(ratio)]);
        legend('Conventional', 'Normalized', 'Emmetropic reference', 'Location', 'Best');
        title('Effect of Varying Anterior Chamber Depth');
      
        grid on; grid minor;
        
    subplot(2,2,2);
        loopSpan = -2:0.1:2;

        refractiveErrors = zeros(length(loopSpan),1);
        dioptricLengths  = zeros(length(loopSpan),1);
        equivalentPowers = zeros(length(loopSpan),1);

        eye = EyeModel;
        lensThickness = eye.getLensThickness();

        for i = 1:length(loopSpan)

            eye = eye.setLensThickness(lensThickness + loopSpan(i));

            refractiveErrors(i) = eye.getRefractiveError();
            dioptricLengths(i)  = eye.getDioptricLength();
            equivalentPowers(i) = eye.getEquivalentPower();

        end

        title('Lens Thickness');
        xlabel('Lens Thickness [mm]');

        yyaxis left;
        plot(lensThickness + loopSpan, refractiveErrors);
        ylabel('Refractive Error [D]');

        yyaxis right;
        ratio = dioptricLengths./equivalentPowers;
        plot(lensThickness + loopSpan, ratio); hold on;
        plot([lensThickness lensThickness], [min(ratio) max(ratio)], 'k'); hold off;
        ylabel('Dioptric Length / Equiv. Power'); ylim([min(ratio) max(ratio)]);

        grid on; grid minor;
        
    subplot(2,2,3);
        loopSpan = -2:0.1:2;
        
        refractiveErrors = zeros(length(loopSpan),1);
        dioptricLengths  = zeros(length(loopSpan),1);
        equivalentPowers = zeros(length(loopSpan),1);
        
        eye = EyeModel;
        vitrealChamberDepth = eye.getVitrealChamberDepth();

        for i = 1:length(loopSpan)
                
            eye = eye.setVitrealChamberDepth(vitrealChamberDepth + loopSpan(i));

            refractiveErrors(i) = eye.getRefractiveError();
            dioptricLengths(i)  = eye.getDioptricLength();
            equivalentPowers(i) = eye.getEquivalentPower();
            
        end
        
        title('Vitreal Chamber Depth');
        xlabel('Vitreal Chamber Depth [mm]');
        
        yyaxis left;
        plot(vitrealChamberDepth + loopSpan, refractiveErrors);
        ylabel('Refractive Error [D]');

        yyaxis right;
        ratio = dioptricLengths./equivalentPowers;
        plot(vitrealChamberDepth + loopSpan, ratio); hold on;
        plot([vitrealChamberDepth vitrealChamberDepth], [min(ratio) max(ratio)], 'k'); hold off;
        ylabel('Dioptric Length / Equiv. Power'); ylim([min(ratio) max(ratio)]);
        
        grid on; grid minor;
        
    subplot(2,2,4);
        loopSpan = 0.8:0.01:1.2;
        
        refractiveErrors = zeros(length(loopSpan),1);
        dioptricLengths  = zeros(length(loopSpan),1);
        equivalentPowers = zeros(length(loopSpan),1);
        
        for i = 1:length(loopSpan)
                
            eye = EyeModel;
            eye = eye.setAxialElongation(loopSpan(i));

            refractiveErrors(i) = eye.getRefractiveError();
            dioptricLengths(i)  = eye.getDioptricLength();
            equivalentPowers(i) = eye.getEquivalentPower();
            
        end
        
        title('Axial Elongation');
        xlabel('Axial Elongation [x]');
        
        yyaxis left;
        plot(loopSpan, refractiveErrors);
        ylabel('Refractive Error [D]');

        yyaxis right;
        ratio = dioptricLengths./equivalentPowers;
        plot(loopSpan, ratio); hold on;
        plot([1.0 1.0], [min(ratio) max(ratio)], 'k'); hold off;
        ylabel('Dioptric Length / Equiv. Power'); ylim([min(ratio) max(ratio)]);
        
        grid on; grid minor;
        
        
    %% Figure 2 - Refractive    
        
    fh = figure(2);
    set(fh, 'Units', 'normalized');
    set(fh, 'Position', [0.1 0.1 0.8 0.8]);
    set(fh, 'Name', 'Refractive');
    
    subplot(2,2,1);
        loopSpan = -0.2:0.01:0.2;
        
        refractiveErrors = zeros(length(loopSpan),1);
        dioptricLengths  = zeros(length(loopSpan),1);
        equivalentPowers = zeros(length(loopSpan),1);
        
        eye = EyeModel;
        aqueousRefractiveIndex = eye.getAqueousRefractiveIndex();

        for i = 1:length(loopSpan)
                
            eye = eye.setAqueousRefractiveIndex(aqueousRefractiveIndex + loopSpan(i));

            refractiveErrors(i) = eye.getRefractiveError();
            dioptricLengths(i)  = eye.getDioptricLength();
            equivalentPowers(i) = eye.getEquivalentPower();
            
        end
        
        title('Aqueous Refractive Index');
        xlabel('Aqueous Refractive Index');
        
        yyaxis left;
        plot(aqueousRefractiveIndex + loopSpan, refractiveErrors);
        ylabel('Refractive Error [D]');

        yyaxis right;
        ratio = dioptricLengths./equivalentPowers;
        plot(aqueousRefractiveIndex + loopSpan, ratio); hold on;
        plot([aqueousRefractiveIndex aqueousRefractiveIndex], [min(ratio) max(ratio)], 'k'); hold off;
        ylabel('Dioptric Length / Equiv. Power'); ylim([min(ratio) max(ratio)]);
        legend('Refractive Error [D]', 'Dioptric Length / Equiv. Power', 'Emmetropic Eye', 'Location', 'Best');
        
        grid on; grid minor;

    subplot(2,2,2);
        loopSpan = -0.2:0.01:0.2;
        
        refractiveErrors = zeros(length(loopSpan),1);
        dioptricLengths  = zeros(length(loopSpan),1);
        equivalentPowers = zeros(length(loopSpan),1);
        
        eye = EyeModel;
        crystallineRefractiveIndex = eye.getCrystallineRefractiveIndex();

        for i = 1:length(loopSpan)
                
            eye = eye.setCrystallineRefractiveIndex(crystallineRefractiveIndex + loopSpan(i));

            refractiveErrors(i) = eye.getRefractiveError();
            dioptricLengths(i)  = eye.getDioptricLength();
            equivalentPowers(i) = eye.getEquivalentPower();
            
        end
        
        title('Crystalline Refractive Index');
        xlabel('Crystalline Refractive Index');
        
        yyaxis left;
        plot(crystallineRefractiveIndex + loopSpan, refractiveErrors);
        ylabel('Refractive Error [D]');

        yyaxis right;
        ratio = dioptricLengths./equivalentPowers;
        plot(crystallineRefractiveIndex + loopSpan, ratio); hold on;
        plot([crystallineRefractiveIndex crystallineRefractiveIndex], [min(ratio) max(ratio)], 'k'); hold off;
        ylabel('Dioptric Length / Equiv. Power'); ylim([min(ratio) max(ratio)]);
        
        grid on; grid minor;
    
    subplot(2,2,3);
        loopSpan = -0.2:0.01:0.2;
        
        refractiveErrors = zeros(length(loopSpan),1);
        dioptricLengths  = zeros(length(loopSpan),1);
        equivalentPowers = zeros(length(loopSpan),1);
        
        eye = EyeModel;
        vitreousRefractiveIndex = eye.getVitreousRefractiveIndex();

        for i = 1:length(loopSpan)
                
            eye = eye.setVitreousRefractiveIndex(vitreousRefractiveIndex + loopSpan(i));

            refractiveErrors(i) = eye.getRefractiveError();
            dioptricLengths(i)  = eye.getDioptricLength();
            equivalentPowers(i) = eye.getEquivalentPower();
            
        end
        
        title('Vitreous Refractive Index');
        xlabel('Vitreous Refractive Index');
        
        yyaxis left;
        plot(vitreousRefractiveIndex + loopSpan, refractiveErrors);
        ylabel('Refractive Error [D]');

        yyaxis right;
        ratio = dioptricLengths./equivalentPowers;
        plot(vitreousRefractiveIndex + loopSpan, ratio); hold on;
        plot([vitreousRefractiveIndex vitreousRefractiveIndex], [min(ratio) max(ratio)], 'k'); hold off;
        ylabel('Dioptric Length / Equiv. Power'); ylim([min(ratio) max(ratio)]);
    
        grid on; grid minor;
        
        
    %% Figure 3 - Radial    
        
    fh = figure(3);
    set(fh, 'Units', 'normalized');
    set(fh, 'Position', [0.1 0.1 0.8 0.8]);
    set(fh, 'Name', 'Radial');
    
    subplot(2,2,1);
        loopSpan = -2:0.1:2;
        
        refractiveErrors = zeros(length(loopSpan),1);
        dioptricLengths  = zeros(length(loopSpan),1);
        equivalentPowers = zeros(length(loopSpan),1);
        
        eye = EyeModel;
        cornealRadius = eye.getCornealRadius();

        for i = 1:length(loopSpan)
                
            eye = eye.setCornealRadius(cornealRadius + loopSpan(i));

            refractiveErrors(i) = eye.getRefractiveError();
            dioptricLengths(i)  = eye.getDioptricLength();
            equivalentPowers(i) = eye.getEquivalentPower();
            
        end
        
        title('Corneal Radius');
        xlabel('Corneal Radius [mm]');
        
        yyaxis left;
        plot(cornealRadius + loopSpan, refractiveErrors);
        ylabel('Refractive Error [D]');

        yyaxis right;
        ratio = dioptricLengths./equivalentPowers;
        plot(cornealRadius + loopSpan, ratio); hold on;
        plot([cornealRadius cornealRadius], [min(ratio) max(ratio)], 'k'); hold off;
        ylabel('Dioptric Length / Equiv. Power'); ylim([min(ratio) max(ratio)]);
        legend('Refractive Error [D]', 'Dioptric Length / Equiv. Power', 'Emmetropic Eye', 'Location', 'Best');
        
        grid on; grid minor;

    subplot(2,2,2);
        loopSpan = -2:0.1:2;
        
        refractiveErrors = zeros(length(loopSpan),1);
        dioptricLengths  = zeros(length(loopSpan),1);
        equivalentPowers = zeros(length(loopSpan),1);
        
        eye = EyeModel;
        anteriorLensRadius = eye.getAnteriorLensRadius();

        for i = 1:length(loopSpan)
                
            eye = eye.setAnteriorLensRadius(anteriorLensRadius + loopSpan(i));

            refractiveErrors(i) = eye.getRefractiveError();
            dioptricLengths(i)  = eye.getDioptricLength();
            equivalentPowers(i) = eye.getEquivalentPower();
            
        end
        
        title('Anterior Lens Radius');
        xlabel('Anterior Lens Radius [mm]');
        
        yyaxis left;
        plot(anteriorLensRadius + loopSpan, refractiveErrors);
        ylabel('Refractive Error [D]');

        yyaxis right;
        ratio = dioptricLengths./equivalentPowers;
        plot(anteriorLensRadius + loopSpan, ratio); hold on;
        plot([anteriorLensRadius anteriorLensRadius], [min(ratio) max(ratio)], 'k'); hold off;
        ylabel('Dioptric Length / Equiv. Power'); ylim([min(ratio) max(ratio)]);
        
        grid on; grid minor;
    
    subplot(2,2,3);
        loopSpan = -2:0.1:2;
        
        refractiveErrors = zeros(length(loopSpan),1);
        dioptricLengths  = zeros(length(loopSpan),1);
        equivalentPowers = zeros(length(loopSpan),1);
        
        eye = EyeModel;
        posteriorLensRadius = eye.getPosteriorLensRadius();

        for i = 1:length(loopSpan)
                
            eye = eye.setPosteriorLensRadius(posteriorLensRadius + loopSpan(i));

            refractiveErrors(i) = eye.getRefractiveError();
            dioptricLengths(i)  = eye.getDioptricLength();
            equivalentPowers(i) = eye.getEquivalentPower();
            
        end
        
        title('Posterior Lens Radius');
        xlabel('Posterior Lens Radius [mm]');
        
        yyaxis left;
        plot(posteriorLensRadius + loopSpan, refractiveErrors);
        ylabel('Refractive Error [D]');

        yyaxis right;
        ratio = dioptricLengths./equivalentPowers;
        plot(posteriorLensRadius + loopSpan, ratio); hold on;
        plot([posteriorLensRadius posteriorLensRadius], [min(ratio) max(ratio)], 'k'); hold off;
        ylabel('Dioptric Length / Equiv. Power'); ylim([min(ratio) max(ratio)]);
    
        grid on; grid minor;
        
    %% Figure 4a - Uniform Scaling (refr. error)
        
    fh = figure(4);
    set(fh, 'Units', 'normalized');
    set(fh, 'Position', [0.1 0.1 0.8 0.8]);
    set(fh, 'Name', 'Uniform Scaling');
    
    subplot(2,2,1);
    
        eye = EyeModel;
        parameter = eye.getAnteriorChamberDepth();
        parameterSpan = [-2, -1, 0, 1, 2];
        legendEntries = cell(length(parameterSpan),1);
        
        loopSpan = 0.5:0.05:2.0;
        
        for p = 1:length(parameterSpan)
            
            refractiveErrors = zeros(length(loopSpan),1);

            eye = EyeModel;
            eye = eye.setAnteriorChamberDepth(parameter + parameterSpan(p));

            for i = 1:length(loopSpan)

                eye = eye.rescale(loopSpan(i));

                refractiveErrors(i) = eye.getRefractiveError();

            end
            
            legendEntries{p} = sprintf('%3.3f mm', parameter + parameterSpan(p));
            plot(loopSpan, refractiveErrors); hold on;

        end
        
        hold off;
        title('Anterior Chamber Depth (parameter)');
        xlabel('Uniform Scale [x]');
        ylabel('Refractive Error [D]');

        legend(legendEntries, 'Location', 'Best');
        
        grid on; grid minor;
        
    subplot(2,2,2);
    
        eye = EyeModel;
        parameter = eye.getLensThickness();
        parameterSpan = [-2, -1, 0, 1, 2];
        legendEntries = cell(length(parameterSpan),1);
        
        loopSpan = 0.5:0.05:2.0;
        
        for p = 1:length(parameterSpan)
            
            refractiveErrors = zeros(length(loopSpan),1);

            eye = EyeModel;
            eye = eye.setLensThickness(parameter + parameterSpan(p));

            for i = 1:length(loopSpan)

                eye = eye.rescale(loopSpan(i));

                refractiveErrors(i) = eye.getRefractiveError();

            end
            
            legendEntries{p} = sprintf('%3.3f mm', parameter + parameterSpan(p));
            plot(loopSpan, refractiveErrors); hold on;

        end
        
        hold off;
        title('Lens Thickness (parameter)');
        xlabel('Uniform Scale [x]');
        ylabel('Refractive Error [D]');

        legend(legendEntries, 'Location', 'Best');
        
        grid on; grid minor;
   
    subplot(2,2,3);
    
        eye = EyeModel;
        parameter = eye.getVitrealChamberDepth();
        parameterSpan = [-2, -1, 0, 1, 2];
        legendEntries = cell(length(parameterSpan),1);
        
        loopSpan = 0.5:0.05:2.0;
        
        for p = 1:length(parameterSpan)
            
            refractiveErrors = zeros(length(loopSpan),1);

            eye = EyeModel;
            eye = eye.setVitrealChamberDepth(parameter + parameterSpan(p));

            for i = 1:length(loopSpan)

                eye = eye.rescale(loopSpan(i));

                refractiveErrors(i) = eye.getRefractiveError();

            end
            
            legendEntries{p} = sprintf('%3.3f mm', parameter + parameterSpan(p));
            plot(loopSpan, refractiveErrors); hold on;

        end
        
        hold off;
        title('Vitreal Chamber Depth (parameter)');
        xlabel('Uniform Scale [x]');
        ylabel('Refractive Error [D]');

        legend(legendEntries, 'Location', 'Best');
        
        grid on; grid minor;
        
    subplot(2,2,4);
        
        parameterSpan = [0.8, 0.9, 1.0, 1.1, 1.2];
        legendEntries = cell(length(parameterSpan),1);
        
        loopSpan = 0.5:0.05:2.0;
        
        for p = 1:length(parameterSpan)
        
            refractiveErrors = zeros(length(loopSpan),1);
            
            eye = EyeModel;
            eye = eye.setAxialElongation(parameterSpan(p));

            for i = 1:length(loopSpan)

                eye = eye.rescale(loopSpan(i));

                refractiveErrors(i) = eye.getRefractiveError();

            end

            legendEntries{p} = sprintf('%3.3f x', parameterSpan(p));
            plot(loopSpan, refractiveErrors); hold on;
            
        end
        
        hold off;
        title('Axial Elongation (parameter)');
        
        xlabel('Uniform Scale [x]');
        ylabel('Refractive Error [D]');

        legend(legendEntries, 'Location', 'Best');
        
        grid on; grid minor;
   
%% Figure 4b - Uniform Scaling (ratio)
        
    fh = figure(5);
    set(fh, 'Units', 'normalized');
    set(fh, 'Position', [0.1 0.1 0.8 0.8]);
    set(fh, 'Name', 'Uniform Scaling');
    
    subplot(2,2,1);
    
        eye = EyeModel;
        parameter = eye.getAnteriorChamberDepth();
        parameterSpan = [-2, -1, 0, 1, 2];
        legendEntries = cell(length(parameterSpan),1);
        
        loopSpan = 0.5:0.05:2.0;
        
        for p = 1:length(parameterSpan)

            dioptricLengths  = zeros(length(loopSpan),1);
            equivalentPowers = zeros(length(loopSpan),1);

            eye = EyeModel;
            eye = eye.setAnteriorChamberDepth(parameter + parameterSpan(p));

            for i = 1:length(loopSpan)

                eye = eye.rescale(loopSpan(i));

                dioptricLengths(i)  = eye.getDioptricLength();
                equivalentPowers(i) = eye.getEquivalentPower();

            end
            
            legendEntries{p} = sprintf('%3.3f mm', parameter + parameterSpan(p));
            ratio = dioptricLengths./equivalentPowers;
            plot(loopSpan, ratio); hold on;

        end
        
        hold off;
        title('Anterior Chamber Depth (parameter)');
        xlabel('Uniform Scale [x]');
        ylabel('Dioptric Length / Equiv. Power');

        legend(legendEntries, 'Location', 'Best');
        
        grid on; grid minor;
        
    subplot(2,2,2);
    
        eye = EyeModel;
        parameter = eye.getLensThickness();
        parameterSpan = [-2, -1, 0, 1, 2];
        legendEntries = cell(length(parameterSpan),1);
        
        loopSpan = 0.5:0.05:2.0;
        
        for p = 1:length(parameterSpan)

            dioptricLengths  = zeros(length(loopSpan),1);
            equivalentPowers = zeros(length(loopSpan),1);

            eye = EyeModel;
            eye = eye.setLensThickness(parameter + parameterSpan(p));

            for i = 1:length(loopSpan)

                eye = eye.rescale(loopSpan(i));

                dioptricLengths(i)  = eye.getDioptricLength();
                equivalentPowers(i) = eye.getEquivalentPower();

            end
            
            legendEntries{p} = sprintf('%3.3f mm', parameter + parameterSpan(p));
            ratio = dioptricLengths./equivalentPowers;
            plot(loopSpan, ratio); hold on;

        end
        
        hold off;
        title('Lens Thickness (parameter)');
        xlabel('Uniform Scale [x]');
        ylabel('Dioptric Length / Equiv. Power');

        legend(legendEntries, 'Location', 'Best');
        
        grid on; grid minor;
   
    subplot(2,2,3);
    
        eye = EyeModel;
        parameter = eye.getVitrealChamberDepth();
        parameterSpan = [-2, -1, 0, 1, 2];
        legendEntries = cell(length(parameterSpan),1);
        
        loopSpan = 0.5:0.05:2.0;
        
        for p = 1:length(parameterSpan)
            
            dioptricLengths  = zeros(length(loopSpan),1);
            equivalentPowers = zeros(length(loopSpan),1);

            eye = EyeModel;
            eye = eye.setVitrealChamberDepth(parameter + parameterSpan(p));

            for i = 1:length(loopSpan)

                eye = eye.rescale(loopSpan(i));

                dioptricLengths(i)  = eye.getDioptricLength();
                equivalentPowers(i) = eye.getEquivalentPower();

            end
            
            legendEntries{p} = sprintf('%3.3f mm', parameter + parameterSpan(p));
            ratio = dioptricLengths./equivalentPowers;
            plot(loopSpan, ratio); hold on;

        end
        
        hold off;
        title('Vitreal Chamber Depth (parameter)');
        xlabel('Uniform Scale [x]');
        ylabel('Dioptric Length / Equiv. Power');

        legend(legendEntries, 'Location', 'Best');
        
        grid on; grid minor;
        
    subplot(2,2,4);
    
        parameterSpan = [0.8, 0.9, 1.0, 1.1, 1.2];
        legendEntries = cell(length(parameterSpan),1);
        
        loopSpan = 0.5:0.05:2.0;
        
        for p = 1:length(parameterSpan)
        
            dioptricLengths  = zeros(length(loopSpan),1);
            equivalentPowers = zeros(length(loopSpan),1);
            
            eye = EyeModel;
            eye = eye.setAxialElongation(parameterSpan(p));

            for i = 1:length(loopSpan)

                eye = eye.rescale(loopSpan(i));

                dioptricLengths(i)  = eye.getDioptricLength();
                equivalentPowers(i) = eye.getEquivalentPower();

            end

            legendEntries{p} = sprintf('%3.3f x', parameterSpan(p));
            ratio = dioptricLengths./equivalentPowers;
            plot(loopSpan, ratio); hold on;
            
        end
        
        hold off;
        title('Axial Elongation (parameter)');
        
        xlabel('Uniform Scale [x]');
        ylabel('Dioptric Length / Equiv. Power');

        legend(legendEntries, 'Location', 'Best');
        
        grid on; grid minor;
