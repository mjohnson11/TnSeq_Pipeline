curr_dir = '/Users/skryazhi/Google Drive File Stream/Team Drives/SKLAB/Projects/TnSeq/Milo/BT_Experiment_Nov_2017/new_s_measurements';
s_column_ix = 9; % column in the data files where the selection coefficient is recorded

%%% v. 3. 2019-05-13. Cleaning up the code
%%% v. 2. 2017-12-28. Excluding Segregant 9 from the analysis

%% Read in the selection coefficients for Milo's BCs

meta.reps = 1;
meta.libs = {};
meta.strains = {};   % list of strains
meta.strainlib = []; % list of library ids for each strain

meta.nstrain = length(meta.strains); % number of strains
meta.nbc = nan(meta.nstrain,1);      % number of unq BCs per strain

bc = cell(meta.nstrain,1); % data containing the barcode information

bc_edge = {};  % First column contains the BC, second column contains the corresponding edge (mutation identifier)

filelist = dir([curr_dir '/*.csv']);

for ifile = 1:length(filelist)
    filename = sprintf('%s/%s', curr_dir, filelist(ifile).name);
    newStr = strsplit(filelist(ifile).name, '_');
    
    if strcmp(newStr{1}, '9')
        continue;
    end
    strain = sprintf('Seg %02d', str2double(newStr{1}));
    assay = newStr{2};
    lib = newStr{3};
    
    istrain = find(strcmp( meta.strains, strain));
    irep = 2;
    if isempty(istrain)
        meta.strains = [meta.strains, strain];
        istrain = length(meta.strains);
        irep = 1;
    end
    
    ixl = find(strcmp( meta.libs, lib));
    if isempty(ixl)
        meta.libs = [meta.libs, lib];
        ixl = length(meta.libs);
    end
    meta.strainlib(istrain, irep) = ixl;
        
    fid = fopen(filename, 'r');
    fs = textscan(fid, '%s%s%d%d%d%d%d%d%f%f%d', 'delimiter', ',', 'headerlines', 1);

    fclose(fid);
    
    if length(bc) >= istrain
        bc{istrain}.bcstr = [bc{istrain}.bcstr; fs{1} fs{2}];
        bc{istrain}.s = [bc{istrain}.s ; fs{s_column_ix}];
        bc{istrain}.lib = [bc{istrain}.lib ; ixl*ones(length(fs{s_column_ix}), 1)];
    else
        bc{istrain}.bcstr = [fs{1} fs{2}];
        bc{istrain}.s = fs{s_column_ix};
        bc{istrain}.lib = ixl*ones(length(fs{s_column_ix}), 1);
    end
    
    bc_edge = [bc_edge; fs{1} fs{2}];
end
clear ifile newStr strain irep assay lib fid fs filelist ix ixl;



%% identifying the mutation id that corresponds to each BC

[meta.strains, ix] = sort(meta.strains);
meta.strainlib = meta.strainlib(ix,:);
bc = bc(ix);
clear ix;

meta.nlib = length(meta.libs); % number of different libraries
meta.nrep = length(meta.reps); % number of replicate measurements per BC
meta.nstrain = length(meta.strains);  % number of strains
meta.strainlibn = size(meta.strainlib,2); % number of libraries per strain

mut.edge = unique( bc_edge(:,2) ); % edge of each mutation
meta.nmut = length(mut.edge);      % number of mutations

for istrain = 1:meta.nstrain
    meta.nbc(istrain) = size( bc{istrain}.s, 1);
    bc{istrain}.mut_ix = zeros(meta.nbc(istrain),1);
    
    if length(unique(bc{istrain}.bcstr(:,1))) ~= length(bc{istrain}.bcstr)
        fprintf('Strain %s has non-unique BCs\n', meta.strains{istrain});
    end
    
    for ilib = 1:meta.strainlibn
        ixl = meta.strainlib(istrain,ilib);
        
        TF = bc{istrain}.lib == ixl;
        
        bc{istrain}.mut_ix(TF) = cellfun(@(edge) find(strcmp(edge, mut.edge)), bc{istrain}.bcstr(TF,2),'UniformOutput', true);                
    end
end

save('sk_data.mat', 'bc', 'mut', 'meta');

clear istrain ilib ixl TF;



%% Populate data for mutations

load('sk_data.mat');

mut.s = cell(meta.nmut, meta.nstrain); % selection coefficient of all barcodes corresponding to each mutation
% mut.s{imut,istrain}(ibc) is the selection coefficient in strain istrain
% of barcode ibc corresponding to mutation imut

for imut = 1:meta.nmut
    for istrain = 1:meta.nstrain
        
        curr_mut_strain_nbc = nnz(bc{istrain}.mut_ix == imut );
        
        mut.s{imut,istrain} = nan( curr_mut_strain_nbc , 1);
        
        TF = bc{istrain}.mut_ix == imut;
        
        nbc = nnz(TF);
        if nbc > 0
            mut.s{imut,istrain} = bc{istrain}.s(TF) * 100;
        end
    end
end
clear curr_mut_strain_nbc imut istrain ilib TF k nbc;

save('sk_data.mat', 'bc', 'mut', 'meta');



%% Fitting data (simple)

load('sk_data.mat');
    
myfit = fit_data_norm_ms( mut.s );
    
save('sk_data.mat', 'bc', 'mut', 'meta', 'myfit');


%% 13-A Plot variation in fitness estimates for all mutations (on both days)

% load('sk_data.mat');

IFPLOT = 1;

batchsize = 5 * 4 * 2^2;

Xvec = 1:meta.nstrain;
xticklabel = [1:8 10:20];

cc = [
    0, 114, 178;    % blue
    213, 94, 0;     % vermillion
    86, 180, 233;   % sky blue
    230 159, 0;     % orange
    204, 121, 167;   % raddish purple
    0, 158, 115;    % bluish green
    240, 228, 66   % yellow
    ]./256;


%%% Specify the following dimensions:
fdim.spwa = 4; % subplotwidth in cm
fdim.spha = 4; % subplotheight in cm

fdim.nx = 10; % number of panels along the horizontal dimension
fdim.ny = 8; % number of panels along the vertical dimension

fdim.xma = [0.5 0.5]; % left right horizontal margin in cm
fdim.yma = [1 0.5]; % bottom top vertical margin cm

fdim.dxa = 0.2; % horizontal distance between panels in cm
fdim.dya = 0.2; % vertical distance between panels in cm

fdim.tickfs = 6;
fdim.labelfs = 8;

%%% These will be computed automatically:
fdim.fw = fdim.spwa * fdim.nx + fdim.dxa * (fdim.nx - 1) + sum(fdim.xma);
fdim.fh = fdim.spha * fdim.ny + fdim.dya * (fdim.ny - 1) + sum(fdim.yma);

fdim.spwr = fdim.spwa / fdim.fw;
fdim.sphr = fdim.spha / fdim.fh;
fdim.xmr = fdim.xma / fdim.fw;
fdim.ymr = fdim.yma / fdim.fh;
fdim.dxr = fdim.dxa / fdim.fw;
fdim.dyr = fdim.dya / fdim.fh;

k = 1;

while k <= meta.nmut
        
    clf;
    
    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [0 0 fdim.fw fdim.fh]);
    
    fdim.spxvec = fdim.xmr(1) + fdim.spwr * ( 0:(fdim.nx-1) ) + fdim.dxr * ( 0:(fdim.nx-1) );
    fdim.spyvec = fdim.ymr(1) + fdim.sphr * ( (fdim.ny-1):-1:0 ) + fdim.dyr * ( (fdim.ny-1):-1:0 );
    
    for imut = k:min(k+batchsize-1, meta.nmut)
        
        iy = floor( (imut-k)/fdim.nx ) + 1;
        ix = mod( imut-k, fdim.nx ) + 1;
        
        subplot('Position', [fdim.spxvec(ix) fdim.spyvec(iy) fdim.spwr fdim.sphr]),
        hold on, box on;
        set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
        
        Ymut = squeeze(mut.s(imut,:));
        
        TF = cellfun(@(x) isempty(x) | all(isnan(x(:))), Ymut); %are they all empty or contain only NaNs?
        
        if all(TF)
            continue;
        end
        
        % Plot ML fit
        % plot(Xvec, myfit.mumat(imut,:), '-', 'Color', 'k', 'LineWidth', 2);
        TF = myfit.pvals(imut,:) < 1e-3;
        plot(Xvec(TF), myfit.mumat(imut,TF), 'd', 'MarkerSize', 8, ...
            'MarkerFaceColor', cc(2,:),  'MarkerEdgeColor', 'none');
        plot(Xvec(~TF), myfit.mumat(imut,~TF), 'd', 'MarkerSize', 8, ...
            'MarkerFaceColor', cc(3,:), 'MarkerEdgeColor', 'none');
        
        for istrain = 1:meta.nstrain
            %         if pval < 1e-3
            %             text( 1.5, myfit.mumat(imut,istrain) + 1, sprintf('P = %.1e', pval),...
            %                 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center' );
            %         end
            
            Ycurr = Ymut{istrain};
            
            % Plot measured values
            Y = Ycurr;
            
            X = Xvec(istrain)*ones(length(Y),1);% + 0.1*rand(length(Y),1);
            
            plot( X', Y', 'o', 'Color', 'k', ...
                'MarkerSize', 4, 'MarkerFaceColor', 'k',...
                'MarkerEdgeColor', 'none');
        end               
        
        % plot( [0.5 meta.nstrain+0.5], [0 0], '-', 'Color', 0.7*[1 1 1]);
        set(gca, 'XLim', [0 meta.nstrain+1], 'YLim', [-15 10],...
            'XTick', 1:meta.nstrain, 'YTick', -15:5:10,...
            'XGrid', 'on', 'YGrid', 'on');
        
        text(10, 9.5, sprintf('Mut #%d', imut), 'FontSize', 12, 'FontName', 'Helvetica',...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

        if ix > 1
            set(gca, 'YTickLabel', {});
        end
        if iy == fdim.ny
            set(gca, 'XTickLabel', xticklabel );
            xtickangle(90);
        else
            set(gca, 'XTickLabel', {});
        end
        
        %         if istrain == 3
        %             text(1.5, 13, sprintf('Mut # %d', imut ),...
        %                 'FontName', 'Helvetica', 'FontSize', 20, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
        %         end

    end
    
    if IFPLOT
        filename = sprintf('figures/muts_%d--%d.eps', k, min(k+batchsize-1, meta.nmut));
        saveas(gcf, filename, 'epsc');
    end
    
    k = k+batchsize;
end

clear fdim cc;
clear ienv istrain ilib it k id irep X Y TF;
clear b bint r rint stats;
clear Xm Ym Xse Yse pval Ymut Ycurr imut dL ilib marker x0




%% Finding statistically significant mutational effects using Benjamini-Hochberg procedure
% load('sk_data.mat');
alpha = 0.001;

myfit.ISBEN = false( size(myfit.mumat) );
myfit.ISDEL = false( size(myfit.mumat) );

mut_type_cnts = zeros(3, meta.nstrain);

for istrain = 1:meta.nstrain
    [pvals, ix] = sort( myfit.pvals(:,istrain), 'ascend');
    
    TF = ~isnan(pvals);
    ntested = nnz(TF);
    pvals = pvals(TF);
    ix = ix(TF);
    ix1 = ix( pvals <= alpha/ntested*(1:ntested)');
    
    TF = false(meta.nmut,1);
    TF(ix1) = true;
    
    myfit.ISBEN(:,istrain) = myfit.mumat(:, istrain) > 0 & TF;
    myfit.ISDEL(:,istrain) = myfit.mumat(:, istrain) < 0 & TF;
    
    TFBEN = myfit.ISBEN(:,istrain);
    TFDEL = myfit.ISDEL(:,istrain);
    
    nben = nnz(myfit.ISBEN(:,istrain));
    ndel = nnz(myfit.ISDEL(:,istrain));
    nneut = ntested - nben - ndel;
    
    mut_type_cnts(:, istrain) = [nneut; ndel; nben];
    
    fprintf( 'Strain %s:\n\t', meta.strains{istrain} );
    fprintf( '%d mut tested\n\t', ntested);
    fprintf( '%d (%.2f%%) ben mut, <s> = %.1f +- %.1f%%\n\t', nben,...
        nben/ntested*100, nanmean( myfit.mumat(TFBEN,istrain) ),...
        nanstd( myfit.mumat(TFBEN,istrain) )  );
    fprintf( '%d (%.2f%%) del mut, <s> = %.1f +- %.1f%%\n\n', ndel,...
        ndel/ntested*100, nanmean( myfit.mumat(TFDEL,istrain) ),...
        nanstd( myfit.mumat(TFDEL,istrain) )  );
end

clear istrain pvals ix TF ntested ix1 nben ndel nneut TFBEN TFDEL;

save('sk_data.mat', 'bc', 'mut', 'meta', 'myfit');



%% Writing results into a text file
filename = 'sk_data.tab';
fid = fopen(filename, 'w');

fprintf(fid, '\t\tEstimated selection coefficient of mutation\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t');
fprintf(fid, 'P-value\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t');
fprintf(fid, 'Is significant at alpha = 0.001 after multiple testing correction?\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n');

fprintf(fid, 'ID\tEDGE\t');
for k=1:3
    for istrain = 1:meta.nstrain
        fprintf(fid, '%s\t', meta.strains{istrain});
    end
end
fprintf(fid, 'Number of strains where significant\n');

% sort mutations by the number of strains in which they are significant:
nSignVec = sum(myfit.ISBEN | myfit.ISDEL, 2);
[tmp, ix] = sort(nSignVec, 'descend');

for imut = ix'
    fprintf(fid, '%d\t%s\t', imut, mut.edge{imut} );
    for istrain = 1:meta.nstrain
        fprintf(fid, '%.2f\t', myfit.mumat(imut, istrain) );
    end
    for istrain = 1:meta.nstrain
        fprintf(fid, '%.2e\t', myfit.pvals(imut, istrain) );
    end
    for istrain = 1:meta.nstrain
        fprintf(fid, '%d\t', myfit.ISBEN(imut, istrain) || myfit.ISDEL(imut, istrain) );
    end
    fprintf(fid, '%d\n', nSignVec(imut));
end

fclose(fid);


