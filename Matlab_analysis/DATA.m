%%                 CoRa ANALYSIS                  %%
%%% DATA: Import output files from  CoRa_Main.jl %%%
%%%       to Matlab to generate figures.         %%%
% Mariana G�mez-Schiavon
% March 2023

clear;
%motifs = {'ATFv1','ATFv2','FADv1','FADv2','BNFv1','BNFv2','FDPv1','FDPv2','FFLv1','BMFv1','BMFv2'};
motifs = {'FADv1'};
save Temp.mat
for m = 1:length(motifs)
    load Temp.mat
    % sim.mm = motifs{m}
    sim.mm = 'FADv1';     % Label for motif file
    sim.ex = 'Fig3';        % Label for parameters file
    sim.pp = 'mY';          % Label for perturbation type
    sim.ax = 'mY';          % Label for condition/range
    sim.an = 'OptCoRa';     % Chose analysis type (Options: ExSSs, CoRams, OptCoRa)
    curves = 1;             % 1 if curves exist in the data  
    npoints = 31;            % number of points of each curve 
    runs = [1];         % list of runs
    clear motifs m
    
    %% Load & parse data
    for run = 1:length(runs)
        data = ['OUT_', sim.an, '_', sim.mm, '_', sim.ex, '_', sim.pp, '_', sim.ax, '_', num2str(runs(run)), '.txt'];
        x = importdata(data,'\t',1);
        if(strcmp(sim.an,'ExSSs'))
            rho.cond   = x.textdata{1};
            rho.values = x.data(:,1);
            for i = 2:(size(x.data,2)-1)
                str = x.textdata{i};
                    [Is,Ie] = regexp(str,"_.+");
                ss.(str(1:3)).(str((Is+1):Ie)) = x.data(:,i);
            end
            rho.pert = x.textdata{end};
            CoRas = x.data(:,end);
            clear Is Ie i x str
        elseif(strcmp(sim.an,'CoRams'))
            for i = 1:size(x.data,2)
                str = x.textdata{i};
                [Is,Ie] = regexp(str,"\d+\.\d+");
                if(Is & Ie)
                    break;
                else
                    p.(str) = x.data(:,i);
                end
            end
            CoRas = x.data(:,i:end);
            rho.name   = sim.pp;
            rhoV = x.textdata(1,i:end);
            rho.values = zeros(size(rhoV));
            for i = 1:length(rhoV)
                rho.values(i) = str2num(rhoV{i});
            end
            clear i rhoV Ie Is x str
        elseif(strcmp(sim.an,'OptCoRa'))
            rho.name = sim.pp;
            x.data_c = x.data(:,1:(end-npoints));
            x.textdata_c = x.textdata(:,1:(end-npoints));
            epsT = x.textdata_c{end-1};
                [Is,Ie] = regexp(epsT,"\d\.\d+");
                epsT  = epsT(Is:Ie);
                clear Is Ie 
            R = x.data_c(:,1);
            if(length(R)==2*max(R))
                I = x.data(:,2);
                CoRai.min = x.data([I == 0],end);
                CoRai.thr = x.data([I == 0],end-1);
                CoRaf.min = x.data([I ~= 0],end);
                CoRaf.thr = x.data([I ~= 0],end-1);
                for i = 3:(size(x.data,2)-2)
                    pi.(x.textdata{i}) = x.data([I == 0],i);
                    pf.(x.textdata{i}) = x.data([I ~= 0],i);
                end
            else
                R0 = 1;
                if curves == 1
                    while(R0<=max(R))
                        CoRas(R0).curv = x.data([R == R0],(end-npoints+1):end);
                        R0 = R0 + 1;
                    end
                end
                x.data = x.data(:,1:(end-npoints));
                x.textdata = x.textdata(:,1:(end-npoints));
                R0 = 1;
                while(R0<=max(R))
                    I = x.data(:,2);
                    CoRas(R0).min = x.data([R == R0],end);
                    CoRas(R0).thr = x.data([R == R0],end-1);
                    for i = 3:(size(x.data,2)-2)
                        ps(R0).(x.textdata{i}) = x.data([R == R0],i);
                    end
                    R0 = R0 + 1;
                end
            end
            clear i x R I R0
        else
           %'ERROR: Undetermined analysis. Options: ExSSs, CoRams, OptCoRa'
        end
        save(cat(2,'DATA_',sim.an,'_',sim.mm,'_',sim.ex,'_',sim.pp,'_',sim.ax,'_',num2str(runs(run)),'.mat'))
        %delete(cat(2,'OUT_',sim.an,'_',sim.mm,'_',sim.ex,'_',sim.pp,'_',sim.ax,'_',run,'.txt'))
    end
end
delete Temp.mat

%% END