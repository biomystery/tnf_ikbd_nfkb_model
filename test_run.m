% BMDM model 
id = struct;
%  Vary i parameters
alter =  [...
 ];
if ~isempty(alter)
    id.inputPid = alter(:,1)';
    id.inputP  = alter(:,2)';
end

% Vary n parameters
alter =       [...
    ];
if ~isempty(alter)
    id.inputvPid = alter(:,1)';
    id.inputvP  = alter(:,2)';
end


doses = [.001, .01, .1, 1, 10, 100, 1000] ; %ng/ml 
mod_colormap = jet(length(doses));
%divergingmap(0:1/(length(doses)-1):1,[12 12 77]/255,[158 ...
%                    4 0]/255);

id.output = {'TNF','IKKK','IkBd','IKK','NFkBn','a20t','a20'}; % output names are in getInit.m
id.DT = 0.05; 
id.sim_time = 60*12;
[n,i] = getRateParams(); % Vary "i" params with id.inputPid/inputP and "n" params with id.inputvPid/inputvP

% Simulate
wt_sim = zeros(length(id.output),round(id.sim_time/id.DT)+1, ...
               length(doses));

for idx_d = 1:length(doses)
    run_id = id;
    run_id.dose = doses(idx_d);
    wt_sim(:,:,idx_d) = getSimData(run_id);
end

% plot
n_output = length(id.output);
figure('Position',[          60         227        1115         433]);
for idx_o = 1:n_output
    subplot(2,4,idx_o)
    hold on
    for idx_d = 1:length(doses)
        plot(0:id.DT:720,wt_sim(idx_o,:,idx_d),'linewidth',1.5,'color',mod_colormap(idx_d,:))
    end
    set(gca,'XLim',[0 720],'Xtick',0:120:720,'Xticklabel',0:2:12)  
    xlabel('Time (hr)')
    title(id.output{idx_o})
end
%subplot(2,4,8) 
%colormap(mod_colormap)
%hcb=colorbar;
%set(hcb,'YTick',log10(doses),'YTicklabel',{'.001','.01','.1','1','10','100','1000'})
%ylabel(hcb,'TNF doses(ng/ml)') 

% save 
print('./figs/dose_sim.png','-dpng')
