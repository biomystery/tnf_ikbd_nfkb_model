%% read exp051916 data 
%  the relative quantification data of IkBs in response to 10 ng/ml
% TNF over time (Western blot of whole cell lysates, normalized to tubulin; 
% 2 biological replicates).

file = './data/160519_TNF time course_summary.xlsx'; 
[tc_data,reps,~]= xlsread(file,1,'C3:J14');
tps=xlsread(file,1,'B4:B14') ; % time points 
[~,species,~]=xlsread(file,1,'C2:J2') ; % time points 
exp051916=struct; 
exp051916.tc_data = tc_data; exp051916.reps=reps; 
exp051916.tps = tps; exp051916.sti = 'TNF'; exp051916.dose=10; 
exp051916.dose_unit='ng/ml'; exp051916.celltype='L929';
exp051916.species=species;
exp051916.species_unique=unique(species);
exp051916.species_unique_no = length(unique(species));

clear tc_data reps tps species;

% plot the exp051916 data 
for i = 1:exp051916.species_unique_no
subplot(2,2,i)
idx =find(strcmp(exp051916.species,exp051916.species_unique{i})) ;
plot(exp051916.tps,exp051916.tc_data(:,idx) ,'o-','linewidth',1.5)
title(exp051916.species_unique{i})
set(gca,'xtick',0:120:1440,'xticklabel',(0:120:1440)/60)
xlabel('Time (h)');ylabel('fold') 
end
print('./figs/exp051916.png','-dpng')

%% simulate once
id = struct;
%  Vary i parameters
alter =  [...
 ];
if ~isempty(alter)
    id.inputPid = alter(:,1)';
    id.inputP  = alter(:,2)';
end

% Vary n parameters
alter =       [
    1 3e-4;
    3 2e-6;
    6 .2;
    %72 2e-7;
    %12 90;
    %4 0.8; 
    %10 0;
    %75 0 ;
    ];
if ~isempty(alter)
    id.inputvPid = alter(:,1)';
    id.inputvP  = alter(:,2)';
end


dose = exp051916.dose ; %ng/ml 

id.output = {'IkBa','IkBan','IkBb','IkBbn','IkBe','IkBen','IkBd','IkBdn'}; % output names are in getInit.m
id.DT = 0.05; 
id.sim_time = exp051916.tps(end);
[n,i] = getRateParams(); % Vary "i" params with id.inputPid/inputP and "n" params with id.inputvPid/inputvP

% Simulate
run_id = id;
run_id.dose = dose;
wt_sim = getSimData(run_id);

%
sim_data = zeros(size(wt_sim,2),4);
for i = 1:4
    sim_data(:,i)= (wt_sim(i*2-1,:) + wt_sim(i*2,:))./(wt_sim(i*2-1,1) + wt_sim(i*2,1));
end


% plot
figure
for i = 1:exp051916.species_unique_no
subplot(2,2,i)
plot(0:id.DT:id.sim_time,sim_data(:,i),'r-','linewidth',1.5)
hold on 
plot(0:id.DT:id.sim_time,wt_sim(i*2-1,:)/(wt_sim(i*2-1,1) + wt_sim(i*2,1)),...
    '--','linewidth',1.5,'color',[255 204 203]/255) % cyto
plot(0:id.DT:id.sim_time,wt_sim(i*2,:)/(wt_sim(i*2-1,1) + wt_sim(i*2,1)),...
    '--','linewidth',1.5,'color',[127 63 62]/255) % nuclear
if(i==1) 
    legend({'tot','cyto','nuclear'}); 
end; 

title(exp051916.species_unique{i})
set(gca,'xtick',0:120:1440,'xticklabel',(0:120:1440)/60)
xlabel('Time (h)');ylabel('fold') 
end
print('./figs/exp051916_sim_basalAE_txnE.png','-dpng')

%% manual calibrations 
% 1. ikba induction is too high 


