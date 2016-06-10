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
plot(exp051916.tps,exp051916.tc_data(:,idx(1)) ,'o-','linewidth',1.5,...
    'color',[30, 144, 255]/255)
hold on 
plot(exp051916.tps,exp051916.tc_data(:,idx(2)) ,'^-','linewidth',1.5,...
    'color',[0,191,255]/255)
if(i==1)
    legend({'rep1','rep2'})
end

title(exp051916.species_unique{i})
set(gca,'xtick',0:120:1440,'xticklabel',(0:120:1440)/60)
xlabel('Time (h)');ylabel('fold') 
end
print('./figs/exp051916.png','-dpng')

%% simulate once
id = struct;
%  Vary i parameters
alter =  [
%    33 1e-7;
%    35 1000; 
%    36 1e-5; 
%    37 0.1; 
 ];
if ~isempty(alter)
    id.inputPid = alter(:,1)';
    id.inputP  = alter(:,2)';
end

% Vary n parameters
alter =       [
%     1 3e-5;
%     %3 2e-6;
     72 3e-8; 
%    4 16;
     5 .3;
     6 .4; %txn e
%     73 .1;



    ];
if ~isempty(alter)
    id.inputvPid = alter(:,1)';
    id.inputvP  = alter(:,2)';
end


dose = exp051916.dose ; %ng/ml 

id.output = {'IkBa','IkBaNFkB','IkBan','IkBaNFkBn',...
    'IkBb','IkBbNFkB','IkBbn','IkBbNFkBn',...
    'IkBe','IkBeNFkB','IkBen','IkBeNFkBn',...
    'IkBd','IkBdNFkB','IkBdn','IkBdNFkBn'...
    'NFkBn','IKK','a20','IkBat','TNF','IkBbt',...
    'IkBet','IkBdt','IkBaNFkBn','IkBdNFkBn',...
    'IKK_off','IKK_i'}; % output names are in getInit.m
id.DT = 0.05; 
id.sim_time = exp051916.tps(end);


% Simulate
run_id = id;
run_id.dose = dose;
[wt_sim,v] = getSimData(run_id);

%
sim_data = zeros(size(wt_sim,2),6);
for i = 1:4
    tmp = (wt_sim(i*4-3,:) + wt_sim(i*4-2,:)+wt_sim(i*4-1,:)+ wt_sim(i*4,:));
    sim_data(:,i)= tmp/tmp(1);
    %sim_data(:,i)= (wt_sim(i*4-1,:) + wt_sim(i*4,:));
end
%
sim_data(:,5:9) = wt_sim(9:13,:)';  

% plot
%figure('position',[ 680         415        1108         563])
figure
tils.species = {exp051916.species_unique{:},'NFkBn','IKK','A20','IkBs','TNF'}; 
for i = 1:4
    subplot(2,2,i)
    idx =find(strcmp(exp051916.species,exp051916.species_unique{i})) ;
    plot(0:id.DT:id.sim_time,sim_data(:,i),'b-','linewidth',1.5)
    hold on
    plot(exp051916.tps,exp051916.tc_data(:,idx(1)) ,'o','linewidth',1.5,...
    'color',[30, 144, 255]/255)
    hold on 
    plot(exp051916.tps,exp051916.tc_data(:,idx(2)) ,'^','linewidth',1.5,...
    'color',[0,191,255]/255)

    %if (i <=4)
    %    tmp = (wt_sim(i*4-3,:) + wt_sim(i*4-2,:)+wt_sim(i*4-1,:)+ wt_sim(i*4,:));
        
  %      plot(0:id.DT:id.sim_time,(wt_sim(i*4-3,:)+wt_sim(i*4-2,:))/tmp(1),...
   %                    '--','linewidth',1.5,'color',[255 204 203]/255) % cyto
        %plot(0:id.DT:id.sim_time,wt_sim(i*2-1,:),...

%        plot(0:id.DT:id.sim_time,(wt_sim(i*4-1,:)+wt_sim(i*4,:))/tmp(1),...
 %           '--','linewidth',1.5,'color',[127 63 62]/255) % nuclear
    %end
%     if(i==5) 
%         plot(0:id.DT:id.sim_time,wt_sim(end-3:end-2,:),...
%             '--','linewidth',1.5,'color',[127 63 62]/255) % nuclear
%     end 
%    
%     if(i==6)
%         plot(0:id.DT:id.sim_time,wt_sim(end-1,:),...
%             'k--','linewidth',1.5) % ikk_off
%         plot(0:id.DT:id.sim_time,wt_sim(end,:),...
%             'b--','linewidth',1.5) % ikk_i
%     end 
%     if(i==8)
%         plot(0:id.DT:id.sim_time,wt_sim(14,:),...
%             'k--','linewidth',1.5) % ikbbt
%         plot(0:id.DT:id.sim_time,wt_sim(15,:),...
%             'g--','linewidth',1.5) % ikbet
%         plot(0:id.DT:id.sim_time,wt_sim(16,:),...
%             'b--','linewidth',1.5) % ikbdt
% 
%     end

    if(i==1)
        legend({'tot','cyto','nuclear'});
    end;
    
    title(tils.species{i})
    set(gca,'xtick',0:120:1440,'xticklabel',(0:120:1440)/60)
    xlabel('Time (h)');ylabel('fold')
end
print('./figs/exp051916_ctxnD_itxnBE.png','-dpng')

%% manual calibrations 
% 1. ikba induction is too high 


