function plotting_param_sens_lam(Diff01,Diff05,Diff1,Diff15,Diff2,Diff100)

cmap = [215,48,39;...
252,141,89;...
254,224,144;...
224,243,248;...
145,191,219;...
69,117,180]/255;

figure
hold on 
plot(linspace(0,72,72)/12,Diff01.uninfected_live,'Color',cmap(1,:),'LineWidth',2)
plot(linspace(0,72,72)/12,Diff05.uninfected_live,'Color',cmap(2,:),'LineWidth',2)
plot(linspace(0,72,72)/12,Diff1.uninfected_live,'Color',cmap(3,:),'LineWidth',2)
plot(linspace(0,72,72)/12,Diff15.uninfected_live,'Color',cmap(4,:),'LineWidth',2)
plot(linspace(0,72,72)/12,Diff2.uninfected_live,'Color',cmap(5,:),'LineWidth',2)
plot(linspace(0,72,72)/12,Diff100.uninfected_live,'Color',cmap(6,:),'LineWidth',2)
xlabel('Time (days)')
ylabel('Cells')
title('Uninfected cells')
set(gca,'FontSize',18)
c = colorbar
colormap(cmap)
set(c,'ticks',[0.0833,0.25,0.41667,0.5835,0.7502,0.9166],'ticklabels',{'\beta{\times}0.1','\beta{\times}0.5','\beta{\times}1','\beta{\times}1.5','\beta{\times}2','\beta{\times}10'})


figure
hold on 
plot(linspace(0,72,72)/12,Diff01.dead,'Color',cmap(1,:),'LineWidth',2)
plot(linspace(0,72,72)/12,Diff05.dead,'Color',cmap(2,:),'LineWidth',2)
plot(linspace(0,72,72)/12,Diff1.dead,'Color',cmap(3,:),'LineWidth',2)
plot(linspace(0,72,72)/12,Diff15.dead,'Color',cmap(4,:),'LineWidth',2)
plot(linspace(0,72,72)/12,Diff2.dead,'Color',cmap(5,:),'LineWidth',2)
plot(linspace(0,72,72)/12,Diff100.dead,'Color',cmap(6,:),'LineWidth',2)
xlabel('Time (days)')
ylabel('Cells')
title('Dead cells')
set(gca,'FontSize',18)
c = colorbar
colormap(cmap)
set(c,'ticks',[0.0833,0.25,0.41667,0.5835,0.7502,0.9166],'ticklabels',{'\beta{\times}0.1','\beta{\times}0.5','\beta{\times}1','\beta{\times}1.5','\beta{\times}2','\beta{\times}10'})

figure
hold on 
plot(linspace(0,72,72)/12,Diff01.extracellular_virus,'Color',cmap(1,:),'LineWidth',2)
plot(linspace(0,72,72)/12,Diff05.extracellular_virus,'Color',cmap(2,:),'LineWidth',2)
plot(linspace(0,72,72)/12,Diff1.extracellular_virus,'Color',cmap(3,:),'LineWidth',2)
plot(linspace(0,72,72)/12,Diff15.extracellular_virus,'Color',cmap(4,:),'LineWidth',2)
plot(linspace(0,72,72)/12,Diff2.extracellular_virus,'Color',cmap(5,:),'LineWidth',2)
plot(linspace(0,72,72)/12,Diff100.extracellular_virus,'Color',cmap(6,:),'LineWidth',2)
xlabel('Time (days)')
ylabel('Cells')
title('Oncolytic virus')
set(gca,'FontSize',18)
c = colorbar
colormap(cmap)
set(c,'ticks',[0.0833,0.25,0.41667,0.5835,0.7502,0.9166],'ticklabels',{'\beta{\times}0.1','\beta{\times}0.5','\beta{\times}1','\beta{\times}1.5','\beta{\times}2','\beta{\times}10'})

end