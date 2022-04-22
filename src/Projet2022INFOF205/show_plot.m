function show_plot(nsteps, xhatE, vxhatE)
%nsteps
x=0:1:nsteps;

figure(1)
%plotvxhatE = 0:0.1:max(vxhatE);
plot(x,vxhatE,'DisplayName','vxhatE');
hold on
%scatter(x,xhatE,'filled','DisplayName','xhatE');
plot(x,xhatE,'DisplayName','xhatE');
hold off
grid
%legend('vxhatE','xhatE')
legend