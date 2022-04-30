function plot_energies(nsteps, xhatE, vxhatE, Epot, Ecin, Etot)
%nsteps
x=0:1:nsteps;

%figure("energies")
figure(2)


plot(x,xhatE,'DisplayName','xhatE');
hold on
plot(x,vxhatE,'DisplayName','vxhatE');
hold on
plot(x,Epot,'DisplayName','Epot');
hold on
plot(x,Ecin,'DisplayName','Ecin');
hold on
plot(x,Etot,'DisplayName','Etot');
hold off
grid
legend