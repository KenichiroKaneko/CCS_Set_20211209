load('vars_2nd')
close all;

ind_b = 1:SENSOR_TPRB.NUM;
ind_f = SENSOR_TPRB.NUM+1:SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM;
BZ2 = BB(ind_b)';
FL2 = BB(ind_f)';
FLC = FC(ind_f);

F = FF+FC;

% figure()
% subplot(3,1,1);hold on;plot(FL);plot(SENSOR_FLXLP.FLXLP,'o-');plot(FL2./(SENSOR_FLXLP.R).^2./2./pi,'o-');legend('FL','SENFL','FL2');
% subplot(3,1,2);hold on;plot(BZ);plot(BZ2);legend('BZ','BZ2');
% subplot(3,1,3);hold on;plot(SENSOR_TPRB.TPRB);plot(BZ2);legend('BZ','BZ2');


figure()
subplot(3,1,1); hold on;
plot(SENSOR_FLXLP.FLXLP)
plot((FL2+FC(ind_f))./factors(2)*2*pi, 'o')

subplot(3,1,2); hold on;
plot(); 
plot(BZ2);
legend('BZ','BZ2');

subplot(3,1,3); hold on;
plot(SENSOR_TPRB.TPRB); plot(BZ2);
legend('BZ','BZ2');