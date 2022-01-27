close all
load('nokaizoudo')
figure()
subplot(1,2,1)
semilogy(param.err)
xlabel('平衡計算の反復回数')
ylabel('誤差')

ax = gca;
ax.FontSize = 18;

load('kaisoudoUPUP')
subplot(1,2,2)
semilogy(param.err)
xlabel('平衡計算の反復回数')
ylabel('誤差')
ax = gca;
ax.FontSize = 18;