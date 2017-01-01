a=3;
figure;
GlobalMins=runABC(a);
semilogy(GlobalMins,'k:');
%semilogy(mean(GlobalMins))
title('griewank函数的适应度值收敛趋势');
xlabel('迭代次数（cycles）');
ylabel('适应度（fitness）');
b=3;
GlobalMins=runABCimprove(b);
GlobalMins1=GlobalMins;
hold on;
semilogy(GlobalMins1,'k-');
legend('原始蜂群算法','本文算法');
