clear all;
close all;

dt = 0.01;
tf = 1;
t = 0:dt:tf;

N = length(t);

% Unit: m/s
j_max = 0.01;
a_max = 0.1;
v_max = 0.6;

j = zeros(1,length(t));
a = zeros(1,length(t));
v = zeros(1,length(t));
p = zeros(1,length(t));

j_sum = zeros(1,length(t));
a_sum = zeros(1,length(t));
v_sum = zeros(1,length(t));
p_sum = zeros(1,length(t));


t01 = sqrt(v_max/j_max);
t01=j_max/a_max;

t1 = t01;
t2 = 2*t01;

for i = 1:N-1
    
    if t(i)>=0 && t(i)<t1
        j(i) = j_max;
        a(i) = j_max*t(i);
        v(i) = 1/2*j_max*t(i)^2;
        p(i) = 1/6*j_max*t(i)^3;
        
        j_sum(i) = j_max;
        a_sum(i+1) = a_sum(i) + j_sum(i)*dt;
        v_sum(i+1) = v_sum(i) + a_sum(i)*dt;
        p_sum(i+1) = p_sum(i) + v_sum(i)*dt;
        
    elseif t(i)>=t1 && t(i)<t2
        j(i) = -j_max;
        a(i) = -j_max*(t(i)-t01)+j_max*t01;
        v(i) = -j_max/2*t(i)*(t(i)-2*t01)+j_max*t01*t(i)-j_max*t01^2;
        p(i) = -j_max/2*(1/3*t(i)^3-t01*t(i)^2)+j_max/2*t01*t(i)^2-j_max*t01^2*t(i)+1/3*j_max*t01^3;
        
        p_t2 = j_max*t01^3;
        
        j_sum(i) = -j_max;
        a_sum(i+1) = a_sum(i) + j_sum(i)*dt;
        v_sum(i+1) = v_sum(i) + a_sum(i)*dt;
        p_sum(i+1) = p_sum(i) + v_sum(i)*dt;
    elseif t(i)>=t2
        j(i) = 0;
        a(i) = 0;
        v(i) = v_max;
        p(i) = p_t2+v(i)*(t(i)-t2);
        
        
        j_sum(i) = 0;
        a_sum(i+1) = 0;
        v_sum(i+1) = v_max;
        p_sum(i+1) = p_sum(i) + v_sum(i)*dt;
        
    end
end

v(end) = v(N-1);
p(end) = p(N-1);

subplot(2,2,1)
plot(t,j);
hold on;
plot(t,j_sum);
title('Jerk')
legend('Numerical','Sum')

subplot(2,2,2)
plot(t,a);
hold on;
plot(t,a_sum);
title('Acc')
legend('Numerical','Sum')

subplot(2,2,3)
plot(t,v);
hold on;
plot(t,v_sum);
title('Vel')
legend('Numerical','Sum')

subplot(2,2,4)
plot(t,p);
hold on;
plot(t,p_sum);
title('Pos')
legend('Numerical','Sum')