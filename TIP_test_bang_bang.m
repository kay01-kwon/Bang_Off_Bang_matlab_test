clear all;
close all;

dt = 0.01;
tf = 100;
t = 0:dt:tf;

N = length(t);

% Unit: m/s
j_max = 0.001;
a_max = 0.01;
v_max = 0.05;

j = zeros(1,length(t));
a = zeros(1,length(t));
v = zeros(1,length(t));
p = zeros(1,length(t));

j_sum = zeros(1,length(t));
a_sum = zeros(1,length(t));
v_sum = zeros(1,length(t));
p_sum = zeros(1,length(t));

goal_distance = 1;

t01 = a_max/j_max;
t02 = v_max/a_max -t01;

if t02 <= 0
    t01 = sqrt(v_max/j_max);
    t1 = t01;
    t02 = t1;
    t2 = 2*t01;
end

a_t1 = j_max*t1;
v_t1 = 1/2*j_max*t1^2;
p_t1 = 1/6*j_max*t1^3;

v_t2 = -j_max*t01 + a_t1*t01;
p_t2 = -j_max/6*t1^3 + a_t1/2*t01^2 + v_t1*t01 + p_t1;

p_stop_start = goal_distance - p_t2;

a_t3 = -j_max*t01;
v_t3 = -j_max/2*t01^2 + v_max;
p_t3 = -j_max/6*t01^3 + v_max*t01 + p_stop_start;


t_stop_start = (goal_distance - 2*p_t2)/v_max + t2;
t3 = t_stop_start + t01;
t4 = t3 + t01;

for i = 1:N
    
    if t(i)>=0 && t(i)<t1
        j(i) = j_max;
        a(i) = j_max*t(i);
        v(i) = 1/2*j_max*t(i)^2;
        p(i) = 1/6*j_max*t(i)^3;
        
    elseif t(i)>=t1 && t(i) < t2
        j(i) = -j_max;
        a(i) = -j_max*(t(i)-t1) + a_t1;
        v(i) = -j_max/2*(t(i)-t1)^2 + a_t1*(t(i)-t1) + v_t1;
        p(i) = -j_max/6*(t(i)-t1)^3 + a_t1/2*(t(i)-t1)^2 + v_t1*(t(i)-t1) + p_t1;
        
        j_sum(i) = -j_max;
        a_sum(i+1) = a_sum(i) + j_sum(i)*dt;
        v_sum(i+1) = v_sum(i) + a_sum(i)*dt;
        p_sum(i+1) = p_sum(i) + v_sum(i)*dt;
        
    elseif t(i)>=t2
        j(i) = 0;
        a(i) = 0;
        v(i) = v_max;
        p(i) = p_t2 + v_max*(t(i)-t2);
        
        j_sum(i) = 0;
        a_sum(i+1) = 0;
        v_sum(i+1) = v_max;
        p_sum(i+1) = p_sum(i) + v_sum(i)*dt;
    end
    
    if t(i) >= t_stop_start && t(i) < t3
        j(i) = -j_max;
        a(i) = -j_max*(t(i)-t_stop_start);
        v(i) = -j_max/2*(t(i)-t_stop_start)^2 + v_max;
        p(i) = -j_max/6*(t(i)-t_stop_start)^3 + v_max*(t(i)-t_stop_start) + p_stop_start;
        
    elseif t(i) >= t3 && t(i) < t4
        j(i) = j_max;
        a(i) = j_max*(t(i)-t3) + a_t3;
        v(i) = j_max/2*(t(i)-t3)^2 + a_t3*(t(i)-t3) + v_t3;
        p(i) = j_max/6*(t(i)-t3)^3 + a_t3/2*(t(i)-t3)^2 + v_t3*(t(i)-t3) + p_t3;
        
    elseif t(i) >=t4
        j(i) = 0;
        a(i) = 0;
        v(i) = 0;
        p(i) = goal_distance;
    end
end


subplot(2,2,1)
plot(t,j);
title('Jerk')

subplot(2,2,2)
plot(t,a);
title('Acc')

subplot(2,2,3)
plot(t,v);
title('Vel')

subplot(2,2,4)
plot(t,p);
title('Pos')