clear all;
close all;

dt = 0.02;
tf = 300;
t = 0:dt:tf;

N = length(t);

% Unit: m/s
j_max = 0.1;
a_max = 0.1;
v_max = 0.6;

j = zeros(1,length(t));
a = zeros(1,length(t));
v = zeros(1,length(t));
p = zeros(1,length(t));
r = zeros(2,length(t));
psi = zeros(1,length(t));

j_sum = zeros(1,length(t));
a_sum = zeros(1,length(t));
v_sum = zeros(1,length(t));
p_sum = zeros(1,length(t));


t01 = a_max/j_max;
t02 = v_max/a_max - t01;

t1 = t01;
t2 = t01+t02;
t3 = t2+t01;

t_vel_max = 2*t01 + t02;
p_t3 = j_max/2*t01^2+(j_max*t01*t02+1/2*j_max*t01^2)*t01+j_max*t01/2*(t01+t02)*t02;


goal_distance = 500;


t_stop_start = t3 + (goal_distance - 2*p_t3)/v_max;
t4 = t_stop_start +t01;
t5 = t4+t02;
t6 = t5+t01;

R_c = 20;
psi_rate = v_max/R_c;

d1 = 20;
d2 = d1 + R_c*pi;
d3 = d2 + d1;
d4 = d3 + R_c*pi;

for i = 1:N
    
    if t(i)>=0 && t(i)<t1
        j(i) = j_max;
        a(i) = j_max*t(i);
        v(i) = 1/2*j_max*t(i)^2;
        p(i) = 1/6*j_max*t(i)^3;
        v_t1 = 1/2*j_max*t1^2;
        p_t1 = 1/6*j_max*t1^3;
        
        
        psi(i) = 0;
        r(1,i) = p(i);
        r(2,i) = 0;
        
        j_sum(i) = j_max;
        a_sum(i+1) = a_sum(i) + j_sum(i)*dt;
        v_sum(i+1) = v_sum(i) + a_sum(i)*dt;
        p_sum(i+1) = p_sum(i) + v_sum(i)*dt;
        
    elseif t(i)>=t1 && t(i)<t2
        j(i) = 0;
        a(i) = a_max;
        v(i) = a_max*(t(i)-t1) + v_t1;
        p(i) = a_max/2*(t(i)-t1)^2 + v_t1*(t(i)-t1) + p_t1;
        
        v_t2 = a_max*(t2-t1) + v_t1;
        p_t2 = a_max/2*(t2-t1)^2 + v_t1*(t2-t1) + p_t1;
        
        psi(i) = 0;
        r(1,i) = p(i);
        r(2,i) = 0;
        
        j_sum(i) = 0;
        a_sum(i+1) = a_sum(i) + j_sum(i)*dt;
        v_sum(i+1) = v_sum(i) + a_sum(i)*dt;
        p_sum(i+1) = p_sum(i) + v_sum(i)*dt;
        
    elseif t(i)>=t2 && t(i) <= t3
        j(i) = -j_max;
        a(i) = -j_max*(t(i)-t2) + a_max;
        v(i) = -j_max/2*(t(i)-t2)^2 + a_max*(t(i)-t2) + v_t2;
        p(i) = -j_max/6*(t(i)-t2)^3 + a_max/2*(t(i)-t2)^2 + v_t2*(t(i)-t2) + p_t2;
        
        v_t3 = -j_max/2*(t3-t2)^2 + a_max*(t3-t2) + v_t2;
        p_t3 = -j_max/6*(t3-t2)^3 + a_max/2*(t3-t2)^2 + v_t2*(t3-t2) + p_t2;
        
        psi(i) = 0;
        r(1,i) = p(i);
        r(2,i) = 0;
        
        j_sum(i) = -j_max;
        a_sum(i+1) = a_sum(i) + j_sum(i)*dt;
        v_sum(i+1) = v_sum(i) + a_sum(i)*dt;
        p_sum(i+1) = p_sum(i) + v_sum(i)*dt;
        
    elseif t(i) > t3
        j(i) = 0;
        a(i) = 0;
        v(i) = v_max;
        p(i) = p_t3+v_max*(t(i)-t3);
        
        p_vel = p_t3+v_max*(t_stop_start-t3);
        p_vel = goal_distance - p_t3;
        
        j_sum(i) = 0;
        a_sum(i+1) = a_sum(i) + j_sum(i)*dt;
        v_sum(i+1) = v_sum(i) + a_sum(i)*dt;
        p_sum(i+1) = p_sum(i) + v_sum(i)*dt;
        
        r(1,i) = p(i);
        r(2,i) = 0;
    end
    
    
    if t(i) >= t_stop_start && t(i) < t4
        j(i) = -j_max;
        a(i) = -j_max*(t(i)-t_stop_start);
        v(i) = v_max - 1/2*j_max*(t(i)-t_stop_start)^2;
        p(i) = p_vel + v_max*(t(i)-t_stop_start) - 1/6*j_max*(t(i)-t_stop_start)^3;
        
        v_t4 = v_max - 1/2*j_max*(t4-t_stop_start)^2;
        p_t4 = p_vel + v_max*t01 - 1/6*j_max*t01^3;
        
    elseif t(i) >= t4 && t(i) <= t5
        j(i) = 0;
        a(i) = -a_max;
        v(i) = -a_max*(t(i)-t4) + v_t4;
        p(i) = -a_max/2*(t(i)-t4)^2 + v_t4*(t(i)-t4) + p_t4;
        
        v_t5 = -a_max*(t5-t4) + v_t4;
        p_t5 = -a_max/2*(t5-t4)^2 + v_t4*(t5-t4) + p_t4;
        
    elseif t(i) >= t5 && t(i) <= t6
        j(i) = j_max;
        a(i) = j_max*(t(i)-t5) - a_max;
        v(i) = j_max/2*(t(i)-t5)^2 - a_max*(t(i)-t5) + v_t5;
        p(i) = j_max/6*(t(i)-t5)^3 - a_max/2*(t(i)-t5)^2 + v_t5*(t(i)-t5) + p_t5;
        
        p_t6 = j_max/6*(t6-t5)^3 - a_max/2*(t6-t5)^2 + v_t5*(t6-t5) + p_t5;
    elseif t(i) > t6
        j(i) = 0;
        a(i) = 0;
        v(i) = 0;
        p(i) = p_t6;
    end
    
    if p(i) >= d1 && p(i) <= d2
        t_ref = (20 - p_t3)/v_max + t3;
        psi(i) = psi_rate*(t(i)-t_ref);
        
        r(1,i) = v_max/psi_rate*sin(psi(i))+20;
        r(2,i) = -v_max/psi_rate*(cos(psi(i))-1);
        
    elseif p(i) >= d2 && p(i) < d3
        psi(i) = pi;
        r(1,i) = 20 - (p(i) - R_c*pi - 20);
        r(2,i) = R_c*2;
        
        r1_end = 20 - (d3 - R_c*pi -20);
        r2_end = R_c*2; 
        
    elseif p(i) >=  d3 && p(i) < d4
        t_ref2 = (d3 - p_t3)/v_max + t3;
        psi(i) = psi_rate*(t(i)-t_ref2) + pi;
        
        r(1,i) =  v_max/psi_rate*sin(psi(i));
        r(2,i) = -v_max/psi_rate*(cos(psi(i))-1);
        
    elseif p(i) >= d4
        t_ref3 = (d4-p_t3)/v_max + t3;
        psi(i) = 2*pi;
        
        r(1,i) = p(i) - d4;
        r(2,i) = 0;
    end
    
end




figure(1)
subplot(2,2,1)
plot(t,j,'r');
hold on;
title('Jerk')
xlabel('s')
ylabel('m/s^3')


subplot(2,2,2)
plot(t,a,'r');
hold on;
title('Acc')
xlabel('s')
ylabel('m/s^2')


subplot(2,2,3)
plot(t,v,'r');
hold on;
title('Vel')
xlabel('s')
ylabel('m/s')

subplot(2,2,4)
plot(t,p,'r');
hold on;
title('Pos')
xlabel('s')
ylabel('m') 


figure(2)
plot(r(1,:),r(2,:))
title('x - y')
xlabel('x (m)')
ylabel('y (m)')
figure(3)
subplot(2,2,1)
plot(t,r(1,:))
subplot(2,2,2)
plot(t,r(2,:))
subplot(2,2,3)
plot(t,psi)