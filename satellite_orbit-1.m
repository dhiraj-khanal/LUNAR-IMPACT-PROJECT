function satellite_orbit

%SATELLITE
%
% launched from position the perigee at t = 0
GM = 3.986012*10^5;
time = 27.7*24*60*60*2;
rpSate = 250+6378.145;
raSate = 35950+6378.145;
theta = 6*pi/180;
O = 200*pi/180;
W = 178*pi/180;

%defining components of initial position (t = 0) for perigee
rpi= rpSate*(cos(O) * cos(W)-sin(O)*sin(W)*cos(theta));
rpj = rpSate*(cos(W) * sin(O)+sin(W)*cos(O)*cos(theta));
rpk = rpSate*(sin(W)*sin(theta));
%i.c for velocity for perigee
c = sqrt((2*GM*raSate/(rpSate+raSate)/rpSate));
vpi = -c*(cos(O) * sin(W)+cos(W)*sin(O)*cos(theta));
vpj = c*(cos(O) * cos(W)*cos(theta)-sin(O)*sin(W));
vpk = c*cos(W)*sin(theta);

%apogee
 %i.c for position/velocity for perigee
rai=-rpi/rpSate*raSate;
raj=-rpj/rpSate*raSate;
rak=-rpk/rpSate*raSate;
ca = sqrt((2*GM*rpSate/(rpSate+raSate)/raSate));
vai = -vpi/c*ca;
vaj = -vpj/c*ca;
vak = -vpk/c*ca;

w0 = [rpi,rpj,rpk,vpi,vpj,vpk]; % Initial conditions for perigee
w1 = [rai,raj,rak,vai,vaj,vak]; % Initial conditions for apogee

optionsSate = odeset('RelTol',0.000001);
[t_values,w_values] = ode45(@(t,w) odefunc(t,w,GM),[0,time],w0, optionsSate);
%plot3(w_values(:,1),w_values(:,2),w_values(:,3),'LineWidth',2,'LineStyle',...
 %   ':','Color',[0 0 1]);

%plot of energy and angular momemtum for satellite

energyS = t_values;
angularmS = t_values;
for i =1:length(t_values)
     rS = sqrt(w_values(i,1)^2 + w_values(i,2)^2+w_values(i,3)^2);
     vmagS = sqrt(w_values(i,4)^2 + w_values(i,5)^2+w_values(i,6)^2);
%TotalE/mass
energyS(i) = -GM/rS + vmagS^2/2;
%mag(TotalH)/mass
angularmS(i) = sqrt((w_values(i,2)*w_values(i,6)-w_values(i,3)*w_values(i,5))^2 ...
    + (w_values(i,1)*w_values(i,6)-w_values(i,3)*w_values(i,4))^2 +...
    (w_values(i,1)*w_values(i,5)-w_values(i,2)*w_values(i,4))^2);
end

%the plot is looks APPROXIMATELY constant
plot(t_values,energyS)
  title('Energy Vs time');
  ylabel('E/m');
  xlabel('Time');
figure
plot(t_values, angularmS)
  title('AngularMomentum Vs time');
  ylabel('Angular Momumtem/Mass');
  xlabel('Time');

%MOON
%
%
%Initial conditions for moon
%date chosen april 7th, 2020
rpMoon = 356906+6378.145;
raMoon = 406692+6378.145;
alpha = 183.325*pi/180; %2h13m18s
delta = 4.0283*pi/180; %+04°01'42"
rMoon = rpMoon*[cos(alpha)*cos(delta), cos(delta)*sin(alpha), sin(delta)];
%mid date april 22
deltaMid = 3.387777*pi/180; %+03°23'16"
alphaMid = 20.2458*pi/180; %01h20m59s
nMid = [cos(alphaMid)*cos(deltaMid), cos(deltaMid)*sin(alphaMid), sin(deltaMid)];
m = cross(nMid,rMoon)/norm(cross(nMid,rMoon));

Vp = sqrt(2*GM*raMoon/rpMoon/(raMoon+rpMoon))*(cross(m,rMoon)/rpMoon);
winc = [rMoon(1),rMoon(2),rMoon(3),Vp(1),Vp(2),Vp(3)]; % Initial conditions

optionsMoon = odeset('Event',@(t,w) detect_impact_point(t,w, O, theta));
[t_vals,sol_vals,t_event,sol_event,index_event] = ode45(@(t,w) odefunc(t,w,GM),[0,time],winc, [optionsMoon, optionsSate]);
distance_between_AP = norm([sol_event(1,1:3)]-[sol_event(2,1:3)])
figure
%plot of moon's trajectory
plot3(sol_vals(:,1),sol_vals(:,2),sol_vals(:,3));
hold on;
%this is time period of the moon was found to be 27.6236 days
time_period_of_moon = (t_event(3) - t_event(1))/86400
time_for_impact_moon= t_event(3)
rmoon1 = sol_event(1,1:3);
plot3(rmoon1(1),rmoon1(2),rmoon1(3),'ro','MarkerSize',11,'MarkerFaceColor','r')
hold on;
rmoon2 = sol_event(2,1:3);
plot3(rmoon2(1),rmoon2(2),rmoon2(3),'ro','MarkerSize',11,'MarkerFaceColor','g')
hold on;

for i = 1:100 
 delta_Vp(i) = 0.67 + i/800;
 [distp(i),traveltimep(i)] = mindisttomoon(delta_Vp(i),rmoon1,[rpi,rpj,rpk],[vpi,vpj,vpk],GM);
end
for i = 1:100 
 delta_Va(i) = 2.5 + i/1000; 
[dista(i),traveltimea(i)] = mindisttomoon(delta_Va(i),rmoon2,[rai,raj,rak],[vai,vaj,vak],GM);
end

figure
plot(distp,delta_Vp)
  title('Perigee');
  ylabel('Additional velocity given');
  xlabel('Distance to the moon');
[distmin_perigee,index_perigee] = min(distp); % Find the min value in the vector dist
critical_delta_V_perigee = delta_Vp(index_perigee) %Print the delta V
time_to_hit_from_perigee = traveltimep(index_perigee) %Print the travel time 
 
figure
plot(dista,delta_Va)
  title('Apogee');
  ylabel('Additional velocity given');
  xlabel('Distance to the moon');
[distmin_apogee,index_apogee] = min(dista); % Find the min value in the vector dist
critical_delta_V_apogee = delta_Va(index_apogee) %Print the delta V
time_to_hit_from_apogee = traveltimea(index_apogee) %Print the travel time 

time_lun = time_for_impact_moon - time_to_hit_from_perigee;
[time_launch, moonVar] = ode45(@(t,w) odefunc(t,w,GM),[0,time_lun],winc, [optionsMoon, optionsSate]);
moonVar(end)
ic = [rpi,rpj,rpk,vpi+vpi/c*critical_delta_V_perigee,vpj+vpj/c*critical_delta_V_perigee,vpk+vpk/c*critical_delta_V_perigee...
    ,moonVar(end,1),moonVar(end,2),moonVar(end,3),moonVar(end,4),moonVar(end,5),moonVar(end,6)];
option = odeset('RelTol',0.000001,'Events',@(t,w) moonColl(t,w));
[tSim,simSol] = ode45(@(t,w) eom(t,w,GM), [0, 27*24*60*60], ic, option);
vel_impact = norm([simSol(end,4), simSol(end,5), simSol(end,6)]- [simSol(end,10), simSol(end,11), simSol(end,12)])
animate_trajectories(simSol(:,7:9), simSol(:,1:3))
end

function [d,time_to_reach_moon] = mindisttomoon(delta_V,r_impact,r0_sat,v0_sat,GM)
% Function to calculate min distance to moon for satellite
% trajectory. r0 and v0 are the position and velocity of the
% satellite when the rocket is fired, and delta_V is the increase in
% its velocity
%
 time = 8*24*60*60.; % 8 days is enough to reach the moon
 y0(1:3) = r0_sat; % Initial position of satellite
 y0(4:6) = v0_sat + delta_V*v0_sat/norm(v0_sat);
 rtol = 0.00000001;
 options = odeset('RelTol',rtol,'Events',@(t,w) min_dist(t,w,r_impact));
% Run the simulation to calculate the trajectory of the satellite
 [times,sols,tevent,sevent,ievent] = ode45(@(t,w) odefunc(t,w,GM),[0,time],y0,options);
 plot3(sols(:,1),sols(:,2),sols(:,3)) % Plot the trajectory
% Find distance to impact point at the event
 rmin = sevent(1,1:3)-r_impact;
 time_to_reach_moon = tevent(1);
 d = norm(rmin); % The min distance
end

function [event_val,stopthecalc,direction] = min_dist(t,w,rImp)
r = w(1:3);
v = w(4:6);
event_val = dot(v,(r-transpose(rImp)));
stopthecalc = 0;
direction = 0;
end

function [event_val,stopthecalc,direction] = detect_impact_point(t,w,O,theta)
 % Formula for unit vector normal to satellite orbit plane
 n = [sin(O)*sin(theta), -cos(O)*sin(theta), cos(theta)];
 % Position vector of moon (this assumes w(1)=x,w(2)=y,w(3)=z)
 r = w(1:3);
 % Detect when r.n=0
 event_val = dot(n,r);
 stopthecalc = 0;
 direction = 0;
end
function [event_val,stopthecalc,direction] = moonColl(t,w)
radius_moon = 1731.1;
rmoon = w(1:3);
rsat = w(7:9);
event_val = norm(rmoon-rsat);
stopthecalc = radius_moon;
direction = -1;
end
 
function dwdt = odefunc(t,w,GM)
   x=w(1); y=w(2); z = w(3); 
   vx=w(4); vy=w(5); vz = w(6);
   r = sqrt(x^2+y^2+z^2);
   dwdt = [vx;vy;vz;-GM*x/r^3;-GM*y/r^3;-GM*z/r^3];
end

function dwdt = eom(t,w,GM)
   x=w(1); y=w(2); z = w(3); 
   vx=w(4); vy=w(5); vz = w(6);
   x1=w(7); y1=w(8); z1 = w(9); 
   vx1=w(10); vy1=w(11); vz1 = w(12);
   rs = sqrt(x^2+y^2+z^2);
   rm = sqrt(x1^2 +y1^2+z1^2);
   dis = sqrt((x-x1)^2+(y-y1)^2+(z-z1)^2);
   dwdt = [vx;vy;vz;-GM*x/rs^3-GM*(abs(x-x1))*0.012298/dis^3;-GM*y/rs^3-GM...
       *(abs(y-y1))*0.012298/dis^3;-GM*z/rs^3-GM*(abs(z-z1))*0.012298...
       /dis^3;vx1;vy1;vz1;-GM*x1/rm^3;-GM*y1/rm^3;-GM*z1/rm^3];
end

function animate_trajectories(moon_pos,sat_pos)
% moon_pos - matrix of position values for moon
% sat_pos - matrix of position values for satellite
 [npoints,~] = size(moon_pos);
 figure
 for i=1:npoints
 clf;
 plot3(moon_pos(:,1),moon_pos(:,2),moon_pos(:,3)); % Plot the full path
 hold on
 plot3(sat_pos(:,1),sat_pos(:,2),sat_pos(:,3)); % Sat trajectory
 plot3(moon_pos(i,1),moon_pos(i,2),moon_pos(i,3),...
 'ko','MarkerFaceColor','y','MarkerSize',25);
 plot3(sat_pos(i,1),sat_pos(i,2),sat_pos(i,3),...
 'ko','MarkerFaceColor','r','MarkerSize',8);
 pause(0.1);
 end
end
