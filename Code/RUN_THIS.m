%RUN THIS FILE FOR ALL THE PLOTS
%NOTE: Run this file only when all the other dependent files exist in the
%same directory along with this file.

% 1D Case
%Plot 1 -- Phase Plane 1
[u, v] = etdrk4_1D(8000, 1.025, 0.2, 0.05, 0.35, 0.35);
x = u(25, :);
y = v(25, :);

figure;
plot(x,y);
xlabel('Prey Density (u)');
ylabel('Predator Density (v)');
title("Phase Plane");

%Plot 2 -- Oscillations 1
[u, v] = etdrk4_1D(1000, 1, 0.2, 0.05, 0.35, 0.35);
t = 0.1:0.1:100;
x = u(25,:);
y = v(25,:);

figure;
hold on;
plot(t, x);
plot(t, y);
xlabel('t');
ylabel('predator/prey density');
title("Predator/Prey Density Oscillations 1");
legend('Prey', 'Predator');
hold off;

%Plot 3 -- Phase Plane 2
[u, v] = etdrk4_1D(7000, 1.5, 0.01, 0.08, 1.5, 0.1);
x = u(25, :);
y = v(25, :);
t = 0.1:0.1:400;

figure;
plot(x,y);
xlabel('Prey Density (u)');
ylabel('Predator Density (v)');
title("Phase Plane");

%Plot 4 -- Oscillations 2
figure;
hold on;
plot(t, x(1:4000));
plot(t, y(1:4000));
xlabel('t');
ylabel('predator/prey density');
title("Predator/Prey Density Oscillations 2");
legend('Prey', 'Predator');
hold off;


%NOTE: Run the 2D Case code only if you have alot of time. The code takes
%too much time to run.

% 2D Case
%[u, v] = etdrk4_2D(5000, 0.2, 0.5, 2, 0.171429, 0.473469);
%x = 0.01:0.01:0.01*50;
%y = 0.01:0.01:0.01*50;

%Plot 5 -- Prey Distribution at t = 20s
%figure;
%imagesc(x, y, transpose(squeeze(u(:,:2000))));
%axis equal tight;
%colorbar;
%xlabel('x');
%ylabel('y');
%title('Prey Distribution at t = 20s');

%Plot 6 -- Predator Distribution at t = 20s
%figure;
%imagesc(x, y, transpose(squeeze(v(:,:,2000))));
%axis equal tight;
%colorbar;
%xlabel('x');
%ylabel('y');
%title('Predator Distribution at t = 20s');

%Plot 7 -- Prey Distribution at t = 50s
%figure;
%imagesc(x, y, transpose(squeeze(u(:,:,5000))));
%axis equal tight;
%colorbar;
%xlabel('x');
%ylabel('y');
%title('Prey Distribution at t = 50s');

%Plot 8 -- Predator Distribution at t = 50s
%figure;
%imagesc(x, y, transpose(squeeze(v(:,:,5000))));
%axis equal tight;
%colorbar;
%xlabel('x');
%ylabel('y');
%title('Predator Distribution at t = 50s');
