clc;
clear all;
close all;

exec_time_mat_size = [0.128889	0.000563307	0.267241	2.99	0.00062471	0.0174659
0.932039	0.00197877	0.537234	6.45094	0.000778633	0.062739
7.29287	0.00760943	1.0781	20.7091	0.00106778	0.423505
57.2809	0.0299144	2.16793	72.9526	0.00148155	2.79163
];

exec_time_rate = [16.3385	0.0307575	0.614106	21.1335	0.00479999	0.694026
44.6883	0.0412871	1.34164	51.2396	0.00289758	1.62222
57.2809	0.0299144	2.16793	72.9526	0.00148155	2.79163
];

exec_time_v_nodes = [57.2809	0.0299144	2.16793	72.9526	0.00148155	2.79163
35.9397	0.0377901	3.75954	42.6659	0.00388177	1.23017
24.0699	0.0302463	5.27402	29.8612	0.00461961	0.751254
17.5419	0.0240446	6.95443	23.2273	0.00519864	0.606409
];

rates = [0.25,0.5,0.75];

figure;
plot([128;256;512;1024], exec_time_mat_size(:,1), 'b--'); hold on;
plot([128;256;512;1024], exec_time_mat_size(:,4), 'r-');
legend('CPU','GPU');
xlabel('Number of columns of H matrix');
ylabel('Execution time (secs)');
title('Construction of G mat');

figure;
plot(rates, exec_time_rate(:,1), 'b--'); hold on;
plot(rates, exec_time_rate(:,4), 'r-');
legend('CPU','GPU');
xlabel('Rate');
ylabel('Execution time (secs)');
title('Construction of G mat');

figure;
plot([128;256;512;1024], exec_time_mat_size(:,2), 'b--'); hold on;
plot([128;256;512;1024], exec_time_mat_size(:,5), 'r-');
legend('CPU','GPU');
xlabel('Number of columns of H matrix');
ylabel('Execution time (secs)');
title('Encoding');

figure;
plot([128;256;512;1024], exec_time_mat_size(:,3), 'b--'); hold on;
plot([128;256;512;1024], exec_time_mat_size(:,6), 'r-');
legend('CPU','GPU');
xlabel('Number of columns of H matrix');
ylabel('Execution time (secs)');
title('Decoding');

figure;
plot(rates, exec_time_rate(:,2), 'b--'); hold on;
plot(rates, exec_time_rate(:,5), 'r-');
legend('CPU','GPU');
xlabel('Rate');
ylabel('Execution time (secs)');
title('Encoding');

figure;
plot(rates, exec_time_rate(:,3), 'b--'); hold on;
plot(rates, exec_time_rate(:,6), 'r-');
legend('CPU','GPU');
xlabel('Rate');
ylabel('Execution time (secs)');
title('Decoding');