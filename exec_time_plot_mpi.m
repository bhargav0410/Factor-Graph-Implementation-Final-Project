clc;
clear all;
close all;

exec_time_rate = [1.21034	0.0020419	0.193653
10.1279	0.00109751	0.226749
13.9003	0.000739846	0.164823
20.7505	0.000592248	0.121223
38.2415	0.00066665	0.104314
44.943	0.000773664	0.100951
		
0.865256	0.00272235	0.127035
5.60817	0.00154382	0.154927
8.13247	0.000839983	0.114
11.6899	0.00122964	0.0854854
17.813	0.000724259	0.0693101
24.4676	0.00101188	0.0785712
		
0.350628	0.00229502	0.0645342
1.6148	0.00132082	0.0660303
4.11123	0.00067701	0.0607128
5.68162	0.000598624	0.0485018
7.52898	0.0012284	0.0430561
10.1986	0.00156201	0.0494491
];

exec_time_size = [1.21034	0.0020419	0.193653
10.1279	0.00109751	0.226749
13.9003	0.000739846	0.164823
20.7505	0.000592248	0.121223
38.2415	0.00066665	0.104314
44.943	0.000773664	0.100951
		
0.865256	0.00272235	0.127035
5.60817	0.00154382	0.154927
8.13247	0.000839983	0.114
11.6899	0.000805	0.0854854
17.813	0.000724259	0.0693101
24.4676	0.00101188	0.0785712
		
0.350628	0.00229502	0.0645342
1.6148	0.00132082	0.0660303
4.11123	0.00067701	0.0607128
5.68162	0.000598624	0.0485018
7.52898	0.0012284	0.0430561
10.1986	0.00156201	0.0494491
];

rates = [0.25;0.5;0.75];
procs = [1;2;4;8;16;32];
figure; hold on;
plot(procs, exec_time_rate(1:6,1), 'b-*');
plot(procs, exec_time_rate(7:12,1), 'r-<');
plot(procs, exec_time_rate(13:18,1), 'k-o');
legend('rate=0.25','rate=0.5','rate=0.75');
xlabel('Number of processor(s)');
ylabel('Execution time (secs)');
title('Construction of G matrix');

%Encoding
figure; hold on;
plot(procs, exec_time_rate(1:6,2), 'b-*');
plot(procs, exec_time_rate(7:12,2), 'r-<');
plot(procs, exec_time_rate(13:18,2), 'k-o');
legend('rate=0.25','rate=0.5','rate=0.75');
xlabel('Number of processor(s)');
ylabel('Execution time (secs)');
title('Encoding');

%Decoding
figure; hold on;
plot(procs, exec_time_rate(1:6,3), 'b-*');
plot(procs, exec_time_rate(7:12,3), 'r-<');
plot(procs, exec_time_rate(13:18,3), 'k-o');
legend('rate=0.25','rate=0.5','rate=0.75');
xlabel('Number of processor(s)');
ylabel('Execution time (secs)');
title('Decoding');

%plots for matrix size
figure; hold on;
plot(procs, exec_time_size(1:6,3), 'b-*');
plot(procs, exec_time_size(7:12,3), 'r-<');
plot(procs, exec_time_size(13:18,3), 'k-o');
legend('size=128','size=512','size=1024');
xlabel('Number of processor(s)');
ylabel('Execution time (secs)');
title('Decoding');

%construction
figure; hold on;
plot(procs, exec_time_size(1:6,1), 'b-*');
plot(procs, exec_time_size(7:12,1), 'r-<');
plot(procs, exec_time_size(13:18,1), 'k-o');
legend('size=128','size=512','size=1024');
xlabel('Number of processor(s)');
ylabel('Execution time (secs)');
title('Construction of G mat');

%encoding
figure; hold on;
plot(procs, exec_time_size(1:6,2), 'b-*');
plot(procs, exec_time_size(7:12,2), 'r-<');
plot(procs, exec_time_size(13:18,2), 'k-o');
legend('size=128','size=512','size=1024');
xlabel('Number of processor(s)');
ylabel('Execution time (secs)');
title('Encoding');