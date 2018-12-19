clc;
clear all;
close all;

exec_time_rate = [1.21034	0.0020419	0.716877
10.1279	0.00109751	1.682831
13.9003	0.000739846	1.344659
20.7505	0.000592248	1.016417
38.2415	0.00066665	0.920012
44.943	0.000773664	0.759943
		
0.865256	0.00272235	0.452484
5.60817	0.00154382	0.880731
8.13247	0.000839983	0.900108
11.6899	0.000805	0.718595
17.813	0.000724259	0.592604
24.4676	0.00101188	0.521696
		
0.350628	0.00229502	0.218838
1.6148	0.00132082	0.593197
4.11123	0.00067701	0.481391
5.68162	0.000598624	0.364992
7.52898	0.0012284	0.305823
10.1986	0.00156201	0.271963
];

exec_time_size = [0.00334924	4.98E-05	0.082966	0.010779	27.773484
0.0473827	5.81E-05	0.171724	0.006237	13.993944
0.0582511	5.32E-05	0.154365	0.003583	7.398815
0.117379	0.000091	0.113423	0.002373	4.090514
0.128736	0.000607844	0.11387	0.002838	2.462502
0.200442	0.000792	0.094001	0.002349	1.604671
				
0.18923	0.000504496	0.348933	0.163467	115.578904
1.6857	0.000295396	0.673532	0.082957	57.757938
1.80153	0.000184531	0.651806	0.043015	30.707095
3.44396	0.000187695	0.527515	0.02558	16.400642
3.82004	0.000681857	0.416307	0.017743	9.063382
6.11223	0.000926438	0.389791	0.017295	6.678328
				
1.21034	0.0020419	0.716877	0.656252	233.13479
10.1279	0.00109751	1.682831	0.328925	119.174823
13.9003	0.000739846	1.344659	0.18026	64.143504
20.7505	0.000592248	1.016417	0.09623	33.268712
38.2415	0.00066665	0.920012	0.085881	25.402884
44.943	0.000773664	0.759943	0.053106	16.854669
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


figure; hold on;
plot(procs, exec_time_size(1:6,4), 'b-*');
plot(procs, exec_time_size(7:12,4), 'r-<');
plot(procs, exec_time_size(13:18,4), 'k-o');
legend('size=128','size=512','size=1024');
xlabel('Number of processor(s)');
ylabel('Execution time (secs)');
title('Encoding (block)');

figure; hold on;
plot(procs, exec_time_size(1:6,5), 'b-*');
plot(procs, exec_time_size(7:12,5), 'r-<');
plot(procs, exec_time_size(13:18,5), 'k-o');
legend('size=128','size=512','size=1024');
xlabel('Number of processor(s)');
ylabel('Execution time (secs)');
title('Decoding (block)');