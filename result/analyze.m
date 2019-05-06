% 首先导入txt数据，示例用output.txt; 再将导入的元胞数据改名为ouput
len = length(output.time_stamp);

% 打印实际位置、 滤波后位置
figure(1);
plot(output.px_state,output.py_state);
hold on;
plot(output.px_ground_truth,output.py_ground_truth);
hold on;

% 打印NIS
NIS_lidar = [];
NIS_radar = [];
len_lidar = 0;
len_radar = 0;
for i= 1:len
    if grp2idx(output.sensor_type(i)) == 1
        len_lidar = len_lidar + 1;
        NIS_lidar(len_lidar) = output.NIS(i);      
    else
        len_radar = len_radar + 1;
        NIS_radar(len_radar) = output.NIS(i); 
    end
end
figure(2);
plot(NIS_lidar);
hold on;
plot(5.99.* ones(len_lidar,1));
hold on;

figure(3);
plot(NIS_radar);
hold on;
plot(7.82.* ones(len_radar,1));
hold on;
