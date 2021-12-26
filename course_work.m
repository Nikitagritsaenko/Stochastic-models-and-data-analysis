clear all
data_dir='C:\Users\Nikita\Desktop\bazhenov\kursach\'
close all
fn='wave_ampl.txt'
FN=strcat(data_dir, fn);
ADC_DATA=dlmread (FN, "," );

% data structure
sample_size=1024;
ch_count=8;
pages_count=10;
level_count=10;
wave_all=reshape(ADC_DATA,sample_size,ch_count,pages_count,level_count);
% /data structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg load interval
addpath('C:\Users\Nikita\Desktop\bazhenov\kursach')
addpath('C:\Users\Nikita\Desktop\bazhenov\kursach\m')
addpath('C:\Users\Nikita\Desktop\bazhenov\kursach\m\polytopes')

% фиксируем один канал
channel_number = 1

% Error model
data_for_level_roi = [];
data_for_level_bg = [];
win1=481;
win2=500;
win=[win1:win2];
winBG=[1:200]+1;
modes = []
for jj=1:10 % levels
  data_for_level_roi = [];
  data_for_level_bg = [];
  for ii=1:10 % runs
    datanow=wave_all(:,channel_number,ii,jj);
    bg_max = max(datanow(winBG));
    bg_min = min(datanow(winBG));
    bg_delta = mean(datanow(winBG));##(bg_max + bg_min)/2;
    bg_d = bg_delta*ones(length(win),1);

    roi = datanow(win)-bg_d;
    data_for_level_roi = [
      data_for_level_roi 
      roi
    ];
    data_for_level_bg= [
      data_for_level_bg 
      datanow(winBG)
    ];
  end

  ERR=std(data_for_level_bg)
  err_const=3*ERR
  estBG=max(data_for_level_bg)-min(data_for_level_bg)

  data_sampleStat=err_const*ones(length(data_for_level_roi),1);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % interval
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  DataStd = midrad(data_for_level_roi, data_sampleStat);
  [mode_single, mu_array, max_mu, mode_ind, c_array, C, multi]= modeIR2(DataStd)

  verbose = 1

  modes = [modes; [inf(mode_single) sup(mode_single)]]

  if verbose > 0
    figure
    x_1 =[]
    for ii=1:length(c_array)
      x_1 = [x_1 inf(c_array(ii))];
    endfor
    plot(x_1,mu_array, '-b', "linewidth", 3);
    xlabel('Data')
    ylabel('\mu')
    title('\mu array')
    
    mode_arr=((modes(jj,1)+modes(jj,2))/2)*ones(2000,1);
    figure
    errorbar (1:length(data_for_level_roi), data_for_level_roi, data_sampleStat,".b");
    xmargin=0.5
    ymargin=err_const/2
    xlim( [1-xmargin length(data_for_level_roi)+xmargin] )
    ylim( [min(inf(DataStd))-ymargin max(sup(DataStd))+ymargin])
    hold on
    plot([1:2000],mode_arr, '-r', "linewidth", 10) 
    title_str=strcat('\it ADCDATA')
    title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
    xlabel('\it channel number');
    ylabel('\it Data');
    set(gca, 'fontsize', 12);
    box off;
    figure_name_out=strcat('ADCData2', '.png')
    print('-dpng', '-r300', figure_name_out), pwd
  end

end

mode_mids = [];
mode_rads = [];
for i = 1: level_count
  mode_mids(i) = (modes(i,2)+modes(i,1))/2;
  mode_rads(i) = (modes(i,2)-modes(i,1))/2;
endfor

figure
xxx = ([[-0.5:0.1:-0.1] [0.1:0.1:0.5]])'
errorbar (xxx, mode_mids, mode_rads,".b");
title_str=strcat('\it Modes for 1st channel')
title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
xlabel('\it V');
ylabel('\it mode');

#ЗЛП
x = ([[-0.5:0.1:-0.1] [0.1:0.1:0.5]])'
y = mode_mids'
eps = mode_rads'
m = size(x)(1)
C = zeros(1, m + 2);
for i = 1:m
  C(i) = 1;
end
A = zeros(2*m, m+2);

for i = 1:m
  A(2 * i - 1, i) = eps(i);
  A(2 * i, i) = eps(i);
  
  A(2 * i - 1, m + 1) = 1;
  A(2 * i, m + 1) = -1;
  
  A(2 * i - 1, m + 2) = x(i);
  A(2 * i, m + 2) = -x(i);
  
end

B = zeros(1, 2*m);
for i = 1:m
  B(2 * i - 1) = y(i);
  B(2 * i) = -y(i);
end

lb = zeros(1, m+2);
for i = 1:m
  lb(i) = 1;
end

lb(m+2) = -inf;

ctype = "";
for i = 1:2 * m
  ctype(i) = 'L';
end

vartype = "";
for i = 1:m + 2
  vartype(i) = 'C';
end

sense = 1

w = glpk(C,A,B,lb,[],ctype,vartype,sense)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(w)-2
scale = max(w(1:n))
for i = 1:n
    eps(i) = eps(i) * scale;
end

X = [ x.^0 x ];                             
lb = [-inf 0];                               
irp_temp = ir_problem(X, y, eps, lb);   
                                               
## График интервальных измерений
figure;
ir_scatter(irp_temp);   
set(gca, 'fontsize', 12)

beta1 = w(n+1)
beta2 = w(n+2)
xx = -0.5:0.01:0.5
yy = beta1 + xx * beta2
hold on
plot(xx, yy, 'r--')
title_str=strcat('\it Modes for 1st channel')
title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
xlabel('\it V');
ylabel('\it mode');

## Графическое представление информационного множества
figure;
ir_plotbeta(irp_temp)
grid on
set(gca, 'fontsize', 12)
xlabel('\beta_1')
ylabel('\beta_2')
title('Information set')

## Вершины информационного множества задачи построения интервальной регрессии
vertices = ir_beta2poly(irp_temp)

## Диаметр и наиболее удаленные вершины информационного множества 
[rhoB, b1, b2] = ir_betadiam(irp_temp)

## Внешние интервальние оценки параметров модели y = beta1 + beta2 * x 
b_int = ir_outer(irp_temp)

## Точечные оценки параметров 
b_maxdiag = (b1 + b2) / 2    # как середина наибольшей диагонали информационного множества

b_gravity = mean(vertices)   # как центр тяжести информационного множества 

b_lsm = (X \ y)'             # методом наименьших квадратов


## Графическое представление внешней интервальной оценки информационного множества
figure('position',[0, 0, 800, 600]);
ir_plotbeta(irp_temp)
hold on
ir_plotrect(b_int,'r-')
grid on
set(gca, 'fontsize', 12)
xlabel('\beta_1')
ylabel('\beta_2')
title('Information set')

## Точечные оценки
plot(b_maxdiag(1), b_maxdiag(2), 'ro')
plot(b_gravity(1), b_gravity(2), 'k+')
plot(b_lsm(1), b_lsm(2), 'gx')
legend("", "", "enclosure", "maxdiag",  "gravity", "lsm")

## Графическое представление коридора совместных зависимостей для модели y = beta1 + beta2 * x
figure;
xlimits = [-1 1];
ir_plotmodelset(irp_temp, xlimits)     # коридор совместных зависимостей

hold on
ir_scatter(irp_temp,'bo')              # интервальные измерения
ir_plotline(b_maxdiag, xlimits, 'r-')   # зависимость с параметрами, оцененными как центр наибольшей диагонали ИМ
ir_plotline(b_gravity, xlim, 'b--')     # зависимость с параметрами, оцененными как центр тяжести ИМ  
ir_plotline(b_lsm, xlim, 'b--')         # зависимость с параметрами, оцененными МНК


grid on
set(gca, 'fontsize', 12)

#[a, b] +- [c, d]
x = 5
a = -0.00538
b =  0.00036
c = 0.9091 * x
d = 0.9246 * x
res = [a+c, b+d]
(res(2) - res(1)) / 2


