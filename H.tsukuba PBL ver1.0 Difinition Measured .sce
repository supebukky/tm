///////////////////////////////////////////////////////////////////////////
//本モデルはTERCにおける高度29.5mにおける顕熱フラックスについて解析したものである。//
///////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//（スキーム読み込み）//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///exec('H.tukuba(vol.2.3 for Semi) [compare of Measured] .sce')
clear//前回までの計算のクリア
tic()//時間計測開始
//A_2009=excel2sci('C:\SC\use scilab clucurate data\tsukuba data 2009 Jan-Apr 1.6m 29.5m new1.0.csv');//ファイル取り込み(筑波公開データ)
//B_2009_Jan_Apr_TERC=evstr(A_2009);//データベクトル化(筑波公開データ)//2009年（9月-12月）データ
//save('C:\SC\use scilab clucurate data\binary\B_2009_Jan_Apr_TERC.sce',B_2009_Jan_Apr_TERC);
load('C:\SC\use scilab clucurate data\binary\B_2009_Jan_Apr_TERC.sce');
load('C:\SC\use scilab clucurate data\binary\B_2008_TERC.sce');

//A_2009=excel2sci('C:\SC\use scilab clucurate data\tsukuba data 2009 Jan-Apr 1.6m 29.5m temperature.csv');//ファイル取り込み(筑波公開データ)
//B_2009_Jan_Apr_TERC_temp=evstr(A_2009);//データベクトル化(筑波公開データ)//2009年（9月-12月）データ
//save('C:\SC\use scilab clucurate data\binary\B_2009_Jan_Apr_TERC_temp.sce',B_2009_Jan_Apr_TERC_temp);
//load('C:\SC\use scilab clucurate data\binary\B_2009_Jan_Apr_TERC_temp.sce');

//A_2008=excel2sci('C:\SC\use scilab clucurate data\tsukuba data 2008 1.6m 12.3m 29.5m temperature ver1.1 .csv');//ファイル取り込み(筑波公開データ)
//B_2008_TERC_temp=evstr(A_2008);//データベクトル化(筑波公開データ)//2009年（9月-12月）データ
//save('C:\SC\use scilab clucurate data\binary\B_2008_TERC_temp.sce',B_2008_TERC_temp);
load('C:\SC\use scilab clucurate data\binary\B_2008_TERC_temp.sce');
load('C:\SC\use scilab clucurate data\binary\B_2009_Jan_Apr_TERC_temp.sce');

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//（ベクトル成分読み込み）///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////（ベクトル成分読み込み:for 2005）17280―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
R1_A_0=48*0+1;
  R1_A=R1_A_0+340;
R1=R1_A_0:R1_A;//データ読み込み用8/1=8257kara
R2_A_0=1;
R2_A=336;
R2=R2_A_0:R2_A;//3120;//データ解析用
R3=2:R2_A+1;//3121;//データ解析用
R4=1:R2_A;//相関係数用
R5=R2_A_0-1;R6=R2_A;//グラフの左端と右端

R7=B_2008_TERC;//年度指定
R8=B_2008_TERC_temp;


//R7=B_2008_May_Aug_TERC;
//R7=B_2008_Sep_Dec_TERC;
R7=B_2009_Jan_Apr_TERC;
R8=B_2009_Jan_Apr_TERC_temp;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//（データ吸い出し）////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////（風速データ）////////////////////////////////////////////////////////////////////
//////[1.6m]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
Horizon_wind_speed_1=R7(R1,1);//（筑波公開データ：1.6ｍ水平風速データ読み込み）
//////[29.5m]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
Horizon_wind_speed_2_NW=R7(R1,2);//（筑波公開データ：29.5ｍNW水平風速データ読み込み）
Horizon_wind_speed_2_SE=R7(R1,3);//（筑波公開データ：29.5ｍSE水平風速データ読み込み）
Horizon_wind_speed_2=sqrt(Horizon_wind_speed_2_NW^2+Horizon_wind_speed_2_SE^2);//（筑波公開データ：29.5ｍ水平風速データ算出）
////（気温データ）////////////////////////////////////////////////////////////////////
Temperture_0=R8(R1,1);//（筑波公開データ：-0.02m気温データ読み込み）
Temperture_0_1=R8(R1,3);//（筑波公開データ：12.6m気温データ読み込み）
//////[1.6m]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
Temperture_1=R7(R1,12);//（筑波公開データ：1.6m気温データ読み込み）
//////[29.5m]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
Temperture_2=R7(R1,13);//（筑波公開データ：29.5m気温データ読み込み）
////（地温データ）////////////////////////////////////////////////////////////////////
Soil_Temperture=R7(R1,14);//（筑波公開データ：-0.02ｍ地温データ読み込み)
////（相対湿度データ）////////////////////////////////////////////////////////////////////
//////[1.6m]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
Relative_humidity_1=R7(R1,15)./100;//（筑波公開データ：1.6ｍ相対湿度データ読み込み）
//////[29.5m]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
Relative_humidity_2=R7(R1,16)./100;//（筑波公開データ：29.5ｍ相対湿度データ読み込み）
////（大気圧データ）////////////////////////////////////////////////////////////////////
Atmospheric_pressure=R7(R1,17);//（筑波公開データ：1.6ｍ地上大気圧データ読み込み）
////（地温データ）////////////////////////////////////////////////////////////////////
Soil_Temperture=R7(R1,14);//（筑波公開データ：-0.02ｍ地温データ読み込み)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//（湿度計算）//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////（湿度関係:1.6m）//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////[Titens]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
//vapor_pressure_SAT=6.11*10.^((7.5.*Temperture_1)./(237.3+Temperture_1));//飽和水蒸気圧[hPa]
//vapor_pressure=vapor_pressure_SAT.*Relative_humidity;//水蒸気圧[hPa]
//////[Goff-Gratch]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
vapor_pressure_SAT_Log_1=-7.90298*(373.15./(Temperture_1+273.15)-1)+5.02808.*log10(373.15./(Temperture_1+273.15))-1.3816.*10.^(-7).*10.^(11.344.*(1-((Temperture_1+273.15)/373.16))-1)+8.1328.*10^(-3)*10.^(-3.49149.*(373.16./(Temperture_1+273.15)-1)-1)+log10(1013.246);
vapor_pressure_SAT_1=10.^(vapor_pressure_SAT_Log_1);//飽和水蒸気圧[hPa]
vapor_pressure_1=vapor_pressure_SAT_1.*Relative_humidity_1;//水蒸気圧[hPa]
//////[比湿計算]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
Absolute_humidety_SAT_1=0.2167*vapor_pressure_SAT_1./(273.15+Temperture_1);//飽和絶対湿度[kg/kg]
Absolute_humidety_T_1=0.2167*vapor_pressure_1./(273.15+Temperture_1);//絶対湿度[kg/kg]
Specific_humidety_SAT_1=0.622*vapor_pressure_SAT_1./(Atmospheric_pressure-0.378*vapor_pressure_SAT_1);//飽和比湿[kg/kg]
Specific_humidety_T_1=0.622*vapor_pressure_1./(Atmospheric_pressure-0.378*vapor_pressure_1);//比湿[kg/kg]
//Specific_humidety_T_1=Absolute_humidety_T_1./p__1//絶対湿度から比湿算出[kg/kg]
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////（湿度関係:29.5m）//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////[Titens]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
//vapor_pressure_SAT_2=6.11*10.^((7.5.*Temperture_2)./(237.3+Temperture_2));//飽和水蒸気圧[hPa]
//vapor_pressure_2=vapor_pressure_SAT_2.*Relative_humidity_2;//水蒸気圧[hPa]
//////[Goff-Gratch]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
vapor_pressure_SAT_Log_2=-7.90298*(373.15./(Temperture_2+273.15)-1)+5.02808.*log10(373.15./(Temperture_2+273.15))-1.3816.*10.^(-7).*10.^(11.344.*(1-((Temperture_2+273.15)/373.16))-1)+8.1328.*10^(-3)*10.^(-3.49149.*(373.16./(Temperture_2+273.15)-1)-1)+log10(1013.246);
vapor_pressure_SAT_2=10.^(vapor_pressure_SAT_Log_2);//飽和水蒸気圧[hPa]
vapor_pressure_2=vapor_pressure_SAT_2.*Relative_humidity_2;//水蒸気圧[hPa]
//////[比湿計算]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
Absolute_humidety_SAT_2=0.2167*vapor_pressure_SAT_2./(273.15+Temperture_2);//飽和絶対湿度[kg/kg]
Absolute_humidety_T_2=0.2167*vapor_pressure_2./(273.15+Temperture_2);//絶対湿度[kg/kg]
Specific_humidety_SAT_2=0.622*vapor_pressure_SAT_2./(Atmospheric_pressure-0.378*vapor_pressure_SAT_2);//飽和比湿[kg/kg]
Specific_humidety_T_2=0.622*vapor_pressure_2./(Atmospheric_pressure-0.378*vapor_pressure_2);//比湿[kg/kg]
//Specific_humidety_T_2=Absolute_humidety_T_2./p__2;//絶対湿度から比湿算出[kg/kg]
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//（諸所の係数）/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////[簡易法]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
g=9.8066;//重力加速度[kg m s-2]
//Cp=1004;//定圧比熱簡易式
//p_=1.293;//空気密度簡易式
////[1.6m]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
l_1=2.5008*10^6-2300.*Temperture_1;//潜熱量1.6m[J/kg]
Cp_1=1004*(1-Specific_humidety_T_1)+1874*(Specific_humidety_T_1); //定圧比熱1.6m[J/kg・K]
p__1=1.293*(273.15./(273.15+Temperture_1)).*(Atmospheric_pressure./1013.25).*(1-0.378.*(vapor_pressure_1./Atmospheric_pressure)); //空気密度1.6m[kg/m^3]
////[1.6m]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
l_1=2.5008*10^6-2300.*Temperture_1;//潜熱量1.6m[J/kg]
Cp_2=1004*(1-Specific_humidety_T_2)+1874*(Specific_humidety_T_2); //定圧比熱1.6m[J/kg・K]
p__2=1.293*(273.15./(273.15+Temperture_2)).*(Atmospheric_pressure./1013.25).*(1-0.378.*(vapor_pressure_2./Atmospheric_pressure)); //空気密度1.6m[kg/m^3]


//////[1.6m]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
Measured_H_1=R7(R1,7).*Cp_1.*p__1;//（筑波大学が実測した顕熱フラックスデータ読み込み）
//////[29.5m]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
Measured_H_2_NW=R7(R1,8);//（筑波大学が実測した顕熱フラックスNWデータ読み込み）
Measured_H_2_SE=R7(R1,9);//（筑波大学が実測した顕熱フラックスSEデータ読み込み）
Measured_H_2=((Measured_H_2_NW+Measured_H_2_SE)./2).*Cp_1.*p__1;//（筑波大学が実測した顕熱フラックスデータ計算）
////[実測潜熱フラックス：熱収支法により算出]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
//////[1.6m]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
//Measured_lE_1(R1,1)=net_radiation(R1,1)-Underground_heat_flux(R1,1)-Measured_H_1(R1,1);//（筑波大学が実測した潜熱フラックスデータ読み込み）
//////[29.5m]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
//Measured_lE_2(R1,1)=net_radiation(R1,1)-Underground_heat_flux(R1,1)-Measured_H_2(R1,1);//（筑波大学が実測した潜熱フラックスデータ読み込み）
////[実測運動量フラックス]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
//////[1.6m]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
Measured_R_stress_1=p__1.*R7(R1,4);//（筑波大学が実測した運動量フラックスデータ読み込み）
//////[29.5m]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
Measured_R_stress_2_NW=R7(R1,5);//（筑波大学が実測した運動量フラックスNWデータ読み込み）
Measured_R_stress_2_SE=R7(R1,6);//（筑波大学が実測した運動量フラックスSEデータ読み込み）
Measured_R_stress_2=((Measured_R_stress_2_SE+Measured_R_stress_2_NW)./2).*p__1;//（筑波大学が実測した運動量フラックスデータ計算）
////[実測摩擦速度]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
//////[1.6m]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
Friction_velocity_1=sqrt(Measured_R_stress_1./p__1);//（運動量フラックスデータから算出した摩擦速度）
//////[29.5m]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
Friction_velocity_2=sqrt(Measured_R_stress_2./p__2);//(R1,1));//（運動量フラックスデータから算出した摩擦速度）
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
sigma_H_S_1_g(1,1)=0;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//（簡易法：Khデータ読み込み）////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////[1.6m]―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
//////（摩擦速度算出）―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
//Us_1=0.4*Horizon_wind_speed_1./(log((1.6-0.05)/(0.05)));//（風速の対数法則逆算：中立仮定　z0=0.05,d=0.05）
Us_1=Friction_velocity_1;//（実測運動量フラックスより摩擦速度算出）
Us_2=Friction_velocity_2;//（実測運動量フラックスより摩擦速度算出）
sigma_H_S_1_g(R2_A_0,1)=0.005;
sigma_H_S_2_g(R2_A_0,1)=0.005;
Equation_H_1(R2_A_0,1)=2;//初期条件（計算顕熱フラックス）
Equation_H_2(R2_A_0,1)=2;//初期条件（計算顕熱フラックス）
zeta_1(1,1)=1;

H_D_grad_1=(Temperture_0_1-Temperture_1)./(12.3-1.6);
H_D_grad_2=(Temperture_2-Temperture_0_1)./(29.5-12.3);



//(スプライシン補間)////////////////////////////////////////////////////////////////////////////////////////////////////

//load('C:\SC\use scilab clucurate data\binary\B_2008_TERC_temp.sce');

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Measured_H_1(m,1)sigma_H_S_1_g(m,1);//
r=1000;

for b=1:4
Potential_Temperature(R2,b)=B_2009_Jan_Apr_TERC_temp(R2,b).*(r./Atmospheric_pressure(R2,1)).^(0.286);
end
//(スプライシン補間)////////////////////////////////////////////////////////////////////////////////////////////////////

//load('C:\SC\use scilab clucurate data\binary\B_2008_TERC_temp.sce');

for s=R2;
n = 2;  // 与えられたデータ数（標本点の数）：n+1 個

x_data = [1.6,12.3,29.5]; // [x1 x2 x3 x4]
y_data = Potential_Temperature(s,2:4); // [y1 y2 y3 y4]

// ddy の算出（ここから）////////////////////////////////////////////////////
for i = 1:n
  h(i) = x_data(i+1) - x_data(i);
end

for i = 1:n-1
  for j = 1:n-1
    if i == j
      A(i,j) = 2*(h(i)+h(i+1));
    elseif i == j-1
      A(i,j) = h(i+1);
    elseif i == j+1
      A(i,j) = h(i);
    else
      A(i,j) = 0;
    end
  end
end

for i = 1:n-1
  z(i+1) = 6*((y_data(i+2) - y_data(i+1))/h(i+1) - (y_data(i+1) - y_data(i))/h(i));
end
b = z(2:n);

p = A\b;   // 連立１次方程式を解く
ddy = [ 0
        p
        0 ];
num = 0;
x = 1.6:0.1:29.5; // 推定値の x 座標

for num=1:size(x,2)
  // x(num) がどの区間にあるかを判別
  for k = 1:n
    if x(num) >= x_data(k) & x(num) <= x_data(k+1)
      i = k;
      break
    end
  end

  // 自然な３次のスプライン補間
  y(num) = (1/6)*(x_data(i+1) - x(num))^3/h(i)*ddy(i) ...
         + (1/6)*(x(num)   - x_data(i))^3/h(i)*ddy(i+1) ...
         + (1/h(i)*y_data(i)   - h(i)/6*ddy(i)  )*(x_data(i+1) - x(num)) ...
         + (1/h(i)*y_data(i+1) - h(i)/6*ddy(i+1))*(x(num)   - x_data(i));
         
//  printf("x = %f, y = %f\n", x(num), y(num));   // 推定値の表示
end
// 推定値の算出（ここまで）///////////////////////////////////////////////////
B_08Su_TERC_temp_grad_1(s,1)=(y(2,1)-y(1,1))./(x(1,2)-x(1,1));
B_08Su_TERC_temp_grad_2(s,1)=(y(279,1)-y(278,1))./(x(1,279)-x(1,278));
end

for s=R2;
k =3;  // 与えられたデータ数
n =2;  // n次関数で近似
//for m=1:10;
//  load('C:\SC\use scilab clucurate data\binary\test00.sce');  
x_data  = [ 1.6   12.3   29.5 ];  // [x1 x2 x3]
y_data  = Potential_Temperature(s,2:4);//test00(m,1:3);  // [y1 y2 y3]
// 行列 X'*X，X'*y の生成（ここから）/////////////////////////////////////////
XX = zeros(n+1,n+1);  // X'*X の初期化
Xy = zeros(n+1,1);    // X'*y の初期化

for n1 = 1:n+1
    num = 2*n + 1 - n1;
    for n2 = 1:n+1
        for i = 1:k
            XX(n1,n2) = XX(n1,n2) + x_data(i)^(num);  // sum_{i=1}^{k}xi^num
        end
        num = num - 1;
    end
end

for n1 = 1:n+1
    num = n + 1 - n1;
    for i = 1:k
        Xy(n1) = Xy(n1) + x_data(i)^(num)*y_data(i);  // sum_{i=1}^{k}xi^num*yi
    end
end
// 行列 X'*X，X'*y の生成（ここまで）/////////////////////////////////////////

// 最小二乗法によるパラメータ決定と結果の表示（ここから）/////////////////////////
a = XX\Xy;  // 連立１次方程式を解く

//for i = 0:n
//  printf("a%d = %f\n",i,a(i+1));
//end
// 最小二乗法によるパラメータ決定と結果の表示（ここまで）/////////////////////////

// グラフの描画のためのデータ生成（ここから）///////////////////////////////////
x = 1.6:0.1:29.5;
y = zeros(1,size(x,2)); // y の初期化

for num = 1:size(x,2)
  for i = 0:n
    y(num) = y(num) + a(i+1)*x(num)^(n-i);  // y = sum_{i=0}^{n}ai*x^(n-i)
  end
end
// グラフの描画のためのデータ生成（ここまで）///////////////////////////////////
B_08Le_TERC_temp_grad_1(s,1)=(y(1,2)-y(1,1))./(x(1,2)-x(1,1));
B_08Le_TERC_temp_grad_2(s,1)=(y(1,279)-y(1,278))./(x(1,279)-x(1,278));
end



for m=R2;
  Start=500; t=1800;
  z(m,1)=1.6;  C(m,1)=(1.6./Start);//0.0032;  Start=500; t=1800;;
  Temperture_1_2(m,1)=Temperture_1(m+1,1);//（Ti＋１読み込み）
  Temperture_1_1(m,1)=Temperture_1(m,1);//（Ti読み込み）
  derta_Temperture_1(m,1)=Temperture_1_2(m,1)-Temperture_1_1(m,1);//(�儺=Ti+1-Ti)
  grad_1(m,1)=B_08Su_TERC_temp_grad_1(m,1);//(Temperture_0_1(m,1)-Temperture_1(m,1))./(12.3-1.6);//sigma_H_S_1_g(m,1);
  grad_abs_1(m,1)=-grad_1(m,1);
    if m==1;
      zeta_1(m,1)=1;
      Boundary_height_1(m,1)=Start;
    else  
      if zeta_1(m-1,1)<0;
        if Measured_H_1(m,1)>0;
          if grad_abs_1(m,1)>0;
            if derta_Temperture_1(m,1)>0;
               Boundary_growth_height_1(m,1)=sqrt(abs((2*1.86./grad_abs_1(m,1)).*(Measured_H_1(m,1)./(Cp_1(m,1).*p__1(m,1))).*t));
            else
               Boundary_growth_height_1(m,1)=0;
            end
             Boundary_height_1(m,1)=Boundary_height_1(m-1,1)+Boundary_growth_height_1(m,1);
          else
             Boundary_height_1(m,1)=Boundary_height_1(m-1,1);
          end 
        else
          Boundary_height_1(m,1)=Start;//Boundary_height_1(m-1,1);
        end
      else
        Boundary_height_1(m,1)=Start;
      end
    end
   if Boundary_height_1(m,1)>1700;
      Boundary_height_1(m,1)=Boundary_height_1(m-1,1);
   else
      Boundary_height_1(m,1)=Boundary_height_1(m,1);
   end
  z_1(m,1)=Boundary_height_1(m,1).*C(m,1);
  Monin_obukhov_length_1(m,1)=-(Us_1(m,1).^3)./((0.4.*(g./(273.15+Temperture_1(m,1))).*((Equation_H_1(m,1))./(Cp_1(m,1).*p__1(m,1)))));//（モーニン・オブコフ長）
  zeta_1(m,1)=z_1(m,1)./Monin_obukhov_length_1(m,1);//（安定度指数：ζ）
  if zeta_1(m,1)>0;
     phi_H_1(m,1)=1+5.3.*((zeta_1(m,1)+zeta_1(m,1).^1.1.*(1+zeta_1(m,1).^1.1).^((1-1.1)./1.1))./(zeta_1(m,1)+(1+zeta_1(m,1).^1.1).^((1)./1.1)));//（不安定の状態：ζ＞０）
     K_H_1(m,1)=0.4*Us_1(m,1).*z_1(m,1)./phi_H_1(m,1); //（安定の状態：乱流拡散係数）
  else
     phi_H_1(m,1)=(1-15.*zeta_1(m,1)).^(-1/2);
     K_H_1(m,1)=0.4*Us_1(m,1).*z_1(m,1)./phi_H_1(m,1); //（不安定の状態：乱流拡散係数）
  end
  derta_t_1(m,1)=1800;//（�冲=ti+1-ti）
  T_H_1(m,1)=derta_Temperture_1(m,1)./derta_t_1(m,1);//(Tt=(Ti+1-Ti)/(ti+1-ti))
  for e=1:m;
    sigma_H_1(e,1)=((abs(sqrt((m+1-(e))*1800))-abs(sqrt((m+1-(e-1))*1800))).*T_H_1(e,1));//sum{sqrt(tN-ti+1)-sqrt(tN-ti).*Tt(e,1)
    sigma_H_S_1(m,1)=sum(sigma_H_1(1:e,1));
  end 
  Equation_H_1(m+1,1)=(-2*Cp_1(m,1).*p__1(m,1).*sqrt(K_H_1(m,1)./%pi).*sigma_H_S_1(m,1));//顕熱フラックス
  sigma_H_S_1_g(m+1,1)=2./(sqrt(%pi*K_H_1(m,1))).*sigma_H_S_1(m,1);
end


for m=R2;
  z(m,1)=29.5;  C(m,1)=(29.5./Start);//0.059;
  z_2(m,1)=Boundary_height_1(m,1).*C(m,1);
  Monin_obukhov_length_2(m,1)=-(Us_2(m,1).^3)./((0.4.*(g./(273.15+Temperture_2(m,1))).*((Equation_H_2(m,1)+1)./(Cp_2(m,1).*p__2(m,1)))));//（モーニン・オブコフ長）
  zeta_2(m+1,1)=z_2(m,1)./Monin_obukhov_length_2(m,1);//（安定度指数：ζ）
  if zeta_2(m,1)>0;
     phi_H_2(m,1)=1+5.3.*((zeta_2(m,1)+zeta_2(m,1).^1.1.*(1+zeta_2(m,1).^1.1).^((1-1.1)./1.1))./(zeta_2(m,1)+(1+zeta_2(m,1).^1.1).^((1)./1.1)));//（不安定の状態：ζ＞０）
     K_H_2(m,1)=0.4*Us_2(m,1).*z_2(m,1);//./phi_H_2(m,1); //（安定の状態：乱流拡散係数）
  else
     phi_H_2(m,1)=(1-15.*zeta_2(m,1)).^(-1/2);
     K_H_2(m,1)=0.4*Us_2(m,1).*z_2(m,1)./phi_H_2(m,1); //（不安定の状態：乱流拡散係数）
  end
  Temperture_2_2(m,1)=Temperture_2(m+1,1);//（Ti＋１読み込み）
  Temperture_2_1(m,1)=Temperture_2(m,1);//（Ti読み込み）
  derta_Temperture_2(m,1)=Temperture_2_2(m,1)-Temperture_2_1(m,1);//(�儺=Ti+1-Ti)
  derta_t_2(m,1)=1800;//（�冲=ti+1-ti）
  T_H_2(m,1)=derta_Temperture_2(m,1)./derta_t_2(m,1);//(Tt=(Ti+1-Ti)/(ti+1-ti))
  for e=1:m;
    sigma_H_2(e,1)=((abs(sqrt((m+1-(e))*1800))-abs(sqrt((m+1-(e-1))*1800))).*T_H_2(e,1));//sum{sqrt(tN-ti+1)-sqrt(tN-ti).*Tt(e,1)
    sigma_H_S_2(m,1)=sum(sigma_H_2(1:e,1));
  end 
  Equation_H_2(m+1,1)=-2*Cp_2(m,1).*p__2(m,1).*sqrt(K_H_2(m,1)./%pi).*sigma_H_S_2(m,1);//顕熱フラックス
  sigma_H_S_2_g(m+1,1)=2./(sqrt(%pi*K_H_2(m,1))).*sigma_H_S_2(m,1);
end

//if Boundary_height_grad(m,1)==0;
//         Boundary_height_grad(m,1)=4;
//      else
//         Boundary_height_grad(m,1)=Boundary_height_grad(m,1);end
//       Boundary_height_grad(m,1)=Boundary_height(m,1)./Start;     
//.*Boundary_height_grad(m,1)


// グラフの描画 ////////////////////////////////////////////////////////////
clf();                                     // グラフィックの消去
xset("font size",3); 

G1=Measured_H_1(R2,1);
G2=Equation_H_1(R2,1);//.*[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1.5,2,2,3,3,4,4,5,5,5.5,5,5,5,4,4,3,3,2,2,1,1,1,1,1,1,1,1,1,1,1,1]';
G3=Temperture_1(R2,1);
plot2d(R2,[G1,G2,G3*7],[1,2,5],,,[R5,-200,R6,400],[5.5,8,2,7],);//グラフ描画
//plot(R2,sigma_H_S_1(R2,1),'black','linewidth',1,'markersize',1);// データの描画（色：赤，印：○，線の太さ：3，印の大きさ：12）
//plot(R2,grad(R2,1),'b','linewidth',1,'markersize',1);// データの描画（色：赤，印：○，線の太さ：3，印の大きさ：12）
//plot(R2,sigma_H_S_1(R2,1),'black','linewidth',1);// 推定値の描画（色：青，線の太さ：3）
//plot(R2,grad(R2,1),'blue','linewidth',1);// 推定値の描画（色：青，線の太さ：3）


xgrid; // グリッドラインの追加         
//title('data of tsukuba:Comparison for H [09,1,1-1,10;K_H=ku*z/a,z=1.6m]')
xlabel('time(2-8,Jan,09) average of 30min');       // 横軸のラベル
ylabel('Sensible heat flux[W/m-2]');               // 縦軸のラベル

graph=gca();  // 軸のハンドル（個々を認識する番号）を取得
L=graph.children.children;
L(1).line_style=1;
//L(2).line_style=3;
//graph.title.font_style = 1;
graph.title.font_size = 3;
//graph.font_style = 1;           // 目盛りのフォント（2: times）
graph.font_size  = 3;           // 目盛りのフォントサイズ（0〜6まで設定可能）
//graph.x_label.font_style = 1;   // 横軸のラベルのフォント（3: times italic）
graph.x_label.font_size  = 3;   // 横軸のラベルのフォントサイズ（0〜6まで設定可能）
//graph.y_label.font_style = 1;   // 縦軸のラベルのフォント（3: times italic）
graph.y_label.font_size  = 3;   // 縦軸のラベルのフォントサイズ（0〜6まで設定可能）

//graph.data_bounds =[R5 -0.2;R6 0.7];   // [xmin ymin; xmax ymax]

legend('Measured_H_1.6m','Equation_H_1.6m Friction verocity','Equation_H_1.6m Wind speed verocity',1);//凡例表示
xgrid; // グリッドラインの追加   
//////////////////////////////////////////////////////////////////////////
////（Hのグラフ）―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
scf();
X=-200:300;
plot2d(Measured_H_1(R2,1)',Equation_H_1(R2,1)',style=-9,rect=[-200,-200,300,300]);//散布図描画
plot2d(X,X);//ｙ＝ｘの線追加
square(-200,-200,300,300);//散布図を正方形にする

xgrid; // グリッドラインの追加         
//title('data of tsukuba:Comparison for H [09,1,1-1,10;K_H=ku*z/a,z=1.6m]')
xlabel('Measured_H_1.6m');       // 横軸のラベル
ylabel('Equation_H_1.6m');               // 縦軸のラベル

graph=gca();  // 軸のハンドル（個々を認識する番号）を取得
L=graph.children.children;
L(1).line_style=1;
//L(2).line_style=3;
//graph.title.font_style = 1;
graph.title.font_size = 3;
//graph.font_style = 1;           // 目盛りのフォント（2: times）
graph.font_size  = 3;           // 目盛りのフォントサイズ（0〜6まで設定可能）
//graph.x_label.font_style = 1;   // 横軸のラベルのフォント（3: times italic）
graph.x_label.font_size  = 3;   // 横軸のラベルのフォントサイズ（0〜6まで設定可能）
//graph.y_label.font_style = 1;   // 縦軸のラベルのフォント（3: times italic）
graph.y_label.font_size  = 3;   // 縦軸のラベルのフォントサイズ（0〜6まで設定可能）

//graph.data_bounds =[R5 -0.2;R6 0.7];   // [xmin ymin; xmax ymax]

xgrid; // グリッドラインの追加   
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

F1_1=Measured_H_1(R4,1);//相関係数計算
G1_1=Equation_H_1(R4,1);//相関係数計算
R_2_H_1=(F1_1'*G1_1/length(F1_1)-mean(F1_1)*mean(G1_1))/stdev(F1_1)/stdev(G1_1)//相関係数計算

RMSD_1=sqrt((sum((Measured_H_1(R2,1)-Equation_H_1(R2,1))^2))./R2_A);

//////////////////////////////////////////////////////////////////////////
//[p,q]=reglin(Measured_H_1(R2,1),Equation_H_1(R2,1));
x = [Measured_H_1(R2,1)];
     
y = [Equation_H_1(R2,1)];

//係数aa,bbを求める．
 n = size(x,1);
 
 sum_x1y1 = 0;
 sum_y1 = 0;
 sum_x2 = 0;
 sum_x1 = 0;
 
 xone = ones(x);
 
  sum_x1 = xone'*x;
  sum_y1 = xone'*y;
  sum_x2 = x'*x;
  sum_x1y1 = x'*y;
    
  Coef_Matrix = inv([sum_x2 sum_x1;sum_x1 n])*[sum_x1y1;sum_y1];
 
  aa_1 = Coef_Matrix(1);//aa=傾き
  bb_1 = Coef_Matrix(2);
  
// グラフ表示
//plot(x,y,'O',x,aa*x+bb)

//f=get("current_figure");
//a=f.children; // the handle on the Axes child
//a.data_bounds=[min(x)-abs(min(x))*0.1 max(x)+abs(max(x))*0.1 min(y)-abs(min(y))*0.1 max(y)+abs(max(y))*0.1]; // xaxis_left xaxis_right yaxis_down yaxis_up 
//a.x_label.text = "x";
//a.y_label.text = "y";
//xgrid(2)
//////////////////////////////////////////////////////////////////////////



// グラフの描画 ////////////////////////////////////////////////////////////
scf();                                     // グラフィックの消去
xset("font size",3); 

G1=Measured_H_2(R2,1);
G2=Equation_H_2(R2,1);//.*[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1.5,2,2,3,3,4,4,5,5,5.5,5,5,5,4,4,3,3,2,2,1,1,1,1,1,1,1,1,1,1,1,1]';
G3=Temperture_2(R2,1);
plot2d(R2,[G1,G2,G3*7],[1,2,3],,,[R5,-200,R6,400],[5.5,8,2,7],);//グラフ描画
//plot(R2,sigma_H_S_1(R2,1),'black','linewidth',1,'markersize',1);// データの描画（色：赤，印：○，線の太さ：3，印の大きさ：12）
//plot(R2,grad(R2,1),'b','linewidth',1,'markersize',1);// データの描画（色：赤，印：○，線の太さ：3，印の大きさ：12）
//plot(R2,sigma_H_S_1(R2,1),'black','linewidth',1);// 推定値の描画（色：青，線の太さ：3）
//plot(R2,grad(R2,1),'blue','linewidth',1);// 推定値の描画（色：青，線の太さ：3）


xgrid; // グリッドラインの追加         
//title('data of tsukuba:Comparison for H [09,1,1-1,10;K_H=ku*z/a,z=1.6m]')
xlabel('time(2-8,Jan,09) average of 30min');       // 横軸のラベル
ylabel('Sensible heat flux[W/m-2]');               // 縦軸のラベル

graph=gca();  // 軸のハンドル（個々を認識する番号）を取得
L=graph.children.children;
L(1).line_style=1;
//L(2).line_style=3;
//graph.title.font_style = 1;
graph.title.font_size = 3;
//graph.font_style = 1;           // 目盛りのフォント（2: times）
graph.font_size  = 3;           // 目盛りのフォントサイズ（0〜6まで設定可能）
//graph.x_label.font_style = 1;   // 横軸のラベルのフォント（3: times italic）
graph.x_label.font_size  = 3;   // 横軸のラベルのフォントサイズ（0〜6まで設定可能）
//graph.y_label.font_style = 1;   // 縦軸のラベルのフォント（3: times italic）
graph.y_label.font_size  = 3;   // 縦軸のラベルのフォントサイズ（0〜6まで設定可能）

//graph.data_bounds =[R5 -0.2;R6 0.7];   // [xmin ymin; xmax ymax]

legend('Measured_H_29.5m','Equation_H_29.5m Friction verocity','Equation_H_29.5m Wind speed verocity',1);//凡例表示
xgrid; // グリッドラインの追加   
//////////////////////////////////////////////////////////////////////////
////（Hのグラフ）―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
scf();
X=-200:400;
plot2d(Measured_H_2(R2,1)',Equation_H_2(R2,1)',style=-9,rect=[-200,-200,400,400]);//散布図描画
plot2d(X,X);//ｙ＝ｘの線追加
square(-200,-200,400,400);//散布図を正方形にする

xgrid; // グリッドラインの追加         
//title('data of tsukuba:Comparison for H [09,1,1-1,10;K_H=ku*z/a,z=1.6m]')
xlabel('Measured_H_29.5m');       // 横軸のラベル
ylabel('Equation_H_29.5m');               // 縦軸のラベル

graph=gca();  // 軸のハンドル（個々を認識する番号）を取得
L=graph.children.children;
L(1).line_style=1;
//L(2).line_style=3;
//graph.title.font_style = 1;
graph.title.font_size = 3;
//graph.font_style = 1;           // 目盛りのフォント（2: times）
graph.font_size  = 3;           // 目盛りのフォントサイズ（0〜6まで設定可能）
//graph.x_label.font_style = 1;   // 横軸のラベルのフォント（3: times italic）
graph.x_label.font_size  = 3;   // 横軸のラベルのフォントサイズ（0〜6まで設定可能）
//graph.y_label.font_style = 1;   // 縦軸のラベルのフォント（3: times italic）
graph.y_label.font_size  = 3;   // 縦軸のラベルのフォントサイズ（0〜6まで設定可能）

//graph.data_bounds =[R5 -0.2;R6 0.7];   // [xmin ymin; xmax ymax]

xgrid; // グリッドラインの追加   
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

F1_1=Measured_H_2(R4,1);//相関係数計算
G1_1=Equation_H_2(R4,1);//相関係数計算
R_2_H_2=(F1_1'*G1_1/length(F1_1)-mean(F1_1)*mean(G1_1))/stdev(F1_1)/stdev(G1_1)//相関係数計算

RMSD_2=sqrt((sum((Measured_H_2(R2,1)-Equation_H_2(R2,1))^2))./R2_A);

//////////////////////////////////////////////////////////////////////////
//[p,q]=reglin(Measured_H_1(R2,1),Equation_H_1(R2,1));
x = [Measured_H_2(R2,1)];
     
y = [Equation_H_2(R2,1)];

//係数aa,bbを求める．
 n = size(x,1);
 
 sum_x1y1 = 0;
 sum_y1 = 0;
 sum_x2 = 0;
 sum_x1 = 0;
 
 xone = ones(x);
 
  sum_x1 = xone'*x;
  sum_y1 = xone'*y;
  sum_x2 = x'*x;
  sum_x1y1 = x'*y;
    
  Coef_Matrix = inv([sum_x2 sum_x1;sum_x1 n])*[sum_x1y1;sum_y1];
 
  aa_2 = Coef_Matrix(1);//aa=傾き
  bb_2 = Coef_Matrix(2);
  
// グラフ表示
//plot(x,y,'O',x,aa*x+bb)

//f=get("current_figure");
//a=f.children; // the handle on the Axes child
//a.data_bounds=[min(x)-abs(min(x))*0.1 max(x)+abs(max(x))*0.1 min(y)-abs(min(y))*0.1 max(y)+abs(max(y))*0.1]; // xaxis_left xaxis_right yaxis_down yaxis_up 
//a.x_label.text = "x";
//a.y_label.text = "y";//
//xgrid(2)
//////////////////////////////////////////////////////////////////////////
scf();
plot2d(R2,[sigma_H_S_1_g(R2,1),B_08Le_TERC_temp_grad_1(R2,1)])
scf();
plot2d(R2,[sigma_H_S_2_g(R2,1),H_D_grad_2(R2,1)])
//plot2d(R2,[B_08Le_TERC_temp_grad_1(R2,1),H_D_grad_2(R2,1)])

G1_1=sigma_H_S_1_g;G1_2=B_08Le_TERC_temp_grad_1;
G2_1=sigma_H_S_2_g;G2_2=H_D_grad_2;

scf();
X=-200:300;
plot2d(X,X);//ｙ＝ｘの線追加
xtitle('data of tsukuba:Comparison for gradient T [1-10,Jan,09 In the daytime;K=ku*z/phi_H z=29.5m]','12.3m-29.5m:real_gradient[dT/dz]','29.5m:half-order_time_method_gradient[dT/dz]');//グラフタイトルと横軸、縦軸
square(-0.5,-0.5,1,1);//散布図を正方形にする
//legend('',1);//凡例表示
xgrid();//グリッド線の表示
plot2d(G1_1(R2,1)',(G1_2(R2,1))',style=-9,rect=[-0.5,-0.5,1,1]);//散布図描画

//for r=1:r_0;
//R9=(16:32)+48*(r-1);
//plot2d(G2_1(R9,1)',(G2_2(R9,1))',style=-9,rect=[-0.05,-0.05,0.05,0.05]);//散布図描画
//end


graph=gca();  // 軸のハンドル（個々を認識する番号）を取得
L=graph.children.children;
L(1).line_style=1;
//L(2).line_style=3;
//graph.title.font_style = 1;
graph.title.font_size = 3;
//graph.font_style = 1;           // 目盛りのフォント（2: times）
graph.font_size  = 3;           // 目盛りのフォントサイズ（0〜6まで設定可能）
//graph.x_label.font_style = 1;   // 横軸のラベルのフォント（3: times italic）
graph.x_label.font_size  = 3;   // 横軸のラベルのフォントサイズ（0〜6まで設定可能）
//graph.y_label.font_style = 1;   // 縦軸のラベルのフォント（3: times italic）
graph.y_label.font_size  = 3;   // 縦軸のラベルのフォントサイズ（0〜6まで設定可能）

//graph.data_bounds =[R5 -0.2;R6 0.7];   // [xmin ymin; xmax ymax]

scf();
X=-200:300;
plot2d(X,X);//ｙ＝ｘの線追加
xtitle('data of tsukuba:Comparison for gradient T [1-10,Jan,09 In the daytime;K=ku*z/phi_H z=29.5m]','12.3m-29.5m:real_gradient[dT/dz]','29.5m:half-order_time_method_gradient[dT/dz]');//グラフタイトルと横軸、縦軸
square(-0.5,-0.5,1,1);//散布図を正方形にする
//legend('',1);//凡例表示
xgrid();//グリッド線の表示
plot2d(G2_1(R2,1)',(G2_2(R2,1))',style=-9,rect=[-0.5,-0.5,1,1]);//散布図描画

//for r=1:r_0;
//R9=(16:32)+48*(r-1);
//plot2d(G2_1(R9,1)',(G2_2(R9,1))',style=-9,rect=[-0.05,-0.05,0.05,0.05]);//散布図描画
//end


graph=gca();  // 軸のハンドル（個々を認識する番号）を取得
L=graph.children.children;
L(1).line_style=1;
//L(2).line_style=3;
//graph.title.font_style = 1;
graph.title.font_size = 3;
//graph.font_style = 1;           // 目盛りのフォント（2: times）
graph.font_size  = 3;           // 目盛りのフォントサイズ（0〜6まで設定可能）
//graph.x_label.font_style = 1;   // 横軸のラベルのフォント（3: times italic）
graph.x_label.font_size  = 3;   // 横軸のラベルのフォントサイズ（0〜6まで設定可能）
//graph.y_label.font_style = 1;   // 縦軸のラベルのフォント（3: times italic）
graph.y_label.font_size  = 3;   // 縦軸のラベルのフォントサイズ（0〜6まで設定可能）

//graph.data_bounds =[R5 -0.2;R6 0.7];   // [xmin ymin; xmax ymax]





R_2_H_1
R_2_H_1^2
RMSD_1
aa_1
bb_1

R_2_H_2
R_2_H_2^2
RMSD_2
aa_2
bb_2


toc()//時間計測終了















