//THIS IS NEW HAYAKAWA CODE(2023/4/5)
//林円盤とQ円盤でのダスト軌道とダスト成長を計算する
//単位系はcgs

#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#define yr 3.16e7 //1yr(s)
#define yrr 1/yr
#define AU 1.50e13 //1AU(cm)
#define AUU 1/AU
#define Msun 1.99e33 //太陽質量(g)
#define Msunn 1/Msun
#define G 6.67e-8 //万有引力定数(cm3/g*s2)
#define GG 1/G

//定数
#define mxstep 5e20;  //step数
#define pi 4.0 * atan(1.0); //π
#define pii 1/pi; //πの逆数
#define kB 1.38e-16;  //ボルツマン定数

#define C1 3/7;
#define C2 1/280;

int main(void){

//宣言部分
double f[11][8],dfdt[11][8];                  
double aa[11],bb[11][10],cc[11];                  
void rkprmt(double aa[11], double bb[11][10],double cc[11]);
rkprmt(aa,bb,cc);   //ルンゲクッタ法に必要な値
int i0,i1,i2,i3,i4,i5,i6,i7;  //パラメータ
int i, k, j;  //ルンゲクッタ

//各パラメータを定義する
double ch[i0] = {1, 2};     //林円盤かQ円盤か
double x_0[i1] = {10,100};   //ダスト初期位置(au)
double ad_0[i2] = {1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.1, 1, 10, 100} ; //ダスト初期半径(cm)　0.01μm~1m(9種類)
double rhop[i3] = {1, 2.5}; //ダスト物質密度(g/cm3) 氷・シリケイト
double M_s[i4] = {1};   //中心星質量(g)　M*=1として規格化している
double alpha[i5] = {0.001};   //乱流パラメータ
double dvfrag[i6] = {100, 1000, 3000};   //ダスト衝突破壊速度(cm/s)
double Q[i7] = {1};   //Toomre Q

double x, y, z, x_0, y_0, z_0, r, rr, r_0;
double sin_cur_0, sin_cur, count_0, count;
double vx_0, vy_0, vz_0, vx, vy, vz, v_0, v;
double istep, t, tt, tt0, dt;
double ad;
double x_pre, y_pre, z_pre, r_pre, rr_pre, sin_pre;
double Fg, Fgx, Fgy, omega, omegaa, vdust, eta, vkep, vgas;
double vrelx, vrely, vrel, dvdg;
double Tdisk, cs, sigg, sigd;
double hg, pg, lg, nyu, R, St, hd, pd, md, nd, M;
double omeg, Cd, Fxx, Fyy, Fx, Fy, F;
double dvB, ts, ts1, ts2, delvg, Dg, Ret, tn, dvtulb, dvdd;
double tgrow, num1, num2, Agl;
double Ftestx, Ftesty;
double dt1, dt2, tyr, wstep, mstep;



for(i0 = 0; i0 < 2; i0++){
  ch = ch[i0];
  for(i1 = 0; i1 < 2; i1++){
    x_0 = x_0[i1];
    for(i2 = 0; i2 < 9; i2++){
      ad_0 = ad_0[i2];
      for(i3 = 0; i3 < 2; i3++){
       rhop = rhop[i3];
        for(i4 = 0; i4 < 1; i4++){
         M_s = M_s[i4];
         for(i5 = 0; i5 < 1; i5++){
           alpha = alpha[i5];
           for(i6 = 0; i6 < 3; i6++){
             dvfrag = dvfrag [i6];
             for(i7 = 0; i7 < 1; i7++){
               Q = Q[i7];

//初期ダスト位置
y_0 = 0;
z_0 = 0;
r_0 = sqrt(x_0*x_0 + y_0*y_0 + z_0*z_0);  //星間距離
//回転数に用いる部分
sin_cur_0 = y_0 / r_0;  //sinθ
count_0 = 0;  //周回数
//初期ダスト速度
vx_0 = 0;
vy_0 = sqrt(1/x_0); //ケプラー速度で発射
vz_0 = 0;
v_0 = sqrt(vx_0*vx_0 + vy_0*vy_0 + vz_0*vz_0);
//
istep = 0;
t = 0;
dt = 1e-8;

for(istep = 0; istep < mxstep; istep++){

tt = t;
t = tt + dt;  //時間(s)

//微分方程式への代入 
f[0][1] = x_0;
f[0][2] = y_0;
f[0][3] = z_0;
f[0][4] = vx_0;
f[0][5] = vy_0;
f[0][6] = vz_0;
f[0][7] = ad_0;

//ルンゲクッタ法による、誤差を小さくする計算部分 
for(i = 0; i <= 10; i++){
 tt = tt0 + aa[i]*dt;
 if(i! = 0){
   for(k = 1; k <= 7; k++){
     f[i][k] = f[0][k];
     }
     for(j = 0; j <= (i-1); j++){
       for(k = 1; k <= 7; k++){
         f[i][k] = f[i][k] + bb[i][j]*dfdt[j][k]*dt;
         }
         }
         }

//微分方程式への代入部分
dfdt[i][1] = f[i][4];
dfdt[i][2] = f[i][5];
dfdt[i][3] = f[i][6];

x = f[i][1];
y = f[i][2];
z = f[i][3];

vx = f[i][4];
vy = f[i][5];
vz = f[i][6];

dfdt[i][7] = f[i][7]; //ダスト成長式
ad = f[i][7];
    
r = sqrt(x*x + y*y + z*z);  //星間距離
rr = 1/r; //星間距離の逆数


//前のステップ(i-1),回転数を計算するための部分
x_pre = f[i-1][1];
y_pre = f[i-1][2];
z_pre = f[i-1][3];
r_pre = sqrt(x_pre*x_pre + y_pre*y_pre + z_pre*z_pre);
rr_pre = 1/r_pre: //１個前の星間半径の逆数

sin_pre = y_pre*rr_pre;
sin_cur = y*rr;
if(sin_pre*sin_cur <= 0){
  count = count + 0.5;   //周回のカウント
}

Fg = rr*rr; //中心力(G,M規格化)
Fgx = Fg*xt*rr;
Fgy = Fg*yt*rr;

omega = sqrt(rr*rr*rr);   //角速度
omegaa = 1/omega; //omegaの逆数
vdust = sqrt(vx*vx + vy*vy);  //ダスト速度
eta = 1.8e-3*sqrt(r);    //ズレ
vkep = sqrt(rr);  //ケプラー速度
vgas = (1 - eta)*vkep;  //ガス速度

//ダストガス相対速度
vrelx = vx - (-vgas*y*rr);
vrely = vy - vgas*x*rr;
vrel = sqrt(vrelx*vrelx + vrely*vrely);
dvdg = vrel;  //ダストガス相対速度

//円盤モデル(円盤温度・音速・円盤ガス面密度・円盤ダスト面密度)
if(ch = 1){
  Tdisk = 280*sqrt(rr); //L = L
  cs = 9.9e4*sqrt(Tdisk*C2);     //音速(cm/s),μ = 2.34,C2 = 1/280
  sigg = 1700*sqrt(rr*rr*rr); //円盤ガス面密度(g/cm2)
  if(0.35 < r && r < 2.7){
      sigd = 7.1*sqrt(rr*rr*rr);  //rock
    }else if(r > 2.7){
      sigd = 30*sqrt(rr*rr*rr); //rock+ice
    } //円盤ダスト面密度(g/cm2)

}else if(ch = 2){
  Tdisk = 150*pow(rr,C1); //C1 = 3/7
  cs = 9.9e4*sqrt(Tdisk*C2);     //音速(cm/s),μ = 2.34,C2 = 1/280
  sigg = cs*omega*pii*GG*yrr/Q;
  sigd = 0.01*sigg;
}//ch = 1なら林円盤、ch = 2ならQ円盤

hg = yr*cs*omegaa; //円盤のスケールハイト(cm)
pg = sigg/(sqrt(2*pi)*hg);  //ガスの質量空間密度(g/cm^3)
lg = 3.9e-9/pg; //平均自由行程(cm)
nyu = 0.353*sqrt(8*pii)*cs*lg;  //動粘性係数(cm2/s)
R = 2*(AU*yrr)*vkep*ad*fabs(dvdg)/nyu; //レイノルズ数
St = pi*rhop*ad/(2*sigg); //Stokes数
hd = hg*pow(1 + ((St/alpha)*(1 + 2*St)/(1 + St)),- 0.5);   //ダストスケールハイト(cm),Arakawa2021
pd = sigd/(sqrt(2*pi)*hd); //g/cm3 ダストの質量空間密度
md = 4*pi*ad*ad*ad*rhop/3; //g ダスト空間質量
nd = pd/md; //ダスト数密度/cm3
M  = (AU*yrr)*fabs(dvdg)/cs; //マッハ数

//補完係数の決定
if(R < 2.0e5){
  omeg = 0.4;
  }else{
    omeg= 0.2;
    }
Cd = pow((0.23*M + pow((24/R + 40/(10 + R)),-1)),-1) + (2 - omeg)*M/(1.6 + M) + omeg; //ガス抵抗係数(無次元)

Fxx = (- 3)*AU*Cd*pg*fabs(vrelx)*vrelx/(8*ad*rhop);
Fyy = (- 3)*AU*Cd*pg*fabs(vrely)*vrely/(8*ad*rhop);
Fx = Fxx; //ガス抵抗力x成分
Fy = Fyy; //ガス抵抗力y成分
F = sqrt(Fx*Fx + Fy*Fy);  //ガス抵抗力(無次元)

//ダストダスト相対速度(多粒子と仮定)
dvB = (yr*AUU)*4*sqrt(kB*Tdisk*pii/md);  //ブラウン運動(全ダスト質量同等)(無次元)

ts = (Msunn)*md*dvdg/F;	//stopping time(無次元)
ts1 = ts;
ts2 = 0.5*ts1;
     
delvg = (yr*AUU)*sqrt(alpha)*cs;	//渦速度
Dg = delvg*delvg/omega;	//拡散係数
Ret = (AU*AU*yrr)*Dg/nyu;	//レイノルズ数(乱流)
tn = pow(Ret,- 0.5)*omega;	//回転時間

if(ts1 < tn){
  dvtulb = delvg*pow(Ret,0.25)*omega*fabs(ts1 - ts2);
  }else if(tn < ts1 && ts1 < 1/omega){
    dvtulb = 1.5*delvg*sqrt(omega*ts1);
    }else if(1/omega < ts1){
      dvtulb = delvg*sqrt((1/(1 + omega*ts1)) + (1/(1 + omega*ts2)));	
      }	//乱流(Ormel and Cuzzi 2007)(無次元)
     
dvdd = dvB + dvtulb;  //ダスト同士の相対速度

tgrow = pii*AUU/(nd*dvdd*ad*ad);    //成長タイムスケール(無次元)
num1 = 1;
num2 = -(log(dvdd/dvfrag)/log(5));
if(num1 < num2){
  Agl = num1;
  }else if(num1 > num2){
    Agl = num2;
    }
     
//運動方程式記述部分
dfdt[i][4] = -Fgx + Fx;
dfdt[i][5] = -Fgy + Fy;
dfdt[i][6] = -rr*rr*rr*z;
dfdt[i][7] = Agl*ad/(3*tgrow); //ダスト進化,Tsukamoto et al.2021	

//force test
Ftestx = Fgx/Fx;
Ftesty = Fgy/Fy;

}//ルンゲクッタによる誤差の終了

for(j = 0; j <= 10; j++){
  for(k = 1; k <= 7; k++){
    f[0][k] = f[0][k] + cc[j]*dfdt[j][k]*dt;
    }
    }

//微分方程式を解いた結果を代入する部分
x = f[0][1];
y = f[0][2];
z = f[0][3];
vx = f[0][4];
vy = f[0][5];
vz = f[0][6];
ad = f[0][7];

//タイムスケールの選定
dt1 = 1e-8*r;  //星間半径で調整　1e-8で0.01μm,1e-6で1μmのダストを計算できる
dt2 = 0.01*tgrow;    //ダスト成長タイムスケール

if(dt1 < dt2){
   dt = dt1;
}else if(dt1 > dt2){
   dt = dt2;
}

tyr = t; //年に変換

//終了条件
	if (r > 500){
		printf("飛んでった\n");
	      	break;
	}
	if (r < 0.35){
		printf("落下\n");
	      	break;
	}
   if (tyr > 100000){
		printf("通常運転\n");
	      	break;
	}

/*
if(tyr > tone){
    printf("一周しました\n");
    break;
}*/

wstep = 1e5*r;             /* stepをいくつ飛ばすか */
mstep = 1e6*r;

//🐷(4/10 14:50)出力パートから


/*計算過程の表示*/
if((istep%mstep)==0){
printf("%ld %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %lf %e %e %e %e %e %e\n",	
		istep,xt,yt,rt,sin_pre,sin_cur,rpt,vxt,vyt,vdust,count,Fx,Fy,Fg,omega,eta,vkep,vgas,vrelx,vrely,dvdg,Tdisk,cs,hg,sigg,pg,lg,nyu,R,sigd,St,hd,pd,md,nd,dvB,dvtulb,dvdd,tgrow,Agl,M,omeg,Cd,Fx,Fy,F,dt1,dt2,dt,tyr,ts1,ts2,delvg,Ret,Dg,tn); 
}
         /* グラフ書き出しのための導出部分 */
        if((istep%wstep)==0)
			{
			fprintf(fp,"%ld %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %lf %e %e %e %e %e %e\n",
            istep,xt,yt,rt,sin_pre,sin_cur,rpt,vxt,vyt,vdust,count,Fx,Fy,Fg,omega,eta,vkep,vgas,vrelx,vrely,dvdg,Tdisk,cs,hg,sigg,pg,lg,nyu,R,sigd,St,hd,pd,md,nd,dvB,dvtulb,dvdd,tgrow,Agl,M,omeg,Cd,Fx,Fy,F,dt1,dt2,dt,tyr,ts1,ts2,delvg,Ret,Dg,tn); 
            }
    } /*itepに対応*/

 printf("rt=%e\n",rt);
 printf("rpt=%e \n",rpt);
 printf("tyr=%lf\n",tyr);
 printf("count=%lf",count);
 /*printf("tone=%lf\n",tone); 一周*/

    long cpu_time;
    double sec,mi;
    /* CPU時間をチェック */
    cpu_time = clock();
    /* 秒に直す */
    sec = (double)cpu_time / CLOCKS_PER_SEC;
    mi = sec/60; 
    printf("計算実時間%f分\n", mi);

    fprintf(fp,"# rt(ダスト位置)=%lf rpt(ダスト半径)=%lf tyr(年)=%lf count(周回数)=%lf mi(計算実時間)=%lf \n",rt,rpt,tyr,count,mi); //計算後パラメータ

    return (0);
}

              }
            }
          }
        }
      }
    }
  }
}


//GPTpart
for (int i = 0; i < 8; i++) {
  char filename[20];
  sprintf(filename, "result%d.txt", i);
  FILE* fp = fopen(filename, "w");
  if (fp == NULL) {
    printf("ファイルを開けませんでした。\n");
    exit(1);
  }

  double a = params[i][0];
  double b = params[i][1];
  double c = params[i][2];

  double result1 = a + b + c;
  double result2 = a * b * c;
  double result3 = (a + b) * c;

  fprintf(fp, "result1: %f\n", result1);
  fprintf(fp, "result2: %f\n", result2);
  fprintf(fp, "result3: %f\n", result3);

  fclose(fp);
}

/*ルンゲクッタ法に用いられる係数の記述*/
void rkprmt(double aa[11], double bb[11][10],double cc[11])
{
    double sqrt21;
    sqrt21=sqrt(21.0);
    aa[0]=0.0;
    aa[1]=0.5;
    aa[2]=aa[1];
    aa[3]=(7.0+sqrt21)/14.0;
    aa[4]=aa[3];
    aa[5]=aa[1];
    aa[6]=(7.0-sqrt21)/14.0;
    aa[7]=aa[6];
    aa[8]=aa[1];
    aa[9]=aa[3];
    aa[10]=1.0;
    bb[1][0]=0.5;
    bb[2][0]=0.25;
    bb[2][1]=0.25;
    bb[3][0]=1.0/7.0;
    bb[3][1]=(-7.0-3.0*sqrt21)/98.0;
    bb[3][2]=(21.0+5.0*sqrt21)/49.0;
    bb[4][0]=(11.0+sqrt21)/84.0;
    bb[4][1]=0.0;
    bb[4][2]=(18.0+4.0*sqrt21)/63.0;
    bb[4][3]=(21.0-sqrt21)/252.0;
    bb[5][0]=(5.0+sqrt21)/48.0;
    bb[5][1]=0.0;
    bb[5][2]=(9.0+sqrt21)/36.0;
    bb[5][3]=(-231.0+14.0*sqrt21)/360.0;
    bb[5][4]=(63.0-7.0*sqrt21)/80.0;
    bb[6][0]=(10.0-sqrt21)/42.0;
    bb[6][1]=0.0;
    bb[6][2]=(-432.0+92.0*sqrt21)/315.0;
    bb[6][3]=(633.0-145.0*sqrt21)/90.0;
    bb[6][4]=(-504.0+115.0*sqrt21)/70.0;
    bb[6][5]=(63.0-13.0*sqrt21)/35.0;
    bb[7][0]=1.0/14.0;
    bb[7][1]=0.0;
    bb[7][2]=0.0;
    bb[7][3]=0.0;
    bb[7][4]=(14.0-3.0*sqrt21)/126.0;
    bb[7][5]=(13.0-3.0*sqrt21)/63.0;
    bb[7][6]=1.0/9.0;
    bb[8][0]=1.0/32.0;
    bb[8][1]=0.0;
    bb[8][2]=0.0;
    bb[8][3]=0.0;
    bb[8][4]=(91.0-21.0*sqrt21)/576.0;
    bb[8][5]=11.0/72.0;
    bb[8][6]=(-385.0-75.0*sqrt21)/1152.0;
    bb[8][7]=(63.0+13.0*sqrt21)/128.0;
    bb[9][0]=1.0/14.0;
    bb[9][1]=0.0;
    bb[9][2]=0.0;
    bb[9][3]=0.0;
    bb[9][4]=1.0/9.0;
    bb[9][5]=(-733.0-147.0*sqrt21)/2205.0;
    bb[9][6]=(515.0+111.0*sqrt21)/504.0;
    bb[9][7]=(-51.0-11.0*sqrt21)/56.0;
    bb[9][8]=(132.0+28.0*sqrt21)/245.0;
    bb[10][0]=0.0;
    bb[10][1]=0.0;
    bb[10][2]=0.0;
    bb[10][3]=0.0;
    bb[10][4]=(-42.0+7.0*sqrt21)/18.0;
    bb[10][5]=(-18.0+28.0*sqrt21)/45.0;
    bb[10][6]=(-273.0-53.0*sqrt21)/72.0;
    bb[10][7]=(301.0+53.0*sqrt21)/72.0;
    bb[10][8]=(28.0-28.0*sqrt21)/45.0;
    bb[10][9]=(49.0-7.0*sqrt21)/18.0;
    cc[0]=0.05;
    cc[1]=0.0;
    cc[2]=0.0;
    cc[3]=0.0;
    cc[4]=0.0;
    cc[5]=0.0;
    cc[6]=0.0;
    cc[7]=49.0/180.0;
    cc[8]=16.0/45.0;
    cc[9]=49.0/180.0;
    cc[10]=0.05;
}
