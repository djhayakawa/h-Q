

#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#define yr 3.16e7 /*1yr(s)*/
#define AU 1.50e13 /*1AU(cm)*/
#define Msun 1.99e33 /*太陽質量(g)*/
#define G 6.67e-8 /*万有引力定数(cm3/g*s2)*/

int main(void)
{
    double x,y,z,xt,yt,zt;
    double v,vx,vy,vz,vxg,vxgd,vyg,vygd,vxt,vyt,vzt; 
    double a1,a2,b1,b2;
    double vgas,vrelx,vrely,vrelz;
    double t,tt,dt0,dt,tt0,o,bs,rr;                    
    double pi,m1,m2,m3,M,a,rH,d;            
    double r,rt,atpow3i,rmin;
    double f[11][8],dfdt[11][8];                  
    double aa[11],bb[11][10],cc[11];                  
    void rkprmt(double aa[11], double bb[11][10],double cc[11]);
    long mxstep,istep,i,j,k,kizami,wstep,mstep,itstep;
    double Fd,Fdx,Fdy,Fdz,Cd,rp,rho,rhogd,rhop,Sigma,Mmh,T,H1,H2,H,Eta,vrel;     
    double jacobi0,jacobi,Emax,c,Vg,e0,e,de,dep,demax,dek; 
    double E,Etot0,Etot,K0,K,U0,U,Eper,tyr;
    double rpt,nd,vdust,tgrow,Agl,omega,num1,num2,dvdg,dvdd,dvB,dvtulb,vkep,eta,pg,pd,md;
    double growth(double dvdd,double vfrag),ranryu();
    double dt1,dt2;
    double teikou(double pg,double dvdg,double cs, double pi,double omega,double R);
    double nyu,cs,lg,R,Fx,Fy,F,omeg;
    double hg,hd,Tdisk,sigg,sigd,Q,c1,Mstar,c2,csd,c3,c4;
    double Fxx,Fyy,Fg,Fgx,Fgy,Ftestx,Ftesty;
    double tone,St,kB,alpha,dvfrag,count,sin_cur,sin_pre,r_pre,x_pre,y_pre,z_pre;
    double ts,ts1,ts2,delvg,Ret,Dg,tn;

    //定数
                                           
    mxstep=5e18;      /* step数 */                     
    rhop=2.5;   /*ダスト内部密度*/   
    rkprmt(aa,bb,cc);        /* ルンゲクッタ法に必要な値 */
    pi=4.0*atan(1.0);    /*π*/
    alpha = 1e-3;  /*乱流定数*/
    kB = 1.38e-16; /*ボルツマン定数*/
    dvfrag = 3000; //衝突破壊速度(cm/s)
    Q = 1 ; //toomre
    a1 = 10;	//初期位置(au)
    c4 = 1e-7;	//ダスト初期半径(cm)

 /*scanf部分*/
    x=a1;
    y=0;
    z=0;
    r=sqrt(x*x+y*y+z*z);       /* 星間距離 */
    sin_cur = y/r; //sinθ
    count = 0;  //周回数

 
    rp=c4;
 /*scanf終わり*/
 /*ファイル出力*/
    FILE *fp;
    if((fp=fopen("./Q_d/10_1.d","w"))==NULL)
    {
        printf("Cannot open output file \n");
        return 0;
    }
 /*出力終わり*/

    vx=0;
    vy=sqrt(1/x);
    vz=0;
    v=sqrt(vx*vx+vy*vy+vz*vz);
    
    istep=0;         /* step数 */
    t=0.0;          /* 時間 */
    dt=1e-9; /*星間半径で調整*/
    tone = sqrt(r*r*r); /*一周時間(年)*/

    fprintf(fp,"# x=%lf vy=%lf rp=%lf rhop=%lf alpha=%lf dvfrag=%lf Q=%lf\n",x,vy,rp,rhop,alpha,dvfrag,Q); //初期パラメータの出力

 fprintf(fp,"# istep,xt,yt,rt,sin_pre,sin_cur,rpt,vxt,vyt,vdust,count,Fx,Fy,Fg,omega,eta,vkep,vgas,vrelx,vrely,dvdg,Tdisk,cs,hg,sigg,pg,lg,nyu,R,sigd,St,hd,pd,md,nd,dvB,dvtulb,dvdd,tgrow,Agl,M,omeg,Cd,Fx,Fy,F,dt1,dt2,dt,tyr,ts1,ts2,delvg,Ret,Dg,tn\n");

 for(istep=0; istep<mxstep; istep++)/*ループ開始*/
    {
        rr = sqrt(x*x + y*y + z*z);/*軌道半径*/
        tt = t;
        t = tt+dt; /*時間(s)*/

        /* 微分方程式への代入 */
        f[0][1]=x;
        f[0][2]=y;
        f[0][3]=z;
        f[0][4]=vx;
        f[0][5]=vy;
        f[0][6]=vz;
        f[0][7]=rp;

        /* ルンゲクッタ法による、誤差を小さくする計算部分 */
        for(i=0; i<=10; i++)
        {
             tt=tt0+aa[i]*dt;
             if(i !=0)
            {
                for(k=1; k<=7; k++)
                {
                    f[i][k]=f[0][k];
                }
                for(j=0; j<=(i-1); j++)
                {
                    for(k=1; k<=7; k++)
                    {
                        f[i][k]=f[i][k]+bb[i][j]*dfdt[j][k]*dt;
                    }
                }
            }

            /*微分方程式への代入部分*/
            dfdt[i][1]=f[i][4];
            dfdt[i][2]=f[i][5];
            dfdt[i][3]=f[i][6];

            xt=f[i][1];
            yt=f[i][2];
            zt=f[i][3];

            vxt=f[i][4];
            vyt=f[i][5];
            vzt=f[i][6];

            dfdt[i][7]=f[i][7];/*ダスト成長部分*/
            rpt=f[i][7];
    
            rt=sqrt(xt*xt+yt*yt+zt*zt);       /* 星間距離 */
            //前のステップ(i-1)
            x_pre = f[i-1][1];
            y_pre = f[i-1][2];
            z_pre = f[i-1][3];
            r_pre = sqrt(x_pre*x_pre+y_pre*y_pre+z_pre*z_pre);
            sin_pre = y_pre/r_pre;
            sin_cur = yt/rt;
            if(sin_pre*sin_cur<=0){
                count =count+0.5;   //周回のカウント
            }


      Fg = 1/(rt*rt); /*中心力*/
      Fgx = Fg*xt/rt;
      Fgy = Fg*yt/rt;

      omega = pow(rt,-3/2);   /*角速度*/
      vdust = sqrt(vxt*vxt+vyt*vyt); /*ダスト速度*/
      eta = 1.8e-3*sqrt(rt);    /*ズレ*/
      vkep = sqrt(1/rt);/*ケプラー速度*/
      vgas = (1-eta)*vkep; /*ガス速度*/

     /*ダストガス相対速度*/
      vrelx = vxt-(-vgas*yt/rt);
      vrely = vyt-vgas*xt/rt;
      vrel = sqrt(vrelx * vrelx + vrely * vrely);
      dvdg = vrel; /*ダストガス相対速度*/


/*円盤モデル*/
    Tdisk = 150*pow(rt,-3/7); /*L=L**/
    cs = 9.9e4*sqrt(Tdisk/280);     /*音速(cm/s),μ=2.34*/
    hg = yr*cs/omega; /*円盤のスケールハイト(cm)*/
    
    sigg = cs*omega/(pi*G*Q*yr); /*ガス面密度(g/cm^2)(Q=1)*/
    pg = sigg/(sqrt(2*pi)*hg); /*ガスの質量空間密度(g/cm^3)*/
     
    lg = 3.9e-9/pg;               /*平均自由行程(cm)*/
    nyu = 0.353*pow(8/pi,1/2)*cs*lg;  /*動粘性係数(cm^2/s)*/
    R = 2*(AU/yr)*vkep*rpt*fabs(dvdg)/nyu; /*レイノルズ数*/

    /*円盤ダスト面密度g/cm2*/
      sigd = 0.01*sigg; /*ダストガス1:100*/
      St = pi*rhop*rpt/(2*sigg); /*Stokes数*/
      hd = hg*pow(1+((St/alpha)*(1+2*St)/(1+St)),-0.5);   //ダストスケールハイト(cm),Arakawa2021
      pd = sigd/(sqrt(2*pi)*hd); /*g/cm3 ダスト質量空間密度*/
      md = 4*pi*pow(rpt,3)*rhop/3; /*g ダスト空間質量*/
      nd = pd/md; /*ダスト数密度/cm3*/


        M  = (AU/yr)*fabs(dvdg)/cs; /*マッハ数*/
   /*補完係数の決定*/
        if(R < 2.0e5){
                omeg = 0.4;
           }else{
                omeg= 0.2;
           }
          Cd = pow((0.23*M + pow((24/R + 40/(10+R)),-1)),-1) + (2-omeg)*M/(1.6+M) + omeg; /*ガス抵抗係数(無次元)*/

      Fxx = (-3)*AU*Cd*pg*fabs(vrelx)*vrelx/(8*rpt*rhop);
      Fyy = (-3)*AU*Cd*pg*fabs(vrely)*vrely/(8*rpt*rhop);
      Fx=Fxx;
      Fy=Fyy;
      F=sqrt(Fx*Fx+Fy*Fy);  /*ガス抵抗力*/


     /*ダストダスト相対速度(多粒子と仮定)*/
     dvB = (yr/AU)*4*sqrt(kB*Tdisk/pi*md);  /*ブラウン運動(全ダスト質量同等)*/
     ts = (1/Msun)*md*dvdg/F;   //stopping time
     ts1 = ts;
     ts2 = 0.5*ts1;

     delvg = (yr/AU)*sqrt(alpha)*cs;    //渦速度
     Dg = delvg*delvg/omega;    //拡散係数
     Ret = (AU*AU/yr)*Dg/nyu;   //レイノルズ数(乱流)
     tn = pow(Ret,-0.5)*omega;  //回転時間

     if(ts1 < tn){
             dvtulb = delvg*pow(Ret,0.25)*omega*fabs(ts1-ts2);
        }else if(tn < ts1 && ts1 < 1/omega){
                dvtulb = 1.5*delvg*sqrt(omega*ts1);
		}else if(1/omega < ts1){
                        dvtulb = delvg*sqrt((1/(1+omega*ts1))+(1/(1+omega*ts2)));
     }  //乱流(Ormel and Cuzzi 2007)
      dvdd = dvB + dvtulb;


      tgrow = 1/(pi*nd*dvdd*rpt*rpt*AU);    /*成長タイムスケール(無次元)*/
      num1 = 1;
      num2 = -(log(dvdd/dvfrag)/log(5));
    if(num1 < num2){
        Agl = num1;
     }else if(num1 > num2){
        Agl = num2;
    }

      
      /* 運動方程式記述部分 */
      dfdt[i][4] = -Fgx + Fx;
      dfdt[i][5] = -Fgy + Fy;
      dfdt[i][6] = -1/(rt*rt*rt)*zt;
      dfdt[i][7] = Agl * rpt / (3 * tgrow); /*ダスト進化*/	

     /*force test*/
      Ftestx = Fgx/Fx;
      Ftesty = Fgy/Fy;
            }

	/* ルンゲクッタによる誤差の終了 */
 
        for(j=0; j<=10; j++)
        {
             for(k=1; k<=7; k++)
            {
                f[0][k]=f[0][k]+cc[j]*dfdt[j][k]*dt;
            }
        }

        /*微分方程式を解いた結果を代入する部分 */
        x=f[0][1];
        y=f[0][2];
        z=f[0][3];
        vx=f[0][4];
        vy=f[0][5];
        vz=f[0][6];
        rp=f[0][7];

/*タイムスケールの選定*/
    dt1 = 1e-9*rt;
    dt2 = 0.01*tgrow;    /*ダスト成長タイムスケール*/

if(dt1 < dt2){
   dt = dt1;
}else if(dt1 > dt2){
   dt = dt2;
}
    tyr = t; /*年に変換*/

/*終了条件*/
	if (rt > 5000){
		printf("飛んでった\n");
	      	break;
	}
	if (rt < 0.35){
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
}
*/
    wstep = 1e6*rt;             /* stepをいくつ飛ばすか */
    mstep = 1e7*rt;
    


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
}/*mainに対応*/



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

