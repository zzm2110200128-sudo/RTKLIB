/*------------------------------------------------------------------------------
* ppp.c : precise point positioning
*
*          Copyright (C) 2010-2020 by T.TAKASU, All rights reserved.
*
* options : -DIERS_MODEL  use IERS tide model
*           -DOUTSTAT_AMB output ambiguity parameters to solution status
*
* references :
*    [1] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*    [2] D.D.McCarthy and G.Petit, IERS Technical Note 32, IERS Conventions
*        2003, November 2003
*    [3] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*        Space Technology Library, 2004
*    [4] J.Kouba, A Guide to using International GNSS Service (IGS) products,
*        May 2009
*    [5] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
*        Code Biases, URA
*    [6] MacMillan et al., Atmospheric gradients and the VLBI terrestrial and
*        celestial reference frames, Geophys. Res. Let., 1997
*    [7] G.Petit and B.Luzum (eds), IERS Technical Note No. 36, IERS
*         Conventions (2010), 2010
*    [8] J.Kouba, A simplified yaw-attitude model for eclipsing GPS satellites,
*        GPS Solutions, 13:1-12, 2009
*    [9] F.Dilssner, GPS IIF-1 satellite antenna phase center and attitude
*        modeling, InsideGNSS, September, 2010
*    [10] F.Dilssner, The GLONASS-M satellite yaw-attitude model, Advances in
*        Space Research, 2010
*    [11] IGS MGEX (http://igs.org/mgex)
*
* version : $Revision:$ $Date:$
* history : 2010/07/20 1.0  new
*                           added api:
*                               tidedisp()
*           2010/12/11 1.1  enable exclusion of eclipsing satellite
*           2012/02/01 1.2  add gps-glonass h/w bias correction
*                           move windupcorr() to rtkcmn.c
*           2013/03/11 1.3  add otl and pole tides corrections
*                           involve iers model with -DIERS_MODEL
*                           change initial variances
*                           suppress acos domain error
*           2013/09/01 1.4  pole tide model by iers 2010
*                           add mode of ionosphere model off
*           2014/05/23 1.5  add output of trop gradient in solution status
*           2014/10/13 1.6  fix bug on P0(a[3]) computation in tide_oload()
*                           fix bug on m2 computation in tide_pole()
*           2015/03/19 1.7  fix bug on ionosphere correction for GLO and BDS
*           2015/05/10 1.8  add function to detect slip by MW-LC jump
*                           fix ppp solution problem with large clock variance
*           2015/06/08 1.9  add precise satellite yaw-models
*                           cope with day-boundary problem of satellite clock
*           2015/07/31 1.10 fix bug on nan-solution without glonass nav-data
*                           pppoutsolsat() -> pppoutstat()
*           2015/11/13 1.11 add L5-receiver-dcb estimation
*                           merge post-residual validation by rnx2rtkp_test
*                           support support option opt->pppopt=-GAP_RESION=nnnn
*           2016/01/22 1.12 delete support for yaw-model bug
*                           add support for ura of ephemeris
*           2018/10/10 1.13 support api change of satexclude()
*           2020/11/30 1.14 use sat2freq() to get carrier frequency
*                           use E1-E5b for Galileo iono-free LC
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define MAX_ITER    8               /* max number of iterations */
#define MAX_STD_FIX 0.15            /* max std-dev (3d) to fix solution */
#define MIN_NSAT_SOL 4              /* min satellite number for solution */
#define THRES_REJECT 4.0            /* reject threshold of posfit-res (sigma) */

#define THRES_MW_JUMP 10.0

/* E1 实验：手机码方差改用实测 C/N0 曲线（0=关闭与基线一致，1=启用）----*/
#define PPP_CODE_CN0_VARERR 0 /* E1 负结果实验已记录（见 doc）；默认关闭，置 1 可复现 */
#define CN0_VAR_A          79.4  /* σ=a*10^(-b*snr) 的系数：按开发集 B 态 IF(Pc) 验后码残差标定 */
#define CN0_VAR_B          0.045 /* 每 dB-Hz 的对数斜率（候选，参数层再扫描） */
#define CN0_VAR_MIN_SIGMA  0.25  /* 码噪声下限候选（m），冻结前做敏感性检查 */
#define CN0_VAR_MAX_SIGMA  12.0  /* 上限，防止 C/N0 外推时方差爆炸 */

#define VAR_POS     SQR(60.0)       /* init variance receiver position (m^2) */
#define VAR_VEL     SQR(10.0)       /* init variance of receiver vel ((m/s)^2) */
#define VAR_ACC     SQR(10.0)       /* init variance of receiver acc ((m/ss)^2) */
#define VAR_CLK     SQR(60.0)       /* init variance receiver clock (m^2) */
#define VAR_ZTD     SQR( 0.6)       /* init variance ztd (m^2) */
#define VAR_GRA     SQR(0.01)       /* init variance gradient (m^2) */
#define VAR_DCB     SQR(30.0)       /* init variance dcb (m^2) */
#define VAR_BIAS    SQR(60.0)       /* init variance phase-bias (m^2) */
#define VAR_IONO    SQR(60.0)       /* init variance iono-delay */
#define VAR_GLO_IFB SQR( 0.6)       /* variance of glonass ifb */

#define ERR_SAAS    0.3             /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */
#define REL_HUMI    0.7             /* relative humidity for saastamoinen model */
#define GAP_RESION  120             /* default gap to reset ionos parameters (ep) */

#define EFACT_GPS_L5 10.0           /* error factor of GPS/QZS L5 */

#define MUDOT_GPS   (0.00836*D2R)   /* average angular velocity GPS (rad/s) */
#define MUDOT_GLO   (0.00888*D2R)   /* average angular velocity GLO (rad/s) */
#define EPS0_GPS    (13.5*D2R)      /* max shadow crossing angle GPS (rad) */
#define EPS0_GLO    (14.2*D2R)      /* max shadow crossing angle GLO (rad) */
#define T_POSTSHADOW 1800.0         /* post-shadow recovery time (s) */
#define QZS_EC_BETA 20.0            /* max beta angle for qzss Ec (deg) */

/* number and index of states */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics?9:3)
#define NC(opt)     (NSYS)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))
#define NI(opt)     ((opt)->ionoopt==IONOOPT_EST?MAXSAT:0)
#define ND(opt)     ((opt)->nf>=3?1:0)
#define NR(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt)+ND(opt))
#define NB(opt)     (NF(opt)*MAXSAT)
#define NX(opt)     (NR(opt)+NB(opt))
#define IC(s,opt)   (NP(opt)+(s))
#define IT(opt)     (NP(opt)+NC(opt))
#define II(s,opt)   (NP(opt)+NC(opt)+NT(opt)+(s)-1)
#define ID(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt))
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)

/* standard deviation of state -----------------------------------------------*/
static double STD(rtk_t *rtk, int i)
{
    if (rtk->sol.stat==SOLQ_FIX) return SQRT(rtk->Pa[i+i*rtk->nx]);
    return SQRT(rtk->P[i+i*rtk->nx]);
}
/* write solution status for PPP ---------------------------------------------*/
extern int pppoutstat(rtk_t *rtk, char *buff, int level)
{
    ssat_t *ssat;
    double tow,pos[3],vel[3],acc[3],*x;
    int i,j,week;
    char id[8],*p=buff;

    if (!rtk->sol.stat) return 0;

    trace(3,"pppoutstat:\n");

    tow=time2gpst(rtk->sol.time,&week);

    x=rtk->sol.stat==SOLQ_FIX?rtk->xa:rtk->x;

    /* receiver position */
    p+=sprintf(p,"$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
               rtk->sol.stat,x[0],x[1],x[2],STD(rtk,0),STD(rtk,1),STD(rtk,2));

    /* receiver velocity and acceleration */
    if (rtk->opt.dynamics) {
        ecef2pos(rtk->sol.rr,pos);
        ecef2enu(pos,rtk->x+3,vel);
        ecef2enu(pos,rtk->x+6,acc);
        p+=sprintf(p,"$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,"
                   "%.4f,%.5f,%.5f,%.5f\n",week,tow,rtk->sol.stat,vel[0],vel[1],
                   vel[2],acc[0],acc[1],acc[2],0.0,0.0,0.0,0.0,0.0,0.0);
    }
    /* receiver clocks */
    i=IC(0,&rtk->opt);
    p+=sprintf(p,"$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
               week,tow,rtk->sol.stat,1,x[i]*1E9/CLIGHT,x[i+1]*1E9/CLIGHT,
               x[i+2]*1E9/CLIGHT,x[i+3]*1E9/CLIGHT,STD(rtk,i)*1E9/CLIGHT,
               STD(rtk,i+1)*1E9/CLIGHT,STD(rtk,i+2)*1E9/CLIGHT,
               STD(rtk,i+2)*1E9/CLIGHT);

    /* tropospheric parameters */
    if (rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG) {
        i=IT(&rtk->opt);
        p+=sprintf(p,"$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n",week,tow,rtk->sol.stat,
                   1,x[i],STD(rtk,i));
    }
    if (rtk->opt.tropopt==TROPOPT_ESTG) {
        i=IT(&rtk->opt);
        p+=sprintf(p,"$TRPG,%d,%.3f,%d,%d,%.5f,%.5f,%.5f,%.5f\n",week,tow,
                   rtk->sol.stat,1,x[i+1],x[i+2],STD(rtk,i+1),STD(rtk,i+2));
    }
    /* ionosphere parameters */
    if (rtk->opt.ionoopt==IONOOPT_EST) {
        for (i=0;i<MAXSAT;i++) {
            ssat=rtk->ssat+i;
            if (!ssat->vs) continue;
            j=II(i+1,&rtk->opt);
            if (rtk->x[j]==0.0) continue;
            satno2id(i+1,id);
            p+=sprintf(p,"$ION,%d,%.3f,%d,%s,%.1f,%.1f,%.4f,%.4f\n",week,tow,
                       rtk->sol.stat,id,rtk->ssat[i].azel[0]*R2D,
                       rtk->ssat[i].azel[1]*R2D,x[j],STD(rtk,j));
        }
    }
    if (level <= 1) return (int)(p-buff);

    /* Write residuals and status */
    for (int i=0;i<MAXSAT;i++) {
        ssat=rtk->ssat+i;
        if (!ssat->vs) continue;
        satno2id(i+1,id);
        for (int j=0;j<NF(&rtk->opt);j++) {
            int k=IB(i+1,j,&rtk->opt);
            p+=sprintf(p,"$SAT,%d,%.3f,%s,%d,%.1f,%.1f,%.4f,%.4f,%d,%.0f,%d,%d,%d,%u,%u,%u,%.2f,%.6f,%.5f\n",
                       week,tow,id,j+1,ssat->azel[0]*R2D,ssat->azel[1]*R2D,
                       ssat->resp[j],ssat->resc[j],ssat->vsat[j],ssat->snr_rover[j],
                       ssat->fix[j],ssat->slip[j]&(LLI_SLIP|LLI_HALFC),ssat->lock[j],ssat->outc[j],
                       ssat->slipc[j],ssat->rejc[j],k<rtk->nx?rtk->x[k]:0,
                       k<rtk->nx?rtk->P[k+k*rtk->nx]:0,ssat->icbias[j]);
        }
    }
    return (int)(p-buff);
}
/* exclude meas of eclipsing satellite (block IIA) ---------------------------*/
static void testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
    double rsun[3],esun[3],r,ang,erpv[5]={0},cosa;
    int i,j;
    const char *type;

    trace(3,"testeclipse:\n");

    /* unit vector of sun direction (ecef) */
    sunmoonpos(gpst2utc(obs[0].time),erpv,rsun,NULL,NULL);
    normv3(rsun,esun);

    for (i=0;i<n;i++) {
        type=nav->pcvs[obs[i].sat-1].type;

        if ((r=norm(rs+i*6,3))<=0.0) continue;

        /* only block IIA */
        if (*type&&!strstr(type,"BLOCK IIA")) continue;

        /* sun-earth-satellite angle */
        cosa=dot3(rs+i*6,esun)/r;
        cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
        ang=acos(cosa);

        /* test eclipse */
        if (ang<PI/2.0||r*sin(ang)>RE_WGS84) continue;

        char tstr[40];
        trace(3,"eclipsing sat excluded %s sat=%2d\n",time2str(obs[0].time,tstr,0),
              obs[i].sat);

        for (j=0;j<3;j++) rs[j+i*6]=0.0;
    }
}
/* nominal yaw-angle ---------------------------------------------------------*/
static double yaw_nominal(double beta, double mu)
{
    if (fabs(beta)<1E-12&&fabs(mu)<1E-12) return PI;
    return atan2(-tan(beta),sin(mu))+PI;
}
/* yaw-angle of satellite ----------------------------------------------------*/
extern int yaw_angle(int sat, const char *type, int opt, double beta, double mu,
                     double *yaw)
{
    (void)sat;
    (void)type;
    (void)opt;
    *yaw=yaw_nominal(beta,mu);
    return 1;
}
/* satellite attitude model --------------------------------------------------*/
static int sat_yaw(gtime_t time, int sat, const char *type, int opt,
                   const double *rs, double *exs, double *eys)
{
    double rsun[3],ri[6],es[3],esun[3],n[3],p[3],en[3],ep[3],ex[3],E,beta,mu;
    double yaw,cosy,siny,erpv[5]={0};
    int i;

    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,NULL);

    /* beta and orbit angle */
    matcpy(ri,rs,6,1);
    ri[3]-=OMGE*ri[1];
    ri[4]+=OMGE*ri[0];
    cross3(ri,ri+3,n);
    cross3(rsun,n,p);
    if (!normv3(rs,es)||!normv3(rsun,esun)||!normv3(n,en)||
        !normv3(p,ep)) return 0;
    beta=PI/2.0-acos(dot3(esun,en));
    E=acos(dot3(es,ep));
    mu=PI/2.0+(dot3(es,esun)<=0?-E:E);
    if      (mu<-PI/2.0) mu+=2.0*PI;
    else if (mu>=PI/2.0) mu-=2.0*PI;

    /* yaw-angle of satellite */
    if (!yaw_angle(sat,type,opt,beta,mu,&yaw)) return 0;

    /* satellite fixed x,y-vector */
    cross3(en,es,ex);
    cosy=cos(yaw);
    siny=sin(yaw);
    for (i=0;i<3;i++) {
        exs[i]=-siny*en[i]+cosy*ex[i];
        eys[i]=-cosy*en[i]-siny*ex[i];
    }
    return 1;
}
/* phase windup model --------------------------------------------------------*/
static int model_phw(gtime_t time, int sat, const char *type, int opt,
                     const double *rs, const double *rr, double *phw)
{
    double exs[3],eys[3],ek[3],exr[3],eyr[3],eks[3],ekr[3],E[9];
    double dr[3],ds[3],drs[3],r[3],pos[3],cosp,ph;
    int i;

    if (opt<=0) return 1; /* no phase windup */

    /* satellite yaw attitude model */
    if (!sat_yaw(time,sat,type,opt,rs,exs,eys)) return 0;

    /* unit vector satellite to receiver */
    for (i=0;i<3;i++) r[i]=rr[i]-rs[i];
    if (!normv3(r,ek)) return 0;

    /* unit vectors of receiver antenna */
    ecef2pos(rr,pos);
    xyz2enu(pos,E);
    exr[0]= E[1]; exr[1]= E[4]; exr[2]= E[7]; /* x = north */
    eyr[0]=-E[0]; eyr[1]=-E[3]; eyr[2]=-E[6]; /* y = west  */

    /* phase windup effect */
    cross3(ek,eys,eks);
    cross3(ek,eyr,ekr);
    for (i=0;i<3;i++) {
        ds[i]=exs[i]-ek[i]*dot3(ek,exs)-eks[i];
        dr[i]=exr[i]-ek[i]*dot3(ek,exr)+ekr[i];
    }
    cosp=dot3(ds,dr)/norm(ds,3)/norm(dr,3);
    if      (cosp<-1.0) cosp=-1.0;
    else if (cosp> 1.0) cosp= 1.0;
    ph=acos(cosp)/2.0/PI;
    cross3(ds,dr,drs);
    if (dot3(ek,drs)<0.0) ph=-ph;

    *phw=ph+floor(*phw-ph+0.5); /* in cycle */
    return 1;
}
/* measurement error variance ------------------------------------------------*/
static double varerr(int sat, int sys, double el, double snr_rover,
                     int f, const prcopt_t *opt, const obsd_t *obs)
{
    (void)sat;
    double a,b,e;
    double snr_max=opt->err[5];
    double fact=1.0;
    double sinel=sin(el),var;
    int frq,code;
    int e1_code_var=0; /* E1 标志：1 表示本观测已使用 C/N0 曲线码方差（已含 IF 域，行尾不再 ×3） */

    frq=f/2;code=f%2; /* 0=phase, 1=code */
    /* increase variance for pseudoranges */
    if (code) fact=opt->eratio[frq];
    if (fact<=0.0) fact=opt->eratio[0];
    /* adjust variances for constellation */
    switch (sys) {
        case SYS_GPS: fact*=EFACT_GPS;break;
        case SYS_GLO: fact*=EFACT_GLO;break;
        case SYS_GAL: fact*=EFACT_GAL;break;
        case SYS_SBS: fact*=EFACT_SBS;break;
        case SYS_QZS: fact*=EFACT_QZS;break;
        case SYS_CMP: fact*=EFACT_CMP;break;
        case SYS_IRN: fact*=EFACT_IRN;break;
        default:      fact*=EFACT_GPS;break;
    }
    if (sys==SYS_GPS||sys==SYS_QZS) {
        if (frq==2) fact*=EFACT_GPS_L5; /* GPS/QZS L5 error factor */
    }
    /* adjust variance for config parameters */
    a=fact*opt->err[1];  /* base term 常数项 */
    b=fact*opt->err[2];  /* el term 高度角项 */
    /* calculate variance */
    var=(a*a+b*b/sinel/sinel);
    if (opt->err[6]>0) {  /* add SNR term */
        e=fact*opt->err[6];
        var+=e*e*(pow(10,0.1*MAX(snr_max-snr_rover,0)));
    }
    if (opt->err[7]>0.0) {   /* add rcvr stdevs term */
        if (code) var+=SQR(opt->err[7]*obs->Pstd[frq]);
        else var+=SQR(opt->err[7]*obs->Lstd[frq]*0.2);
    }
#if PPP_CODE_CN0_VARERR
    /* E1: 手机伪距方差按实测 C/N0 曲线直接给出（IF/Pc 域，σ 已含 IF 组合
     * 噪声放大，因此行尾的 IFLC ×3 必须跳过，否则重复放大）。仅作用于
     * IFLC + code 观测；C/N0 缺失或异常（<=0 或 >=snr_max）时回退上面的
     * 原公式（含 ×3），与 A 基线行为一致。相位观测不受影响。*/
    if (code&&opt->ionoopt==IONOOPT_IFLC&&snr_rover>0.0&&snr_rover<snr_max) {
        double sigma=CN0_VAR_A*pow(10.0,-CN0_VAR_B*snr_rover); /* 候选曲线 σ(C/N0) */
        sigma=MAX(CN0_VAR_MIN_SIGMA,MIN(CN0_VAR_MAX_SIGMA,sigma)); /* 钳制上下限 */
        var=SQR(sigma);  /* 直接把 IF 域码方差给滤波器 */
        e1_code_var=1;
    }
#endif
    /* FIXME: the scaling factor is not 3 for other signals/constellations than GPS L1/L2 */
    var*=(opt->ionoopt==IONOOPT_IFLC&&!e1_code_var)?SQR(3.0):1.0;
    return var;
}
/* initialize state and covariance -------------------------------------------*/
static inline void initx(rtk_t *rtk, double xi, double var, int i)
{
    int j;
    rtk->x[i]=xi;
    for (j=0;j<rtk->nx;j++) rtk->P[i+j*rtk->nx]=0.0;
    for (j=0;j<rtk->nx;j++) rtk->P[j+i*rtk->nx]=0.0;
    rtk->P[i+i*rtk->nx]=var;
}
/* geometry-free phase measurement -------------------------------------------*/
static double gfmeas(const obsd_t *obs, const nav_t *nav, int f2)
{
    double freq1,freq2;

    freq1=sat2freq(obs->sat,obs->code[0],nav);
    freq2=sat2freq(obs->sat,obs->code[f2],nav);
    if (freq1==0.0||freq2==0.0||obs->L[0]==0.0||obs->L[f2]==0.0) return 0.0;
    return (obs->L[0]/freq1-obs->L[f2]/freq2)*CLIGHT;
}
/* Melbourne-Wubbena linear combination --------------------------------------*/
static double mwmeas(const obsd_t *obs, const nav_t *nav, int f2)
{
    double freq1,freq2;

    freq1=sat2freq(obs->sat,obs->code[0],nav);
    freq2=sat2freq(obs->sat,obs->code[f2],nav);

    if (freq1==0.0||freq2==0.0||obs->L[0]==0.0||obs->L[f2]==0.0||
        obs->P[0]==0.0||obs->P[f2]==0.0) return 0.0;
    return (obs->L[0]-obs->L[f2])*CLIGHT/(freq1-freq2)-
           (freq1*obs->P[0]+freq2*obs->P[f2])/(freq1+freq2);
}
/* antenna corrected measurements --------------------------------------------*/
static void corr_meas(
    const obsd_t *obs,  /* 一颗卫星的原始观测 */
    const nav_t *nav,   /* 频率和码偏差等导航数据 */
    const double *azel, /* 该卫星的方位角和高度角 */
    const prcopt_t *opt,/* PPP 处理选项 */
    const double *dantr,/* 接收机天线改正，单位 m */
    const double *dants,/* 卫星天线改正，单位 m */
    double phw,         /* 相位缠绕改正，单位周 */
    double *L,          /* 输出：各频率改正后载波相位，单位 m */
    double *P,          /* 输出：各频率改正后伪距，单位 m */
    double *Lc,         /* 输出：无电离层组合载波相位，单位 m */
    double *Pc)         /* 输出：无电离层组合伪距，单位 m */
{
    double freq[NFREQ]={0},C1,C2; /* 各观测频率及无电离层组合的两个系数 */
    int i,ix=0,frq2,sys=satsys(obs->sat,NULL); /* 频率下标、组合第二频率和卫星系统；ix 当前未使用 */

    for (i=0;i<opt->nf;i++) { /* 逐频率改正原始观测 */
        L[i]=P[i]=0.0;         /* 先置为无效值，只有通过全部检查才写入 */
        /* skip if low SNR or missing observations */
        freq[i]=sat2freq(obs->sat,obs->code[i],nav); /* 由卫星系统和观测码确定真实载波频率 */
        if (freq[i]==0.0||obs->L[i]==0.0||obs->P[i]==0.0) continue;
        if (testsnr(0,0,azel[1],obs->SNR[i],&opt->snrmask)) continue;

        /* antenna phase center and phase windup correction */
        L[i]=obs->L[i]*CLIGHT/freq[i]-dants[i]-dantr[i]-phw*CLIGHT/freq[i]; /* 周数×波长，并减天线/相位缠绕改正 */
        P[i]=obs->P[i]               -dants[i]-dantr[i]; /* 伪距本来就是 m，只减天线改正 */
        double P_nobias = P[i]; /* 保存应用码偏差前的值，仅供日志对比 */
        if (opt->sateph==EPHOPT_SSRAPC||opt->sateph==EPHOPT_SSRCOM) {
            /* apply SSR correction */
            P[i]-=nav->ssr[obs->sat-1].cbias[obs->code[i]-1]; /* SSR 模式使用实时 SSR 码偏差 */
        }
        else {   /* apply code bias corrections from file */
            P[i]-=code2bias(nav,sys,obs->sat,obs->code[i],1); /* 否则使用导航文件中的绝对码偏差 */
        }
        trace(4,"sys=%d sat=%d frq=%d, P: %.3f->%.3f, dt=%.3f\n",sys,obs->sat,i,P_nobias,P[i],(P[i]-P_nobias)/(1E-9*CLIGHT));
    }
    /* choose freqs for iono-free LC */
    *Lc=*Pc=0.0;                                    /* 默认组合观测无效 */
    frq2=seliflc(opt->nf,satsys(obs->sat,NULL));    /* 按系统选择与第一频率配对的第二频率 */
    if (freq[0]==0.0||freq[frq2]==0.0) return;      /* 任一组合频率无效则不能形成 IFLC */
    C1= SQR(freq[0])/(SQR(freq[0])-SQR(freq[frq2])); /* IFLC 第一频率系数 f1^2/(f1^2-f2^2) */
    C2=-SQR(freq[frq2])/(SQR(freq[0])-SQR(freq[frq2])); /* IFLC 第二频率系数 -f2^2/(f1^2-f2^2) */

    if (L[0]!=0.0&&L[frq2]!=0.0) *Lc=C1*L[0]+C2*L[frq2]; /* 组合相位，一阶电离层项被消除 */
    if (P[0]!=0.0&&P[frq2]!=0.0) *Pc=C1*P[0]+C2*P[frq2]; /* 组合伪距，一阶电离层项被消除 */
    trace(4,"corr_meas: sat=%d f2=%d, Lc=%.3f Pc=%.3f\n",obs->sat,frq2,*Lc,*Pc);
}
/* detect cycle slip by LLI --------------------------------------------------*/
static void detslp_ll(rtk_t *rtk, const obsd_t *obs, int n)
{
    int i,j,nf=rtk->opt.nf;

    trace(3,"detslp_ll: n=%d\n",n);

    if (nf > NFREQ) nf = NFREQ; // Quieten compiler warnings on slip[] write.
    for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<nf;j++) {
        if (obs[i].L[j]==0.0||!(obs[i].LLI[j]&(LLI_SLIP|LLI_HALFC))) continue;

        trace(3,"detslp_ll: slip detected sat=%2d f=%d\n",obs[i].sat,j+1);

        rtk->ssat[obs[i].sat-1].slip[j]=LLI_SLIP;
    }
}
/* detect cycle slip by geometry free phase jump -----------------------------*/
static void detslp_gf(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    double gf0,gf1;
    int i,k,sat;

    trace(4,"detslp_gf: n=%d\n",n);

    if (rtk->opt.thresslip==0) return;  /* return if check disabled */
    for (i=0;i<n&&i<MAXOBS;i++) {
        sat=obs[i].sat;
        for (k=1;k<rtk->opt.nf;k++) {
            /* skip check if slip already detected */
            if (rtk->ssat[sat-1].slip[k]&LLI_SLIP) continue;
            /* calc SD geomotry free LC of phase between freq0 and freqk */
            if ((gf1=gfmeas(obs+i,nav,k))==0.0) continue;

            gf0=rtk->ssat[sat-1].gf[k-1];    /* retrieve previous gf */
            rtk->ssat[sat-1].gf[k-1]=gf1;    /* save current gf for next epoch */

            if (gf0!=0.0&&fabs(gf1-gf0)>rtk->opt.thresslip) {
                rtk->ssat[sat-1].slip[0]|=LLI_SLIP;
                rtk->ssat[sat-1].slip[k]|=LLI_SLIP;
                trace(3,"slip detected GF jump (sat=%2d L1-L%d dGF=%.3f)\n",
                    sat,k+1,gf0-gf1);
            }
        }
    }
}
/* detect slip by Melbourne-Wubbena linear combination jump ------------------*/
static void detslp_mw(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    double mw0,mw1;
    int i,j,k,sat;

    trace(4,"detslp_mw: n=%d\n",n);

    for (i=0;i<n&&i<MAXOBS;i++) {
        sat=obs[i].sat;
        for (k=1;k<rtk->opt.nf;k++) {
            /* skip check if slip already detected */
            if (rtk->ssat[sat-1].slip[k]&LLI_SLIP) continue;
            /* calc MW LC of phase between freq0 and freqk */
            if ((mw1=mwmeas(obs+i,nav,k))==0.0) continue;

            mw0=rtk->ssat[sat-1].mw[k-1];    /* retrieve previous mw */
            rtk->ssat[sat-1].mw[k-1]=mw1;    /* save current mw for next epoch */

            if (mw0!=0.0&&fabs(mw1-mw0)>THRES_MW_JUMP) {
                rtk->ssat[sat-1].slip[0]|=LLI_SLIP;
                rtk->ssat[sat-1].slip[k]|=LLI_SLIP;
                trace(3,"slip detected MW jump (sat=%2d L1-L%d dMW=%.3f)\n",
                    sat,k+1,mw0-mw1);
            }
        }
    }
}
/* temporal update of position -----------------------------------------------*/
static void udpos_ppp(rtk_t *rtk) /* 根据定位模式，初始化或预测 rtk 中的位置、速度和加速度状态 */
{
    double *F;            /* 状态转移矩阵：描述位置、速度、加速度怎样随时间传播 */
    double *P;            /* 从 rtk->P 中取出的有效状态协方差子矩阵 */
    double *FP;           /* 计算 F×P 时使用的临时矩阵 */
    double *x;            /* 从 rtk->x 中取出的有效状态子向量 */
    double *xp;           /* 经过 F×x 得到的预测状态子向量,得到[新位置,新速度,新加速度]*/
    double pos[3];        /* 接收机的大地纬度、经度和高度，用于坐标系转换 */
    double Q[9]={0};      /* 当地坐标系中的 3×3 加速度过程噪声协方差，初始清零 */
    double Qv[9];         /* 将 Q 转到地心地固坐标系后得到的 3×3 协方差 */
    double var=0.0;       /* X、Y、Z 三个位置状态方差的平均值 */
    int i,j;              /* 循环计数器 */
    int *ix;              /* 保存参与本次预测的有效状态在完整状态向量中的下标 */
    int nx;               /* 参与本次预测的有效状态数量 */

    trace(3,"udpos_ppp:\n"); /* 在调试日志中记录程序进入了位置状态更新函数 */

  /* 1.fixed mode 固定坐标模式： 使用用户给定坐标，几乎完全不允许变化   */
    if (rtk->opt.mode==PMODE_PPP_FIXED) { /* 固定坐标模式：接收机坐标已知，不在 PPP 中估计位置 */
        for (i=0;i<3;i++)                 /* i=0、1、2 分别对应 ECEF 的 X、Y、Z 坐标状态 */
            initx(rtk,rtk->opt.ru[i],1E-8,i); /* 用用户给定坐标初始化状态，并赋予极小方差 */
        return;                           /* 固定坐标已设置完成，无需执行后面的位置预测 */
    }
    /* initialize position for first epoch */
    if (norm(rtk->x,3)<=0.0) {       /* X、Y、Z 状态全为 0：认为这是尚未初始化的第一个历元 */
        for (i=0;i<3;i++)            /* 初始化 X、Y、Z 三个位置状态 */
            initx(rtk,rtk->sol.rr[i],VAR_POS,i); /* 初值取单点定位坐标，初始标准差为 60 m */
        if (rtk->opt.dynamics) {      /* 如果启用了速度和加速度动力学模型 */
            for (i=3;i<6;i++)        /* 初始化 X、Y、Z 三个速度状态 */
                initx(rtk,rtk->sol.rr[i],VAR_VEL,i); /* 初值取单点定位速度，初始标准差为 10 m/s */
            for (i=6;i<9;i++)        /* 初始化 X、Y、Z 三个加速度状态 */
                initx(rtk,1E-6,VAR_ACC,i); /* 用极小非零值作为初值，初始标准差为 10 m/s^2 */
        }
    }
    /* 2.static ppp mode 静态PPP：保留上一历元PPP位置，只缓慢增加方差  */
    if (rtk->opt.mode==PMODE_PPP_STATIC) { /* 静态 PPP：位置值保持不变，只更新其不确定程度 */
        for (i=0;i<3;i++) {                /* 依次处理 X、Y、Z 三个位置状态 */
            rtk->P[i*(1+rtk->nx)]+=SQR(rtk->opt.prn[5])*fabs(rtk->tt); /* 位置方差增加：过程噪声平方×历元间隔 */
        }
        return;                            /* 静态模式无需执行后面的运动状态预测 */
    }
    /* 3.kinematic mode without dynamics 动态PPP但未启用动力学：每个历元用当前单点定位坐标重新初始化位置  */
    if (!rtk->opt.dynamics) {       /* 动态定位但未启用速度、加速度模型，无法从上一历元预测位置 */
        for (i=0;i<3;i++) {         /* 依次重新初始化 X、Y、Z 三个位置状态 */
            initx(rtk,rtk->sol.rr[i],VAR_POS,i); /* 每个历元都用当前单点定位坐标作为 PPP 位置初值 */
        }
        return;                     /* 初始化完成，不执行后面的动力学状态转移 */
    }
    /* check variance of estimated position */
    for (i=0;i<3;i++)                      /* 依次读取 X、Y、Z 的协方差对角元素 */
        var+=rtk->P[i+i*rtk->nx];          /* 将三个位置状态的方差累加到 var */
    var/=3.0;                              /* 除以 3，得到三维位置方差的平均值 */

    if (var>VAR_POS) {                    /* 平均位置方差超过 60^2 m^2：认为当前动力学状态已不可靠 */
        /* reset position with large variance */
        for (i=0;i<3;i++)                 /* 重置 X、Y、Z 位置状态 */
            initx(rtk,rtk->sol.rr[i],VAR_POS,i); /* 位置重新取当前单点定位坐标 */
        for (i=3;i<6;i++)                 /* 重置 X、Y、Z 速度状态 */
            initx(rtk,rtk->sol.rr[i],VAR_VEL,i); /* 速度重新取当前单点定位速度 */
        for (i=6;i<9;i++)                 /* 重置 X、Y、Z 加速度状态 */
            initx(rtk,1E-6,VAR_ACC,i);     /* 加速度重新设为接近 0 的非零初值 */
        trace(2,"reset rtk position due to large variance: var=%.3f\n",var); /* 记录重置原因和平均方差 */
        return;                            /* 重置已经代替本次预测，直接结束位置更新 */
    }
    /* generate valid state index */
    ix=imat(rtk->nx,1);                 /* 最多按完整状态数分配一个整数下标数组 */
    for (i=nx=0;i<rtk->nx;i++) {       /* i 遍历完整状态；nx 从 0 开始统计有效状态数 */
        if (i<9||(rtk->x[i]!=0.0&&rtk->P[i+i*rtk->nx]>0.0)) /* 前 9 个运动状态总是保留；其他状态须已初始化且方差有效 */
            ix[nx++]=i;                 /* 记录该状态在完整向量中的下标，然后有效状态数加 1 */
    }
    /* state transition of position/velocity/acceleration */
    F=eye(nx);           /* 生成 nx×nx 单位矩阵，作为状态转移矩阵 F 的初始形式 */
    P=mat(nx,nx);        /* 为有效状态的 nx×nx 协方差子矩阵分配空间 */
    FP=mat(nx,nx);       /* 为中间乘积 F×P 分配 nx×nx 空间 */
    x=mat(nx,1);         /* 为预测前的有效状态子向量分配 nx×1 空间 */
    xp=mat(nx,1);        /* 为预测后的有效状态子向量分配 nx×1 空间 */

    for (i=0;i<6;i++) {             /* i=0..2 设置速度对位置的影响；i=3..5 设置加速度对速度的影响 */
        F[i+(i+3)*nx]=rtk->tt;      /* 对应系数均为历元时间间隔 dt，使新状态加入“变化率×dt” */
    }
    /* include accel terms if filter is converged */
    if (var<rtk->opt.thresar[1]) {       /* 平均位置方差低于阈值：认为位置已较稳定，可使用加速度项 */
        for (i=0;i<3;i++) {              /* 分别设置 X、Y、Z 加速度对对应位置的影响 */
            F[i+(i+6)*nx]=SQR(rtk->tt)/2.0; /* 系数为 dt^2/2，对应位移公式 0.5×加速度×dt^2 */
        }
    }
    else trace(3,"pos var too high for accel term: %.4f,%.4f\n", var,rtk->opt.thresar[1]); /* 方差过大时跳过位置加速度项并记录日志 */
    for (i=0;i<nx;i++) {                 /* i 遍历下标表中的每一个有效状态 */
        x[i]=rtk->x[ix[i]];              /* 按 ix[i] 从完整状态 rtk->x 取出第 i 个有效状态 */
        for (j=0;j<nx;j++) {             /* j 遍历每一个有效状态，提取状态间协方差 */
            P[i+j*nx]=rtk->P[ix[i]+ix[j]*rtk->nx]; /* 从完整 rtk->P 取出对应行列，组成协方差子矩阵 P */
        }
    }
    /* x=F*x, P=F*P*F+Q */
    matmul("NN",nx,1,nx,F,x,xp);   /* 计算 xp=F×x，得到传播到当前历元的预测状态 */
    matmul("NN",nx,nx,nx,F,P,FP);  /* 第一步计算中间矩阵 FP=F×P */
    matmul("NT",nx,nx,nx,FP,F,P);  /* 第二步计算 P=FP×F^T，即 P=F×P×F^T */

    for (i=0;i<nx;i++) {                 /* i 遍历每一个有效预测状态 */
        rtk->x[ix[i]]=xp[i];             /* 按 ix[i] 把预测状态写回完整状态向量 rtk->x */
        for (j=0;j<nx;j++) {             /* j 遍历有效状态协方差的每一列 */
            rtk->P[ix[i]+ix[j]*rtk->nx]=P[i+j*nx]; /* 将协方差子矩阵对应元素写回完整 rtk->P */
        }
    }
    /* process noise added to only acceleration */
    Q[0]=Q[4]=SQR(rtk->opt.prn[3])*fabs(rtk->tt); /* 设置东、北两个水平方向的加速度过程噪声方差 */
    Q[8]=SQR(rtk->opt.prn[4])*fabs(rtk->tt);      /* 设置天向（垂直方向）的加速度过程噪声方差 */
    ecef2pos(rtk->x,pos);                         /* 将预测后的 ECEF 坐标转成纬度、经度和高度 */
    covecef(pos,Q,Qv);                            /* 将当地东、北、天协方差 Q 转成 ECEF 协方差 Qv */
    for (i=0;i<3;i++)                             /* i 遍历 X、Y、Z 加速度协方差的行 */
        for (j=0;j<3;j++) {                       /* j 遍历 X、Y、Z 加速度协方差的列 */
            rtk->P[i+6+(j+6)*rtk->nx]+=Qv[i+j*3]; /* 把过程噪声加到完整协方差中第 6～8 项的加速度块 */
        }
    free(ix); /* 释放有效状态下标数组 */
    free(F);  /* 释放状态转移矩阵 */
    free(P);  /* 释放有效状态协方差子矩阵 */
    free(FP); /* 释放矩阵乘法中间结果 */
    free(x);  /* 释放预测前的有效状态子向量 */
    free(xp); /* 释放预测后的有效状态子向量 */
}
/* temporal update of clock --------------------------------------------------*/
static void udclk_ppp(rtk_t *rtk) /* 每个历元重新初始化接收机钟差和系统间钟差状态 */
{
    double dtr; /* 当前正在处理的接收机钟差临时值，单位为秒 */
    int i;      /* 卫星系统钟差状态的循环下标 */

    trace(3,"udclk_ppp:\n"); /* 在调试日志中记录程序进入了钟差状态更新函数 */

    /* initialize every epoch for clock (white noise) */
    for (i=0;i<NSYS;i++) {                         /* 逐个初始化本版本支持的卫星系统钟差状态 */
        if (rtk->opt.sateph==EPHOPT_PREC) {        /* 使用精密星历/钟差产品时，其时间基准按 GPST 处理 */
            /* time of prec ephemeris is based gpst */
            /* neglect receiver inter-system bias  */
            dtr=rtk->sol.dtr[0];                   /* 各系统先使用单点定位得到的公共接收机钟差 */
        }
        else {
            dtr=i==0?rtk->sol.dtr[0]:rtk->sol.dtr[0]+rtk->sol.dtr[i]; /* GPS 用公共钟差；其他系统再加相对 GPS 的系统间偏差 */
        }
        initx(rtk,CLIGHT*dtr,VAR_CLK,IC(i,&rtk->opt)); /* 秒乘光速转成米，并写入第 i 个钟差状态 */
    }
}
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop_ppp(rtk_t *rtk) /* 初始化或预测 PPP 的对流层延迟及其水平梯度状态 */
{
    double pos[3];                    /* 接收机的大地纬度、经度和高度 */
    double azel[]={0.0,PI/2.0};       /* 初始化模型使用的天顶方向：方位角 0、高度角 90 度 */
    double ztd;                       /* 天顶方向的对流层总延迟，单位为米 */
    double var;                       /* 对流层延迟初值的方差 */
    int i=IT(&rtk->opt);              /* 对流层状态在完整 PPP 状态向量中的起始下标 */
    int j;                            /* 初始化或更新水平梯度状态时使用的循环下标 */

    trace(3,"udtrop_ppp:\n"); /* 在调试日志中记录程序进入了对流层状态更新函数 */

    if (rtk->x[i]==0.0) {                    /* ZTD 状态值为 0：认为对流层状态尚未初始化 */
        ecef2pos(rtk->sol.rr,pos);            /* 将当前接收机 ECEF 坐标转换为纬度、经度和高度 */
        ztd=sbstropcorr(rtk->sol.time,pos,azel,&var); /* 按当前时间、位置和天顶方向计算 ZTD 初值及方差 */
        initx(rtk,ztd,var,i);                 /* 把 ZTD 初值和方差写入第 i 个 PPP 状态 */

        if (rtk->opt.tropopt>=TROPOPT_ESTG) { /* 若配置还要求估计南北、东西两个水平梯度 */
            for (j=i+1;j<i+3;j++)             /* 梯度状态紧跟在 ZTD 状态后，共两个 */
                initx(rtk,1E-6,VAR_GRA,j);     /* 用接近 0 的非零值和预设方差初始化梯度 */
        }
    }
    else {                                    /* ZTD 已经初始化：保留状态值，只传播不确定性 */
        rtk->P[i+i*rtk->nx]+=SQR(rtk->opt.prn[2])*fabs(rtk->tt); /* 给 ZTD 方差增加过程噪声 */

        if (rtk->opt.tropopt>=TROPOPT_ESTG) { /* 如果同时估计两个水平梯度 */
            for (j=i+1;j<i+3;j++) {           /* 依次更新南北、东西梯度的不确定性 */
                rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[2]*0.1)*fabs(rtk->tt); /* 梯度过程噪声取 ZTD 的 0.1 倍 */
            }
        }
    }
}
/* temporal update of ionospheric parameters ---------------------------------*/
static void udiono_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav) /* 初始化或预测每颗卫星的垂直电离层延迟状态 */
{
    double freq1,freq2;             /* 第一频率和选定第二频率的载波频率，单位 Hz */
    double ion;                     /* 当前卫星的垂直 L1 等效电离层延迟，单位 m */
    double sinel;                   /* 卫星高度角的正弦值，用于缩放过程噪声 */
    double pos[3];                  /* 接收机纬度、经度和高度 */
    double *azel;                   /* 指向当前卫星方位角、高度角的指针 */
    double var=VAR_IONO;            /* 电离层状态初始方差 */
    char *p;                        /* 在 PPP 扩展选项字符串中查找参数的位置 */
    int i,j;                        /* 观测/卫星循环下标和状态下标 */
    int f2;                         /* 与第一频率组成双频组合的第二频率下标 */
    int gap_resion=GAP_RESION;      /* 电离层状态允许连续缺测的默认历元数 */
    int sat;                        /* 当前观测对应的 RTKLIB 内部卫星编号 */

    trace(3,"udiono_ppp:\n"); /* 记录进入电离层状态更新函数 */

    if ((p=strstr(rtk->opt.pppopt,"-GAP_RESION="))) { /* 检查用户是否在扩展选项中覆盖缺测阈值 */
        sscanf(p,"-GAP_RESION=%d",&gap_resion);       /* 读取用户设置的最大连续缺测历元数 */
    }
    /* reset ionosphere delay estimate if outage too long */
    for (i=0;i<MAXSAT;i++) {                         /* 遍历所有可能的卫星 */
        j=II(i+1,&rtk->opt);                         /* 计算该卫星电离层状态在 rtk->x 中的下标 */
        if (rtk->x[j]!=0.0&&(int)rtk->ssat[i].outc[0]>gap_resion) { /* 状态存在但连续缺测过久 */
            rtk->x[j]=0.0;                           /* 将状态置 0，等待卫星恢复后重新初始化 */
        }
    }
    /* reset ionosphere states if VTEC corrections just became available */
    if (rtk->vtec_used==0&&nav->vtec.nlay>0) {       /* 此前未用 VTEC，而当前导航数据首次提供了 VTEC */
        for (i=0;i<MAXSAT;i++) {                     /* 遍历全部卫星电离层状态 */
            rtk->x[II(i+1,&rtk->opt)]=0.0;           /* 清除旧初值，使其下一步改用 VTEC 重新初始化 */
        }
        rtk->vtec_used=1;                            /* 记住本次以后已经使用过 VTEC 系数 */
    }
    for (i=0;i<n;i++) {                              /* 遍历当前历元的每条卫星观测 */
        sat=obs[i].sat;                               /* 取得当前观测的卫星编号 */
        j=II(sat,&rtk->opt);                          /* 取得该卫星电离层状态下标 */
        if (rtk->x[j]==0.0&&(int)rtk->ssat[i].outc[0]<=gap_resion) { /* 状态尚未初始化且缺测未超限 */
            /* initialize ionosphere delay estimates if zero */
            ecef2pos(rtk->sol.rr,pos);                /* 获取接收机纬度、经度和高度 */
            azel=rtk->ssat[sat-1].azel;               /* 获取当前卫星方位角和高度角 */
            f2=seliflc(rtk->opt.nf,satsys(sat,NULL)); /* 按系统和频率数选择第二频率 */
            if (testsnr(0,0,azel[1],obs[i].SNR[0],&rtk->opt.snrmask)) continue; /* 第一频率信噪比不合格则跳过 */
            freq1=sat2freq(sat,obs[i].code[0],nav);   /* 根据卫星和观测码取得第一载波频率 */
            freq2=sat2freq(sat,obs[i].code[f2],nav);  /* 取得选定第二载波频率 */
            if (freq1==0.0) continue;                 /* 无法得到第一频率时不能初始化 */
            if (nav->vtec.nlay>0) {  /* use VTEC if corrections available */
                ionvtec(obs[i].time,nav,pos,azel,freq1,&ion,&var); /* 优先由 VTEC 模型得到斜向电离层延迟和方差 */
                if (var==0.0) continue;                         /* VTEC 结果无有效方差则跳过 */
            } else {
                if (obs[i].P[0]==0.0||obs[i].P[f2]==0.0||freq2==0.0||
                        testsnr(0,f2,azel[1],obs[i].SNR[f2],&rtk->opt.snrmask)) {
                    continue;                         /* 缺少双频伪距/频率或第二频率信噪比差时无法初始化 */
                }
                /* use pseudorange difference adjusted by freq for initial estimate */
                int sys=satsys(sat,NULL);             /* 取得当前卫星所属系统 */
                double P0_corr=obs[i].P[0];           /* 第一频率伪距的可改正副本 */
                double Pf_corr=obs[i].P[f2];          /* 第二频率伪距的可改正副本 */
                if (rtk->opt.sateph==EPHOPT_SSRAPC||rtk->opt.sateph==EPHOPT_SSRCOM) {
                    /* apply SSR correction */
                    P0_corr-=nav->ssr[obs->sat-1].cbias[obs[i].code[0]-1]; /* 从第一频率伪距减去 SSR 码偏差 */
                    Pf_corr-=nav->ssr[obs->sat-1].cbias[obs[i].code[f2]-1]; /* 从第二频率伪距减去 SSR 码偏差 */
                }
                else {   /* apply code bias corrections from file */
                    P0_corr-=code2bias(nav,sys,sat,obs[i].code[0],1);  /* 应用文件提供的第一频率码偏差 */
                    Pf_corr-=code2bias(nav,sys,sat,obs[i].code[f2],1); /* 应用文件提供的第二频率码偏差 */
                }
                ion=(P0_corr-Pf_corr)/(SQR(FREQL1/freq1)-SQR(FREQL1/freq2)); /* 用双频伪距差估计斜向 L1 等效电离层延迟 */
                trace(3,"P1=%.3f P2=%.3f frq1=%.1f frq2=%.1f\n",obs[i].P[0],obs[i].P[f2],freq1,freq2);
                var=VAR_IONO;                         /* 双频伪距法使用预设电离层初始方差 */
            }
            /* adjust delay estimate by path length */
            ion/=ionmapf(pos,azel);                    /* 除以映射函数，将斜向延迟换成垂直延迟状态 */
            initx(rtk,ion,var,j);                       /* 写入该卫星电离层状态及初始方差 */
            trace(3,"ion init: sat=%d ion=%.4f var=%.1f\n",sat,ion,var); /* 记录初始化结果 */
        }
        else { /* temporal update */
            sinel=sin(MAX(rtk->ssat[sat-1].azel[1],5.0*D2R)); /* 高度角最低按 5 度计算，避免过程噪声无限放大 */
            /* update variance of delay state */
            rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[1]/sinel)*fabs(rtk->tt); /* 保留电离层值并按高度角增加过程噪声方差 */
        }
    }
}
/* temporal update of L5-receiver-dcb parameters -----------------------------*/
static void uddcb_ppp(rtk_t *rtk) /* 初始化三频 PPP 使用的接收机第 3 频率码偏差状态 */
{
    int i=ID(&rtk->opt); /* ID() 返回该 DCB 状态在完整状态向量 rtk->x 中的下标 */

    trace(3,"uddcb_ppp:\n"); /* 写入调试日志 */

    if (rtk->x[i]==0.0) {          /* 值为 0 表示这个状态还没有初始化 */
        initx(rtk,1E-6,VAR_DCB,i); /* 以接近 0 m 的值和 VAR_DCB 方差建立该状态 */
    }
}
/* temporal update of phase biases -------------------------------------------*/
static void udbias_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav) /* 更新载波相位偏差（模糊度，单位 m） */
{
    double L[NFREQ],P[NFREQ],Lc,Pc; /* 改正后的各频点相位/伪距，以及无电离层组合相位/伪距 */
    double bias[MAXOBS],offset=0.0; /* 各观测的模糊度初值；offset 为公共相位－码跳变平均量 */
    double freq1,freq2,ion;         /* 载波频率和由双频伪距估计的电离层延迟 */
    double dantr[NFREQ]={0},dants[NFREQ]={0}; /* 此处初始化模糊度时暂不重复加入收、发天线改正 */
    int i,j,k,f,sat;               /* 循环下标、状态下标、频率下标和卫星编号 */
    int slip[MAXOBS]={0};           /* 当前历元各观测在本频率上是否发生周跳 */
    int clk_jump=0;                 /* 是否处在启用日界跳变处理的 GPS 日边界 */

    trace(3,"udbias  : n=%d\n",n); /* 记录本历元观测数 */

    /* handle day-boundary clock jump */
    if (rtk->opt.posopt[5]) { /* 只有打开 posopt[5] 才检查接收机钟在日界处的跳变 */
        clk_jump=ROUND(time2gpst(obs[0].time,NULL)*10)%864000==0; /* 0.1 s 计数恰好落在 86400 s 边界 */
    }
    for (i=0;i<MAXSAT;i++) for (j=0;j<rtk->opt.nf;j++) { /* 先清除上一历元的临时周跳标志 */
        rtk->ssat[i].slip[j]=0;                           /* 后面三种探测方法会重新置位 */
    }
    /* detect cycle slip by LLI */
    detslp_ll(rtk,obs,n); /* 用接收机给出的 LLI 标志探测周跳 */

    /* detect cycle slip by geometry-free phase jump */
    detslp_gf(rtk,obs,n,nav); /* 用无几何相位组合的历元跳变探测周跳 */

    /* detect slip by Melbourne-Wubbena linear combination jump */
    detslp_mw(rtk,obs,n,nav); /* 用 Melbourne-Wubbena 组合的历元跳变探测周跳 */

    for (f=0;f<NF(&rtk->opt);f++) { /* 逐个参与 PPP 的频率/组合处理模糊度 */
        offset=0;                    /* 每个频率重新累计公共跳变量 */
        /* reset phase-bias if expire obs outage counter */
        for (i=0;i<MAXSAT;i++) { /* 对所有卫星增加一次连续缺测计数 */
            if (++rtk->ssat[i].outc[f]>(uint32_t)rtk->opt.maxout||
                rtk->opt.modear==ARMODE_INST||clk_jump) { /* 缺测过久、瞬时 AR 或日界跳变都使旧模糊度失效 */
                initx(rtk,0.0,0.0,IB(i+1,f,&rtk->opt));   /* 将对应模糊度状态及方差清零，等待重建 */
            }
        }
        for (i=k=0;i<n&&i<MAXOBS;i++) { /* 第一遍：计算各卫星模糊度初值并统计公共跳变 */
            sat=obs[i].sat;              /* 当前卫星编号 */
            j=IB(sat,f,&rtk->opt);       /* 当前卫星、当前频率的模糊度状态下标 */
            corr_meas(obs+i,nav,rtk->ssat[sat-1].azel,&rtk->opt,dantr,dants,
                      0.0,L,P,&Lc,&Pc); /* 得到以 m 为单位的相位、伪距和无电离层组合 */

            bias[i]=0.0; /* 0 表示本观测暂时无法形成有效模糊度初值 */

            if (rtk->opt.ionoopt==IONOOPT_IFLC) { /* 无电离层组合模式 */
                bias[i]=Lc-Pc;                    /* 相位－伪距，近似得到组合模糊度（m） */
                int f2=seliflc(rtk->opt.nf,rtk->ssat[sat-1].sys); /* 选出组成 IFLC 的第二频率 */
                slip[i]=rtk->ssat[sat-1].slip[0]||rtk->ssat[sat-1].slip[f2]; /* 任一组合频率周跳即需重置 */
            }
            else if (L[f]!=0.0&&P[f]!=0.0) { /* 非组合模式且本频率相位、伪距都有效 */
                freq1=sat2freq(sat,obs[i].code[0],nav); /* 基准频率 */
                freq2=sat2freq(sat,obs[i].code[f],nav); /* 当前频率 */
                slip[i]=rtk->ssat[sat-1].slip[f];       /* 读取当前频率周跳结果 */
                if (f==0||obs[i].P[0]==0.0||obs[i].P[f]==0.0||freq1==0.0||freq2==0.0)
                    ion=0; /* 第一频率或数据不足时，不由伪距差估计电离层项 */
                else
                    ion=(obs[i].P[0]-obs[i].P[f])/(1.0-SQR(freq1/freq2)); /* 双频伪距差估计当前频率电离层量 */
                bias[i]=L[f]-P[f]+2.0*ion*SQR(freq1/freq2); /* 相位与伪距的电离层符号相反，所以补偿 2I */
            }
            if (rtk->x[j]==0.0||slip[i]||bias[i]==0.0) continue; /* 只用连续、已有状态的卫星估计公共跳变 */

            offset+=bias[i]-rtk->x[j]; /* 新算初值减去滤波器旧值 */
            k++;                       /* 有效参与平均的卫星数 */
        }
        /* correct phase-code jump to ensure phase-code coherence */
        if (k>=2&&fabs(offset/k)>0.0005*CLIGHT) { /* 至少两星且公共差值超过约 0.5 ms 对应距离 */
            for (i=0;i<MAXSAT;i++) {             /* 将同一个公共跳变量加到所有有效模糊度 */
                j=IB(i+1,f,&rtk->opt);            /* 第 i+1 颗卫星、本频率模糊度下标 */
                if (rtk->x[j]!=0.0) rtk->x[j]+=offset/k; /* 保持相位与伪距钟的一致性 */
            }
            char tstr[40];
            trace(2,"phase-code jump corrected: %s n=%2d dt=%12.9fs\n",
                  time2str(rtk->sol.time,tstr,0),k,offset/k/CLIGHT);
        }
        for (i=0;i<n&&i<MAXOBS;i++) { /* 第二遍：预测旧模糊度方差，必要时重新初始化 */
            sat=obs[i].sat;            /* 当前卫星编号 */
            j=IB(sat,f,&rtk->opt);     /* 当前模糊度状态下标 */

            rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[0])*fabs(rtk->tt); /* 随时间给模糊度对角方差加入随机游走噪声 */

            if (bias[i]==0.0||(rtk->x[j]!=0.0&&!slip[i])) continue; /* 无初值或旧状态仍连续时无需重建 */

            /* reinitialize phase-bias if detecting cycle slip */
            initx(rtk,bias[i],VAR_BIAS,IB(sat,f,&rtk->opt)); /* 首次见星或周跳后，以新 bias 和初始方差重建 */
            trace(3,"init bias: sat=%d frq=%d\n", sat,f);

            /* reset fix flags */
            for (k=0;k<MAXSAT;k++) rtk->ambc[sat-1].flags[k]=0; /* 模糊度重置后，旧整数固定关系全部作废 */

            trace(3,"udbias_ppp: sat=%2d bias=%.3f\n",sat,bias[i]);
        }
    }
}
/* temporal update of states --------------------------------------------------*/
static void udstate_ppp(
    rtk_t *rtk,          /* 输入并更新：保存上一历元状态，本函数将其预测到当前历元 */
    const obsd_t *obs,   /* 输入：当前历元的卫星观测值 */
    int n,               /* 输入：当前历元的观测记录数量 */
    const nav_t *nav)    /* 输入：星历、钟差、偏差和电离层等导航数据 */
{
    trace(3,"udstate_ppp: n=%d\n",n); /* 在调试日志中记录本次状态更新的观测记录数量 */

    /* temporal update of position */
    udpos_ppp(rtk); /* 首先初始化或预测接收机的位置、速度和加速度状态 */

    /* temporal update of clock */
    udclk_ppp(rtk); /* 为各卫星系统初始化当前历元接收机钟差状态 */

    /* temporal update of tropospheric parameters */
    if (rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG) { /* 仅估计 ZTD/ZTD+梯度时需要状态预测 */
        udtrop_ppp(rtk); /* 初始化或保留对流层状态，并给方差加入过程噪声 */
    }
    /* temporal update of ionospheric parameters */
    if (rtk->opt.ionoopt==IONOOPT_EST) { /* 仅非组合且把电离层作为未知数时执行 */
        udiono_ppp(rtk,obs,n,nav);       /* 初始化/预测每颗卫星的垂直电离层延迟 */
    }
    /* temporal update of L5-receiver-dcb parameters */
    if (rtk->opt.nf>=3) { /* 第 3 频率伪距进入滤波时才需要接收机 DCB 状态 */
        uddcb_ppp(rtk);    /* 初始化第 3 频率接收机码偏差 */
    }
    /* temporal update of phase-bias */
    udbias_ppp(rtk,obs,n,nav); /* 最后更新各卫星、各频率的载波相位偏差/模糊度 */
}
/* satellite antenna phase center variation ----------------------------------*/
static void satantpcv(const double *rs, const double *rr, const pcv_t *pcv,
                      double *dant)
{
    double ru[3],rz[3],eu[3],ez[3],nadir,cosa;
    int i;

    for (i=0;i<3;i++) {
        ru[i]=rr[i]-rs[i];
        rz[i]=-rs[i];
    }
    if (!normv3(ru,eu)||!normv3(rz,ez)) return;

    cosa=dot3(eu,ez);
    cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
    nadir=acos(cosa);

    antmodel_s(pcv,nadir,dant);
}
/* precise tropospheric model ------------------------------------------------*/
static double trop_model_prec(gtime_t time, const double *pos,
                              const double *azel, const double *x, double *dtdx,
                              double *var)
{
    const double zazel[]={0.0,PI/2.0};
    double zhd,m_h,m_w,cotz,grad_n,grad_e;

    /* zenith hydrostatic delay */
    zhd=tropmodel(time,pos,zazel,0.0);

    /* mapping function */
    m_h=tropmapf(time,pos,azel,&m_w);

    if (azel[1]>0.0) {

        /* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
        cotz=1.0/tan(azel[1]);
        grad_n=m_w*cotz*cos(azel[0]);
        grad_e=m_w*cotz*sin(azel[0]);
        m_w+=grad_n*x[1]+grad_e*x[2];
        dtdx[1]=grad_n*(x[0]-zhd);
        dtdx[2]=grad_e*(x[0]-zhd);
    }
    dtdx[0]=m_w;
    *var=SQR(0.01);
    return m_h*zhd+m_w*(x[0]-zhd);
}
/* tropospheric model ---------------------------------------------------------*/
static int model_trop(gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, const double *x, double *dtdx,
                      const nav_t *nav, double *dtrp, double *var)
{
    (void)nav;
    double trp[3]={0};

    if (opt->tropopt==TROPOPT_SAAS) {
        *dtrp=tropmodel(time,pos,azel,REL_HUMI);
        *var=SQR(ERR_SAAS);
        return 1;
    }
    if (opt->tropopt==TROPOPT_SBAS) {
        *dtrp=sbstropcorr(time,pos,azel,var);
        return 1;
    }
    if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
        matcpy(trp,x+IT(opt),opt->tropopt==TROPOPT_EST?1:3,1);
        *dtrp=trop_model_prec(time,pos,azel,trp,dtdx,var);
        return 1;
    }
    return 0;
}
/* ionospheric model ---------------------------------------------------------*/
static int model_iono(gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, int sat, const double *x,
                      const nav_t *nav, double *dion, double *var)
{
    if (opt->ionoopt==IONOOPT_SBAS) {
        return sbsioncorr(time,nav,pos,azel,dion,var);
    }
    if (opt->ionoopt==IONOOPT_TEC) {
        return iontec(time,nav,pos,azel,1,dion,var);
    }
    if (opt->ionoopt==IONOOPT_BRDC) {
        *dion=ionmodel(time,nav->ion_gps,pos,azel);
        *var=SQR(*dion*ERR_BRDCI);
        return 1;
    }
    if (opt->ionoopt==IONOOPT_EST) {
        /* Estimated delay is a vertical delay, apply the mapping function. */
        *dion=x[II(sat,opt)]*ionmapf(pos,azel);
        *var=0.0;
        return 1;
    }
    if (opt->ionoopt==IONOOPT_IFLC) {
        *dion=*var=0.0;
        return 1;
    }
    return 0;
}
/* phase and code residuals --------------------------------------------------*/
static int ppp_res(
    int post,              /* -1/0：滤波前残差；>0：滤波后残差及粗差复检 */
    const obsd_t *obs,     /* 当前历元观测数组 */
    int n,                 /* 观测数组中的卫星数 */
    const double *rs,      /* 卫星位置、速度，每颗卫星占 6 个 double */
    const double *dts,     /* 卫星钟差、钟漂，每颗卫星占 2 个 double，钟差单位 s */
    const double *var_rs,  /* 各卫星星历和钟差误差方差，单位 m^2 */
    const int *svh,        /* 各卫星健康状态 */
    const double *dr,      /* 地球潮汐等造成的接收机坐标改正，单位 m */
    int *exc,              /* 输入/输出：卫星排除标志 */
    const nav_t *nav,      /* 星历、钟差、天线、偏差等导航数据 */
    const double *x,       /* 本次计算残差所使用的 PPP 状态向量 */
    rtk_t *rtk,            /* PPP 工作区，用于选项和每颗卫星的状态记录 */
    double *v,             /* 输出：有效观测的残差向量 */
    double *H,             /* 输出：设计矩阵，按“状态数 × 残差数”列优先存储 */
    double *R,             /* 输出：观测误差协方差矩阵 */
    double *azel)          /* 输出：每颗卫星的方位角、高度角 */
{
    prcopt_t *opt=&rtk->opt; /* 处理选项的简写指针 */
    double y;                /* 当前参与计算的载波相位或伪距观测值，单位 m */
    double r,cdtr,bias;      /* 几何距离、接收机钟差距离、载波相位偏差，单位 m */
    double rr[3],pos[3],e[3]; /* 潮汐改正后 ECEF 坐标、大地坐标、接收机到卫星单位向量 */
    double dtdx[3];          /* 对流层延迟对 ZTD/梯度状态的偏导数 */
    double L[NFREQ],P[NFREQ],Lc,Pc; /* 改正后的相位/伪距及其无电离层组合，单位 m */
    double var[MAXOBS*2*NFREQ]; /* 每一个有效残差的总方差 */
    double dtrp=0.0,dion=0.0;   /* 对流层延迟和 L1 等效斜向电离层延迟，单位 m */
    double vart=0.0,vari=0.0;   /* 对流层、电离层改正模型的方差 */
    double dcb,freq;             /* 第 3 频率接收机码偏差和当前观测载波频率 */
    double dantr[NFREQ]={0},dants[NFREQ]={0}; /* 接收机、卫星天线改正，单位 m */
    double ve[MAXOBS*2*NFREQ]={0},vmax=0; /* 超限的后验残差列表及其中最大值 */
    char str[40];                /* 用于调试输出的历元时间字符串 */
    int ne=0;                    /* 超限后验残差的数量 */
    int obsi[MAXOBS*2*NFREQ]={0},frqi[MAXOBS*2*NFREQ]; /* 超限残差对应的观测和频率/类型下标 */
    int maxobs,maxfrq,rej;       /* 最大粗差对应的观测、频率/类型及其列表下标 */
    int i,j,k,sat,sys;           /* 循环下标、卫星编号和卫星系统 */
    int nv=0,nx=rtk->nx,stat=1; /* 当前有效残差数、完整状态数和后验检查状态 */
    int frq,code;                /* 频率下标；code=0 为相位，code=1 为伪距 */

    time2str(obs[0].time,str,2); /* 当前历元时间转成字符串，仅供日志使用 */

    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) rtk->ssat[i].vsat[j]=0; /* 清除上一轮“相位有效”标志 */

    for (i=0;i<3;i++) rr[i]=x[i]+dr[i]; /* 状态中的接收机坐标加上潮汐位移改正 */
    ecef2pos(rr,pos);                    /* ECEF 坐标转纬度、经度和高程，供高度角/大气模型使用 */

    for (i=0;i<n&&i<MAXOBS;i++) { /* 逐颗卫星构造相位和伪距残差 */
        sat=obs[i].sat;            /* RTKLIB 内部卫星编号 */

        /* line-of-sight vector from receiver to satellite */
        if ((r=geodist(rs+i*6,rr,e))<=0.0|| /* 计算几何距离 r 和视线单位向量 e */
            satazel(pos,e,azel+i*2)<opt->elmin) { /* 计算方位/高度角并检查截止高度角 */
            exc[i]=1;                            /* 几何无效或高度角太低，排除整颗卫星 */
            continue;
        }
        if (!(sys=satsys(sat,NULL))||!rtk->ssat[sat-1].vs||
            satexclude(sat,var_rs[i],svh[i],opt)||exc[i]) { /* 检查系统、卫星位置、健康和用户排除设置 */
            exc[i]=1;                                  /* 任一检查失败就排除该卫星 */
            continue;
        }
        /* tropospheric and ionospheric model */
        if (!model_trop(obs[i].time,pos,azel+i*2,opt,x,dtdx,nav,&dtrp,&vart)|| /* 算对流层延迟、方差和状态偏导 */
            !model_iono(obs[i].time,pos,azel+i*2,opt,sat,x,nav,&dion,&vari)) { /* 算电离层延迟及方差 */
            continue;
        }
        /* satellite and receiver antenna model */
        if (opt->posopt[0]) satantpcv(rs+i*6,rr,nav->pcvs+sat-1,dants); /* 可选：卫星天线相位中心变化改正 */
        antmodel(opt->pcvr,opt->antdel[0],azel+i*2,opt->posopt[1],dantr); /* 接收机天线相位中心改正 */

        /* phase windup model */
        if (!model_phw(rtk->sol.time,sat,nav->pcvs[sat-1].type,
                       opt->posopt[2]?2:0,rs+i*6,rr,&rtk->ssat[sat-1].phw)) {
            continue;
        }
        /* corrected phase and code measurements */
        corr_meas(obs+i,nav,azel+i*2,&rtk->opt,dantr,dants,
                  rtk->ssat[sat-1].phw,L,P,&Lc,&Pc); /* 把原始观测改正并统一为 m */

        /* stack phase and code residuals {L1,P1,L2,P2,...} */
        for (j=0;j<2*NF(opt);j++) { /* 顺序为 L1、P1、L2、P2…… */
            double C=0.0;           /* 当前频率电离层系数：相位为负，伪距为正 */

            dcb=bias=0.0; /* 每条观测先清空可选 DCB 和相位偏差项 */
            code=j%2;     /* 偶数 j 是载波相位，奇数 j 是伪距 */
            frq=j/2;      /* j=0/1 对应频率 0，j=2/3 对应频率 1 */

            if (opt->ionoopt==IONOOPT_IFLC) { /* 无电离层组合模式只使用 Lc/Pc */
                if ((y=code==0?Lc:Pc)==0.0) continue; /* 缺少组合观测就跳过 */
            }
            else {
                if ((y=code==0?L[frq]:P[frq])==0.0) continue; /* 非组合模式取当前频率 L/P */

                if ((freq=sat2freq(sat,obs[i].code[frq],nav))==0.0) continue; /* 观测码无法映射到频率则跳过 */
                /* The iono paths have already applied a slant factor. */
                C=SQR(FREQL1/freq)*(code==0?-1.0:1.0); /* I_f=(f_L1/f)^2 I_L1，相位含 -I、伪距含 +I */
            }
            if (H) {
                for (k=0;k<nx;k++) H[k+nx*nv]=0.0; /* 先清空第 nv 条残差对应的一整列偏导 */
                for (k=0;k<3;k++) H[k+nx*nv]=-e[k]; /* 距离对接收机 XYZ 的偏导为视线向量负值 */
            }

            /* receiver clock */
            switch (sys) { /* GPS 使用基准钟；其他系统使用各自的系统间钟差状态 */
                case SYS_GLO: k=1; break;
                case SYS_GAL: k=2; break;
                case SYS_CMP: k=3; break;
                case SYS_IRN: k=4; break;
                default:      k=0; break;
            }
            cdtr=x[IC(k,opt)]; /* 从状态向量取对应系统的接收机钟差，单位已经是 m */
            if (H) {
                H[IC(k,opt)+nx*nv]=1.0; /* 预测距离对接收机钟差的偏导为 +1 */

                if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                    for (k=0;k<(opt->tropopt>=TROPOPT_ESTG?3:1);k++) {
                        H[IT(opt)+k+nx*nv]=dtdx[k]; /* 填入对 ZTD（以及南北/东西梯度）的偏导 */
                    }
                }
            }
            if (opt->ionoopt==IONOOPT_EST) {
                if (rtk->x[II(sat,opt)]==0.0) continue; /* 该卫星电离层状态尚未初始化，不能使用 */
                /* The vertical iono delay is estimated, but the residual is
                 * in the direction of the slant, so apply the slant factor
                 * mapping function. */
                if (H) H[II(sat,opt)+nx*nv]=C*ionmapf(pos,azel+i*2); /* 垂直状态映射到斜路径后再乘频率/类型系数 */
            }
            if (frq==2&&code==1) { /* L5-receiver-dcb */
                dcb+=rtk->x[ID(opt)];              /* 第 3 频率伪距加入接收机 DCB 状态 */
                if (H) H[ID(opt)+nx*nv]=1.0;       /* 对该 DCB 状态的偏导为 +1 */
            }
            if (code==0) { /* phase bias */
                if ((bias=x[IB(sat,frq,opt)])==0.0) continue; /* 相位必须已有本星本频率模糊度状态 */
                if (H) H[IB(sat,frq,opt)+nx*nv]=1.0;          /* 相位预测值对模糊度的偏导为 +1 */
            }
            /* residual */
            double res=y-(r+cdtr-CLIGHT*dts[i*2]+dtrp+C*dion+dcb+bias); /* v=实测－(几何+接收机钟－卫星钟+大气+偏差) */
            if (v) v[nv]=res; /* 把这条有效残差按紧凑顺序放入滤波器输入向量 */

            if (code==0) rtk->ssat[sat-1].resc[frq]=res;  /* carrier phase */
            else         rtk->ssat[sat-1].resp[frq]=res;  /* pseudorange */

            /* variance */
            var[nv]=varerr(sat,sys,azel[1+i*2], /* 观测自身噪声：系统、相位/码、频率、高度角和 SNR 等共同决定 */
                           rtk->ssat[sat-1].snr_rover[frq],
                           j,opt,obs+i);
            var[nv] +=vart+SQR(C)*vari+var_rs[i]; /* 再加对流层、电离层和卫星轨道/钟差方差 */
            if (sys==SYS_GLO&&code==1) var[nv]+=VAR_GLO_IFB; /* GLONASS 伪距额外考虑频间偏差不确定度 */

            trace(3,"%s post=%2d sat=%2d %s%d res=%9.4f sig=%9.4f el=%4.1f\n",
                  str,post,sat,code?"P":"L",frq+1,res,sqrt(var[nv]),azel[1+i*2]*R2D);

            /* reject satellite by pre-fit residuals */
            double maxinno = (post==-1?1000:1)*opt->maxinno[code]; /* 初始粗定位阶段放宽 1000 倍，避免过早删星 */
            if (post<=0&&opt->maxinno[code]>0.0&&fabs(res)>maxinno) {
                trace(2,"outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n",
                      post,str,sat,code?"P":"L",frq+1,res,azel[1+i*2]*R2D);
                exc[i]=1; rtk->ssat[sat-1].rejc[frq]++; /* 先验残差过大：立即排除该卫星并累计拒绝次数 */
                continue;
            }
            /* record large post-fit residuals */
            if (post>0&&fabs(res)>sqrt(var[nv])*THRES_REJECT) {
                obsi[ne]=i; frqi[ne]=j; ve[ne]=res; ne++; /* 记录超限后验残差，稍后只剔除最大的一条 */
            }
            if (code==0) rtk->ssat[sat-1].vsat[frq]=1; /* 相位残差成功入列，标记该星该频率有效 */
            nv++; /* 有效残差数加一；下一条观测使用 v/H/var 的下一个位置 */
        }
    }
    /* reject satellite with large and max post-fit residual */
    if (post>0&&ne>0) {
        vmax=ve[0]; maxobs=obsi[0]; maxfrq=frqi[0]; rej=0; /* 先假设第一条是最大粗差 */
        for (j=1;j<ne;j++) {
            if (fabs(vmax)>=fabs(ve[j])) continue;
            vmax=ve[j]; maxobs=obsi[j]; maxfrq=frqi[j]; rej=j; /* 更新绝对值最大的后验残差 */
        }
        sat=obs[maxobs].sat;
        trace(2,"outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n",
              post,str,sat,maxfrq%2?"P":"L",maxfrq/2+1,vmax,azel[1+maxobs*2]*R2D);
        exc[maxobs]=1; rtk->ssat[sat-1].rejc[maxfrq%2]++; stat=0; /* 排除其整颗卫星，并要求外层重新滤波 */
        ve[rej]=0; /* 当前局部列表中清掉已选中的粗差值 */
    }
    if (R) {
        for (j=0;j<nv;j++) for (i=0;i<nv;i++) R[i+j*nv]=0.0; /* 本实现假定各观测误差互不相关，先清零非对角项 */
        for (i=0;i<nv;i++) R[i+i*nv]=var[i];                  /* 把每条残差方差填到 R 的主对角线 */
    }
    return post>0?stat:nv; /* 后验阶段返回是否通过；先验阶段返回可交给 filter() 的残差数 */
}
/* number of estimated states ------------------------------------------------*/
extern int pppnx(const prcopt_t *opt)
{
    return NX(opt);
}
/* update solution status ----------------------------------------------------*/
static void update_stat(rtk_t *rtk, const obsd_t *obs, int n, int stat)
{
    const prcopt_t *opt=&rtk->opt;
    int i,j;

    /* test # of valid satellites */
    rtk->sol.ns=0;
    for (i=0;i<n&&i<MAXOBS;i++) {
        for (j=0;j<opt->nf;j++) {
            if (!rtk->ssat[obs[i].sat-1].vsat[j]) continue;
            rtk->ssat[obs[i].sat-1].lock[j]++;
            rtk->ssat[obs[i].sat-1].outc[j]=0;
            if (j==0) rtk->sol.ns++;
        }
    }
    rtk->sol.stat=rtk->sol.ns<MIN_NSAT_SOL?SOLQ_NONE:stat;

    if (rtk->sol.stat==SOLQ_FIX) {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->xa[i];
            rtk->sol.qr[i]=(float)rtk->Pa[i+i*rtk->na];
        }
        rtk->sol.qr[3]=(float)rtk->Pa[1];
        rtk->sol.qr[4]=(float)rtk->Pa[1+2*rtk->na];
        rtk->sol.qr[5]=(float)rtk->Pa[2];
    }
    else {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->x[i];
            rtk->sol.qr[i]=(float)rtk->P[i+i*rtk->nx];
        }
        rtk->sol.qr[3]=(float)rtk->P[1];
        rtk->sol.qr[4]=(float)rtk->P[2+rtk->nx];
        rtk->sol.qr[5]=(float)rtk->P[2];

        if (rtk->opt.dynamics) { /* velocity and covariance */
            for (i=3;i<6;i++) {
                rtk->sol.rr[i]=rtk->x[i];
                rtk->sol.qv[i-3]=(float)rtk->P[i+i*rtk->nx];
            }
            rtk->sol.qv[3]=(float)rtk->P[4+3*rtk->nx];
            rtk->sol.qv[4]=(float)rtk->P[5+4*rtk->nx];
            rtk->sol.qv[5]=(float)rtk->P[5+3*rtk->nx];
        }
    }
    rtk->sol.dtr[0]=rtk->x[IC(0,opt)]/CLIGHT; /* GPS */
    rtk->sol.dtr[1]=(rtk->x[IC(1,opt)]-rtk->x[IC(0,opt)])/CLIGHT; /* GLO-GPS */
    rtk->sol.dtr[2]=(rtk->x[IC(2,opt)]-rtk->x[IC(0,opt)])/CLIGHT; /* GAL-GPS */
    rtk->sol.dtr[3]=(rtk->x[IC(3,opt)]-rtk->x[IC(0,opt)])/CLIGHT; /* BDS-GPS */

    for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<opt->nf;j++) {
        rtk->ssat[obs[i].sat-1].snr_rover[j]=obs[i].SNR[j];
        rtk->ssat[obs[i].sat-1].snr_base[j] =0;
    }
    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) {
        if (rtk->ssat[i].slip[j]&(LLI_SLIP|LLI_HALFC)) rtk->ssat[i].slipc[j]++;
        if (rtk->ssat[i].fix[j]==2&&stat!=SOLQ_FIX) rtk->ssat[i].fix[j]=1;
    }
}
/* test hold ambiguity -------------------------------------------------------*/
static int test_hold_amb(rtk_t *rtk)
{
    int i,j,stat=0;

    /* no fix-and-hold mode */
    if (rtk->opt.modear!=ARMODE_FIXHOLD) return 0;

    /* reset # of continuous fixed if new ambiguity introduced */
    for (i=0;i<MAXSAT;i++) {
        if (rtk->ssat[i].fix[0]!=2&&rtk->ssat[i].fix[1]!=2) continue;
        for (j=0;j<MAXSAT;j++) {
            if (rtk->ssat[j].fix[0]!=2&&rtk->ssat[j].fix[1]!=2) continue;
            if (!rtk->ambc[j].flags[i]||!rtk->ambc[i].flags[j]) stat=1;
            rtk->ambc[j].flags[i]=rtk->ambc[i].flags[j]=1;
        }
    }
    if (stat) {
        rtk->nfix=0;
        return 0;
    }
    /* test # of continuous fixed */
    return ++rtk->nfix>=rtk->opt.minfix;
}
/* precise point positioning -------------------------------------------------*/
extern void pppos(       /* 一句话理解 pppos()：接收观测与导航数据，进行一个历元的 PPP 定位 */
    rtk_t *rtk,          /* PPP 工作区：保存配置、待估状态、协方差和定位结果 */
    const obsd_t *obs,   /* 当前历元的卫星观测值：伪距、载波相位、信噪比等 */
    int n,               /* obs 数组中观测记录的数量 */
    const nav_t *nav)    /* 导航数据：星历、卫星钟差及其他改正数据 */
{
    const prcopt_t *opt=&rtk->opt; /* 用 opt 简写 rtk->opt；这里只读取定位配置 */
    double *rs;           /* 每颗卫星的位置和速度，后面由 satposs() 算出 */
    double *dts;          /* 每颗卫星的钟差和钟漂，后面由 satposs() 算出 */
    double *var;          /* 每颗卫星位置和钟差计算结果的方差（不确定程度） */
    double *v;            /* 观测残差：实际观测值与当前模型计算值之差 */
    double *H;            /* 设计矩阵：描述各个状态量会怎样影响观测值 */
    double *R;            /* 观测噪声协方差矩阵：描述各项观测值的可信程度 */
    double *azel;         /* 每颗卫星的方位角和高度角：表示卫星在天空中的方向 */
    double *xp;           /* 状态向量 rtk->x 的临时副本：滤波器先在这里试算更新 */
    double *Pp;           /* 协方差矩阵 rtk->P 的临时副本：表示试算状态的不确定程度 */
    double dr[3]={0};     /* 地球潮汐引起的接收机三维位移改正量，初始都为 0 */
    double std[3];        /* 最终三维位置状态的标准差，用来检查固定解是否足够可靠 */
    char str[40];         /* 把当前观测时间转成文字后暂存在这里，主要用于日志输出 */
    int i,j;              /* 循环计数器 */
    int nv;               /* 当前建立出的有效残差（观测方程）数量 */
    int info;             /* filter() 的返回码：0 表示正常，非 0 表示滤波失败 */
    int svh[MAXOBS];      /* 每条观测对应卫星的健康状态标记 */
    int exc[MAXOBS]={0};  /* 卫星排除标记：0=暂不排除，1=排除；初始全部为 0 */
    int stat=SOLQ_SINGLE; /* 当前解状态先设为单点解，PPP 成功后再升级为 PPP 解 */

    time2str(obs[0].time,str,2); /* 把当前历元时间转成保留 2 位小数的文字，写入 str */
    trace(3,"pppos   : time=%s nx=%d n=%d\n",str,rtk->nx,n); /* 输出时间、状态数和观测数到调试日志 */

    rs=mat(6,n);          /* 分配 6×n 空间：每颗卫星存 3 个位置和 3 个速度分量 */
    dts=mat(2,n);         /* 分配 2×n 空间：每颗卫星存钟差和钟漂 */
    var=mat(1,n);         /* 分配 1×n 空间：每颗卫星存一个位置/钟差方差 */
    azel=zeros(2,n);      /* 分配并清零 2×n 空间：每颗卫星存方位角和高度角 */

    for (i=0;i<MAXSAT;i++)       /* i 遍历程序支持的所有卫星编号 */
        for (j=0;j<opt->nf;j++)  /* j 遍历当前配置使用的所有频率 */
            rtk->ssat[i].fix[j]=0; /* 清除该卫星、该频率上一次留下的模糊度固定标记，不是删除状态向量中的模糊度估计值 */
    for (i=0;i<n&&i<MAXOBS;i++)       /* i 遍历当前历元的观测记录，但不超过 MAXOBS */
        for (j=0;j<opt->nf;j++) {     /* j 遍历当前配置使用的所有频率 */
            rtk->ssat[obs[i].sat-1].snr_rover[j]=obs[i].SNR[j]; /* 保存流动站该频率的信噪比 */
            rtk->ssat[obs[i].sat-1].snr_base[j] =0;             /* PPP 没有基准站观测，故设为 0 */
        }

    /* temporal update of ekf states */
    udstate_ppp(rtk,obs,n,nav); /* 将上一历元的 PPP 状态预测到当前历元，并初始化新出现的状态 */

    /* satellite positions and clocks  计算卫星位置和钟差 */
    satposs(
        obs[0].time,      /* 输入：当前历元的观测时间 */
        obs,              /* 输入：当前历元的卫星观测记录 */
        n,                /* 输入：观测记录数量 */
        nav,              /* 输入：计算卫星轨道和钟差所需的导航数据 */
        rtk->opt.sateph,  /* 输入：选择使用广播星历、精密星历还是 SSR 星历 */
        rs,               /* 输出：各卫星的位置和速度 */
        dts,              /* 输出：各卫星的钟差和钟漂 */
        var,              /* 输出：各卫星位置和钟差结果的方差 */
        svh);             /* 输出：各卫星的健康状态 */

    /* exclude measurements of eclipsing satellite (block IIA)  卫星地影检查 */
    if (rtk->opt.posopt[3]) {       /* 如果配置中启用了卫星地影检查 */
        testeclipse(obs,n,nav,rs);  /* 检查处于地影期的 Block IIA 卫星，并使其不参与后续定位 */
    }
    /* earth tides correction  地球潮汐改正 */
    if (opt->tidecorr) {            /* 如果配置中启用了地球潮汐改正 */
        tidedisp(gpst2utc(obs[0].time),rtk->x,opt->tidecorr,&nav->erp,
                 opt->odisp[0],dr); /* 计算潮汐造成的测站三维位移，并把结果写入 dr */
    }
    nv=n*rtk->opt.nf*2+MAXSAT+3; /* 暂存观测方程容量上限：每频率含相位和伪距，并额外预留空间=n×频率数×2*/
    xp=mat(rtk->nx,1);           /* 分配临时状态向量：共 rtk->nx 个状态 */
    Pp=zeros(rtk->nx,rtk->nx);   /* 分配并清零临时状态协方差矩阵 */
    v=mat(nv,1);                 /* 按容量上限分配残差向量 */
    H=mat(rtk->nx,nv);           /* 按“状态数×容量上限”分配设计矩阵 */
    R=mat(nv,nv);                /* 按“容量上限×容量上限”分配观测噪声协方差矩阵 */

    for (i=0;i<MAX_ITER;i++) { /* 最多尝试 MAX_ITER（本版本为 8）次滤波和残差检查 */

        matcpy(xp,rtk->x,rtk->nx,1);           /* 把正式预测状态 rtk->x 复制到临时状态 xp */
        matcpy(Pp,rtk->P,rtk->nx,rtk->nx);     /* 把正式协方差 rtk->P 复制到临时协方差 Pp */

      /* prefit residuals    验前残差 → filter更新 → 验后残差检查
                                                      ↓
                                               合格：正式保存结果
                                              不合格：剔除异常观测后重算
         * NOTE: use different limit for pre-fit residuals in first iteration
         *       by using argument post = -1
         */
        if (!(nv=ppp_res(i==0?-1:0,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,v,H,R,azel))) { /* 建立验前残差及 H、R；nv=0 表示没有有效观测 */
            trace(2,"%s ppp (%d) no valid obs data\n",str,i+1); /* 在日志中记录时间和第几次尝试没有有效观测 */
            break; /* 退出最多 8 次的滤波迭代循环，不再继续本历元的 PPP 更新，但仍可处理下一个历元。*/
        }
        /* measurement update of ekf states */
        if ((info=filter(xp,Pp,H,v,R,rtk->nx,nv))) { /* 用 H、v、R 更新 xp、Pp；返回 0 成功，非 0 失败 */
            trace(2,"%s ppp (%d) filter error info=%d\n",str,i+1,info); /* 记录滤波失败的尝试次数和错误码 */
            break; /* 滤波失败，退出本历元的迭代循环 */
        }
        /* postfit residuals */
        if (ppp_res(i+1,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,NULL,NULL,NULL,azel)) { /* 用滤波后的 xp 计算验后残差；通过检查时返回真 */
            matcpy(rtk->x,xp,rtk->nx,1);       /* 验后检查通过：把临时状态 xp 写回正式状态 rtk->x */
            matcpy(rtk->P,Pp,rtk->nx,rtk->nx); /* 把临时协方差 Pp 写回正式协方差 rtk->P */
            stat=SOLQ_PPP;                     /* 将本历元解状态标记为 PPP 浮点解 */
            break;                             /* 已得到可接受结果，结束迭代 */
        }
    }
    if (i>=MAX_ITER) { /* 循环自然结束且 i 达到上限，表示 8 次尝试都没有得到可接受结果 */
        trace(2,"%s ppp (%d) iteration overflows\n",str,i); /* 在日志中记录本历元 PPP 迭代超限 */
    }
    if (stat==SOLQ_PPP) { /* 只有先得到可接受的 PPP 浮点解，才尝试固定载波相位模糊度 */

        if (ppp_ar(rtk,obs,n,exc,nav,azel,xp,Pp)&& /* 尝试生成固定解；但本仓库该函数目前固定返回 0 */
            ppp_res(9,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,NULL,NULL,NULL,azel)) { /* 若固定成功，再检查固定解的验后残差 */

            matcpy(rtk->xa,xp,rtk->nx,1);       /* 保存通过残差检查的固定候选状态 */
            matcpy(rtk->Pa,Pp,rtk->nx,rtk->nx); /* 保存固定候选状态的协方差矩阵 */

            for (i=0;i<3;i++)                   /* 依次处理接收机 X、Y、Z 三个位置状态 */
                std[i]=sqrt(Pp[i+i*rtk->nx]);   /* 协方差对角元素开平方，得到该方向的位置标准差 */
            if (norm(std,3)<MAX_STD_FIX)         /* 三维综合标准差小于 0.15 m，才认可固定解 */
                stat=SOLQ_FIX;                   /* 将本历元解状态正式标记为固定解 */
        }
        else {                    /* 模糊度固定失败，或固定候选解未通过残差检查 */
            rtk->nfix=0;          /* 连续固定成功次数清零 */
        }
        /* update solution status */
        update_stat(rtk,obs,n,stat); /* 将本历元的解类型、坐标和卫星状态整理到输出结果中 */

        if (stat==SOLQ_FIX&&test_hold_amb(rtk)) { /* 当前为固定解，且 fix-and-hold 的连续固定条件已满足 */
            matcpy(rtk->x,xp,rtk->nx,1);         /* 用固定状态覆盖浮点状态，作为后续历元的起点 */
            matcpy(rtk->P,Pp,rtk->nx,rtk->nx);   /* 同时用固定解协方差覆盖浮点协方差 */
            trace(2,"%s hold ambiguity\n",str);  /* 记录本历元开始保持模糊度 */
            rtk->nfix=0;                         /* 本次保持完成后，重新累计连续固定次数 */
        }
    }
    free(rs);   /* 释放卫星位置和速度数组 */
    free(dts);  /* 释放卫星钟差和钟漂数组 */
    free(var);  /* 释放卫星位置和钟差方差数组 */
    free(azel); /* 释放卫星方位角和高度角数组 */
    free(xp);   /* 释放临时状态向量 */
    free(Pp);   /* 释放临时状态协方差矩阵 */
    free(v);    /* 释放残差向量 */
    free(H);    /* 释放设计矩阵 */
    free(R);    /* 释放观测噪声协方差矩阵 */
}
