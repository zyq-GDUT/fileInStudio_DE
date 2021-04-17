/*该遗传算法中,适应度只取正数,同时,
每个个体的适应度都和目标函数的函数值相同
测试函数的变量取值范围从文件“D://bound.txt”中读取
算法的输出结果输出到文件“D://delog.txt"

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<string.h>

#define POPSIZE 10   /*个体的数量*/
#define ARRAYSIZE POPSIZE+2 /*种群数组的大小,多出来的两个位置,用来保存最好个体和最差个体*/
#define MAXGENS 3000 /*最大的迭代次数*/
#define NVARS 30      /*函数的维数,即函数未知量个数*/
#define dim   NVARS    /*测试函数中的函数维数*/
#define CR  0.9     /*变量发生交叉的概率*/
#define FUNCTIONNUM 23
#define TRAILNUM     10
#define TRUE 1
#define FALSE 0
#define F 0.5     /*缩放因子F一般取[0,2]的常量*/
#define LAMBDA F  /*λ值与缩放因子取值一致*/
#define BEST_INDEX POPSIZE
#define WORST_INDEX POPSIZE+1
#define DELOG "DEreport\\DE_rand_1_bin.txt"
#define STRATEGY_NUM 10 /*10个差分策略*/



int generation;/*当前的代数*/
FILE *delog ;    /*文件指针*/
double best_val[TRAILNUM+1];//十次测试的最优值
double worst_val[TRAILNUM+1];//十次测试的最差值
double mean;//用来保存十次测试最优值的平均值
double stddev;//用来保存十次测试最优值的标准差
double accu_best_gen;//用来累加十次中最优个体的迭代次数

struct genotype /*基因型,代表一个个体*/
{
    double gene[NVARS]; //用来保存每个基因，即每个变量
    double fitness; //当前个体的适应度
    double upper[NVARS];    //变量的取值范围中的上界，不同下标的值代表对应的变量的上界
    double lower[NVARS];    //变量的取值范围中的下界，不同下标的值代表对应的变量的下界
    int generation;  //得到最优个体时的迭达次数
};
struct genotype population[ARRAYSIZE]; //用来保存所有个体，最后两个位置用来保存最好和最差个体
struct genotype newpopulation[ARRAYSIZE];//用来保存所有新生成的个体
struct genotype mutatepopulation[ARRAYSIZE];
struct genotype crosspopulation[ARRAYSIZE];

/*用来保存测试函数的变量取值的上界和下界*/
struct boundtype
{
    double lbound;
    double ubound;
    int    nvars;//当前测试函数的维数
};
struct boundtype  bound[FUNCTIONNUM+1];

/*不同差分策略的测试报告打印在不同文件*/
struct printfpath
{
    char* strategyname;
    char* path;
};
struct printfpath path[STRATEGY_NUM+1];
/************************************测试函数**********************************/
double f1(double x[])
{
	double rt1 = 0;
	int i;
	for(i=0;i<dim;i++){
		rt1 += (x[i]*x[i]);
	}
	return rt1;
}

double f2(double x[])
{
	double rt1, rt2;
	int i;
	rt1 = 0 , rt2 = 1;
	for(i=0;i<dim;i++){
		if(x[i]>0){
			rt1 += x[i];
			rt2 *= x[i];
		}else{
			rt1 -= x[i];
			rt2 *= -x[i];
		}
	}
	rt1 += rt2;
	return rt1;
}

double f3(double x[])
{
	double rt1, rt2;
	int i,j;
	rt1 = 0;
	for(i=0;i<dim;i++){
		rt2 = 0;
		for(j=0;j<=i;j++){
			rt2 += x[j];
		}
		rt1 += (rt2*rt2);
	}
	return rt1;
}

double f4(double x[])
{
	double rt1;
	int i;

	rt1 = 0;
	for(i=0;i<dim;i++){
		if(x[i]>0 && x[i] > rt1){
			rt1 = x[i];
		}
		if(x[i]<0&&-x[i] > rt1){
			rt1 = -x[i];
		}
	}
	return rt1;
}

double f5(double x[])
{	double rt1, rt2;
	int i;
	rt1 = 0,rt2 = 0;
	for(i=0;i<dim-1;i++){
		rt1 += 100*((x[i+1]-x[i]*x[i])*(x[i+1]-x[i]*x[i])) + (x[i]-1)*(x[i]-1);
	}
	return rt1;
}

double f6(double x[])
{
	double rt1;
	int i;
	rt1 = 0;
	for(i=0;i<dim;i++){
		rt1 += (x[i]+0.5)*(x[i]+0.5);
	}
	return rt1;
}

double f7(double x[])
{	double rt1;
	int i;
	rt1 = 0;
	for(i=0;i<dim;i++){
		rt1 += (i+1)*x[i]*x[i]*x[i]*x[i];
	}
	rt1 += rand()*1.0/RAND_MAX;
	return rt1;
}

double f8(double x[])
{	double rt1, rt2;
	int i;
	rt1 = 0;
	for(i=0;i<dim;i++){
		if(x[i]<0)rt2 = -x[i];
		else rt2 = x[i];
		rt1 += -x[i]*sin(sqrt(rt2));
	}
	return rt1;
}

double f9(double x[])
{
	double rt1;
	int i;
	double pi;
	rt1 = 0;
	pi = acos(-1.0);
	for(i=0;i<dim;i++){
		rt1 += x[i]*x[i]-10*cos(2*pi*x[i])+10;
	}
	return rt1;
}

double f10(double x[])
{
	double rt1, rt2;
	int i;
	double pi;
	rt1 = 0;
	rt2 = 0;
	pi = acos(-1.0);
	for(i=0;i<dim;i++){
		rt1 += x[i]*x[i];
		rt2 +=cos(2*pi*x[i]);
	}
	rt1 = -20*exp(-0.2*sqrt(rt1/dim))-exp(rt2/dim)+20+exp(1);
	return rt1;
}

double f11(double x[])
{
	double rt1, rt2;
	int i;
	double pi;
	rt1 = 0;
	rt2 = 1;
	pi = acos(-1.0);
	for(i=0;i<dim;i++){
		rt1 += x[i]*x[i];
		rt2 *=cos(x[i]/sqrt((double)i+1));
	}
	rt1 = 1.0/4000*rt1-rt2+1;
	return rt1;
}

double f12(double x[])
{	double rt1, rt2;
	int i;
	double pi,u,ya,yb,y1;
	rt1 = 0;
	rt2 = 0;
	pi = acos(-1.0);
	for(i=0;i<dim;i++){
		if(x[i]>10)
			u = 100*(x[i]-10)*(x[i]-10)*(x[i]-10)*(x[i]-10);
		if(x[i]>=-10&&x[i]<=10)
			u = 0;
		if(x[i]<-10)
			u = 100*(-x[i]-10)*(-x[i]-10)*(-x[i]-10)*(-x[i]-10);
		rt2 += u;
		ya = yb;
		yb = 1+ 1.0/4*(x[i]+1);
		if(i!=0){
			double t = sin(pi*yb);
			rt1 += (ya-1)*(ya-1)*(1+10*t*t);
		}else{
			y1 = yb;
		}
	}
	rt1 = pi/dim*(10*sin(pi*y1)*sin(pi*y1) + rt1 + (yb-1)*(yb-1)) + rt2;
	return rt1;
}

double f13(double x[])
{
	double rt1, rt2;
	int i;
	double pi,u,t;
	rt1 = 0;
	rt2 = 0;
	pi = acos(-1.0);
	for(i=0;i<dim;i++){
		if(x[i]>5)
			u = 100*(x[i]-5)*(x[i]-5)*(x[i]-5)*(x[i]-5);
		if(x[i]>=-5&&x[i]<=5)
			u = 0;
		if(x[i]<-5)
			u = 100*(-x[i]-5)*(-x[i]-5)*(-x[i]-5)*(-x[i]-5);
		rt2 += u;
		if(i!=0){
			t = sin(3*pi*x[i]);
			rt1 += (x[i-1]-1)*(x[i-1]-1)*(1+t*t);
		}
	}
	rt1 = 0.1*(sin(3*pi*x[0])*sin(3*pi*x[0]) + rt1 + (x[dim-1]-1)*(x[dim-1]-1)*(1+ sin(2*pi* x[dim-1])*sin(2*pi* x[dim-1]))) +rt2;
	return rt1;
}

double f14(double x[])
{
	double rt1;
	int j;
	rt1 = 0;
	double a[2][25]={{-32,-16,0,16,32,-32,-16,0,16,32,-32,-16,0,16,32,-32,-16,0,16,32,-32,-16,0,16,32},
			{-32,-32,-32,-32,-32,-16,-16,-16,-16,-16,0,0,0,0,0,16,16,16,16,16,32,32,32,32,32}};
	for(j=0;j<25;j++){
		rt1 += 1.0/((j+1) + pow(x[0] - a[0][j],6) + pow(x[1] - a[1][j],6));
	}
	rt1 = 1.0/(1.0/500 + rt1);
	return rt1;
}

double f15(double x[])
{
	double rt1;
	int i;
	rt1 = 0;
	double a[] = {0.1957,0.1947,0.1735,0.1600,0.0844,0.0627,0.0456,0.0342,0.0323,0.0235,0.0246};
	double b[] = {0.25,0.5,1,2,4,6,8,10,12,14,16};
	for(i=0;i<11;i++){
		double b1 = 1.0/b[i];
		rt1 += (a[i] - x[0]*(b1*b1+b1*x[1])/(b1*b1 + b1*x[2] +x[3]))*
				(a[i] - x[0]*(b1*b1+b1*x[1])/(b1*b1 + b1*x[2] +x[3]));
	}
	return rt1;
}

double f16(double x[])
{
	double rt1;
	rt1 = 0;
	rt1 = 4*x[0]*x[0] -2.1*pow(x[0],4) +1.0/3*pow(x[0],6) +x[0]*x[1]-4*x[1]*x[1] +4*pow(x[1],4);
	return rt1;
}

double f17(double x[])
{
	double pi,rt1;
	rt1 = 0;
	 pi = acos(-1.0);
	rt1 = ((x[1]+5) - 5.1/(4*pi*pi)*x[0]*x[0] + 5.0/pi*x[0]-6)*
			((x[1]+5) - 5.1/(4*pi*pi)*x[0]*x[0] + 5.0/pi*x[0]-6) +10 * (1-1.0/(8*pi))*cos(x[0]) +10;
	return rt1;
}

double f18(double x[])
{
	double rt1;
	rt1 = 0;
	rt1 =  (1+(x[0]+x[1]+1)*(x[0]+x[1]+1)*(19-14*x[0]+3*x[0]*x[0]-14*x[1]
		+6*x[0]*x[1]+3*x[1]*x[1]))*(30+(2*x[0]-3*x[1])*(2*x[0]-3*x[1])*(18-32*x[0]
		+12*x[0]*x[0]+48*x[1]-36*x[0]*x[1]+27*x[1]*x[1]));
	return rt1;
}

double f21(double x[])
{
	double rt1;
	int i;
	double a[10][4]={{4,4,4,4},{1,1,1,1},{8,8,8,8},{6,6,6,6},{3,7,3,7},
					{2,9,2,9},{5,5,3,3},{8,1,8,1},{6,2,6,2},{7,3.6,7,3.6}};
	double c[10] = {0.1,0.2,0.2,0.4,0.4,0.6,0.3,0.7,0.5,0.5};
	rt1 = 0;
	for(i=0;i<5;i++){
		rt1 += 1.0/((x[0]-a[i][0])*(x[0]-a[i][0]) + (x[1]-a[i][1])*(x[1]-a[i][1])+
				(x[2]-a[i][2])*(x[2]-a[i][2]) + (x[3]-a[i][3])*(x[3]-a[i][3]) +c[i]);
	}
	rt1 = -rt1;
	return rt1;
}

double f22(double x[])
{
	double rt1;
	int i;

	double a[10][4]={{4,4,4,4},{1,1,1,1},{8,8,8,8},{6,6,6,6},{3,7,3,7},
					{2,9,2,9},{5,5,3,3},{8,1,8,1},{6,2,6,2},{7,3.6,7,3.6}};
	double c[10] = {0.1,0.2,0.2,0.4,0.4,0.6,0.3,0.7,0.5,0.5};
	rt1 = 0;
	for(i=0;i<7;i++){
		rt1 += 1.0/((x[0]-a[i][0])*(x[0]-a[i][0]) + (x[1]-a[i][1])*(x[1]-a[i][1])+
				(x[2]-a[i][2])*(x[2]-a[i][2]) + (x[3]-a[i][3])*(x[3]-a[i][3]) +c[i]);
	}
	rt1 = -rt1;
	return rt1;
}

double f23(double x[])
{
	double rt1;
	int i;
	double a[10][4]={{4,4,4,4},{1,1,1,1},{8,8,8,8},{6,6,6,6},{3,7,3,7},
					{2,9,2,9},{5,5,3,3},{8,1,8,1},{6,2,6,2},{7,3.6,7,3.6}};
	double c[10] = {0.1,0.2,0.2,0.4,0.4,0.6,0.3,0.7,0.5,0.5};
	rt1 = 0;
	for(i=0;i<10;i++){
		rt1 += 1.0/((x[0]-a[i][0])*(x[0]-a[i][0]) + (x[1]-a[i][1])*(x[1]-a[i][1])+
				(x[2]-a[i][2])*(x[2]-a[i][2]) + (x[3]-a[i][3])*(x[3]-a[i][3]) +c[i]);
	}
	rt1 = -rt1;
	return rt1;
}

/*选择测试函数*/
double CalFitness(double x[],int fun)
/******************************************
***输入：一组有效解************************
***输出：该解的适应值**********************
***备注：修里改适应值函数在这进行**********
******************************************/
{

	switch(fun){
	case 1:	//function 1
		return f1(x);
		break;
	case 2: //function 2
		return f2(x);
		break;
	case 3: //function 3
		return f3(x);
		break;
	case 4: //function 4
		return f4(x);
		break;
	case 5:		//function 5
		return f5(x);
		break;
	case 6:		//function 6
		return f6(x);
		break;
	case 7:		//function 7
		return f7(x);
		break;
	case 8:		//function 8
		return f8(x);
		break;
	case 9:		//function 9
		return f9(x);
		break;
	case 10:	//function 10
		return f10(x);
		break;
	case 11:	//function 11
		return f11(x);
		break;
	case 12:	//function 12
		return f12(x);
		break;
	case 13:	//function 13
		return f13(x);
		break;
	case 14:	//function 14
		return f14(x);
		break;
	case 15:	//function 15
		return f15(x);
		break;
	case 16:	//function 16
		return f16(x);
		break;
	case 17:	//function 17
		return f17(x);
		break;
	case 18:	//function 18
		return f18(x);
		break;
	case 21:	//function 21
		return f21(x);
		break;
	case 22:	//function 22
		return f22(x);
		break;
	case 23:	//function 23
		return f23(x);
		break;
	default:
		printf("Error!No such function!\n");
		return -1000;
	}
}


double randval(double,double);//生成随机值的函数
void initialize(int function);//初始化函数
void initializeBound(void);//初始化函数的定义域
void copyGene(double *a,double *b,int num);//将数组a复制到数组b
void printHead(int function,int strategy);//在文件中输出测试报告的头部信息
void printf_trailNo_bestVal_bestGen(int trailNo,int founction,int strategy);//
void save_multipletrail_best_worst(int trailNo);//保存每次测试的最好个体和最差个体的适应度
void printf_mean_stddev_worst_avggen(int strategy);//打印平均值,标准差,最差个体的适应度，获取最优值的平均代数
void  DE(int trail,int function,int strategy);//单次测试函数
void  multiple_DE(int num ,int function,int strategy);//用DE算法测试第function个函数num次
void test_DE(int num,int strategy);//每个函数测试num次
void get_mean_stdDev();//计算最优适应度的平均值和标准差
void exp_cross(int function);

void initialize_strategy_path(){

    path[1].strategyname="DE/rand/1/bin";
    path[1].path="DEreport\\DE_rand_1_bin.txt";

    path[2].strategyname="DE/rand/1/exp";
    path[2].path="DEreport\\DE_rand_1_exp.txt";

    path[3].strategyname="DE/best/1/bin";
    path[3].path="DEreport\\DE_best_1_bin.txt";

    path[4].strategyname="DE/best/1/exp";
    path[4].path="DEreport\\DE_best_1_exp.txt";

    path[5].strategyname="DE/rand/2/bin";
    path[5].path="DEreport\\DE_rand_2_bin.txt";

    path[6].strategyname="DE/rand/2/exp";
    path[6].path="DEreport\\DE_rand_2_exp.txt";

    path[7].strategyname="DE/best/2/bin";
    path[7].path="DEreport\\DE_best_2_bin.txt";

    path[8].strategyname="DE/best/2/exp";
    path[8].path="DEreport\\DE_best_2_exp.txt";

    path[9].strategyname="DE/rand-to-best/1/bin";
    path[9].path="DEreport\\DE_rand-to-best_1_bin.txt";

    path[10].strategyname="DE/rand-to-best/1/exp";
    path[10].path="DEreport\\DE_rand-to-best_1_exp.txt";

    FILE *fp;
    int i;
    //如果D:/DEreport文件夹下不存在如下文件则新建，注意D盘下要自行创建DEreport文件夹
    for(i=1;i<=STRATEGY_NUM;i++){
        if((fp=fopen(path[i].path,"r"))==NULL) //判断文件是否存在
        {
        fp=fopen(path[i].path,"w"); //如果文件不存在，新建一个文件
        }
        fclose(fp); //关闭文件
    }
}


/***五种变异策略***/

void mutate_rand_1(int function){
    int i,j;
    int rand1,rand2,rand3;
    for(i=0;i<POPSIZE;i++){
      rand1=(int)rand()%POPSIZE;
      rand2=(int)rand()%POPSIZE;
      rand3=(int)rand()%POPSIZE;
      //i与几个随机数都不能相同
      while(rand1==rand2||rand1==rand3||rand2==rand3||rand1==i||rand2==i||rand3==i)
     {
       rand1=(int)rand()%POPSIZE;
       rand2=(int)rand()%POPSIZE;
       rand3=(int)rand()%POPSIZE;
    //   printf("输出随机数：%d,%d,%d,%d\n",i,rand1,rand2,rand3);
      }
      for(j=0;j<NVARS;j++){

        mutatepopulation[i].gene[j]=population[rand1].gene[j]+F*(population[rand2].gene[j]-population[rand3].gene[j]);

         //变异后的值不能超出定义域
         if(mutatepopulation[i].gene[j]>bound[function].ubound){

            mutatepopulation[i].gene[j]=bound[function].ubound;
      }

         if(mutatepopulation[i].gene[j]<bound[function].lbound){

            mutatepopulation[i].gene[j]=bound[function].lbound;
         }
      }


    }


}

void mutate_best_1(int function){

    int i,j;
    int rand1,rand2;
    for(i=0;i<POPSIZE;i++){
      rand1=(int)rand()%POPSIZE;
      rand2=(int)rand()%POPSIZE;
      //i与几个随机数都不能相同
      while(rand1==rand2||rand1==i||rand2==i)
     {
       rand1=(int)rand()%POPSIZE;
       rand2=(int)rand()%POPSIZE;
    //   printf("输出随机数：%d,%d,%d,%d\n",i,rand1,rand2,rand3);
      }
      for(j=0;j<NVARS;j++){

        mutatepopulation[i].gene[j]=population[BEST_INDEX].gene[j]+F*(population[rand1].gene[j]-population[rand2].gene[j]);

         //变异后的值不能超出定义域
         if(mutatepopulation[i].gene[j]>bound[function].ubound){

            mutatepopulation[i].gene[j]=bound[function].ubound;
      }

         if(mutatepopulation[i].gene[j]<bound[function].lbound){

            mutatepopulation[i].gene[j]=bound[function].lbound;
         }
      }

    }

}
void mutate_rand_2(int function){
    int i,j;
    int rand1,rand2,rand3,rand4,rand5;
    for(i=0;i<POPSIZE;i++){
      rand1=(int)rand()%POPSIZE;
      rand2=(int)rand()%POPSIZE;
      rand3=(int)rand()%POPSIZE;
      rand4=(int)rand()%POPSIZE;
      rand5=(int)rand()%POPSIZE;
      //i与几个随机数都不能相同,
      while(rand1==rand2||rand1==rand3||rand2==rand3||rand1==i||rand2==i||rand3==i||rand1==rand4||rand1==rand5||rand2==rand4||rand2==rand5||rand3==rand4||rand3==rand5||rand4==rand5||i==rand4||i==rand5)
     {
       rand1=(int)rand()%POPSIZE;
       rand2=(int)rand()%POPSIZE;
       rand3=(int)rand()%POPSIZE;
       rand4=(int)rand()%POPSIZE;
       rand5=(int)rand()%POPSIZE;
    //   printf("输出随机数：%d,%d,%d,%d\n",i,rand1,rand2,rand3);
      }
      for(j=0;j<NVARS;j++){

        mutatepopulation[i].gene[j]=population[rand1].gene[j]+F*(population[rand2].gene[j]+population[rand3].gene[j]-population[rand4].gene[j]-population[rand5].gene[j]);

         //变异后的值不能超出定义域
         if(mutatepopulation[i].gene[j]>bound[function].ubound){

            mutatepopulation[i].gene[j]=bound[function].ubound;
      }

         if(mutatepopulation[i].gene[j]<bound[function].lbound){

            mutatepopulation[i].gene[j]=bound[function].lbound;
         }
      }


    }

}

void mutate_best_2(int function){
    int i,j;
    int rand1,rand2,rand3,rand4;
    for(i=0;i<POPSIZE;i++){
      rand1=(int)rand()%POPSIZE;
      rand2=(int)rand()%POPSIZE;
      rand3=(int)rand()%POPSIZE;
      rand4=(int)rand()%POPSIZE;

      //i与几个随机数都不能相同,
      while(rand1==rand2||rand1==rand3||rand2==rand3||rand1==i||rand2==i||rand3==i||rand1==rand4||rand2==rand4||rand3==rand4||i==rand4)
     {
       rand1=(int)rand()%POPSIZE;
       rand2=(int)rand()%POPSIZE;
       rand3=(int)rand()%POPSIZE;
       rand4=(int)rand()%POPSIZE;

    //   printf("输出随机数：%d,%d,%d,%d\n",i,rand1,rand2,rand3);
      }
      for(j=0;j<NVARS;j++){

        mutatepopulation[i].gene[j]=population[BEST_INDEX].gene[j]+F*(population[rand1].gene[j]+population[rand2].gene[j]-population[rand3].gene[j]-population[rand4].gene[j]);

         //变异后的值不能超出定义域
         if(mutatepopulation[i].gene[j]>bound[function].ubound){

            mutatepopulation[i].gene[j]=bound[function].ubound;
      }

         if(mutatepopulation[i].gene[j]<bound[function].lbound){

            mutatepopulation[i].gene[j]=bound[function].lbound;
         }
      }

    }

}


void mutate_randToBest_1(int function){

    int i,j;
    int rand1,rand2;
    for(i=0;i<POPSIZE;i++){
      rand1=(int)rand()%POPSIZE;
      rand2=(int)rand()%POPSIZE;

      //i与几个随机数都不能相同
      while(rand1==rand2||rand1==i||rand2==i)
     {
       rand1=(int)rand()%POPSIZE;
       rand2=(int)rand()%POPSIZE;
    //   printf("输出随机数：%d,%d,%d,%d\n",i,rand1,rand2,rand3);
      }
      for(j=0;j<NVARS;j++){

        mutatepopulation[i].gene[j]=population[i].gene[j]+LAMBDA*(population[BEST_INDEX].gene[j]-population[i].gene[j])+F*(population[rand1].gene[j]-population[rand2].gene[j]);

         //变异后的值不能超出定义域
         if(mutatepopulation[i].gene[j]>bound[function].ubound){

            mutatepopulation[i].gene[j]=bound[function].ubound;
      }

         if(mutatepopulation[i].gene[j]<bound[function].lbound){

            mutatepopulation[i].gene[j]=bound[function].lbound;
         }
      }

    }

}



/***两种交叉策略***/
/*满足二项式分布的交叉操作*/
void bin_cross(int function){
   int i,j,randj;
   double randp;
   for(i=0;i<POPSIZE;i++){
     randj=(int)(rand()%NVARS);
    for(j=0;j<NVARS;j++){
           randp=((double)(rand()%1001)/1000);
//           printf("输出概率p=%lf\n",randp);
        if(randp<=CR||randj==j){
            crosspopulation[i].gene[j]=mutatepopulation[i].gene[j];
        }else{
            crosspopulation[i].gene[j]=population[i].gene[j];
        }

    }
    //计算个体的适应度
    crosspopulation[i].fitness=CalFitness(crosspopulation[i].gene,function);
   }

}

/**综合起来，十种差分更新策略**/
/* DE/rand/1/bin */
void DE_rand_1_bin(int function){
    mutate_rand_1(function);
    bin_cross(function);
}

/* DE/rand/1/exp */
void DE_rand_1_exp(int function){
    mutate_rand_1(function);
    exp_cross(function);
}


/* DE/best/1/bin */
void DE_best_1_bin(int function){

    mutate_best_1(function);
    bin_cross(function);
}

/* DE/best/1/exp */
void DE_best_1_exp(int function){

    mutate_best_1(function);
    exp_cross(function);
}

/* DE/rand/2/bin */
void DE_rand_2_bin(int function){

   mutate_rand_2(function);
   bin_cross(function);

}

/* DE/rand/2/exp */
void DE_rand_2_exp(int function){

  mutate_rand_2(function);
  exp_cross(function);
}

/* DE/best/2/bin */
void DE_best_2_bin(int function){

   mutate_best_2(function);
   bin_cross(function);

}

/* DE/best/2/exp */
void DE_best_2_exp(int function){

   mutate_best_2(function);
   exp_cross(function);
}

/* DE/rand-to-best/1/bin */
void DE_randToBest_1_bin(int function){

    mutate_randToBest_1(function);
    bin_cross(function);

}
/* DE/rand-to-best/1/exp */
void DE_randToBest_1_exp(int function){

    mutate_randToBest_1(function);
    exp_cross(function);

}


int DE_strategy(int function,int strategy){

    switch(strategy){
    case 1:
        DE_rand_1_bin(function);
        break;
    case 2:
        DE_rand_1_exp(function);
        break;
    case 3:
        DE_best_1_bin(function);
        break;
    case 4:
        DE_best_1_exp(function);
        break;
    case 5:
        DE_rand_2_bin(function);
        break;
    case 6:
        DE_rand_2_exp(function);
        break;
    case 7:
        DE_best_2_bin(function);
        break;
    case 8:
        DE_best_2_exp(function);
        break;
    case 9:
        DE_randToBest_1_bin(function);
        break;
    case 10:
        DE_randToBest_1_exp(function);
        break;
    default:
		printf("没有该演化策略\n");
		return -1000;


    }


}



/*满足指数分布的交叉操作*/
void exp_cross(int function){
   int i,j,k,randj;
   double randp;
   for(i=0;i<POPSIZE;i++){
     randj=(int)(rand()%NVARS);
     for(j=0;j<=randj;j++){
        if(j<randj){
          crosspopulation[i].gene[j]=population[i].gene[j];
        }else{
          crosspopulation[i].gene[j]=mutatepopulation[i].gene[j];
        }

     }
    for(j=randj+1;j<NVARS;j++){
        randp=((double)(rand()%1001)/1000);
//      printf("输出概率p=%lf\n",randp);
        if(randp<=CR){
            crosspopulation[i].gene[j]=mutatepopulation[i].gene[j];
        }else{
            //当randp>CR时，结束交叉操作
            for(k=j;k<NVARS;k++)
             crosspopulation[i].gene[k]=population[i].gene[k];

             j=NVARS;//退出外出循环
        }

    }
    //计算个体的适应度
    crosspopulation[i].fitness=CalFitness(crosspopulation[i].gene,function);
   }

}



/**选择操作**/
/*从交叉后的种群和目标种群中，选择适应度更高的个体*/

void select_individual(int generation){
    int i;
    for(i=0;i<POPSIZE;i++){
        if(crosspopulation[i].fitness>population[i].fitness)
         population[i].fitness=crosspopulation[i].fitness;
         copyGene(crosspopulation[i].gene,population[i].gene,NVARS);
         keep_best_worst(i,generation);
    }

}


/*在种群population中的最后两个位置分别保存最好个体和最差个体*/
void keep_best_worst(int i,int generation){

    if(population[i].fitness>population[BEST_INDEX].fitness){
       population[BEST_INDEX].fitness=population[i].fitness;
       copyGene(population[i].gene,population[BEST_INDEX].gene,NVARS);
       population[BEST_INDEX].generation=generation;
    }
    if(population[i].fitness<population[WORST_INDEX].fitness)
        population[WORST_INDEX].fitness=population[i].fitness;
        copyGene(population[i].gene,population[WORST_INDEX].gene,NVARS);
}









/*初始化函数中，每个变量的值都是在对应定义域的随机值*/
void initialize(int function){


int i,j;
double lbound,ubound;

   for(i=0;i<POPSIZE;i++)
    {
     //初始化变量
      for(j=0;j<NVARS;j++)
      {

            population[i].gene[j]=randval(bound[function].lbound,bound[function].ubound);

      }
     //初始化个体的适应度
    population[i].fitness= CalFitness(population[i].gene,function);
    //初始化适应度最好和最差个体
    if(i==0){
         population[BEST_INDEX].fitness=population[i].fitness;
         population[WORST_INDEX].fitness=population[i].fitness;
         copyGene( population[i].gene, population[BEST_INDEX].gene,NVARS);
         copyGene( population[i].gene, population[WORST_INDEX].gene,NVARS);
    }else{
      keep_best_worst(i,1);
    }

    }


}


/*从文件中获取测试函数变量的上下界*/
void initializeBound(void)
{
    FILE *infile;
    int i,j;
    double lbound,ubound;

    if((infile=fopen("bound.txt","r"))==NULL)
    {
        //fprintf(delog,"\nCannot open input file!\n");
        printf("打开文件失败\n");
        exit(1);
    }
     //printf("打开文件成功\n");
   /*从文件中获取上界下界并初始化参数*/
   for(i=1;i<=FUNCTIONNUM;i++)
    {

       fscanf(infile,"%lf",&bound[i].lbound);
       fscanf(infile,"%lf",&bound[i].ubound);
       fscanf(infile,"%d",&bound[i].nvars);

    }
    fclose(infile);

}


/*通过上界下界产生随机值*/
double randval(double low,double high)
{
    double val;
    val=((double)(rand()%1001)/1000)*(high-low)+low;
    return(val);
}





//void evaluate_population(int founction,struct genotype * population){
//
//int mem;
//int i;
//double x[NVARS+1];
//for(mem=0;mem<POPSIZE;mem++)
//  {
//     population[mem].fitness=CalFitness(population[mem].gene,founction);
//
//  }
//
//}
//




/*将数组a复制到数组b*/
void copyGene(double *a,double *b,int num){
    int i;
    for(i=0;i<num;i++){
        b[i]=a[i];
    }
}




/*打印报告的头部*/
void printHead(int function,int strategy)
{


    if((delog=fopen(path[strategy].path,"a"))==NULL)
    {
        exit(1);
    }
    //fseek(galog,0,SEEK_END);
    fprintf(delog,"\nFno          [low,up]                  nvars        pop       CR         F     TRIALNUM     maxGen       \tDEStrategy\n");
//    fprintf(galog,"  %d         [%6.3f,%6.3f]            %d           %d     %6.1f      %6.1f        %d              %d\n",i,bound[i].lbound,bound[i].ubound,bound[i].nvars,POPSIZE,PXOVER,PMUTATION ,TRAILNUM,MAXGENS);
    fprintf(delog,"  %d         [%6.3f,%6.3f]   ",function,bound[function].lbound,bound[function].ubound);
   if(function==1||function==3||function==4||function==6||function==8||function==11)
    {
    fprintf(delog,"     %d           %d     %6.1f      %6.1f        %d              %d\t\t%s\n",bound[function].nvars,POPSIZE,CR,F ,TRAILNUM,MAXGENS,path[strategy].strategyname);
    }
   else
    {
        fprintf(delog,"        %d           %d     %6.1f      %6.1f        %d              %d\t\t%s\n",bound[function].nvars,POPSIZE,CR,F,TRAILNUM,MAXGENS,path[strategy].strategyname);

    }
    fprintf(delog,"\n");
    fprintf(delog,"trial          \t\tbest_Val\t\t\t\t              best_Gen\n");
    fclose(delog);
}


/*打印测试序号，最好适应度，最好适应度对应的代数*/

void printf_trailNo_bestVal_bestGen(int trailNo,int founction,int strategy)
{

    if((delog=fopen(path[strategy].path,"a"))==NULL)
    {
        exit(1);
    }

     double best_fitness;

         best_fitness=population[BEST_INDEX].fitness;


//   //将每次测试的最优值保存起来
//    best_val[trailNo]=best_fitness;

       //fprintf(delog,"trial         best_Val              best_Gen\n");
        fprintf(delog,"%d         \t\t%lf\t\t\t              %d\n",trailNo,best_fitness,population[BEST_INDEX].generation);
        fclose(delog);

}


/*打印平均值,标准差,最差个体的适应度，获取最优值的平均代数*/
void printf_mean_stddev_worst_avggen(int strategy)
{
     if((delog=fopen(path[strategy].path,"a"))==NULL)
    {


        exit(1);
    }
        fprintf(delog,"best:%6.3lf \tworst:%6.3lf \t\t mean:%6.3lf  \tstdev:%6.3f \tavgGen:%lf\n",best_val[0],worst_val[0],mean,stddev,(accu_best_gen/TRAILNUM));
        fprintf(delog,"*********************************************************************************************************************");
        fclose(delog);
}
//保存每次测试的最好个体和最差个体的适应度
void save_multipletrail_best_worst(int trailNo)
{
    best_val[trailNo]=population[BEST_INDEX].fitness;
    worst_val[trailNo]=population[WORST_INDEX].fitness;
    if(trailNo==1){
        best_val[0]=population[BEST_INDEX].fitness;
        worst_val[0]=population[WORST_INDEX].fitness;
    }
    //0号位置的值，是多次测试中最优值的最优值和最差值的最差值
    if(best_val[trailNo]>best_val[0]){
        best_val[0]=best_val[trailNo];
    }
     if(worst_val[trailNo]<worst_val[0]){
        worst_val[0]=worst_val[trailNo];
    }
}

//计算最优适应度的平均值和标准差
void get_mean_stdDev()
{
    int i;

    double sum=0.0;/*总的适应度*/
    double std=0.0;

    double fitness;

    for(i=1;i<=TRAILNUM;i++)
    {
        sum+=best_val[i];

    }
    mean=sum/TRAILNUM;

    for(i=1;i<=TRAILNUM;i++)
    {
        std+=pow(best_val[i]-mean,2);
    }
    std/=TRAILNUM-1;
    stddev=sqrt(std);

}
/*清空文件filename*/
int  clean_file(char *filename){
    //清空文件galog.txt再操作
   if((delog=fopen(filename,"w"))==NULL)
     {
        return 0;
     }
     fclose(delog);
     return 1;
}



/*单次测试函数*/
void  DE(int trail,int function,int strategy)
{
    int i;
    generation=1;
    initialize(function);
    while(generation<MAXGENS)
    {
        generation++;
        DE_strategy(function,strategy);
        select_individual(generation);

    }
    printf_trailNo_bestVal_bestGen(trail,function,strategy);

//    printf("当前测试次数：%d 最好适应度：%lf   获得最好适应度的代数%d\n",trail,population[BEST_INDEX].fitness,population[BEST_INDEX].generation);
    accu_best_gen+=population[BEST_INDEX].generation;

    printf("当前的测试序号：%d\n",trail);
    if(trail==TRAILNUM){
         printf("\n****第%d个函数已经测试完成****\n",function);
    }
//    printf("最差的个体的适应度是：%6.3lf\n",population[POPSIZE+1].fitness);

}

/*用函数function测试DE算法num次,用的差分策略是strategy*/
void  multiple_DE(int num ,int function,int strategy){

   int i;
   accu_best_gen=0;
   printHead(function,strategy);
   for(i=1;i<=num;i++)
    {
         DE(i,function,strategy);
         save_multipletrail_best_worst(i);
    }

    //计算十次最优值的标准差和平均值
    get_mean_stdDev();


}




/*每个函数测试num次*/
void test_DE(int num,int strategy){
    int i;
    if(!clean_file(path[strategy].path)){
        printf("文件%s打开失败\n",DELOG);
        return 0;
    }
    initializeBound();//初始化各个函数的定义域
    for(i=1;i<=FUNCTIONNUM;i++){
        if(i==19||i==20)
        continue;
        multiple_DE(num,i,strategy);
        printf_mean_stddev_worst_avggen(strategy);

    }

}
/*测试前num个差分策略*/
int test_DE_strategy(int num){
    int i;
    if(num>STRATEGY_NUM){
        printf("差分策略个数有误，最多为%d\n",STRATEGY_NUM);
        return 0;
    }
    initialize_strategy_path();
    for(i=1;i<=num;i++){
        test_DE(TRAILNUM,i);
    }


}



/*测试函数，打印种群中的个体信息*/
void printfPopulation(struct genotype* population,int pop_num,int var_num){
int i;
int j;
for(i=0;i<pop_num;i++)
  {
    printf("个体%d：",i+1);
    for(j=0;j<var_num;j++)
        printf("变量%d为：%lf \t",j+1,population[i].gene[j]);

    //printf("相对适应度为:%lf",population[i].rfitness);
     printf("\n适应度为:%lf",population[i].fitness);
    printf("\n");
  }
}

//测试函数，打印最差个体的数组
void printfWorst(){
int i;
int j;
for(i=1;i<=TRAILNUM;i++)
  {
    printf("第%d次测试的：",i);
//    for(j=0;j<NVARS;j++)
//        printf("变量%d为：%lf \t",j+1,population[i].gene[j]);

    //printf("相对适应度为:%lf",population[i].rfitness);
     printf("最差的适应度为:%lf\n",worst_val[i]);
    printf("\n");
  }
}
//测试函数，打印最好个体的数组
void printfBest(){
int i;
int j;
for(i=1;i<=TRAILNUM;i++)
  {
    printf("第%d次测试的：",i);
//    for(j=0;j<NVARS;j++)
//        printf("变量%d为：%lf \t",j+1,population[i].gene[j]);

    //printf("相对适应度为:%lf",population[i].rfitness);
     printf("最好的适应度为:%lf\n",best_val[i]);
    printf("\n");
  }
}

//测试函数，打印population数组中的gene数组

void printfGene(struct genotype * population){

    double *gene;
    int i=0;
    gene=population->gene;
    for(int i=0;i<NVARS;i++){
        printf("输出变量%d的值为%lf\n",i,gene[i]);

    }

}




int main()
{

    srand((unsigned)time(NULL));
//    initializeBound();
//    initialize(1);
//    mutate_rand_1(1);
//    bin_cross(1);

//    printf("输出population种群\n");
//    printf_pop(population,2);
//    printf("\n输出mutatepopulation种群\n");
//    printf_pop(mutatepopulation,2);
//    printf("\n输出crosspopulation种群\n");
//    printf_pop(crosspopulation,2);

//    evaluate_population(1,population);
//    evaluate_population(1,crosspopulation);
//    evaluate_population(1,mutatepopulation);
//
//    multiple_DE(1,1);
//    printf("输出population种群\n");
//    printfPopulation(population,ARRAYSIZE,0);
//    printf("\n输出crosspopulation种群\n");
//    printfPopulation(crosspopulation,POPSIZE,0);
//    select_individual(1);
//    printf("选择后输出population种群\n");
//    printfPopulation(population,ARRAYSIZE,0);
//    printf("输出mutatepopulation种群\n");
//    printfPopulation(mutatepopulation,2,NVARS);

// int boolean=clean_file(DELOG);
//printf("%d",boolean);

//  test_DE(TRAILNUM,1);
// initializeBound();
// initialize(1);
// mutate_rand_1(1);
// exp_cross(1);
//    printf("输出population种群\n");
//    printfPopulation(population,1,NVARS);
//    printf("输出mutatepopulation种群\n");
//    printfPopulation(mutatepopulation,1,NVARS);
//    printf("\n输出crosspopulation种群\n");
//    printfPopulation(crosspopulation,1,NVARS);
//
test_DE_strategy(STRATEGY_NUM);

    return 0;
}
