/*���Ŵ��㷨��,��Ӧ��ֻȡ����,ͬʱ,
ÿ���������Ӧ�ȶ���Ŀ�꺯���ĺ���ֵ��ͬ
���Ժ����ı���ȡֵ��Χ���ļ���D://bound.txt���ж�ȡ
�㷨��������������ļ���D://delog.txt"

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<string.h>

#define POPSIZE 10   /*���������*/
#define ARRAYSIZE POPSIZE+2 /*��Ⱥ����Ĵ�С,�����������λ��,����������ø����������*/
#define MAXGENS 3000 /*���ĵ�������*/
#define NVARS 30      /*������ά��,������δ֪������*/
#define dim   NVARS    /*���Ժ����еĺ���ά��*/
#define CR  0.9     /*������������ĸ���*/
#define FUNCTIONNUM 23
#define TRAILNUM     10
#define TRUE 1
#define FALSE 0
#define F 0.5     /*��������Fһ��ȡ[0,2]�ĳ���*/
#define LAMBDA F  /*��ֵ����������ȡֵһ��*/
#define BEST_INDEX POPSIZE
#define WORST_INDEX POPSIZE+1
#define DELOG "DEreport\\DE_rand_1_bin.txt"
#define STRATEGY_NUM 10 /*10����ֲ���*/



int generation;/*��ǰ�Ĵ���*/
FILE *delog ;    /*�ļ�ָ��*/
double best_val[TRAILNUM+1];//ʮ�β��Ե�����ֵ
double worst_val[TRAILNUM+1];//ʮ�β��Ե����ֵ
double mean;//��������ʮ�β�������ֵ��ƽ��ֵ
double stddev;//��������ʮ�β�������ֵ�ı�׼��
double accu_best_gen;//�����ۼ�ʮ�������Ÿ���ĵ�������

struct genotype /*������,����һ������*/
{
    double gene[NVARS]; //��������ÿ�����򣬼�ÿ������
    double fitness; //��ǰ�������Ӧ��
    double upper[NVARS];    //������ȡֵ��Χ�е��Ͻ磬��ͬ�±��ֵ�����Ӧ�ı������Ͻ�
    double lower[NVARS];    //������ȡֵ��Χ�е��½磬��ͬ�±��ֵ�����Ӧ�ı������½�
    int generation;  //�õ����Ÿ���ʱ�ĵ������
};
struct genotype population[ARRAYSIZE]; //�����������и��壬�������λ������������ú�������
struct genotype newpopulation[ARRAYSIZE];//�����������������ɵĸ���
struct genotype mutatepopulation[ARRAYSIZE];
struct genotype crosspopulation[ARRAYSIZE];

/*����������Ժ����ı���ȡֵ���Ͻ���½�*/
struct boundtype
{
    double lbound;
    double ubound;
    int    nvars;//��ǰ���Ժ�����ά��
};
struct boundtype  bound[FUNCTIONNUM+1];

/*��ͬ��ֲ��ԵĲ��Ա����ӡ�ڲ�ͬ�ļ�*/
struct printfpath
{
    char* strategyname;
    char* path;
};
struct printfpath path[STRATEGY_NUM+1];
/************************************���Ժ���**********************************/
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

/*ѡ����Ժ���*/
double CalFitness(double x[],int fun)
/******************************************
***���룺һ����Ч��************************
***������ý����Ӧֵ**********************
***��ע���������Ӧֵ�����������**********
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


double randval(double,double);//�������ֵ�ĺ���
void initialize(int function);//��ʼ������
void initializeBound(void);//��ʼ�������Ķ�����
void copyGene(double *a,double *b,int num);//������a���Ƶ�����b
void printHead(int function,int strategy);//���ļ���������Ա����ͷ����Ϣ
void printf_trailNo_bestVal_bestGen(int trailNo,int founction,int strategy);//
void save_multipletrail_best_worst(int trailNo);//����ÿ�β��Ե���ø�������������Ӧ��
void printf_mean_stddev_worst_avggen(int strategy);//��ӡƽ��ֵ,��׼��,���������Ӧ�ȣ���ȡ����ֵ��ƽ������
void  DE(int trail,int function,int strategy);//���β��Ժ���
void  multiple_DE(int num ,int function,int strategy);//��DE�㷨���Ե�function������num��
void test_DE(int num,int strategy);//ÿ����������num��
void get_mean_stdDev();//����������Ӧ�ȵ�ƽ��ֵ�ͱ�׼��
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
    //���D:/DEreport�ļ����²����������ļ����½���ע��D����Ҫ���д���DEreport�ļ���
    for(i=1;i<=STRATEGY_NUM;i++){
        if((fp=fopen(path[i].path,"r"))==NULL) //�ж��ļ��Ƿ����
        {
        fp=fopen(path[i].path,"w"); //����ļ������ڣ��½�һ���ļ�
        }
        fclose(fp); //�ر��ļ�
    }
}


/***���ֱ������***/

void mutate_rand_1(int function){
    int i,j;
    int rand1,rand2,rand3;
    for(i=0;i<POPSIZE;i++){
      rand1=(int)rand()%POPSIZE;
      rand2=(int)rand()%POPSIZE;
      rand3=(int)rand()%POPSIZE;
      //i�뼸���������������ͬ
      while(rand1==rand2||rand1==rand3||rand2==rand3||rand1==i||rand2==i||rand3==i)
     {
       rand1=(int)rand()%POPSIZE;
       rand2=(int)rand()%POPSIZE;
       rand3=(int)rand()%POPSIZE;
    //   printf("����������%d,%d,%d,%d\n",i,rand1,rand2,rand3);
      }
      for(j=0;j<NVARS;j++){

        mutatepopulation[i].gene[j]=population[rand1].gene[j]+F*(population[rand2].gene[j]-population[rand3].gene[j]);

         //������ֵ���ܳ���������
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
      //i�뼸���������������ͬ
      while(rand1==rand2||rand1==i||rand2==i)
     {
       rand1=(int)rand()%POPSIZE;
       rand2=(int)rand()%POPSIZE;
    //   printf("����������%d,%d,%d,%d\n",i,rand1,rand2,rand3);
      }
      for(j=0;j<NVARS;j++){

        mutatepopulation[i].gene[j]=population[BEST_INDEX].gene[j]+F*(population[rand1].gene[j]-population[rand2].gene[j]);

         //������ֵ���ܳ���������
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
      //i�뼸���������������ͬ,
      while(rand1==rand2||rand1==rand3||rand2==rand3||rand1==i||rand2==i||rand3==i||rand1==rand4||rand1==rand5||rand2==rand4||rand2==rand5||rand3==rand4||rand3==rand5||rand4==rand5||i==rand4||i==rand5)
     {
       rand1=(int)rand()%POPSIZE;
       rand2=(int)rand()%POPSIZE;
       rand3=(int)rand()%POPSIZE;
       rand4=(int)rand()%POPSIZE;
       rand5=(int)rand()%POPSIZE;
    //   printf("����������%d,%d,%d,%d\n",i,rand1,rand2,rand3);
      }
      for(j=0;j<NVARS;j++){

        mutatepopulation[i].gene[j]=population[rand1].gene[j]+F*(population[rand2].gene[j]+population[rand3].gene[j]-population[rand4].gene[j]-population[rand5].gene[j]);

         //������ֵ���ܳ���������
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

      //i�뼸���������������ͬ,
      while(rand1==rand2||rand1==rand3||rand2==rand3||rand1==i||rand2==i||rand3==i||rand1==rand4||rand2==rand4||rand3==rand4||i==rand4)
     {
       rand1=(int)rand()%POPSIZE;
       rand2=(int)rand()%POPSIZE;
       rand3=(int)rand()%POPSIZE;
       rand4=(int)rand()%POPSIZE;

    //   printf("����������%d,%d,%d,%d\n",i,rand1,rand2,rand3);
      }
      for(j=0;j<NVARS;j++){

        mutatepopulation[i].gene[j]=population[BEST_INDEX].gene[j]+F*(population[rand1].gene[j]+population[rand2].gene[j]-population[rand3].gene[j]-population[rand4].gene[j]);

         //������ֵ���ܳ���������
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

      //i�뼸���������������ͬ
      while(rand1==rand2||rand1==i||rand2==i)
     {
       rand1=(int)rand()%POPSIZE;
       rand2=(int)rand()%POPSIZE;
    //   printf("����������%d,%d,%d,%d\n",i,rand1,rand2,rand3);
      }
      for(j=0;j<NVARS;j++){

        mutatepopulation[i].gene[j]=population[i].gene[j]+LAMBDA*(population[BEST_INDEX].gene[j]-population[i].gene[j])+F*(population[rand1].gene[j]-population[rand2].gene[j]);

         //������ֵ���ܳ���������
         if(mutatepopulation[i].gene[j]>bound[function].ubound){

            mutatepopulation[i].gene[j]=bound[function].ubound;
      }

         if(mutatepopulation[i].gene[j]<bound[function].lbound){

            mutatepopulation[i].gene[j]=bound[function].lbound;
         }
      }

    }

}



/***���ֽ������***/
/*�������ʽ�ֲ��Ľ������*/
void bin_cross(int function){
   int i,j,randj;
   double randp;
   for(i=0;i<POPSIZE;i++){
     randj=(int)(rand()%NVARS);
    for(j=0;j<NVARS;j++){
           randp=((double)(rand()%1001)/1000);
//           printf("�������p=%lf\n",randp);
        if(randp<=CR||randj==j){
            crosspopulation[i].gene[j]=mutatepopulation[i].gene[j];
        }else{
            crosspopulation[i].gene[j]=population[i].gene[j];
        }

    }
    //����������Ӧ��
    crosspopulation[i].fitness=CalFitness(crosspopulation[i].gene,function);
   }

}

/**�ۺ�������ʮ�ֲ�ָ��²���**/
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
		printf("û�и��ݻ�����\n");
		return -1000;


    }


}



/*����ָ���ֲ��Ľ������*/
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
//      printf("�������p=%lf\n",randp);
        if(randp<=CR){
            crosspopulation[i].gene[j]=mutatepopulation[i].gene[j];
        }else{
            //��randp>CRʱ�������������
            for(k=j;k<NVARS;k++)
             crosspopulation[i].gene[k]=population[i].gene[k];

             j=NVARS;//�˳����ѭ��
        }

    }
    //����������Ӧ��
    crosspopulation[i].fitness=CalFitness(crosspopulation[i].gene,function);
   }

}



/**ѡ�����**/
/*�ӽ�������Ⱥ��Ŀ����Ⱥ�У�ѡ����Ӧ�ȸ��ߵĸ���*/

void select_individual(int generation){
    int i;
    for(i=0;i<POPSIZE;i++){
        if(crosspopulation[i].fitness>population[i].fitness)
         population[i].fitness=crosspopulation[i].fitness;
         copyGene(crosspopulation[i].gene,population[i].gene,NVARS);
         keep_best_worst(i,generation);
    }

}


/*����Ⱥpopulation�е��������λ�÷ֱ𱣴���ø����������*/
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









/*��ʼ�������У�ÿ��������ֵ�����ڶ�Ӧ����������ֵ*/
void initialize(int function){


int i,j;
double lbound,ubound;

   for(i=0;i<POPSIZE;i++)
    {
     //��ʼ������
      for(j=0;j<NVARS;j++)
      {

            population[i].gene[j]=randval(bound[function].lbound,bound[function].ubound);

      }
     //��ʼ���������Ӧ��
    population[i].fitness= CalFitness(population[i].gene,function);
    //��ʼ����Ӧ����ú�������
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


/*���ļ��л�ȡ���Ժ������������½�*/
void initializeBound(void)
{
    FILE *infile;
    int i,j;
    double lbound,ubound;

    if((infile=fopen("bound.txt","r"))==NULL)
    {
        //fprintf(delog,"\nCannot open input file!\n");
        printf("���ļ�ʧ��\n");
        exit(1);
    }
     //printf("���ļ��ɹ�\n");
   /*���ļ��л�ȡ�Ͻ��½粢��ʼ������*/
   for(i=1;i<=FUNCTIONNUM;i++)
    {

       fscanf(infile,"%lf",&bound[i].lbound);
       fscanf(infile,"%lf",&bound[i].ubound);
       fscanf(infile,"%d",&bound[i].nvars);

    }
    fclose(infile);

}


/*ͨ���Ͻ��½�������ֵ*/
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




/*������a���Ƶ�����b*/
void copyGene(double *a,double *b,int num){
    int i;
    for(i=0;i<num;i++){
        b[i]=a[i];
    }
}




/*��ӡ�����ͷ��*/
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


/*��ӡ������ţ������Ӧ�ȣ������Ӧ�ȶ�Ӧ�Ĵ���*/

void printf_trailNo_bestVal_bestGen(int trailNo,int founction,int strategy)
{

    if((delog=fopen(path[strategy].path,"a"))==NULL)
    {
        exit(1);
    }

     double best_fitness;

         best_fitness=population[BEST_INDEX].fitness;


//   //��ÿ�β��Ե�����ֵ��������
//    best_val[trailNo]=best_fitness;

       //fprintf(delog,"trial         best_Val              best_Gen\n");
        fprintf(delog,"%d         \t\t%lf\t\t\t              %d\n",trailNo,best_fitness,population[BEST_INDEX].generation);
        fclose(delog);

}


/*��ӡƽ��ֵ,��׼��,���������Ӧ�ȣ���ȡ����ֵ��ƽ������*/
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
//����ÿ�β��Ե���ø�������������Ӧ��
void save_multipletrail_best_worst(int trailNo)
{
    best_val[trailNo]=population[BEST_INDEX].fitness;
    worst_val[trailNo]=population[WORST_INDEX].fitness;
    if(trailNo==1){
        best_val[0]=population[BEST_INDEX].fitness;
        worst_val[0]=population[WORST_INDEX].fitness;
    }
    //0��λ�õ�ֵ���Ƕ�β���������ֵ������ֵ�����ֵ�����ֵ
    if(best_val[trailNo]>best_val[0]){
        best_val[0]=best_val[trailNo];
    }
     if(worst_val[trailNo]<worst_val[0]){
        worst_val[0]=worst_val[trailNo];
    }
}

//����������Ӧ�ȵ�ƽ��ֵ�ͱ�׼��
void get_mean_stdDev()
{
    int i;

    double sum=0.0;/*�ܵ���Ӧ��*/
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
/*����ļ�filename*/
int  clean_file(char *filename){
    //����ļ�galog.txt�ٲ���
   if((delog=fopen(filename,"w"))==NULL)
     {
        return 0;
     }
     fclose(delog);
     return 1;
}



/*���β��Ժ���*/
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

//    printf("��ǰ���Դ�����%d �����Ӧ�ȣ�%lf   ��������Ӧ�ȵĴ���%d\n",trail,population[BEST_INDEX].fitness,population[BEST_INDEX].generation);
    accu_best_gen+=population[BEST_INDEX].generation;

    printf("��ǰ�Ĳ�����ţ�%d\n",trail);
    if(trail==TRAILNUM){
         printf("\n****��%d�������Ѿ��������****\n",function);
    }
//    printf("���ĸ������Ӧ���ǣ�%6.3lf\n",population[POPSIZE+1].fitness);

}

/*�ú���function����DE�㷨num��,�õĲ�ֲ�����strategy*/
void  multiple_DE(int num ,int function,int strategy){

   int i;
   accu_best_gen=0;
   printHead(function,strategy);
   for(i=1;i<=num;i++)
    {
         DE(i,function,strategy);
         save_multipletrail_best_worst(i);
    }

    //����ʮ������ֵ�ı�׼���ƽ��ֵ
    get_mean_stdDev();


}




/*ÿ����������num��*/
void test_DE(int num,int strategy){
    int i;
    if(!clean_file(path[strategy].path)){
        printf("�ļ�%s��ʧ��\n",DELOG);
        return 0;
    }
    initializeBound();//��ʼ�����������Ķ�����
    for(i=1;i<=FUNCTIONNUM;i++){
        if(i==19||i==20)
        continue;
        multiple_DE(num,i,strategy);
        printf_mean_stddev_worst_avggen(strategy);

    }

}
/*����ǰnum����ֲ���*/
int test_DE_strategy(int num){
    int i;
    if(num>STRATEGY_NUM){
        printf("��ֲ��Ը����������Ϊ%d\n",STRATEGY_NUM);
        return 0;
    }
    initialize_strategy_path();
    for(i=1;i<=num;i++){
        test_DE(TRAILNUM,i);
    }


}



/*���Ժ�������ӡ��Ⱥ�еĸ�����Ϣ*/
void printfPopulation(struct genotype* population,int pop_num,int var_num){
int i;
int j;
for(i=0;i<pop_num;i++)
  {
    printf("����%d��",i+1);
    for(j=0;j<var_num;j++)
        printf("����%dΪ��%lf \t",j+1,population[i].gene[j]);

    //printf("�����Ӧ��Ϊ:%lf",population[i].rfitness);
     printf("\n��Ӧ��Ϊ:%lf",population[i].fitness);
    printf("\n");
  }
}

//���Ժ�������ӡ�����������
void printfWorst(){
int i;
int j;
for(i=1;i<=TRAILNUM;i++)
  {
    printf("��%d�β��Եģ�",i);
//    for(j=0;j<NVARS;j++)
//        printf("����%dΪ��%lf \t",j+1,population[i].gene[j]);

    //printf("�����Ӧ��Ϊ:%lf",population[i].rfitness);
     printf("������Ӧ��Ϊ:%lf\n",worst_val[i]);
    printf("\n");
  }
}
//���Ժ�������ӡ��ø��������
void printfBest(){
int i;
int j;
for(i=1;i<=TRAILNUM;i++)
  {
    printf("��%d�β��Եģ�",i);
//    for(j=0;j<NVARS;j++)
//        printf("����%dΪ��%lf \t",j+1,population[i].gene[j]);

    //printf("�����Ӧ��Ϊ:%lf",population[i].rfitness);
     printf("��õ���Ӧ��Ϊ:%lf\n",best_val[i]);
    printf("\n");
  }
}

//���Ժ�������ӡpopulation�����е�gene����

void printfGene(struct genotype * population){

    double *gene;
    int i=0;
    gene=population->gene;
    for(int i=0;i<NVARS;i++){
        printf("�������%d��ֵΪ%lf\n",i,gene[i]);

    }

}




int main()
{

    srand((unsigned)time(NULL));
//    initializeBound();
//    initialize(1);
//    mutate_rand_1(1);
//    bin_cross(1);

//    printf("���population��Ⱥ\n");
//    printf_pop(population,2);
//    printf("\n���mutatepopulation��Ⱥ\n");
//    printf_pop(mutatepopulation,2);
//    printf("\n���crosspopulation��Ⱥ\n");
//    printf_pop(crosspopulation,2);

//    evaluate_population(1,population);
//    evaluate_population(1,crosspopulation);
//    evaluate_population(1,mutatepopulation);
//
//    multiple_DE(1,1);
//    printf("���population��Ⱥ\n");
//    printfPopulation(population,ARRAYSIZE,0);
//    printf("\n���crosspopulation��Ⱥ\n");
//    printfPopulation(crosspopulation,POPSIZE,0);
//    select_individual(1);
//    printf("ѡ������population��Ⱥ\n");
//    printfPopulation(population,ARRAYSIZE,0);
//    printf("���mutatepopulation��Ⱥ\n");
//    printfPopulation(mutatepopulation,2,NVARS);

// int boolean=clean_file(DELOG);
//printf("%d",boolean);

//  test_DE(TRAILNUM,1);
// initializeBound();
// initialize(1);
// mutate_rand_1(1);
// exp_cross(1);
//    printf("���population��Ⱥ\n");
//    printfPopulation(population,1,NVARS);
//    printf("���mutatepopulation��Ⱥ\n");
//    printfPopulation(mutatepopulation,1,NVARS);
//    printf("\n���crosspopulation��Ⱥ\n");
//    printfPopulation(crosspopulation,1,NVARS);
//
test_DE_strategy(STRATEGY_NUM);

    return 0;
}
