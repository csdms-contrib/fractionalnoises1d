#ifdef HAVE_MALLOC_H
# include<malloc.h>
#endif
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define FREE_ARG char*
#define NR_END 1

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
        return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(idum)
int *idum;
{
        static int inext,inextp;
        static long ma[56];
        static int iff=0;
        long mj,mk;
        int i,ii,k;

        if (*idum < 0 || iff == 0) {
                iff=1;
                mj=MSEED-(*idum < 0 ? -*idum : *idum);
                mj %= MBIG;
                ma[55]=mj;
                mk=1;
                for (i=1;i<=54;i++) {
                        ii=(21*i) % 55;
                        ma[ii]=mk;
                        mk=mj-mk;
                        if (mk < MZ) mk += MBIG;
                        mj=ma[ii];
                }
                for (k=1;k<=4;k++)
                        for (i=1;i<=55;i++) {
                                ma[i] -= ma[1+(i+30) % 55];
                                if (ma[i] < MZ) ma[i] += MBIG;
                        }
                inext=0;
                inextp=31;
                *idum=1;
        }
        if (++inext == 56) inext=1;
        if (++inextp == 56) inextp=1;
        mj=ma[inext]-ma[inextp];
        if (mj < MZ) mj += MBIG;
        ma[inext]=mj;
        return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

float gasdev(idum)
int *idum;
{
        static int iset=0;
        static float gset;
        float fac,r,v1,v2;
        float ran3();

        if  (iset == 0) {
                do {
                        v1=2.0*ran3(idum)-1.0;
                        v2=2.0*ran3(idum)-1.0;
                        r=v1*v1+v2*v2;
                } while (r >= 1.0);
                fac=sqrt(-2.0*log(r)/r);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset=0;
                return gset;
        }
}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(float data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
#undef SWAP

void realft(float data[], unsigned long n, int isign)
{
	void four1(float data[], unsigned long nn, int isign);
	unsigned long i,i1,i2,i3,i4,np3;
	float c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592653589793/(double) (n>>1);
	if (isign == 1) {
		c2 = -0.5;
		four1(data,n>>1,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) {
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	} else {
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		four1(data,n>>1,-1);
	}
}

main()
{    float *precyear,beta,sum,sddevcomp,mean,sddev;
     int length,year,format; 
     int idum;
     char outputfile[30];
     FILE *fp;   
     
     printf("\nName of output file: ");
     scanf("%s",outputfile);
     printf("\nLength of desired noise: (factor of two, please) ");
     scanf("%d",&length);
     printf("\nGaussian or Log-normal dataset? (1 for gaussian 2 for l-n): ");
     scanf("%d",&format);
     printf("\nMean: (if log-normal this mean is that of the final, log transformed data) ");
     scanf("%f",&mean);
     printf("\nStandard Dev: ");
     scanf("%f",&sddev);
     printf("\nBeta: ");
     scanf("%f",&beta);
     printf("\nRandom number seed: (any negative integer) ");
     scanf("%d",&idum);      
     fp=fopen(outputfile,"w");
     precyear=vector(1,2*length);
     /*we contruct a vector of length two times the input length so that we
     may cut off the ends to eliminate the periodicity introduced by filtering*/
     for (year=1;year<=2*length;year++) 
       precyear[year]=gasdev(&idum);    
     realft(precyear,2*length,1);
     precyear[1]=0.0;
     precyear[2]=precyear[2]/pow(0.5,beta/2.0);
     for (year=3;year<=2*length;year++) 
       precyear[year]=precyear[year]/pow((year/2)/(float)(2*length),beta/2.0);
     realft(precyear,2*length,-1);
     sum=0.0;
     for (year=length/2;year<=length+length/2;year++) 
        {sum+=precyear[year];}
     sum=sum/length; 
     sddevcomp=0.0;
     for (year=length/2;year<=length+length/2;year++)
       sddevcomp+=(precyear[year]-sum)*(precyear[year]-sum); 
     sddevcomp=sddevcomp/length;
     /*below we rescale the amplitudes to restore the desired moments*/
     for (year=length/2;year<=length+length/2;year++)
      {precyear[year]=(precyear[year]-sum)/sqrt(sddevcomp);
       precyear[year]=(precyear[year]*sddev)+mean;
       if (format==2) if (mean!=0.0) precyear[year]=
        sqrt(log(1+sddev*sddev/(mean*mean)))*precyear[year]+
         log(mean/sqrt(1+sddev*sddev/(mean*mean)));}
     if (format==2) for (year=length/2;year<=length+length/2;year++)
      precyear[year]=exp(precyear[year]);
     for (year=length/2;year<=length+length/2;year++)
      fprintf(fp,"%d %f\n",year-length/2,precyear[year]+(float)(year-length/2)/(length));
     fclose(fp);
}   
