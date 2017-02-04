/*
Sheffield Auditory Model (First developed at Univ. of Sheffield).
Adapted by Nicoleta Roman, DeLiang Wang & Guy J. Brown (2003) to perform binaural feature extraction.
Useful Outputs:
  
   namefile_itd, namefile_iid, namefile_ratio = to store ITD/IID and a priori relative strength estimates for individual mixtures.
   output_itd, output_iid, output_ratio = to store ITD/IID and a priori relative strength estimates for all the training data.
  
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ieee754.h>
#include "common.h"

typedef struct
{
	double cf, bw, criticalRate, z, gain, expCoeff;
	double midEarCoeff;
	double p0, p1, p2, p3, p4;
	double q0, q1, q2, q3, q4;
	double u0, u1;
	double v0, v1;
	double c,q,w;
	double rate;
	double buffer[MAX_BUFFER_SIZE];
	int ptr;
} channel;


/* function prototypes */

void initHairCells(void);
/* initialise the meddis hair cell parameters */

void initOuterMiddleEar(void);
/* set parameters of equal-loudness functions */

float DBtoAmplitude(float dB);
/* convert dB to amplitude */

float loudnessLevelInPhons(float dB, float freq);
/* compute loudness level */

void initChannels(int lowerCF, int upperCF, int numChannels);
/* initialise filterbank channels */

float updateCochlea(channel *c, float sigval, int tim,int chan);
/* process one sample of the input though the cochlea */

int msToSamples(float ms);
/* converts time in ms to samples at srate */

float compenergy(channel*c);
/* compute energy*/ 

float compute_IID(int chan );
/* compute IID as the energy ratio between right and left ear*/ 

float compute_ITD(int chan);
/* compute cross correlation functions and find maximum peak in a 2PI interval around the target reference */


/* global variables */

channel cochlea_left[MAX_CHANNEL];
channel cochlea_right[MAX_CHANNEL];

channel cochlea_t[MAX_CHANNEL];
channel cochlea_n[MAX_CHANNEL];

float CCF[MAX_CHANNEL];

float t, dt, twoPi, twoPiT, vmin, k, l;
short verboseOutput=FALSE;
float ymdt, xdt, ydt, lplusrdt, rdt, gdt, hdt;


/* outer/middle ear stuff */

float f[MIDDLE_EAR_SIZE];
float af[MIDDLE_EAR_SIZE];
float bf[MIDDLE_EAR_SIZE];
float tf[MIDDLE_EAR_SIZE];


/* function declarations */

void initHairCells(void)
{
	ymdt=MED_Y*MED_M*dt;
	xdt=MED_X*dt;
	ydt=MED_Y*dt;
	lplusrdt=(MED_L+MED_R)*dt;
	rdt=MED_R*dt;
	gdt=MED_G*dt;
	hdt=MED_H; /* should be multiplied by ts really */
}

void initOuterMiddleEar(void)
/*
parameters of equal-loudness functions from BS3383,"Normal equal-loudness level
contours for pure tones under free-field listening conditions", table 1.
f is the tone frequency
af and bf are frequency-dependent coefficients
 tf is the threshold sound pressure level of the tone, in dBs
*/
{
	f[0]=20.0;     af[0]=2.347;  bf[0]=0.00561;   tf[0]=74.3;
	f[1]=25.0;     af[1]=2.190;  bf[1]=0.00527;   tf[1]=65.0;
	f[2]=31.5;     af[2]=2.050;  bf[2]=0.00481;   tf[2]=56.3;
	f[3]=40.0;     af[3]=1.879;  bf[3]=0.00404;   tf[3]=48.4;
	f[4]=50.0;     af[4]=1.724;  bf[4]=0.00383;   tf[4]=41.7;
	f[5]=63.0;     af[5]=1.579;  bf[5]=0.00286;   tf[5]=35.5;
	f[6]=80.0;     af[6]=1.512;  bf[6]=0.00259;   tf[6]=29.8;
	f[7]=100.0;    af[7]=1.466;  bf[7]=0.00257;   tf[7]=25.1;
	f[8]=125.0;    af[8]=1.426;  bf[8]=0.00256;   tf[8]=20.7;
	f[9]=160.0;    af[9]=1.394;  bf[9]=0.00255;   tf[9]=16.8;
	f[10]=200.0;   af[10]=1.372; bf[10]=0.00254;  tf[10]=13.8;
	f[11]=250.0;   af[11]=1.344; bf[11]=0.00248;  tf[11]=11.2;
	f[12]=315.0;   af[12]=1.304; bf[12]=0.00229;  tf[12]=8.9;
	f[13]=400.0;   af[13]=1.256; bf[13]=0.00201;  tf[13]=7.2;
	f[14]=500.0;   af[14]=1.203; bf[14]=0.00162;  tf[14]=6.0;
	f[15]=630.0;   af[15]=1.135; bf[15]=0.00111;  tf[15]=5.0;
	f[16]=800.0;   af[16]=1.062; bf[16]=0.00052;  tf[16]=4.4;
	f[17]=1000.0;  af[17]=1.000; bf[17]=0.00000;  tf[17]=4.2;
	f[18]=1250.0;  af[18]=0.967; bf[18]=-0.00039; tf[18]=3.7;
	f[19]=1600.0;  af[19]=0.943; bf[19]=-0.00067; tf[19]=2.6;
	f[20]=2000.0;  af[20]=0.932; bf[20]=-0.00092; tf[20]=1.0;
	f[21]=2500.0;  af[21]=0.933; bf[21]=-0.00105; tf[21]=-1.2;
	f[22]=3150.0;  af[22]=0.937; bf[22]=-0.00104; tf[22]=-3.6;
	f[23]=4000.0;  af[23]=0.952; bf[23]=-0.00088; tf[23]=-3.9;
	f[24]=5000.0;  af[24]=0.974; bf[24]=-0.00055; tf[24]=-1.1;
	f[25]=6300.0;  af[25]=1.027; bf[25]=0.00000;  tf[25]=6.6;
	f[26]=8000.0;  af[26]=1.135; bf[26]=0.00089;  tf[26]=15.3;
	f[27]=10000.0; af[27]=1.266; bf[27]=0.00211;  tf[27]=16.4;
	f[28]=12500.0; af[28]=1.501; bf[28]=0.00488;  tf[28]=11.6;
}

float DBtoAmplitude(float dB)
{
	return pow(10.0,(dB/20.0));
}

float loudnessLevelInPhons(float dB, float freq)
/*
Uses linear interpolation of the look-up tables to compute the loudness level,
in phons, of a pure tone of frequency freq using the reference curve for sound
 pressure level dB.
The equation is taken from section 4 of BS3383.
*/
{
	int i=0;
	float afy, bfy, tfy;

	if ((freq<20.0) | (freq>12500.0))
	{
		fprintf(stderr,"Can't compute a outer/middle ear gain for that frequency\n");
		exit(0);
	}

	while (f[i] < freq)
		i++;
	
	afy = af[i-1]+(freq-f[i-1])*(af[i]-af[i-1])/(f[i]-f[i-1]);
	bfy = bf[i-1]+(freq-f[i-1])*(bf[i]-bf[i-1])/(f[i]-f[i-1]);
	tfy = tf[i-1]+(freq-f[i-1])*(tf[i]-tf[i-1])/(f[i]-f[i-1]);
	
	return 4.2+afy*(dB-tfy)/(1.0+bfy*(dB-tfy));
}

void initChannels(int lowerCF, int upperCF, int numChannels)
{
	float lowerERB, upperERB, spaceERB, kt;
	channel c;
	int chan,i;

	dt = 1.0/(float)SAMPLING_FREQUENCY;
	twoPi = 2.0*M_PI;
	twoPiT = 2.0*M_PI*dt;

	lowerERB = hzToERBrate(lowerCF);
	upperERB = hzToERBrate(upperCF);

	if (numChannels > 1)
		spaceERB = (upperERB-lowerERB)/(numChannels-1);
	else
		spaceERB = 0.0;

	for (chan=0; chan<numChannels; chan++)
	{
		c.criticalRate = lowerERB+chan*spaceERB;
		c.cf = ERBrateToHz(c.criticalRate);

		CCF[chan]=c.cf;

		c.midEarCoeff=DBtoAmplitude(loudnessLevelInPhons(DB,c.cf)-DB);
		c.bw = erb(c.cf)*BW_CORRECTION;
		c.z = exp(-twoPiT*c.bw);
		c.expCoeff = c.cf*twoPiT;
		c.gain = c.midEarCoeff*sqr(sqr(2*M_PI*c.bw*dt))/3.0;

		if (verboseOutput)
		{
			fprintf(stderr,"cf=%1.4f:\n",c.cf);
			fprintf(stderr,"criticalRate=%1.4f  midEarCoeff=%1.4f  bw=%1.4f  gain=%1.4f\n",c.criticalRate,c.midEarCoeff,c.bw,c.gain);
		}

		c.p0 = 0.0; c.p1 = 0.0; c.p2 = 0.0; c.p3 = 0.0; c.p4 = 0.0;
		c.q0 = 0.0; c.q1 = 0.0; c.q2 = 0.0; c.q3 = 0.0; c.q4 = 0.0;
		c.u0 = 0.0; c.u1 = 0.0;
		c.v0 = 0.0; c.v1 = 0.0;

		kt = MED_G*MED_A/(MED_A+MED_B);
		c.c = MED_M*MED_Y*kt/(MED_L*kt+MED_Y*(MED_L+MED_R));
		c.q = c.c*(MED_L+MED_R)/kt;
		c.w = c.c*MED_R/MED_X;
		c.ptr = MAX_BUFFER_SIZE-1;

		for(i=0;i<MAX_BUFFER_SIZE;i++)
		{
			c.buffer[i]=0.0;
		}

		cochlea_left[chan]=c;
		cochlea_right[chan]=c;
		cochlea_t[chan]=c;
		cochlea_n[chan]=c;
	}
}


float updateCochlea(channel *c, float sigval, int tim, int chan)
{
	double zz, bm, hc, pow, amp;
	double replenish, eject, reuptakeandloss, reuptake, reprocess, kt;
	float filtered_value;

	zz = c->z;
	c->p0 = sigval*cos(c->expCoeff*tim)+zz*(4*c->p1-zz*(6*c->p2-zz*(4*c->p3-zz*c->p4)));
	c->q0 =-sigval*sin(c->expCoeff*tim)+zz*(4*c->q1-zz*(6*c->q2-zz*(4*c->q3-zz*c->q4)));
	c->u0 = zz*(c->p1+zz*(4*c->p2+zz*c->p3));
	c->v0 = zz*(c->q1+zz*(4*c->q2+zz*c->q3));
	bm = (c->u0*cos(c->expCoeff*tim)-c->v0*sin(c->expCoeff*tim))*c->gain;
	pow = sqr(c->u0)+sqr(c->v0);
	amp = sqrt(pow)*c->gain;

	/* hair cell */

	if ((bm+MED_A)>0.0)
		kt=gdt*(bm+MED_A)/(bm+MED_A+MED_B);
	else
		kt=0.0;

	if (c->q<MED_M )
		replenish=ymdt-ydt*c->q;
	else
		replenish=0.0;

	eject = kt*c->q;
	reuptakeandloss = lplusrdt*c->c;
	reuptake = rdt*c->c;
	reprocess = xdt*c->w;
	c->q += replenish-eject+reprocess;

	if (c->q<0.0)
		c->q=0.0;
	
	c->c+=eject-reuptakeandloss;
	
	if (c->c<0.0)
		c->c=0.0;

	c->w += reuptake-reprocess;
	if (c->w<0.0)
		c->w=0.0;

	hc = hdt*c->c;

	/* filter coefficients */

	c->p4 = c->p3; c->p3 = c->p2; c->p2 = c->p1; c->p1 = c->p0;
	c->q4 = c->q3; c->q3 = c->q2; c->q2 = c->q1; c->q1 = c->q0;
	c->u1 = c->u0; c->v1 = c->v0;

	c->ptr++;
	if (c->ptr==MAX_BUFFER_SIZE)
		c->ptr=0;

	if (bm<0)
		bm=0;

	bm=sqrt(bm);
	c->buffer[c->ptr]=bm;

	if(isnan(c->buffer[c->ptr]))
		printf("eroare\n");

	return (c->buffer[c->ptr]);
}

double getBufferVal(channel *c, int i)
{
	int idx;
	idx=(c->ptr+i+1);
	if (idx>=MAX_BUFFER_SIZE) idx-=MAX_BUFFER_SIZE;
	return c->buffer[idx];
}

int msToSamples(float ms)
{
	return (int)((float)SAMPLING_FREQUENCY*ms/1000.0);
}



/* compute cross-correlation and energy */

float compenergy(channel* c)
{
	int win;
	float s, e=0;

	for(win=0;win<MAX_BUFFER_SIZE-1;win++)
	{
		s=getBufferVal(c,win);
		e+=(s*s);
	}

	return (e/(float)MAX_BUFFER_SIZE);
}

   
float compute_IID(int chan )
{
	int win;
	float s_left,s_right,e_left=0, e_right=0;
	float eps=0.0000000001;

	for(win=0;win<MAX_BUFFER_SIZE-1;win++)
	{
		s_left = getBufferVal(&cochlea_left[chan],win);
		s_right = getBufferVal(&cochlea_right[chan],win); 

		e_left += (s_left*s_left); 
		e_right += (s_right*s_right);
	}
	return (20*log10(e_right/(e_left+eps)+eps));
}


int index1[MAX_CHANNEL],index2[MAX_CHANNEL];

float compute_ITD(int chan)
{
	float s,s1,s2;
	int delay,win,i,j;
	int ind,idx1,idx2;
	float w[MAX_BUFFER_SIZE];
	int n=MAX_WINDOW;
	float xval;
	float sigma;
	FILE*f;
	int d,I,I1;
	float max_value,max_ind;
	float avg1,avg2,d1,d2;

	float cov[2*MAX_DISP+1];

	for (ind=0; ind<n; ind++)
		w[ind]=0.5-0.5*cos(2.0*ind*M_PI/(n-1.0));

	for(delay=-MAX_DISP;delay<=MAX_DISP;delay++)
	{
		ind=delay+MAX_DISP;
		cov[ind]=0.0;
		idx1=MAX_BUFFER_SIZE-1-MAX_DISP;
		idx2=MAX_BUFFER_SIZE-1-MAX_DISP-delay;
		avg1=0;avg2=0;d1=0;d2=0;

		for(win=0;win<n;win++)
		{
			s1 = getBufferVal(&cochlea_left[chan],idx1-win);
			s2 = getBufferVal(&cochlea_right[chan],idx2-win);

			avg1 += s1;
			avg2 += s2;
		}

		avg1 /= n;
		avg2 /= n;

		for(win=0;win<n;win++)
		{
			s1 = getBufferVal(&cochlea_left[chan],idx1-win);
			s2 = getBufferVal(&cochlea_right[chan],idx2-win);

			d1 += (s1-avg1)*(s1-avg1);
			d2 += (s2-avg2)*(s2-avg2);
		}

		d1 = sqrt(d1);
		d2 = sqrt(d2);

		for(win=0;win<n;win++)
		{
			s1 = getBufferVal(&cochlea_left[chan],idx1-win);
			s2 = getBufferVal(&cochlea_right[chan],idx2-win);
			s = (s1-avg1)*(s2-avg2);

			cov[ind] += s;
		}
		cov[ind] /= ((float)(d1*d2));
	}

//////////////////////////
///find peak location

	max_value =- 100;
    max_ind = 0;
    for(i=index1[chan];i<=index2[chan];i++)
    {
    	if(cov[i]>max_value)
    	{
    		max_ind = i;
    		max_value = cov[i];
    	}
    }
    return (max_ind+1);
}

/*------------------------------------------------------*/
/* Main program */
/*------------------------------------------------------*/

int main (int argc, char **argv)
{
	int numChannels=128;
	int chan,tim;
	int lowerCF=80;
	int upperCF=5000;
	float sigVal_left,sigVal_right,sigVal_t,sigVal_n;
	float haircell_left,haircell_right;
	FILE *ofp_left,*ofp_right,*ofp_t,*ofp_n;
	int frame,win,delay;
	FILE*f,*itd_file,*iid_file,*ratio_file,*output_ITD,*output_IID,*output_Ratio;
	float itd;
	float xval1,xval2,ratio,iid,eps=0.000000000001;
	char s[50];
	float loc[MAX_CHANNEL],pi;


///////////////////////
///channel initialization

	initOuterMiddleEar();
	initChannels(lowerCF,upperCF,numChannels);
	initHairCells();

////////////////////////
///reference ITD for target
		
	sprintf(s,"%s/%s",TARGET_DIR,INPUTFILE_ITD);
	f = fopen(s,"r");

	for (chan=0;chan<MAX_CHANNEL;chan++)
	{
		fscanf(f,"%f ",&xval1);loc[chan]=xval1;
	}

	fclose(f);

	for (chan=0;chan<MAX_CHANNEL;chan++)
	{
		pi = 44100.0/(2.0*CCF[chan]);
		index1[chan] = (int)(loc[chan]-pi);

		if (index1[chan]<0)
			index1[chan] = 0;

		index2[chan] = (int)(loc[chan]+pi);
		
		if (index2[chan]>2*MAX_DISP)
			index2[chan]=2*MAX_DISP;
	}
////////////////////////////
//// binaural input (mixture)

	sprintf(s,"%s/%s",argv[1],INPUTFILE_LEFT);
	ofp_left = fopen(s,"r");
	printf("%s\n",s);

	sprintf(s,"%s/%s",argv[1],INPUTFILE_RIGHT);
	ofp_right = fopen(s,"r");


////////////////////////////////
///target / interference input at left ear


	sprintf(s,"%s/%s",argv[1],INPUTFILE_TARGET);
	ofp_t = fopen(s,"r");

	sprintf(s,"%s/%s",argv[1],INPUTFILE_NOISE);
	ofp_n = fopen(s,"r");



////////////////////////
///output files

	sprintf(s,"%s/%s",argv[1],INPUTFILE_ITD);
	itd_file = fopen(s,"w");

	sprintf(s,"%s/%s",argv[1],INPUTFILE_IID);
	iid_file = fopen(s,"w");

	if (COMPUTE_R)
	{
		sprintf(s,"%s/%s",argv[1],INPUTFILE_R);
		ratio_file = fopen(s,"w");
	}

//////////////////////////
//output files for training

	if (TRAIN)
	{
		sprintf(s,"%s",TRAINFILE_ITD);
		output_ITD = fopen(s,"a");

		sprintf(s,"%s",TRAINFILE_IID);
		output_IID = fopen(s,"a");

		sprintf(s,"%s",TRAINFILE_R);
		output_Ratio = fopen(s,"a");
	}
   
   
 /********************/
	tim = 0;
	frame = 0;

	while(!feof(ofp_left))
	{
		fscanf(ofp_left,"%f\n",&sigVal_left);
		fscanf(ofp_right,"%f\n",&sigVal_right);

		if(COMPUTE_R)
		{
			fscanf(ofp_t,"%f\n",&sigVal_t);
			fscanf(ofp_n,"%f\n",&sigVal_n);
		}

		for(chan=0;chan<MAX_CHANNEL;chan++)
		{
			updateCochlea(&cochlea_left[chan],sigVal_left,tim,chan);
			updateCochlea(&cochlea_right[chan],sigVal_right,tim,chan);

			if (COMPUTE_R)
			{
				updateCochlea(&cochlea_t[chan],sigVal_t,tim,chan);
				updateCochlea(&cochlea_n[chan],sigVal_n,tim,chan);
			}

			if(fmod(tim,OFFSET)==0 && tim>0)
			{
				itd = compute_ITD(chan);
				iid = compute_IID(chan);

				fprintf(itd_file,"%f ",itd);
				fprintf(iid_file,"%f ",iid);

				if (COMPUTE_R)
				{
					xval1 = compenergy(&cochlea_t[chan]);
					xval2 = compenergy(&cochlea_n[chan]);

					ratio = sqrt(xval1)/(sqrt(xval1)+sqrt(xval2)+eps);
					fprintf(ratio_file,"%f ",ratio);
				}

				if (TRAIN)
				{
					fprintf(output_ITD,"%f ",itd);
					fprintf(output_IID,"%f ",iid);
					fprintf(output_Ratio,"%f ",ratio);
				}
       //////////////////
       /////////////////
			}
		}

		tim++;

	}

	fclose(output_ITD);
	fclose(output_IID);
	fclose(output_Ratio);

	fclose(itd_file);
	fclose(iid_file);
	fclose(ratio_file);

	fclose(ofp_left);
	fclose(ofp_right);
	fclose(ofp_n);
	fclose(ofp_t);
}