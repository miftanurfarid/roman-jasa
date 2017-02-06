/* Part of speech segregation program for the 1999 IEEE Trans. Neural Net. */
/* paper by D.L. Wang and G.J. Brown.                                      */
/* This part performs resynthesis based on segregated streams (binary      */
/* masks). It takes a segregated stream (a binary mask) and the original   */
/* mixture and produces the corresponding segregated signal.               */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
//#include <ieeefp.h>
#include "common.h"


typedef struct
{
	double cf, bw, criticalRate, z, gain, expCoeff;
	double p0, p1, p2, p3, p4;
	double q0, q1, q2, q3, q4;
	double u0, u1;
	double v0, v1;
} channel;

/* function prototypes */

void help(void);
/* UNIX help */

void blip(void);
/* write dots on stderr */

void initChannels(int lowerCF, int upperCF, int numChannels);
/* initialise filterbank channels */

double updateCochlea(channel *c, float sigval, int tim);
/* process one sample of the input though the cochlea */

int msToSamples(float ms);
/* converts time in ms to samples at srate */

/* global variables */

channel cochlea[MAX_CHANNEL];
float t, dt, twoPi, twoPiT;
float  mask[MAX_FRAME][MAX_CHANNEL];

/* function declarations */

void help(void)
{
	fprintf(stderr,"-l int  lowest filter centre frequency (Hz) (500)\n");
	fprintf(stderr,"-u int  highest filter centre frequency (Hz) (2000)\n");
	fprintf(stderr,"-n int  number of channels (32)\n");
	fprintf(stderr,"-a string name of left input file\n");
	fprintf(stderr,"-b string name of right input file\n");
	fprintf(stderr,"-d float buffer decay time in ms (20.0)\n");
	fprintf(stderr,"-v bool verbose output (FALSE)\n");
}

int readSignal(float s[],char name[50])
{
	int sample = 0;
	int i;
	float x;
	FILE*f;

	f = fopen(name,"r");

	printf("%s\n",name);

	for (i=0; i<MAX_SIGNAL; i++)
		s[i] = 0.0;

	while ((fscanf(f,"%f",&x)!=EOF) && sample<MAX_SIGNAL)
	{
		s[sample]=(float)x;
		sample++;
	}

	return sample;
}

void blip(void)
{
  static int count = 0;

  fprintf(stderr,".");
  count += 1;

  if (count>32)
  	count = 0;
}

float DBtoAmplitude(float dB)
{
	return pow(10.0,(dB/20.0));
}

void initChannels(int lowerCF, int upperCF, int numChannels)
{
	float lowerERB, upperERB, spaceERB;
	channel c;
	int chan;

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
		c.bw = erb(c.cf)*BW_CORRECTION;
		c.z = exp(-twoPiT*c.bw);
		c.expCoeff = c.cf*twoPiT;
		c.gain = sqr(sqr(2*M_PI*c.bw*dt))/3.0;
		c.p0 = 0.0; c.p1 = 0.0; c.p2 = 0.0; c.p3 = 0.0; c.p4 = 0.0;
		c.q0 = 0.0; c.q1 = 0.0; c.q2 = 0.0; c.q3 = 0.0; c.q4 = 0.0;
		c.u0 = 0.0; c.u1 = 0.0;
		c.v0 = 0.0; c.v1 = 0.0;
		cochlea[chan] = c;
	}
}

double updateCochlea(channel *c, float sigval, int tim)
{
	double zz, bm;

	zz = c->z;
	c->p0 = sigval*cos(c->expCoeff*tim)+zz*(4*c->p1-zz*(6*c->p2-zz*(4*c->p3-zz*c->p4)));
	c->q0 =-sigval*sin(c->expCoeff*tim)+zz*(4*c->q1-zz*(6*c->q2-zz*(4*c->q3-zz*c->q4)));
	c->u0 = zz*(c->p1+zz*(4*c->p2+zz*c->p3));
	c->v0 = zz*(c->q1+zz*(4*c->q2+zz*c->q3));
	bm = (c->u0*cos(c->expCoeff*tim)-c->v0*sin(c->expCoeff*tim))*c->gain;

	/* filter coefficients */

	c->p4 = c->p3; c->p3 = c->p2; c->p2 = c->p1; c->p1 = c->p0;
	c->q4 = c->q3; c->q3 = c->q2; c->q2 = c->q1; c->q1 = c->q0;
	c->u1 = c->u0; c->v1 = c->v0;

	return(bm);
}

int msToSamples(float ms)
{
	return (int)((float)SAMPLING_FREQUENCY*ms/1000.0);
}

int numFrames;

void readUnaryMask(int maxFrame, int maxChan)
{
	int frame, chan;

	for (frame=0; frame<maxFrame; frame++)
	{
		for (chan=0; chan<maxChan; chan++)
			mask[frame][chan]=1;
	}
	fprintf(stderr,"ok.\n");
}

void readMask(char fname[], int maxFrame, int maxChan)
{
	int frame, chan;
	FILE *ifp;
	float x;

	fprintf(stderr,"reading mask from file %s...",fname);
	ifp = fopen(fname,"r");

	if (ifp==NULL)
	{
		fprintf(stderr,"Cannot open file %s\n",fname);
		exit(0);
	}

	fscanf(ifp,"%d\n",&maxFrame);

	numFrames = maxFrame;
	fprintf(stderr,"MAxframe=%d\n",maxFrame);

	for (frame=0; frame<maxFrame; frame++)
	{
		for (chan=0; chan<maxChan; chan++)
			fscanf(ifp,"%f ",&mask[frame][chan]);

		fscanf(ifp,"\n");
	}

	fprintf(stderr,"ok.\n");
}

void core(char fname[50],char name[50],char name_out[50], int unary)
{
	int numChannels=128;
	int lonnels=128;
	int lowerCF=80;
	int upperCF=5000;
	float* signal;
	float* resynth;
	float* nerve;
	float* nerve2;
	float* w;
	int numSamples;
	//int numFrames;
	int chan;
	int frame, i;
	float tmp;

	FILE * fout;

	signal = malloc(MAX_SIGNAL*sizeof(float));
	resynth = malloc(MAX_SIGNAL*sizeof(float));
	nerve = malloc(MAX_SIGNAL*sizeof(float));
	nerve2 = malloc(MAX_SIGNAL*sizeof(float));
	w = malloc(MAX_SIGNAL*sizeof(float));

	fout = fopen(name_out,"w");

	fprintf(stderr,"Numframes=%d\n",numFrames);
	initChannels(lowerCF,upperCF,numChannels);

	/* read the signal and mask */

	numSamples=readSignal(signal,name);

	fprintf(stderr,"read %d samples.\n",numSamples);

	numFrames=(int)(numSamples)/OFFSET;

	fprintf(stderr,"computed : %d samples.\n",numSamples);
	fprintf(stderr,"max frames=%d\n",numFrames);

	if (unary==1)
	{
		readUnaryMask(numFrames,numChannels);
	}
	else
	{
		readMask(fname,numFrames,numChannels);
	}

	for (i=0; i<numSamples; i++)
		resynth[i]=0.0;

		/* do the filtering */

	for (chan=0; chan<numChannels; chan++)
	{
		for (i=0; i<MAX_SIGNAL; i++)
			w[i]=0.0;

		for (frame=-1; frame<numFrames-1; frame++)
		{
			for (i=0; i<MAX_WINDOW/2; i++)
			{
				if ((frame*OFFSET+i)>=0)
					w[frame*OFFSET+i]=w[frame*OFFSET+i]+mask[frame+1][chan]*0.5*(1.0+cos(i*M_PI/(HWIN)+M_PI));
			}

			for (i=MAX_WINDOW/2; i<MAX_WINDOW; i++)
		      w[frame*OFFSET+i]=w[frame*OFFSET+i]+mask[frame+1][chan]*0.5*(1.0+cos((i-HWIN)*M_PI/(HWIN)));
		}

		for (i=0; i<numSamples; i++)
			nerve[numSamples-i-1]=(float)updateCochlea(&cochlea[chan],signal[i],i);

		for (i=0; i<numSamples; i++)
			nerve2[numSamples-i-1]=(float)updateCochlea(&cochlea[chan],nerve[i],i);

		for (i=0; i<numSamples; i++)
			resynth[i]+=w[i]*nerve2[i];
	}

/* write */
	for (i=0; i<numSamples; i++)
		fprintf(fout,"%f\n",resynth[i]);
}

/*------------------------------------------------------*/
/* Main program */
/*------------------------------------------------------*/

int main (int argc, char **argv)
{
	char fname[50],name[50],name_out[50];
	
	sprintf(fname,"%s/%s",TEST_DIR,INPUTFILE_MASK);
	sprintf(name,"%s/%s",TEST_DIR,INPUTFILE_MIXT);
	sprintf(name_out,"%s/%s",TEST_DIR,INPUTFILE_RESYNTH);

	core(fname,name,name_out,0);
}
