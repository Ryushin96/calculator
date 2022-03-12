#include <stdio.h>

#include <stdlib.h>

#include <math.h>

#include "complex.h"

#include "fftsg.h"

/*400G �^�C�~���O�W�b�^�̃v���O����*/

/*?????v???O???????Adevd.dat???l???????????B?t?H?[?g????????aipt.dat???f?[?^???Aaitd.dat???????????B

hpc?????A?f?q????????30mm????????????????????????*/

     complex cmplxadd(),cmplxmul(),cmplxconjg(),cmplxmai(),cmplxsub(),cmplxexp(),cmplxdiv(),cmplxmul2(),cmplxmul3();

     double cmplxreal();




    /*64416684?o?C?g?K?v*/

	/*???f?????????Z*/

	complex cmplxadd(complex c1,complex c2)

    {

    complex csum;



    csum.re=c1.re+c2.re;

    csum.im=c1.im+c2.im;



    return(csum);

    }

	/*???f?????|???Z*/

	complex cmplxmul(complex c3,complex c4)

    {

    complex cmul;



    cmul.re=(c3.re*c4.re)-(c3.im*c4.im);

    cmul.im=(c3.re*c4.im)+(c3.im*c4.re);



    return(cmul);

    }

	/*???f????????*/

	complex cmplxconjg(complex c5)

    {

    complex conjg;



    conjg.re=c5.re;

    conjg.im=(-1)*(c5.im);



    return(conjg);

    }

		/*???f?????}?C?i?X??*/

	complex cmplxmai(complex c6)

  /*  complex c6;*/

    {

    complex mainas;



    mainas.re=(-1)*(c6.re);

	mainas.im=(-1)*(c6.im);



    return(mainas);

    }

	

	/*???f??????*/

	double cmplxreal(complex c7)

  /*  complex c7;*/

    {

    double real;



    real=c7.re;



    return(real);

    }

	/*???f?????????Z*/

	complex cmplxsub(complex c8,complex c9)

  /*  complex c8,c9;*/

    {

    complex csub;



    csub.re=c8.re-c9.re;

    csub.im=c8.im-c9.im;



    return(csub);

    }

	/*e?????f???*/

	complex cmplxexp(complex c10)

  /*  complex c10;*/

    {

    complex cexp;



    cexp.re=exp(c10.re)*cos(c10.im);

    cexp.im=exp(c10.re)*sin(c10.im);



    return(cexp);

    }



	/*???f?????????Z*/

	complex cmplxdiv(complex c11,double c12)

  /*  complex c11;

	double c12;*/

    {

    complex cdiv;



    cdiv.re=c11.re/c12;

    cdiv.im=c11.im/c12;



    return(cdiv);

    }

	

		/*???f?????|???Z2*/

	complex cmplxmul2(complex c13,double c14)

  /*  complex c13;

	double c14;*/

    {

    complex cmul2;



    cmul2.re=c13.re*c14;

    cmul2.im=c13.im*c14;



    return(cmul2);

    }

	

			/*???f?????|???Z3*/

	complex cmplxmul3(double c15,complex c16)

  /*  complex c16;

	double c15;*/

    {

    complex cmul3;



    cmul3.re=c16.re*c15;

    cmul3.im=c16.im*c15;



    return(cmul3);

    }



	complex tocomplex(double c17,double c18)

{

    complex toc;

    toc.re=c17;

    toc.im=c18;

	return(toc);



}



double ran1(int *idum)

{

	int j;

	int k;

	static int iy=0;

	static int iv[32];

	double temp;



	if(*idum <=0 || !iy){

		if (-(*idum)<1) *idum=1;

		else *idum= -(*idum);

        for (j=39;j>=0;j--){

			k=(*idum)/127773;

			*idum=16807*(*idum-k*127773)-2836*k;

			if (*idum < 0) *idum += 2147483647;

			if (j < 32) iv[j] = *idum;

		}

		iy=iv[0];

	}

	k=(*idum)/127773;

	*idum=16807*(*idum-k*127773)-2836*k;

	if (*idum < 0) *idum += 2147483647;

	j=iy/(1+(2147483646)/32);

	iy=iv[j];

	iv[j] = *idum;

	if((temp=(1.0/2147483647)*iy) > (1.0-1.2e-7)) return (1.0-1.2e-7);

	else return temp;

}



double gasdev(int *idum)

/*????0?A???U1?????K???z?????????????B???l????????????ran1(idum)???g???B*/

{

	double ran1(int *idum);

	static int iset=0;

	static double gset;

	double fac,rsq,v1,v2;



	if (iset==0){

		do{

			v1=2.0*ran1(idum)-1.0;

			v2=2.0*ran1(idum)-1.0;

			rsq=v1*v1+v2*v2;

		} while (rsq >= 1.0 || rsq == 0.0);

		fac=sqrt(-2.0*log(rsq)/rsq);

		gset=v1*fac;

		iset=1;

		return v2*fac;

	} else {

		iset=0;

		return gset;

	}

}



main(void){





    complex oioo,ooio,oooo;

/*	<?p???X???????g?`=at(V/m)  ???????p???X?s?[?N?d?E?U??=at0(V/m)  ?????p???X??=t0(ps)>*/

    complex ast[2][100001],ait[2][100001],apt[2][100001],aht[2][100001];

	complex astout[200001],aitout[200001],aptout[200001],ahtout[200001];





	double ast0,ts0,ts,apt0,tp0,tp;

	int jfulls,jsmux,jts,jsch;

	int jfullp,jpmux,jtp,jpch;




/* <?p???X???X?y?N?g??=aw(V/m)>*/

    complex asw[2][20001],aiw[2][20001],apw[2][20001],ahw[2][20001];



/*<?o???p???X?p???[=aopow>*/

    double asopow[20001],aiopow[20001];

	double apopow[20001],ahopow[20001];



/*<?o???p???X????=aoarg>*/

    double asoarg[20001],aioarg[20001];

	double apoarg[20001],ahoarg[20001];



/*<?r?[???f????=s0  ???????p???X?s?[?N?p???[=pkp0  ?p???[?????W??=pow>*/

    double s0,pkps0,pkpp0,pows,powi,powp,powh,pkpin;



/*<?g??=lamdai>*/

    double lamdas[20001],lamdai[20001];

	double lamdap[20001],lamdah[20001];

	double clamdas,clamdai,clamdap,clamdah;



/*<???`??????=n>*/

    double ns[20001],ni[20001];

	double np[20001],nh[20001];

	double cns,cni,cnp,cnh;



/*<?????????`?W??=d33  |??2|=kay2>*/

    double d33,kay2;



/*<?^???????g??=k  ???g?H?`??????=b>*/

    double ks[20001],ki[20001];

	double kp[20001],kh[20001];

	double bs[20001],bi[20001];

	double bp[20001],bh[20001];

	double cks,cki,ckp,ckh,cbs,cbi,cbp,cbh;



/*<???g??=f>*/

    double fs[20001],fi[20001];

	double fp[20001],fh[20001];

	double cfs,cfi,cfp,cfh,df;



/*<?????`?W??=m>*/

    double ms[20001],mi[20001];

	double mp[20001],mpph[20001],msih[20001];

	double cms,cmi,cmp,cmpph,cmsih;



/*<?f?o?C?X??=l  ?Z?N?V??????=dsec  ?Z?N?V??????=nsec  ???]??????=nrev>*/

    double l,dsec,chirp,refl;

	int nsec,nrev;



/*<?f?o?C?X???xt0(K)>*/

    double t0;



/*<z???W????=z(i)>*/

    double z[100001];

	double z0;



/*<z???W????????=dz(i)>*/

    double dz[100001];



/*<?|???v???Q???x=vgp(m/s)>*/

    double vgp,vgh;



/* <????????=tps(ps)>*/

    double tps,tsops,tiops,tpops,thops,dt;

	int itpout,ithout;





/*<???K?????p???????????`>*/

    double rtn;

	int ix,ii,nat,ire,iim,ijit,njit,jit;

	double rm,sd,anore,anoim;

	double x[100001],jitter[1001],are[1000001],aim[1000001];

	double xxx,sumx,avex,x01,p01;



/*<FFT?p?????????`>*/

    int n;

/*	int i;*/

	int ip[65536];

    double w[65536];





	int nfft;

/*	int iopt,ier;*/

/*	double wk[131072];*/

	double niii;



	double f[262144],ff[262144];

	double fso[262144],fpo[262144],fio[262144],fho[262144];

/*	double fsoy[131072],fpoy[131072],fioy[131072],fhoy[131072];*/

/*	complex cso[131072],cpo[131072],cio[131072],cho[131072];*/





/* <do?p????>*/

    int i,j,js,ji,jp,jh,ls,li,lp,lsih,lpph;

	int iii,jjj,kkk,iw,io,jo,ippp,iipp,iiip;

	int	ko;

	int ijs,iji,ijp,ijh,soo,ioo,poo,hoo;

	int ijss,ijii,ijpp,ijhh,inoise;



/*<?`???????????p????>*/

    double zout,zl,zno,ddz,zlno,zlre,zlreno;

	complex snonli,inonli,pnonli,sihnonli,pphnonli;

	double dbs,dbi,dbp,dbsih,dbpph;

	int izno,iz;



/*<?s?[?N?p???[???p???X???????p????>*/

    int ima,imae,itmi,itma;

	double maxpwi,pwi,maxe,pwie,pwiie,dit,tipsz,pwp,maxpwp,pwpe,pwpie;



/*<?~????=pi  ?^???????????x=c0  ?^?????U?d??=ep0>*/

    double pi,c0,ep0;





    FILE *devd,*astd,*aswd,*aitd,*aiwd,*aptd,*nyuumark,*syutumark;

	FILE *apwd,*ahtd,*ahwd;

	FILE *jouhou;









	/*<?o???t?@?C?????m??>*/

	devd = fopen("devd.dat", "w");

    astd = fopen("astd.dat", "w");

	aswd = fopen("aswd.dat", "w");

	aitd = fopen("aitd.dat", "w");

	aiwd = fopen("aiwd.dat", "w");

	aptd = fopen("aptd.dat", "w");

	apwd = fopen("apwd.dat", "w");

	ahtd = fopen("ahtd.dat", "w");

	ahwd = fopen("ahwd.dat", "w");

	nyuumark = fopen("nyuumark.txt", "w");

	syutumark = fopen("syutumark.txt", "w");

	jouhou = fopen("jouhou.txt","w");









	pi=3.14159265358979;

	c0=2.998e8;

	ep0=8.8542e-12;













	oioo.re=1.0;

	oioo.im=0.0;

	ooio.re=0.0;

	ooio.im=1.0;

	oooo.re=0.0;

	oooo.im=0.0;


//���w��jouhou///////////////////////////////////////////////////////////////////////

//400Gjitter

//���̓p���[�ƃN���b�N�p���[�̐ݒ�
	/*PcL^2=pcl*/
	/*PinL^2=pinl*/
	double pcl,pinl;
	pcl=1.5;
	pinl=0.5;

//�{���w��
	int kara,made,honsuu;

	kara=1;
	made=16;
	honsuu=made-kara+1;

        ts0=0.5;            //////////////////////////////////��ts0=0.25///////////////////////////////////////

//�I�t�Z�b�g
	double toff,ptoff,tfwhm;
	ptoff = 0.1;
        tfwhm=2.08125*ts0;  /////////////////////////////////////��tfwhm=0.83333*ts0////////////////////////////////////////////
        toff=ptoff*tfwhm;

//�W�b�^���x�ݒ�
	double yokozure;
	yokozure = 0.50;     //////////////////////

//�m�C�Y���x
	double zatu;
	zatu = 0.0;


//���o��
	fprintf(jouhou,"400Gjitter\n");
	fprintf(jouhou,"PcL^2=%.3f\tPinL^2=%.3f\n",pcl,pinl);
	fprintf(jouhou,"jit\t%d~%d\t honsuu\t%d\n",kara,made,honsuu);
	fprintf(jouhou,"zatu%.4f\n",zatu);
        fprintf(jouhou,"ptoff=%.3f\n",ptoff);
        fprintf(jouhou,"tfwhm=%.4f\n",tfwhm);
	fprintf(jouhou,"yokozure%.4f\n",yokozure);
	fprintf(jouhou,"toff=%.3f\n",toff);

	fclose(jouhou);



//�쐻�덷�p������`
    for(i=0;i<=100000;i++){

	x[i]=0.0;
    }



/*<????????=ix  ????????????=n>*/

    ix=-55;

    n=100000;



/*<?????l=rm  ?W??????=sd>*/

    rm=0.0;

    sd=1.0;



	 for(i=0;i<=100000;i++){

	x[i]=gasdev(&ix);

	 }



/////////////////////////////////////////////////////////////////////////////////////////////


    for(jit=kara;jit<=made;jit++){





	for(i=0;i<=100000;i++){

	ast[0][i].re=0.0;

	ast[0][i].im=0.0;

    ait[0][i].re=0.0;

	ait[0][i].im=0.0;

    apt[0][i].re=0.0;

	apt[0][i].im=0.0;

    aht[0][i].re=0.0;

	aht[0][i].im=0.0;

	ast[1][i].re=0.0;

	ast[1][i].im=0.0;

    ait[1][i].re=0.0;

	ait[1][i].im=0.0;

    apt[1][i].re=0.0;

	apt[1][i].im=0.0;

    aht[1][i].re=0.0;

	aht[1][i].im=0.0;

    }

	

	for(i=0;i<=200000;i++){

    astout[i].re=0.0;

	astout[i].im=0.0;

    aitout[i].re=0.0;

	aitout[i].im=0.0;

    aptout[i].re=0.0;

	aptout[i].im=0.0;

    ahtout[i].re=0.0;

	ahtout[i].im=0.0;

    }

	

	for(i=0;i<=20000;i++){

    asw[0][i].re=0.0;

	aiw[0][i].re=0.0;

    apw[0][i].re=0.0;

    ahw[0][i].re=0.0;

	asw[1][i].re=0.0;

	aiw[1][i].re=0.0;

    apw[1][i].re=0.0;

    ahw[1][i].re=0.0;

	asw[0][i].im=0.0;

	aiw[0][i].im=0.0;

    apw[0][i].im=0.0;

    ahw[0][i].im=0.0;

	asw[1][i].im=0.0;

	aiw[1][i].im=0.0;

    apw[1][i].im=0.0;

    ahw[1][i].im=0.0;

    asopow[i]=0.0;

    aiopow[i]=0.0;

    apopow[i]=0.0;

    ahopow[i]=0.0;

    asoarg[i]=0.0;

    aioarg[i]=0.0;

    apoarg[i]=0.0;

    ahoarg[i]=0.0;

    lamdas[i]=0.0;

    lamdai[i]=0.0;

    lamdap[i]=0.0;

    lamdah[i]=0.0;

    ns[i]=0.0;

    ni[i]=0.0;

    np[i]=0.0;

    nh[i]=0.0;

    ks[i]=0.0;

    ki[i]=0.0;

    kp[i]=0.0;

    kh[i]=0.0;

    bs[i]=0.0;

    bi[i]=0.0;

    bp[i]=0.0;

    bh[i]=0.0;

    fs[i]=0.0;

    fi[i]=0.0;

    fp[i]=0.0;

    fh[i]=0.0;

    ms[i]=0.0;

    mi[i]=0.0;

    mp[i]=0.0;

    mpph[i]=0.0;

    msih[i]=0.0;

	}

	

    for(i=0;i<=100000;i++){

    z[i]=0.0;

    dz[i]=0.0;

//	x[i]=0.0;

	}



	for(i=0;i<=262143;i++){

    f[i]=0.0;

    ff[i]=0.0;

    fso[i]=0.0;

    fio[i]=0.0;

    fpo[i]=0.0;

    fho[i]=0.0;

	}

	for(i=0;i<=65535;i++){

	ip[i]=0;

	}

	for(i=0;i<=65535;i++){

    w[i]=0.0;

	}

	zout=0.0;

	zl=0.0;

	zno=0.0;

	ddz=0.0;

	zlno=0.0;

	zlre=0.0;

	zlreno=0.0;

	snonli.re=0.0;

	inonli.re=0.0;

	pnonli.re=0.0;

	sihnonli.re=0.0;

	pphnonli.re=0.0;

	snonli.im=0.0;

	inonli.im=0.0;

	pnonli.im=0.0;

	sihnonli.im=0.0;

	pphnonli.im=0.0;

	dbs=0.0;

	dbi=0.0;

	dbp=0.0;

	dbsih=0.0;

	dbpph=0.0;

	izno=0;

	iz=0;

	ima=0;

	imae=0;

	itmi=0;

	itma=0;

	maxpwi=0.0;

	pwi=0.0;

	maxpwp=0.0;

	pwp=0.0;

	maxe=0.0;

	pwie=0.0;

	pwiie=0.0;

	pwpe=0.0;

	pwpie=0.0;

	dit=0.0;

	tipsz=0.0;

	tps=0.0;

	tsops=0.0;

	tiops=0.0;

	tpops=0.0;

	thops=0.0;

	dt=0.0;

	itpout=0;

	ithout=0;

	z0=0.0;

    vgp=0.0;

	vgh=0.0;

	iipp=0;



/*	<PPLN???p?????[?^??????>*/

	t0=363.0;

	/*l???????@t0?????x*/

	l=12.55e-3;          ///////////////////////////////////��l=5.02e-3//////////////////////////////////////////

/*	refl=l-1.0e-4;*/

/* <?????`???[?v????????????????????=chirp>*/

    chirp=0.0e-5;

	

	fprintf(devd,"'t0(K)='\t%le\t\n",t0);

	fprintf(devd,"l(m)=\t%le\t\n",l);

    fprintf(devd,"chirp(%)=\t%le\t\n",chirp*100.0);



	clamdas=1520.0e-9;

    clamdap=1550.0e-9;

    clamdah=clamdap*clamdap/(clamdap+clamdap);

    clamdai=1.0/((1.0/clamdah)-(1.0/clamdas));

    cfs=c0/clamdas;

    cfi=c0/clamdai;

    cfp=c0/clamdap;

    cfh=c0/clamdah;

	fprintf(devd,"'clamdas(m)='\t%le\t\n",clamdas);

	fprintf(devd,"'clamdai(m)='\t%le\t\n",clamdai);

    fprintf(devd,"'clamdap(m)='\t%le\t\n",clamdap);

	fprintf(devd,"'clamdah(m)='\t%le\t\n",clamdah);

	

/*<???S???g?????????`??????>*/

	cns=sqrt(4.913+(0.1173+1.65*1.0e-8*t0*t0)/((clamdas*1.0e6)*(clamdas*1.0e6)-(0.212+2.7*1.0e-8*t0*t0)*(0.212+2.7*1.0e-8*t0*t0))-2.78*1.0e-2*(clamdas*1.0e6)*(clamdas*1.0e6));

	cni=sqrt(4.913+(0.1173+1.65*1.0e-8*t0*t0)/((clamdai*1.0e6)*(clamdai*1.0e6)-(0.212+2.7*1.0e-8*t0*t0)*(0.212+2.7*1.0e-8*t0*t0))-2.78*1.0e-2*(clamdai*1.0e6)*(clamdai*1.0e6));

	cnp=sqrt(4.913+(0.1173+1.65*1.0e-8*t0*t0)/((clamdap*1.0e6)*(clamdap*1.0e6)-(0.212+2.7*1.0e-8*t0*t0)*(0.212+2.7*1.0e-8*t0*t0))-2.78*1.0e-2*(clamdap*1.0e6)*(clamdap*1.0e6));

	cnh=sqrt(4.913+(0.1173+1.65*1.0e-8*t0*t0)/((clamdah*1.0e6)*(clamdah*1.0e6)-(0.212+2.7*1.0e-8*t0*t0)*(0.212+2.7*1.0e-8*t0*t0))-2.78*1.0e-2*(clamdah*1.0e6)*(clamdah*1.0e6));

	

	cks=2.0*pi/clamdas;

    cki=2.0*pi/clamdai;

    ckp=2.0*pi/clamdap;

    ckh=2.0*pi/clamdah;

    cbs=cns*cks;

    cbi=cni*cki;

    cbp=cnp*ckp;

    cbh=cnh*ckh;

	

    d33=(pi/2.0)*16.5e-12;

    kay2=2.0*d33;

	

    cms=cks*cks*kay2/2.0/cbs;

    cmi=cki*cki*kay2/2.0/cbi;

    cmp=ckp*ckp*kay2/2.0/cbp;

    cmpph=ckh*ckh*kay2/2.0/cbh;

    cmsih=cmpph;

	

    dsec=fabs(pi/(cbh-2.0*cbp));

    nsec=(int)(l/dsec);

    nrev=nsec/2;

    nsec=2*nrev;

    fprintf(devd,"'dsec(m)='\t%le\t\n",dsec);

	fprintf(devd,"'nsec='\t%i\t\n",nsec);

    fprintf(devd,"'drev='\t%i\t\n",nrev);

	

/* <<???K????????>>

  <???????W?????????R?q?[?????X??(=dsec)????????????=xxx??????>*/	

	xxx=0.00;



/*<???K?????????T?u???[?`?????????l=rtn(=0)???\>*/

	  ijit=400;

      njit=1000;

      

	

	/*write(6,*)rtn*/

       nat=1000000;

      ire=jit;

      iim=jit+100;

	  for(i=0;i<=1000;i++){

		  jitter[i]=gasdev(&ijit);

	  }

	  for(i=0;i<=1000000;i++){

      are[i]=gasdev(&ire);

     /* write(6,*)rtn*/

      aim[i]=gasdev(&iim);

     /* write(6,*)rtn*/

	  }

/*<???????z??????>*/

      sumx=0.0;

    for(ii=1;ii<=nsec;ii++){

    sumx=sumx+x[ii];

    }

	avex=sumx/nsec;

/*<?????O?????????z???????l=avex???\>*/

    

	/*  write(6,*)avex */

	

	for(ii=1;ii<=nsec;ii++){

        x[ii]=x[ii]-avex;

	}

	sumx=0.0;

	for(ii=1;ii<=nsec;ii++){

       sumx=sumx+x[ii];

	}

	avex=sumx/nsec;

/*<???????????????z???????l=avex???\>*/

    

	/*write(6,*)avex*/

	

	

/* <????????=dz(i)????????(MKS?P???n???g?p)>*/

    for(ii=1;ii<=nsec;ii++){

       dz[ii]=x[ii]*dsec*xxx;

	}

	

/* <<z???W????=z(i)??????(MKS?P???n???g?p)>>

     <?????????????>*/

	z[0]=0.0;

	for(ii=1;ii<=nsec;ii++){

       z[ii]=z[ii-1]+dsec+chirp*dsec*(ii-1)-chirp*dsec*(nsec/2.0);

	}

	

	fprintf(devd,"'dsecin(m)'\t%le\t\n",z[1]-z[0]);

	fprintf(devd,"'dsecave(m)'\t%le\t\n",dsec);

	fprintf(devd,"'dsecout(m)'\t%le\t\n",z[nsec]-z[nsec-1]);

	fprintf(devd,"'error='\t%le\t\n",0.0);

	

//	<duty cycle error?????????????>

/*      z[0]=0.0;

      for(ii=1;ii<=nsec;ii++){

        z[ii]=ii*dsec+dz[ii];

      }

	  fprintf(devd,"'duty cycle error='\t%le\t\n",xxx);
/*


     /* <domain length error?????????????>*/
/*
      z[0]=0.0;

      for(ii=1;ii<=nsec;ii++){

        z[ii]=z[ii-1]+dsec+dz[ii];

      }



	  fprintf(devd,"'domain length error='\t%le\t\n",xxx);
*/


  /*    s0=2.646e-12;*/

	

	

	s0=8.0e-12;



	/*PcL^2=pcl*/

	pkps0=pcl/1.0;      /////////////////////pkps0=pcl/0.25///////////////////////////////////////

	/*PinL^2=pinl*/

	pkpp0=pinl/1.0;    /////////////////////pkpp0=pinl=0.25///////////////////////////////




	fprintf(devd,"'s0(m2)='\t%le\t\n",s0);

	fprintf(devd,"'pkps0(W)='\t%le\t\n",pkps0);

	fprintf(devd,"'pkpp0(W)='\t%le\t\n",pkpp0);

	pows=ep0*c0*s0*cns/2.0;

    powi=ep0*c0*s0*cni/2.0;

    powp=ep0*c0*s0*cnp/2.0;

    powh=ep0*c0*s0*cnh/2.0;

    ast0=sqrt(pkps0/pows);

    apt0=sqrt(pkpp0/powp);



/* <<?????M?????p???X???????g?`=ast(1,Ts)(V/m)>>

   <?????M?????p???X??=ts0(ps)> */

//	ts0=0.5;    /*100G*/

/*<1(ps)????????=dt>*/

	dt=400.0;       ///////////////////////////////////////dt=800.0////////////////////////////////////////////

	fprintf(devd,"'ts0(ps)='\t%le\t\n",ts0);



	jfulls=(int)(6*ts0*dt);

/*	for(jsmux=0;jsmux<=160;jsmux++){*/

	for(jsmux=1;jsmux<=1;jsmux++){

	       x01=x[jsmux*5+100];

		   if(x01>=0.0){

		     p01=1.0;

		   }

		   if(x01<0.0){

		   /*<?????_???p?^?[???p????=p01>*/

		     /*p01=0.0*/

			      p01=1.0;

		   }

           for(jts=-1*jfulls;jts<=jfulls;jts++){

		   	jsch=(int)(5*ts0*dt*jsmux);

            ts=jts;

			ast[0][jsch+jts+jfulls].re=cmplxadd(ast[0][jsch+jts+jfulls],cmplxmul2(oioo,p01*ast0*exp(-1.0*(ts*ts)/(2.0*(ts0*dt)*(ts0*dt))))).re;

			ast[0][jsch+jts+jfulls].im=cmplxadd(ast[0][jsch+jts+jfulls],cmplxmul2(oioo,p01*ast0*exp(-1.0*(ts*ts)/(2.0*(ts0*dt)*(ts0*dt))))).im;





			}

	   	 }





/*<<?????|???v???p???X???????g?`=apt(1,Tp)(V/m)>>

      <?????|???v???p???X??=tp0(ps)>*/

	tp0=0.5;          ///////////////////////////////////tp0=0.125//////////////////////////////////////

	fprintf(devd,"'tp0(ps)='\t%le\t\n",tp0);



  /*  fclose(devd);*/



	jfullp=(int)(6*tp0*dt);

/*	for(jpmux=0;jpmux<=160;jpmux++){*/

	for(jpmux=1;jpmux<=1;jpmux++){

	  for(jtp=-1*jfullp;jtp<=jfullp;jtp++){

		jpch=(int)(yokozure*jitter[jit]*tp0*dt+(5-ptoff*1.665)*ts0*dt);

        tp=jtp;

		apt[0][jpch+jtp+jfullp].re=cmplxadd(apt[0][jpch+jtp+jfullp],cmplxmul2(oioo,apt0*exp(-1.0*(tp*tp)/(2.0*(tp0*dt)*(tp0*dt))))).re;

		apt[0][jpch+jtp+jfullp].im=cmplxadd(apt[0][jpch+jtp+jfullp],cmplxmul2(oioo,apt0*exp(-1.0*(tp*tp)/(2.0*(tp0*dt)*(tp0*dt))))).im;



	    }

	}



    for(inoise=0;inoise<=100000;inoise++){

        anore=are[inoise];

        anoim=aim[inoise];

        apt[0][inoise].re=cmplxadd(apt[0][inoise],cmplxadd(cmplxmul2(oioo,anore*apt0*zatu),cmplxmul2(ooio,anoim*apt0*zatu*0.83))).re;

		apt[0][inoise].im=cmplxadd(apt[0][inoise],cmplxadd(cmplxmul2(oioo,anore*apt0*zatu),cmplxmul2(ooio,anoim*apt0*zatu*0.83))).im;

    }





	z0=0.0;

    tps=0.0;

	for(i=0;i<=6000;i=i+1){

	fprintf(astd,"%le\t%le\t%le\t%le\t%le\n",z0,tps,cmplxmul(cmplxmul2(ast[0][i],pows),cmplxconjg(ast[0][i])).re,ast[0][i].re,ast[0][i].im);

/*	fprintf(aitd,"%le\t%le\t%le\t%le\t%le\n",z0,tps,cmplxmul(cmplxmul2(ait[0][i],powi),cmplxconjg(ait[0][i])).re,ait[0][i].re,ait[0][i].im);

	fprintf(aptd,"%le\t%le\t%le\t%le\t%le\n",z0,tps,cmplxmul(cmplxmul2(apt[0][i],powp),cmplxconjg(apt[0][i])).re,apt[0][i].re,apt[0][i].im);

	fprintf(ahtd,"%le\t%le\t%le\t%le\t%le\n",z0,tps,cmplxmul(cmplxmul2(aht[0][i],powh),cmplxconjg(aht[0][i])).re,aht[0][i].re,aht[0][i].im);

  */  tps=tps+(1.0/dt)*1.0;

	}

	/*write(201,*)

      write(301,*)

      write(401,*)

      write(501,*)*/

/*<?????M?????p???X???????g?`=ast(1,Ts)(V/m)??

>*/

	nfft=131072;

    for(i=0;i<=100000;i++){

        f[i*2]=ast[0][i].re;

		f[i*2+1]=ast[0][i].im;

    }





	for(i=0;i<=65535;i++){

	ip[i]=0;

    w[i]=0.0;

	}



  /******  call cfft(f,nfft,1,c,wk,ier)*/

   cdft(2*nfft,-1,f,ip,w);

  /******    write(6,*)ier*/



/*<?????M?????p???X???X?y?N?g???g?`=asw(1,Ws)(V/m)??????>*/

    for(i=0;i<=4096;i++){

        asw[0][i+10000].re=f[i*2]/nfft;

		asw[0][i+10000].im=f[i*2+1]/nfft;

    }

	for(j=-4096;j<=-1;j++){

        asw[0][j+10000].re=f[2*nfft+j*2]/nfft;

		asw[0][j+10000].im=f[2*nfft+j*2+1]/nfft;

    }



/*<?????|???v???p???X???????g?`=apt(1,Tp)(V/m)??FFT>*/

	for(i=0;i<=100000;i++){

        ff[i*2]=apt[0][i].re;

		ff[i*2+1]=apt[0][i].im;

    }





	for(i=0;i<=65535;i++){

	ip[i]=0;

    w[i]=0.0;

	}



	/*call cfft(f,nfft,1,c,wk,ier)*/

    cdft(2*nfft,-1,ff,ip,w);

   /*****  write(6,*)ier */



/* <?????|???v???p???X???X?y?N?g???g?`=apw(1,Wp)(V/m)??????>*/



	  for(i=0;i<=4096;i++){

        apw[0][i+10000].re=ff[i*2]/nfft;

		apw[0][i+10000].im=ff[i*2+1]/nfft;

    }

	for(j=-4096;j<=-1;j++){

        apw[0][j+10000].re=ff[2*nfft+j*2]/nfft;

		apw[0][j+10000].im=ff[2*nfft+j*2+1]/nfft;

    }

      df=1.0/(nfft*(1.0e-12/dt));

      for (i=-5000;i<= 5000;i++){

        fs[i+10000]=cfs+df*i;

        fi[i+10000]=cfi+df*i;

        fp[i+10000]=cfp+df*i;

        fh[i+10000]=cfh+df*i;

      }



      for (i=-5000;i<= 5000;i++){

        lamdas[i+10000]=c0/fs[i+10000];

        lamdai[i+10000]=c0/fi[i+10000];

        lamdap[i+10000]=c0/fp[i+10000];

        lamdah[i+10000]=c0/fh[i+10000];

      }



      z0=0.0;



 /*     for (i=-4096;i<= 4096;i++){

		fprintf(aswd,"%le\t%le\t%le\t%le\t%le\t%le\n",z0,fs[i+10000],lamdas[i+10000],cmplxmul(cmplxmul2(asw[0][i+10000],pows),cmplxconjg(asw[0][i+10000])).re,asw[0][i+10000].re,asw[0][i+10000].im);

		fprintf(aiwd,"%le\t%le\t%le\t%le\t%le\t%le\n",z0,fi[i+10000],lamdai[i+10000],cmplxmul(cmplxmul2(aiw[0][i+10000],powi),cmplxconjg(aiw[0][i+10000])).re,aiw[0][i+10000].re,aiw[0][i+10000].im);

		fprintf(apwd,"%le\t%le\t%le\t%le\t%le\t%le\n",z0,fp[i+10000],lamdap[i+10000],cmplxmul(cmplxmul2(apw[0][i+10000],powp),cmplxconjg(apw[0][i+10000])).re,apw[0][i+10000].re,apw[0][i+10000].im);

        fprintf(ahwd,"%le\t%le\t%le\t%le\t%le\t%le\n",z0,fh[i+10000],lamdah[i+10000],cmplxmul(cmplxmul2(ahw[0][i+10000],powh),cmplxconjg(ahw[0][i+10000])).re,ahw[0][i+10000].re,ahw[0][i+10000].im);

	  }*/

      /*write(202,*)

      write(302,*)

      write(402,*)

      write(502,*)*/



      for (i=-5000;i<= 5000;i++){

        ns[i+10000]=sqrt(4.913+(0.1173+1.65*1.0e-8*t0*t0)/((lamdas[i+10000]*1.0e6)*(lamdas[i+10000]*1.0e6)-(0.212+2.7*1.0e-8*t0*t0)*(0.212+2.7*1.0e-8*t0*t0))-2.78*1.0e-2*(lamdas[i+10000]*1.0e6)*(lamdas[i+10000]*1.0e6));

        ni[i+10000]=sqrt(4.913+(0.1173+1.65*1.0e-8*t0*t0)/((lamdai[i+10000]*1.0e6)*(lamdai[i+10000]*1.0e6)-(0.212+2.7*1.0e-8*t0*t0)*(0.212+2.7*1.0e-8*t0*t0))-2.78*1.0e-2*(lamdai[i+10000]*1.0e6)*(lamdai[i+10000]*1.0e6));

		np[i+10000]=sqrt(4.913+(0.1173+1.65*1.0e-8*t0*t0)/((lamdap[i+10000]*1.0e6)*(lamdap[i+10000]*1.0e6)-(0.212+2.7*1.0e-8*t0*t0)*(0.212+2.7*1.0e-8*t0*t0))-2.78*1.0e-2*(lamdap[i+10000]*1.0e6)*(lamdap[i+10000]*1.0e6));

        nh[i+10000]=sqrt(4.913+(0.1173+1.65*1.0e-8*t0*t0)/((lamdah[i+10000]*1.0e6)*(lamdah[i+10000]*1.0e6)-(0.212+2.7*1.0e-8*t0*t0)*(0.212+2.7*1.0e-8*t0*t0))-2.78*1.0e-2*(lamdah[i+10000]*1.0e6)*(lamdah[i+10000]*1.0e6));

      }



      for (i=-5000;i<= 5000;i++){

        ks[i+10000]=2.0*pi/lamdas[i+10000];

        ki[i+10000]=2.0*pi/lamdai[i+10000];

        kp[i+10000]=2.0*pi/lamdap[i+10000];

        kh[i+10000]=2.0*pi/lamdah[i+10000];

      }



      for (i=-5000;i<= 5000;i++){

        bs[i+10000]=ns[i+10000]*ks[i+10000];

        bi[i+10000]=ni[i+10000]*ki[i+10000];

        bp[i+10000]=np[i+10000]*kp[i+10000];

        bh[i+10000]=nh[i+10000]*kh[i+10000];

      }



      for (i=-5000;i<= 5000;i++){

        ms[i+10000]=ks[i+10000]*ks[i+10000]*kay2/2.0/bs[i+10000];

        mi[i+10000]=ki[i+10000]*ki[i+10000]*kay2/2.0/bi[i+10000];

        mp[i+10000]=kp[i+10000]*kp[i+10000]*kay2/2.0/bp[i+10000];

        mpph[i+10000]=kh[i+10000]*kh[i+10000]*kay2/2.0/bh[i+10000];

        msih[i+10000]=kh[i+10000]*kh[i+10000]*kay2/2.0/bh[i+10000];

      }



      vgp=2.0*pi*((fp[10001]-fp[10000])/(bp[10001]-bp[10000]));

      vgh=2.0*pi*((fh[10001]-fh[10000])/(bh[10001]-bh[10000]));



      /*zout=2.0e-4;*//* 1.0e-3  */

      zl=0.0;

      zno=3.125e4;/*�v�Z�̕�����*/          /////////////////////////////zno=1.25e4//////////////////////////////

      izno=(int)(zno);

      ddz=l/zno;

	  zout=ddz;

      iz=0;

      zlno=z[iz];

      zlre=z[iz+1];

      zlreno=z[iz+2];

      snonli.re=0.0,snonli.im=0.0;

      inonli.re=0.0,inonli.im=0.0;

      pnonli.re=0.0,pnonli.im=0.0;

      sihnonli.re=0.0,sihnonli.im=0.0;

      pphnonli.re=0.0,pphnonli.im=0.0;



      for (iipp=1;iipp<=25488;iipp++){





        if(zlno<=zl && zl<zlre){

          for (js= -256;js<=256;js++){

            for (ls= -1024;ls<=1024;ls++){

              dbs=bh[js+ls+10000]-bs[js+10000]-bi[ls+10000];

              if(dbs<0.0||dbs>0.0){

                snonli.re=cmplxadd(snonli,cmplxmul3((ms[js+10000]/dbs),cmplxmul(ahw[0][js+ls+10000],cmplxmul(cmplxconjg(aiw[0][ls+10000]),cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbs),zl)),cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbs),(zl+ddz)))))))).re;

				snonli.im=cmplxadd(snonli,cmplxmul3((ms[js+10000]/dbs),cmplxmul(ahw[0][js+ls+10000],cmplxmul(cmplxconjg(aiw[0][ls+10000]),cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbs),zl)),cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbs),(zl+ddz)))))))).im;

              }

            }

            asw[1][js+10000].re=cmplxadd(asw[0][js+10000],snonli).re;

            asw[1][js+10000].im=cmplxadd(asw[0][js+10000],snonli).im;

            snonli.re=0.0;

			snonli.im=0.0;

          }



          for(ji= -256;ji<= 256;ji++){

            for(li= -1024;li<=1024;li++){

              dbi=bh[ji+li+10000]-bs[li+10000]-bi[ji+10000];

              if(dbi<0.0||dbi>0.0){

				inonli.re=cmplxadd(inonli,cmplxmul3((mi[ji+10000]/dbi),cmplxmul(ahw[0][ji+li+10000],cmplxmul(cmplxconjg(asw[0][li+10000]),cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbi),zl)),cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbi),(zl+ddz)))))))).re;

				inonli.im=cmplxadd(inonli,cmplxmul3((mi[ji+10000]/dbi),cmplxmul(ahw[0][ji+li+10000],cmplxmul(cmplxconjg(asw[0][li+10000]),cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbi),zl)),cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbi),(zl+ddz)))))))).im;

              }

            }

            aiw[1][ji+10000].re=cmplxadd(aiw[0][ji+10000],inonli).re;

			aiw[1][ji+10000].im=cmplxadd(aiw[0][ji+10000],inonli).im;

            inonli.re=0.0;

			inonli.im=0.0;

          }



          for(jp= -256;jp<= 256;jp++){

            for(lp= -1024;lp<=1024;lp++){

              dbp=bh[jp+lp+10000]-bp[lp+10000]-bp[jp+10000];

              if(dbp<0.0||dbp>0.0){

                pnonli.re=cmplxadd(pnonli,cmplxmul3((mp[jp+10000]/dbp),cmplxmul(ahw[0][jp+lp+10000],cmplxmul(cmplxconjg(apw[0][lp+10000]),cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbp),zl)),cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbp),(zl+ddz)))))))).re;

				pnonli.im=cmplxadd(pnonli,cmplxmul3((mp[jp+10000]/dbp),cmplxmul(ahw[0][jp+lp+10000],cmplxmul(cmplxconjg(apw[0][lp+10000]),cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbp),zl)),cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbp),(zl+ddz)))))))).im;

              }

            }

            apw[1][jp+10000].re=cmplxadd(apw[0][jp+10000],pnonli).re;

			apw[1][jp+10000].im=cmplxadd(apw[0][jp+10000],pnonli).im;

            pnonli.re=0.0;

	        pnonli.im=0.0;

          }



          for(jh= -512;jh<=512;jh++){

            for(lsih= -1024;lsih<=1024;lsih++){

              dbsih=bh[jh+10000]-bs[jh+lsih+10000]-bi[-1*lsih+10000];

              if(dbsih<0.0||dbsih>0.0){

			    sihnonli.re=cmplxadd(sihnonli,cmplxmul3((msih[jh+10000]/dbsih),cmplxmul(asw[0][jh+lsih+10000],cmplxmul(aiw[0][-1*lsih+10000],cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(ooio,dbsih),zl)),cmplxexp(cmplxmul2(cmplxmul2(ooio,dbsih),(zl+ddz)))))))).re;

				sihnonli.im=cmplxadd(sihnonli,cmplxmul3((msih[jh+10000]/dbsih),cmplxmul(asw[0][jh+lsih+10000],cmplxmul(aiw[0][-1*lsih+10000],cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(ooio,dbsih),zl)),cmplxexp(cmplxmul2(cmplxmul2(ooio,dbsih),(zl+ddz)))))))).im;

              }

            }

            for(lpph= -1024;lpph<= 1024;lpph++){

              dbpph=bh[jh+10000]-bp[jh+lpph+10000]-bp[-1*lpph+10000];

              if(dbpph<0.0||dbpph>0.0){

			    pphnonli.re=cmplxadd(pphnonli,cmplxmul3((mpph[jh+10000]/dbpph),cmplxmul(apw[0][jh+lpph+10000],cmplxmul(apw[0][-1*lpph+10000],cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(ooio,dbpph),zl)),cmplxexp(cmplxmul2(cmplxmul2(ooio,dbpph),(zl+ddz)))))))).re;

				pphnonli.im=cmplxadd(pphnonli,cmplxmul3((mpph[jh+10000]/dbpph),cmplxmul(apw[0][jh+lpph+10000],cmplxmul(apw[0][-1*lpph+10000],cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(ooio,dbpph),zl)),cmplxexp(cmplxmul2(cmplxmul2(ooio,dbpph),(zl+ddz)))))))).im;

              }

            }

            ahw[1][jh+10000].re=cmplxsub(ahw[0][jh+10000],cmplxadd(sihnonli,cmplxdiv(pphnonli,2.0))).re;

			ahw[1][jh+10000].im=cmplxsub(ahw[0][jh+10000],cmplxadd(sihnonli,cmplxdiv(pphnonli,2.0))).im;

            sihnonli.re=0.0;

	        sihnonli.im=0.0;

            pphnonli.re=0.0;

			pphnonli.im=0.0;

          }



          for(iii= -10000;iii<=10000;iii++){

            asw[0][iii+10000].re=0.0;

            asw[0][iii+10000].im=0.0;

            aiw[0][iii+10000].re=0.0;

            aiw[0][iii+10000].im=0.0;

			apw[0][iii+10000].re=0.0;

            apw[0][iii+10000].im=0.0;

			ahw[0][iii+10000].re=0.0;

            ahw[0][iii+10000].im=0.0;

          }

          for(jjj= -10000;jjj<=10000;jjj++){

            asw[0][jjj+10000].re=asw[1][jjj+10000].re;

			asw[0][jjj+10000].im=asw[1][jjj+10000].im;

            aiw[0][jjj+10000].re=aiw[1][jjj+10000].re;

			aiw[0][jjj+10000].im=aiw[1][jjj+10000].im;

            apw[0][jjj+10000].re=apw[1][jjj+10000].re;

			apw[0][jjj+10000].im=apw[1][jjj+10000].im;

            ahw[0][jjj+10000].re=ahw[1][jjj+10000].re;

			ahw[0][jjj+10000].im=ahw[1][jjj+10000].im;

          }

          for(kkk= -10000;kkk<=10000;kkk++){

            asw[1][kkk+10000].re=0.0;

            asw[1][kkk+10000].im=0.0;

            aiw[1][kkk+10000].re=0.0;

            aiw[1][kkk+10000].im=0.0;

			apw[1][kkk+10000].re=0.0;

            apw[1][kkk+10000].im=0.0;

			ahw[1][kkk+10000].re=0.0;

            ahw[1][kkk+10000].im=0.0;

          }

          zl=zl+ddz;

		  }





        if(zlre<=zl && zl<zlreno){

          for(js= -256;js<=256;js++){

            for(ls= -1024;ls<=1024;ls++){

              dbs=bh[js+ls+10000]-bs[js+10000]-bi[ls+10000];

              if(dbs<0.0||dbs>0.0){

                snonli.re=cmplxadd(snonli,cmplxmul3((ms[js+10000]/dbs),cmplxmul(ahw[0][js+ls+10000],cmplxmul(cmplxconjg(aiw[0][ls+10000]),cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbs),zl)),cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbs),(zl+ddz)))))))).re;

				snonli.im=cmplxadd(snonli,cmplxmul3((ms[js+10000]/dbs),cmplxmul(ahw[0][js+ls+10000],cmplxmul(cmplxconjg(aiw[0][ls+10000]),cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbs),zl)),cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbs),(zl+ddz)))))))).im;

			  }

            }

            asw[1][js+10000].re=cmplxsub(asw[0][js+10000],snonli).re;

			asw[1][js+10000].im=cmplxsub(asw[0][js+10000],snonli).im;

            snonli.re=0.0;

            snonli.im=0.0;

          }



          for(ji=-256;ji<=256;ji++){

            for(li=-1024;li<=1024;li++){

              dbi=bh[ji+li+10000]-bs[li+10000]-bi[ji+10000];

              if(dbi<0.0||dbi>0.0){

				inonli.re=cmplxadd(inonli,cmplxmul3((mi[ji+10000]/dbi),cmplxmul(ahw[0][ji+li+10000],cmplxmul(cmplxconjg(asw[0][li+10000]),cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbi),zl)),cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbi),(zl+ddz)))))))).re;

				inonli.im=cmplxadd(inonli,cmplxmul3((mi[ji+10000]/dbi),cmplxmul(ahw[0][ji+li+10000],cmplxmul(cmplxconjg(asw[0][li+10000]),cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbi),zl)),cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbi),(zl+ddz)))))))).im;

              }

            }

	        aiw[1][ji+10000].re=cmplxsub(aiw[0][ji+10000],inonli).re;

			aiw[1][ji+10000].im=cmplxsub(aiw[0][ji+10000],inonli).im;

            inonli.re=0.0;

            inonli.im=0.0;

          }

          for(jp= -256;jp<=256;jp++){

            for(lp= -1024;lp<=1024;lp++){

              dbp=bh[jp+lp+10000]-bp[lp+10000]-bp[jp+10000];

              if(dbp<0.0||dbp>0.0){

			    pnonli.re=cmplxadd(pnonli,cmplxmul3((mp[jp+10000]/dbp),cmplxmul(ahw[0][jp+lp+10000],cmplxmul(cmplxconjg(apw[0][lp+10000]),cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbp),zl)),cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbp),(zl+ddz)))))))).re;

				pnonli.im=cmplxadd(pnonli,cmplxmul3((mp[jp+10000]/dbp),cmplxmul(ahw[0][jp+lp+10000],cmplxmul(cmplxconjg(apw[0][lp+10000]),cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbp),zl)),cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),dbp),(zl+ddz)))))))).im;

              }

            }

	        apw[1][jp+10000].re=cmplxsub(apw[0][jp+10000],pnonli).re;

			apw[1][jp+10000].im=cmplxsub(apw[0][jp+10000],pnonli).im;

            pnonli.re=0.0;

            pnonli.im=0.0;

          }



          for(jh= -512;jh<=512;jh++){

            for(lsih= -1024;lsih<=1024;lsih++){

              dbsih=bh[jh+10000]-bs[jh+lsih+10000]-bi[-1*lsih+10000];

              if(dbsih<0.0||dbsih>0.0){

                sihnonli.re=cmplxadd(sihnonli,cmplxmul3((msih[jh+10000]/dbsih),cmplxmul(asw[0][jh+lsih+10000],cmplxmul(aiw[0][-1*lsih+10000],cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(ooio,dbsih),zl)),cmplxexp(cmplxmul2(cmplxmul2(ooio,dbsih),(zl+ddz)))))))).re;

				sihnonli.im=cmplxadd(sihnonli,cmplxmul3((msih[jh+10000]/dbsih),cmplxmul(asw[0][jh+lsih+10000],cmplxmul(aiw[0][-1*lsih+10000],cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(ooio,dbsih),zl)),cmplxexp(cmplxmul2(cmplxmul2(ooio,dbsih),(zl+ddz)))))))).im;

              }

            }

            for(lpph= -1024;lpph<=1024;lpph++){

              dbpph=bh[jh+10000]-bp[jh+lpph+10000]-bp[-1*lpph+10000];

              if(dbpph<0.0||dbpph>0.0){

				pphnonli.re=cmplxadd(pphnonli,cmplxmul3((mpph[jh+10000]/dbpph),cmplxmul(apw[0][jh+lpph+10000],cmplxmul(apw[0][-1*lpph+10000],cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(ooio,dbpph),zl)),cmplxexp(cmplxmul2(cmplxmul2(ooio,dbpph),(zl+ddz)))))))).re;

				pphnonli.im=cmplxadd(pphnonli,cmplxmul3((mpph[jh+10000]/dbpph),cmplxmul(apw[0][jh+lpph+10000],cmplxmul(apw[0][-1*lpph+10000],cmplxsub(cmplxexp(cmplxmul2(cmplxmul2(ooio,dbpph),zl)),cmplxexp(cmplxmul2(cmplxmul2(ooio,dbpph),(zl+ddz)))))))).im;

              }

            }

			niii=2.0;

	        ahw[1][jh+10000].re=cmplxadd(ahw[0][jh+10000],cmplxadd(sihnonli,cmplxdiv(pphnonli,niii))).re;

            ahw[1][jh+10000].im=cmplxadd(ahw[0][jh+10000],cmplxadd(sihnonli,cmplxdiv(pphnonli,niii))).im;

            sihnonli.re=0.0;

            sihnonli.im=0.0;

            pphnonli.re=0.0;

            pphnonli.im=0.0;

          }



          for(iii= -10000;iii<=10000;iii++){

            asw[0][iii+10000].re=0.0;

            asw[0][iii+10000].im=0.0;

			aiw[0][iii+10000].re=0.0;

            aiw[0][iii+10000].im=0.0;

			apw[0][iii+10000].re=0.0;

            apw[0][iii+10000].im=0.0;

			ahw[0][iii+10000].re=0.0;

            ahw[0][iii+10000].im=0.0;

          }

          for(jjj= -10000;jjj<=10000;jjj++){

			asw[0][jjj+10000].re=asw[1][jjj+10000].re;

            asw[0][jjj+10000].im=asw[1][jjj+10000].im;

			aiw[0][jjj+10000].re=aiw[1][jjj+10000].re;

            aiw[0][jjj+10000].im=aiw[1][jjj+10000].im;

			apw[0][jjj+10000].re=apw[1][jjj+10000].re;

            apw[0][jjj+10000].im=apw[1][jjj+10000].im;

			ahw[0][jjj+10000].re=ahw[1][jjj+10000].re;

            ahw[0][jjj+10000].im=ahw[1][jjj+10000].im;

          }



          for(kkk= -10000;kkk<=10000;kkk++){

            asw[1][kkk+10000].re=0.0;

            asw[1][kkk+10000].im=0.0;

			aiw[1][kkk+10000].re=0.0;

            aiw[1][kkk+10000].im=0.0;

			apw[1][kkk+10000].re=0.0;

            apw[1][kkk+10000].im=0.0;

			ahw[1][kkk+10000].re=0.0;

            ahw[1][kkk+10000].im=0.0;

          }

          zl=zl+ddz;

        }



        if(zlreno<=zl){

          iz=iz+2;

          zlno=z[iz];

          zlre=z[iz+1];

          zlreno=z[iz+2];

        }



        if(zl<=4.904e-7){

          for(iw= -4096;iw<=4096;iw++){

  /*          fprintf(aswd,"%le\t%le\t%le\t%le\t%le\t%le\n",zout*(1.0e3),fs[iw+10000],lamdas[iw+10000],pows*(cmplxmul(asw[0][iw+10000],cmplxconjg(asw[0][iw+10000]))).re,asw[0][iw+10000].re,asw[0][iw+10000].im);

            fprintf(aiwd,"%le\t%le\t%le\t%le\t%le\t%le\n",zout*(1.0e3),fi[iw+10000],lamdai[iw+10000],powi*(cmplxmul(aiw[0][iw+10000],cmplxconjg(aiw[0][iw+10000]))).re,aiw[0][iw+10000].re,aiw[0][iw+10000].im);

*/            fprintf(apwd,"%le\t%le\t%le\t%le\t%le\t%le\n",zout*(1.0e3),fp[iw+10000],lamdap[iw+10000],powp*(cmplxmul(apw[0][iw+10000],cmplxconjg(apw[0][iw+10000]))).re,apw[0][iw+10000].re,apw[0][iw+10000].im);

/*            fprintf(ahwd,"%le\t%le\t%le\t%le\t%le\t%le\n",zout*(1.0e3),fh[iw+10000],lamdah[iw+10000],powh*(cmplxmul(ahw[0][iw+10000],cmplxconjg(ahw[0][iw+10000]))).re,ahw[0][iw+10000].re,ahw[0][iw+10000].im);

*/		  }



  /*        write(202,*)

          write(302,*)

          write(402,*)

          write(502,*)  */







		  for(i=0;i<=262143;i++){

              fso[i]=0.0;

              fio[i]=0.0;

              fpo[i]=0.0;

              fho[i]=0.0;

		  }





          for(io=0;io<=4096;io++){

    /*        fso[io*2]=asw[0][io+10000].re;

			fso[io*2+1]=asw[0][io+10000].im;

            fio[io*2]=aiw[0][io+10000].re;

            fio[io*2+1]=aiw[0][io+10000].im;

			fpo[io*2]=apw[0][io+10000].re;

			fpo[io*2+1]=apw[0][io+10000].im;

			fho[io*2]=ahw[0][io+10000].re;

		    fho[io*2+1]=ahw[0][io+10000].im;*/

            fso[io*2]=cmplxmul(asw[0][io+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bs[io+10000]),zout)),(double)1)).re;

			fso[io*2+1]=cmplxmul(asw[0][io+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bs[io+10000]),zout)),(double)1)).im;

            fio[io*2]=cmplxmul(aiw[0][io+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bi[io+10000]),zout)),(double)1)).re;

			fio[io*2+1]=cmplxmul(aiw[0][io+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bi[io+10000]),zout)),(double)1)).im;

			fpo[io*2]=cmplxmul(apw[0][io+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bp[io+10000]),zout)),(double)1)).re;

			fpo[io*2+1]=cmplxmul(apw[0][io+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bp[io+10000]),zout)),(double)1)).im;

			fho[io*2]=cmplxmul(ahw[0][io+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bh[io+10000]),zout)),(double)1)).re;

			fho[io*2+1]=cmplxmul(ahw[0][io+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bh[io+10000]),zout)),(double)1)).im;

          }



          for(jo=-4096;jo<=-1;jo++){

/*			fso[2*nfft+jo*2]=asw[0][jo+10000].re;

			fso[2*nfft+jo*2+1]=asw[0][jo+10000].im;

            fio[2*nfft+jo*2]=aiw[0][jo+10000].re;

            fio[2*nfft+jo*2+1]=aiw[0][jo+10000].im;

			fpo[2*nfft+jo*2]=apw[0][jo+10000].re;

			fpo[2*nfft+jo*2+1]=apw[0][jo+10000].im;

			fho[2*nfft+jo*2]=ahw[0][jo+10000].re;

		    fho[2*nfft+jo*2+1]=ahw[0][jo+10000].im;*/

            fso[2*nfft+jo*2]=cmplxmul(asw[0][jo+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bs[jo+10000]),zout)),(double)1)).re;

			fso[2*nfft+jo*2+1]=cmplxmul(asw[0][jo+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bs[jo+10000]),zout)),(double)1)).im;

			fio[2*nfft+jo*2]=cmplxmul(aiw[0][jo+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bi[jo+10000]),zout)),(double)1)).re;

			fio[2*nfft+jo*2+1]=cmplxmul(aiw[0][jo+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bi[jo+10000]),zout)),(double)1)).im;

			fpo[2*nfft+jo*2]=cmplxmul(apw[0][jo+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bp[jo+10000]),zout)),(double)1)).re;

			fpo[2*nfft+jo*2+1]=cmplxmul(apw[0][jo+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bp[jo+10000]),zout)),(double)1)).im;

			fho[2*nfft+jo*2]=cmplxmul(ahw[0][jo+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bh[jo+10000]),zout)),(double)1)).re;

			fho[2*nfft+jo*2+1]=cmplxmul(ahw[0][jo+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bh[jo+10000]),zout)),(double)1)).im;

          }









          itpout=(int)((zout/vgp)*1.0e12*dt);

          ithout=(int)((zout/vgh)*1.0e12*dt);







	      for(i=0;i<=65535;i++){

		  ip[i]=0;

          w[i]=0.0;

		  }

			  cdft(2*nfft,1,fso,ip,w);



		  for(ijs=0;ijs<=100000;ijs++){

            astout[ijs-itpout+100000].re=fso[ijs*2];

            astout[ijs-itpout+100000].im=fso[ijs*2+1];

          }

          for(ijss=0;ijss<=100000;ijss++){

            ast[1][ijss].re=astout[ijss+100000].re;

			ast[1][ijss].im=astout[ijss+100000].im;

          }



		   cdft(2*nfft,1,fio,ip,w);





		  for(iji=0;iji<=100000;iji++){

            aitout[iji-itpout+100000].re=fio[iji*2];

            aitout[iji-itpout+100000].im=fio[iji*2+1];

          }

          for(ijii=0;ijii<=100000;ijii++){

            ait[1][ijii].re=aitout[ijii+100000].re;

			ait[1][ijii].im=aitout[ijii+100000].im;

          }



           cdft(2*nfft,1,fpo,ip,w);



		  for(ijp=0;ijp<=100000;ijp++){

            aptout[ijp-itpout+100000].re=fpo[ijp*2];

            aptout[ijp-itpout+100000].im=fpo[ijp*2+1];

          }

          for(ijpp=0;ijpp<=100000;ijpp++){

            apt[1][ijpp].re=aptout[ijpp+100000].re;

			apt[1][ijpp].im=aptout[ijpp+100000].im;

          }





          cdft(2*nfft,1,fho,ip,w);



          for(ijh=0;ijh<=100000;ijh++){

            ahtout[ijh-ithout+100000].re=fho[ijh*2];

            ahtout[ijh-ithout+100000].im=fho[ijh*2+1];

          }

          for(ijhh=0;ijhh<=100000;ijhh++){

            aht[1][ijhh].re=ahtout[ijhh+100000].re;

			aht[1][ijhh].im=ahtout[ijhh+100000].im;

          }



/*		  printf("%le  ",zl);*/



		  tsops=0.0;
if(0.0<zout&&zout<1.0e-3){         /////////////////////ifadd commentin//////////////////////

          for(soo=0;soo<=6000;soo=soo+1){

            fprintf(astd,"%le\t%le\t%le\t%le\t%le\t%le\n",zout*(1.0e3),tsops,cmplxmul3(pows,cmplxmul(ast[1][soo],cmplxconjg(ast[1][soo]))).re,ast[1][soo].re,ast[1][soo].im,zl);

			tsops=tsops+(1.0/dt)*1.0;

        }


}

          maxpwi=0.0;
	  maxpwp=0.0;

          for(ima=0;ima<=10000;ima++){

            pwi=cmplxmul3(powi,cmplxmul(ait[1][ima],cmplxconjg(ait[1][ima]))).re;
	    pwp=cmplxmul3(powp,cmplxmul(apt[1][ima],cmplxconjg(apt[1][ima]))).re;

	        if(maxpwi<pwi){
              maxpwi=pwi;
            }

	        if(maxpwp<pwp){
              maxpwp=pwp;
            }

          }



          maxe=maxpwi*exp(-1.0);

		  for(imae=0;imae<=100000;imae++){

            pwie=cmplxmul3(powi,cmplxmul(ait[1][imae],cmplxconjg(ait[1][imae]))).re;

            pwiie=cmplxmul3(powi,cmplxmul(ait[1][imae+1],cmplxconjg(ait[1][imae+1]))).re;

            if(pwie<=maxe&&maxe<=pwiie){

              itmi=imae;

            }

            if(pwie>=maxe&&maxe>=pwiie){

              itma=imae;

            }

          }

          dit=itma-itmi;



          tipsz=dit/dt;







          tiops=0.0;

	if(9.5e-3<zout&&zout<10.5e-3){       //////////////ifnakami///////////////
         for(ioo=0;ioo<=6000;ioo=ioo+1){

			fprintf(aitd,"%le\t%le\t%le\t%le\t%le\t%le\t%le\n",zout*(1.0e3),tiops,cmplxmul3(powi,cmplxmul(ait[1][ioo],cmplxconjg(ait[1][ioo]))).re,ait[1][ioo].re,ait[1][ioo].im,maxpwi,tipsz);
			fprintf(syutumark,"%le\t%le\t%le\n",tiops,cmplxmul3(powi,cmplxmul(ait[1][ioo],cmplxconjg(ait[1][ioo]))).re,maxpwi);

            tiops=tiops+(1.0/dt)*1.0;

          }
	}


          tpops=0.0;
if(0.0<zout&&zout<1.0e-3){
          for(poo=0;poo<=6000;poo=poo+1){

			fprintf(aptd,"%le\t%le\t%le\t%le\t%le\n",zout*(1.0e3),tpops,cmplxmul3(powp,cmplxmul(apt[1][poo],cmplxconjg(apt[1][poo]))).re,apt[1][poo].re,apt[1][poo].im);
			fprintf(nyuumark,"%le\t%le\t%le\n",tpops,cmplxmul3(powp,cmplxmul(apt[1][poo],cmplxconjg(apt[1][poo]))).re,maxpwp);

            tpops=tpops+(1.0/dt)*1.0;

		  }
}




          thops=0.0;

     /*     if(refl<=zl){*/

 /*         for(hoo=0;hoo<=100000;hoo=hoo+1){

			fprintf(ahtd,"%le\t%le\t%le\t%le\t%le\n",zout*(1.0e3),thops,cmplxmul3(powh,cmplxmul(aht[1][hoo],cmplxconjg(aht[1][hoo]))).re,aht[1][hoo].re,aht[1][hoo].im);

            thops=thops+(1.0/dt)*1.0;

          }*/

		/*  }*/



          zout=1.0e-3;/*1.0e-3*/



        }





        if(zout<=zl){

          for(iw= -4096;iw<=4096;iw++){

  /*          fprintf(aswd,"%le\t%le\t%le\t%le\t%le\t%le\n",zout*(1.0e3),fs[iw+10000],lamdas[iw+10000],pows*(cmplxmul(asw[0][iw+10000],cmplxconjg(asw[0][iw+10000]))).re,asw[0][iw+10000].re,asw[0][iw+10000].im);

            fprintf(aiwd,"%le\t%le\t%le\t%le\t%le\t%le\n",zout*(1.0e3),fi[iw+10000],lamdai[iw+10000],powi*(cmplxmul(aiw[0][iw+10000],cmplxconjg(aiw[0][iw+10000]))).re,aiw[0][iw+10000].re,aiw[0][iw+10000].im);
*/
/*            fprintf(apwd,"%le\t%le\t%le\t%le\t%le\t%le\n",zout*(1.0e3),fp[iw+10000],lamdap[iw+10000],powp*(cmplxmul(apw[0][iw+10000],cmplxconjg(apw[0][iw+10000]))).re,apw[0][iw+10000].re,apw[0][iw+10000].im);

          fprintf(ahwd,"%le\t%le\t%le\t%le\t%le\t%le\n",zout*(1.0e3),fh[iw+10000],lamdah[iw+10000],powh*(cmplxmul(ahw[0][iw+10000],cmplxconjg(ahw[0][iw+10000]))).re,ahw[0][iw+10000].re,ahw[0][iw+10000].im);

*/		  }



  /*        write(202,*)

          write(302,*)

          write(402,*)

          write(502,*)  */







		  for(i=0;i<=262143;i++){

              fso[i]=0.0;

              fio[i]=0.0;

              fpo[i]=0.0;

              fho[i]=0.0;

		  }



          for(io=0;io<=4096;io++){

    /*        fso[io*2]=asw[0][io+10000].re;

			fso[io*2+1]=asw[0][io+10000].im;

            fio[io*2]=aiw[0][io+10000].re;

            fio[io*2+1]=aiw[0][io+10000].im;

			fpo[io*2]=apw[0][io+10000].re;

			fpo[io*2+1]=apw[0][io+10000].im;

			fho[io*2]=ahw[0][io+10000].re;

		    fho[io*2+1]=ahw[0][io+10000].im;*/

            fso[io*2]=cmplxmul(asw[0][io+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bs[io+10000]),zout)),(double)1)).re;

			fso[io*2+1]=cmplxmul(asw[0][io+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bs[io+10000]),zout)),(double)1)).im;

            fio[io*2]=cmplxmul(aiw[0][io+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bi[io+10000]),zout)),(double)1)).re;

			fio[io*2+1]=cmplxmul(aiw[0][io+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bi[io+10000]),zout)),(double)1)).im;

			fpo[io*2]=cmplxmul(apw[0][io+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bp[io+10000]),zout)),(double)1)).re;

			fpo[io*2+1]=cmplxmul(apw[0][io+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bp[io+10000]),zout)),(double)1)).im;

			fho[io*2]=cmplxmul(ahw[0][io+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bh[io+10000]),zout)),(double)1)).re;

			fho[io*2+1]=cmplxmul(ahw[0][io+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bh[io+10000]),zout)),(double)1)).im;

          }



          for(jo=-4096;jo<=-1;jo++){

/*			fso[2*nfft+jo*2]=asw[0][jo+10000].re;

			fso[2*nfft+jo*2+1]=asw[0][jo+10000].im;

            fio[2*nfft+jo*2]=aiw[0][jo+10000].re;

            fio[2*nfft+jo*2+1]=aiw[0][jo+10000].im;

			fpo[2*nfft+jo*2]=apw[0][jo+10000].re;

			fpo[2*nfft+jo*2+1]=apw[0][jo+10000].im;

			fho[2*nfft+jo*2]=ahw[0][jo+10000].re;

		    fho[2*nfft+jo*2+1]=ahw[0][jo+10000].im;*/

            fso[2*nfft+jo*2]=cmplxmul(asw[0][jo+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bs[jo+10000]),zout)),(double)1)).re;

			fso[2*nfft+jo*2+1]=cmplxmul(asw[0][jo+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bs[jo+10000]),zout)),(double)1)).im;

			fio[2*nfft+jo*2]=cmplxmul(aiw[0][jo+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bi[jo+10000]),zout)),(double)1)).re;

			fio[2*nfft+jo*2+1]=cmplxmul(aiw[0][jo+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bi[jo+10000]),zout)),(double)1)).im;

			fpo[2*nfft+jo*2]=cmplxmul(apw[0][jo+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bp[jo+10000]),zout)),(double)1)).re;

			fpo[2*nfft+jo*2+1]=cmplxmul(apw[0][jo+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bp[jo+10000]),zout)),(double)1)).im;

			fho[2*nfft+jo*2]=cmplxmul(ahw[0][jo+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bh[jo+10000]),zout)),(double)1)).re;

			fho[2*nfft+jo*2+1]=cmplxmul(ahw[0][jo+10000],cmplxmul2(cmplxexp(cmplxmul2(cmplxmul2(cmplxmai(ooio),bh[jo+10000]),zout)),(double)1)).im;

          }



          itpout=(int)((zout/vgp)*1.0e12*dt);

          ithout=(int)((zout/vgh)*1.0e12*dt);



	      for(i=0;i<=65535;i++){

		  ip[i]=0;

          w[i]=0.0;

		  }



			  cdft(2*nfft,1,fso,ip,w);



		  for(ijs=0;ijs<=100000;ijs++){

            astout[ijs-itpout+100000].re=fso[ijs*2];

            astout[ijs-itpout+100000].im=fso[ijs*2+1];

          }

          for(ijss=0;ijss<=100000;ijss++){

            ast[1][ijss].re=astout[ijss+100000].re;

			ast[1][ijss].im=astout[ijss+100000].im;

          }



		   cdft(2*nfft,1,fio,ip,w);



		  for(iji=0;iji<=100000;iji++){

            aitout[iji-itpout+100000].re=fio[iji*2];

            aitout[iji-itpout+100000].im=fio[iji*2+1];

          }

          for(ijii=0;ijii<=100000;ijii++){

            ait[1][ijii].re=aitout[ijii+100000].re;

			ait[1][ijii].im=aitout[ijii+100000].im;

          }



          cdft(2*nfft,1,fpo,ip,w);



		  for(ijp=0;ijp<=100000;ijp++){

            aptout[ijp-itpout+100000].re=fpo[ijp*2];

            aptout[ijp-itpout+100000].im=fpo[ijp*2+1];

          }

          for(ijpp=0;ijpp<=100000;ijpp++){

            apt[1][ijpp].re=aptout[ijpp+100000].re;

			apt[1][ijpp].im=aptout[ijpp+100000].im;

          }



          cdft(2*nfft,1,fho,ip,w);



          for(ijh=0;ijh<=100000;ijh++){

            ahtout[ijh-ithout+100000].re=fho[ijh*2];

            ahtout[ijh-ithout+100000].im=fho[ijh*2+1];

          }

          for(ijhh=0;ijhh<=100000;ijhh++){

            aht[1][ijhh].re=ahtout[ijhh+100000].re;

			aht[1][ijhh].im=ahtout[ijhh+100000].im;

          }



/*		  printf("%le  ",zl);*/



		  tsops=0.0;

if(0.0<zout&&zout<1.0e-3){           ///////////////ifadd commentin soo///////////////
          for(soo=0;soo<=6000;soo=soo+1){

            fprintf(astd,"%le\t%le\t%le\t%le\t%le\t%le\n",zout*(1.0e3),tsops,cmplxmul3(pows,cmplxmul(ast[1][soo],cmplxconjg(ast[1][soo]))).re,ast[1][soo].re,ast[1][soo].im,zl);

			tsops=tsops+(1.0/dt)*1.0;

        }


}

          maxpwi=0.0;

          for(ima=0;ima<=100000;ima++){

            pwi=cmplxmul3(powi,cmplxmul(ait[1][ima],cmplxconjg(ait[1][ima]))).re;

	        if(maxpwi<pwi){

              maxpwi=pwi;

            }

          }



          maxe=maxpwi*exp(-1.0);

		  for(imae=0;imae<=100000;imae++){

            pwie=cmplxmul3(powi,cmplxmul(ait[1][imae],cmplxconjg(ait[1][imae]))).re;

            pwiie=cmplxmul3(powi,cmplxmul(ait[1][imae+1],cmplxconjg(ait[1][imae+1]))).re;

            if(pwie<=maxe&&maxe<=pwiie){

              itmi=imae;

            }

            if(pwie>=maxe&&maxe>=pwiie){

              itma=imae;

            }

          }

          dit=itma-itmi;



          tipsz=dit/dt;







          tiops=0.0;

		  if(9.5e-3<zout&&zout<10.5e-3){       ////////////ifnakami///////////////
          for(ioo=0;ioo<=6000;ioo=ioo+1){

			fprintf(aitd,"%le\t%le\t%le\t%le\t%le\t%le\t%le\n",zout*(1.0e3),tiops,cmplxmul3(powi,cmplxmul(ait[1][ioo],cmplxconjg(ait[1][ioo]))).re,ait[1][ioo].re,ait[1][ioo].im,maxpwi,tipsz);
			fprintf(syutumark,"%le\t%le\t%le\n",tiops,cmplxmul3(powi,cmplxmul(ait[1][ioo],cmplxconjg(ait[1][ioo]))).re,maxpwi);

            tiops=tiops+(1.0/dt)*1.0;

          }
		  }


 /*         tpops=0.0;

          for(poo=0;poo<=10000;poo=poo+1){

			fprintf(aptd,"%le\t%le\t%le\t%le\t%le\n",zout*(1.0e3),tpops,cmplxmul3(powp,cmplxmul(apt[1][poo],cmplxconjg(apt[1][poo]))).re,apt[1][poo].re,apt[1][poo].im);

            tpops=tpops+(1.0/dt)*1.0;

		  }*/





          thops=0.0;

          if(refl<=zl){              //////////////commentin///////////////

          for(hoo=0;hoo<=6000;hoo=hoo+1){    ///////////////////commentin hoo//////////////////////

			fprintf(ahtd,"%le\t%le\t%le\t%le\t%le\n",zout*(1.0e3),thops,cmplxmul3(powh,cmplxmul(aht[1][hoo],cmplxconjg(aht[1][hoo]))).re,aht[1][hoo].re,aht[1][hoo].im);

            thops=thops+(1.0/dt)*1.0;

          }

		  }





          zout=zout+1.0e-3;/*1.0e-3*/

        }





	}



	/*  fclose(devd);



	  fclose(aswd);

	  fclose(aiwd);

	  fclose(apwd);

	  fclose(ahwd);



	 

	  fclose(astd);

	  fclose(aitd);

	  fclose(aptd);

	  fclose(ahtd);*/

	  



	  

      

 /**     write(6,*)'The End'*/

/*

 *****    1000 format(1X,5e18.8e3)

 *****    2000 format(1X,6e18.8e3)

 *****    3000 format(1X,3e15.6e3)



      end

*/

 }

 

          fclose(devd);

	  fclose(aswd);

	  fclose(aiwd);

	  fclose(apwd);

	  fclose(ahwd);



	 

	  fclose(astd);

	  fclose(aitd);

	  fclose(aptd);

	  fclose(ahtd);

	  fclose(nyuumark);

	  fclose(syutumark);

	  

 

 

     printf("The End");

	return 0;



}

