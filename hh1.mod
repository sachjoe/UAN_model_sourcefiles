TITLE TRPV1 channel
COMMENT
Unit check passed only with units off
Author: Satchithananthi Aruljothi
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 
ENDCOMMENT

NEURON {
	SUFFIX hh1
	USEION  Ca READ eCa WRITE iCa VALENCE 2
	RANGE gmax, g, iCa,ntauamp, cap,ctauamp, gc, capson, gcmax, tempon,ica1,ica2
	GLOBAL ninf, ntau, cinf, ctau
	THREADSAFE
	}

UNITS {
(molar) = (1/liter)
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
	(mS) = (millimho)
	(nM) =	(nanomolar)
}

PARAMETER {
	gmax =0.001	(mho/cm2)	<0,1e9>
	gcmax (mho/cm2)
	eCa	        (mV)
	celsius = 22 (degC)
	ntauamp = 200
	ctauamp= 200
}

ASSIGNED {
	v	(mV)
	g	(mho/cm2)
	gc	(mho/cm2)
	iCa	(mA/cm2)
	ica1	(mA/cm2)
	ica2	(mA/cm2)
	ninf   (1)
	ntau  (ms)
	:ntauamp (1)
	cinf   (1)
	ctau  (ms)
	:ctauamp (1)
	:celsius (degC)
	cap(nM)
	capson(1)
	tempon(1)
}



STATE {

	n   c
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gmax*n
	gc= gcmax*c
	ica1=(g*(v-eCa))
	:ica2=capson*(gc*(v-eca))
	ica2=0
	iCa= ica1 + ica2


}

INITIAL {
	rates(v)
	n = ninf
	c = cinf
}

DERIVATIVE states {
	rates(v)

	n' = (ninf-n)/ntau
	c' = (cinf-c)/ctau
}

PROCEDURE rates(v(mV)) {
	LOCAL n_alpha, n_beta

UNITSOFF
	n_alpha = 1.61e37*(exp(-(208000/(8.314*(celsius+273.15)))))*(exp((0.5*0.71*9.65e4*v*1e-3)/(8.314*(celsius+273.15))))
	n_beta =  9.67e5*(exp(-(23200/(8.314*(celsius+273.15)))))*(exp(-((1-0.5)*0.71*9.65e4*v*1e-3)/(8.314*(celsius+273.15))))
	:gmax= -0.0181+(0.01588*exp(0.01456*celsius))

	ntau = (ntauamp/(n_alpha+n_beta))
	ctau = (ctauamp/(n_alpha+n_beta))

	ninf = (n_alpha/(n_alpha+n_beta))
	cinf =  (n_alpha/(n_alpha+n_beta))
}
UNITSON
