TITLE Ih.mod   
:Kouranova et al 2008
:Kouranova-2008-Hyperpolarization-activated cyclic nucgated channel mRNA and protein expression in large versus

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 

NEURON {
        SUFFIX HcnKov
		USEION h READ eh WRITE ih VALENCE 1
        RANGE gbarfast, gbarslow, g, ih
        RANGE mtauf, mtausl, minf
		THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
	   gbarfast = 0.01105 (S/cm2) : 2 nS /pi*24^2
	   gbarslow = 0.0221 (S/cm2) : 2 nS /pi*24^2
       eh = -30 (mV)
}
 
STATE {
        mf msl
}
 
ASSIGNED {
	v (mV)
	celsius (degC)
	
	g (S/cm2)
	ih (mA/cm2)
     
    minf
	mtauf (ms)
	mtausl (ms)
}
 

BREAKPOINT {
        SOLVE states METHOD cnexp
        g = gbarfast*mf+gbarslow*msl
		ih = g*(v - eh)
}
 

INITIAL {
	rates(v)
	mf = minf
	msl = minf
	
}


DERIVATIVE states {  
        rates(v)
        mf'  =  (minf-mf)/mtauf
		msl' =  (minf-msl)/mtausl
}


PROCEDURE rates(v(mV)) {  
		LOCAL q10
UNITSOFF
		:q10 = 1.8^((celsius - 22(degC))/10(degC)): old q10
		q10 = 3^((celsius - 22(degC))/10(degC))
		minf = 1/(1+exp((v+87.2)/9.7)) : Kouronova 2008

		if (v < -70){
			mtauf = (1/q10)*(250 + 12*exp((v+240)/50))
		}
		else{
			mtauf = (1/q10)*(140 + 50*exp((v+25)/-20))
		}
		
		if (v < -70){
			mtausl = (1/q10)*(2500 + 100*exp((v+240)/50))
		}
		else{
			mtausl = (1/q10)*(300 + 542*exp((v+25)/-20))
		}
}
 
 
UNITSON
