import time,math
from datetime import datetime


pi =3.141592653589793238462643
tpi = 2 * 3.141592653589793238462643
degs = 180.0/3.141592653589793238462643
rads = 3.141592653589793238462643/180.0

def C2K(temp_c):
  return temp_c + 273.16

def C2F(temp_c):
  return ((temp_c * 9.0)/5.0) + 32.0
      
def F2C(temp_f):
  return ((temp_f - 32.0) * 5.0) / 9.0

def dewpoint(celsius,humidity):
	RATIO = 373.15 / (273.15 + celsius)
	RHS = -7.90298 * (RATIO - 1)
	RHS += 5.02808 * math.log10(RATIO)
	RHS += -1.3816e-7 * (math.pow(10, (11.344 * (1 - 1/RATIO ))) - 1) 
	RHS += 8.1328e-3 * (math.pow(10, (-3.49149 * (RATIO - 1))) - 1) 
	RHS += math.log10(1013.246)
	VP = pow(10, RHS - 3) * humidity
	T = math.log(VP/0.61078)   
	return (241.88 * T) / (17.558 - T)

 
  
def fix(n):
   if n >= 0.0:      
     return floor(n)
   else:
      return ceil(n)

def normal(angle):
   rev = angle / 360.0
   rev = rev - fix(rev)
   if rev < 0.0:
      rev += 1.0
   return (rev*360.0)

def ora():
   return datetime.now().strftime('%Y-%m-%d %H:%M:%S')

def data_oggi():
   return datetime.now().strftime('%Y-%m-%d ')

def es(ta):
	# Hardy, R.; ITS-90 Formulations for Vapor Pressure, Frostpoint
	# Temperature, Temperature and Enhancement Factors in the
	# Range -100 to 100 C; 
	# Proceedings of Third International Symposium on Humidity and Moisture
	# edited by National Physical Laboratory (NPL), London, 1998, pp. 214-221
	# http:#www.thunderscientific.com/tech_info/reflibrary/its90formulas.pdf
	# (retrieved 2008-10-01)
 
	g= [ -2.8365744E3,-6.028076559E3,1.954263612E1,-2.737830188E-2,1.6261698E-5,7.0229056E-10,-1.8680009E-13,2.7150305]
				  
	tk = ta+273.15; # air temp in K
	es = g[7]*math.log10(tk)
	for x in range(0, 7):
		es = es+g[x]*math.pow(tk,float(x-2)); 
	es = math.exp(es)*0.01; # convert Pa to hPa
	return es


def dewpt(t,hpa):
	p=math.log(100.0*hpa)
	if t>=0.0:
		dpt=(-60.45)+7.0322*p+.37*p*p
	else:
		dpt=(-35.957)-1.8726*p+1.1689*p*p
	return (dpt)



def degsat(t, rh, pa):
	pws = rh/100*6.10* math.pow(2.718281828,( 17.27*t/( 237.7+t)))
	mu=rh*(1.0-pws/pa)/(1.0-rh*pws/pa)
	return mu




def heatindex( t, rh):
	if t<27.0: 
		return t
	else:
		tf= t * (9.0 / 5.0)+ 32.0
		tf2=math.pow(tf, 2.0)
		ur2=math.pow(rh, 2.0)
	hif= -42.379 + 2.04901523 * tf + 10.1433127 * rh - 0.22475541 *tf * rh- 6.83783 * 0.001* tf2 - 5.481717 * 0.01* ur2 +1.22874 * 0.001* tf2* rh+ 8.5282 * 0.0001* tf * ur2 -1.99 * 0.000001* tf2* ur2
	hif=((5.0 / 9.0) * (hif - 32.0))
	return hif




def p_vap( t, rh):
	t += 273.15
	if t>273.15:
		pvap=(math.exp(-5800.2206/t+1.3914993-.048640239*t+(.41764768e-4)*math.pow(t,
				2.0)-(.14452093e-7)*math.pow(t, 3.0)+6.5459673*math.log(t))/1000.0)
	else:
		pvap=(math.exp(-5674.5359/t+6.3925247-(.9677843e-2)*t+(.62215701e-6)*math.pow(
				t, 2.0)+(.20747825e-8)*math.pow(t, 3.0)-(.9484024e-12)*math.pow(t, 4.0)
				+4.1635019*math.log(t))/1000.0)
	pvap=pvap*(rh/100)			
	return pvap		

# Base metabolism in function of T   
           
def metabolism(t):
	return (-3.0909 * t + 203.64); 


def frostime(t,wind):
	if wind> 100.1 or wind < 0.0:
		return 999.9
	elif t > -10.0 or t < -60.0:
		return 999.9
	else:
		ft=(((-24.5*((0.667*wind)+4.8)+2111)*(math.pow((-t-4.8),-1.668)))/60)
	return ft


def p_saturazione(t):
	t = t+273.15
	p_saturazione=(math.exp(-5674.5359/t+6.3925247-(.9677843e-2)*t+(.62215701e-6)*math.pow(t, 2.0)+(.20747825e-8)*math.pow(t, 3.0)-(.9484024e-12)*math.pow(t, 4.0)+4.1635019*math.log(t))/1000.0)
	if t>273.15:
		p_saturazione=(math.exp(-5800.2206/t+1.3914993-.048640239*t+(.41764768e-4)*math.pow(t,2.0)-(.14452093e-7)*math.pow(t, 3.0)+6.5459673*math.log(t))/1000.0)
	return p_saturazione



def pmv_hoppe_iso( t, rh, wind, mtrad, iclo):
	eta = 0.01; # Mechanical efficiency
	age = 35.0; # Age
	mbody = 75.0; # Weigth in kg
	ht = 1.75; # Heigth in m
	tcl = 30.005
	MAX_LOOP = 200
	MAX_LOOP_HALF = MAX_LOOP / 2
	tcl_eps = 0.05
	eps = 0.97
	sigm = 5.67e-8
	adu = 0.203 * math.pow(mbody, 0.425) * math.pow(ht, 0.725)
	metbf = 3.19 * math.pow(mbody, (3.0/4.0)) * (1.0 + 0.004* (30.0- age)+0.018 * ((ht * 100.0/ math.pow(mbody, (1.0/3.0))) - 42.1)) 
	metbm = 3.45 * math.pow(mbody, (3.0/4.0)) * (1.0 + 0.004* (30.0- age)+0.010 * ((ht * 100.0/ math.pow(mbody, (1.0/3.0))) - 43.4)) 
	vpa = (rh / 100) * 6.105 * math.pow(2.718281828, ( 17.27*t / ( 237.7 + t ) ))
	fcl = 1.0 + iclo * 0.15
	metb = metabolism(t)
	metf = metbf + metb
	metm = metbm + metb
	metb = (metf+metm)/2.0
	h = metb * (1.0 - eta)
	aef = 0.71 * fcl * adu
	p1 = 35.7 - 0.032 * (metb / (adu * 1.16))* (1 - eta)
	tcl1 = tcl
	for x in range(0, MAX_LOOP_HALF):
		if x < MAX_LOOP_HALF:
			hc = 12.06 * math.sqrtf(wind)
			abhc = 0.0
		else:
			hc = 2.38 * math.pow(fabsf(tcl1 - t), 4.0)
			abhc = 0.6 * fabsf(math.pow((tcl1 - t), -0.75))
		
		tcl2 = p1 - 0.155 * iclo * (3.94 * 0.00000001* fcl *(math.pow((tcl1 + 273.2),4.0)- math.pow((mtrad+ 273.2), 4.0))+fcl * hc* (tcl1 - t))
		diff = math.abs(tcl1 - tcl2)
		if diff < tcl_eps:
			break
		abtcl = -0.155 * iclo * (4.0 * 3.94* 0.00000001* fcl *math.pow((tcl1+ 273.2),3.0) + fcl * hc- fcl *(tcl1 - t)* abhc)- 1.0
		tcl1 = tcl1 - (tcl2 - tcl1) / abtcl
		difhc = (12.06 * math.sqrt(wind)) - (2.38 * (math.pow(math.abs(t - tcl1), 0.25)))
		if difhc > 0.0 and i == MAX_LOOP_HALF:
			break
	tsk = 35.7 - (0.028 * h / adu)
	esw = 0.42 * adu * (h / adu - 58.08)
	if esw  < 0.0: 
		esw= 0.0
	rsum = aef * eps * sigm * (math.pow((tcl1 + 273.2), 4.0) - math.pow((mtrad + 273.2),4.0))
	csum = adu * fcl * hc * (tcl1 - t)
	erel = 0.0023 * metb * (44.0 - 0.75 * vpa)
	eres = 0.0014 * metb * (34.0 - t)
	ed = 0.406 * adu * (1.92 * tsk - 25.3- 0.75 * vpa)
	load = (h - ed - erel - eres - esw - rsum - csum) / adu
	ts = (0.303 * math.expf(-0.036 * (metb / adu)) + 0.028)
	pmv= ts * load
	return pmv



def poda(t,rh,p):
	vpa = (rh / 100) * 6.105 * math.pow(2.718281828, ( 17.27*t / ( 237.7 + t ) ))
	poda = 80.51 * p / (t + 273.15) * (1.0 - vpa / p)
	return poda



def thom(t, p_hPa):
	tw = (t*(0.45+(0.006*t*math.sqrt(p_hPa/1060.0))));  
	thom = 0.4*(t+tw)+4.8
	return thom



def utci(t,rh,wind,tmrt):
	utci_v=-99.9
	if t<-50.0 or t>50.0:
		return utci_v
	if tmrt<t-30.0 or tmrt>t+70.0:
		return utci_v
	if wind<0.5 or wind>30.0:
		return utci_v
	if rh<=0.0 or rh>=100.0:
		return utci_v
	va=wind
	ta=t
	e=es(ta)
	pa=e*(rh/100)
	pa=(pa/10)
	dtm=(tmrt-ta)
	utci_v = ta+(6.07562052E-01)+(-2.27712343E-02)*ta+(8.06470249E-04)*ta*ta+-1.54271372E-04 *ta*ta*ta+(-3.24651735E-06)*ta*ta*ta*ta+(7.32602852E-08)*ta*ta*ta*ta*ta+(1.35959073E-09)*ta*ta*ta*ta*ta*ta+(-2.25836520E+00)*va+(8.80326035E-02)*ta*va+(2.16844454E-03)*ta*ta*va+(-1.53347087E-05)*ta*ta*ta*va+(-5.72983704E-07)*ta*ta*ta*ta*va+(-2.55090145E-09)*ta*ta*ta*ta*ta*va+(-7.51269505E-01)*va*va+(-4.08350271E-03)*ta*va*va+(-5.21670675E-05)*ta*ta*va*va+(1.94544667E-06)*ta*ta*ta*va*va+(1.14099531E-08)*ta*ta*ta*ta*va*va+(1.58137256E-01)*va*va*va+(-6.57263143E-05)*ta*va*va*va+(2.22697524E-07)*ta*ta*va*va*va+(-4.16117031E-08)*ta*ta*ta*va*va*va+(-1.27762753E-02)*va*va*va*va+(9.66891875E-06)*ta*va*va*va*va+(2.52785852E-09)*ta*ta*va*va*va*va+(4.56306672E-04)*va*va*va*va*va+(-1.74202546E-07)*ta*va*va*va*va*va+(-5.91491269E-06)*va*va*va*va*va*va+(3.98374029E-01)*dtm+(1.83945314E-04)*ta*dtm+(-1.73754510E-04)*ta*ta*dtm+(-7.60781159E-07)*ta*ta*ta*dtm+(3.77830287E-08)*ta*ta*ta*ta*dtm+(5.43079673E-10)*ta*ta*ta*ta*ta*dtm+(-2.00518269E-02)*va*dtm+(8.92859837E-04)*ta*va*dtm+(3.45433048E-06)*ta*ta*va*dtm+(-3.77925774E-07)*ta*ta*ta*va*dtm+(-1.69699377E-09)*ta*ta*ta*ta*va*dtm+(1.69992415E-04)*va*va*dtm+(-4.99204314E-05)*ta*va*va*dtm+(2.47417178E-07)*ta*ta*va*va*dtm+(1.07596466E-08)*ta*ta*ta*va*va*dtm+(8.49242932E-05)*va*va*va*dtm+(1.35191328E-06)*ta*va*va*va*dtm+(-6.21531254E-09)*ta*ta*va*va*va*dtm+(-4.99410301E-06)*va*va*va*va*dtm+(-1.89489258E-08)*ta*va*va*va*va*dtm+(8.15300114E-08)*va*va*va*va*va*dtm+(7.55043090E-04)*dtm*dtm+(-5.65095215E-05)*ta*dtm*dtm+(-4.52166564E-07)*ta*ta*dtm*dtm+(2.46688878E-08)*ta*ta*ta*dtm*dtm+(2.42674348E-10)*ta*ta*ta*ta*dtm*dtm+(1.54547250E-04)*va*dtm*dtm+(5.24110970E-06)*ta*va*dtm*dtm+(-8.75874982E-08)*ta*ta*va*dtm*dtm+(-1.50743064E-09)*ta*ta*ta*va*dtm*dtm+(-1.56236307E-05)*va*va*dtm*dtm+(-1.33895614E-07)*ta*va*va*dtm*dtm+(2.49709824E-09)*ta*ta*va*va*dtm*dtm+(6.51711721E-07)*va*va*va*dtm*dtm+(1.94960053E-09)*ta*va*va*va*dtm*dtm+(-1.00361113E-08)*va*va*va*va*dtm*dtm+(-1.21206673E-05)*dtm*dtm*dtm+(-2.18203660E-07)*ta*dtm*dtm*dtm+(7.51269482E-09)*ta*ta*dtm*dtm*dtm+(9.79063848E-11)*ta*ta*ta*dtm*dtm*dtm+(1.25006734E-06)*va*dtm*dtm*dtm+(-1.81584736E-09)*ta*va*dtm*dtm*dtm+(-3.52197671E-10)*ta*ta*va*dtm*dtm*dtm+(-3.36514630E-08)*va*va*dtm*dtm*dtm+(1.35908359E-10)*ta*va*va*dtm*dtm*dtm+(4.17032620E-10)*va*va*va*dtm*dtm*dtm+(-1.30369025E-09)*dtm*dtm*dtm*dtm+(4.13908461E-10)*ta*dtm*dtm*dtm*dtm+(9.22652254E-12)*ta*ta*dtm*dtm*dtm*dtm+(-5.08220384E-09)*va*dtm*dtm*dtm*dtm+(-2.24730961E-11)*ta*va*dtm*dtm*dtm*dtm+(1.17139133E-10)*va*va*dtm*dtm*dtm*dtm+(6.62154879E-10)*dtm*dtm*dtm*dtm*dtm+(4.03863260E-13)*ta*dtm*dtm*dtm*dtm*dtm+(1.95087203E-12)*va*dtm*dtm*dtm*dtm*dtm+(-4.73602469E-12)*dtm*dtm*dtm*dtm*dtm*dtm+(5.12733497E+00)*pa+(-3.12788561E-01)*ta*pa+(-1.96701861E-02)*ta*ta*pa+(9.99690870E-04)*ta*ta*ta*pa+(9.51738512E-06)*ta*ta*ta*ta*pa+(-4.66426341E-07)*ta*ta*ta*ta*ta*pa+(5.48050612E-01)*va*pa+(-3.30552823E-03)*ta*va*pa+(-1.64119440E-03)*ta*ta*va*pa+(-5.16670694E-06)*ta*ta*ta*va*pa+(9.52692432E-07)*ta*ta*ta*ta*va*pa+(-4.29223622E-02)*va*va*pa+(5.00845667E-03)*ta*va*va*pa+(1.00601257E-06)*ta*ta*va*va*pa+(-1.81748644E-06)*ta*ta*ta*va*va*pa+(-1.25813502E-03)*va*va*va*pa+(-1.79330391E-04)*ta*va*va*va*pa+(2.34994441E-06)*ta*ta*va*va*va*pa+(1.29735808E-04)*va*va*va*va*pa+(1.29064870E-06)*ta*va*va*va*va*pa+(-2.28558686E-06)*va*va*va*va*va*pa+(-3.69476348E-02)*dtm*pa+(1.62325322E-03)*ta*dtm*pa+(-3.14279680E-05)*ta*ta*dtm*pa+(2.59835559E-06)*ta*ta*ta*dtm*pa+(-4.77136523E-08)*ta*ta*ta*ta*dtm*pa+(8.64203390E-03)*va*dtm*pa+(-6.87405181E-04)*ta*va*dtm*pa+(-9.13863872E-06)*ta*ta*va*dtm*pa+(5.15916806E-07)*ta*ta*ta*va*dtm*pa+(-3.59217476E-05)*va*va*dtm*pa+(3.28696511E-05)*ta*va*va*dtm*pa+(-7.10542454E-07)*ta*ta*va*va*dtm*pa+(-1.24382300E-05)*va*va*va*dtm*pa+(-7.38584400E-09)*ta*va*va*va*dtm*pa+(2.20609296E-07)*va*va*va*va*dtm*pa+(-7.32469180E-04)*dtm*dtm*pa+(-1.87381964E-05)*ta*dtm*dtm*pa+(4.80925239E-06)*ta*ta*dtm*dtm*pa+(-8.75492040E-08)*ta*ta*ta*dtm*dtm*pa+(2.77862930E-05)*va*dtm*dtm*pa+(-5.06004592E-06)*ta*va*dtm*dtm*pa+(1.14325367E-07)*ta*ta*va*dtm*dtm*pa+(2.53016723E-06)*va*va*dtm*dtm*pa+(-1.72857035E-08)*ta*va*va*dtm*dtm*pa+(-3.95079398E-08)*va*va*va*dtm*dtm*pa+(-3.59413173E-07)*dtm*dtm*dtm*pa+(7.04388046E-07)*ta*dtm*dtm*dtm*pa+(-1.89309167E-08)*ta*ta*dtm*dtm*dtm*pa+(-4.79768731E-07)*va*dtm*dtm*dtm*pa+(7.96079978E-09)*ta*va*dtm*dtm*dtm*pa+(1.62897058E-09)*va*va*dtm*dtm*dtm*pa+(3.94367674E-08)*dtm*dtm*dtm*dtm*pa+(-1.18566247E-09)*ta*dtm*dtm*dtm*dtm*pa+(3.34678041E-10)*va*dtm*dtm*dtm*dtm*pa+(-1.15606447E-10)*dtm*dtm*dtm*dtm*dtm*pa+(-2.80626406E+00)*pa*pa+(5.48712484E-01)*ta*pa*pa+(-3.99428410E-03)*ta*ta*pa*pa+(-9.54009191E-04)*ta*ta*ta*pa*pa+(1.93090978E-05)*ta*ta*ta*ta*pa*pa+(-3.08806365E-01)*va*pa*pa+(1.16952364E-02)*ta*va*pa*pa+(4.95271903E-04)*ta*ta*va*pa*pa+(-1.90710882E-05)*ta*ta*ta*va*pa*pa+(2.10787756E-03)*va*va*pa*pa+(-6.98445738E-04)*ta*va*va*pa*pa+(2.30109073E-05)*ta*ta*va*va*pa*pa+(4.17856590E-04)*va*va*va*pa*pa+(-1.27043871E-05)*ta*va*va*va*pa*pa+(-3.04620472E-06)*va*va*va*va*pa*pa+(5.14507424E-02)*dtm*pa*pa+(-4.32510997E-03)*ta*dtm*pa*pa+(8.99281156E-05)*ta*ta*dtm*pa*pa+(-7.14663943E-07)*ta*ta*ta*dtm*pa*pa+(-2.66016305E-04)*va*dtm*pa*pa+(2.63789586E-04)*ta*va*dtm*pa*pa+(-7.01199003E-06)*ta*ta*va*dtm*pa*pa+(-1.06823306E-04)*va*va*dtm*pa*pa+(3.61341136E-06)*ta*va*va*dtm*pa*pa+(2.29748967E-07)*va*va*va*dtm*pa*pa+(3.04788893E-04)*dtm*dtm*pa*pa+(-6.42070836E-05)*ta*dtm*dtm*pa*pa+(1.16257971E-06)*ta*ta*dtm*dtm*pa*pa+(7.68023384E-06)*va*dtm*dtm*pa*pa+(-5.47446896E-07)*ta*va*dtm*dtm*pa*pa+(-3.59937910E-08)*va*va*dtm*dtm*pa*pa+(-4.36497725E-06)*dtm*dtm*dtm*pa*pa+(1.68737969E-07)*ta*dtm*dtm*dtm*pa*pa+(2.67489271E-08)*va*dtm*dtm*dtm*pa*pa+(3.23926897E-09)*dtm*dtm*dtm*dtm*pa*pa+(-3.53874123E-02)*pa*pa*pa+(-2.21201190E-01)*ta*pa*pa*pa+(1.55126038E-02)*ta*ta*pa*pa*pa+(-2.63917279E-04)*ta*ta*ta*pa*pa*pa+(4.53433455E-02)*va*pa*pa*pa+(-4.32943862E-03)*ta*va*pa*pa*pa+(1.45389826E-04)*ta*ta*va*pa*pa*pa+(2.17508610E-04)*va*va*pa*pa*pa+(-6.66724702E-05)*ta*va*va*pa*pa*pa+(3.33217140E-05)*va*va*va*pa*pa*pa+(-2.26921615E-03)*dtm*pa*pa*pa+(3.80261982E-04)*ta*dtm*pa*pa*pa+(-5.45314314E-09)*ta*ta*dtm*pa*pa*pa+(-7.96355448E-04)*va*dtm*pa*pa*pa+(2.53458034E-05)*ta*va*dtm*pa*pa*pa+(-6.31223658E-06)*va*va*dtm*pa*pa*pa+(3.02122035E-04)*dtm*dtm*pa*pa*pa+(-4.77403547E-06)*ta*dtm*dtm*pa*pa*pa+(1.73825715E-06)*va*dtm*dtm*pa*pa*pa+(-4.09087898E-07)*dtm*dtm*dtm*pa*pa*pa+(6.14155345E-01)*pa*pa*pa*pa+(-6.16755931E-02)*ta*pa*pa*pa*pa+(1.33374846E-03)*ta*ta*pa*pa*pa*pa+(3.55375387E-03)*va*pa*pa*pa*pa+(-5.13027851E-04)*ta*va*pa*pa*pa*pa+(1.02449757E-04)*va*va*pa*pa*pa*pa+(-1.48526421E-03)*dtm*pa*pa*pa*pa+(-4.11469183E-05)*ta*dtm*pa*pa*pa*pa+(-6.80434415E-06)*va*dtm*pa*pa*pa*pa+(-9.77675906E-06)*dtm*dtm*pa*pa*pa*pa+(8.82773108E-02)*pa*pa*pa*pa*pa+(-3.01859306E-03)*ta*pa*pa*pa*pa*pa+(1.04452989E-03)*va*pa*pa*pa*pa*pa+(2.47090539E-04)*dtm*pa*pa*pa*pa*pa+(1.48348065E-03)*pa*pa*pa*pa*pa*pa
	return utci_v


def wbgt( t, rh, wind):
	wbgt = 999.9
	if rh> 100.1 or rh < 0.0:
		return wbgt
	elif t > 100.0 or t < -100.0:
		return wbgt
	else:
		e = (rh/100.0)*(6.105*math.exp((t*17.27)/(237.7+t)))
		wbgt = (0.567*t)+(0.393*e)+3.94;
	return wbgt



def proj(sunelev):
	if sunelev < 0.0:
		return 0.0
	return 0.308 * cos(rads * (sunelev* (0.998- (math.pow(sunelev, 2.0) / 50000.0))))


def temprad(t,rh,rshort,rdiffuse,sunelev,albedo):
	e = (rh/100.0)*(6.105*math.exp((t*17.27)/(237.7+t)))
	sig = 5.67e-8
	emiair = 0.66 + 0.039 * math.sqrtf(e)
	tsk = t + 273.12
	ratio=0.0429*sin(sunelev*rads)+0.345*cos(sunelev*rads)
	proj=0.308 * cos(rads * (sunelev* (0.998- (math.pow(sunelev, 2.0) / 50000.0))))
	temprad= math.pow(273.16 - (emiair * math.pow(tsk, 4.0) + (1-albedo) * (rdiffuse) / (sig* 0.97)+(1-albedo) * proj * ratio* ((rshort-rdiffuse)/(sig*0.97))),0.25)- 273.16
	return temprad


def humidex( t,rh):
	humidex = 999.9
	if rh > 100.1 or rh < 0.0:
		return humidex
	elif t > 100.0 or t < -100.0:
		return humidex
	else:
		e = (rh/100.0)*(6.105*math.exp((t*17.27)/(237.7+t)))
	humidex = t+(0.5555*(e-10.0))
	return humidex


def hi( t, rh):
	hi = 999.9
	if rh> 100.1 or rh < 0.0:
		return hi
	elif t > 100.0 or t < -100.0:
		return hi
	else:
		hi = -8.784695+(1.61139411*t)+(2.338549*rh)-(0.14611605*t*rh)-(1.2308094*math.pow(10,-2.0)*math.pow(t,2.0))-(1.6424828*math.pow(10,-2.0)*math.pow(rh,2.0))+(2.211732*math.pow(10,-3.0)*math.pow(t,2.0)*rh)+(7.2546*math.pow(10,-4.0)*t*math.pow(rh,2.0))-(3.582*math.pow(10,-6.0)*math.pow(rh,2.0))
	return hi


def net(t,rh,wind):
	net = 999.9
	if rh > 100.1 or rh < 0.0:
		return net
	elif wind > 130.0 or wind < 0.0:
		return 999.9
	elif t > 100.0 or t < -100.0:
		return 999.9
	else:
		net = 37-((37-t)/(0.68-(0.0014*rh)+(1/(1.76+(1.4*(math.pow(wind,0.75)))))))-(0.29*t*(1.0-(0.01*rh)))
	return net



def ssi( t,rh):
	ssi = 999.9
	if rh > 100.1 or rh < 0.0:
		return ssi
	elif t > 100.0 or t < -100.0:
		return ssi
	else:
		ssi = ((1.98*((((9.0/5.0)*t)+32.0)-(0.55-0.0055*rh)*((((9.0/5.0)*t)+32.0)-58.0))-56.83)-32.0)/1.8
	return  ssi


def steadman_outdoor_sun( t, rh, wind, rshort, sunelev):
	steadman_outdoor_sun=-999.9
	ee = (rh/1000.0)*(6.105*math.exp((t*17.27)/(237.7+t)))
	q_glob = 0.56*(0.386-(0.0032*sunelev))*rshort + 0.224*(0.1*rshort)+ 0.028*rshort - 150.0*(0.38-0.16*(math.pow(ee,0.5)));  
	if rh > 100.1 or rh < 0.0:
		return steadman_outdoor_sun
	elif t > 100.0 or t < -100.0:
		return steadman_outdoor_sun 	  
	if q_glob > 0.0:
		steadman_outdoor_sun = t+3.48*(ee)-0.7*wind +0.7*q_glob/(wind+10.0)-4.25
	return steadman_outdoor_sun


def steadman_indoor(t, rh):
	steadman_indoor = 999.9
	if rh > 100.1 or rh < 0.0:
		return steadman_indoor
	elif t > 100.0 or t < -100.0:
		return steadman_indoor
	else:
		e = (rh/100.0)*(6.105*math.exp((t*17.27)/(237.7+t)))
	steadman_indoor = -2.56+(0.89*t)+(0.382*e);  
	return steadman_indoor


def steadman_outdoor_shade( t, rh, wind):
	steadman_outdoor_shade = 999.9
	if rh > 100.1 or rh < 0.0:
		return steadman_outdoor_shade
	elif wind > 130.0 or wind < 0.0:
		return steadman_outdoor_shade
	elif t > 100.0 or t < -100.0:
		return steadman_outdoor_shade
	else:
		e = (rh/100.0)*(6.105*math.exp((t*17.27)/(237.7+t)))
	steadman_outdoor_shade = t+(0.33*e)-(0.7*wind)-4.0
	return steadman_outdoor_shade



def ppd( pmv):
	ppd = 100.0 - 95.0 * math.expf(-0.2179 * math.pow(pmv, 2.0)) - 0.03353* math.pow(pmv, 4.0)
	return ppd

def new_windchill(t,wind):
	nsw = 13.12+(0.6215*t)-(11.37*math.pow(wind,0.16))+(0.3965*t*math.pow(wind,0.16)) 
	return nsw


def t_apparent_aus( t, rh, wind):
	e = (rh/100.0)*(6.105*math.exp((t*17.27)/(237.7+t)))
	t_app = t +0.33*e-0.70*wind-4
	return t_app

def clomin( t, rh, wind, trad):
	MAX_ITER = 40
	PMV_GOOD = 0.5
	pmv = -1.0
	clomin = 0.1
	for x in range(0,MAX_ITER):
		pmv = pmv_hoppe_iso(t, rh, wind, trad, clomin)
		if pmv > PMV_GOOD:
			break
		clomin += 0.1
	return clomin


def clomax( t, rh, wind, trad):
	MAX_ITER = 40
	PMV_GOOD = 0.5
	pmv = 1.0
	clomax = 5
	for x in range(0,MAX_ITER):
		pmv = pmv_hoppe_iso(t, rh, wind, trad, clomax)
		if pmv < PMV_GOOD:
			break
		clomax -= 0.1
	return clomax

def rdiffuse(radteoric, rshort):
	if radteoric<=0.0:
		rdiffuse=0.0
	else:
		kg=rshort / radteoric
		kd=(0.365 - 0.25 * kg) * sin(tpi * kg)
		if (kd * radteoric)>0.0:
			rdiffuse=kd * radteoric
		else:
			rdiffuse=0.0
		if (rdiffuse > rshort):
			rdiffuse=rshort
	return rdiffuse


def ta_comfort(rh, iclo, wind, M, tskin):
	t = 40.0; # Initial estimation value
	Ra = 1 / 9
	fcl = 1+.31 * iclo
	Tsk =35.7 - .0285 * M - 3
	if H == 999.9:
		Tsk = 35.7 - .0285 * M
	Rst = iclo * .155 + Ra / fcl
	WS = .0052 * (M - 58)
	if WS>.7:
		WS = .7
	corr = math.exp(.043 - .398 * wind + .066 * wind * wind- .378 * WS +.094 * WS * WS)
	Rdyn = Rst * corr
	for x in range(0,1000):
		# Calculation of Convective heat loss from the skin
		t=t-0.1
		C = (Tsk - t) / Rdyn
		# Calculation of radiation heat exchange
		hr = 5.67E-08 * .97 * .77 * (math.exp(4 * math.log(Tsk +273.15)) - math.exp(4 * math.log(t + 273.15))) / (Tsk - t)
		hc = 8.7 * math.exp (.6 * math.log(wind))
		if wind<1:
			hc = 3.5+5.2 * wind
		Fcl = 1 / ((hc + hr) * iclo * .155 + 1 / fcl)
		R = hr * Fcl * (Tsk - t)
		# Calculation of Evaporative Heat Loss from the Skin
		Psk = .1333 * math.exp(18.6686 - 4030.183 / (Tsk + 235))
		Pa = rh * .1333 * math.exp(18.6686 - 4030.183 / (t +235)) / 100
		Im = .38 * (4.9 - 6.5 * corr + 2.6 * corr * corr)
		if Im>.9: 
			Im=.9;	
		Retdyn = (Rdyn / Im )/ 16.65
		w = .001 * M
		E = w * (Psk - Pa) / Retdyn
		# Calculation of Convective Heat Loss from Respiration
		mres = 2.58 * .000001 * M
		Tex = 29 +.2 * t
		Cres = 1007 * mres * (Tex - t) / 1.8
		# Calculation of Evaporative Heat Loss from Respiration
		Wa = .622 * Pa / (101.325 - Pa)
		Pex = .1333 * math.exp(18.6686 - 4030.183 / (Tex + 235))
		Wex = .622 * Pex / (101.325 - Pex)
		Eres = 2423000 * mres * (Wex - Wa) / 1.8	
		# Calculation of heat debt or heat storage
		S =1
		if H == -999.0:
			S = 0
		else: 
			S = 40 / H
		balance = M - C - R - E - Cres - Eres - S
		if  balance < 0:
			break	
	return t-0.1



def p_local(press, topo, temp):
	temp=temp+273
	L=-0.0065;# temperature lapse L = -0.0065 K/m
	R_cost=287.05 ;#gas constant for dry air, J/(kg*degK) = 287.05
	T0=temp-(L/2)*topo;# sea level standard T0 = 288.15 K
	p_local = press*math.exp(-topo*9.81/(R_cost*T0))
	return p_local 
 
def utci_class_10(t,tmrt,wind,rh): 
	utci_v=utci(t,tmrt,wind,rh)
	if t<-90:
		return -99.9
	if t>90:
		return -99.9
	utci_c=-99.9
	if utci_v > 46.0:
		utci_c=10.0
	elif utci_v>38.0 and utci_v<=46.0:
		utci_c=9.0
	elif utci_v>32.0 and utci_v<=38.0:
		utci_c=8.0
	elif utci_v>26.0 and utci_v<=32.0:
		utci_c=7.0
	elif utci_v >9.0 and utci_v<=26.0:
		utci_c=6.0
	elif utci_v>0.0 and utci_v<=9.0:
		utci_c=5.0
	elif utci_v>-13.0 and utci_v<=0:
		utci_c=4.0
	elif utci_v>-27.0 and utci_v<=-13.0:
		utci_c=3.0
	elif utci_v>-40.0 and utci_v<=-27.0:
		utci_c=2.0
	elif utci_v<=-40.0:
		utci_c=1.0
	return utci_c
	
def utci_class_7(t,tmrt,wind,rh): 
	utci_v=utci(t,tmrt,wind,rh)
	if t<-90:
		return -99.9
	if t>90:
		return -99.9
	utci_c=-99.9
	if utci_v > 46.0:
		utci_c=7.0
	elif utci_v>38.0 and utci_v<=46.0:
		utci_c=7.0
	elif utci_v>32.0 and utci_v<=38.0:
		utci_c=7.0
	elif utci_v>26.0 and utci_v<=32.0:
		utci_c=6.0
	elif utci_v >16.0 and utci_v<=26.0:
		utci_c=5.0
	elif utci_v>0.0 and utci_v<=16.0:
		utci_c=4.0
	elif utci_v>-13.0 and utci_v<=0:
		utci_c=3.0
	elif utci_v>-27.0 and utci_v<=-13.0:
		utci_c=2.0
	elif utci_v>-40.0 and utci_v<=-27.0:
		utci_c=1.0
	elif utci_v<=-40.0:
		utci_c=1.0
	return utci_c
	
