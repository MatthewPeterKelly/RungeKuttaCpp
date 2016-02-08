#include <iostream>
#include <fstream>
using namespace std;

#include <cmath>

#include "integrator.h"


void printState(std::ofstream& file, double t, double z[], int nDim) {
	file << t ;
	for (int j = 0; j < nDim; j++) {
		file << ", " << z[j];
	}
	file << "\n";
}


/* Takes a simple euler step for the system */
void eulerStep(DynFun dynFun, double t0, double t1, double z0[], double z1[], int nDim) {

	double dt = t1 - t0;
	double *dz;
	dz = new double[nDim];

	dynFun(t0, z0, dz);

	for (int i = 0; i < nDim; i++) {
		z1[i] = z0[i] + dt * dz[i];
	}

	delete [] dz;

}


/* Time step using the mid-point method */
void midPointStep(DynFun dynFun, double tLow, double tUpp, double zLow[], double zUpp[], int nDim) {

	double dt = tUpp - tLow;

	/// Stage 0
	double t0 = tLow;
	double *z0 = zLow;
	double *f0; f0 = new double[nDim];
	dynFun(t0, z0, f0);

	/// Stage 1
	double t1 = t0 + 0.5 * dt;
	double *z1; z1 = new double[nDim];
	double *f1; f1 = new double[nDim];
	for (int i = 0; i < nDim; i++) {
		z1[i] = z0[i] + 0.5 * dt * f0[i];
	}
	dynFun(t1, z1, f1);

	/// Collect Stages:
	for (int i = 0; i < nDim; i++) {
		zUpp[i] = zLow[i] + dt * f1[i];
	}

	delete [] f0;
	delete [] f1;
	delete [] z1;

}


/* Time step using 4th-order Runge Kutta */
void rungeKuttaStep(DynFun dynFun, double tLow, double tUpp, double zLow[], double zUpp[], int nDim) {

	double dt = tUpp - tLow;

	/// Stage 0
	double t0 = tLow;
	double *z0 = zLow;
	double *f0; f0 = new double[nDim];
	dynFun(t0, z0, f0);

	/// Stage 1
	double t1 = t0 + 0.5 * dt;
	double *z1; z1 = new double[nDim];
	double *f1; f1 = new double[nDim];
	for (int i = 0; i < nDim; i++) {
		z1[i] = zLow[i] + 0.5 * dt * f0[i];
	}
	dynFun(t1, z1, f1);

	/// Stage 2
	double t2 = t1;
	double *z2; z2 = new double[nDim];
	double *f2; f2 = new double[nDim];
	for (int i = 0; i < nDim; i++) {
		z2[i] = zLow[i] + 0.5 * dt * f1[i];
	}
	dynFun(t2, z2, f2);

	/// Stage 3
	double t3 = tLow + dt;
	double *z3; z3 = new double[nDim];
	double *f3; f3 = new double[nDim];
	for (int i = 0; i < nDim; i++) {
		z3[i] = zLow[i] + dt * f2[i];
	}
	dynFun(t3, z3, f3);

	/// Collect Stages:
	for (int i = 0; i < nDim; i++) {
		zUpp[i] = zLow[i] + (dt / 6) * (f0[i] + 2.0 * f1[i] + 2.0 * f2[i] + f3[i]);
	}

	delete [] f0;
	delete [] f1;
	delete [] z1;
	delete [] f2;
	delete [] z2;
	delete [] f3;
	delete [] z3;

}



/* Time step using 8th-order Runge-Kutta-Fehlberg. Data for method from:
http://sce.uhcl.edu/rungekutta/rk108.txt
Paper:  "A tenth-order Runge-Kutta method with error estimate"
By Feagin*/

static double RK_A[] = {
	0.000000000000000000000000000000000000000000000000000000000000,
	0.100000000000000000000000000000000000000000000000000000000000,
	0.539357840802981787532485197881302436857273449701009015505500,
	0.809036761204472681298727796821953655285910174551513523258250,
	0.309036761204472681298727796821953655285910174551513523258250,
	0.981074190219795268254879548310562080489056746118724882027805,
	0.833333333333333333333333333333333333333333333333333333333333,
	0.354017365856802376329264185948796742115824053807373968324184,
	0.882527661964732346425501486979669075182867844268052119663791,
	0.642615758240322548157075497020439535959501736363212695909875,
	0.357384241759677451842924502979560464040498263636787304090125,
	0.117472338035267653574498513020330924817132155731947880336209,
	0.833333333333333333333333333333333333333333333333333333333333,
	0.309036761204472681298727796821953655285910174551513523258250,
	0.539357840802981787532485197881302436857273449701009015505500,
	0.100000000000000000000000000000000000000000000000000000000000,
	1.00000000000000000000000000000000000000000000000000000000000
};

static double RK_C[] = {
	0.0333333333333333333333333333333333333333333333333333333333333,
	0.0250000000000000000000000000000000000000000000000000000000000,
	0.0333333333333333333333333333333333333333333333333333333333333,
	0.000000000000000000000000000000000000000000000000000000000000,
	0.0500000000000000000000000000000000000000000000000000000000000,
	0.000000000000000000000000000000000000000000000000000000000000,
	0.0400000000000000000000000000000000000000000000000000000000000,
	0.000000000000000000000000000000000000000000000000000000000000,
	0.189237478148923490158306404106012326238162346948625830327194,
	0.277429188517743176508360262560654340428504319718040836339472,
	0.277429188517743176508360262560654340428504319718040836339472,
	0.189237478148923490158306404106012326238162346948625830327194,
	-0.0400000000000000000000000000000000000000000000000000000000000,
	-0.0500000000000000000000000000000000000000000000000000000000000,
	-0.0333333333333333333333333333333333333333333333333333333333333,
	-0.0250000000000000000000000000000000000000000000000000000000000,
	0.0333333333333333333333333333333333333333333333333333333333333
};

static double RK_B[17][16] = {0.0};
static bool flagInitRK = false;
void initializeRungeKuttaConstants(void) {
	RK_B[1][0] = 0.100000000000000000000000000000000000000000000000000000000000;
	RK_B[2][0] = -0.915176561375291440520015019275342154318951387664369720564660;
	RK_B[2][1] = 1.45453440217827322805250021715664459117622483736537873607016;
	RK_B[3][0] = 0.202259190301118170324681949205488413821477543637878380814562;
	RK_B[3][1] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[3][2] = 0.606777570903354510974045847616465241464432630913635142443687;
	RK_B[4][0] = 0.184024714708643575149100693471120664216774047979591417844635;
	RK_B[4][1] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[4][2] = 0.197966831227192369068141770510388793370637287463360401555746;
	RK_B[4][3] = -0.0729547847313632629185146671595558023015011608914382961421311;
	RK_B[5][0] = 0.0879007340206681337319777094132125475918886824944548534041378;
	RK_B[5][1] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[5][2] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[5][3] = 0.410459702520260645318174895920453426088035325902848695210406;
	RK_B[5][4] = 0.482713753678866489204726942976896106809132737721421333413261;
	RK_B[6][0] = 0.0859700504902460302188480225945808401411132615636600222593880;
	RK_B[6][1] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[6][2] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[6][3] = 0.330885963040722183948884057658753173648240154838402033448632;
	RK_B[6][4] = 0.489662957309450192844507011135898201178015478433790097210790;
	RK_B[6][5] = -0.0731856375070850736789057580558988816340355615025188195854775;
	RK_B[7][0] = 0.120930449125333720660378854927668953958938996999703678812621;
	RK_B[7][1] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[7][2] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[7][3] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[7][4] = 0.260124675758295622809007617838335174368108756484693361887839;
	RK_B[7][5] = 0.0325402621549091330158899334391231259332716675992700000776101;
	RK_B[7][6] = -0.0595780211817361001560122202563305121444953672762930724538856;
	RK_B[8][0] = 0.110854379580391483508936171010218441909425780168656559807038;
	RK_B[8][1] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[8][2] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[8][3] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[8][4] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[8][5] = -0.0605761488255005587620924953655516875526344415354339234619466;
	RK_B[8][6] = 0.321763705601778390100898799049878904081404368603077129251110;
	RK_B[8][7] = 0.510485725608063031577759012285123416744672137031752354067590;
	RK_B[9][0] = 0.112054414752879004829715002761802363003717611158172229329393;
	RK_B[9][1] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[9][2] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[9][3] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[9][4] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[9][5] = -0.144942775902865915672349828340980777181668499748506838876185;
	RK_B[9][6] = -0.333269719096256706589705211415746871709467423992115497968724;
	RK_B[9][7] = 0.499269229556880061353316843969978567860276816592673201240332;
	RK_B[9][8] = 0.509504608929686104236098690045386253986643232352989602185060;
	RK_B[10][0] = 0.113976783964185986138004186736901163890724752541486831640341;
	RK_B[10][1] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[10][2] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[10][3] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[10][4] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[10][5] = -0.0768813364203356938586214289120895270821349023390922987406384;
	RK_B[10][6] = 0.239527360324390649107711455271882373019741311201004119339563;
	RK_B[10][7] = 0.397774662368094639047830462488952104564716416343454639902613;
	RK_B[10][8] = 0.0107558956873607455550609147441477450257136782823280838547024;
	RK_B[10][9] = -0.327769124164018874147061087350233395378262992392394071906457;
	RK_B[11][0] = 0.0798314528280196046351426864486400322758737630423413945356284;
	RK_B[11][1] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[11][2] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[11][3] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[11][4] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[11][5] = -0.0520329686800603076514949887612959068721311443881683526937298;
	RK_B[11][6] = -0.0576954146168548881732784355283433509066159287152968723021864;
	RK_B[11][7] = 0.194781915712104164976306262147382871156142921354409364738090;
	RK_B[11][8] = 0.145384923188325069727524825977071194859203467568236523866582;
	RK_B[11][9] = -0.0782942710351670777553986729725692447252077047239160551335016;
	RK_B[11][10] = -0.114503299361098912184303164290554670970133218405658122674674;
	RK_B[12][0] = 0.985115610164857280120041500306517278413646677314195559520529;
	RK_B[12][1] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[12][2] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[12][3] = 0.330885963040722183948884057658753173648240154838402033448632;
	RK_B[12][4] = 0.489662957309450192844507011135898201178015478433790097210790;
	RK_B[12][5] = -1.37896486574843567582112720930751902353904327148559471526397;
	RK_B[12][6] = -0.861164195027635666673916999665534573351026060987427093314412;
	RK_B[12][7] = 5.78428813637537220022999785486578436006872789689499172601856;
	RK_B[12][8] = 3.28807761985103566890460615937314805477268252903342356581925;
	RK_B[12][9] = -2.38633905093136384013422325215527866148401465975954104585807;
	RK_B[12][10] = -3.25479342483643918654589367587788726747711504674780680269911;
	RK_B[12][11] = -2.16343541686422982353954211300054820889678036420109999154887;
	RK_B[13][0] = 0.895080295771632891049613132336585138148156279241561345991710;
	RK_B[13][1] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[13][2] = 0.197966831227192369068141770510388793370637287463360401555746;
	RK_B[13][3] = -0.0729547847313632629185146671595558023015011608914382961421311;
	RK_B[13][4] = 0.0000000000000000000000000000000000000000000000000000000000000;
	RK_B[13][5] = -0.851236239662007619739049371445966793289359722875702227166105;
	RK_B[13][6] = 0.398320112318533301719718614174373643336480918103773904231856;
	RK_B[13][7] = 3.63937263181035606029412920047090044132027387893977804176229;
	RK_B[13][8] = 1.54822877039830322365301663075174564919981736348973496313065;
	RK_B[13][9] = -2.12221714704053716026062427460427261025318461146260124401561;
	RK_B[13][10] = -1.58350398545326172713384349625753212757269188934434237975291;
	RK_B[13][11] = -1.71561608285936264922031819751349098912615880827551992973034;
	RK_B[13][12] = -0.0244036405750127452135415444412216875465593598370910566069132;
	RK_B[14][0] = -0.915176561375291440520015019275342154318951387664369720564660;
	RK_B[14][1] = 1.45453440217827322805250021715664459117622483736537873607016;
	RK_B[14][2] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[14][3] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[14][4] = -0.777333643644968233538931228575302137803351053629547286334469;
	RK_B[14][5] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[14][6] = -0.0910895662155176069593203555807484200111889091770101799647985;
	RK_B[14][7] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[14][8] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[14][9] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[14][10] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[14][11] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[14][12] = 0.0910895662155176069593203555807484200111889091770101799647985;
	RK_B[14][13] = 0.777333643644968233538931228575302137803351053629547286334469;
	RK_B[15][0] = 0.100000000000000000000000000000000000000000000000000000000000;
	RK_B[15][1] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[15][2] = -0.157178665799771163367058998273128921867183754126709419409654;
	RK_B[15][3] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[15][4] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[15][5] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[15][6] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[15][7] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[15][8] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[15][9] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[15][10] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[15][11] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[15][12] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[15][13] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[15][14] = 0.157178665799771163367058998273128921867183754126709419409654;
	RK_B[16][0] = 0.181781300700095283888472062582262379650443831463199521664945;
	RK_B[16][1] = 0.675000000000000000000000000000000000000000000000000000000000;
	RK_B[16][2] = 0.342758159847189839942220553413850871742338734703958919937260;
	RK_B[16][3] = 0.000000000000000000000000000000000000000000000000000000000000;
	RK_B[16][4] = 0.259111214548322744512977076191767379267783684543182428778156;
	RK_B[16][5] = -0.358278966717952089048961276721979397739750634673268802484271;
	RK_B[16][6] = -1.04594895940883306095050068756409905131588123172378489286080;
	RK_B[16][7] = 0.930327845415626983292300564432428777137601651182965794680397;
	RK_B[16][8] = 1.77950959431708102446142106794824453926275743243327790536000;
	RK_B[16][9] = 0.100000000000000000000000000000000000000000000000000000000000;
	RK_B[16][10] = -0.282547569539044081612477785222287276408489375976211189952877;
	RK_B[16][11] = -0.159327350119972549169261984373485859278031542127551931461821;
	RK_B[16][12] = -0.145515894647001510860991961081084111308650130578626404945571;
	RK_B[16][13] = -0.259111214548322744512977076191767379267783684543182428778156;
	RK_B[16][14] = -0.342758159847189839942220553413850871742338734703958919937260;
	RK_B[16][15] = -0.675000000000000000000000000000000000000000000000000000000000;
}

void rk8Step(DynFun dynFun, double tLow, double tUpp, double zLow[], double zUpp[], int nDim) {

	if (!flagInitRK) {
		initializeRungeKuttaConstants();
		flagInitRK = true;
	}

	double dt = tUpp - tLow;

	/// Allocate memory:
	double t[17];
	double** z = new double*[17];
	for (int i = 0; i < 17; i++) {
		z[i] = new double[nDim];
	}
	double** f = new double*[17];
	for (int i = 0; i < 17; i++) {
		f[i] = new double[nDim];
	}

	/// Populate time grid:
	for (int iStage = 0; iStage < 17; iStage++) {
		t[iStage] = tLow + dt * RK_A[iStage];
	}

	/// Initial State:
	int iStage = 0;
	for (int iDim = 0; iDim < nDim; iDim++) {
		z[iDim][iStage] = zLow[iDim];
	}

	/// Dynamics at initial point:
	dynFun(t[0], z[0], f[0]);

	/// March through each stage:
	double sum;
	for (int iStage = 1; iStage < 17; iStage++) {
		for (int iDim = 0; iDim < nDim; iDim++) {
			sum = 0.0;
			for (int j = 0; j < iStage; j++) {
				sum = sum + RK_B[iStage][j] * f[iStage - 1][iDim];
			}
			z[iStage][iDim] = z[iStage - 1][iDim] + dt * sum;
		}
		dynFun(t[iStage], z[iStage], f[iStage]);
	}

	/// Compute the final estimate:
	for (int iDim = 0; iDim < nDim; iDim++) {
		sum = 0.0;
		for (int iStage = 0; iStage < 17; iStage++) {
			sum = sum + RK_C[iStage] * f[iStage][iDim];
		}
		zUpp[iDim] = zLow[iDim] + dt * sum;
	}

	/// Release memory:
	for (int i = 0; i < 17; i++) {
		delete [] z[i];
	}
	delete [] z;
	for (int i = 0; i < 17; i++) {
		delete [] f[i];
	}
	delete [] f;
}

/* Runs several time steps using euler integration */
void simulate(DynFun dynFun, double t0, double t1, double z0[], double z1[],
              int nDim, int nStep, IntegrationMethod method)
{
	double dt, tLow, tUpp;
	double *zLow;
	double *zUpp;

	/// Allocate memory:
	zLow = new double[nDim];
	zUpp = new double[nDim];

	/// File IO stuff:
	ofstream logFile;
	logFile.open("logFile.csv");

	/// Initial conditions
	tLow = t0;
	for (int i = 0; i < nDim; i++) {
		zLow[i] = z0[i];
	}

	/// March forward in time:
	dt = (t1 - t0) / ((double) nStep);
	for (int i = 0; i < nStep; i++) {
		tUpp = tLow + dt;
		switch (method) {
		case Euler:
			eulerStep(dynFun, tLow, tUpp, zLow, zUpp, nDim); break;
		case MidPoint:
			midPointStep(dynFun, tLow, tUpp, zLow, zUpp, nDim); break;
		case RungeKutta:
			rungeKuttaStep(dynFun, tLow, tUpp, zLow, zUpp, nDim); break;
		case RK8:
			rk8Step(dynFun, tLow, tUpp, zLow, zUpp, nDim); break;
		}

		/// Print the state of the simulation:
		printState(logFile, tLow, zLow, nDim);

		/// Advance temp variables:
		tLow = tUpp;
		for (int j = 0; j < nDim; j++) {
			zLow[j] = zUpp[j];
		}
	}
	printState(logFile, tLow, zLow, nDim);

	delete [] zLow;
	delete [] zUpp;

	logFile.close();

}