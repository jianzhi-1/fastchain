#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <fstream>
#include <cstring>

#include "vector3d.h"

using namespace std;
#define ldd long double
typedef pair<double, double> pi;

double currenttime = 0.0;

const static int N = 1000;
FILE *fx, *fy, *ftheta, *fxdot, *fydot, *fthetadot;

double g = 9.81; //gravity value
int n = 3; //number of sticks
double deltat = 0.0005; //small time interval
double ddt = 0.0000001;
vector<pi> x[N], y[N], theta[N], xdot[N], ydot[N], thetadot[N]; //respective x y theta coordinates and dots //pair<time, xdot>
vector<double> cxdot[N], cydot[N], cthetadot[N]; //temporary storage of velocities
double h0 = 0.25; //distance from cg of bottom rod
double lA = 0.37; //long string
double la = 0.13; //short string
double dh = (lA + la)/2; //subsequent differences in change of height of the rods
double w = 0.24; //distance between 2 end points of a string on rod;
double theta0 = asin(((lA - la)/2)/w);
double epsilon = 0.00001; //a very small constant;

double ll = 0.3; //length of 2 rods
double m = 0.012; //mass of the rods
double inertia = (1.0/12.0)*m*ll*ll;  //inertia of the rods
double e = 0.85; //coefficient of restitution;

double hitfloor[50];
long double pii = 3.14159265359;


void fileoutput(int index){
	for (int i = 0; i < (int)x[index].size(); i++){
		fprintf(fx, "%lf %lf\n", x[index][i].first, x[index][i].second);
	}
	for (int i = 0; i < (int)y[index].size(); i++){
		fprintf(fy, "%lf %lf\n", y[index][i].first, y[index][i].second);
	}
	for (int i = 0; i < (int)theta[index].size(); i++){
		fprintf(ftheta, "%lf %lf\n", theta[index][i].first, theta[index][i].second);
	}
	for (int i = 0; i < (int)xdot[index].size(); i++){
		fprintf(fxdot, "%lf %lf\n", xdot[index][i].first, xdot[index][i].second);
	}
	for (int i = 0; i < (int)ydot[index].size(); i++){
		fprintf(fydot, "%lf %lf\n", ydot[index][i].first, ydot[index][i].second);
	}
	for (int i = 0; i < (int)thetadot[index].size(); i++){
		fprintf(fthetadot, "%lf %lf\n", thetadot[index][i].first, thetadot[index][i].second);
	}
}

 //Output Function
 //q == 1 if print out x, 2 if print y, 3 if print theta
 //q == 4 if print xdot, 5 if print ydot, 6 if print thetado
void printing(int q, int index){ 
	if (q == 1){
		for (int i = 0; i < (int)x[index].size(); i++){
			printf("%lf %lf\n", x[index][i].first, x[index][i].second);
		}
	} else if (q == 2){
		for (int i = 0; i < (int)y[index].size(); i++){
			printf("%lf %lf\n", y[index][i].first, y[index][i].second);
		}
	} else if (q == 3){
		for (int i = 0; i < (int)theta[index].size(); i++){
			printf("%lf %lf\n", theta[index][i].first, theta[index][i].second);
		}
	} else if (q == 4){
		for (int i = 0; i < (int)xdot[index].size(); i++){
			printf("%lf %lf\n", xdot[index][i].first, xdot[index][i].second);
		}
	} else if (q == 5){
		for (int i = 0; i < (int)ydot[index].size(); i++){
			printf("%lf %lf\n", ydot[index][i].first, ydot[index][i].second);
		}
	} else if (q == 6){
		for (int i = 0; i < (int)thetadot[index].size(); i++){
			printf("%lf %lf\n", thetadot[index][i].first, thetadot[index][i].second);
		}
	} else {
		printf("Wrong command\n");
	}
	
}

//retuens x to the power of y     OK
double power(double x, int y){
	double xx = 1.0;
	for (int i = 0; i < y; i++){
		xx *= x;
	}
	return xx;
}

//returns the square of distance between point 1 and point 2
//OK
double dist(pi p1, pi p2){ 
	return power(p1.first - p2.first, 2) + power(p1.second - p2.second, 2);
}

//returns the absolute value of x OK
double abso(double x){
	if (x >= 0.0000) return x;
	return -x;
}

//returns x to the range [-pi, pi] OK
double modpi(double x){
	while (x >= pii){
		x -= 2.0*pii;
	}
	
	while (x <= -pii){
		x += 2.0*pii;
	}
	return x;
}

//returns the coordinate of the point at w/2 from the center of rod
//index is the index of the rod, b is 0 if looking at left side, otherwise right
//OK
pi gcw(int index, int b){ 
	if (b == 0){
		return make_pair(x[index].back().second - (w/2.0)*cos(theta[index].back().second), y[index].back().second - (w/2.0)*sin(theta[index].back().second));
	} else {
		return make_pair(x[index].back().second + (w/2.0)*cos(theta[index].back().second), y[index].back().second + (w/2.0)*sin(theta[index].back().second));
	}
}

//returns the coordinate of the point at l/2 from the center of the rod
//index is the index of the rod, b is 0 if looking at the left side, otherwise right, used for the l
//OK
pi gcl(int index, int b){
	if (b == 0){
		return make_pair(x[index].back().second - (ll/2.0)*cos(theta[index].back().second), y[index].back().second - (ll/2.0)*sin(theta[index].back().second));
	} else {
		return make_pair(x[index].back().second + (ll/2.0)*cos(theta[index].back().second), y[index].back().second + (ll/2.0)*sin(theta[index].back().second));
	}
}



//setting up all the base variables OK
void init(){  
	currenttime = 0.0;
	for (int i = 0; i < n; i++){
		x[i].push_back(make_pair(currenttime, 0.0));
		y[i].push_back(make_pair(currenttime, h0 + i*dh));
		if (i % 2){
			theta[i].push_back(make_pair(currenttime, -theta0));
		} else {
			theta[i].push_back(make_pair(currenttime, theta0));
		}
		xdot[i].push_back(make_pair(currenttime, 0.0));
		ydot[i].push_back(make_pair(currenttime, 0.0));
		thetadot[i].push_back(make_pair(currenttime, 0.0));
	}
}

//simulates what will happen normally if dt passes for all the sticks above and including index
//OK
void updatenormal(int index){ 
	currenttime += deltat; //increment by deltat
	for (int i = index; i < n; i++){
		ydot[i].push_back(make_pair(currenttime, ydot[i].back().second - deltat*g));
		y[i].push_back(make_pair(currenttime, y[i].back().second + deltat*ydot[i].back().second));
		x[i].push_back(make_pair(currenttime, x[i].back().second + deltat*xdot[i].back().second));
		theta[i].push_back(make_pair(currenttime, modpi(theta[i].back().second + deltat*thetadot[i].back().second)));
	}
}



//returns the vector of the string between rods indexed index and index + 1
//b is 0 if wish to get left string, otherwise right
//OK
v3 sv(int index, int b){ 
	if (index >= n || index < 0){
		printf("Wrong parameters given to function sv!\n");
		return v3(0.0, 0.0, 0.0);
	}
	pi rodi = gcw(index, b);
	pi rodii = gcw(index + 1, b);
	return v3(rodii.first - rodi.first, rodii.second - rodi.second, 0);
}

//Computing function when a string is tight
//for the left string
void computell(int index){
	//printf("Left string is tight\n");
	if (index == n - 1){
		printf("Computing left string of invalid index %d\n", index);
		return;
	}
	v3 strvector = sv(index, 0);
	double theta1 = acos(dot(strvector, v3(1.0, 0.0, 0.0))/magnitude(strvector));
	double wp = w/2.0;
	
	double xdot1 = xdot[index].back().second;
	double xdot2 = xdot[index + 1].back().second;
	double ydot1 = ydot[index].back().second;
	double ydot2 = ydot[index + 1].back().second;
	double thetadot1 = thetadot[index].back().second;
	double thetadot2 = thetadot[index + 1].back().second;
	
	double phi1 = theta[index].back().second;
	double phi2 = theta[index + 1].back().second;
	
	//double xfdot1 = (-inertia*x1dot*power(cos(phi1 + theta1), 2) + 2.0*inertia*xdot2*cos(phi1 + theta1)*cos(phi2 - theta2) + inertia*xdot1*power(cos(phi2 - theta2, 2)) + 2*inertia*thetadot1*wp*cos(phi1 + theta1)*)
	//double xfdot1 = (-inertia*xdot1*power(cos(phi1 + theta1), 2) + 2.0*inertia*xdot2*cos(phi1 + theta1)*cos(phi2 - theta2) + inertia*xdot1*power(cos(phi2 - theta2), 2) + 2.0*inertia*thetadot1*wp*cos(phi1 + theta1)*sin(theta1) + m*power(wp, 2)*xdot1*power(sin(theta1), 2) - 2.0*inertia*ydot1*cos(phi1 + theta1)*sin(phi1 + theta1) + inertia*xdot1*power(sin(phi1 + theta1), 2) - 2.0*inertia*ydot2*cos(phi1 + theta1)*sin(phi2 - theta2) + inertia*xdot1*power(sin(phi2 - theta2), 2) - 2.0*inertia*thetadot2*wp*cos(phi1 + theta1)*sin(theta2) + m*power(wp, 2)*xdot1*power(sin(theta2), 2))/(inertia*power(cos(phi1 + theta1), 2) + inertia*power(cos(phi2 - theta2), 2) + m*power(wp, 2)*power(sin(theta1), 2) + inertia*power(sin(phi1 + theta1), 2) + inertia*power(sin(phi2 - theta2), 2) + m*power(wp, 2)*power(sin(theta2), 2));
	double xfdot1 = (2.0*inertia*xdot2*power(cos(theta1), 2) + 2.0*inertia*thetadot2*wp*cos(theta1)*sin(phi2 - theta1) + m*xdot1*power(wp*sin(phi2 - theta1), 2) - 2.0*inertia*thetadot1*wp*cos(theta1)*sin(phi1 - theta1) + m*power(wp*sin(phi1 - theta1), 2)*xdot1 - 2.0*inertia*ydot1*cos(theta1)*sin(theta1) + 2.0*inertia*ydot2*cos(theta1)*sin(theta1) + 2.0*inertia*xdot1*power(sin(theta1), 2))/(2.0*inertia*power(cos(theta1), 2) + m*power(wp*sin(phi2 - theta1), 2) + m*power(wp*sin(phi1 - theta1), 2) + 2.0*inertia*power(sin(theta1), 2));
	//double yfdot1 = (2.0*inertia*ydot1*power(cos(theta1), 2) + m*power(wp*sin(phi2 - theta1), 2)*ydot1 + m*power(wp*sin(phi1 - theta1), 2)*ydot1 - 2.0*inertia*xdot1*cos(theta1)*sin(theta1) + 2.0*inertia*xdot2*cos(theta1)*sin(theta1) + 2.0*inertia*thetadot2*wp*sin(phi2 - theta1)*sin(theta1) - 2.0*inertia*thetadot1*wp*sin(phi1 - theta1)*sin(theta1) + 2.0*inertia*ydot2*power(sin(theta1), 2))/(2.0*inertia*power(cos(theta1), 2) + m*power(wp*sin(phi2 - theta1), 2) + m*power(wp*sin(phi1 - theta1), 2) + 2.0*inertia*power(sin(theta1), 2));
	double yfdot1 = (2*inertia*ydot1*power(cos(theta1),2) + m*power(wp,2)*ydot1*power(sin(phi2 - theta1),2) + m*power(wp,2)*ydot1*power(sin(phi1 - theta1),2) - 2*inertia*xdot1*cos(theta1)*sin(theta1) + 2*inertia*xdot2*cos(theta1)*sin(theta1) + 2*inertia*thetadot2*wp*sin(phi2 - theta1)*sin(theta1) - 2*inertia*thetadot1*wp*sin(phi1 - theta1)*sin(theta1) + 2*inertia*ydot2*power(sin(theta1),2))/(2*inertia*power(cos(theta1),2) + m*power(wp,2)*power(sin(phi2 - theta1),2) + m*power(wp,2)*power(sin(phi1 - theta1),2) + 2*inertia*power(sin(theta1),2));
	double xfdot2 = (2.0*inertia*xdot1*power(cos(theta1), 2) - 2.0*inertia*thetadot2*wp*cos(theta1)*sin(phi2 - theta1) + m*power(wp*sin(phi2 - theta1), 2)*xdot2 + 2.0*inertia*thetadot1*wp*cos(theta1)*sin(phi1 - theta1) + m*power(wp, 2)*xdot2*power(sin(phi1 - theta1), 2) + 2.0*inertia*ydot1*cos(theta1)*sin(theta1) - 2.0*inertia*ydot2*cos(theta1)*sin(theta1) + 2.0*inertia*xdot2*power(sin(theta1), 2))/(2.0*inertia*power(cos(theta1), 2) + m*power(wp*sin(phi2 - theta1), 2) + m*power(wp*sin(phi1 - theta1), 2) + 2.0*inertia*power(sin(theta1), 2));
   	double yfdot2 = (2*inertia*ydot2*power(cos(theta1),2) + m*power(wp,2)*ydot2*power(sin(phi2 - theta1),2) + m*power(wp,2)*ydot2*power(sin(phi1 - theta1),2) + 2*inertia*xdot1*cos(theta1)*sin(theta1) - 2*inertia*xdot2*cos(theta1)*sin(theta1) - 2*inertia*thetadot2*wp*sin(phi2 - theta1)*sin(theta1) + 2*inertia*thetadot1*wp*sin(phi1 - theta1)*sin(theta1) + 2*inertia*ydot1*power(sin(theta1),2))/(2*inertia*power(cos(theta1),2) + m*power(wp,2)*power(sin(phi2 - theta1),2) + m*power(wp,2)*power(sin(phi1 - theta1),2) + 2*inertia*power(sin(theta1),2));
   	double thetafdot1 = (2*inertia*thetadot1*power(cos(theta1),2) + m*thetadot1*power(wp,2)*power(sin(phi2 - theta1),2) - 2*m*wp*xdot1*cos(theta1)*sin(phi1 - theta1) + 2*m*wp*xdot2*cos(theta1)*sin(phi1 - theta1) + 2*m*thetadot2*power(wp,2)*sin(phi2 - theta1)*sin(phi1 - theta1) - m*thetadot1*power(wp,2)*power(sin(phi1 - theta1),2) - 2*m*wp*ydot1*sin(phi1 - theta1)*sin(theta1) + 2*m*wp*ydot2*sin(phi1 - theta1)*sin(theta1) + 2*inertia*thetadot1*power(sin(theta1),2))/(2*inertia*power(cos(theta1),2) + m*power(wp,2)*power(sin(phi2 - theta1),2) + m*power(wp,2)*power(sin(phi1 - theta1),2) + 2*inertia*power(sin(theta1),2));
   	double thetafdot2 = (2*inertia*thetadot2*power(cos(theta1),2) + 2*m*wp*xdot1*cos(theta1)*sin(phi2 - theta1) - 2*m*wp*xdot2*cos(theta1)*sin(phi2 - theta1) - m*thetadot2*power(wp,2)*power(sin(phi2 - theta1),2) + 2*m*thetadot1*power(wp,2)*sin(phi2 - theta1)*sin(phi1 - theta1) + m*thetadot2*power(wp,2)*power(sin(phi1 - theta1),2) + 2*m*wp*ydot1*sin(phi2 - theta1)*sin(theta1) - 2*m*wp*ydot2*sin(phi2 - theta1)*sin(theta1) + 2*inertia*thetadot2*power(sin(theta1),2))/(2*inertia*power(cos(theta1),2) + m*power(wp,2)*power(sin(phi2 - theta1),2) + m*power(wp,2)*power(sin(phi1 - theta1),2) + 2*inertia*power(sin(theta1),2));
	
	cxdot[index].push_back(xfdot1);
	cydot[index].push_back(yfdot1);
	cthetadot[index].push_back(thetafdot1);
	cxdot[index + 1].push_back(xfdot2);
	cydot[index + 1].push_back(yfdot2);
	cthetadot[index + 1].push_back(thetafdot2);
}

//Computing function when a string is tight
//for the right string
void computerr(int index){
	//printf("Right string is tight\n");
	if (index == n - 1){
		printf("Computing right string of invalid index %d\n", index);
		return;
	}
	v3 strvector = sv(index, 1);
	double theta1 = acos(dot(strvector, v3(1.0, 0.0, 0.0))/magnitude(strvector));
	double wp = w/2.0;
	
	double xdot1 = xdot[index].back().second;
	double xdot2 = xdot[index + 1].back().second;
	double ydot1 = ydot[index].back().second;
	double ydot2 = ydot[index + 1].back().second;
	double thetadot1 = thetadot[index].back().second;
	double thetadot2 = thetadot[index + 1].back().second;
	
	double phi1 = theta[index].back().second;
	double phi2 = theta[index + 1].back().second;
	
	double xfdot1 = (2*inertia*xdot2*power(cos(theta1),2) - 2*inertia*thetadot2*wp*cos(theta1)*sin(phi2 - theta1) + m*power(wp,2)*xdot1*power(sin(phi2 - theta1),2) + 2*inertia*thetadot1*wp*cos(theta1)*sin(phi1 - theta1) + m*power(wp,2)*xdot1*power(sin(phi1 - theta1),2) - 2*inertia*ydot1*cos(theta1)*sin(theta1) + 2*inertia*ydot2*cos(theta1)*sin(theta1) + 2*inertia*xdot1*power(sin(theta1),2))/(2*inertia*power(cos(theta1),2) + m*power(wp,2)*power(sin(phi2 - theta1),2) + m*power(wp,2)*power(sin(phi1 - theta1),2) + 2*inertia*power(sin(theta1),2));
    double yfdot1 = (2*inertia*ydot1*power(cos(theta1),2) + m*power(wp,2)*ydot1*power(sin(phi2 - theta1),2) + m*power(wp,2)*ydot1*power(sin(phi1 - theta1),2) - 2*inertia*xdot1*cos(theta1)*sin(theta1) + 2*inertia*xdot2*cos(theta1)*sin(theta1) - 2*inertia*thetadot2*wp*sin(phi2 - theta1)*sin(theta1) + 2*inertia*thetadot1*wp*sin(phi1 - theta1)*sin(theta1) + 2*inertia*ydot2*power(sin(theta1),2))/(2*inertia*power(cos(theta1),2) + m*power(wp,2)*power(sin(phi2 - theta1),2) + m*power(wp,2)*power(sin(phi1 - theta1),2) + 2*inertia*power(sin(theta1),2));
    double xfdot2 = (2*inertia*xdot1*power(cos(theta1),2) + 2*inertia*thetadot2*wp*cos(theta1)*sin(phi2 - theta1) + m*power(wp,2)*xdot2*power(sin(phi2 - theta1),2) - 2*inertia*thetadot1*wp*cos(theta1)*sin(phi1 - theta1) + m*power(wp,2)*xdot2*power(sin(phi1 - theta1),2) + 2*inertia*ydot1*cos(theta1)*sin(theta1) - 2*inertia*ydot2*cos(theta1)*sin(theta1) + 2*inertia*xdot2*power(sin(theta1),2))/(2*inertia*power(cos(theta1),2) + m*power(wp,2)*power(sin(phi2 - theta1),2) + m*power(wp,2)*power(sin(phi1 - theta1),2) + 2*inertia*power(sin(theta1),2));
    double yfdot2 = (2*inertia*ydot2*power(cos(theta1),2) + m*power(wp,2)*ydot2*power(sin(phi2 - theta1),2) + m*power(wp,2)*ydot2*power(sin(phi1 - theta1),2) + 2*inertia*xdot1*cos(theta1)*sin(theta1) - 2*inertia*xdot2*cos(theta1)*sin(theta1) + 2*inertia*thetadot2*wp*sin(phi2 - theta1)*sin(theta1) - 2*inertia*thetadot1*wp*sin(phi1 - theta1)*sin(theta1) + 2*inertia*ydot1*power(sin(theta1),2))/(2*inertia*power(cos(theta1),2) + m*power(wp,2)*power(sin(phi2 - theta1),2) + m*power(wp,2)*power(sin(phi1 - theta1),2) + 2*inertia*power(sin(theta1),2));
    double thetafdot1 = (2*inertia*thetadot1*power(cos(theta1),2) + m*thetadot1*power(wp,2)*power(sin(phi2 - theta1),2) + 2*m*wp*xdot1*cos(theta1)*sin(phi1 - theta1) - 2*m*wp*xdot2*cos(theta1)*sin(phi1 - theta1) + 2*m*thetadot2*power(wp,2)*sin(phi2 - theta1)*sin(phi1 - theta1) - m*thetadot1*power(wp,2)*power(sin(phi1 - theta1),2) + 2*m*wp*ydot1*sin(phi1 - theta1)*sin(theta1) - 2*m*wp*ydot2*sin(phi1 - theta1)*sin(theta1) + 2*inertia*thetadot1*power(sin(theta1),2))/(2*inertia*power(cos(theta1),2) + m*power(wp,2)*power(sin(phi2 - theta1),2) + m*power(wp,2)*power(sin(phi1 - theta1),2) + 2*inertia*power(sin(theta1),2));
    double thetafdot2 = (2*inertia*thetadot2*power(cos(theta1),2) - 2*m*wp*xdot1*cos(theta1)*sin(phi2 - theta1) + 2*m*wp*xdot2*cos(theta1)*sin(phi2 - theta1) - m*thetadot2*power(wp,2)*power(sin(phi2 - theta1),2) + 2*m*thetadot1*power(wp,2)*sin(phi2 - theta1)*sin(phi1 - theta1) + m*thetadot2*power(wp,2)*power(sin(phi1 - theta1),2) - 2*m*wp*ydot1*sin(phi2 - theta1)*sin(theta1) + 2*m*wp*ydot2*sin(phi2 - theta1)*sin(theta1) + 2*inertia*thetadot2*power(sin(theta1),2))/(2*inertia*power(cos(theta1),2) + m*power(wp,2)*power(sin(phi2 - theta1),2) + m*power(wp,2)*power(sin(phi1 - theta1),2) + 2*inertia*power(sin(theta1),2));

	cxdot[index].push_back(xfdot1);
	cydot[index].push_back(yfdot1);
	cthetadot[index].push_back(thetafdot1);
	cxdot[index + 1].push_back(xfdot2);
	cydot[index + 1].push_back(yfdot2);
	cthetadot[index + 1].push_back(thetafdot2);
}

//returns the radius vector of rod indexed index
//b = 0 if radius towards left end, b = 1 if radius towards right end
//OK
v3 radii(int index, int b){
	if (b == 0){  //radius towards left
		return v3(- (w/2.0)*cos(theta[index].back().second), - (w/2.0)*sin(theta[index].back().second), 0);
	} else {
		return v3( (w/2.0)*cos(theta[index].back().second),  (w/2.0)*sin(theta[index].back().second), 0);
	}
}


double velocityalongleft(int index, v3 str){
	v3 ori = v3(xdot[index].back().second, ydot[index].back().second, 0.0);
	//v3 angular = cross(v3(0, 0, thetadot[index].back()), v3(gcw(index, 0).first, gcw(index, 0).second, 0));
	v3 angular = cross(v3(0, 0, thetadot[index].back().second), radii(index, 0));
	v3 resultant = vadd(ori, angular);
	double valong = dot(resultant, str)/magnitude(str);
	return valong;
}

double velocityalongright(int index, v3 str){
	v3 ori = v3(xdot[index].back().second, ydot[index].back().second, 0);
	//v3 angular = cross(v3(0, 0, thetadot[index].back()), v3(gcw(index, 1).first, gcw(index, 1).second, 0));
	v3 angular = cross(v3(0, 0, thetadot[index].back().second), radii(index, 1));
	v3 resultant = vadd(ori, angular);
	double valong = dot(resultant, str)/magnitude(str);
	return valong;
}

//checks if any strings are right and performs relevant operations
//index is the index of the rod
//OK
//returns true if a string is tight
bool checkstring(int index){ 
	bool tighten = false;
	//check left string
	for (int i = index; i < n - 1; i++){
		double str; //stores the string length
		if (i % 2){str = la;} else {str = lA;}
		if (dist(gcw(i, 0), gcw(i + 1, 0)) > power(str, 2) + epsilon){ //if true then string is tight, otherwise string is not tight
			//printf("Left String is Tight\n");
			tighten = true;
			//if (velocityalongleft(i + 1, sv(i, 0)) < velocityalongleft(i, sv(i, 0))) continue; //no problem
			computell(i);
		}
	}
	
	//check right string
	for (int i = index; i < n - 1; i++){
		double str; //stores the string length
		if (i % 2){str = lA;} else {str = la;}
		if (dist(gcw(i, 1), gcw(i + 1, 1)) > power(str, 2) + epsilon){ //if true then string is tight
			//printf("Right String is Tight\n");
			tighten = true;
			//if (velocityalongright(i + 1, sv(i, 1)) < velocityalongright(i, sv(i, 1))) continue;
			computerr(i);
		}
	}
	return tighten;
}

void checkstring2(int index){ 

	double str; //stores the string length
	if (index % 2){str = la;} else {str = lA;}
	if (dist(gcw(index, 0), gcw(index + 1, 0)) > power(str, 2) + epsilon){ //if true then string is tight, otherwise string is not tight
		//printf("Left String is Tight\n");
		//if (velocityalongleft(i + 1, sv(i, 0)) < velocityalongleft(i, sv(i, 0))) continue; //no problem
		computell(index);
	}
	
	if (index % 2){str = lA;} else {str = la;}
	if (dist(gcw(index, 1), gcw(index + 1, 1)) > power(str, 2) + epsilon){ //if true then string is tight
		//printf("Right String is Tight\n");
		//if (velocityalongright(i + 1, sv(i, 1)) < velocityalongright(i, sv(i, 1))) continue;
		computerr(index);
	}
}

//updates the velocity after multiple strings have become tight
//only updates rods in the range [index, n)
//OK
void updatevelocity(int index){ 
	currenttime += ddt;
	for (int i = index; i < n; i++){
		//update the x velocity
		double totall = 0.0; //totall stores all the deltas, addition of all impulses
		while (cxdot[i].size() > 0){
			totall += (cxdot[i].back() - xdot[i].back().second);
			cxdot[i].pop_back();
		}
		xdot[i].push_back(make_pair(currenttime, totall + xdot[i].back().second));
		cxdot[i].clear();
		
		totall = 0.0;
		while (cydot[i].size() > 0){
			totall += (cydot[i].back() - ydot[i].back().second);
			cydot[i].pop_back();
		}
		ydot[i].push_back(make_pair(currenttime, totall + ydot[i].back().second));
		cydot[i].clear();
		
		totall = 0.0;
		while (cthetadot[i].size() > 0){
			totall += (cthetadot[i].back() - thetadot[i].back().second);
			cthetadot[i].pop_back();
		}
		thetadot[i].push_back(make_pair(currenttime, totall + thetadot[i].back().second));
		cthetadot[i].clear();
	}
}

//Update function for when the rod hits the ground with left end first
void hitgroundupdateleft(int index){ 
	currenttime += ddt;
	//printf("Rod hit ground on the left end\n");
	//double newthetadot = (ll*thetadot[index].back() + 6.0*(1 + e)*ydot[index].back()*cos(theta[index].back()) - 3.0*e*ll*thetadot[index].back()*power(cos(theta[index].back()), 2))/((1.0 + 3.0*power(cos(theta[index].back()), 2))*ll);
	//double newydot = (6.0*ydot[index].back()*power(cos(theta[index].back()), 2) + ll*(1.0 + e)*thetadot[index].back()*cos(theta[index].back()) - 2.0*e*ydot[index].back())/(2.0*(1.0 + 3.0*power(cos(theta[index].back()), 2)));
	double newthetadot = (ll*thetadot[index].back().second + 6.0*(1 + e)*ydot[index].back().second*cos(theta[index].back().second) - 3.0*e*ll*thetadot[index].back().second*power(cos(theta[index].back().second), 2))/((1.0 + 3.0*power(cos(theta[index].back().second), 2))*ll);
	double newydot = (6.0*ydot[index].back().second*power(cos(theta[index].back().second), 2) + ll*(1.0 + e)*thetadot[index].back().second*cos(theta[index].back().second) - 2.0*e*ydot[index].back().second)/(2.0*(1.0 + 3.0*power(cos(theta[index].back().second), 2)));	
	ydot[index].push_back(make_pair(currenttime, newydot));
	thetadot[index].push_back(make_pair(currenttime, newthetadot));
}

//Update function for when the rod hits the ground with right end first
void hitgroundupdateright(int index){
	currenttime += ddt;
	//printf("Rod hit ground on the right end\n");
	//double newthetadot = (ll*thetadot[index].back() - 6.0*(1 + e)*ydot[index].back()*cos(theta[index].back()) - 3.0*e*ll*thetadot[index].back()*power(cos(theta[index].back()), 2))/((1.0 + 3.0*power(cos(theta[index].back()), 2))*ll);
	//double newydot = (6.0*ydot[index].back()*power(cos(theta[index].back()), 2) - ll*(1.0 + e)*thetadot[index].back()*cos(theta[index].back()) - 2.0*e*ydot[index].back())/(2.0*(1.0 + 3.0*power(cos(theta[index].back()), 2)));
	double newthetadot = (ll*thetadot[index].back().second - 6.0*(1 + e)*ydot[index].back().second*cos(theta[index].back().second) - 3.0*e*ll*thetadot[index].back().second*power(cos(theta[index].back().second), 2))/((1.0 + 3.0*power(cos(theta[index].back().second), 2))*ll);
	double newydot = (6.0*ydot[index].back().second*power(cos(theta[index].back().second), 2) - ll*(1.0 + e)*thetadot[index].back().second*cos(theta[index].back().second) - 2.0*e*ydot[index].back().second)/(2.0*(1.0 + 3.0*power(cos(theta[index].back().second), 2)));
	ydot[index].push_back(make_pair(currenttime, newydot));
	thetadot[index].push_back(make_pair(currenttime, newthetadot));	
}

//clear all the arrays    OK
void clearall(){
	for (int i = 0; i < n; i++){
		x[i].clear();
		y[i].clear();
		theta[i].clear();
		xdot[i].clear();
		ydot[i].clear();
		thetadot[i].clear();
	}
}

string tostring(int x){
	if (x == 0) return "0";
	string a = "";
	while (x != 0){
		a = (char)((x%10) + 48) + a;
		x -= (x%10);
		x /= 10;
	}
	return a;
}

int main(){
	FILE *fx[N], *fy[N], *ftheta[N], *fxdot[N], *fydot[N], *fthetadot[N];
	FILE *hitf;
	for (int i = 0; i < n; i++){
		fx[i] = fopen(("x[" + tostring(i) + "].txt").c_str(), "w");
		fy[i] = fopen(("y[" + tostring(i) + "].txt").c_str(), "w");
		ftheta[i] = fopen(("theta[" + tostring(i) + "].txt").c_str(), "w");
		fxdot[i] = fopen(("xdot[" + tostring(i) + "].txt").c_str(), "w");
		fydot[i] = fopen(("ydot[" + tostring(i) + "].txt").c_str(), "w");
		fthetadot[i] = fopen(("thetadot[" + tostring(i) + "].txt").c_str(), "w");
	}
	hitf = fopen("hitfloor.txt", "w");

	
	clearall();
	init();
	int bottom = 0;   //now considering the bottom most rod first;
	while (bottom < n){
		
		while (y[bottom].back().second >= epsilon + (ll/2)*abso(sin(theta[bottom].back().second))){
			updatenormal(bottom);
			for (int i = bottom; i < n - 1; i++){
				checkstring2(i);
				updatevelocity(i);
				updatenormal(bottom);
			}
			//checkstring(bottom);
			//updatevelocity(bottom);
		}
		
		//now the bottom most rod has hit the floor
		if (bottom == 0){
			//printf("%lf\n", currenttime);
		}
		hitfloor[bottom] = currenttime;
		fprintf(hitf, "%d %lf\n", bottom, currenttime);
		if (theta[bottom].back().second > 0.0){  //means left end touching floor
			hitgroundupdateleft(bottom);
			updatenormal(bottom);
			for (int i = bottom; i < n - 1; i++){
				checkstring2(i);
				updatevelocity(i);
				updatenormal(bottom);
			}
			while (y[bottom].back().second >= epsilon + (ll/2.0)*abso(sin(theta[bottom].back().second))){
				if (bottom != n - 1 && y[bottom].back().second > y[bottom + 1].back().second){
					break;
				}
				for (int i = bottom; i < n - 1; i++){
					checkstring2(i);
					updatevelocity(i);
					updatenormal(bottom);
				}
				updatenormal(bottom);
				//checkstring(bottom);
				//updatevelocity(bottom);
			}
		} else {
			hitgroundupdateright(bottom);
			updatenormal(bottom);
			for (int i = bottom; i < n - 1; i++){
				checkstring2(i);
				updatevelocity(i);
				updatenormal(bottom);
			}
			while (y[bottom].back().second >= epsilon + (ll/2.0)*abso(sin(theta[bottom].back().second))){
				if (bottom != n - 1 && y[bottom].back().second > y[bottom + 1].back().second){
					break;
				}
				for (int i = bottom; i < n - 1; i++){
					checkstring2(i);
					updatevelocity(i);
					updatenormal(bottom);
				}
				updatenormal(bottom);
				//checkstring(bottom);
				//updatevelocity(bottom);
			}
		}
		//printf("The %dth stick is done\n", bottom + 1);
		bottom++; //move onto the next stick
		
	}
	for (int j = 0; j < n; j++){
		int indian = j;
		for (int i = 0; i < (int)x[indian].size(); i++){
			fprintf(fx[j], "%lf %lf\n", x[indian][i].first, x[indian][i].second);
		}
		for (int i = 0; i < (int)y[indian].size(); i++){
			fprintf(fy[j], "%lf %lf\n", y[indian][i].first, y[indian][i].second);
		}
		for (int i = 0; i < (int)theta[indian].size(); i++){
			fprintf(ftheta[j], "%lf %lf\n", theta[indian][i].first, theta[indian][i].second);
		}
		for (int i = 0; i < (int)xdot[indian].size(); i++){
			fprintf(fxdot[j], "%lf %lf\n", xdot[indian][i].first, xdot[indian][i].second);
		}
		for (int i = 0; i < (int)ydot[indian].size(); i++){
			fprintf(fydot[j], "%lf %lf\n", ydot[indian][i].first, ydot[indian][i].second);
		}
		for (int i = 0; i < (int)thetadot[indian].size(); i++){
			fprintf(fthetadot[j], "%lf %lf\n", thetadot[indian][i].first, thetadot[indian][i].second);
		}
	}
	
	fclose(hitf);
	for (int i = 0; i < n; i++){
		fclose(fx[i]);
		fclose(fy[i]);
		fclose(ftheta[i]);
		fclose(fxdot[i]);
		fclose(fydot[i]);
		fclose(fthetadot[i]);	
	}
	//do whatever command you need
	//printing(5, n - 1);
	//printing(2, n - 1);
	//printf("%lf %lf\n", h0, sqrt(y[n - 1][0]*2.0/g) - (y[n - 1].size()*deltat));

	
}
