/************
/ Programme that checks proposed enhanced Bianchi's model
************************/
#include <math.h>
#include <iostream>
#include <cmath>
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <vector>


using namespace std;

const double Wo=16;    // the minimal size of the CW
const double mDASH=0;  // The maximal size of CW is 2^(mDASH)*W 
const double m=7; 	// or 4 to check The maximum  retransmission number

const double ACK=14*8;	// ACK size 14 Bytes
const double RH = 6e6; // Header rates (valid for PCLP; MAC and ACK)

const double TACK=ACK/RH; // Duration of an acknowlege

const double Tslot=9*pow(10,-6);// Duration of an emty slot

const double TH= (28*8)/(6e6)+20*pow(10,-6);//Duration of a PHY head + a MAC head

const double PS=1448*8; // in Bites to express the MAC Payload size.

const double delta=1*pow(10,-6);//Propagation delay, 1µs

const double DIFS= 28*pow(10,-6);// Distributed Inter Frame Space
const double SIFS= 10*pow(10,-6);// Short Inter Frame Space

const double BER=pow(10,-5);

class Ju_model{
public:
  //double BER;  		// The BER
  double Lp;   //size of the frame including PHY and MAC headers.
  double taux; // transmission probabillity
  double tauxR; // upper bound value for taux
  double tauxL; // lower bound value for taux
  double tauxM;	// The mean value
  
  double diff;
  double eps;   // The acceptable error
  		
  
  double n6;  // nbre of users operating on 11Mbps
  double n9;  // nbre of users operating on 5.5Mbps
  double n12; // nbre of users operating on 2Mbps
  double n18; // nbre of users operating on 1Mbps
  double n24;
  double n36;
  double n48; 
  double n54;
  double n;         // total nbre of users whithin the coverage area of an AP

  vector<double> N; // vector of type double contains the Clients number according to the Data Rates
  double pc;         // collision probability
  double pe;        /* erroneous received packet probability
		     pe expresses the the reception quality exp when BER = 10^-5 
		     pe equals to 0.08*/
  vector<double> Ptr_x;		// vector containing Probability that there is at least one transmission for each set
  double Ptr;			// value of the Probability that there is at least one transmission
  vector<double> Ps_x;		// probability that the transmission is succesful 
  
  vector<double> T_x; //
  vector<double> Tr_x;	// duration of unsuccesful transmission leading to retransmission
  vector<double> Ts_x; // duration whithin it the channel is occupied by a successful transmission 

  vector<double> x;   // vector of different data rates 1, 2 , 5.5 and 11 Mbps
  
  double pr;        // retransmission probability (pe+pc)
  vector<double> S_x; // vector of the different throughput according to the corresponding data rates
  double S; 		// the total throughput
  double S6;		// the throughput of the set of stations operating on 6Mbps 
  double S9;		// the throughputof the set of stations operating on 9Mbps
  double S12;		// the throughputof the set of stations operating on 12Mbps
  double S18;		// the throughputof the set of stations operating on 18Mbps
  double S24;		// the throughputof the set of stations operating on 24Mbps
  double S36;		// the throughputof the set of stations operating on 36Mbps
  double S48;		// the throughputof the set of stations operating on 48Mbps
  double S54;		// the throughputof the set of stations operating on 54Mbps
  
  
  void Print();
  void Initialize();  // Initialize some paprameters to do computation  
  void Compute_taux_collisionPrb(); // Compute 'taux' and 'commision probability'
  void Compute_throughput();
  double Compute_sumation();
};



/*Ju_model::Ju_model(){ 
 * 
}

Ju_model::~Ju_model(){
}*/

void Ju_model::Print(){

/*
cout<<" The total number of uers   "<< n6 + n9 + n12 + n18 + n24 + n36 + n48 + n54 << endl;
cout <<"  The number of users operating on 6Mbps:  "<< n6 <<endl;

cout <<"  The number of users operating on 9Mbps:  "<< n9 <<endl;

cout <<"  The number of users operating on 12Mbps:  "<< n12 <<endl;

cout <<"  The number of users operating on 18Mbps:  "<< n18 <<endl;

cout <<"  The number of users operating on 24Mbps:  "<< n24 <<endl;


cout <<"  The number of users operating on 36Mbps:  "<< n36 <<endl;

cout <<"  The number of users operating on 48Mbps:  "<< n48 <<endl;

cout <<"  The number of users operating on 54Mbps:  "<< n54 <<endl;
*/


N.push_back(n6);
N.push_back(n9);
N.push_back(n12);
N.push_back(n18);
N.push_back(n24);
N.push_back(n36);
N.push_back(n48);
N.push_back(n54);

}

void Ju_model::Initialize()
{
 taux = 0.0;
 tauxL = 0.0;
 tauxR = 1.2;
 tauxM = 0.0001;
 
 eps=1e-4;
 Lp = 100; // to see signification later
 pe=0.24;//0.2; // Initially we consider an ideal channel (without error probability) 
  
}
void Ju_model::Compute_taux_collisionPrb(){
	n=n6 + n9 + n12 + n18 + n24 + n36 + n48 + n54; // the sum of all users
	
bool cont = true;
    while (cont){
  	tauxM=(tauxR+tauxL)/2; 
	//cout <<tauxM << " The tauxM value"<<endl;
	
	pr=1-pow((1-tauxM),(n-1))+pe;		//Prob de collision p = X
	//cout <<" display pr "<< pr  <<endl;
	//sleep(2);

	if(pr==1){pr=0.99999;} // the taux equation is not define for pr = 1.
				// When the total number of users is more than 42 ==> pr = 1
	if(pr==0.5){pr=0.499999;} // to avoid dividing by zeros
	//taux=(2*(1-2*pr))/((1-2*pr)*(Wo+1)+pr*Wo*(1-pow((2*pr),m))); //Prob de transmission taux
	
	taux= ((1-pow(pr,m+1) )/(1-pr))* ( (2*(1-2*pr)*(1-pr)) / ((Wo)*(1-pow(2*pr,m+1))*(1-pr)+(1-2*pr)*(1-pow((pr),m+1))) ) ; //Prob de transmission taux 
        //cout <<" display taux "<< taux<<endl;
	//sleep(1);
	diff=(tauxM-taux);
	//cout <<" The difference is " <<diff << endl;
	
	if (diff >= 0){
			tauxR=tauxM;
			
			if (abs(diff) <= eps){
				//<<"La solution est (p,taux) "
			cout << n <<"    "<< pr << "    " << taux;
				//taux=y;
				//p=X;
				break;	
			}	
			
	}

	if (diff < 0){
			tauxL=tauxM;
			
			if (abs(diff) <= eps){
				//<<"La solution est (p,taux) "
			cout << n <<"    "<< pr << "    " << taux;
				//taux=y;
				//p=X;
				break;	
			}
			
	}
    }
}


void Ju_model::Compute_throughput(){
  // the allowed data rates: "6", "9", "12", "18", "24", "36", "48", and "54" Mbit/s
  x.push_back(6*1e6); x.push_back(9*1e6); x.push_back(12*1e6); x.push_back(18*1e6); x.push_back(24*1e6);x.push_back(36*1e6); x.push_back(48*1e6); x.push_back(54*1e6);   // declare a vector ans intitiate it.
  
  //cout <<" The size of N is   "<< N.size() <<endl;
  
  for(unsigned int i=0; i<N.size() ;i++){
    double ptr=1-pow((1-taux),N[i]);	// transmission prb for users of data rate i
    Ptr_x.push_back(ptr);
    
    //cout << "  ptr  " << ptr <<endl;
    double psx=0;
    if(N[i] != 0){
      psx= ( N[i]*taux* pow((1-taux),(n-1)) * (1-pe) ) / (1-pow((1-taux),N[i]));
    }     
    Ps_x.push_back(psx);
    
    //cout << "  psx  " << psx <<endl;
    double t= PS/x[i]; // duration of the MAC payload
    T_x.push_back(t);  // Save this duration in the vector T_x
    
    //cout << "  T_x  " << t <<endl;
    double ts= DIFS + TH + t + delta + SIFS + TACK + delta;  // time for succesful transmission
    Ts_x.push_back(ts);
    
    //cout << "  Ts_x  " << ts <<endl;
    double tr= DIFS + TH + t + delta;
    Tr_x.push_back(tr);
    
    //cout << "  Tr_x  " << tr <<endl;
    
  }
  //exit(0);
  
/*  for(unsigned int i=0; i<Ptr_x.size() ;i++){
    cout <<"   Transmission prb for ransmission of data rate " << i+1<< "  "<< Ptr_x[i]<<endl;
    cout <<"   Ts_x prb for ransmission of data rate " << i+1<< "  " << Ts_x[i]<<endl;
    cout <<"   Tr_x prb for ransmission of data rate " << i+1<< "  " << Tr_x[i]<<endl;
  }
 */ 
  //cout <<" I finish jut here  "<< N.size() <<endl;
  //exit(0);
  
  Ptr=1-(pow(1-taux,n)); // computation of the transmission prob among all stations 
  double sum = Compute_sumation(); // computation of the stwo summation at the denominator
  
  // comutation of the throughput of each zone
  for(unsigned int i=0; i<Ptr_x.size() ;i++){
    double sx= (Ps_x[i]*Ptr_x[i]*T_x[i])/((1-Ptr)*Tslot+sum);
    S_x.push_back(sx);
  }
  
  // computation of the overall throughput
  S=0; // initialization of the value of throughput
  for(unsigned int i=0; i<Ptr_x.size();i++ ){
   S=S+ S_x[i]*x[i]; // the summution over all values 
  }
  
}

double Ju_model::Compute_sumation(){
  double s=0, s1=0, s2=0; //s3=0, s4=0, s5=0, s6=0, s7=0, s7=0
  for(unsigned int i=0;i<Ptr_x.size();i++){
 
    s1 = s1+Ps_x[i]*Ptr_x[i]*Ts_x[i];
    s2 = s2+(1-Ps_x[i])*Ptr_x[i]*Tr_x[i];
  }
  s=s1+s2;
    
  return s;
}

int main (int argc, char **argv){
  
  Ju_model JM;
/*  cout<< " Enter the number of users operating at 6Mbps :  ";
  cin>>JM.n6;	// enter the 

  cout<< " Enter the number of users operating at 9Mbps :  ";
  cin>>JM.n9;	// enter the 
  
  cout<< " Enter the number of users operating at 12Mbps :  ";
  cin>>JM.n12;	// enter the 
  
  cout<< " Enter the number of users operating at 18Mbps :  ";
  cin>>JM.n18;	// enter the 
  
  cout<< " Enter the number of users operating at 24Mbps :  ";
  cin>>JM.n24;	// enter the
  
  cout<< " Enter the number of users operating at 36Mbps :  ";
  cin>>JM.n36;	// enter the 
  
  cout<< " Enter the number of users operating at 48Mbps :  ";
  cin>>JM.n48;	// enter the 
  
  cout<< " Enter the number of users operating at 54Mbps :  ";
  cin>>JM.n54;	// enter the 

*/  
JM.n6=atoi(argv[1]);	
JM.n9=atoi(argv[2]);	// enter the 
JM.n12=atoi(argv[3]);	// enter the 
JM.n18=atoi(argv[4]);	// enter the   
JM.n24=atoi(argv[5]);	// enter the
JM.n36=atoi(argv[6]);	// enter the
JM.n48=atoi(argv[7]);	// enter the 
JM.n54=atoi(argv[8]);
  
  
  JM.Print();
  JM.Initialize();
  
  //cout<< " Enter probability of error : pe ";
  //cin>>JM.pe;
  
  JM.Compute_taux_collisionPrb();
  JM.Compute_throughput();
  cout <<"  "<< JM.S<<endl;
  
  
  
}
