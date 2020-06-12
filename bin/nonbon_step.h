extern "C" void nonbon_step_( double xl[],double xr[],  
          int iacir[], int iacil[],int nonr[], int nonl[],  int& nonp,
          double& enon, double& epote, int& cartstateHandle);
double nonbon_step_cartstate(const int cartstatehandle, const int recType, const int ligType, const double dist);
extern "C" bool   use_steppot_(int & cartstatehandle);