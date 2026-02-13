#ifndef AprimeGenerator_h
#define AprimeGenerator_h

#include <G4ThreeVector.hh>
#include <random>
#include <cmath>


class AprimeGenerator {
private: 

  int fVerbose_level=0; 
  
  double fBeam_energy = 2200;                 /// beam energy, MeV
  double fMass_Aprime;                        /// Rest mass of A', MeV/c^2
  double fE_electron;                         /// electron/positrion energy in A' rest frame, MeV
  
  /// ---- geometrical / momentum acceptance of both detectors
  double fRange_L_x_sv[2] = {-85., +85.}; 
  double fRange_L_y_sv[2] = {-70., +80.}; 
  
  double fRange_R_x_sv[2] = {-85., +85.}; 
  double fRange_R_y_sv[2] = {-80., +70.}; 
  
  double fRange_dp  = +0.06; 
  
  double fSpectrometer_p0 = 1108;         /// Central momentum of spectrometer, MeV/c 
  
  long long int fN_updates=0; 
    
    //uniform random number generator
    std::uniform_real_distribution<double> fRandGen_uniform; 

    //normally-distributed number generator
    std::normal_distribution<double> fRandGen_gaus; 
    
    //pseudorandom engine
    std::mt19937 fTwister; 

    /// @return Uniformly-distributed random number on the interval [0, 1)
    inline double Rand() { return fRandGen_uniform(fTwister); }

    /// @return Normally-distributed random number with mean=0 and RMS=1 
    inline double Gaus() { return fRandGen_gaus(fTwister); }

    struct State_t {

        G4ThreeVector vertex;              /// react vertex, HCS (m)
        double x=1.;                  /// E_A / E_0 
        double theta=0.;              /// Angle between A' and beam, lab frame 
        double phi=0.;                /// azimuthal angle of A' around beam
        
        G4ThreeVector P_electron;          /// momentum of electron in A' rest frame
        G4ThreeVector Pe_boosted;          /// momentum of electorn in the lab-frame
        G4ThreeVector Pp_boosted;          /// momentum of positron in the lab-frame
        int foil_num=4;               /// the number of the production foil we're generating from 

        double amplitude=0.;          /// 'Amplitude' of the PDF(x,theta) at current values of x,theta
        bool in_acceptance=false;     /// is this phase-point in the acceptance?  
    }; 

    State_t fState; /// current state of the system

    /// @param  i_foil the index of production foil to use
    /// @return z-position of production foil
    double Get_foil_z(int i_foil) { return -238. + ((double)i_foil)*55.; }

    /// perform random rotation of given G4ThreeVector, with given RMS rotation magnitude. preserve magnitude of vector.  
    G4ThreeVector Random_rotation(const G4ThreeVector& v, double RMS_rotation_mag); 

    //given the ratio x := Ea / E0 (A' energy to intial beam energy, in Lab), and 'theta_A' (angle between beam and A' in Lab), 
    //this returns the rel. differential cross section. 
    double A_production_amplitude(double x, double theta_A); 

  /// @param is_RHRS right(true) or left(false) arm
  /// @param R react-vertex (HCS, meters)
  /// @param dR direction, HCS (any units)
  /// @param x_sv projection of track onto sieve (mm)
  /// @param y_sv projection of track onto sieve (mm)
  void Project_onto_sieve(bool is_RHRS, const G4ThreeVector& R, const G4ThreeVector& dR, double& x_sv, double& y_sv);  

  /// reset state to a sensible point
  void ResetState(); 
  
public: 
  
  AprimeGenerator(double _m_A=140., double _beam_E=2200., double _spec_p0=1108.); 

  ~AprimeGenerator() {};

  void SetRange_x_sv(bool is_RHRS, double min, double max);
  void SetRange_y_sv(bool is_RHRS, double min, double max);
    
  void SetRange_p0(double _dp) { fRange_dp=_dp; }

  void Set_p0(double _p0) { fSpectrometer_p0=_p0; }  
  void Set_beamEnergy(double _E) { fBeam_energy=_E; }

  /// perform metropolis update 
  void Update(); 

  void Set_AprimeMass(double mA) { fMass_Aprime=mA; ResetState(); }
  
  /// @return true if current state is inside acceptance
  inline bool InsideAcceptance() const { return fState.in_acceptance; }

  /// @return event vertex, in HCS (m)
  inline G4ThreeVector GetVertex() const { return fState.vertex; }; 

  /// @return 3-momentum of electron, in HCS (MeV/c)
  inline G4ThreeVector GetPe() const { return fState.Pe_boosted; }; 

  /// @return 3-momentum of positron, in HCS (MeV/c)
  inline G4ThreeVector GetPp() const { return fState.Pp_boosted; }; 

  inline long long int Get_NUpdates() const { return fN_updates; }
  
  void SetVerbosity(int _v) { fVerbose_level=_v; }
  
};



#endif
