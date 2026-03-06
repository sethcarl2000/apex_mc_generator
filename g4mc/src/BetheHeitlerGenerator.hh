#ifndef BetheHeitlerGenerator_h
#define BetheHeitlerGenerator_h

#include <G4ThreeVector.hh>
#include <random>
#include <cmath>


class BetheHeitlerGenerator {
private: 

  int fVerbose_level=0; 
  double fBeam_energy = 2200.;                 /// beam energy, MeV
  static constexpr double fMass_e  = 0.501;   /// e+/e- rest mass, MeV/c^2 
  double fE_electron;                         /// electron/positrion energy in A' rest frame, MeV
    
  /// ---- geometrical / momentum acceptance of both detectors
  double fRange_L_x_sv[2] = {-85.e-3, +85.e-3}; 
  double fRange_L_y_sv[2] = {-70.e-3, +80.e-3}; 

  double fRange_R_x_sv[2] = {-85.e-3, +85.e-3};  
  double fRange_R_y_sv[2] = {-80.e-3, +70.e-3}; 

  double fRange_dp  = +0.06; 

  double fSpectrometer_p0 = 1104.;            /// Central momentum of spectrometer, MeV/c 

  double fRange_invariantMass[2] = {110., 270}; 
    
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

    G4ThreeVector vertex;        /// react vertex, HCS (m)
        
    G4ThreeVector Pp;            /// 3-momentum of the positron, HCS (MeV)
    G4ThreeVector Pe;            /// 3-momentum of the electron, HCS (MeV)

    int foil_num=4;              /// the number of the production foil we're generating fro 

    double amplitude=0.;         /// 'Amplitude' of the PDF(x,theta,P1) at current values of x,theta
    bool in_acceptance=false;    /// is this phase-point in the acceptance?  

    G4ThreeVector P_electron;    /// 3-momentum of electron in rest frame of beam electron
    double   inv_mass;            /// invariant mass of virtual photon 
  }; 

  State_t fState; /// current state of the system

  /// @param  i_foil the index of production foil to use
  /// @return z-position of production foil
  double Get_foil_z(int i_foil) { return -238.8  + ((double)i_foil)*55.; }

  /// perform random rotation of given G4ThreeVector, with given RMS rotation magnitude. preserve magnitude of vector.  
  G4ThreeVector Random_rotation(const G4ThreeVector& v, double RMS_rotation_mag); 

  /// @param is_RHRS right(true) or left(false) arm
  /// @param R react-vertex (HCS, meters)
  /// @param dR direction, HCS (any units)
  /// @param x_sv projection of track onto sieve (mm)
  /// @param y_sv projection of track onto sieve (mm)
  void Project_onto_sieve(bool is_RHRS, const G4ThreeVector& R, const G4ThreeVector& dR, double& x_sv, double& y_sv);  

  /// reset state to a sensible point
  //void ResetState(); 

  long long int fN_updates=0; 
  
public: 
  
  BetheHeitlerGenerator(double _beam_E=2200., double _spec_p0=1108.); 

  ~BetheHeitlerGenerator() {};


  void SetRange_invariantMass(double min, double max);

  void SetRange_x_sv(bool is_RHRS, double min, double max);
  void SetRange_y_sv(bool is_RHRS, double min, double max);
    
  void SetRange_p0(double _dp) { fRange_dp=_dp; }

  void Set_p0(double _p0) { fSpectrometer_p0=_p0; }  
  void Set_beamEnergy(double _E) { fBeam_energy=_E; }

  /// perform metropolis update 
  void Update(); 
  
  /// @return true if current state is inside acceptance
  inline bool InsideAcceptance() const { return fState.in_acceptance; }

  /// @return event vertex, in HCS (m)
  inline G4ThreeVector GetVertex() const { return fState.vertex; }; 

  /// @return 3-momentum of electron, in HCS (MeV/c)
  inline G4ThreeVector GetPe() const { return fState.Pe; }; 

  /// @return 3-momentum of positron, in HCS (MeV/c)
  inline G4ThreeVector GetPp() const { return fState.Pp; }; 

  /// @return current invariant mass
  inline double GetM() const { return fState.inv_mass; }; 
  
  inline long long int Get_NUpdates() const { return fN_updates; }
  
  void SetVerbosity(int _v) { fVerbose_level=_v; }

  /// @brief Compute the differential Bethe-Heitler production amplitude
  /// @param Pp 3-momentum of positron, (lab frame)
  /// @param Pe 3-momentum of electron, (lab frame)
  /// @param k energy of photon (lab frame)
  double BetheHeitler_pairprod_amplitude(const G4ThreeVector& Pp, const G4ThreeVector& Pm, const double k) const; 

};


#endif
