/** @file Atom.h
 *  @brief An atom is descriped parsed from PQR file.
 *
 *  An atom is described as a point charge (implicit solvent model) to compute
 *  electrostatic potential energy for every point in the protein.
 *
 *  @author Ilkay Sakalli
 */

#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <iomanip>
#include <vcg/space/point3.h>

using namespace vcg;
using namespace std;

/**
 * @class Atom
 *
 * @brief An atom is defined by this class.
 *
 * An atom has an atom number, a name, is assigned to a conformation, a residue with 
 * a residue number, it has a charge, a van-der-Waals radius and is assigned to a
 * specific segment. Although not all informations are necessary, they are provided.
 *
 * @author Ilkay Sakalli
 *
 */ 

class Atom {
public:
  /** @brief Constructor
   */
  Atom(){};

  /** @brief Overloaded constructor.
   *         
   *  All information of an atom can be assigned via construction.
   *  
   *  @param int Atom number
   *  @param string Atom name
   *  @param char Conformation ID
   *  @param string Residue name
   *  @param char Chain ID
   *  @param int Residue ID
   *  @param Point3 Coordinates of atom
   *  @param float Charge
   *  @param float Van-der-Waals (vdW) radius
   *  @param string Segment name
   */    
 Atom(int _atomNumber, string _atomName, char _confID, string _residueName, char _chainID, int _residueNumber, Point3<float> _coord, float _charge, float _radius, string _segName):
  atomNumber(_atomNumber),
    atomName(_atomName),
    confID(_confID),
    residueName(_residueName),
    chainID(_chainID),
    residueNumber(_residueNumber),
    coord(_coord),
    charge(_charge),
    radius(_radius),
    segName(_segName)
    {
      
    }
  
  /** @brief Prints out all atom information.
   *         
   *  @return Void.
   */    
  void print(){
    cout << "Atom number : " << atomNumber << endl;
    cout << "Atom name   : " << atomName << endl;
    cout << "Conform. nr : " << confID << endl;
    cout << "Residue name: " << residueName << endl;
    cout << "Chain nr    : " << chainID << endl;
    cout << "Residue nr  : " << residueNumber << endl;
    cout << "Coordinates : (" << coord.X() << ", " << coord.Y() << ", " << coord.Z() << ")" << endl;
    cout << "Charge      : " << charge << endl;
    cout << "Radius      : " << radius << endl;
    cout << "Segment name: " << segName << endl;
  }
  
  /** @brief Returns the coordinates of the atom.
   *         
   *  @return Point3<float> Returns x-, y-, z-coordinates of the current atom.
   */    
  Point3<float> getCoord(){
    return coord;
  }
  
  /** @brief Returns the radius of the atom.
   *         
   *  @return float Radius of atom in Angstroems.
   */    
  float getRadius(){
    return radius;
  }
  
  /** @brief Returns the residue number the atom is assigned to.
   *         
   *  @return int Returns the residue number.
   */    
  int getResidueNumber(){
    return residueNumber;
  }
  
  /** @brief Returns the residue name the atom is assigned to.
   *         
   *  @return string Returns the residue name the atom is assigned to.
   */    
  string getResidueName(){
    return residueName;
  }
  
  /** @brief Sets the residue name.
   *         
   *  @param string Residue name the atom should be assigned to.
   *  @return Void.
   */    
  void setResidueName(string _residueName){
    residueName = _residueName;
  }
  
  /** @brief Returns the chain ID the atom is assigned to.
   *         
   *  @return string Returns the chain ID the atom is assigned to.
   */    
  char getChainID(){
    return chainID;
  }
  
  /** @brief Returns the atom name the atom is assigned to.
   *         
   *  @return string Returns the atom name the atom is assigned to.
   */    
  string getAtomName(){
    return atomName;
  }
  
  /** @brief Sets the atom charge.
   *         
   *  @param string Charge which should be assigned to the atom.
   *  @return Void.
   */    
  void setCharge(float _charge){
    charge = _charge;
  }
  
  /** @brief Returns the chare of the atom.
   *         
   *  @return float Returns the charge of the atom.
   */    
  float getCharge(){
    return charge;
  }
  

  /** @brief Returns the atom ID.
   *         
   *  @return int Returns the atom ID.
   */    
  int getAtomNumber(){
    return atomNumber;
  }
  

  /** @brief Prints out the atom information as a line in the PQR file
   *         
   *  @return string Formatted PQR line representation of the current atom.
   */    
  string pqrLine(){
    stringstream ss;
    ss << "ATOM" << setw(7) << atomNumber << setw(5) << atomName << setw(1)
       << confID << setw(3) << residueName << setw(2) << chainID << setw(4)
       << residueNumber << setw(12) << coord.X() << setw(8) << coord.Y()
       << setw(8) << coord.Z() << setw(6) << charge << setw(6)
       << radius << setw(7) << segName;
    return ss.str();
  }
  
 private:
  /// ID of current Atom
  int atomNumber;
  /// Name of current Atom
  string atomName;
  /// ID of conformation
  char confID;
  /// Name of residue where current atom belongs to
  string residueName;
  /// ID of chain where current atom belongs to
  char chainID;
  /// Number of residue where current atom belongs to
  int residueNumber;
  /// Coordinate of current atom
  Point3<float> coord;
  /// Charge of current atom in Coulomb
  float charge;
  /// Van-der-Waals radius of current atom in Angstroems
  float radius;
  /// Name of segment where current atom belongs to
  string segName;
};

#endif
