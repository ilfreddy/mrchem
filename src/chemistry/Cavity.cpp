#include "Cavity.h"
#include <cmath> 
#include <vector>
#include <array>

namespace mrchem {

Cavity::Cavity(std::vector<std::array<double, 3>> coord, std::vector<double> R, double slope){
  this->pos = coord;
  this->R = R;
  this->d = slope; 
}


Cavity::Cavity(std::vector<std::array<double, 3>> coord, std::vector<double> R, double slope, double eps_i, double eps_o){
  this->pos = coord;
  this->R = R;
  this->d = slope;
  this->e_i = e_i;
  this->e_o = e_o;

}

void Cavity::eval_epsilon(bool argument, bool implement){
  this->b = argument;
  this->is_linear = implement;
}


double Cavity::evalf(const double *r) const {
  double C = 1.0;
  double s, O;  
 
  for(int i = 0; i < pos.size(); i++){
    s = std::sqrt(std::pow(pos[i][0] - r[0], 2) + std::pow(pos[i][1] - r[1], 2) + std::pow(pos[i][2] - r[2], 2)) - R[i];
    O = 0.5 * (1 + std::erf(s/d));
    C *= 1 - (1 - O);
  }
  C = 1 - C;

  if(b == false){
    return C;

  }else if(b == true){

    if(is_linear == true){
      return 1/(e_o + C*(e_i - e_o));

    }else{
      return (1/e_i)*std::exp(log(e_i/e_o)*(1 - C));

    }
     

  }
}



} //namespace mrchem
