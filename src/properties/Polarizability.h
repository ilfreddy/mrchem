/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#pragma once

#include "mrchem.h"

namespace mrchem {
class Polarizability final {
public:
    explicit Polarizability(double w = 0.0, mrcpp::Coord<3> o = {0.0, 0.0, 0.0}, bool v = false)
        : frequency(w)
        , origin(o)
        , velocity(v)
        , tensor(DoubleMatrix::Zero(3, 3)) {}
    ~Polarizability() { }

    double getFrequency() const { return this->frequency; }
    bool getVelocityGauge() const { return this->velocity; }
    const mrcpp::Coord<3> getGaugeOrigin() const { return this->origin; }
    DoubleMatrix &get() { return this->tensor; }

    friend std::ostream& operator<<(std::ostream &o, const Polarizability &pol) {
        
        double w_au = 0.0;  // Only static Polarizability
        double isoPolarau = pol.tensor.trace()/3.0;
        
        Eigen::IOFormat clean_format(10, 0, ", ", "\n", "[", "]");
        double isoPolarsi = isoPolarau * 0.0; // Luca: FIX THIS
        
        int oldPrec = mrcpp::Printer::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                   Polarizability tensor                    "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Frequency:       (au)       " << std::setw(30) << w_au      <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Isotropic average(au)       " << std::setw(30) << isoPolarau<<std::endl;
        o<<" Isotropic average(SI)       TO BE FIXED                    "<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"-------------------------- Tensor --------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<                pol.tensor.format(clean_format)               <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                                                            "<<std::endl;
        mrcpp::Printer::setPrecision(oldPrec);
        return o;
    }
private:
    bool velocity;
    double frequency;
    mrcpp::Coord<3> origin;
    DoubleMatrix tensor;
};
    
} //namespace mrchem
