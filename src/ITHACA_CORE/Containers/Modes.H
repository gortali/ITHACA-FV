/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------
License
    This file is part of ITHACA-FV
    ITHACA-FV is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.
Class
    EigenFunctions
Description
    Container that contains a list of modes with additional operations
SourceFiles
    Modes.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the Modes class.
/// \dir
/// Directory containing the header and source files for all the containers defined in ITHACA-FV


#ifndef Modes_H
#define Modes_H
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <Eigen/Eigen>
#include <unsupported/Eigen/SparseExtra>
#include <unsupported/Eigen/CXX11/Tensor>
#pragma GCC diagnostic pop
#include "fvCFD.H"
#include "Foam2Eigen.H"
#include "ITHACAutilities.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
  Class Modes Declaration
\*---------------------------------------------------------------------------*/

//--------------------------------------------------------------------------
///
/// @brief      Implementation of a container class derived from PtrList. It is meant to store POD modes.
///
/// @details This a child class of the PtrList class, it is meant to store modes and contains
/// some additional operations to project the modes and to reconstruct a solution from the coefficients.
///
/// @tparam     T    Type of the Geometric Field used inside the List, it can be vector or scalar
///
template<class T>
class Modes : public PtrList<GeometricField<T, fvPatchField, volMesh>>
{
    public:

        /// List of Matrices that contains the internalField and the additional matrices for the boundary patches.
        List<Eigen::MatrixXd> EigenModes;

        /// Number of patches
        int NBC;

        /// Method that convert a PtrList of modes into Eigen matrices filling the EigenModes object
        List<Eigen::MatrixXd> toEigen();

        //--------------------------------------------------------------------------
        /// @brief      Function that returns the Modes object as a standard PtrList
        ///
        /// @return     The PtrList containing the modes
        ///
        PtrList<GeometricField<T, fvPatchField, volMesh>>& toPtrList()
        {
            return static_cast<PtrList<GeometricField<T, fvPatchField, volMesh>>&>(*this);
        }

        //--------------------------------------------------------------------------
        /// @brief     A function that project an FvMatrix (OpenFoam linear System) on the modes
        ///
        /// @details   Used to project and FvMatrix that is a standard linear system of OpenFoam
        /// onto the modes defined inside the Modes container. The output is a list of matrices that
        /// contains in the first element the reduced matrix and in the second element the source term.
        ///
        /// @param      Af   The OpenFOAM fvMatrix linear system, it can be vector or scalar matrix
        ///
        /// @return     A list of Eigen Matrices of dimension 2. The first element of the list is the
        /// reduced matrix of the linear system, the second element is the reduced source term of the linear
        /// system.
        ///
        List<Eigen::MatrixXd> project(fvMatrix<T>& Af);

        //----------------------------------------------------------------------
        /// @brief      Function to reconstruct the solution starting from the
        ///             coefficients
        ///
        /// @param[in]  Coeff  The coefficients of the POD expansion
        /// @param[in]  Name   The name of the field you want to return
        ///
        /// @return     A volScalarField or volScalarField given by the reconstruction.
        ///
        GeometricField<T, fvPatchField, volMesh> reconstruct(Eigen::MatrixXd& Coeff,
                word Name);
};

#ifdef NoRepository
#   include "ModesTemplates.C"
#endif
#endif


