/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | FickianTransportFoam 
   \\    /   O peration     |
    \\  /    A nd           | 
     \\/     M anipulation  | 2024, Aalto University, Finland
-------------------------------------------------------------------------------
License
    This file is part of FickianTransportFoam library, derived from OpenFOAM.

    https://github.com/Aalto-CFD/FickianTransportFoam

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "basicFickianTransportModel.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvmSup.H"
#include "fvmDiv.H"
#include "surfaceInterpolate.H"
#include "Function2Evaluate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
basicFickianTransportModel<BasicThermophysicalTransportModel>::basicFickianTransportModel
(
    const word& type,
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo,
    bool constantLewis
)
:
    BasicThermophysicalTransportModel
    (
        type,
        momentumTransport,
        thermo
    ),
    constantLewis_(constantLewis),
    Jc_
     (
         IOobject
         (
             thermo.phasePropertyName("Jc"),
             this->thermo().T().mesh().time().name(),
             this->thermo().T().mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         this->thermo().T().mesh(),
         dimensionedScalar(dimMass/dimArea/dimTime, 0)
     ),

    implicitFlux_(true),
 
    DFuncs_(this->thermo().species().size()),
    
    Le_(this->thermo().species().size())
    
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
bool basicFickianTransportModel<BasicThermophysicalTransportModel>::read()
{
    if
    (
        BasicThermophysicalTransportModel::read()
    )
    {
        const speciesTable& species = this->thermo().species();

	if(constantLewis_) {
	   Info << "Setting Lewis numbers " << endl;
	   dictionary LewisNumberDict(this->coeffDict().subDict("Le"));

           const PtrList<volScalarField>& Y = this->thermo().Y();
           for (label i=0; i<species.size(); i++)
           {
              if (LewisNumberDict.found(Y[i].member()))
              {
                Le_[i] = readScalar(LewisNumberDict.lookup(Y[i].member()));
                Info<<"Lewis number of specie "<<Y[i].name()<<" is: "<<Le_[i]<<endl;
              }
              else {
           	Info<<"Setting Lewis number of specie "<<Y[i].name()<<" to default value 1.0 "<<endl;
           	Le_[i] = 1;
              }
           }
	}
	else 
	{
            const dictionary& Ddict = this->coeffDict().subDict("D");

            // Read the array of specie binary mass diffusion coefficient
            // functions
            forAll(species, i)
            {
                DFuncs_[i].setSize(species.size());

                forAll(species, j)
                {
                    if (j >= i)
                    {
                        const word nameij(species[i] + '-' + species[j]);
                        const word nameji(species[j] + '-' + species[i]);

                        word Dname;

                        if (Ddict.found(nameij) && Ddict.found(nameji))
                        {
                            if (i != j)
                            {
                                WarningInFunction
                                    << "Binary mass diffusion coefficients "
                                       "for Both " << nameij
                                    << " and " << nameji << " provided, using "
                                    << nameij << endl;
                            }

                            Dname = nameij;
                        }
                        else if (Ddict.found(nameij))
                        {
                            Dname = nameij;
                        }
                        else if (Ddict.found(nameji))
                        {
                            Dname = nameji;
                        }
                        else
                        {
                            FatalIOErrorInFunction(Ddict)
                                << "Binary mass diffusion coefficient for pair "
                                << nameij << " or " << nameji << " not provided"
                                << exit(FatalIOError);
                        }

                        DFuncs_[i].set
                        (
                            j,
                            Function2<scalar>::New(Dname,                                
                                 dimPressure,
                                 dimTemperature,
                                 dimKinematicViscosity,
                                 Ddict).ptr()
                        );
                    }
                }
            }
        }
        
	implicitFlux_ = this->coeffDict().lookupOrDefault("implicitHeatFlux",true);

	if (implicitFlux_ && this->thermo().he().name() == "e")
	{
	    WarningInFunction
	        << "The implicit heat-flux formulation is derived from "
	           "dh_s = c_p dT and is not strictly consistent with the "
	           "sensible internal energy equation. Falling back to the "
	           "explicit formulation (implicitHeatFlux false)." << endl;
	    implicitFlux_ = false;
	}

	Info << "Selecting "<< (implicitFlux_ ? "implicit" : "explicit") << " formulation for the heat flux" << endl;
	
        return true;
    }
    else
    {
        return false;
    }
}



template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField> basicFickianTransportModel<BasicThermophysicalTransportModel>::q() const
{
    tmp<surfaceScalarField> tmpq
    (
        surfaceScalarField::New
        (
            IOobject::groupName
            (
                "q",
                this->momentumTransport().alphaRhoPhi().group()
            ),
           -fvc::interpolate(this->alpha()*this->kappaEff())
           *fvc::snGrad(this->thermo().T())
        )
    );

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();
    
    if (Y.size())
    {
        surfaceScalarField sumJh
        (
            surfaceScalarField::New
            (
                "sumJh",
                Y[0].mesh(),
                dimensionedScalar(dimMass/dimArea/dimTime*dimEnergy/dimMass, 0)
            )
        );

        forAll(Y, i)
        {

                const volScalarField hi(this->thermo().hsi(i, p, T));

                const surfaceScalarField ji(BasicThermophysicalTransportModel::j(Y[i]));
                sumJh += ji*fvc::interpolate(hi);
            
        }



        tmpq.ref() += sumJh;
    }
    tmpq.ref() -= Jc_*fvc::interpolate(this->thermo().he());
    return tmpq;
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix> basicFickianTransportModel<BasicThermophysicalTransportModel>::divq
(
    volScalarField& he
) const
{
     correctJc();
     
     tmp<fvScalarMatrix> tmpDivq
     (
        implicitFlux_ ? -fvm::laplacian(this->alpha()*this->alphaEff(), he) :
        fvm::Su 
         (
             -fvc::laplacian(this->alpha()*this->kappaEff(), this->thermo().T()),
             he
         )
     );
     
     const PtrList<volScalarField>& Y = this->thermo().Y();
     const volScalarField& p = this->thermo().p();
     const volScalarField& T = this->thermo().T();

     if (!implicitFlux_) 
     {
     tmpDivq.ref() -=
         correction(fvm::laplacian(this->alpha()*this->alphaEff(), he));
     }

 
     surfaceScalarField sumJh
     (
         surfaceScalarField::New
         (
             "sumJh",
             he.mesh(),
             dimensionedScalar(dimMass/dimArea/dimTime*he.dimensions(), 0)
         )
     );
 
     forAll(Y, i)
     {

             const volScalarField hi(this->thermo().hsi(i, p, T));
 
             const surfaceScalarField ji(BasicThermophysicalTransportModel::j(Y[i]));
 
             sumJh += ji*fvc::interpolate(hi);
             if (implicitFlux_) {
                 sumJh += fvc::interpolate(this->alpha()*this->alphaEff()*hi)*fvc::snGrad(Y[i]);
            
             }
         
     }
 
 
     tmpDivq.ref() += fvc::div(sumJh*he.mesh().magSf());
     tmpDivq.ref() -= fvm::div(Jc_*he.mesh().magSf(),he,"div(phic,Yi_h)");
     return tmpDivq;
 }

template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField> basicFickianTransportModel<BasicThermophysicalTransportModel>::j
(
    const volScalarField& Yi
) const
{
       return BasicThermophysicalTransportModel::j(Yi) - Jc_*fvc::interpolate(Yi); 
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix> basicFickianTransportModel<BasicThermophysicalTransportModel>::divj
(
    volScalarField& Yi
) const
{
    const volScalarField& T = this->thermo().T();
    return BasicThermophysicalTransportModel::divj(Yi) - fvm::div(Jc_*T.mesh().magSf(),Yi,"div(phic,Yi_h)");  
}


template<class BasicThermophysicalTransportModel>
void basicFickianTransportModel<BasicThermophysicalTransportModel>::updateDm() const
{
    if (constantLewis_) return correctJc();
    
    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();
    const volScalarField Wm(this->thermo().W());
    Dm_.setSize(Y.size()); 
    volScalarField Dji
    (
        volScalarField::New
        (
            "Dji",
            T.mesh(),
            dimKinematicViscosity
        )
    );

    forAll(Y, i)
    {
        // Calculate first the denominators of the mixture-averaged
        // diffusion coefficients
        Dm_.set
        (
            i,
            volScalarField::New
            (
                "Dm_" + Y[i].name(),
                T.mesh(),
                dimensionedScalar(dimMoles/dimMass/dimKinematicViscosity, 0)
            )
        );
        for(label j = 0; j < i; j++)
        {

            Dji = evaluate(DFuncs_[j][i], dimKinematicViscosity, p, T);
            dimensionedScalar Wi = this->thermo().Wi(i);
            dimensionedScalar Wj = this->thermo().Wi(j);
            Dm_[i] += Y[j]/Dji*(1/Wj + (1/Wi - 1/Wj)*Y[i]);
            Dm_[j] += Y[i]/Dji*(1/Wi + (1/Wj - 1/Wi)*Y[j]);
        }
    }

    forAll(Dm_, i)
    {
        // At the limit Yk = 1, use the self-diffusion coefficients
        Dji = evaluate(DFuncs_[i][i], dimKinematicViscosity, p, T);
        Dm_.set(i, max(1-Y[i],small)/max(Dm_[i]*Wm,small/Dji));
    }       
    correctJc();
}

 template<class BasicThermophysicalTransportModel>
 void basicFickianTransportModel<BasicThermophysicalTransportModel>::predict()
 {
     BasicThermophysicalTransportModel::predict();
     updateDm();
 }
  
  
 template<class BasicThermophysicalTransportModel>
 bool basicFickianTransportModel<BasicThermophysicalTransportModel>::movePoints()
 {
     return true;
 }
  
  
 template<class BasicThermophysicalTransportModel>
 void basicFickianTransportModel<BasicThermophysicalTransportModel>::topoChange
 (
     const polyTopoChangeMap& map
 )
 {
     // Delete the cached Dm, will be re-created in predict
     Dm_.clear();
 }
  
  
template<class BasicThermophysicalTransportModel>
void basicFickianTransportModel<BasicThermophysicalTransportModel>::mapMesh
 (
     const polyMeshMap& map
 )
 {
     // Delete the cached Dm, will be re-created in predict
     Dm_.clear();
 }
  
  
template<class BasicThermophysicalTransportModel>
void basicFickianTransportModel<BasicThermophysicalTransportModel>::distribute
 (
     const polyDistributionMap& map
 )
 {
     // Delete the cached Dm, will be re-created in predict
     Dm_.clear();
 }

 template<class BasicThermophysicalTransportModel>
 const PtrList<volScalarField>&
 basicFickianTransportModel<BasicThermophysicalTransportModel>::Dm() const
 {
     if (!Dm_.size())
     {
         updateDm();
     }
  
     return Dm_;
 }
   
template<class BasicThermophysicalTransportModel>
void basicFickianTransportModel<BasicThermophysicalTransportModel>::correctJc() const
{
   
   const PtrList<volScalarField>& Y = this->thermo().Y();
   Jc_ *= scalar(0);
   forAll(Y, i)
   {
     Jc_ -= fvc::interpolate(this->alpha()*this->DEff(Y[i]))*fvc::snGrad(Y[i]);
   }  

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
