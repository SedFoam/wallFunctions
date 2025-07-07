/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "wilcoxOmegaWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "nutkWallFunctionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void wilcoxOmegaWallFunctionFvPatchScalarField::checkType()
{
    if (not isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}

void wilcoxOmegaWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os)
const
{
    os.writeKeyword("kn") << kn_ << token::END_STATEMENT << nl;
}


void wilcoxOmegaWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<wilcoxOmegaWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            wilcoxOmegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            opf.master() = master;
        }
    }
}


wilcoxOmegaWallFunctionFvPatchScalarField&
wilcoxOmegaWallFunctionFvPatchScalarField::omegaPatch(const label patchi)
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const wilcoxOmegaWallFunctionFvPatchScalarField& opf =
        refCast<const wilcoxOmegaWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<wilcoxOmegaWallFunctionFvPatchScalarField&>(opf);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

wilcoxOmegaWallFunctionFvPatchScalarField::
wilcoxOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    kn_(1e-6),
    omega_(),
    initialised_(false),
    master_(-1)
{
    checkType();
}


wilcoxOmegaWallFunctionFvPatchScalarField::
wilcoxOmegaWallFunctionFvPatchScalarField
(
    const wilcoxOmegaWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    kn_(ptf.kn_),
    omega_(),
    initialised_(false),
    master_(-1)
{
    checkType();
}


wilcoxOmegaWallFunctionFvPatchScalarField::
wilcoxOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    kn_(dict.lookupOrDefault<scalar>("kn", 1e-6)),
    omega_(),
    initialised_(false),
    master_(-1)
{
    checkType();

    // apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


wilcoxOmegaWallFunctionFvPatchScalarField::
wilcoxOmegaWallFunctionFvPatchScalarField
(
    const wilcoxOmegaWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedValueFvPatchField<scalar>(owfpsf),
    kn_(owfpsf.kn_),
    omega_(),
    initialised_(false),
    master_(-1)
{
    checkType();
}


wilcoxOmegaWallFunctionFvPatchScalarField::
wilcoxOmegaWallFunctionFvPatchScalarField
(
    const wilcoxOmegaWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(owfpsf, iF),
    kn_(owfpsf.kn_),
    omega_(),
    initialised_(false),
    master_(-1)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalarField& wilcoxOmegaWallFunctionFvPatchScalarField::omega(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            omega_ = 0.0;
        }

        return omega_;
    }

    return omegaPatch(master_).omega(init);
}


void wilcoxOmegaWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    scalarField newOmegaValues(this->patchInternalField());

    const label patchi = patch().index();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const tmp<volScalarField> tnut = turbModel.nut();
    const scalarField& nutw = tnutw();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    // Set omega
    forAll(nutw, facei)
    {
        scalar utau = max(1.0e-6, sqrt((nutw[facei] + nuw[facei])*
                      magGradUw[facei]));
        scalar knplus = kn_*utau/nuw[facei];
        scalar SR = 0e0;
        //Info<<utau<<endl;
        //Info<<knplus<<endl;
        if (knplus<=25e0)
          {
            SR = pow(50e0/knplus, 2);
          }
        else
          {
            SR = 100e0/knplus;
          }

        newOmegaValues[facei] = pow(utau, 2) * SR / nuw[facei];
    }
    this->operator==(newOmegaValues);

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void wilcoxOmegaWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeLocalEntries(os);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    wilcoxOmegaWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
