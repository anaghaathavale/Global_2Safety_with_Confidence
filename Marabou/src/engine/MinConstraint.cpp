/*********************                                                        */
/*! \file MinConstraint.cpp
 ** \verbatim
 ** Top contributors (to current version):
 **   Haoze Wu, Derek Huang, Guy Katz, Shantanu Thakoor
 ** This file is part of the Marabou project.
 ** Copyright (c) 2017-2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved. See the file COPYING in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** See the description of the class in MaxConstraint.h.
 **/

#include "MinConstraint.h"

#include "Debug.h"
#include "FloatUtils.h"
#include "ITableau.h"
#include "InputQuery.h"
#include "List.h"
#include "MStringf.h"
#include "MarabouError.h"
#include "PiecewiseLinearCaseSplit.h"
#include "Statistics.h"


#include <algorithm>

#ifdef _WIN32
#undef max
#undef min
#endif

MinConstraint::MinConstraint( unsigned f, const Set<unsigned> &elements )
    : PiecewiseLinearConstraint( elements.size() )
    , _f( f )
    , _elements( elements )
    , _initialElements( elements )
    , _obsolete( false )
    , _minUpperBound( FloatUtils::infinity() )
    , _haveFeasibleEliminatedPhases( false )
    , _minValueOfEliminatedPhases( FloatUtils::infinity() ) //doubt
{
}

MinConstraint::MinConstraint( const String &serializedMin )
{
    String constraintType = serializedMin.substring( 0, 3 );
    ASSERT( constraintType == String( "min" ) );

    // Remove the constraint type in serialized form
    String serializedValues = serializedMin.substring( 4, serializedMin.length() - 4 );
    List<String> values = serializedValues.tokenize( "," );

    auto valuesIter = values.begin();
    unsigned f = atoi( valuesIter->ascii() );
    ++valuesIter;

    Set<unsigned> elements;
    for ( ; *valuesIter != "e"; ++valuesIter )
        elements.insert( atoi( valuesIter->ascii() ) );

    // Save flag and values indicating eliminated variables
    ++valuesIter;
    bool eliminatedVariableFromString = ( atoi( valuesIter->ascii() ) == 1 );

    ++valuesIter;
    double minValueOfEliminatedFromString = FloatUtils::infinity();

    if ( eliminatedVariableFromString )
        minValueOfEliminatedFromString = std::stod( valuesIter->ascii() );

    *(this) = MinConstraint( f, elements );
    _haveFeasibleEliminatedPhases = eliminatedVariableFromString;
    _minValueOfEliminatedPhases = minValueOfEliminatedFromString;
}

MinConstraint::~MinConstraint()
{
    _elements.clear();
}

PiecewiseLinearFunctionType MinConstraint::getType() const
{
    return PiecewiseLinearFunctionType::MIN;
}

PiecewiseLinearConstraint *MinConstraint::duplicateConstraint() const
{
    MinConstraint *clone = new MinConstraint( _f, _elements );
    *clone = *this;
    this->initializeDuplicateCDOs( clone );
    return clone;
}

void MinConstraint::restoreState( const PiecewiseLinearConstraint *state )
{
    const MinConstraint *min = dynamic_cast<const MinConstraint *>( state );

    CVC4::context::CDO<bool> *activeStatus = _cdConstraintActive;
    CVC4::context::CDO<PhaseStatus> *phaseStatus = _cdPhaseStatus;
    CVC4::context::CDList<PhaseStatus> *infeasibleCases = _cdInfeasibleCases;
    *this = *min;
    _cdConstraintActive = activeStatus;
    _cdPhaseStatus = phaseStatus;
    _cdInfeasibleCases = infeasibleCases;
}

void MinConstraint::registerAsWatcher( ITableau *tableau )
{
    for ( unsigned element : _elements )
    {
        tableau->registerToWatchVariable( this, element );
        if ( _elementToAux.exists( element ) )
            tableau->registerToWatchVariable( this, _elementToAux[element] );
    }

    if ( !_elements.exists( _f ) )
        tableau->registerToWatchVariable( this, _f );
}

void MinConstraint::unregisterAsWatcher( ITableau *tableau )
{
    for ( unsigned element : _initialElements )
    {
        tableau->unregisterToWatchVariable( this, element );
        if ( _elementToAux.exists( element ) )
            tableau->unregisterToWatchVariable( this, _elementToAux[element] );
    }

    if ( !_initialElements.exists( _f ) )
        tableau->unregisterToWatchVariable( this, _f );
}

void MinConstraint::notifyUpperBound( unsigned variable, double value )
{
    if ( _statistics )
        _statistics->incLongAttribute(
            Statistics::NUM_BOUND_NOTIFICATIONS_TO_PL_CONSTRAINTS );

    if ( _boundManager == nullptr )
    {
        if ( existsUpperBound( variable ) &&
             !FloatUtils::lt( value, getUpperBound( variable ) ) )
            return;

        setUpperBound( variable, value );
    }

    /*
      See if we can eliminate any cases.
    */
    if ( _auxToElement.exists( variable ) && FloatUtils::isPositive( value ) )
    {
        // The case that this variable is an aux variable.
        // We can eliminate the corresponding case if the aux variable is
        // positive.
        eliminateCase( _auxToElement[variable] );
    }
    else if ( variable == _f || _elements.exists( variable ) )
    {
        // If the variable is either in _element or _f, there will be eliminated
        // cases when the new lowerBound is greater than the _maxLowerBound.
        if ( FloatUtils::lt( value, _minUpperBound ) )
        {
            _minUpperBound = value;
            List<unsigned> toRemove;
            for ( auto element : _elements )
            {
                if ( element == variable )
                    continue;
                if ( existsLowerBound( element ) &&
                     FloatUtils::gt( getLowerBound( element ), value ) )
                {
                    toRemove.append( element );
                }
            }
            for ( unsigned removeVar : toRemove )
                eliminateCase( removeVar );

            _haveFeasibleEliminatedPhases =
                FloatUtils::gte( _minUpperBound, _minValueOfEliminatedPhases ); //doubt
        }
    }

    if ( phaseFixed() )
        _phaseStatus = ( _haveFeasibleEliminatedPhases
                             ? MIN_PHASE_ELIMINATED
                             : variableToPhase( *_elements.begin() ) );

    if ( isActive() && _boundManager )
    {
        // TODO: optimize this. Don't need to recompute ALL possible bounds,
        // Can focus only on the newly learned bound and possible consequences.
        List<Tightening> tightenings;
        getEntailedTightenings( tightenings );
        for ( const auto &tightening : tightenings )
        {
            if ( tightening._type == Tightening::LB )
                _boundManager->tightenLowerBound( tightening._variable,
                                                  tightening._value );
            else if ( tightening._type == Tightening::UB )
                _boundManager->tightenUpperBound( tightening._variable,
                                                  tightening._value );
        }
    }
}

void MinConstraint::notifyLowerBound( unsigned variable, double value )
{
    if ( _statistics )
        _statistics->incLongAttribute(
            Statistics::NUM_BOUND_NOTIFICATIONS_TO_PL_CONSTRAINTS );

    if ( _boundManager == nullptr )
    {
        if ( existsLowerBound( variable ) &&
             !FloatUtils::gt( value, getLowerBound( variable ) ) )
            return;

        setLowerBound( variable, value );
    }
    /*
      See if we can eliminate any cases.
    */
    if ( _auxToElement.exists( variable ) )
    {
        /* The case that this variable is an aux variable.
           If the aux variable is 0. This means f - auxToElement[aux] = aux = 0.
           We can eliminate all other cases and set
           _haveFeasibleEliminatedPhases to false;
        */
        if ( FloatUtils::isZero( value ) )
        {
            unsigned currentElement = _auxToElement[variable];
            Set<unsigned> toRemove = _elements;
            toRemove.erase( currentElement );
            for ( const auto &element : toRemove )
                eliminateCase( element );
            _haveFeasibleEliminatedPhases = false;
        }
    }
    else if ( _elements.exists( variable ) )
    {
        // If the upper bound of this variable is below the _maxLowerBound,
        // this case is infeasible.
        if ( FloatUtils::gt( value, _minUpperBound ) )
            eliminateCase( variable );
    }

    if ( phaseFixed() )
        _phaseStatus = ( _haveFeasibleEliminatedPhases
                             ? MIN_PHASE_ELIMINATED
                             : variableToPhase( *_elements.begin() ) ); //doubt

    // There is no need to recompute the max lower bound and max index here.

    if ( isActive() && _boundManager )
    {
        // TODO: optimize this. Don't need to recompute ALL possible bounds,
        // Can focus only on the newly learned bound and possible consequences.
        List<Tightening> tightenings;
        getEntailedTightenings( tightenings );
        for ( const auto &tightening : tightenings )
        {
            if ( tightening._type == Tightening::LB )
                _boundManager->tightenLowerBound( tightening._variable,
                                                  tightening._value );
            else if ( tightening._type == Tightening::UB )
                _boundManager->tightenUpperBound( tightening._variable,
                                                  tightening._value );
        }
    }
}

void MinConstraint::getEntailedTightenings( List<Tightening> &tightenings ) const
{
    // Lower and upper bounds for the f variable
    double fLB = existsLowerBound( _f ) ? getLowerBound( _f ) : FloatUtils::negativeInfinity();
    double fUB = existsUpperBound( _f ) ? getUpperBound( _f ) : FloatUtils::infinity();

    // Compute the maximal bounds (lower and upper) for the elements
    double minElementLB = FloatUtils::infinity();
    double minElementUB = FloatUtils::infinity();

    for ( const auto &element : _elements )
    {
        if ( existsUpperBound( element ) )
            minElementUB = FloatUtils::min( getUpperBound( element ), minElementUB );

        if ( !existsLowerBound( element ) )
            minElementLB = FloatUtils::negativeInfinity();
        else
            minElementLB = FloatUtils::min( getLowerBound( element ), minElementLB );
    }

    minElementLB = FloatUtils::min( _minValueOfEliminatedPhases, minElementLB );
    minElementUB = FloatUtils::min( _minValueOfEliminatedPhases, minElementUB );

    // f_UB and maxElementUB need to be equal. If not, the lower of the two wins.
    if ( FloatUtils::areDisequal( fLB, minElementLB ) )
    {
        if ( FloatUtils::lt( fLB, minElementLB ) )
        {
            tightenings.append( Tightening( _f, minElementLB, Tightening::LB ) );
        }
        else
        {
            // f_UB <= maxElementUB
            for ( const auto &element : _elements )
            {
                if ( !existsLowerBound( element ) ||
                     FloatUtils::lt( getLowerBound( element ), fLB ) )
                    tightenings.append
                        ( Tightening( element, fLB, Tightening::LB ) );
            }
        }
    }

    // fLB cannot be smaller than maxElementLB
    if ( FloatUtils::gt( fUB, minElementUB ) )
        tightenings.append( Tightening( _f, minElementUB, Tightening::UB ) );

    // TODO: bound tightening for aux vars.
}

bool MinConstraint::participatingVariable( unsigned variable ) const
{
    return ( variable == _f ) || _elements.exists( variable ) ||
        _auxToElement.exists( variable );
}

List<unsigned> MinConstraint::getParticipatingVariables() const
{
    List<unsigned> result;

    for ( auto element : _elements )
        result.append( element );

    for ( auto pair : _auxToElement )
        result.append( pair.first );

    if ( !_elements.exists( _f ) )
        result.append( _f );

    return result;
}

List<unsigned> MinConstraint::getElements() const
{
    List<unsigned> result;
    for ( auto element : _elements )
        result.append( element );
    return result;
}

unsigned MinConstraint::getF() const
{
    return _f;
}

bool MinConstraint::satisfied() const
{
    DEBUG({
            if ( !( existsAssignment( _f ) ) )
                throw MarabouError
                    ( MarabouError::PARTICIPATING_VARIABLE_MISSING_ASSIGNMENT,
                      Stringf( "f(x%u) assignment missing.", _f ).ascii() );
            for ( const auto &element : _elements )
                if ( !( existsAssignment( element ) ) )
                    throw MarabouError
                        ( MarabouError::PARTICIPATING_VARIABLE_MISSING_ASSIGNMENT,
                          Stringf( "input(x%u) assignment missing.",
                                   element ).ascii() );
        });

    double fValue = getAssignment( _f );
    double minValue = _minValueOfEliminatedPhases;  //doubt
    for ( const auto &element : _elements )
    {
        double currentValue = getAssignment( element );
        if ( FloatUtils::lt( currentValue, minValue ) )
            minValue = currentValue;
    }
    return FloatUtils::areEqual( minValue, fValue );
}

bool MinConstraint::isCaseInfeasible( unsigned variable ) const
{
    return PiecewiseLinearConstraint::isCaseInfeasible( variableToPhase( variable ) );
}

List<PiecewiseLinearConstraint::Fix> MinConstraint::getPossibleFixes() const
{
    // Reluplex does not currently work with Gurobi.
    ASSERT( _gurobi == NULL );
    return List<PiecewiseLinearConstraint::Fix>();
}

List<PiecewiseLinearConstraint::Fix> MinConstraint::getSmartFixes( ITableau * ) const
{
    // Reluplex does not currently work with Gurobi.
    ASSERT( _gurobi == NULL );
    return getPossibleFixes();
}

List<PhaseStatus> MinConstraint::getAllCases() const
{
    List<PhaseStatus> cases;
    for ( unsigned element : _elements )
        cases.append( variableToPhase( element ) );

    if ( _haveFeasibleEliminatedPhases )
        cases.append( MIN_PHASE_ELIMINATED );

    return cases;
}

List<PiecewiseLinearCaseSplit> MinConstraint::getCaseSplits() const
{
    List<PiecewiseLinearCaseSplit> splits;

    for ( const auto &phase : getAllCases() )
        splits.append( getCaseSplit( phase ) );

    return splits;
}

bool MinConstraint::phaseFixed() const
{
    return ( ( _elements.size() == 1 && !_haveFeasibleEliminatedPhases ) ||
             ( _elements.size() == 0 && _haveFeasibleEliminatedPhases ) );
}

bool MinConstraint::isImplication() const
{
    return _elements.exists( _f ) || numFeasibleCases() == 1u;
}

PiecewiseLinearCaseSplit MinConstraint::getImpliedCaseSplit() const
{
    ASSERT( phaseFixed() );

    PhaseStatus phase = getPhaseStatus();

    ASSERT( phase != PHASE_NOT_FIXED );

    return getCaseSplit( phase );
}

PiecewiseLinearCaseSplit MinConstraint::getValidCaseSplit() const
{
    return getImpliedCaseSplit();
}

PiecewiseLinearCaseSplit MinConstraint::getCaseSplit( PhaseStatus phase ) const
{
    if ( phase == MIN_PHASE_ELIMINATED )
    {
        PiecewiseLinearCaseSplit eliminatedPhase;
        eliminatedPhase.storeBoundTightening
            ( Tightening( _f, _minValueOfEliminatedPhases,
                          Tightening::LB ) );
        eliminatedPhase.storeBoundTightening
            ( Tightening( _f, _minValueOfEliminatedPhases,
                          Tightening::UB ) );
        return eliminatedPhase;
    }
    else
    {
        unsigned argMax = phaseToVariable( phase );
        PiecewiseLinearCaseSplit maxPhase;

        if ( argMax != _f )
        {
            // We had f - argMax = aux and
            maxPhase.storeBoundTightening( Tightening( _elementToAux[argMax], 0,
                                                           Tightening::UB ) );
        }
        return maxPhase;
    }
}

void MinConstraint::updateVariableIndex( unsigned oldIndex, unsigned newIndex )
{
    // Variable re-indexing can only occur in preprocessing before Gurobi is
    // registered.
    ASSERT( _gurobi == NULL );

    _lowerBounds[newIndex] = _lowerBounds[oldIndex];
    _upperBounds[newIndex] = _upperBounds[oldIndex];

    if ( oldIndex == _f )
        _f = newIndex;
    else if ( _elements.exists( oldIndex ) )
    {
        _elements.erase( oldIndex );
        _elements.insert( newIndex );

        unsigned auxVar = _elementToAux[oldIndex];
        _elementToAux.erase( oldIndex );
        _elementToAux[newIndex] = auxVar;
        _auxToElement[auxVar] = newIndex;

        if ( _phaseStatus == variableToPhase( oldIndex ) )
            _phaseStatus = variableToPhase( newIndex ) ;
    }
    else
    {
        ASSERT( _auxToElement.exists( oldIndex ) );
        unsigned element = _auxToElement[oldIndex];
        _elementToAux[element] = newIndex;
        _auxToElement.erase( oldIndex );
        _auxToElement[newIndex] = element;
    }
}

bool MinConstraint::constraintObsolete() const
{
    return _obsolete;
}

void MinConstraint::eliminateVariable( unsigned var, double value )
{
    if ( var == _f )
    {
        _obsolete = true;
        return;
    }
    else if ( _elements.exists( var ) )
    {
        eliminateCase( var );

        _minUpperBound = FloatUtils::min( value, _minUpperBound );

        _minValueOfEliminatedPhases = FloatUtils::min
            ( value, _minValueOfEliminatedPhases );

        _haveFeasibleEliminatedPhases =
            FloatUtils::lte( _minValueOfEliminatedPhases, _minUpperBound ); //doubt
    }
    else if ( _auxToElement.exists( var ) )
    {
        // Aux var
        unsigned currentElement = _auxToElement[var];
        if ( FloatUtils::isZero( value ) )
        {
            _obsolete = true;
        }
        else
        {
            // The case corresponding to aux is infeasible.
            eliminateCase( currentElement );
        }
    }

    if ( phaseFixed() )
        _phaseStatus = ( _haveFeasibleEliminatedPhases ?
                         MIN_PHASE_ELIMINATED :
                         variableToPhase( *_elements.begin() ) );

    if ( _elements.size() == 0 )
        _obsolete = true;
}

void MinConstraint::transformToUseAuxVariables( InputQuery &inputQuery )
{
    if ( _auxToElement.size() > 0 )
        return;

    bool fInInput = false;
    for ( const auto &element : _elements )
    {
        // If element is equal to _f, skip this step.
        // The reason is to avoid adding equations like `1.00x00 -1.00x00 -1.00x01 = 0.00`.
        if ( element == _f )
        {
            fInInput = true;
            continue;
        }
        // Create an aux variable
        unsigned auxVariable = inputQuery.getNumberOfVariables();
        inputQuery.setNumberOfVariables( auxVariable + 1 );

        // f >= element, or f - element - aux = 0, for non-negative aux
        Equation equation( Equation::EQ );
        equation.addAddend( 1.0, element );
        equation.addAddend( -1.0, _f );
        equation.addAddend( -1.0, auxVariable );
        equation.setScalar( 0 );
        inputQuery.addEquation( equation );

        // Set the bounds for the aux variable
        inputQuery.setLowerBound( auxVariable, 0 );

        _elementToAux[element] = auxVariable;
        _auxToElement[auxVariable] = element;
    }
    if ( fInInput )
    {
        List<unsigned> toRemove;
        for ( const auto &element : _elements )
        {
            if ( element == _f )
                continue;
            toRemove.append( element );
        }
        for ( const auto &element : toRemove )
            eliminateCase( element );
        _obsolete = true;
    }
}

void MinConstraint::getCostFunctionComponent( LinearExpression &cost,
                                              PhaseStatus phase ) const
{
    // If the constraint is not active or is fixed, it contributes nothing
    if( !isActive() || phaseFixed() )
        return;

    if ( phase == MIN_PHASE_ELIMINATED )
    {
        // The cost term corresponding to this phase is f - maxValueOfEliminated.
        if ( !cost._addends.exists( _f ) )
            cost._addends[_f] = 0;
        cost._addends[_f] = cost._addends[_f] + 1;
        cost._constant = _minValueOfEliminatedPhases - cost._constant ;
    }
    else
    {
        unsigned element = phaseToVariable( phase );
        unsigned aux = _elementToAux[element];
        if ( !cost._addends.exists( aux ) )
            cost._addends[aux] = 0;
        cost._addends[aux] += 1;
    }
}

PhaseStatus MinConstraint::getPhaseStatusInAssignment( const Map<unsigned, double>
                                                       &assignment ) const
{
    auto byAssignment = [&](const unsigned& a, const unsigned& b) {
                            return assignment[a] < assignment[b];
                        };
    unsigned smallestVariable =  *std::min_element( _elements.begin(),
                                                   _elements.end(),
                                                   byAssignment );
    double value = assignment[smallestVariable];
    if ( _haveFeasibleEliminatedPhases &&
         FloatUtils::gt( value, _minValueOfEliminatedPhases ) )
        return MIN_PHASE_ELIMINATED;
    else
        return variableToPhase( smallestVariable );
}

String MinConstraint::serializeToString() const
{
    // Output format: max,f,element_1,element_2,element_3,...
    Stringf output = Stringf( "min,%u", _f );
    for ( const auto &element : _elements )
        output += Stringf( ",%u", element );

    // Special delimiter ",e" represents elimination flag and variables
    output += Stringf( ",e" );
    output += Stringf( ",%u", _haveFeasibleEliminatedPhases ? 1 : 0 );
    if ( _haveFeasibleEliminatedPhases )
        output += Stringf( ",%f", _minValueOfEliminatedPhases );
    else
        // Will be ignored in any case
        output += Stringf( ",%u", 0 );

    return output;
}

void MinConstraint::eliminateCase( unsigned variable )
{
    if ( _cdInfeasibleCases )
    {
        markInfeasible( variableToPhase( variable ) );
    }
    else
    {
        _elements.erase( variable );
        if ( _elementToAux.exists( variable ) )
        {
            unsigned aux = _elementToAux[variable];
            _elementToAux.erase( variable );
            _auxToElement.erase( aux );
        }
    }
}

bool MinConstraint::haveOutOfBoundVariables() const
{
    double fValue = getAssignment( _f );
    if ( FloatUtils::gt( getLowerBound( _f ), fValue ) ||
         FloatUtils::lt( getUpperBound( _f ), fValue ) )
        return true;

    for ( const auto &element : _elements )
    {
        double value = getAssignment( element );
        if ( FloatUtils::gt( getLowerBound( element ), value ) ||
             FloatUtils::lt( getUpperBound( element ), value ) )
        return true;
        unsigned aux = _elementToAux[element];
        double auxValue = getAssignment( aux );
        if ( FloatUtils::gt( getLowerBound( aux ), auxValue ) ||
             FloatUtils::lt( getUpperBound( aux ), auxValue ) )
            return true;

    }
    return false;
}