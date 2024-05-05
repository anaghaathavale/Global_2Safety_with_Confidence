
/*********************                                                        */
/*! \file AcasParser.cpp
 ** \verbatim
 ** Top contributors (to current version):
 **   Guy Katz, Andrew Wu
 ** This file is part of the Marabou project.
 ** Copyright (c) 2017-2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved. See the file COPYING in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** [[ Add lengthier description here ]]

**/

#include "AcasParser.h"
#include "FloatUtils.h"
#include "InputParserError.h"
#include "InputQuery.h"
#include "MString.h"
#include "ReluConstraint.h"

AcasParser::NodeIndex::NodeIndex( unsigned layer, unsigned node )
    : _layer( layer )
    , _node( node )
{
}

bool AcasParser::NodeIndex::operator<( const NodeIndex &other ) const
{
    if ( _layer != other._layer )
        return _layer < other._layer;

    return _node < other._node;
}

AcasParser::AcasParser( const String &path )
    : _acasNeuralNetwork( path )
{
}

void AcasParser::generateQuery( InputQuery &inputQuery )
{
    // First encode the actual network
    // _acasNeuralNetwork doesn't count the input layer, so add 1
    unsigned numberOfLayers = _acasNeuralNetwork.getNumLayers() + 1;
    unsigned inputLayerSize = _acasNeuralNetwork.getLayerSize( 0 );
    unsigned outputLayerSize = _acasNeuralNetwork.getLayerSize( numberOfLayers - 1 );

    unsigned numberOfInternalNodes = 0;
    for ( unsigned i = 1; i < numberOfLayers - 1; ++i )
        numberOfInternalNodes += _acasNeuralNetwork.getLayerSize( i );

    printf( "Number of layers: %u. Input layer size: %u. Output layer size: %u. Number of ReLUs: %u\n",
            numberOfLayers, inputLayerSize, outputLayerSize, numberOfInternalNodes );

    // The total number of variables required for the encoding is computed as follows:
    //   1. Each input node appears once
    //   2. Each internal node has a B variable and an F variable
    //   3. Each output node appears once
    unsigned numberOfVariables = inputLayerSize + ( 2 * numberOfInternalNodes ) + outputLayerSize;
    printf( "Total number of variables: %u\n", numberOfVariables );

    inputQuery.setNumberOfVariables( numberOfVariables );

    // Next, we want to map each node to its corresponding
    // variables. We group variables according to this order: f's from
    // layer i, b's from layer i+1, and repeat.
    unsigned currentIndex = 0;
    for ( unsigned i = 1; i < numberOfLayers; ++i )
    {
        unsigned previousLayerSize = _acasNeuralNetwork.getLayerSize( i - 1 );
        unsigned currentLayerSize = _acasNeuralNetwork.getLayerSize( i );

        // First add the F variables from layer i-1
        for ( unsigned j = 0; j < previousLayerSize; ++j )
        {
            _nodeToF[NodeIndex( i - 1, j )] = (currentIndex);
            ++currentIndex;
        }

        // Now add the B variables from layer i
        for ( unsigned j = 0; j < currentLayerSize; ++j )
        {
            _nodeToB[NodeIndex( i, j )] = (currentIndex);
            ++currentIndex;
        }
    }

    // Now we set the variable bounds. Input bounds are
    // given as part of the network. B variables are
    // unbounded, and F variables are non-negative.
    for ( unsigned i = 0; i < inputLayerSize; ++i )
    {
        double min, max;
        _acasNeuralNetwork.getInputRange( i, min, max );

        inputQuery.setLowerBound( _nodeToF[NodeIndex(0, i)], min );
        inputQuery.setUpperBound( _nodeToF[NodeIndex(0, i)], max );
    }

    for ( const auto &fNode : _nodeToF )
    {
        // Be careful not to override the bounds for the input layer
        if ( fNode.first._layer != 0 )
        {
            inputQuery.setLowerBound( fNode.second, 0.0 );
            inputQuery.setUpperBound( fNode.second, FloatUtils::infinity() );
        }
    }

    for ( const auto &fNode : _nodeToB )
    {
        inputQuery.setLowerBound( fNode.second, FloatUtils::negativeInfinity() );
        inputQuery.setUpperBound( fNode.second, FloatUtils::infinity() );
    }

    // Next come the actual equations
    for ( unsigned layer = 0; layer < numberOfLayers - 1; ++layer )
    {
        unsigned targetLayerSize = _acasNeuralNetwork.getLayerSize( layer + 1 );
        for ( unsigned target = 0; target < targetLayerSize; ++target )
        {
            // This will represent the equation:
            //   sum - b + fs = -bias
            Equation equation;

            // The b variable
            unsigned bVar = _nodeToB[NodeIndex(layer + 1, target)];
            equation.addAddend( -1.0, bVar );

            // The f variables from the previous layer
            for ( unsigned source = 0; source < _acasNeuralNetwork.getLayerSize( layer ); ++source )
            {
                unsigned fVar = _nodeToF[NodeIndex(layer, source)];
                equation.addAddend( _acasNeuralNetwork.getWeight( layer, source, target ), fVar );
            }

            // The bias
            equation.setScalar( -_acasNeuralNetwork.getBias( layer + 1, target ) );

            // Add the equation to the input query
            inputQuery.addEquation( equation );
        }
    }

    // Add the ReLU constraints
    for ( unsigned i = 1; i < numberOfLayers - 1; ++i )
    {
        unsigned currentLayerSize = _acasNeuralNetwork.getLayerSize( i );

        for ( unsigned j = 0; j < currentLayerSize; ++j )
        {
            unsigned b = _nodeToB[NodeIndex(i, j)];
            unsigned f = _nodeToF[NodeIndex(i, j)];
            PiecewiseLinearConstraint *relu = new ReluConstraint( b, f );

            inputQuery.addPiecewiseLinearConstraint( relu );
        }
    }

    // Mark the input and output variables
    for ( unsigned i = 0; i < inputLayerSize; ++i )
        inputQuery.markInputVariable( _nodeToF[NodeIndex( 0, (i) )], (i) );

    for ( unsigned i = 0; i < outputLayerSize; ++i )
        inputQuery.markOutputVariable( _nodeToB[NodeIndex( numberOfLayers - 1, (i) )], (i) );
}

unsigned AcasParser::getNumInputVaribales() const
{
    return _acasNeuralNetwork.getLayerSize( 0 );
}

unsigned AcasParser::getNumOutputVariables() const
{
    return _acasNeuralNetwork.getLayerSize( _acasNeuralNetwork.getNumLayers() );
}

unsigned AcasParser::getInputVariable( unsigned index ) const
{
    return getFVariable( 0, index );
}

unsigned AcasParser::getOutputVariable( unsigned index ) const
{
    return getBVariable( _acasNeuralNetwork.getNumLayers(), index );
}

unsigned AcasParser::getBVariable( unsigned layer, unsigned index ) const
{
    NodeIndex nodeIndex( layer, index );
    if ( !_nodeToB.exists( nodeIndex ) )
        throw InputParserError( InputParserError::VARIABLE_INDEX_OUT_OF_RANGE );

    return _nodeToB.get( nodeIndex );
}

unsigned AcasParser::getFVariable( unsigned layer, unsigned index ) const
{
    NodeIndex nodeIndex( layer, index );
    if ( !_nodeToF.exists( nodeIndex ) )
        throw InputParserError( InputParserError::VARIABLE_INDEX_OUT_OF_RANGE );

    return _nodeToF.get( nodeIndex );
}

void AcasParser::evaluate( const Vector<double> &inputs, Vector<double> &outputs ) const
{
    _acasNeuralNetwork.evaluate( inputs, outputs,
                                 _acasNeuralNetwork.getLayerSize( _acasNeuralNetwork.getNumLayers() ) );
}




void AcasParser::generateQueryMod( InputQuery &inputQuery, unsigned alpha )
{
    // First encode the actual network
    // _acasNeuralNetwork doesn't count the input layer, so add 1
    unsigned numberOfLayers = _acasNeuralNetwork.getNumLayers() + 1;
    unsigned inputLayerSize = _acasNeuralNetwork.getLayerSize( 0 );
    unsigned outputLayerSize = _acasNeuralNetwork.getLayerSize( numberOfLayers - 1 );
    //unsigned alpha = 22;
    unsigned numberOfInternalNodes = 0;
    for ( unsigned i = 1; i < numberOfLayers - 1; ++i )
        numberOfInternalNodes += _acasNeuralNetwork.getLayerSize( i );

    printf( "Number of layers: %u. Input layer size: %u. Output layer size: %u. Number of ReLUs: %u\n",
            numberOfLayers, inputLayerSize, outputLayerSize, numberOfInternalNodes );

    // The total number of variables required for the encoding is computed as follows:
    //   1. Each input node appears once
    //   2. Each internal node has a B variable and an F variable
    //   3. Each output node appears once
    unsigned numberOfVariables = inputLayerSize + ( 2 * numberOfInternalNodes ) + outputLayerSize;
    printf( "Total number of variables: %u\n", numberOfVariables );

    inputQuery.setNumberOfVariables( (numberOfVariables) );

    // Next, we want to map each node to its corresponding
    // variables. We group variables according to this order: f's from
    // layer i, b's from layer i+1, and repeat.
    unsigned currentIndex = 0;
    for ( unsigned i = 1; i < numberOfLayers; ++i )
    {
        unsigned previousLayerSize = 2*(_acasNeuralNetwork.getLayerSize( i - 1 ));
        unsigned currentLayerSize = 2*(_acasNeuralNetwork.getLayerSize( i ));

        // First add the F variables from layer i-1
        for ( unsigned j = 0; j < previousLayerSize/2; ++j )
        {
            _nodeToF[NodeIndex( i - 1, j )] = (currentIndex);
            ++currentIndex;
        }
        for ( unsigned j = previousLayerSize/2; j < previousLayerSize; ++j )
        {
            _dupnodeToF[NodeIndex( i - 1, j )] = (currentIndex);
            ++currentIndex;
        }


        // Now add the B variables from layer i
        for ( unsigned j = 0; j < currentLayerSize/2; ++j )
        {
            _nodeToB[NodeIndex( i, j )] = (currentIndex+alpha);
            ++currentIndex;
        }
        for ( unsigned j = currentLayerSize/2; j < currentLayerSize; ++j )
        {
            _dupnodeToB[NodeIndex( i, j )] = (currentIndex+alpha);
            ++currentIndex;
        }
    }

    // Now we set the variable bounds. Input bounds are
    // given as part of the network. B variables are
    // unbounded, and F variables are non-negative.
    for ( unsigned i = 0; i < inputLayerSize; ++i )
    {
        double min, max;
        _acasNeuralNetwork.getInputRange( i, min, max );

        inputQuery.setLowerBound( _nodeToF[NodeIndex(0, i)], min );
        inputQuery.setUpperBound( _nodeToF[NodeIndex(0, i)], max );
    }

    for ( const auto &fNode : _nodeToF )
    {
        // Be careful not to override the bounds for the input layer
        if ( fNode.first._layer != 0 )
        {
            inputQuery.setLowerBound( fNode.second, 0.0 );
            inputQuery.setUpperBound( fNode.second, FloatUtils::infinity() );
        }
    }

    for ( const auto &fNode : _nodeToB )
    {
        inputQuery.setLowerBound( fNode.second, FloatUtils::negativeInfinity() );
        inputQuery.setUpperBound( fNode.second, FloatUtils::infinity() );
    }

    // Next come the actual equations
    for ( unsigned layer = 0; layer < numberOfLayers - 1; ++layer )
    {
        unsigned targetLayerSize = _acasNeuralNetwork.getLayerSize( layer + 1 );
        for ( unsigned target = 0; target < targetLayerSize; ++target )
        {
            // This will represent the equation:
            //   sum - b + fs = -bias
            Equation equation;

            // The b variable
            unsigned bVar = _nodeToB[NodeIndex(layer + 1, target)];
            equation.addAddend( -1.0, bVar );

            // The f variables from the previous layer
            for ( unsigned source = 0; source < _acasNeuralNetwork.getLayerSize( layer ); ++source )
            {
                unsigned fVar = _nodeToF[NodeIndex(layer, source)];
                equation.addAddend( _acasNeuralNetwork.getWeight( layer, source, target ), fVar );
            }

            // The bias
            equation.setScalar( -_acasNeuralNetwork.getBias( layer + 1, target ) );

            // Add the equation to the input query
            inputQuery.addEquation( equation );
        }
    }

    // Add the ReLU constraints
    for ( unsigned i = 1; i < numberOfLayers - 1; ++i )
    {
        unsigned currentLayerSize = _acasNeuralNetwork.getLayerSize( i );

        for ( unsigned j = 0; j < currentLayerSize; ++j )
        {
            unsigned b = _nodeToB[NodeIndex(i, j)];
            unsigned f = _nodeToF[NodeIndex(i, j)];
            PiecewiseLinearConstraint *relu = new ReluConstraint( b, f );

            inputQuery.addPiecewiseLinearConstraint( relu );
        }
    }

    // Mark the input and output variables
    for ( unsigned i = 0; i < inputLayerSize; ++i )
        inputQuery.markInputVariable( _nodeToF[NodeIndex( 0, (i) )], (i) );

    for ( unsigned i = 0; i < outputLayerSize; ++i )
        inputQuery.markOutputVariable( _nodeToB[NodeIndex( numberOfLayers - 1, (i) )], (i) );
}










void AcasParser::generateQueryAnagha( InputQuery &inputQuery )
{
    // First encode the actual network
    // _acasNeuralNetwork doesn't count the input layer, so add 1
    unsigned numberOfLayers = _acasNeuralNetwork.getNumLayers() + 1;
    inputLayerSize = (2*(_acasNeuralNetwork.getLayerSize( 0 )));
    outputLayerSize = (2*(_acasNeuralNetwork.getLayerSize( numberOfLayers - 1 )));
    //unsigned alpha = 22;
    unsigned numberOfInternalNodes = 0;
    for ( unsigned i = 1; i < numberOfLayers - 1; ++i )
        numberOfInternalNodes += 2*(_acasNeuralNetwork.getLayerSize( i ));

    printf( "Number of layers: %u. Input layer size: %u. Output layer size: %u. Number of ReLUs: %u\n",
            numberOfLayers, inputLayerSize, outputLayerSize, numberOfInternalNodes );

    // The total number of variables required for the encoding is computed as follows:
    //   1. Each input node appears once
    //   2. Each internal node has a B variable and an F variable
    //   3. Each output node appears once
    unsigned numberOfVariables = inputLayerSize + (2*numberOfInternalNodes) + outputLayerSize;
    printf( "Total number of variables: %u\n", numberOfVariables );

    inputQuery.setNumberOfVariables((numberOfVariables) );

    // Next, we want to map each node to its corresponding
    // variables. We group variables according to this order: f's from
    // layer i, b's from layer i+1, and repeat.
    unsigned currentIndex = 0;
    //unsigned _dupcurrentIndex = numberOfVariables/2;

    for ( unsigned i = 1; i < numberOfLayers; ++i )
    {
        unsigned previousLayerSize = 2*(_acasNeuralNetwork.getLayerSize( i - 1 ));
        printf("previousLayerSize = %u", previousLayerSize);
        unsigned currentLayerSize = 2*(_acasNeuralNetwork.getLayerSize( i ));
        printf("currentLayerSize = %u", currentLayerSize);
        // First add the F variables from layer i-1
        for ( unsigned j = 0; j < previousLayerSize; ++j )
        {
            _nodeToF[NodeIndex( i - 1, j )] = (currentIndex);
            ++currentIndex;
        }
        /*for ( unsigned j = previousLayerSize/2; j < previousLayerSize; ++j )
        {
            _dupnodeToF[NodeIndex( i - 1, j )] = (currentIndex);
            ++_dupcurrentIndex;
        }*/

        // Now add the B variables from layer i
        for ( unsigned j = 0; j < currentLayerSize; ++j )
        {
            _nodeToB[NodeIndex( i, j )] = (currentIndex);
            ++currentIndex;
        }
        /*for ( unsigned j = currentLayerSize/2; j < currentLayerSize; ++j )
        {
            _dupnodeToB[NodeIndex( i, j )] = (currentIndex);
            ++_dupcurrentIndex;
        }*/
    }

    // Now we set the variable bounds. Input bounds are
    // given as part of the network. B variables are
    // unbounded, and F variables are non-negative.
    for ( unsigned i = 0; i < inputLayerSize/2; ++i )
    {
        double min, max;
        _acasNeuralNetwork.getInputRange( i, min, max );

        inputQuery.setLowerBound( _nodeToF[NodeIndex(0, i)], min );
        inputQuery.setUpperBound( _nodeToF[NodeIndex(0, i)], max );

        inputQuery.setLowerBound( _nodeToF[NodeIndex(0, (i+(inputLayerSize/2)))], min );
        inputQuery.setUpperBound( _nodeToF[NodeIndex(0, (i+(inputLayerSize/2)))], max );
    }
    
    for ( const auto &fNode : _nodeToF )
    {
        // Be careful not to override the bounds for the input layer
        if ( fNode.first._layer != 0 )
        {
            inputQuery.setLowerBound( fNode.second, 0.0 );
            inputQuery.setUpperBound( fNode.second, FloatUtils::infinity() );
        }
    }
   
    for ( const auto &fNode : _nodeToB )
    {
        inputQuery.setLowerBound( fNode.second, FloatUtils::negativeInfinity() );
        inputQuery.setUpperBound( fNode.second, FloatUtils::infinity() );
    }
    
    // Next come the actual equations
    for ( unsigned layer = 0; layer < numberOfLayers - 1; ++layer )
    {
        unsigned targetLayerSize = _acasNeuralNetwork.getLayerSize( layer + 1 );
        for ( unsigned target = 0; target < targetLayerSize; ++target )
        {
            // This will represent the equation:
            //   sum - b + fs = -bias
            Equation equation;

            // The b variable
            unsigned bVar = _nodeToB[NodeIndex(layer + 1, target)];
            equation.addAddend( -1.0, bVar );

            // The f variables from the previous layer
            for ( unsigned source = 0; source < _acasNeuralNetwork.getLayerSize( layer ); ++source )
            {
                unsigned fVar = _nodeToF[NodeIndex(layer, source)];
                equation.addAddend( _acasNeuralNetwork.getWeight( layer, source, target ), fVar );
            }

            // The bias
            equation.setScalar( -_acasNeuralNetwork.getBias( layer + 1, target ) );

            // Add the equation to the input query
            inputQuery.addEquation( equation );
        }
    }

    for ( unsigned layer = 0; layer < numberOfLayers - 1; ++layer )
    {
        unsigned targetLayerSize = _acasNeuralNetwork.getLayerSize( layer + 1 );
        for ( unsigned target = targetLayerSize; target < 2*targetLayerSize; ++target )
        {
            // This will represent the equation:
            //   sum - b + fs = -bias
            Equation equation;

            // The b variable
            unsigned bVar = _nodeToB[NodeIndex(layer + 1, target)];
            equation.addAddend( -1.0, bVar );

            // The f variables from the previous layer
            for ( unsigned source = _acasNeuralNetwork.getLayerSize( layer ); source < 2*(_acasNeuralNetwork.getLayerSize( layer )); ++source )
            {
                unsigned fVar = _nodeToF[NodeIndex(layer, (source))];
                equation.addAddend( _acasNeuralNetwork.getWeight (layer, (source - _acasNeuralNetwork.getLayerSize( layer )), (target - _acasNeuralNetwork.getLayerSize( layer + 1 ))), fVar );
            }

            // The bias
            equation.setScalar( -_acasNeuralNetwork.getBias( layer + 1, (target - _acasNeuralNetwork.getLayerSize( layer + 1 )) ) );

            // Add the equation to the input query
            inputQuery.addEquation( equation );
        }
    }

    // Add the ReLU constraints
    for ( unsigned i = 1; i < numberOfLayers - 1; ++i )
    {
        unsigned currentLayerSize = _acasNeuralNetwork.getLayerSize( i );

        for ( unsigned j = 0; j < currentLayerSize; ++j )
        {
            unsigned b = _nodeToB[NodeIndex(i, j)];
            unsigned f = _nodeToF[NodeIndex(i, j)];
            PiecewiseLinearConstraint *relu = new ReluConstraint( b, f );
            inputQuery.addPiecewiseLinearConstraint( relu );
            
        }
    }

    for ( unsigned i = 1; i < numberOfLayers - 1; ++i )
    {
        unsigned currentLayerSize = (2*_acasNeuralNetwork.getLayerSize( i ));

        for ( unsigned j = currentLayerSize/2; j < currentLayerSize; ++j )
        {
            unsigned dupb = _nodeToB[NodeIndex(i, j)];
            unsigned dupf = _nodeToF[NodeIndex(i, j)];
            printf("dupb = %u, dupf = %u", dupb, dupf);
            PiecewiseLinearConstraint *duprelu = new ReluConstraint( dupb, dupf );
            inputQuery.addPiecewiseLinearConstraint( duprelu );
        }
    }

    // Mark the input and output variables
    for ( unsigned i = 0; i < inputLayerSize; ++i )
    {    
        inputQuery.markInputVariable( _nodeToF[NodeIndex( 0, (i) )], (i) );
        //inputQuery.markInputVariable( _dupnodeToF[NodeIndex( 0, (i) )], (i) );

    }
    /*for ( unsigned i = inputLayerSize/2; i < inputLayerSize; ++i )
    {    
        inputQuery.markInputVariable( _dupnodeToF[NodeIndex( 0, (i) )], (i) );
    }*/

    for ( unsigned i = 0; i < outputLayerSize; ++i )
    {
        inputQuery.markOutputVariable( _nodeToB[NodeIndex( numberOfLayers - 1, (i) )], (i) );
        //inputQuery.markOutputVariable( _dupnodeToB[NodeIndex( numberOfLayers - 1, (i) )], (i) );
    }
    /*for ( unsigned i = outputLayerSize/2; i < outputLayerSize; ++i )
    {
        inputQuery.markOutputVariable( _dupnodeToB[NodeIndex( numberOfLayers - 1, (i) )], (i) );
    }*/



    //_acasNeuralNetwork.incrementLayer();
    //unsigned int layerindex = _acasNeuralNetwork.getNumLayers();
  
}



  
    /*List<unsigned> outlist1;
    List<unsigned> outlist2;
    Set<unsigned> outSet1;
    Set<unsigned> outSet2;    
    unsigned counter = 0;
    for ( const auto &pair : inputQuery._outputIndexToVariable)
    {
        if(counter < outputLayerSize/2)
        {
            outlist1.append( pair.second );
            ++counter;
        }    
        else
        {
            outlist2.append(pair.second);
            ++counter;
        }
    }

    for (const auto &outVar1 : outlist1)
    {
        Set<unsigned> maxSet;
        for (const auto &outVar2 : outlist1)
        {
            if (&outVar1 != &outVar2)
            {
                maxSet.insert(outVar2);
                
            }
        }
        inputQuery.incrementNumberOfVariables();
        unsigned var = inputQuery.getNumberOfVariables();
        //inputQuery.setNumberOfVariables( (var+1) );
        //inputQuery.setLowerBound(var+1, 0.0);
        MaxConstraint *max = new MaxConstraint(var, maxSet);
        inputQuery.addPiecewiseLinearConstraint(max);
        printf("%u=max(",var);
        maxSet.print();
        printf(")\t");

        inputQuery.incrementNumberOfVariables();
        unsigned var2 = inputQuery.getNumberOfVariables();
        inputQuery.setNumberOfVariables( var2 );
        Equation equation2;
        equation2.addAddend(1, var2);
        equation2.addAddend(-1, outVar1);
        equation2.addAddend(1, var);
        equation2.setScalar(0);
        inputQuery.addEquation(equation2);
        
        outSet1.insert(var2);
        inputQuery.incrementNumberOfVariables();
        unsigned var_1 = inputQuery.getNumberOfVariables(); // var_1 softmax ka output hai
        inputQuery.setNumberOfVariables( var_1 );
        //NLR::NetworkLevelReasoner *nlr = inputQuery._networkLevelReasoner;;
        //nlr->setNeuronVariable( NLR::NeuronIndex( layerindex, var_1 ), var_1 );
        //nlr->setNeuronVariable( var_1, var_1 );
        inputQuery.markOutputVariable(var_1, var_1);
        //inputQuery.setLowerBound(var_1,0.0);
        //printf("%f",inputQuery.getLowerBound(var_1));
        SigmoidConstraint *sig = new SigmoidConstraint(var2,var_1);
        inputQuery.addTranscendentalConstraint(sig);    
    }



    
    inputQuery.incrementNumberOfVariables();
    unsigned var5 = inputQuery.getNumberOfVariables();
    //inputQuery.setLowerBound(var5, 0.0);
    //inputQuery.setLowerBound(var5, 1.0);
    MaxConstraint *maxout = new MaxConstraint(var5, outSet1);
    inputQuery.addPiecewiseLinearConstraint(maxout);

    inputQuery.incrementNumberOfVariables();
    unsigned conf = inputQuery.getNumberOfVariables();
    SigmoidConstraint *sig1 = new SigmoidConstraint(var5,conf);
    inputQuery.addTranscendentalConstraint(sig1);
    printf("\nConfidence1 = x%u(",conf);

    printf("\nmaxOut%u=max(",var5);
    outSet1.print();
    printf(")\t");
    
    for (const auto &outVar1 : outlist2)
    {
        Set<unsigned> maxSet;
        for (const auto &outVar2 : outlist2)
        {
            if (&outVar1 != &outVar2)
            {
                maxSet.insert(outVar2);
                
            }
        }
        inputQuery.incrementNumberOfVariables();
        unsigned var = inputQuery.getNumberOfVariables();
        MaxConstraint *max = new MaxConstraint(var, maxSet);
        inputQuery.addPiecewiseLinearConstraint(max);
        printf("%u=max(",var);
        maxSet.print();
        printf(")\t");

        inputQuery.incrementNumberOfVariables();
        unsigned var2 = inputQuery.getNumberOfVariables();
        Equation equation2;
        equation2.addAddend(1, var2);
        equation2.addAddend(-1, outVar1);
        equation2.addAddend(1, var);
        equation2.setScalar(0);
        inputQuery.addEquation(equation2);
        
        //outSet2.insert(var2);
        inputQuery.incrementNumberOfVariables();
        unsigned var_1 = inputQuery.getNumberOfVariables(); // var_1 softmax ka output hai
        //nlr->setNeuronVariable( NeuronIndex index, unsigned variable )
        outSet2.insert(var_1);
        SigmoidConstraint *sig = new SigmoidConstraint(var2,var_1);
        inputQuery.addTranscendentalConstraint(sig);

    }
    inputQuery.incrementNumberOfVariables();
    unsigned var6 = inputQuery.getNumberOfVariables();
    MaxConstraint *maxout2 = new MaxConstraint(var6, outSet2);
    inputQuery.addPiecewiseLinearConstraint(maxout2);

    inputQuery.incrementNumberOfVariables();
    
    unsigned conf2 = inputQuery.getNumberOfVariables();
    SigmoidConstraint *sig2 = new SigmoidConstraint(var6, conf2);
    inputQuery.addTranscendentalConstraint(sig2);
    printf("\nConfidence2 = x%u(",conf2);
    
    inputQuery.incrementNumberOfVariables();
    printf("\nmaxOut%u=max(",var6);
    outSet2.print();
    printf(")\t");
    

    
   

    unsigned counterX = 0;
    for ( const auto &pair : inputQuery._inputIndexToVariable)
    {
        if(counterX < inputLayerSize/2)
        {
            Equation equation4(Equation::LE);
            equation4.addAddend(1, (pair.second));
            equation4.addAddend(-1, (pair.second+(inputLayerSize/2)));
            equation4.setScalar(0.001);
            inputQuery.addEquation(equation4);
            printf("\n\ninput diff:\n");
            equation4.dump();
            printf("\n\n");
            printf("\n\nProperty Anagha: \nx%u - x%u\n", pair.second, (pair.second+(inputLayerSize/2)));
            ++counterX;
        }    
        else
        {
            ++counterX;
        } 
    }
   
    //printf("\n\nConfidence > k: \nx%u > 0.9999 and x%u > 0.9999 \n", conf, conf2);
    Equation equation5(Equation::GE);
    equation5.addAddend(1, var5);
    equation5.setScalar(0.6);
    //inputQuery.addEquation(equation5);
    printf("\n\nConf1:\n");
    //equation5.dump();

    Equation equation6(Equation::GE);
    equation6.addAddend(1, var6);
    equation6.setScalar(0.7);
    //inputQuery.addEquation(equation6);
    printf("\n\nConf2:\n");
    //equation6.dump();



    //std::list<unsigned>Super::iterator it;
    List<unsigned int>::iterator it = outlist1.begin();
    List<unsigned int>::iterator it2 = outlist2.begin();
    for(; ((it != outlist1.end()) && (it2 != outlist2.end())); ++it, ++it2)
    {
        Equation equation7(Equation::LE);
        equation7.addAddend(1, *it);
        equation7.addAddend(-1, *it2);
        equation7.setScalar(0.1);
        inputQuery.addEquation(equation7);

        printf("\n\nOutput difference %u - %u\n*", *it, *it2);
    }

*/



    



//
// Local Variables:
// compile-command: "make -C ../.. "
// tags-file-name: "../../TAGS"
// c-basic-offset: 4
// End:
//
