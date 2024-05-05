/*********************                                                        */
/*! \file Marabou.h
 ** \verbatim
 ** Top contributors (to current version):
 **   Guy Katz
 ** This file is part of the Marabou project.
 ** Copyright (c) 2017-2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved. See the file COPYING in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** [[ Add lengthier description here ]]

 **/

#ifndef __Marabou_h__
#define __Marabou_h__

#include "AcasParser.h"
#include "Engine.h"
#include "InputQuery.h"

class Marabou
{
public:
    Marabou();
    ~Marabou();

    /*
      Entry point of this class
    */
    void run();
    unsigned var5;
    unsigned var_6;

private:
    InputQuery _inputQuery;
    InputQuery _dupInputQuery;
    InputQuery iq;

    /*
      Extract the options and input files (network and property), and
      use them to generate the input query
    */
    void prepareInputQuery();
    void extractSplittingThreshold();

    unsigned sigmoid_reduced1(unsigned var2);
    unsigned sigmoid_reduced2(unsigned var2);
    unsigned sigmoid_reduced9(unsigned var2);
    unsigned sigmoid_reduced12(unsigned var2);
    unsigned sigmoid_reduced19(unsigned var2);
    unsigned sigmoid_anagha_final(unsigned aaa);
    /*
      Invoke the engine to solve the input query
    */
    void solveQuery();

    /*
      Display the results
    */
    void displayResults( unsigned long long microSecondsElapsed ) const;

    /*
      Export assignment as per Options
     */
    void exportAssignment() const;

    /*
      Import assignment for debugging as per Options
     */
    void importDebuggingSolution();

    /*
      ACAS network parser
    */
    AcasParser *_acasParser;

    /*
      The solver
    */
    Engine _engine;
};

#endif // __Marabou_h__

//
// Local Variables:
// compile-command: "make -C ../.. "
// tags-file-name: "../../TAGS"
// c-basic-offset: 4
// End:
//
