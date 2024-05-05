/*********************                                                        */
/*! \file Marabou.cpp
 ** \verbatim
 ** Top contributors (to current version):
 **   Guy Katz
 ** This file is part of the Marabou project.
 ** Copyright (c) 2017-2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved. See the file COPYING in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief [[ Add one-line brief description here ]]
 **
 ** [[ Add lengthier description here ]]
 **/

#include "AcasParser.h"
#include "AutoFile.h"
#include "GlobalConfiguration.h"
#include "File.h"
#include "MStringf.h"
#include "Marabou.h"
#include "Options.h"
#include "PropertyParser.h"
#include "MarabouError.h"
#include "QueryLoader.h"
#include "Equation.h"

#ifdef _WIN32
#undef ERROR
#endif

Marabou::Marabou()
    : _acasParser(NULL), _engine()
{
}

Marabou::~Marabou()
{
    if (_acasParser)
    {
        delete _acasParser;
        _acasParser = NULL;
    }
}

void Marabou::run()
{
    struct timespec start = TimeUtils::sampleMicro();

    prepareInputQuery();
    solveQuery();

    struct timespec end = TimeUtils::sampleMicro();

    unsigned long long totalElapsed = TimeUtils::timePassed(start, end);
    displayResults(totalElapsed);

    if (Options::get()->getBool(Options::EXPORT_ASSIGNMENT))
        exportAssignment();
}

void Marabou::prepareInputQuery()
{
    String inputQueryFilePath = Options::get()->getString(Options::INPUT_QUERY_FILE_PATH);
    if (inputQueryFilePath.length() > 0)
    {
        /*
          Step 1: extract the query
        */
        if (!File::exists(inputQueryFilePath))
        {
            printf("Error: the specified inputQuery file (%s) doesn't exist!\n", inputQueryFilePath.ascii());
            throw MarabouError(MarabouError::FILE_DOESNT_EXIST, inputQueryFilePath.ascii());
        }

        printf("InputQuery: %s\n", inputQueryFilePath.ascii());
        _inputQuery = QueryLoader::loadQuery(inputQueryFilePath);

        printf("\n\ntest1\n\n");
        _inputQuery.dumpEquations();
        // anagha: add max() and sigmoid() to _inputQuery here

        //_inputQuery.addTranscendentalConstraint( new SigmoidConstraint( 5, 7 ) );

        _inputQuery.constructNetworkLevelReasoner();
    }
    else
    {
        /*
          Step 1: extract the network
        */
        String networkFilePath = Options::get()->getString(Options::INPUT_FILE_PATH);
        if (!File::exists(networkFilePath))
        {
            printf("Error: the specified network file (%s) doesn't exist!\n", networkFilePath.ascii());
            throw MarabouError(MarabouError::FILE_DOESNT_EXIST, networkFilePath.ascii());
        }
        printf("Network: %s\n", networkFilePath.ascii());

        // For now, assume the network is given in ACAS format
        _acasParser = new AcasParser(networkFilePath);
        _acasParser->generateQuery(_inputQuery);
        printf("\n\n\n\n**********Input query before any modifications begins**********");
        _inputQuery.dump();
        printf("**********Input query before any modifications ends**********\n\n\n\n");
        _inputQuery.constructNetworkLevelReasoner();

        NLR::NetworkLevelReasoner *nlr = _inputQuery.getNetworkLevelReasoner();


        _acasParser->generateQueryAnagha(iq);
        iq.constructNetworkLevelReasoner();

        //NLR::NetworkLevelReasoner *nlriq = iq.getNetworkLevelReasoner();
        printf("\n\n\niq dump\n");
        iq.dump();
        printf("\n\n\n\n");




        

        //List<Equation> &equations( _dupInputQuery.getEquations() );
        //for ( auto &equation : equations )
        //equation.updateVariableIndex( 0, 20);
    
        unsigned last_layer_index = ((nlr->getNumberOfLayers()) - 1);
        printf("Index of last layer: %u", last_layer_index);

        List<unsigned> outlist = _inputQuery.outputVars();
        Set<unsigned> outSet;

        for (const auto &outVar1 : outlist)
        {
            Set<unsigned> maxSet;
            for (const auto &outVar2 : outlist)
            {
                if (&outVar1 != &outVar2)
                {
                    maxSet.insert(outVar2);
                    
                }
            }
            _inputQuery.incrementNumberOfVariables();
            unsigned var = _inputQuery.getNumberOfVariables();
            Equation equation1;
            equation1.addAddend(1, var);
            equation1.addAddend(-1, outVar1);
            equation1.setScalar(0);
            /*_inputQuery.addEquation(equation1);*/

            MaxConstraint *max = new MaxConstraint(var, maxSet);
            _inputQuery.addPiecewiseLinearConstraint(max);
            printf("%u=max(",var);
            maxSet.print();
            printf(")\t");


            _inputQuery.incrementNumberOfVariables();
            unsigned var2 = _inputQuery.getNumberOfVariables();
            
            Equation equation2;
            equation2.addAddend(1, var2);
            equation2.addAddend(-1, outVar1);
            equation2.addAddend(1, var);
            equation2.setScalar(0);
            _inputQuery.addEquation(equation2);

           

           /*_inputQuery.incrementNumberOfVariables();
            unsigned var4 = _inputQuery.getNumberOfVariables();
            Equation equation4;
            equation4.addAddend(1, var4);
            //equation1.addAddend(-1, outVar1);
            equation4.setScalar(0);
            _inputQuery.addEquation(equation4);*/
            //_inputQuery.setLowerBound(var3,0);
            //_inputQuery.setUpperBound(var3,1.0);
            //printf("here is the lower bound%f",_inputQuery.getLowerBound(var3));
            //printf("the upper bound%f",_inputQuery.getUpperBound(var3));
            //printf("\t\t%u\t\t", var3);

            outSet.insert(var2);
            SigmoidConstraint *sig = new SigmoidConstraint(var2,outVar1);
            _inputQuery.addTranscendentalConstraint(sig);


            //_inputQuery.addTranscendentalConstraint( new SigmoidConstraint( 10, var ) );
            // printf("\t\t%u\t\t", outVar1);
            // maxSet.print();
        }
        _inputQuery.incrementNumberOfVariables();
        unsigned var5 = _inputQuery.getNumberOfVariables();
        MaxConstraint *maxout = new MaxConstraint(var5, outSet);
        _inputQuery.addPiecewiseLinearConstraint(maxout);


        printf("/nmaxOut%u=max(",var5);
        outSet.print();
        printf(")\t");


        _inputQuery.dump();

        //unsigned alpha = (_inputQuery.getNumberOfVariables() + 1); //index to add to every var in duplicate query
        //unsigned orig_no_of_vars = _inputQuery.getNumberOfVariables();
        unsigned alpha = (_inputQuery.getNumberOfVariables() + 1 );
        printf("\nalpha%u", alpha);
        _acasParser->generateQueryMod(_dupInputQuery, alpha);
        _dupInputQuery.constructNetworkLevelReasoner();
        NLR::NetworkLevelReasoner *dupnlr = _dupInputQuery.getNetworkLevelReasoner();

        printf("\n***dup query begins***\n");
        

        unsigned duplast_layer_index = ((dupnlr->getNumberOfLayers()) - 1);
        printf("Index of last layer: %u", duplast_layer_index);

        List<unsigned> dupoutlist = _dupInputQuery.outputVars();
        //_dupInputQuery.setNumberOfVariables(2*orig_no_of_vars);
        Set<unsigned> dupoutSet;
        for (const auto &dupoutVar1 : dupoutlist)
        {
            Set<unsigned> dupmaxSet;
            for (const auto &dupoutVar2 : dupoutlist)
            {
                if (&dupoutVar1 != &dupoutVar2)
                {
                    dupmaxSet.insert(dupoutVar2);
                    
                }
            }
            _dupInputQuery.incrementNumberOfVariables();
            unsigned dupvar = _dupInputQuery.getNumberOfVariables();
            /*Equation dupequation1;
            dupequation1.addAddend(1, dupvar);
            dupequation1.addAddend(-1, dupoutVar1);
            dupequation1.setScalar(0);
            _dupInputQuery.addEquation(dupequation1);*/

            MaxConstraint *dupmax = new MaxConstraint(dupvar, dupmaxSet);
            _dupInputQuery.addPiecewiseLinearConstraint(dupmax);
            printf("%u=max(",dupvar);
            dupmaxSet.print();
            printf(")\t");


            _dupInputQuery.incrementNumberOfVariables();
            
            unsigned dupvar2 = _dupInputQuery.getNumberOfVariables();
             
            Equation dupequation2;
            dupequation2.addAddend(1, dupvar2);
            dupequation2.addAddend(-1, dupoutVar1);
            dupequation2.addAddend(1, dupvar);
            dupequation2.setScalar(0);
            _dupInputQuery.addEquation(dupequation2);

            /*_dupInputQuery.incrementNumberOfVariables();
            unsigned dupvar4 = _dupInputQuery.getNumberOfVariables();
            Equation dupequation4;
            dupequation4.addAddend(1, dupvar4);
            dupequation4.setScalar(0);
            _dupInputQuery.addEquation(dupequation4);*/
            
            dupoutSet.insert(dupvar2);
            SigmoidConstraint *dupsig = new SigmoidConstraint(dupvar2,dupoutVar1);
            _dupInputQuery.addTranscendentalConstraint(dupsig);

            //_dupInputQuery.incrementNumberOfVariables();
            //unsigned dupvar5 = _dupInputQuery.getNumberOfVariables();
            

        }
        
        _dupInputQuery.incrementNumberOfVariables();
        unsigned dupvar5 = _dupInputQuery.getNumberOfVariables();
        MaxConstraint *dupmaxout = new MaxConstraint(dupvar5, dupoutSet);
        _dupInputQuery.addPiecewiseLinearConstraint(dupmaxout);


        printf("DupmaxOut%u=max(",dupvar5);
        dupoutSet.print();
        printf(")\t");

        _dupInputQuery.dump();
        printf("\nno of dup vars = %u\n", _dupInputQuery.getNumberOfVariables());
        printf("\n***dup query ends***\n");







        //printf("\t\tNo of vars%u", _inputQuery.getNumberOfVariables());
        //_inputQuery.incrementNumberOfVariables();
        //printf("\t\tNo of vars%u", _inputQuery.getNumberOfVariables());



        // double ub = _inputQuery.getUpperBound(5);
        // double lb = _inputQuery.getLowerBound(5);

        //_inputQuery.dumpOutputVars();
        //_inputQuery.dumpInputVars();
        //_inputQuery.dumpEquations();
        //printf("***********************************************************************************");
        //_inputQuery.dump();
        



        /*for duplicating query*/
        /*unsigned x = _inputQuery.getNumberOfVariables();
        _inputQuery.setNumberOfVariables(2*x);
        unsigned y = _inputQuery.getNumberOfVariables();
        printf("\n\n\n*******************%u\n*****************",y);
        List<unsigned> dup_input = _inputQuery.getInputVariables();
        for(const auto &ip : dup_input)
        {

        }*/


        /*_inputQuery.incrementNumberOfVariables();
        unsigned var5 = _inputQuery.getNumberOfVariables();
        _inputQuery.markInputVariable(var5, var5);*/


        /*unsigned var5;
        unsigned ii = 1;
        unsigned oo = _inputQuery.getNumInputVariables();
        while(ii <= oo)
        {
            _inputQuery.incrementNumberOfVariables();
            var5 = _inputQuery.getNumberOfVariables();
            _inputQuery.markInputVariable(var5, var5);
            ii++;
        }*/
        //_inputQuery.setLowerBound(22,-0.5);

        //_inputQuery.getNumOutputVariables();
        //_inputQuery.markOutputVariable(var5, var5);

        /*List<unsigned> dup_output = _inputQuery.getOutputVariables();
        printf("****no of output vars%u****",_inputQuery.getNumOutputVariables());
        unsigned var6;
        unsigned i = 1;
        unsigned o = _inputQuery.getNumOutputVariables();
        while(i <= o)
        {
            _inputQuery.incrementNumberOfVariables();
            var6 = _inputQuery.getNumberOfVariables();
            _inputQuery.markOutputVariable(var6, var6);

            i++;
            
        }*/

        /*unsigned var7;
        for (const auto &oq : _inputQuery._outputIndexToVariable)
        {
            _inputQuery.incrementNumberOfVariables();
            var7 = _inputQuery.getNumberOfVariables();
            _inputQuery.markOutputVariable(var7, var7);
            double var8 = _inputQuery.getLowerBound(oq.first);
            printf("%f", var8);
            //_inputQuery.setLowerBound(var7, var8);

            //oq.second
        }*/




        /*for duplicating query*/


        /*
          Step 2: extract the property in question
        */
        String propertyFilePath = Options::get()->getString(Options::PROPERTY_FILE_PATH);
        if (propertyFilePath != "")
        {
            printf("Property: %s\n", propertyFilePath.ascii()); // called
//            PropertyParser().parse(propertyFilePath, _inputQuery);
        }
        else
            printf("Property: None\n");

        printf("\n");
    }

    if (Options::get()->getBool(Options::DEBUG_ASSIGNMENT))
        importDebuggingSolution();

    String queryDumpFilePath = Options::get()->getString(Options::QUERY_DUMP_FILE);
    if (queryDumpFilePath.length() > 0)
    {
        _inputQuery.saveQuery(queryDumpFilePath);
        printf("\nInput query successfully dumped to file\n");
        exit(0);
    }
}

void Marabou::importDebuggingSolution()
{
    String fileName = Options::get()->getString(Options::IMPORT_ASSIGNMENT_FILE_PATH);
    AutoFile input(fileName);

    if (!IFile::exists(fileName))
    {
        throw MarabouError(MarabouError::FILE_DOES_NOT_EXIST, Stringf("File %s not found.\n", fileName.ascii()).ascii());
    }

    input->open(IFile::MODE_READ);

    unsigned numVars = atoi(input->readLine().trim().ascii());
    ASSERT(numVars == _inputQuery.getNumberOfVariables());

    unsigned var;
    double value;
    String line;

    // Import each assignment
    for (unsigned i = 0; i < numVars; ++i)
    {
        line = input->readLine();
        List<String> tokens = line.tokenize(",");
        auto it = tokens.begin();
        var = atoi(it->ascii());
        ASSERT(var == i);
        it++;
        value = atof(it->ascii());
        it++;
        ASSERT(it == tokens.end());
        _inputQuery.storeDebuggingSolution(var, value);
    }

    input->close();
}

void Marabou::exportAssignment() const
{
    String assignmentFileName = "assignment.txt";
    AutoFile exportFile(assignmentFileName);
    exportFile->open(IFile::MODE_WRITE_TRUNCATE);

    unsigned numberOfVariables = _inputQuery.getNumberOfVariables();
    // Number of Variables
    exportFile->write(Stringf("%u\n", numberOfVariables));

    // Export each assignment
    for (unsigned var = 0; var < numberOfVariables; ++var)
        exportFile->write(Stringf("%u, %f\n", var, _inputQuery.getSolutionValue(var)));

    exportFile->close();
}

void Marabou::solveQuery()
{
    if (_engine.processInputQuery(_inputQuery))
        _engine.solve(Options::get()->getInt(Options::TIMEOUT));

    if (_engine.getExitCode() == Engine::SAT)
        _engine.extractSolution(_inputQuery);
}

void Marabou::displayResults(unsigned long long microSecondsElapsed) const
{
    Engine::ExitCode result = _engine.getExitCode();
    String resultString;

    if (result == Engine::UNSAT)
    {
        resultString = "unsat";
        printf("unsat\n");
    }
    else if (result == Engine::SAT)
    {
        resultString = "sat";
        printf("sat\n");

        printf("Input assignment:\n");
        for (unsigned i = 0; i < _inputQuery.getNumInputVariables(); ++i)
            printf("\tx%u = %lf\n", i, _inputQuery.getSolutionValue(_inputQuery.inputVariableByIndex(i)));

        if (_inputQuery._networkLevelReasoner)
        {
            double *input = new double[_inputQuery.getNumInputVariables()];
            for (unsigned i = 0; i < _inputQuery.getNumInputVariables(); ++i)
                input[i] = _inputQuery.getSolutionValue(_inputQuery.inputVariableByIndex(i));

            NLR::NetworkLevelReasoner *nlr = _inputQuery._networkLevelReasoner;
            NLR::Layer *lastLayer = nlr->getLayer(nlr->getNumberOfLayers() - 1);
            double *output = new double[lastLayer->getSize()];

            nlr->evaluate(input, output);

            printf("\n");
            printf("Output:\n");
            for (unsigned i = 0; i < lastLayer->getSize(); ++i)
                printf("\ty%u = %lf\n", i, output[i]);
            printf("\n");
            delete[] input;
            delete[] output;
        }
        else
        {
            printf("\n");
            printf("Output:\n");
            for (unsigned i = 0; i < _inputQuery.getNumOutputVariables(); ++i)
                printf("\ty%u = %lf\n", i, _inputQuery.getSolutionValue(_inputQuery.outputVariableByIndex(i)));
            printf("\n");
        }
    }
    else if (result == Engine::TIMEOUT)
    {
        resultString = "TIMEOUT";
        printf("Timeout\n");
    }
    else if (result == Engine::ERROR)
    {
        resultString = "ERROR";
        printf("Error\n");
    }
    else
    {
        resultString = "UNKNOWN";
        printf("UNKNOWN EXIT CODE! (this should not happen)");
    }

    // Create a summary file, if requested
    String summaryFilePath = Options::get()->getString(Options::SUMMARY_FILE);
    if (summaryFilePath != "")
    {
        File summaryFile(summaryFilePath);
        summaryFile.open(File::MODE_WRITE_TRUNCATE);

        // Field #1: result
        summaryFile.write(resultString);

        // Field #2: total elapsed time
        summaryFile.write(Stringf(" %u ", microSecondsElapsed / 1000000)); // In seconds

        // Field #3: number of visited tree states
        summaryFile.write(Stringf("%u ",
                                  _engine.getStatistics()->getUnsignedAttribute(Statistics::NUM_VISITED_TREE_STATES)));

        // Field #4: average pivot time in micro seconds
        summaryFile.write(Stringf("%u",
                                  _engine.getStatistics()->getAveragePivotTimeInMicro()));

        summaryFile.write("\n");
    }
}

//
// Local Variables:
// compile-command: "make -C ../.. "
// tags-file-name: "../../TAGS"
// c-basic-offset: 4
// End:
//

