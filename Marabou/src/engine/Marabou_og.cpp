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
#include "Tableau.h"
#include "File.h"
#include "MStringf.h"
#include "Marabou.h"
#include "Options.h"
#include "PropertyParser.h"
#include "MarabouError.h"
#include "QueryLoader.h"
#include "Equation.h"
#include "MinConstraint.h"

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



unsigned Marabou::sigmoid_anagha_final(unsigned var2)
{
    Set<unsigned int> set3, set4, minSet3, minSet4;

    unsigned q2_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q2_ + 1);
    _inputQuery.setLowerBound(q2_, -100.0);
    _inputQuery.setUpperBound(q2_, 100.0);
    Equation eq2_;
    eq2_.addAddend(1, q2_);
    eq2_.addAddend(-0.002543865904564067, var2);
    eq2_.setScalar(0.01805060314480649);
    _inputQuery.addEquation(eq2_);
    set3.insert(q2_);
    unsigned q3_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q3_ + 1);
    _inputQuery.setLowerBound(q3_, -100.0);
    _inputQuery.setUpperBound(q3_, 100.0);
    Equation eq3_;
    eq3_.addAddend(1, q3_);
    eq3_.addAddend(-0.008671798219180675, var2);
    eq3_.setScalar(0.05003576561700055);
    _inputQuery.addEquation(eq3_);
    set3.insert(q3_);
    unsigned q4_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q4_ + 1);
    _inputQuery.setLowerBound(q4_, -100.0);
    _inputQuery.setUpperBound(q4_, 100.0);
    Equation eq4_;
    eq4_.addAddend(1, q4_);
    eq4_.addAddend(-0.0160918859304203, var2);
    eq4_.setScalar(0.08240230673319072);
    _inputQuery.addEquation(eq4_);
    set3.insert(q4_);
    unsigned q5_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q5_ + 1);
    _inputQuery.setLowerBound(q5_, -100.0);
    _inputQuery.setUpperBound(q5_, 100.0);
    Equation eq5_;
    eq5_.addAddend(1, q5_);
    eq5_.addAddend(-0.02434297406286057, var2);
    eq5_.setScalar(0.1143649740425666);
    _inputQuery.addEquation(eq5_);
    set3.insert(q5_);
    unsigned q6_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q6_ + 1);
    _inputQuery.setLowerBound(q6_, -100.0);
    _inputQuery.setUpperBound(q6_, 100.0);
    Equation eq6_;
    eq6_.addAddend(1, q6_);
    eq6_.addAddend(-0.03648513005820122, var2);
    eq6_.setScalar(0.15616581066719448);
    _inputQuery.addEquation(eq6_);
    set3.insert(q6_);
    
    unsigned q7_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q7_ + 1);
    _inputQuery.setLowerBound(q7_, -100.0);
    _inputQuery.setUpperBound(q7_, 100.0);
    Equation eq7_;
    eq7_.addAddend(1, q7_);
    eq7_.addAddend(-0.053930892473068974, var2);
    eq7_.setScalar(0.20870845479245959);
    _inputQuery.addEquation(eq7_);
    set3.insert(q7_);
    unsigned q8_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q8_ + 1);
    _inputQuery.setLowerBound(q8_, -100.0);
    _inputQuery.setUpperBound(q8_, 100.0);
    Equation eq8_;
    eq8_.addAddend(1, q8_);
    eq8_.addAddend(-0.0712222309563578, var2);
    eq8_.setScalar(0.25402251637804774);
    _inputQuery.addEquation(eq8_);
    set3.insert(q8_);
    unsigned q9_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q9_ + 1);
    _inputQuery.setLowerBound(q9_, -100.0);
    _inputQuery.setUpperBound(q9_, 100.0);
    Equation eq9_;
    eq9_.addAddend(1, q9_);
    eq9_.addAddend(-0.08517649902476777, var2);
    eq9_.setScalar(0.28712509843219053);
    _inputQuery.addEquation(eq9_);
    set3.insert(q9_);
    unsigned q10_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q10_ + 1);
    _inputQuery.setLowerBound(q10_, -100.0);
    _inputQuery.setUpperBound(q10_, 100.0);
    Equation eq10_;
    eq10_.addAddend(1, q10_);
    eq10_.addAddend(-0.10105713280442323, var2);
    eq10_.setScalar(0.32137391451541486);
    _inputQuery.addEquation(eq10_);
    set3.insert(q10_);
    unsigned q11_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q11_ + 1);
    _inputQuery.setLowerBound(q11_, -100.0);
    _inputQuery.setUpperBound(q11_, 100.0);
    Equation eq11_;
    eq11_.addAddend(1, q11_);
    eq11_.addAddend(-0.11877226533845178, var2);
    eq11_.setScalar(0.35576112149866723);
    _inputQuery.addEquation(eq11_);
    set3.insert(q11_);
    unsigned q12_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q12_ + 1);
    _inputQuery.setLowerBound(q12_, -100.0);
    _inputQuery.setUpperBound(q12_, 100.0);
    Equation eq12_;
    eq12_.addAddend(1, q12_);
    eq12_.addAddend(-0.13805388170120814, var2);
    eq12_.setScalar(0.38903490752059966);
    _inputQuery.addEquation(eq12_);
    set3.insert(q12_);
    unsigned q13_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q13_ + 1);
    _inputQuery.setLowerBound(q13_, -100.0);
    _inputQuery.setUpperBound(q13_, 100.0);
    Equation eq13_;
    eq13_.addAddend(1, q13_);
    eq13_.addAddend(-0.15841276484258954, var2);
    eq13_.setScalar(0.4197833653601326);
    _inputQuery.addEquation(eq13_);
    set3.insert(q13_);
    unsigned q14_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q14_ + 1);
    _inputQuery.setLowerBound(q14_, -100.0);
    _inputQuery.setUpperBound(q14_, 100.0);
    Equation eq14_;
    eq14_.addAddend(1, q14_);
    eq14_.addAddend(-0.17911003633714492, var2);
    eq14_.setScalar(0.44658826147689556);
    _inputQuery.addEquation(eq14_);
    set3.insert(q14_);
    unsigned q15_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q15_ + 1);
    _inputQuery.setLowerBound(q15_, -100.0);
    _inputQuery.setUpperBound(q15_, 100.0);
    Equation eq15_;
    eq15_.addAddend(1, q15_);
    eq15_.addAddend(-0.1991616624642956, var2);
    eq15_.setScalar(0.46824484877026984);
    _inputQuery.addEquation(eq15_);
    set3.insert(q15_);
    unsigned q16_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q16_ + 1);
    _inputQuery.setLowerBound(q16_, -100.0);
    _inputQuery.setUpperBound(q16_, 100.0);
    Equation eq16_;
    eq16_.addAddend(1, q16_);
    eq16_.addAddend(-0.217391566491753, var2);
    eq16_.setScalar(0.48401859038912826);
    _inputQuery.addEquation(eq16_);
    set3.insert(q16_);
    unsigned q17_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q17_ + 1);
    _inputQuery.setLowerBound(q17_, -100.0);
    _inputQuery.setUpperBound(q17_, 100.0);
    Equation eq17_;
    eq17_.addAddend(1, q17_);
    eq17_.addAddend(-0.23254124615360464, var2);
    eq17_.setScalar(0.49388113246220744);
    _inputQuery.addEquation(eq17_);
    set3.insert(q17_);
    unsigned q18_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q18_ + 1);
    _inputQuery.setLowerBound(q18_, -100.0);
    _inputQuery.setUpperBound(q18_, 100.0);
    Equation eq18_;
    eq18_.addAddend(1, q18_);
    eq18_.addAddend(-0.2465459176469087, var2);
    eq18_.setScalar(0.4996724297871519);
    _inputQuery.addEquation(eq18_);
    set3.insert(q18_);


    unsigned q101_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q101_ + 1);
    _inputQuery.setLowerBound(q101_, -100.0);
    _inputQuery.setUpperBound(q101_, 100.0);
    Equation eq101_;
    eq101_.addAddend(1, q101_);
    eq101_.addAddend(-0.2465459176469087, var2);
    eq101_.setScalar(0.5003275702128485);
    _inputQuery.addEquation(eq101_);
    unsigned q102_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q102_ + 1);
    _inputQuery.setLowerBound(q102_, -100.0);
    _inputQuery.setUpperBound(q102_, 100.0);
    Equation eq102_;
    eq102_.addAddend(1, q102_);
    eq102_.addAddend(-0.2251459482999449, var2);
    eq102_.setScalar(0.5101081756873811);
    _inputQuery.addEquation(eq102_);
    unsigned q103_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q103_ + 1);
    _inputQuery.setLowerBound(q103_, -100.0);
    _inputQuery.setUpperBound(q103_, 100.0);
    Equation eq103_;
    eq103_.addAddend(1, q103_);
    eq103_.addAddend(-0.1991616624642951, var2);
    eq103_.setScalar(0.5317551512297305);
    _inputQuery.addEquation(eq103_);
    unsigned q104_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q104_ + 1);
    _inputQuery.setLowerBound(q104_, -100.0);
    _inputQuery.setUpperBound(q104_, 100.0);
    Equation eq104_;
    eq104_.addAddend(1, q104_);
    eq104_.addAddend(-0.17911003633714528, var2);
    eq104_.setScalar(0.5534117385231038);
    _inputQuery.addEquation(eq104_);
    unsigned q105_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q105_ + 1);
    _inputQuery.setLowerBound(q105_, -100.0);
    _inputQuery.setUpperBound(q105_, 100.0);
    Equation eq105_;
    eq105_.addAddend(1, q105_);
    eq105_.addAddend(-0.15841276484258868, var2);
    eq105_.setScalar(0.5802166346398685);
    _inputQuery.addEquation(eq105_);
    unsigned q106_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q106_ + 1);
    _inputQuery.setLowerBound(q106_, -100.0);
    _inputQuery.setUpperBound(q106_, 100.0);
    Equation eq106_;
    eq106_.addAddend(1, q106_);
    eq106_.addAddend(-0.13805388170120766, var2);
    eq106_.setScalar(0.6109650924794012);
    _inputQuery.addEquation(eq106_);
    unsigned q107_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q107_ + 1);
    _inputQuery.setLowerBound(q107_, -100.0);
    _inputQuery.setUpperBound(q107_, 100.0);
    Equation eq107_;
    eq107_.addAddend(1, q107_);
    eq107_.addAddend(-0.11877226533845185, var2);
    eq107_.setScalar(0.6442388785013329);
    _inputQuery.addEquation(eq107_);
    unsigned q108_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q108_ + 1);
    _inputQuery.setLowerBound(q108_, -100.0);
    _inputQuery.setUpperBound(q108_, 100.0);
    Equation eq108_;
    eq108_.addAddend(1, q108_);
    eq108_.addAddend(-0.10983042222425944, var2);
    eq108_.setScalar(0.6606306476130153);
    _inputQuery.addEquation(eq108_);
    unsigned q109_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q109_ + 1);
    _inputQuery.setLowerBound(q109_, -100.0);
    _inputQuery.setUpperBound(q109_, 100.0);
    Equation eq109_;
    eq109_.addAddend(1, q109_);
    eq109_.addAddend(-0.040470623980445174, var2);
    eq109_.setScalar(0.8313480269902461);
    _inputQuery.addEquation(eq109_);
    unsigned q110_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q110_ + 1);
    _inputQuery.setLowerBound(q110_, -100.0);
    _inputQuery.setUpperBound(q110_, 100.0);
    Equation eq110_;
    eq110_.addAddend(1, q110_);
    eq110_.addAddend(-0.03229299282395175, var2);
    eq110_.setScalar(0.8578503415691678);
    _inputQuery.addEquation(eq110_);

    unsigned q111_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q111_ + 1);
    _inputQuery.setLowerBound(q111_, -100.0);
    _inputQuery.setUpperBound(q111_, 100.0);
    Equation eq111_;
    eq111_.addAddend(1, q111_);
    eq111_.addAddend(-0.025668092163887983, var2);
    eq111_.setScalar(0.8809372759083436);
    _inputQuery.addEquation(eq111_);
    unsigned q112_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q112_ + 1);
    _inputQuery.setLowerBound(q112_, -100.0);
    _inputQuery.setUpperBound(q112_, 100.0);
    Equation eq112_;
    eq112_.addAddend(1, q112_);
    eq112_.addAddend(-0.02033951725545943, var2);
    eq112_.setScalar(0.9008072241656719);
    _inputQuery.addEquation(eq112_);
    unsigned q113_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q113_ + 1);
    _inputQuery.setLowerBound(q113_, -100.0);
    _inputQuery.setUpperBound(q113_, 100.0);
    Equation eq113_;
    eq113_.addAddend(1, q113_);
    eq113_.addAddend(-0.016077813562956105, var2);
    eq113_.setScalar(0.9177391169334725);
    _inputQuery.addEquation(eq113_);
    unsigned q114_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q114_ + 1);
    _inputQuery.setLowerBound(q114_, -100.0);
    _inputQuery.setUpperBound(q114_, 100.0);
    Equation eq114_;
    eq114_.addAddend(1, q114_);
    eq114_.addAddend(-0.012684539525458829, var2);
    eq114_.setScalar(0.932049037105814);
    _inputQuery.addEquation(eq114_);
    unsigned q115_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q115_ + 1);
    _inputQuery.setLowerBound(q115_, -100.0);
    _inputQuery.setUpperBound(q115_, 100.0);
    Equation eq115_;
    eq115_.addAddend(1, q115_);
    eq115_.addAddend(-0.00890273748673438, var2);
    eq115_.setScalar(0.9490553646042065);
    _inputQuery.addEquation(eq115_);
    unsigned q116_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q116_ + 1);
    _inputQuery.setLowerBound(q116_, -100.0);
    _inputQuery.setUpperBound(q116_, 100.0);
    Equation eq116_;
    eq116_.addAddend(1, q116_);
    eq116_.addAddend(-0.005501670887814405, var2);
    eq116_.setScalar(0.9658501156231741);
    _inputQuery.addEquation(eq116_);
    unsigned q117_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q117_ + 1);
    _inputQuery.setLowerBound(q117_, -100.0);
    _inputQuery.setUpperBound(q117_, 100.0);
    Equation eq117_;
    eq117_.addAddend(1, q117_);
    eq117_.addAddend(-0.003390707677573013, var2);
    eq117_.setScalar(0.9773047373082717);
    _inputQuery.addEquation(eq117_);
    unsigned q118_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q118_ + 1);
    _inputQuery.setLowerBound(q118_, -100.0);
    _inputQuery.setUpperBound(q118_, 100.0);
    Equation eq118_;
    eq118_.addAddend(1, q118_);
    eq118_.addAddend(-0.0016654340025765723, var2);
    eq118_.setScalar(0.9876154660829327);
    _inputQuery.addEquation(eq118_);

    unsigned negative_q118_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q118_+1);
    _inputQuery.setLowerBound(negative_q118_,-100.0);
    _inputQuery.setUpperBound(negative_q118_,100.0);
    Equation negative_eq118_;
    negative_eq118_.addAddend(1,q118_);
    negative_eq118_.addAddend(1,negative_q118_);
    negative_eq118_.setScalar(0);
    _inputQuery.addEquation(negative_eq118_);
    set4.insert(negative_q118_);

    unsigned negative_q117_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q117_+1);
    _inputQuery.setLowerBound(negative_q117_,-100.0);
    _inputQuery.setUpperBound(negative_q117_,100.0);
    Equation negative_eq117_;
    negative_eq117_.addAddend(1,q117_);
    negative_eq117_.addAddend(1,negative_q117_);
    negative_eq117_.setScalar(0);
    _inputQuery.addEquation(negative_eq117_);
    set4.insert(negative_q117_);

    unsigned negative_q116_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q116_+1);
    _inputQuery.setLowerBound(negative_q116_,-100.0);
    _inputQuery.setUpperBound(negative_q116_,100.0);
    Equation negative_eq116_;
    negative_eq116_.addAddend(1,q116_);
    negative_eq116_.addAddend(1,negative_q116_);
    negative_eq116_.setScalar(0);
    _inputQuery.addEquation(negative_eq116_);
    set4.insert(negative_q116_);

    unsigned negative_q115_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q115_+1);
    _inputQuery.setLowerBound(negative_q115_,-100.0);
    _inputQuery.setUpperBound(negative_q115_,100.0);
    Equation negative_eq115_;
    negative_eq115_.addAddend(1,q115_);
    negative_eq115_.addAddend(1,negative_q115_);
    negative_eq115_.setScalar(0);
    _inputQuery.addEquation(negative_eq115_);
    set4.insert(negative_q115_);

    unsigned negative_q114_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q114_+1);
    _inputQuery.setLowerBound(negative_q114_,-100.0);
    _inputQuery.setUpperBound(negative_q114_,100.0);
    Equation negative_eq114_;
    negative_eq114_.addAddend(1,q114_);
    negative_eq114_.addAddend(1,negative_q114_);
    negative_eq114_.setScalar(0);
    _inputQuery.addEquation(negative_eq114_);
    set4.insert(negative_q114_);
    
    unsigned negative_q113_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q113_+1);
    _inputQuery.setLowerBound(negative_q113_,-100.0);
    _inputQuery.setUpperBound(negative_q113_,100.0);
    Equation negative_eq113_;
    negative_eq113_.addAddend(1,q113_);
    negative_eq113_.addAddend(1,negative_q113_);
    negative_eq113_.setScalar(0);
    _inputQuery.addEquation(negative_eq113_);
    set4.insert(negative_q113_);

    unsigned negative_q112_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q112_+1);
    _inputQuery.setLowerBound(negative_q112_,-100.0);
    _inputQuery.setUpperBound(negative_q112_,100.0);
    Equation negative_eq112_;
    negative_eq112_.addAddend(1,q112_);
    negative_eq112_.addAddend(1,negative_q112_);
    negative_eq112_.setScalar(0);
    _inputQuery.addEquation(negative_eq112_);
    set4.insert(negative_q112_);

    unsigned negative_q111_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q111_+1);
    _inputQuery.setLowerBound(negative_q111_,-100.0);
    _inputQuery.setUpperBound(negative_q111_,100.0);
    Equation negative_eq111_;
    negative_eq111_.addAddend(1,q111_);
    negative_eq111_.addAddend(1,negative_q111_);
    negative_eq111_.setScalar(0);
    _inputQuery.addEquation(negative_eq111_);
    set4.insert(negative_q111_);

    unsigned negative_q110_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q110_+1);
    _inputQuery.setLowerBound(negative_q110_,-100.0);
    _inputQuery.setUpperBound(negative_q110_,100.0);
    Equation negative_eq110_;
    negative_eq110_.addAddend(1,q110_);
    negative_eq110_.addAddend(1,negative_q110_);
    negative_eq110_.setScalar(0);
    _inputQuery.addEquation(negative_eq110_);
    set4.insert(negative_q110_);

    unsigned negative_q109_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q109_+1);
    _inputQuery.setLowerBound(negative_q109_,-100.0);
    _inputQuery.setUpperBound(negative_q109_,100.0);
    Equation negative_eq109_;
    negative_eq109_.addAddend(1,q109_);
    negative_eq109_.addAddend(1,negative_q109_);
    negative_eq109_.setScalar(0);
    _inputQuery.addEquation(negative_eq109_);
    set4.insert(negative_q109_);

    unsigned negative_q108_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q108_+1);
    _inputQuery.setLowerBound(negative_q108_,-100.0);
    _inputQuery.setUpperBound(negative_q108_,100.0);
    Equation negative_eq108_;
    negative_eq108_.addAddend(1,q108_);
    negative_eq108_.addAddend(1,negative_q108_);
    negative_eq108_.setScalar(0);
    _inputQuery.addEquation(negative_eq108_);
    set4.insert(negative_q108_);
    unsigned negative_q107_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q107_+1);
    _inputQuery.setLowerBound(negative_q107_,-100.0);
    _inputQuery.setUpperBound(negative_q107_,100.0);
    Equation negative_eq107_;
    negative_eq107_.addAddend(1,q107_);
    negative_eq107_.addAddend(1,negative_q107_);
    negative_eq107_.setScalar(0);
    _inputQuery.addEquation(negative_eq107_);
    set4.insert(negative_q107_);
    unsigned negative_q106_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q106_+1);
    _inputQuery.setLowerBound(negative_q106_,-100.0);
    _inputQuery.setUpperBound(negative_q106_,100.0);
    Equation negative_eq106_;
    negative_eq106_.addAddend(1,q106_);
    negative_eq106_.addAddend(1,negative_q106_);
    negative_eq106_.setScalar(0);
    _inputQuery.addEquation(negative_eq106_);
    set4.insert(negative_q106_);
    unsigned negative_q105_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q105_+1);
    _inputQuery.setLowerBound(negative_q105_,-100.0);
    _inputQuery.setUpperBound(negative_q105_,100.0);
    Equation negative_eq105_;
    negative_eq105_.addAddend(1,q105_);
    negative_eq105_.addAddend(1,negative_q105_);
    negative_eq105_.setScalar(0);
    _inputQuery.addEquation(negative_eq105_);
    set4.insert(negative_q105_);
    unsigned negative_q104_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q104_+1);
    _inputQuery.setLowerBound(negative_q104_,-100.0);
    _inputQuery.setUpperBound(negative_q104_,100.0);
    Equation negative_eq104_;
    negative_eq104_.addAddend(1,q104_);
    negative_eq104_.addAddend(1,negative_q104_);
    negative_eq104_.setScalar(0);
    _inputQuery.addEquation(negative_eq104_);
    set4.insert(negative_q104_);
    unsigned negative_q103_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q103_+1);
    _inputQuery.setLowerBound(negative_q103_,-100.0);
    _inputQuery.setUpperBound(negative_q103_,100.0);
    Equation negative_eq103_;
    negative_eq103_.addAddend(1,q103_);
    negative_eq103_.addAddend(1,negative_q103_);
    negative_eq103_.setScalar(0);
    _inputQuery.addEquation(negative_eq103_);
    set4.insert(negative_q103_);
    unsigned negative_q102_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q102_+1);
    _inputQuery.setLowerBound(negative_q102_,-100.0);
    _inputQuery.setUpperBound(negative_q102_,100.0);
    Equation negative_eq102_;
    negative_eq102_.addAddend(1,q103_);
    negative_eq102_.addAddend(1,negative_q102_);
    negative_eq102_.setScalar(0);
    _inputQuery.addEquation(negative_eq102_);
    set4.insert(negative_q102_);    
    unsigned negative_q101_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q101_+1);
    _inputQuery.setLowerBound(negative_q101_,-100.0);
    _inputQuery.setUpperBound(negative_q101_,100.0);
    Equation negative_eq101_;
    negative_eq101_.addAddend(1,q101_);
    negative_eq101_.addAddend(1,negative_q101_);
    negative_eq101_.setScalar(0);
    _inputQuery.addEquation(negative_eq101_);
    set4.insert(negative_q101_);

    unsigned q_at_x_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(q_at_x_+1);
    //0.70448665747663167881338441808708
    _inputQuery.setLowerBound(q_at_x_, -1.0);
    _inputQuery.setUpperBound(q_at_x_, -1.0);
    set4.insert(q_at_x_);


    unsigned negative_q_at_x_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_q_at_x_+1);
    _inputQuery.setLowerBound(negative_q_at_x_, -0.5);
    _inputQuery.setUpperBound(negative_q_at_x_, -0.5);
    printf("\nq_at_x_ = %d negative_q_at_x_ = %d", q_at_x_, negative_q_at_x_);

    unsigned zero_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(zero_+1);
    _inputQuery.setLowerBound(zero_, 0.0);
    _inputQuery.setUpperBound(zero_, 0.0);
    set3.insert(zero_);


    unsigned fmax3 = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(fmax3+1);
    _inputQuery.setLowerBound(fmax3,-100.0);
    _inputQuery.setUpperBound(fmax3,100.0);
    MaxConstraint *m3 = new MaxConstraint(fmax3, set3);
    _inputQuery.addPiecewiseLinearConstraint(m3);

    unsigned fmax4 = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(fmax4+1);
    _inputQuery.setLowerBound(fmax4,-100.0);
    _inputQuery.setUpperBound(fmax4,100.0);
    MaxConstraint *m4 = new MaxConstraint(fmax4, set4);
    _inputQuery.addPiecewiseLinearConstraint(m4);

    unsigned negative_fmax4 = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(negative_fmax4+1);
    _inputQuery.setLowerBound(negative_fmax4,-10.0);
    _inputQuery.setUpperBound(negative_fmax4,10.0);
    printf ("\n fmax3 = %d \tfmax4 = %d", fmax3,fmax4);
    Equation eq37_;
    eq37_.addAddend(1,fmax4);
    eq37_.addAddend(1,negative_fmax4);
    eq37_.setScalar(0);
    _inputQuery.addEquation(eq37_);

    unsigned negative_fmax3 = _inputQuery.getNumberOfVariables();
    printf("\tnegative_fmax1 = %u\t", negative_fmax3);
    _inputQuery.setNumberOfVariables(negative_fmax3+1);
    _inputQuery.setLowerBound(negative_fmax3,-100.0);
    _inputQuery.setUpperBound(negative_fmax3,100.0);
    Equation eq36_;
    eq36_.addAddend(1,fmax3);
    eq36_.addAddend(1,negative_fmax3);
    eq36_.setScalar(0);
    _inputQuery.addEquation(eq36_);

    unsigned variable_min3 = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(variable_min3+1);
    _inputQuery.setLowerBound(variable_min3, -100.0);
    _inputQuery.setUpperBound(variable_min3, 100.0);

    unsigned variable_min4 = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(variable_min4+1);
    _inputQuery.setLowerBound(variable_min4, -100.0);
    _inputQuery.setUpperBound(variable_min4, 100.0);
     
    unsigned point5 = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(point5+1);
    _inputQuery.setLowerBound(point5, 0.5);
    _inputQuery.setUpperBound(point5, 0.5);

    minSet3.insert(negative_q_at_x_);
    minSet3.insert(negative_fmax3);

    minSet4.insert(negative_fmax4);
    minSet4.insert(point5);

    unsigned neg_variable_min3 = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(neg_variable_min3+1);
    _inputQuery.setLowerBound(neg_variable_min3,-100.0);
    _inputQuery.setUpperBound(neg_variable_min3,100.0);

    unsigned neg_variable_min4 = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(neg_variable_min4+1);
    _inputQuery.setLowerBound(neg_variable_min4,-100.0);
    _inputQuery.setUpperBound(neg_variable_min4,100.0);

    MaxConstraint *min3 = new MaxConstraint(variable_min3, minSet3);
    _inputQuery.addPiecewiseLinearConstraint(min3);
    MaxConstraint *min4 = new MaxConstraint(variable_min4, minSet4);
    _inputQuery.addPiecewiseLinearConstraint(min4);

    printf("\n\nneg_variable_min3 = %d  neg_variable_min4 = %d\n\n", neg_variable_min3, neg_variable_min4);

    Equation eq38_;
    eq38_.addAddend(1,variable_min3);
    eq38_.addAddend(1,neg_variable_min3);
    eq38_.setScalar(0);
    _inputQuery.addEquation(eq38_);

    Equation eq39_;
    eq39_.addAddend(1,variable_min4);
    eq39_.addAddend(1,neg_variable_min4);
    eq39_.setScalar(0);
    _inputQuery.addEquation(eq39_);

    unsigned answer_ = _inputQuery.getNumberOfVariables();
    _inputQuery.setNumberOfVariables(answer_+1);
    Equation eq40_;
    eq40_.addAddend(1,answer_);
    eq40_.addAddend(-1,neg_variable_min3);
    eq40_.addAddend(-1,variable_min4);
    eq40_.setScalar(-0.5);
    _inputQuery.addEquation(eq40_);
    printf("\t\tanswer = %d\t\t", answer_);

    return answer_;
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
        _acasParser->generateQueryAnagha(_inputQuery);
        _inputQuery.constructNetworkLevelReasoner();
        
        
        /*List<unsigned> outlist, outlist1;
        outlist = _inputQuery.getOutputVariables();
        outlist.
        unsigned num_output_vars = _inputQuery.getNumOutputVariables();
        for(unsigned i=0; i<num_output_vars/2; i++)
        {
            outlist1 = outlist;
        }*/
        
        //_inputQuery.setLowerBound(12,1.0);
        //_inputQuery.setUpperBound(12,1.0);
        /*_inputQuery.setLowerBound(13,10.0);
        _inputQuery.setUpperBound(13,10.0);
        _inputQuery.setLowerBound(14,-10.0);
        _inputQuery.setUpperBound(14,-10.0);*/
        /*_inputQuery.setLowerBound(15,-10.0);
        _inputQuery.setUpperBound(15,-10.0);
        _inputQuery.setLowerBound(16,10.0);
        _inputQuery.setUpperBound(16,10.0);
        _inputQuery.setUpperBound(17,-10.0);
        _inputQuery.setLowerBound(17,-10.0);*/
        _inputQuery.setLowerBound(0,0.0);
        _inputQuery.setUpperBound(0,5.0);
        _inputQuery.setLowerBound(1,0.0);
        _inputQuery.setUpperBound(1,5.0);
        _inputQuery.setLowerBound(2,0.0);
        _inputQuery.setUpperBound(2,5.0);
        _inputQuery.setLowerBound(3,0.0);
        _inputQuery.setUpperBound(3,5.0);
      
        List<unsigned> outlist1;
        List<unsigned> outlist2;
        Set<unsigned> outSet1;
        Set<unsigned> outSet2;  
        //Set<unsigned> maxSet1, maxSet2;
        unsigned outputLayerSize = _inputQuery.getNumOutputVariables();  
        unsigned counter = 0;
        Map<unsigned,unsigned> id1;
        unsigned idCounter1 = 1;
        Map<unsigned,unsigned> id2;
        unsigned idCounter2 = 1;
        for ( const auto &pair : _inputQuery._outputIndexToVariable)
        {
            if(counter < outputLayerSize/2)
            {
                outlist1.append( pair.second );
                outSet1.insert( pair.second );
                id1.insert(pair.second, idCounter1);
                printf("\n<outlist1 = %d, idCounter1 = %d>", pair.second, idCounter1);
                ++idCounter1;
                ++counter;
            }    
            else
            {
                outlist2.append(pair.second);
                outSet2.insert(pair.second);

                id2.insert(pair.second, idCounter2);
                printf("\n<outlist2 = %d, idCounter2 = %d>", pair.second, idCounter2);
                ++idCounter2;
                ++counter;
            }
        }

        Map<unsigned, unsigned> map1;
        Set<unsigned> confSet1, confSet2;
        for (const auto &outVar1 : outlist1)
        {
            Set<unsigned> maxSet1;
            for (const auto &outVar2 : outlist1)
            {
                if (&outVar1 != &outVar2)
                {
                    maxSet1.insert(outVar2); 
                }
            }
            unsigned max_var1 = _inputQuery.getNumberOfVariables();
            _inputQuery.setNumberOfVariables(max_var1+1);
            MaxConstraint *max1 = new MaxConstraint(max_var1, maxSet1);
            _inputQuery.addPiecewiseLinearConstraint(max1);
            printf("\t\tmax_var1 = %d\t\t", max_var1);

            unsigned var2 = _inputQuery.getNumberOfVariables();
            _inputQuery.setNumberOfVariables(var2+1);
            Equation equation2;
            equation2.addAddend(1, var2);
            equation2.addAddend(-1, outVar1);
            equation2.addAddend(1, max_var1);
            equation2.setScalar(0);
            _inputQuery.addEquation(equation2);
            
            //outSet1.insert(var2);
            //_inputQuery.markOutputVariable(var_1, var_1);
            //SigmoidConstraint *sig = new SigmoidConstraint(var2,var_1);
            printf("\n\noutVar1 = %d max_var1 = %d\n\n", outVar1, max_var1);
            printf("\t\tvar2(before sigmoid fn) = %d\t\t", var2);
            unsigned conf1 = _inputQuery.getNumberOfVariables();
            _inputQuery.setNumberOfVariables(conf1+1);
            _inputQuery.setUpperBound(conf1,-1.0);
            _inputQuery.setLowerBound(conf1,1.0);
            conf1 = sigmoid_anagha_final(var2); 
            confSet1.insert(conf1);   
            map1.insert(conf1, id1.get(outVar1));
            

            printf("\nconfidence = %d", conf1);
           
        }
        Map<unsigned, unsigned> map2;
        for (const auto &outVar1 : outlist2)
        {
            Set<unsigned> maxSet2;
            for (const auto &outVar2 : outlist2)
            {
                if (&outVar1 != &outVar2)
                {
                    maxSet2.insert(outVar2); 
                }
            }
            unsigned max_var2 = _inputQuery.getNumberOfVariables();
            _inputQuery.setNumberOfVariables(max_var2+1);
            MaxConstraint *max2 = new MaxConstraint(max_var2, maxSet2);
            _inputQuery.addPiecewiseLinearConstraint(max2);
            printf("\t\tmax_var2 = %d\t\t", max_var2);
            unsigned var3 = _inputQuery.getNumberOfVariables();
            _inputQuery.setNumberOfVariables(var3+1);
            //_inputQuery.setUpperBound(var3,-1000.0);
            //_inputQuery.setLowerBound(var3,1000.0);
            Equation equation3;
            equation3.addAddend(1, var3);
            equation3.addAddend(-1, outVar1);
            equation3.addAddend(1, max_var2);
            equation3.setScalar(0);
            _inputQuery.addEquation(equation3);
            
            //outSet1.insert(var2);
            //_inputQuery.markOutputVariable(var_1, var_1);
            //SigmoidConstraint *sig = new SigmoidConstraint(var2,var_1);
            unsigned conf2 = _inputQuery.getNumberOfVariables();
            _inputQuery.setNumberOfVariables(conf2+1);
            _inputQuery.setUpperBound(conf2,-1.0);
            _inputQuery.setLowerBound(conf2,1.0);
            printf("\t\tvar3(before sigmoid fn) = %d\t\t", var3);
            conf2 = sigmoid_anagha_final(var3);    
            confSet2.insert(conf2);   
            map2.insert(conf2, id2.get(outVar1));

            printf("\nconfidence = %d", conf2);
            
        }
        unsigned max_conf1 = _inputQuery.getNumberOfVariables();
        _inputQuery.setNumberOfVariables(max_conf1+1);
        MaxConstraint *maxConfidence1 = new MaxConstraint(max_conf1, confSet1);
        _inputQuery.addPiecewiseLinearConstraint(maxConfidence1);
        printf("\nMaxConfidence1 = %d", max_conf1);
        unsigned max_conf2 = _inputQuery.getNumberOfVariables();
        _inputQuery.setNumberOfVariables(max_conf2+1);
        MaxConstraint *maxConfidence2 = new MaxConstraint(max_conf2, confSet2);
        _inputQuery.addPiecewiseLinearConstraint(maxConfidence2);
        printf("\nMaxConfidence2 = %d", max_conf2);

        unsigned counterX = 0;
        for ( const auto &pair : _inputQuery._inputIndexToVariable)
        {
            if(counterX < (_inputQuery.getNumInputVariables())/2)
            {
                Equation equation4(Equation::LE);
                equation4.addAddend(1, (pair.second));
                equation4.addAddend(-1, (pair.second+(_inputQuery.getNumInputVariables()/2)));
                equation4.setScalar(0.0001);
                _inputQuery.addEquation(equation4);
                printf("\n\ninput diff:\n");
                equation4.dump();
                printf("\n\n");
                printf("\n\nProperty Anagha: \nx%u - x%u\n", pair.second, (pair.second+(_inputQuery.getNumInputVariables()/2)));
                ++counterX;
            }    
            else
            {
                ++counterX;
            } 
        }

        //printf("\n\nConfidence > k: \nx%u > 0.9999 and x%u > 0.9999 \n", conf, conf2);
        Equation equation5(Equation::GE);
        equation5.addAddend(1, max_conf1);
        equation5.setScalar(0.9);
        _inputQuery.addEquation(equation5);
        

        Equation equation6(Equation::GE);
        equation6.addAddend(1, max_conf2);
        equation6.setScalar(0.9);
        _inputQuery.addEquation(equation6);
        

/*        List<unsigned int>::iterator it = outlist1.begin();
        List<unsigned int>::iterator it2 = outlist2.begin();
        for(; ((it != outlist1.end()) && (it2 != outlist2.end())); ++it, ++it2)
        {
            Equation equation7(Equation::GE);
            equation7.addAddend(1, *it);
            equation7.addAddend(-1, *it2);
            equation7.setScalar(0.00001);
            _inputQuery.addEquation(equation7);

            printf("\n\nOutput difference %u - %u\n*", *it, *it2);
        }
*/

	unsigned out1 = _inputQuery.getNumberOfVariables();
        _inputQuery.setNumberOfVariables(out1+1);
//	unsigned xxx = outlist1.front();
//        _inputQuery.getUpperBound(outlist1.front());
	_inputQuery.setUpperBound(out1, 100);
        _inputQuery.setLowerBound(out1, -100);
//	_inputQuery.setUpperBound(out1, _inputQuery.getUpperBound(outlist1.front()));
//        _inputQuery.setLowerBound(out1, _inputQuery.getLowerBound(outlist1.front()));
	MaxConstraint *out_1 = new MaxConstraint(out1, outSet1);
        _inputQuery.addPiecewiseLinearConstraint(out_1);

        unsigned out2 = _inputQuery.getNumberOfVariables();
        _inputQuery.setNumberOfVariables(out2+1);
//	_inputQuery.setUpperBound(out2, _inputQuery.getUpperBound(outlist2.front()));
//        _inputQuery.setLowerBound(out2, _inputQuery.getLowerBound(outlist2.front()));
	 _inputQuery.setUpperBound(out2, 100);
        _inputQuery.setLowerBound(out2, -100);  

      MaxConstraint *out_2 = new MaxConstraint(out2, outSet2);
        _inputQuery.addPiecewiseLinearConstraint(out_2);

        Equation equation7(Equation::GE);
        equation7.addAddend(1, out1);
        equation7.addAddend(-1, out2);
        equation7.setScalar(0.01);
        _inputQuery.addEquation(equation7); 

/*        
        unsigned var7 = _inputQuery.getNumberOfVariables();
        _inputQuery.setNumberOfVariables(var7+1);
        _inputQuery.setUpperBound(var7,-3.0);
        _inputQuery.setLowerBound(var7,-3.0);
        printf("sig(-3) = %d", var7);
        //sigmoid_anagha_test(var7);
        //sigmoid_anagha_test2(var7);
        sigmoid_anagha_final(var7);

        unsigned var5 = _inputQuery.getNumberOfVariables();
        _inputQuery.setNumberOfVariables(var5+1);
        _inputQuery.setUpperBound(var5,-1.0);
        _inputQuery.setLowerBound(var5,-1.0);
        printf("sig(-1) = %d", var5);
        //sigmoid_anagha_test(var5);
        //sigmoid_anagha_test2(var5);
        sigmoid_anagha_final(var5);

        unsigned var3 = _inputQuery.getNumberOfVariables();
        _inputQuery.setNumberOfVariables(var3+1);
        _inputQuery.setUpperBound(var3,0.0);
        _inputQuery.setLowerBound(var3,0.0);
        printf("sig(0) = %d", var3);
        //sigmoid_anagha_test(var3);
        //sigmoid_anagha_test2(var3);
        sigmoid_anagha_final(var3);

        unsigned var4 = _inputQuery.getNumberOfVariables();
        _inputQuery.setNumberOfVariables(var4+1);
        _inputQuery.setUpperBound(var4,1.0);
        _inputQuery.setLowerBound(var4,1.0);
        printf("sig(1) = %d", var4);
        //sigmoid_anagha_test(var4);
        //sigmoid_anagha_test2(var4);
        sigmoid_anagha_final(var4);
 
        unsigned var8 = _inputQuery.getNumberOfVariables();
        _inputQuery.setNumberOfVariables(var8+1);
        _inputQuery.setUpperBound(var8,2.0);
        _inputQuery.setLowerBound(var8,2.0);
        printf("sig(2) = %d", var8);
        //sigmoid_anagha_test(var8);
        //sigmoid_anagha_test2(var8);
        sigmoid_anagha_final(var8);

        unsigned var6 = _inputQuery.getNumberOfVariables();
        _inputQuery.setNumberOfVariables(var6+1);
        _inputQuery.setUpperBound(var6,3.0);
        _inputQuery.setLowerBound(var6,3.0);
        printf("sig(3) = %d", var6);
        //sigmoid_anagha_test(var6);
        //sigmoid_anagha_test2(var6);
        sigmoid_anagha_final(var6);

        unsigned var9 = _inputQuery.getNumberOfVariables();
        _inputQuery.setNumberOfVariables(var9+1);
        _inputQuery.setUpperBound(var9,4.0);
        _inputQuery.setLowerBound(var9,4.0);
        printf("sig(4) = %d", var9);
        //sigmoid_anagha_test(var9);
        //sigmoid_anagha_test2(var9);
        sigmoid_anagha_final(var9);

        unsigned var10 = _inputQuery.getNumberOfVariables();
        _inputQuery.setNumberOfVariables(var10+1);
        _inputQuery.setUpperBound(var10,5.0);
        _inputQuery.setLowerBound(var10,5.0);
        printf("sig(5) = %d", var10);
        //sigmoid_anagha_test(var10);
        //sigmoid_anagha_test2(var10);
        sigmoid_anagha_final(var10);
*/

         


        /*
          Step 2: extract the property in question
        */
        String propertyFilePath = Options::get()->getString(Options::PROPERTY_FILE_PATH);
        if (propertyFilePath != "")
        {
            printf("Property: %s\n", propertyFilePath.ascii()); // called
            //PropertyParser().parse(propertyFilePath, _inputQuery);
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

            
            for(unsigned i = 0; i < _inputQuery.getNumberOfVariables(); i++)
            {
                double xy = _inputQuery.getSolutionValue(i);
                printf("\nSolution value x%u = %f", i, xy);   
            }

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
