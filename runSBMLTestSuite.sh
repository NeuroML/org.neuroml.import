#!/bin/bash


if [ ! -d "sbml-test-suite" ]; then

    echo ""
    echo "A copy of the SBML Test Suite needs to be checked out from GitHub. Try:"
    echo ""
    echo "    git clone https://github.com/sbmlteam/sbml-test-suite.git"
    echo ""
    exit 1

fi

java -Xmx400M -classpath target/*jar-with-dependencies.jar org.neuroml.importer.sbml.SBMLImporter -runSBMLTestSuite


if [ "$?" -ne "0" ]; then
    echo ""
    echo "To run the SBML Test Suite, you need to create a jar in the target directory with all of the dependencies of org.neuroml.import included. Try:"
    echo ""
    echo "    mvn assembly:assembly -DdescriptorId=jar-with-dependencies"
    echo ""
    echo "You need to have org.neuroml.model, org.neuroml.export, etc. installed locally too."
    echo "A useful way to do this is with https://github.com/NeuroML/jNeuroML/blob/master/getNeuroML.py"
    echo ""
    exit 1
fi