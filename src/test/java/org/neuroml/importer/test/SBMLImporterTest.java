package org.neuroml.importer.test;

import java.io.File;
import java.io.IOException;

import javax.xml.stream.XMLStreamException;

import org.lemsml.jlems.core.expression.ParseError;
import org.lemsml.jlems.core.logging.E;
import org.lemsml.jlems.core.logging.MinimalMessageHandler;
import org.lemsml.jlems.core.run.ConnectionError;
import org.lemsml.jlems.core.run.RuntimeError;
import org.lemsml.jlems.core.sim.ContentError;
import org.lemsml.jlems.core.sim.ParseException;
import org.lemsml.jlems.core.type.BuildException;
import org.lemsml.jlems.core.type.Lems;
import org.lemsml.jlems.core.xml.XMLException;
import org.lemsml.jlems.io.util.FileUtil;
import org.lemsml.jlems.io.xmlio.XMLSerializer;
import org.neuroml.export.Utils;
import org.neuroml.importer.sbml.SBMLImporter;
import org.sbml.jsbml.SBMLException;

import junit.framework.TestCase;

public class SBMLImporterTest extends TestCase {

    String exampleSrcDir = "src/test/resources";
    /*
    if (sbmlFile.getName().indexOf("Izh")>=0) len = 140;
    if (sbmlFile.getName().indexOf("0039")>=0) len = 50;
    if (sbmlFile.getName().indexOf("00118")>=0) len = 50;
    if (sbmlFile.getName().indexOf("00184")>=0) len = 1000;*/

	public void testSimple()  throws ContentError, ParseError, IOException, ParseException, BuildException, XMLException, XMLStreamException, SBMLException, org.sbml.jsbml.text.parser.ParseException, ConnectionError, RuntimeError {

        File sbmlFile = new File(exampleSrcDir+"/Simple3Species.xml");
        convertSBMLtoLEMSFile(sbmlFile, 8);
	}

	public void testIzhikevich()  throws ContentError, ParseError, IOException, ParseException, BuildException, XMLException, XMLStreamException, SBMLException, org.sbml.jsbml.text.parser.ParseException, ConnectionError, RuntimeError {

        File sbmlFile = new File(exampleSrcDir+"/Izhikevich.xml");
        convertSBMLtoLEMSFile(sbmlFile, 150);
	}

	public void testBio39()  throws ContentError, ParseError, IOException, ParseException, BuildException, XMLException, XMLStreamException, SBMLException, org.sbml.jsbml.text.parser.ParseException, ConnectionError, RuntimeError {

        File sbmlFile = new File(exampleSrcDir+"/BIOMD0000000039.xml");
        convertSBMLtoLEMSFile(sbmlFile, 50);
	}
	
		
	public void testBio184()  throws ContentError, ParseError, IOException, ParseException, BuildException, XMLException, XMLStreamException, SBMLException, org.sbml.jsbml.text.parser.ParseException, ConnectionError, RuntimeError {

        File sbmlFile = new File(exampleSrcDir+"/BIOMD0000000184.xml");
        convertSBMLtoLEMSFile(sbmlFile, 1000);
	}
	
	/*
	public void testBio118()  throws ContentError, ParseError, IOException, ParseException, BuildException, XMLException, XMLStreamException, SBMLException, org.sbml.jsbml.text.parser.ParseException {

        File sbmlFile = new File(exampleSrcDir+"/BIOMD0000000118.xml");
        convertSBMLtoLEMSFile(sbmlFile, 50);
	}*/

	public void convertSBMLtoLEMSFile(File sbmlFile, float simDuration) throws ContentError, ParseError, IOException, ParseException, BuildException, XMLException, XMLStreamException, SBMLException, org.sbml.jsbml.text.parser.ParseException, ConnectionError, RuntimeError {

    	MinimalMessageHandler.setVeryMinimal(true);
		
        System.out.println( "Converting "+sbmlFile.getAbsolutePath()+" from SBML to LEMS" );
        float dt = (float)(simDuration/10000.0);

        Lems lems = SBMLImporter.convertSBMLToLEMS(sbmlFile, simDuration, dt);
        lems.resolve();
        String lemsString  = XMLSerializer.serialize(lems);
    
        //E.info("Created: \n"+lemsString);
        //E.info("Info: \n"+lems.textSummary());

        File testFile = new File(sbmlFile.getParentFile(), sbmlFile.getName().replaceAll(".xml", "")+"_LEMS.xml");

        FileUtil.writeStringToFile(lemsString, testFile);

        E.info("Written to: "+ testFile.getCanonicalPath());

        E.info("Loading LEMS file from: "+ testFile.getAbsolutePath());

		Lems lems2 = Utils.readLemsNeuroMLFile(testFile).getLems();
		lems2.resolve();
	}
		
		

}
