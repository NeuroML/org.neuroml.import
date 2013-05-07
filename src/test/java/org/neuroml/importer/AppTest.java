package org.neuroml.importer;



import java.io.File;
import java.io.IOException;

import org.lemsml.jlems.core.xml.XMLElementReader;
import org.lemsml.jlems.io.util.FileUtil;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;


/**
 * Unit test for simple App.
 */
public class AppTest extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AppTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
    }

    public void testVersions() throws IOException
    {
    	System.out.println("Running a test on version usage, making all references to versions are: v"+Main.ORG_NEUROML_IMPORT_VERSION+"...");

    	String jnmlPom = FileUtil.readStringFromFile(new File("pom.xml"));

    	XMLElementReader xer = new XMLElementReader(jnmlPom);
    	assertEquals(Main.ORG_NEUROML_IMPORT_VERSION, xer.getRootElement().getElement("version").getBody());
    	
    }

    public void testApp()
    {

    	System.out.println("Completed the Import test...");
    	
    }
    
    
    
}
