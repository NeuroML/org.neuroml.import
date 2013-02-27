package org.neuroml.importer;


import java.io.File;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import org.lemsml.jlems.type.Dimension;
import org.neuroml.model.IzhikevichCell;
import org.neuroml.model.NeuroMLDocument;

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

    public void testApp()
    {

    	System.out.println("Completed the Import test...");
    	
    }
    
    
    
}
