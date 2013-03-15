package org.neuroml.importer.sbml;



import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import javax.xml.stream.XMLStreamException;

import org.lemsml.jlems.core.expression.ParseError;
import org.lemsml.jlems.core.logging.E;
import org.lemsml.jlems.core.sim.ContentError;
import org.lemsml.jlems.core.sim.Sim;
import org.lemsml.jlems.core.type.BuildException;
import org.lemsml.jlems.core.type.Component;
import org.lemsml.jlems.core.type.ComponentType;
import org.lemsml.jlems.core.type.Constant;
import org.lemsml.jlems.core.type.Dimension;
import org.lemsml.jlems.core.type.Exposure;
import org.lemsml.jlems.core.type.Lems;
import org.lemsml.jlems.core.type.Target;
import org.lemsml.jlems.core.type.dynamics.DerivedVariable;
import org.lemsml.jlems.core.type.dynamics.Dynamics;
import org.lemsml.jlems.core.type.dynamics.OnCondition;
import org.lemsml.jlems.core.type.dynamics.OnStart;
import org.lemsml.jlems.core.type.dynamics.StateAssignment;
import org.lemsml.jlems.core.type.dynamics.StateVariable;
import org.lemsml.jlems.core.type.dynamics.TimeDerivative;
import org.lemsml.jlems.core.xml.XMLException;
import org.lemsml.jlems.io.logging.DefaultLogger;
import org.lemsml.jlems.io.out.FileResultWriterFactory;
import org.lemsml.jlems.io.util.FileUtil;
import org.lemsml.jlems.io.xmlio.XMLSerializer;
import org.lemsml.jlems.viz.datadisplay.SwingDataViewerFactory;
import org.lemsml.jlems.viz.plot.ColorUtil;
import org.neuroml.export.Utils;
import org.sbml.jsbml.ASTNode;
import org.sbml.jsbml.AssignmentRule;
import org.sbml.jsbml.Compartment;
import org.sbml.jsbml.Event;
import org.sbml.jsbml.EventAssignment;
import org.sbml.jsbml.FunctionDefinition;
import org.sbml.jsbml.JSBML;
import org.sbml.jsbml.KineticLaw;
import org.sbml.jsbml.LocalParameter;
import org.sbml.jsbml.Model;
import org.sbml.jsbml.Parameter;
import org.sbml.jsbml.RateRule;
import org.sbml.jsbml.Reaction;
import org.sbml.jsbml.Rule;
import org.sbml.jsbml.SBMLDocument;
import org.sbml.jsbml.SBMLException;
import org.sbml.jsbml.SBMLReader;
import org.sbml.jsbml.Species;
import org.sbml.jsbml.SpeciesReference;
import org.sbml.jsbml.text.parser.ParseException;

public class SBMLImporter  {


    public SBMLImporter() {
        E.info("Created new SBMLImporter...");
    }

    public static String tscaleName = "tscale";

    public static Lems convertSBMLToLEMS(File sbmlFile, float simDuration, float simDt) throws ContentError, XMLStreamException, ParseError, org.lemsml.jlems.core.sim.ParseException, BuildException, XMLException, IOException {
    	return convertSBMLToLEMS(sbmlFile, simDuration, simDt, null);
    }

    @SuppressWarnings("deprecation")
	public static Lems convertSBMLToLEMS(File sbmlFile, float simDuration, float simDt, File dirForResults) throws ContentError, XMLStreamException, ParseError, org.lemsml.jlems.core.sim.ParseException, BuildException, XMLException, IOException {

    	E.setDebug(false);
        SBMLReader sr = new SBMLReader();

        SBMLDocument doc = sr.readSBML(sbmlFile);

        Model model = doc.getModel();

        E.info("Read in SBML from "+sbmlFile.getAbsolutePath());

        //TODO: update!!!
        File f = new File("../NeuroML2/NeuroML2CoreTypes/Simulation.xml");
        //File f = new File("NeuroML2CoreTypes/Cells.xml");

        E.info("Loading LEMS file from: "+ f.getAbsolutePath());


		Lems lems = Utils.loadLemsFile(f);

        boolean addModel = true;

        ComponentType ct = new ComponentType(model.getId());
        if (addModel)
            lems.addComponentType(ct);

        Component comp = new Component(model.getId()+"_0", ct);
        if (addModel)
            lems.addComponent(comp);

        Dimension noDim = new Dimension(Dimension.NO_DIMENSION);

        Dynamics b = new Dynamics();
        ct.dynamicses.add(b);

        OnStart os = new OnStart();

        for(Compartment c: model.getListOfCompartments()){
            Dimension compDim = noDim;
            String size = c.getSize()+"";
            if (!c.isSetSize())
            	size="1";
            
            E.info("Adding: "+c+" (size = "+size+" (set? "+c.isSetSize()+"), constant = "+c.isConstant()+")");

            if (c.isConstant()){
                Constant constComp = new Constant(c.getId(), compDim, size);
                ct.constants.add(constComp);
            } else {
                Exposure ex = new Exposure(c.getId(), compDim);
                ct.exposures.add(ex);

                StateVariable sv = new StateVariable(c.getId(), compDim, ex);
                b.stateVariables.add(sv);

                StateAssignment sa = new StateAssignment(c.getId(), c.getSize()+"");
                os.stateAssignments.add(sa);
            }


        }

        for(Parameter p: model.getListOfParameters()) {
            //org.neuroml.lems.type.Lems
            E.info("Adding: "+p);
            Dimension paramDim = noDim;

            if (!p.isConstant()) {
                Exposure ex = new Exposure(p.getId(), paramDim);
                ct.exposures.add(ex);

                if (p.isSetValue()){

                    StateVariable sv = new StateVariable(p.getId(), paramDim, ex);
                    b.stateVariables.add(sv);

                    StateAssignment sa = new StateAssignment(p.getId(), p.getValue()+"");
                    os.stateAssignments.add(sa);
                }
            }
            else {
            	org.lemsml.jlems.core.type.Parameter lp = new org.lemsml.jlems.core.type.Parameter(p.getId(), paramDim);
                lp.name = p.getId();
                ct.parameters.add(lp);
                comp.setParameter(p.getId(), p.getValue()+"");
            }
        }
        
        ArrayList<FunctionDefinition> functions = new ArrayList<FunctionDefinition>();

        for (FunctionDefinition fd: model.getListOfFunctionDefinitions()){
            functions.add(fd);
        }

        E.info("functions: "+functions);


        for(Species s: model.getListOfSpecies()) {
            Dimension speciesDim = noDim;
            Exposure ex = new Exposure(s.getId(), speciesDim);
            ct.exposures.add(ex);
            StateVariable sv = new StateVariable(s.getId(), speciesDim, ex);
            b.stateVariables.add(sv);

            if (s.isSetValue()){
                StateAssignment sa = new StateAssignment(s.getId(), s.getValue()+"");
                os.stateAssignments.add(sa);
                E.info("Init sa: "+sa.getValueExpression());
            }

        }

        Constant timeScale = new Constant(tscaleName, lems.dimensions.getByName("per_time"), "1per_s");
        ct.constants.add(timeScale);


        for (Rule r: model.getListOfRules()){

            String formula = r.getFormula();
            
            try {
                formula = replaceFunctionDefinitions(formula, functions);
            }
            catch (Exception ex) {
                throw new ContentError("Problem substituting function definitions in SBML", ex);
            }

            E.info("Adding rule: "+r);

            if (r.isRate()){
                RateRule rr = (RateRule)r;


                TimeDerivative td = new TimeDerivative(rr.getVariable(), timeScale.getName() +" * ("+formula +")");
                b.timeDerivatives.add(td);
            }
            if (r.isAssignment()){
                Dimension speciesDim = noDim;
                AssignmentRule ar = (AssignmentRule)r;
                DerivedVariable dv = new DerivedVariable(ar.getVariable(), speciesDim,  formula, ar.getVariable());
                b.derivedVariables.add(dv);

                //TimeDerivative td = new TimeDerivative(rr.getVariable(), timeScale.getName() +" * ("+ rr.getFormula()+")");
                //b.timeDerivatives.add(td);
            }
        }

        for (Event e: model.getListOfEvents()){
            String testFormula = e.getTrigger().getFormula();
            E.info("Adding event: "+e+", test: "+testFormula);

            try {
                testFormula = replaceFunctionDefinitions(testFormula, functions);
                String test = replaceOperatorsAndTime(testFormula);
                E.info("Test tidied to: "+test);

                OnCondition oc = new OnCondition(test);
                b.onConditions.add(oc);
                for (EventAssignment ea: e.getListOfEventAssignments() ){
                    String formula = ea.getFormula();

                    formula = replaceFunctionDefinitions(formula, functions);


                    StateAssignment sa = new StateAssignment(ea.getVariable(), formula);
                    oc.stateAssignments.add(sa);
                }
            } catch (Exception ex) {
                throw new ContentError("Problem substituting function definitions in SBML", ex);
            }

        }
        if (os.stateAssignments.size()>0)
            b.onStarts.add(os);

        HashMap<String, StringBuilder> rates = new HashMap<String, StringBuilder>();

        for (Reaction reaction: model.getListOfReactions()){
            KineticLaw kl = reaction.getKineticLaw();
            HashMap<String, String> toReplace = new HashMap<String, String>();

            for(LocalParameter p: reaction.getKineticLaw().getListOfLocalParameters()) {
                //org.neuroml.lems.type.Lems
                String localName = p.getId()+"_"+reaction.getId();
                toReplace.put(p.getId(), localName);
                E.info("Adding: "+localName);
                Dimension paramDim = noDim;

                org.lemsml.jlems.core.type.Parameter lp = new org.lemsml.jlems.core.type.Parameter(localName, paramDim);
                lp.name = localName;
                ct.parameters.add(lp);
                comp.setParameter(localName, p.getValue()+"");
            }

            String formula = "("+kl.getFormula()+")";

            try{
                formula = replaceFunctionDefinitions(formula, functions);
            }
            catch (Exception ex){
                throw new ContentError("Problem substituting function definitions in SBML", ex);
            }
            
            //
            for(String old: toReplace.keySet()){
                String new_ = toReplace.get(old);
                String[] pres = new String[]{"\\(","\\+","-","\\*","/","\\^"};
                String[] posts = new String[]{"\\)","\\+","-","\\*","/","\\^"};

                for(String pre: pres){
                    for(String post: posts){
                        String o = pre+old+post;
                        String n = pre+" "+new_+" "+post;
                        formula = formula.replaceAll(o, n);
                        //E.info("Replacing "+o+" with "+n+": "+formula);
                    }
                }
                /*
                //formula = formula.replaceAll("/("+old, "/( "+new_);
                formula = formula.replaceAll("/+"+old, "+ "+new_);
                formula = formula.replaceAll("-"+old, "- "+new_);
                formula = formula.replaceAll("/*"+old, "* "+new_);
                formula = formula.replaceAll("//"+old, "/ "+new_);
                //formula = formula.replaceAll(new_+")", new_+" )");
                formula = formula.replaceAll(new_+"+", new_+" +");
                formula = formula.replaceAll(new_+"-", new_+" -");
                formula = formula.replaceAll(new_+"/*", new_+" *");
                formula = formula.replaceAll(new_+"//", new_+" /");*/
            }

            for (SpeciesReference product: reaction.getListOfProducts()){
                String s = product.getSpecies();
                Species sp = model.getListOfSpecies().get(s);

                if (!sp.getBoundaryCondition()) {
                    if (rates.get(s)==null) rates.put(s, new StringBuilder("0"));

                    StringBuilder sb = rates.get(s);

                    if (product.isSetStoichiometry() && product.getStoichiometry()!=1){
                        sb.append(" + ("+formula+" * "+product.getStoichiometry()+")");
                    } else {
                        sb.append(" + "+formula+"");
                    }
                    	
                }

            }
            for (SpeciesReference reactant: reaction.getListOfReactants()){

                String s = reactant.getSpecies();
                Species sp = model.getListOfSpecies().get(s);

                if (!sp.getBoundaryCondition()) {
                    if (rates.get(s)==null) rates.put(s, new StringBuilder("0"));

                    StringBuilder sb = rates.get(s);
                    sb.append(" - "+formula+"");
                    
                    if (reactant.isSetStoichiometry() && reactant.getStoichiometry()!=1){
                        sb.append(" * "+reactant.getStoichiometry()+"");
                    }
                }

            }

        }
        
        System.out.println(">>>> "+rates);

        for(String s: rates.keySet()){
        	Species sp = model.getSpecies(s);
            TimeDerivative td = new TimeDerivative(s, timeScale.getName() +" * ("+rates.get(s).toString()+") / "+sp.getCompartment());

            System.out.println(">>>> TimeDerivative "+td.getValueExpression());
            b.timeDerivatives.add(td);
        }

        if (simDuration>0 && simDt>0){

            Target dr = new Target();
            Component sim1 = new Component("sim1", lems.getComponentTypeByName("Simulation"));
            sim1.setParameter("length", simDuration+"s");
            sim1.setParameter("step", simDt+"s");
            sim1.setParameter("target", comp.getID());
            //sim1.setParameter("report",comp.getID()+"_report.txt");
            ////dr.timesFile = "examples/"+model.getId()+"_time.dat";

            Component disp1 = new Component("disp1", lems.getComponentTypeByName("Display"));
            disp1.setParameter("timeScale", "1s");
            disp1.setParameter("title", "Simulation of SBML model: "+ model.getId()+" from file: "+ sbmlFile);

            disp1.setParameter("xmin", "0");
            disp1.setParameter("xmax", simDuration+"");
            disp1.setParameter("ymin", "10");
            disp1.setParameter("ymax", "-10");

            sim1.addToChildren("displays", disp1);
            
            Component outF = new Component("outputFile1", lems.getComponentTypeByName("OutputFile"));
            String path = ".";
            if (dirForResults!=null)
            {
            	path = dirForResults.getAbsolutePath()+ System.getProperty("file.separator");
            }
            outF.setParameter("fileName", model.getId()+".dat");
            outF.setParameter("path", path);

            sim1.addToChildren("outputs", outF);

            int count = 1;
            /*
            for(Parameter p: model.getListOfParameters()) {

                if (!p.isConstant()){
                    Component lineCpt = new Component("lp_"+p.getId(), lems.getComponentTypeByName("Line"));
                    lineCpt.setParameter("scale", "1");
                    lineCpt.setParameter("quantity", p.getId());
                    Color c = ColorUtil.getSequentialColour(count);
                    String rgb = Integer.toHexString(c.getRGB());
                    rgb = rgb.substring(2, rgb.length());

                    lineCpt.setParameter("color", "#"+rgb);
                    lineCpt.setParameter("timeScale", "1s");

                    disp1.addToChildren("lines", lineCpt);

                    Component outputColumn = new Component("op_"+p.getId(), lems.getComponentTypeByName("OutputColumn"));
                    outputColumn.setParameter("quantity", p.getId());
                    outF.addToChildren("outputColumn", outputColumn);
                    
                    count++;
                }
            }*/

            for(Species s: model.getListOfSpecies()) {

                Component lineCpt = new Component("ls_"+s.getId(), lems.getComponentTypeByName("Line"));
                lineCpt.setParameter("scale", "1");
                lineCpt.setParameter("quantity", s.getId());
                Color c = ColorUtil.getSequentialColour(count);
                String rgb = Integer.toHexString(c.getRGB());
                rgb = rgb.substring(2, rgb.length());

                lineCpt.setParameter("color", "#"+rgb);
                lineCpt.setParameter("timeScale", "1s");

                disp1.addToChildren("lines", lineCpt);

                Component outputColumn = new Component("os_"+s.getId(), lems.getComponentTypeByName("OutputColumn"));
                outputColumn.setParameter("quantity", s.getId());
                outF.addToChildren("outputColumn", outputColumn);
                
                count++;
            }

            for(Compartment c: model.getListOfCompartments()) {

                if (!c.isConstant()){
                    Component lineCpt = new Component("lc_"+c.getId(), lems.getComponentTypeByName("Line"));
                    lineCpt.setParameter("scale", "1");
                    lineCpt.setParameter("quantity", c.getId());
                    Color col = ColorUtil.getSequentialColour(count);
                    String rgb = Integer.toHexString(col.getRGB());
                    rgb = rgb.substring(2, rgb.length());

                    lineCpt.setParameter("color", "#"+rgb);
                    lineCpt.setParameter("timeScale", "1s");

                    disp1.addToChildren("lines", lineCpt);

                    Component outputColumn = new Component("oc_"+c.getId(), lems.getComponentTypeByName("OutputColumn"));
                    outputColumn.setParameter("quantity", c.getId());
                    outF.addToChildren("outputColumn", outputColumn);
                    
                    count++;
                }
            }

            if (addModel)
                lems.addComponent(sim1);

            dr.component = sim1.getID();

            if (addModel)
                lems.targets.add(dr);
        }

        return lems;

    }

    private static String replaceFunctionDefinitions(String formula, ArrayList<FunctionDefinition> functions) throws SBMLException, ParseException{

        E.info("Substituting function defs in: "+formula);
        for (FunctionDefinition fd: functions){

            ASTNode exp = fd.getBody();
            E.info("-- Function def is: "+fd.getId()+"(...) = "+ exp);
            int count=0;
            String origFormula = new String(formula);

            while (formula.contains(fd.getId())){
                count++;
                if (count>20) throw new ParseException("Problem with formula: "+origFormula);
                
                int start = formula.indexOf(fd.getId());
                int end = formula.indexOf(")", start);
                String orig = formula.substring(start, end+1);
                E.info("orig: "+orig);

                String[] args = orig.substring(orig.indexOf("(")+1,orig.indexOf(")")).split(",");


                for(int i=0;i<args.length;i++){
                    String arg = args[i];
                    ASTNode a = fd.getArgument(i);
                    System.out.println("replacing "+a+" by "+arg);
                    ASTNode newA = JSBML.parseFormula(arg);
                    exp.replaceArgument(a.toString(), newA);
                }
                String newExp = "( "+ ASTNode.formulaToString(exp)+" )";

                orig = orig.replace("(", "\\(").replace(")", "\\)");
                E.info("formula was: "+ formula+", replacing "+ orig+" with "+newExp+", at "+formula.indexOf(orig)+", id: "+fd.getId());
                formula = formula.replaceAll(orig, newExp);
                if (formula.startsWith("( ") && formula.endsWith(" )"))
                    formula = formula.substring(2, formula.length()-2);
                E.info("formula is: "+ formula);
            }

        }
        return formula;
    }

    private static void runTest(File sbmlFile, float simDuration, float simDt) throws Exception
    {
        Lems lems = convertSBMLToLEMS(sbmlFile, simDuration, simDt, sbmlFile.getParentFile());

        //E.info("Generated: "+ lems.textSummary(true));
        //E.info("Generated: "+ lems.getComponentTypeByName("case00053").summary());
        lems.resolve();
        String lemsString  = XMLSerializer.serialize(lems);

        File testFile = new File(sbmlFile.getParent(), sbmlFile.getName().replaceAll(".xml", "")+"_SBML.xml");

        FileUtil.writeStringToFile(lemsString, testFile);

        E.info("Written to: "+ testFile.getCanonicalPath());

        E.info("Loading LEMS file from: "+ testFile.getAbsolutePath());

		Sim sim = Utils.loadLemsFileToSim(testFile);

        sim.readModel();
		sim.build();
		sim.run();
        E.info("Ran file at: "+ testFile.getCanonicalPath());
		
 

    }

    public static void main(String[] args) throws Exception
    {
        E.setDebug(true);

        
    	FileResultWriterFactory.initialize();
		DefaultLogger.initialize();
        /*
        SBMLImporter si = new SBMLImporter();

        SBMLDocument doc0 = new SBMLDocument(2, 4);

        E.info("Created doc: "+ doc0);

        new SBMLWriter().write(doc0, new BufferedOutputStream(System.out),
                                "TesterP", "Test2");*/

        //File sbmlFile = new File("exportImportUtils/SBML/Simple3Species.xml");
        String srcDir = "src/test/resources";
        File sbmlFile = new File(srcDir+"/Izhikevich.xml");

        sbmlFile = new File(srcDir+"/BIOMD0000000118.xml");
        
        sbmlFile = new File(srcDir+"/BIOMD0000000184.xml");
        
        sbmlFile = new File(srcDir+"/Simple3Species.xml");
        
        File sbmlFileDir = new File(srcDir+"/sbmlTestSuite/cases/semantic/");
            if (sbmlFileDir.exists()){
            sbmlFile = sbmlFileDir;
        }

        boolean sbmlTestSuite = sbmlFile.getAbsolutePath().indexOf("sbmlTestSuite")>=0;



        float len = 10;
        if (sbmlFile.getName().indexOf("Izh")>=0) len = 140;
        if (sbmlFile.getName().indexOf("0039")>=0) len = 50;
        if (sbmlFile.getName().indexOf("00118")>=0) len = 50;
        if (sbmlFile.getName().indexOf("00184")>=0) len = 1000;
        float dt = (float)(len/20000.0);

        HashMap<Integer, String> problematic = new HashMap<Integer, String>();
        //////problematic.put(21, "Incorrectly reads component size when units=volume ");
        ///problematic.put(22, "Incorrectly reads component size when units=volume ");
        //problematic.put(27, "Incorrectly reads component size when units=volume & init assignment");
        //problematic.put(86, "Complex args for function");

        StringBuilder errors = new StringBuilder();

        if (!sbmlTestSuite){
        	SwingDataViewerFactory.initialize();
            runTest(sbmlFile, len, dt);
        }
        else {
            int successful = 0;
            int completed = 0;
            int failed = 0;
            int notFound = 0;
            int skipped = 0;

            int numToStart = 1;
            int numToStop = 21;
            //numToStart = 18;
            //numToStop = 200;
            numToStop = 1123;
            //numToStop = 500;
            
            int numLemsPoints = 10000;
            float tolerance = 0.01f;

            if ((numToStop-numToStart)<=10)
        		SwingDataViewerFactory.initialize();
        		
            boolean exitOnError = false;
            //boolean exitOnMismatch = true;
            boolean exitOnMismatch = false;
            boolean skipFuncDefinitions = false;
            boolean skipUnitDefinitions = false;
            //String version = "l2v4";
            String version = "l3v1";
            
            for(int i=numToStart;i<=numToStop;i++){

                if (!problematic.keySet().contains(i)){

                    String testCase= "0"+i;

                    while(testCase.length()<5) testCase ="0"+testCase;


                    sbmlFile = new File(srcDir+"/sbmlTestSuite/cases/semantic/"+testCase+"/"+testCase+"-sbml-"+version+".xml");
                    if (!sbmlFile.exists()){
                        System.out.println("   ----  File not found: "+sbmlFile.getAbsolutePath()+"!!   ---- \n\n");
                        notFound++;
                    }
                    else
                    {
                        //File parent
                        File propsFile = new File(sbmlFile.getParentFile(), testCase+"-settings.txt");
                        String info = FileUtil.readStringFromFile(propsFile);
                        String duration = info.substring(info.indexOf("duration: ")+9, info.indexOf("steps:")-1).trim();
                        len = Float.parseFloat(duration);
                        dt = (float)(len/10000.0);


                        System.out.println("\n\n---------------------------------\n\nSBML test: "+testCase+" going to be run!");
                        try{
                            System.gc();
                            runTest(sbmlFile, len, dt);

                            System.out.println("SBML test: "+testCase+" completed!");
                            HashMap<String, float[]> targets = new HashMap<String, float[]>();
                            //HashMap<String, float[]> results = new HashMap<String, float[]>();

                            File targetFile = new File(sbmlFile.getParentFile(), testCase+"-results.csv");

                            System.out.println("Loading target data from "+targetFile.getCanonicalPath());
                            BufferedReader reader = new BufferedReader(new FileReader(targetFile));

                            String line = reader.readLine();
                            String[] dataNames = line.split(",");
                            int count = 0;
                            while ((line=reader.readLine()) != null) {
                                count++;
                            }
                            for(int d=0;d<dataNames.length;d++){

                                if (dataNames[d].indexOf("`")>=0)
                                    dataNames[d] = dataNames[d].substring(dataNames[d].indexOf("`") +1);

                                targets.put(dataNames[d], new float[count]);
                            }
                            E.info("Creating data arrays of size: "+count+" for "+targets.keySet());
                            reader.close();
                            
                            reader = new BufferedReader(new FileReader(targetFile));
                            line = reader.readLine();

                            int index = 0;
                            while ((line=reader.readLine()) != null) {
                                String[] data = line.split(",");
                                for(int d=0;d<data.length;d++){
                                    float[] dataArray = targets.get(dataNames[d]);
                                    //E.info("Adding index "+index+" to "+dataNames[d]+": "+ data[d]);
                                    dataArray[index] = Float.parseFloat(data[d]);
                                }
                                index++;
                            }
                            E.info(targets.toString());


                            SBMLReader sr = new SBMLReader();

                            SBMLDocument doc = sr.readSBML(sbmlFile);

                            Model model = doc.getModel();

                            E.info("Model: "+model.getListOfUnitDefinitions().size());
                            
                            if (model.getListOfFunctionDefinitions().size()>0 && skipFuncDefinitions) {
                            	String infoMessage = "Skipping: "+testCase+" due to function definitions";
	                            E.info(infoMessage);
	                            errors.append(testCase+": "+ infoMessage+"\n");
	                            skipped++;
                            } else if (model.getListOfUnitDefinitions().size()>0 && skipUnitDefinitions) {
                            	String infoMessage = "Skipping: "+testCase+" due to unit definitions";
	                            E.info(infoMessage);
	                            errors.append(testCase+": "+ infoMessage+"\n");
	                            skipped++;
                            } else {

	                            File resultFile = new File(sbmlFile.getParentFile(), "case"+testCase+".dat");
	                            ArrayList<String> cols = new ArrayList<String>();
	                            cols.add("time");
	
	                            for(Species s: model.getListOfSpecies()) {
	
	                                cols.add(s.getId());
	                            }
	                            for(Compartment c: model.getListOfCompartments()) {
	                            	if (!c.isConstant())
	                            		cols.add(c.getId());
	                            }
	
	                            E.info("Checking columns: "+cols+" in "+resultFile);
	
	                            reader.close();
	                            reader = new BufferedReader(new FileReader(resultFile));
	
	                            int lines = 0;
	                            while ((line=reader.readLine()) != null) {
	                            	lines++;
	                            }
	                            
	                            float[][] data = new float[cols.size()][lines];
	                            
	                            int lineNum = 0;
	                            reader.close();
	                            reader = new BufferedReader(new FileReader(resultFile));
	                            
	                            while ((line=reader.readLine()) != null) {
	                            	String[] words = line.split("\\s");
	                            	for (int c=0;c<cols.size();c++){
	    	                            float factor = 1;
	    	                            //if (cols.get(c).equals("time"))
	    	                            //    factor = 1000f;
	                            		data[c][lineNum] = Float.parseFloat(words[c])* factor;
	                            	}
	                            	lineNum++;
	                            }
	
	                            float[] timeTarg = targets.get("time");
	
	
	                            int ir = 0;
	                            boolean match = true;
	
	                            for (int it=0;it<timeTarg.length;it++){
	                                float tt = timeTarg[it];
	                                float tr = data[0][ir];
	
	                                while (tr < tt && ir<lines-1){
	                                    ir++;
	                                    tr = data[0][ir];
	                                }
	                                //E.info("------------ Testing target time point "+tt+" ("+it+" of "+timeTarg.length+") against sim data time point "+tr+" ("+ir+" of "+lines+")");
	
	
	                                for(int c=0;c<cols.size();c++){
	                                	String dataName = cols.get(c);
	                                    if(!dataName.equals("time")){
	                                        float[] dataTarg = targets.get(dataName);
	                                        float t = dataTarg[it];
	                                        float r = data[c][ir];
	                                        //E.info("--- Comparing val for "+dataName+" ("+c+") simulated: "+r+" against target "+t);
	                                        if (t!=0 && r!=0){
	                                            float diff = Math.abs((t-r)/t);
	                                            if (diff <=tolerance)
	                                            {
	                                                //E.info("---   Points match: Comparing val for "+dataName+" simulated: "+r+" against target "+t);
	                                            }
	                                            else
	                                            {
	                                                E.info("---   Points DON'T match at time "+tt+": Comparing val for "+dataName+" simulated: "+r+" against target "+t+ ", diff aimed at: "+tolerance+", real diff: "+ diff);
	                                                match = false;
	                                            }
	                                        }
	                                    }
	                                }
	                            }
	                            completed++;
	
	
	                            if (match)
	                            {
	                                successful++;
	                                E.info("\n    Success of test: "+testCase+"!!\n    So far: "+successful+" successful out of "+(successful+failed)+" ("+completed+" completed)\n");
	                            }
	                            else
	                            {
	                                E.info("Failure of test: "+testCase+"!\n    So far: "+successful+" successful out of "+(successful+failed)+" ("+completed+" completed)\n");
	                                if (exitOnMismatch) System.exit(-1);
	                                failed++;
	                            }
                            }

                            reader.close();
                            
                        } catch(Exception e){
                            System.out.println("\n\nSBML test: "+testCase+" failed!!\n");
                            e.printStackTrace();
                            errors.append(testCase+": "+ e.getMessage()+"\n");
                            if (exitOnError) System.exit(-1);

                            failed++;
                        }
                        
                    }
                }

            }

            E.info("\nAll finished!\n"
                    + "  Number successful:   "+successful+" out of "+(successful+failed)+" ("+(float)successful*100 / (successful+failed)+" %)\n"
                    + "  Number completed:     "+completed+"\n"
                    + "  Number failed:        "+failed+"\n"
                    + "  Number skipped:       "+skipped+"\n"
                    + "  Number not found:     "+notFound+"\n\n"
                    + "  Number LEMS points:   "+numLemsPoints+"\n"
                    + "  Tolerance:            "+tolerance);

            FileUtil.writeStringToFile(errors.toString(), new File("Errors.txt"));
        }


    }

    public static String replaceOperatorsAndTime(String formula)
    {
        String opsReplaced = formula.replaceAll("<=", ".leq.").replaceAll(">=", ".geq.").replaceAll("<", ".lt.").replaceAll(">", ".gt.").replaceAll("=", ".eq.");
        String timeReplaced = opsReplaced.replaceAll("time ", "(t * "+tscaleName+") ").replaceAll(" time", " (t * "+tscaleName+")");
        return timeReplaced;
    }
}


















