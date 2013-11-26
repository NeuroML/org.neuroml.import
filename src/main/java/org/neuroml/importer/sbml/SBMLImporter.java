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
import org.lemsml.jlems.core.expression.ParseTree;
import org.lemsml.jlems.core.expression.Parser;
import org.lemsml.jlems.core.logging.E;
import org.lemsml.jlems.core.run.ConnectionError;
import org.lemsml.jlems.core.run.RuntimeError;
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
import org.lemsml.jlems.core.type.Unit;
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
import org.lemsml.jlems.io.util.JUtil;
import org.lemsml.jlems.io.xmlio.XMLSerializer;
import org.lemsml.jlems.viz.datadisplay.SwingDataViewerFactory;
import org.lemsml.jlems.viz.plot.ColorUtil;
import org.neuroml.export.Utils;
import org.neuroml.importer.Main;
import org.neuroml.model.util.NeuroML2Validator;
import org.neuroml.model.util.NeuroMLElements;
import org.sbml.jsbml.ASTNode;
import org.sbml.jsbml.AssignmentRule;
import org.sbml.jsbml.Compartment;
import org.sbml.jsbml.Event;
import org.sbml.jsbml.EventAssignment;
import org.sbml.jsbml.FunctionDefinition;
import org.sbml.jsbml.InitialAssignment;
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
import org.sbml.jsbml.UnitDefinition;
import org.sbml.jsbml.text.parser.ParseException;

public class SBMLImporter  {

	static String DIM_SUFFIX = "_dimension";
	static String UNIT_SUFFIX = "_unit";
	static String INIT_PREFIX = "init_";

	//static boolean useUnits = true;
	static boolean useUnits = false;

    static final Dimension noDim = new Dimension(Dimension.NO_DIMENSION);
    
    static final Unit noUnit = new Unit(Unit.NO_UNIT, "", noDim);
    
    
    static Parser lemsExpressionParser = new Parser();
	
	
    public SBMLImporter() {
        E.info("Created new SBMLImporter...");
    }

    public static String tscaleName = "tscale";

    public static Lems convertSBMLToLEMS(File sbmlFile, float simDuration, float simDt) throws ContentError, XMLStreamException, ParseError, org.lemsml.jlems.core.sim.ParseException, BuildException, XMLException, IOException, SBMLException, ParseException, ConnectionError, RuntimeError, UnsupportedSBMLFeature {
    	return convertSBMLToLEMS(sbmlFile, simDuration, simDt, null);
    }

    private static Unit getUnit(String unit, HashMap<String, org.lemsml.jlems.core.type.Unit> units) {
    	if (unit==null || !units.containsKey(unit))
    		return noUnit;
    	return units.get(unit);
    }
    
    private static Dimension getDim(String dim, HashMap<String, Dimension> dims) {
    	if (dim==null || !dims.containsKey(dim))
    		return noDim;
    	return dims.get(dim);
    }
    

    @SuppressWarnings("deprecation")
	public static Lems convertSBMLToLEMS(File sbmlFile, float simDuration, float simDt, File dirForResults) throws ContentError, XMLStreamException, ParseError, org.lemsml.jlems.core.sim.ParseException, BuildException, XMLException, IOException, SBMLException, ParseException, ConnectionError, RuntimeError, UnsupportedSBMLFeature {

    	E.setDebug(false);
        SBMLReader sr = new SBMLReader();

        SBMLDocument doc = sr.readSBML(sbmlFile);

        Model model = doc.getModel();
        
        HashMap<String, Dimension> dims = new HashMap<String, Dimension>();
        HashMap<String, org.lemsml.jlems.core.type.Unit> units = new HashMap<String, org.lemsml.jlems.core.type.Unit>();


        E.info("Read in SBML from "+sbmlFile.getAbsolutePath());
        

    	NeuroML2Validator nmlv = new NeuroML2Validator();
    	
		String content = JUtil.getRelativeResource(nmlv.getClass(), "/NeuroML2CoreTypes/Simulation.xml");

        Sim sim = Utils.readLemsNeuroMLFile(content);
        
        //sim.build();
		Lems lems = sim.getLems();
		
		E.info("Loaded LEMS: "+lems.toString());
		

        boolean addModel = true;

        ComponentType ct = new ComponentType(model.getId());
        if (addModel)
            lems.addComponentType(ct);

        Component comp = new Component(model.getId()+"_0", ct);
        if (addModel)
            lems.addComponent(comp);


        Dynamics dyn = new Dynamics();
        ct.dynamicses.add(dyn);

        OnStart os = new OnStart();
        
        ArrayList<String> timeAliases = new ArrayList<String>();

        
        for(UnitDefinition ud: model.getListOfUnitDefinitions()){
            Dimension newDim = new Dimension(ud.getId()+DIM_SUFFIX);
            Unit newUnit = new Unit(ud.getId()+UNIT_SUFFIX, ud.getId()+UNIT_SUFFIX, newDim);
            
            for (org.sbml.jsbml.Unit u: ud.getListOfUnits()){
            	String kind = u.getKind().getName().toLowerCase();
            	if (kind.equals("ampere")) {
            		newDim.setI(1);
            	} else if (u.getKind().equals("farad")) {
            		newDim.setM(-1);
            		newDim.setL(-2);
            		newDim.setT(4);
            		newDim.setI(2);
            	} else if (kind.equals("litre")) {
            		newDim.setM(3);
            		newUnit.setPower(-3);
            	} else if (kind.equals("metre")) {
            		newDim.setM(1);
            		newUnit.setPower(1);
            	} else if (kind.equals("mole")) {
            		newDim.setN(1);
            	} else if (kind.equals("second")) {
            		newDim.setT(1);
            	} else {
            		//TODO: Add all unit kinds from section 4.4.2 in SBML specs: http://sbml.org/Documents/Specifications
            		System.err.print("Add more unit definitions! Missing: "+kind);
            		if (useUnits)
            			System.exit(1);
            		
            	}

            }
            lems.addDimension(newDim);
            lems.addUnit(newUnit);

            if (!useUnits) {
            	newDim = noDim;
            	newUnit = noUnit;
            } 

            dims.put(ud.getId(), newDim);
            units.put(ud.getId(), newUnit);
         
        }

        Constant timeScale = new Constant(tscaleName, lems.dimensions.getByName("per_time"), "1per_s");
        ct.constants.add(timeScale);

        for(Compartment c: model.getListOfCompartments()){
            Dimension compDim = getDim(c.getUnits(), dims);
            
            Unit compUnit = getUnit(c.getUnits(), units);
  
            String size = c.getSize()+"";
            if (!c.isSetSize())
            	size="1";
            
            E.info("Adding: "+c+" (size = "+size+" (set? "+c.isSetSize()+"), constant = "+c.isConstant()+", units = "+c.getUnits()+" ("+compDim+"))");
            
            size = size+" "+compUnit.getSymbol();
            
            boolean isInitAss = false;
            
            for(InitialAssignment ia: model.getListOfInitialAssignments()) {
        		if ( ia.getVariable().equals(c.getId()) ) {
        			isInitAss = true;
        		}
            }
            

            if (c.isConstant() && !isInitAss){
                Constant constComp = new Constant(c.getId(), compDim, size);
                ct.constants.add(constComp);
            } else {
                Exposure ex = new Exposure(c.getId(), compDim);
                ct.exposures.add(ex);

                //StateVariable sv = new StateVariable(c.getId(), compDim, ex);
                ///dyn.stateVariables.add(sv);
                
                boolean isStateVar = false;

                for(Rule r: model.getListOfRules()) {
                	if (r.isRate()) {
                		RateRule rr = (RateRule)r;
                		if ( rr.getVariable().equals(c.getId()) ) {
                            isStateVar = true;
                		}
                	}
                }
                for(Event e: model.getListOfEvents()) {
                	for (EventAssignment ea: e.getListOfEventAssignments()) {
                		if ( ea.getVariable().equals(c.getId()) ) {
                            isStateVar = true;
                		}
                	}
                }
                
                for(InitialAssignment ia: model.getListOfInitialAssignments()) {
            		if ( ia.getVariable().equals(c.getId()) ) {
                        isStateVar = true;
            		}
                }
                E.info("Adding: "+c+" isStateVar: "+isStateVar);
                
                if (isStateVar) {
                	StateVariable sv = new StateVariable(c.getId(), compDim, ex);
                    dyn.stateVariables.add(sv);
                    StateAssignment sa = new StateAssignment(c.getId(), c.getSize()+" "+compUnit.getSymbol());
                    os.stateAssignments.add(sa);
                }

            }


        }

        for(Parameter p: model.getListOfParameters()) {
            //org.neuroml.lems.type.Lems
            E.info("Adding: "+p);
            Dimension paramDim = getDim(p.getUnits(), dims);

            boolean isInitAss = false;
            
            for(InitialAssignment ia: model.getListOfInitialAssignments()) {
        		if ( ia.getVariable().equals(p.getId()) ) {
        			isInitAss = true;
        		}
            }

            if (!p.isConstant() || isInitAss) {
                Exposure ex = new Exposure(p.getId(), paramDim);
                ct.exposures.add(ex);

                boolean isStateVar = false;
                boolean hasRateRule = false;
                boolean hasAssRule = false;

                for(Rule r: model.getListOfRules()) {
                	if (r.isRate()) {
                		RateRule rr = (RateRule)r;
                		if ( rr.getVariable().equals(p.getId()) ) {
                            isStateVar = true;
                            hasRateRule = true;
                		}
                	}
                    else if (r.isAssignment()) {
                		AssignmentRule ar = (AssignmentRule)r;
                		if ( ar.getVariable().equals(p.getId()) ) {
                            hasAssRule = true;
                		}
                	}
                }
                for(Event e: model.getListOfEvents()) {
                	for (EventAssignment ea: e.getListOfEventAssignments()) {
                		if ( ea.getVariable().equals(p.getId()) ) {
                            isStateVar = true;
                		}
                	}
                }
                
                if (isInitAss)
                    isStateVar = true;
                
                E.info("  ---- Param: "+p.getId()+", isStateVar: "+isStateVar+", hasRateRule: "+hasRateRule+", isInitAss: "+isInitAss);
                
                if (isStateVar) {
                    StateVariable sv = new StateVariable(p.getId(), paramDim, ex);
                    dyn.stateVariables.add(sv);
                    E.info("  -- Adding state var: "+sv);
                    if (!hasRateRule) {
                        TimeDerivative td = new TimeDerivative(p.getId(), timeScale.getName() +" * 0");
                        dyn.timeDerivatives.add(td);
                    }
                }

                if (p.isSetValue()){
                	if (isStateVar) {
	                    E.info("Setting init param: "+p.getId() +" = "+p.getValue());
	                    StateAssignment sa = new StateAssignment(p.getId(), p.getValue()+"");
	                    os.stateAssignments.add(sa);
                	} else if (!hasAssRule) {
	                    E.info("! p.isSetValue()");
                		//E.info("Problem with "+sbmlFile+"\n");
                        //System.exit(-1);

                    	org.lemsml.jlems.core.type.Parameter lp = new org.lemsml.jlems.core.type.Parameter(p.getId(), paramDim);
                        lp.name = p.getId();
                        ct.parameters.add(lp);
                        comp.setParameter(p.getId(), p.getValue()+"");
                		
                	}
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
        
        
        HashMap<String, String> initVals = new HashMap<String, String>();

        for(Species s: model.getListOfSpecies()) {
            Dimension speciesDim = getDim(s.getSubstanceUnits(), dims);
            Unit speciesUnit = getUnit(s.getSubstanceUnits(), units);
            
            Exposure ex = new Exposure(s.getId(), speciesDim);
            ct.exposures.add(ex);
            StateVariable sv = new StateVariable(s.getId(), speciesDim, ex);
            dyn.stateVariables.add(sv);

            if (s.isSetValue()){

            	//TODO: check the need for this!!
            	String initConst = INIT_PREFIX+s.getId();
            	initVals.put(s.getId(), initConst);
                Constant constComp = new Constant(initConst, speciesDim, s.getValue()+" "+speciesUnit.getSymbol());
                ct.constants.add(constComp);
                
                StateAssignment sa = new StateAssignment(s.getId(), initConst);
                
                os.stateAssignments.add(sa);
                E.info("Init sa: "+sa.getValueExpression());
            }

        }
        

        for (InitialAssignment ia: model.getListOfInitialAssignments()) {
        	String var = ia.getVariable();
        	

            E.info("ct.constants: "+ct.constants);
            
        	String formula = handleFormula(ia.getMath(), functions, timeAliases);
            formula = replaceInFormula(formula, initVals);
            
            if (ct.constants.getByName(var)!=null)  {  // e.g. a compartment which is constant=true... case 00027
            	Constant c = ct.constants.getByName(var);
            	c.value = formula;
	            E.info("Resetting constant: "+var +" to "+c.value);
            } else {
	            StateAssignment sa = new StateAssignment(var, formula);
	            os.stateAssignments.add(sa);
	            E.info("InitialAssignment: "+var +" = "+sa.getValueExpression());
            }
        	
        }


        for (Rule r: model.getListOfRules()){
            
            String formula = handleFormula(r.getMath(), functions, timeAliases);
   

            E.info("Adding rule: "+r);

            if (r.isRate()){
                RateRule rr = (RateRule)r;
                boolean isSpecies = model.getSpecies(rr.getVariable())!=null;
                String scaling = "";
                if (isSpecies)
                    scaling = " * "+model.getSpecies(rr.getVariable()).getCompartment();
                TimeDerivative td = new TimeDerivative(rr.getVariable(), timeScale.getName()+ scaling +" * ("+formula +")");
                
                dyn.timeDerivatives.add(td);
                
            } else if (r.isAssignment()){
                
                AssignmentRule ar = (AssignmentRule)r;
                
                boolean isParameter = ct.parameters.getByName(ar.getVariable())!=null;
                if (isParameter) {
                    E.info("Resetting parameter: "+ar.getVariable()+" to: "+formula);
                    comp.setParameter(ar.getVariable(), formula);
                } else {
                
                    Dimension speciesDim = noDim;

                    DerivedVariable dv = new DerivedVariable(ar.getVariable(), speciesDim,  formula, ar.getVariable());
                    E.info("DerivedVariable: "+dv);

                    dyn.derivedVariables.add(dv);

                    formula = replaceInFormula(formula, initVals);

                    StateAssignment sa = null;
                    //E.info("Sas: "+os.getStateAssignments());
                    if (os.getStateAssignments()!=null) {
                        for (StateAssignment sa0: os.getStateAssignments()) {
                            //E.info("Sa: "+sa0);
                            if (sa0.getVariable().equals(ar.getVariable())) {
                                E.info("Replacing function for init: "+sa0.getVariable()+" = "+sa0.getValueExpression()
                                        +" with "+formula);
                                sa = sa0;
                            }
                        }
                    }

                    if (sa!=null) {
                        sa.setValue(formula);
                    } else {
                        ///sa = new StateAssignment(ar.getVariable());
                        ///os.stateAssignments.add(sa);
                        ////sa.setValue(formula);
                    }
                }

            } else if (r.isAlgebraic()){
            	throw new UnsupportedSBMLFeature("Algebraic rules are not yet supported in the SBML -> LEMS converter!");
            }
        }

        for (Event e: model.getListOfEvents()){

            String test = handleFormula(e.getTrigger().getMath(), functions, timeAliases, true);
            E.info("Adding event: "+e+", test: "+test);
            
            //test = replaceOperators(test);
            E.info("Test tidied to: "+test);

            OnCondition oc = new OnCondition(test);
            dyn.onConditions.add(oc);
            for (EventAssignment ea: e.getListOfEventAssignments() ){

                String formula = handleFormula(ea.getMath(), functions, timeAliases);

                StateAssignment sa = new StateAssignment(ea.getVariable(), formula);
                oc.stateAssignments.add(sa);
            }
  

        }
        if (os.stateAssignments.size()>0)
            dyn.onStarts.add(os);

        for (StateAssignment sa: os.getStateAssignments()) {
        	E.info("OnStarts: "+sa.variable+" = "+sa.value);
        }

        HashMap<String, StringBuilder> speciesTotalRates = new HashMap<String, StringBuilder>();

        for (Reaction reaction: model.getListOfReactions()){
            KineticLaw kl = reaction.getKineticLaw();
            HashMap<String, String> localParams = new HashMap<String, String>();
            HashMap<String, String> speciesScales = new HashMap<String, String>();

            for(LocalParameter p: reaction.getKineticLaw().getListOfLocalParameters()) {
                //org.neuroml.lems.type.Lems
                String localName = p.getId()+"_"+reaction.getId();
                localParams.put(p.getId(), localName);
                E.info("Adding: "+localName);
                Dimension paramDim = noDim;

                org.lemsml.jlems.core.type.Parameter lp = new org.lemsml.jlems.core.type.Parameter(localName, paramDim);
                lp.name = localName;
                ct.parameters.add(lp);
                comp.setParameter(localName, p.getValue()+"");
            }
            
            for (Species s: model.getListOfSpecies()) {

            	speciesScales.put(s.getId(), "("+s.getId()+"/"+s.getCompartment()+")");
            }

            String formula = handleFormula(kl.getMath(), functions, timeAliases);
            E.info("formula: "+formula+", derived units: "+kl.getDerivedUnits()+" ud "+kl.containsUndeclaredUnits());
      
            
            //
            formula = replaceInFormula(formula, localParams);
            formula = replaceInFormula(formula, speciesScales);
        
            E.info("formula now: "+formula);
            String rateDvName = "rate__"+reaction.getId();
            DerivedVariable dv = new DerivedVariable(rateDvName, noDim, formula);
            E.info("DerivedVariable: "+dv);

            dyn.derivedVariables.add(dv);
            
            for (SpeciesReference product: reaction.getListOfProducts()){
                String s = product.getSpecies();
                if (product.getId().length()>0)
                {
                	throw new UnsupportedSBMLFeature("stoichiometryMath/using id in speciesReference not supported!");
                }
                Species sp = model.getListOfSpecies().get(s);

                if (!sp.getBoundaryCondition()) {
                    if (speciesTotalRates.get(s)==null) speciesTotalRates.put(s, new StringBuilder(""));

                    StringBuilder sb = speciesTotalRates.get(s);
                    String pre = " + ";
                    if (sb.length()==0) pre = "";

                    if (product.isSetStoichiometry() && product.getStoichiometry()!=1){
                        sb.append(pre+"("+rateDvName+" * "+product.getStoichiometry()+")");
                    } else {
                        sb.append(pre+rateDvName+"");
                    }
                    	
                }

            }
            for (SpeciesReference reactant: reaction.getListOfReactants()){

                String s = reactant.getSpecies();
                if (reactant.getId().length()>0)
                {
                	throw new UnsupportedSBMLFeature("stoichiometryMath/using id in speciesReference not supported!");
                }
                Species sp = model.getListOfSpecies().get(s);

                if (!sp.getBoundaryCondition()) {
                    if (speciesTotalRates.get(s)==null) speciesTotalRates.put(s, new StringBuilder(""));

                    StringBuilder sb = speciesTotalRates.get(s);
                    if (sb.length()==0) {
                    	sb.append("-1*"+rateDvName+"");
                    } else {
                    	sb.append(" - "+rateDvName+"");
                    }
                    
                    if (reactant.isSetStoichiometry() && reactant.getStoichiometry()!=1){

                        sb.append(" * "+reactant.getStoichiometry()+"");
                    }
                }

            }

        }
        
        E.info(">>>> speciesTotalRates: "+speciesTotalRates);

        for(String s: speciesTotalRates.keySet()){
        	//Species sp = model.getSpecies(s);
            TimeDerivative td = new TimeDerivative(s, timeScale.getName() +" * ("+speciesTotalRates.get(s).toString()+") ");

            E.info(">>>> TimeDerivative "+td.getValueExpression());
            dyn.timeDerivatives.add(td);
        }
        
        for (String timeVarName: timeAliases) {
            StateVariable sv = new StateVariable(timeVarName, noDim);
            E.info(">>>> StateVariable "+sv);
            dyn.stateVariables.add(sv);

            TimeDerivative td = new TimeDerivative(sv.getName(), "1");
            dyn.timeDerivatives.add(td);
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
            disp1.setParameter("ymin", "0");
            disp1.setParameter("ymax", "1");

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
            

            for(Species s: model.getListOfSpecies()) {

                Component lineCpt = new Component(s.getId()+"__S", lems.getComponentTypeByName("Line"));
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
                    Component lineCpt = new Component(c.getId()+"__C", lems.getComponentTypeByName("Line"));
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

            for(Parameter p: model.getListOfParameters()) {

                if (!p.isConstant()){
                    Component lineCpt = new Component(p.getId()+"__P", lems.getComponentTypeByName("Line"));
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
            }

            if (addModel)
                lems.addComponent(sim1);

            dr.component = sim1.getID();

            if (addModel)
                lems.targets.add(dr);
        }
		E.info("Loaded LEMS: "+lems.toString());
        

        return lems;

    }

    private static String replaceInFormula(String formula, HashMap<String, String> oldVsNew) {

        //E.info("Replacing in: "+formula+" with: "+oldVsNew);
    	for(String old: oldVsNew.keySet()) {
    		String new_ = oldVsNew.get(old);
    		formula = replaceInFormula(formula, old, new_);
    	}
        E.info("Now: "+formula);
    	return formula;
    }
    
    private static String replaceInFormula(String formula, String oldVal, String newVal) {
    	formula = " "+formula+" ";
    	String[] pres = new String[]{"\\(","\\+","-","\\*","/","\\^", " "};
        String[] posts = new String[]{"\\)","\\+","-","\\*","/","\\^", " "};

        for(String pre: pres){
            for(String post: posts){
                String o = pre+oldVal+post;
                String n = pre+" "+newVal+" "+post;
	                //E.info("Replacing "+o+" with "+n+": "+formula);
                //if (formula.indexOf(o)>=0) {
	                formula = formula.replaceAll(o, n);
                //}
            }
        }
        return formula.trim();
    }
    
    
    
    private static String astNodeToString(ASTNode ast, boolean condition) throws XMLStreamException, ParseError, ContentError
    {
        //String mml = JSBML.writeMathMLToString(ast);
        //E.info("MathML: "+mml);
        /*String mml2 = mml.replaceAll("cn type=\"integer\"", "cn type=\"real\"");
        E.info("MathML2: "+mml2);
        ASTNode ast2 = JSBML.readMathMLFromString(mml2);*/
        
        String expression = JSBML.formulaToString(ast);
        //E.info("expression: "+expression);
          
        return expression;
        
    }

    private static String handleFormula(ASTNode ast, ArrayList<FunctionDefinition> functions, ArrayList<String> timeAliases) throws SBMLException, ParseException, UnsupportedSBMLFeature, ParseError, ContentError, XMLStreamException {
        return handleFormula(ast, functions, timeAliases, false);
    }

    private static String handleFormula(ASTNode ast, ArrayList<FunctionDefinition> functions, ArrayList<String> timeAliases, boolean condition) throws SBMLException, ParseException, UnsupportedSBMLFeature, ParseError, ContentError, XMLStreamException {
        
        String formula = astNodeToString(ast,condition);
        String mml = JSBML.writeMathMLToString(ast);
        //E.info("MathML: "+mml);
        if (mml.indexOf("csymbol")>=0 && mml.indexOf("http://www.sbml.org/sbml/symbols/time")>=0) {
            int start = mml.indexOf("http://www.sbml.org/sbml/symbols/time");
            String varName = mml.substring(mml.indexOf(">", start)+1, mml.indexOf("<", start)).trim();
            
            if (!timeAliases.contains(varName)) 
                timeAliases.add(varName);
        }
        
    	checkFormula(formula);
    	String formula0 = replaceTime(formula);
    	String formula1 = replaceFactorial(formula0);
    	String formula2 = replaceFunctionDefinitions(formula1, functions);
        
        if (condition)
        {
            return replaceOperators(formula2);
        }
        else
        {
            ParseTree pt = lemsExpressionParser.parseExpression(formula2);
    
            return pt.toExpression();
        }
    }
    
    private static void checkFormula(String formula) throws UnsupportedSBMLFeature {
    	
    	E.info("Checking: "+formula);
        if (formula.indexOf("piecewise")>=0)
        {
        	E.info("Fail!");
        	throw new UnsupportedSBMLFeature("Piecewise function expressions are not yet supported in the SBML -> LEMS converter!");
        }
    }

    private static String replaceFunctionDefinitions(String formula, ArrayList<FunctionDefinition> functions) throws SBMLException, ParseException, UnsupportedSBMLFeature{

        E.info("---- Substituting function defs in: "+formula);
        
        for (FunctionDefinition fd: functions){

            int count=0;
            String origFormula = new String(formula);
            

            while (formula.contains(fd.getId())){
                E.info("Working on formula: "+ formula);

                ASTNode exp0 = fd.getBody();
                ASTNode exp = (ASTNode)exp0.clone();
                E.info("-- Function def is: "+fd.getId()+"(...) = "+ exp);
                count++;
                if (count>20) throw new ParseException("Problem with formula: "+origFormula);
                
                int start = formula.lastIndexOf(fd.getId());
                //int end = formula.indexOf(")", start);
                int depth = 0;
                int end = -1;
                int index = start+fd.getId().length()+1;
                while (end<0 && index<formula.length())
                {
                	char c = formula.charAt(index);
                	//E.info(c+"");
                	if (c=='(') depth++;
                	if (c==')')
                	{
                		if (depth==0)
                			end=index;
                		else
                			depth--;
                	}
                	index++;
                }
                String orig = formula.substring(start, end+1);
                E.info("orig: "+orig);

                String[] args = orig.substring(orig.indexOf("(")+1,orig.lastIndexOf(")")).split(",");


                for(int i=0;i<args.length;i++){
                    String arg = args[i].trim();
                    ASTNode a = fd.getArgument(i);
                    E.info("replacing <"+a+"> by <"+arg+"> in "+exp);
                    ASTNode newA = JSBML.parseFormula(arg);
                    //E.info(exp.toFormula());
                    exp.replaceArgument(a.toString(), newA);
                    //E.info(exp.toFormula());
                }
                String newExp = "( "+ ASTNode.formulaToString(exp)+" )";

                orig = orig.replace("(", "\\(").replace(")", "\\)").replace("*", "\\*");
                E.info("formula was: "+ formula+", replacing "+ orig+" with "+newExp+", at "+formula.indexOf(orig)+", id: "+fd.getId());
                formula = formula.replaceAll(orig, newExp);
                if (formula.startsWith("( ") && formula.endsWith(" )"))
                    formula = formula.substring(2, formula.length()-2);
                E.info("formula is: "+ formula);
            }

        }
        return formula;
    }
    
    public static File convertSBMLToLEMSFile(File sbmlFile, float simDuration, float simDt, boolean comment) throws SBMLException, ContentError, XMLStreamException, ParseError, org.lemsml.jlems.core.sim.ParseException, BuildException, XMLException, IOException, ParseException, ConnectionError, RuntimeError, UnsupportedSBMLFeature
    {
        Lems lems = convertSBMLToLEMS(sbmlFile, simDuration, simDt, sbmlFile.getParentFile());
        //E.info("Generated: "+ lems.textSummary(true));
        //E.info("Generated ct: "+ lems.getComponentTypeByName("case00471"));
        lems.resolve();
        String lemsString  = XMLSerializer.serialize(lems);
        
        if (comment)
        {
        	String commentString = "    This LEMS file has been generated by org.neuroml.importer.sbml.SBMLImporter (see https://github.com/NeuroML/org.neuroml.import)\n" +
        			"         org.neuroml.import  v"+Main.ORG_NEUROML_IMPORT_VERSION+"\n" +
                    "         org.neuroml.model   v"+NeuroMLElements.ORG_NEUROML_MODEL_VERSION+"\n" +
                    "         jLEMS               v"+org.lemsml.jlems.io.Main.VERSION;
        	
        	lemsString = lemsString.replace("<Lems>", "<Lems>\n<!--\n"+commentString+"\n-->\n");
        }
        
        File lemsFile = new File(sbmlFile.getParent(), sbmlFile.getName().replaceAll(".xml", "")+"_LEMS.xml");

        FileUtil.writeStringToFile(lemsString, lemsFile);

        E.info("Written to: "+ lemsFile.getCanonicalPath());
        
        return lemsFile;
    }

    private static void runTest(File sbmlFile, float simDuration, float simDt) throws SBMLException, ContentError, XMLStreamException, ParseError, org.lemsml.jlems.core.sim.ParseException, BuildException, XMLException, IOException, ParseException, ConnectionError, RuntimeError, UnsupportedSBMLFeature
    {
    	File testFile = convertSBMLToLEMSFile(sbmlFile, simDuration, simDt, true);

        E.info("Loading LEMS file from: "+ testFile.getAbsolutePath());

		Sim sim = Utils.readLemsNeuroMLFile(testFile);
		sim.build();
		sim.run();
        E.info("Ran file at: "+ testFile.getCanonicalPath());
		
 

    }

    public static void main(String[] args) throws Exception
    {
        E.setDebug(true);

        String[] exprs = {"(4)!","((4)!+4)","( (5 + (4)!) +4)", "sin(g)", "((ceil(p1*S1))!*p2^(-1))"};
        for (String expr: exprs){
        	//E.info("------------------------");
        	String rep = replaceFactorial(expr);
        	//E.info("Expr: ["+expr+"] -> ["+rep+"]");
        }
        //System.exit(0);
        
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
        
        sbmlFile = new File(srcDir+"/BIOMD0000000184.xml");
        
        sbmlFile = new File(srcDir+"/Simple3Species.xml");

        sbmlFile = new File(srcDir+"/BIOMD0000000118.xml");
        
        sbmlFile = new File(srcDir+"/BIOMD0000000138_SBML-L3V1.xml");
        
        
        boolean overrideLocalSTS = true;
        overrideLocalSTS = false;
        
        File sbmlFileDir = new File("sbmlTestSuite/cases/semantic/");
            if (sbmlFileDir.exists() && !overrideLocalSTS){
            sbmlFile = sbmlFileDir;
        }

        
        boolean sbmlTestSuite = sbmlFile.getAbsolutePath().indexOf("sbmlTestSuite")>=0;


        float len = 10;
        if (sbmlFile.getName().indexOf("Izh")>=0) len = 140;
        if (sbmlFile.getName().indexOf("0039")>=0) len = 50;
        if (sbmlFile.getName().indexOf("00118")>=0) len = 500;
        if (sbmlFile.getName().indexOf("00184")>=0) len = 1000;
        if (sbmlFile.getName().indexOf("00138")>=0) len = 3000;
        float dt = (float)(len/20000.0);

        HashMap<Integer, String> problematic = new HashMap<Integer, String>();
        //////problematic.put(21, "Incorrectly reads component size when units=volume ");
        ///problematic.put(22, "Incorrectly reads component size when units=volume ");
        //problematic.put(27, "Incorrectly reads component size when units=volume & init assignment");
        //problematic.put(86, "Complex args for function");

        StringBuilder errors = new StringBuilder();

        if (overrideLocalSTS || !sbmlTestSuite){
        	SwingDataViewerFactory.initialize();
            runTest(sbmlFile, len, dt);
        }
        else {
            int successful = 0;
            int completed = 0;
            int failedMismatch = 0;
            int failedError = 0;
            int unsupported = 0;
            
            int notFound = 0;
            int skipped = 0;

            int numToStart = 1; 
            int numToStop = 200;
            numToStop = 1180;
            //numToStop = 500;
            //numToStart = 1065;
            //numToStop = 700;
            
            int numLemsPoints = 30000;
            float tolerance = 0.01f;

            if ((numToStop-numToStart)<=10)
        		SwingDataViewerFactory.initialize();

            boolean exitOnError = false;
            //exitOnError = true;
            boolean exitOnMismatch = true;
            exitOnMismatch = false;

            boolean skipFuncDefinitions = false;
            
            boolean skipUnitDefinitions = false;

            boolean skipDelays = true;
            
            //String version = "l2v4";
            String version = "l3v1";
            
            for(int i=numToStart;i<=numToStop;i++){

                if (!problematic.keySet().contains(i)){

                    String testCase= "0"+i;

                    while(testCase.length()<5) testCase ="0"+testCase;


                    sbmlFile = new File("sbmlTestSuite/cases/semantic/"+testCase+"/"+testCase+"-sbml-"+version+".xml");
                    if (!sbmlFile.exists()){
                        E.info("   ----  File not found: "+sbmlFile.getAbsolutePath()+"!!   ---- \n\n");
                        notFound++;
                    }
                    else
                    {
                        //File parent
                        File propsFile = new File(sbmlFile.getParentFile(), testCase+"-settings.txt");
                        String info = FileUtil.readStringFromFile(propsFile);
                        String duration = info.substring(info.indexOf("duration: ")+9, info.indexOf("steps:")-1).trim();
                        len = Float.parseFloat(duration);
                        dt = (float)(len/numLemsPoints);


                        E.info("\n\n---------------------------------\n\nSBML test: "+testCase+" going to be run!");
                        try{
                            System.gc();
                            runTest(sbmlFile, len, dt);

                            E.info("SBML test: "+testCase+" completed!");
                            HashMap<String, float[]> targets = new HashMap<String, float[]>();
                            //HashMap<String, float[]> results = new HashMap<String, float[]>();

                            File targetFile = new File(sbmlFile.getParentFile(), testCase+"-results.csv");

                            E.info("Loading target data from "+targetFile.getCanonicalPath());
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

                            boolean containsDelays = false;
                            for(Event e: model.getListOfEvents()) {
                            	if (e.getDelay()!=null)
                            		containsDelays = true;
                            }
                            
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
                            } else if (containsDelays && skipDelays) {
                            	String infoMessage = "Skipping: "+testCase+" due to Delay in Event";
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
	                            /*
	                            for(Compartment c: model.getListOfCompartments()) {
	                            	if (!c.isConstant())
	                            		cols.add(c.getId());
	                            }*/
	
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
                                
                                float minFactor = 1e-12f;
                                float[] maxAsbValues = new float[cols.size()];
	
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
                                            maxAsbValues[c] = Math.max(maxAsbValues[c], Math.abs(t));
	                                        //E.info("--- Comparing val for "+dataName+" ("+c+") simulated: "+r+" against target "+t+", max so far: "+maxAsbValues[c]);
	                                        if (t!=0 && r!=0){
	                                            float diffAbs = Math.abs(t-r);
	                                            float diffFract = diffAbs/t;
	                                            if (diffFract <=tolerance)
	                                            {
	                                                //E.info("---   Points match: Comparing val for "+dataName+" simulated: "+r+" against target "+t);
	                                            }
                                                else if ( (maxAsbValues[c]*diffAbs) < minFactor)
	                                            {
	                                                E.info("---   Points DON'T match at time "+tt+" BUT ALLOWING because much smaller than max absolute value ("+maxAsbValues[c]+"): Comparing val for "+dataName+" simulated: "+r+" against target "+t+ ", diff aimed at: "+tolerance+", real diff: "+ diffFract);
	                                                //match = false;
	                                            }
	                                            else
	                                            {
	                                                E.info("---   Points DON'T match at time "+tt+": Comparing val for "+dataName+" simulated: "+r+" against target "+t+ ", diff aimed at: "+tolerance+", real diff: "+ diffFract);
	                                                match = false;
	                                            }
	                                        }
	                                    }
	                                }
	                            }
	                            completed++;
	

	                            int tried = successful+failedMismatch+failedError;
	                            if (match)
	                            {
	                                successful++;
	                                E.info("\n    Success of test: "+testCase+"!!\n    So far: "+successful+" successful out of "+tried+" ("+completed+" completed)\n");
	                            }
	                            else
	                            {
	                                E.info("Failure of test: "+testCase+"!\n    So far: "+successful+" successful out of "+tried+" ("+completed+" completed)\n");
	                                if (exitOnMismatch) System.exit(-1);
	                                failedMismatch++;
	                            }
                            }

                            reader.close();
                            
                        } catch(UnsupportedSBMLFeature e){
                            E.info("\n\nSBML test: "+testCase+" can't be run due to unsupported feature!!\n");
                            e.printStackTrace();
                            errors.append(testCase+": "+ e.getMessage()+"\n");

                            unsupported++;
                        } catch(Exception e){
                            E.info("\n\nSBML test: "+testCase+" failed!!\n");
                            e.printStackTrace();
                            errors.append(testCase+": "+ e.getMessage()+"\n");
                            if (exitOnError) 
                            	System.exit(-1);

                            failedError++;
                        }
                        
                    }
                }

            }
            
            int tried = successful+failedMismatch+failedError;
            
            E.info("\nAll finished!\n"
                    + "  Number successful:          "+successful+" out of "+tried+" ("+(float)successful*100 / (tried)+" %)\n"
                    + "  Number completed:           "+completed+"\n"
                    + "  Number failed (mismatch):   "+failedMismatch+"\n"
                    + "  Number failed (exception):  "+failedError+"\n"
                    + "  Number unsupported:         "+unsupported+"\n"
                    + "  Number skipped:             "+skipped+"\n"
                    + "  Number not found:           "+notFound+"\n\n"
                    + "  Range:                      "+numToStart+" -> "+numToStop+"\n\n"
                    + "  Number LEMS points:         "+numLemsPoints+"\n"
                    + "  Tolerance:                  "+tolerance);

            FileUtil.writeStringToFile(errors.toString(), new File("Errors.txt"));
        }


    }

    public static String replaceTime(String formula)
    {
        String timeReplaced = replaceInFormula(formula, "time", "(t * "+tscaleName+")");
        return timeReplaced;
    }
    public static String replaceOperators(String formula)
    {
        String opsReplaced = formula.replaceAll("<=", ".leq.").replaceAll(">=", ".geq.").replaceAll("<", ".lt.").replaceAll(">", ".gt.").replaceAll("=", ".eq.");
    
        return opsReplaced;
    }
    

    public static String replaceFactorial(String formula){
        while (formula.contains(")!")){
        	int index = formula.indexOf(")!");
        	int indexStart = -1;
        	int depth = 0;
        	for (int i =index+1; i<formula.length();i++) {
        		char c = formula.charAt(i);
        		if (c == ')') depth++;
        		if (c == '(') depth--;
        	}
        	int depth2 = 0;
        	for (int i =0; i<index;i++) {
        		char c = formula.charAt(i);
        		if (c == '('){
        			if (depth2==depth)
        			{
        				indexStart = i;
        			}
        			depth2++;
        		}
        		if (c == ')') depth2--;
        		
        	}
        	//E.info("Index: "+index+", depth: "+depth+", depth2: "+depth2+", indexStart: "+indexStart);
        	formula = formula.replaceFirst("\\)!", ")");
        	formula = formula.substring(0,indexStart)+"factorial("+formula.substring(indexStart+1);
        }
        return formula;
    }
}


















