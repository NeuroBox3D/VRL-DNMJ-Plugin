/**
 * VRL-DNMJ-Plugin
 * DNMJ.java
 * 
 * @date 2014-06-18
 * @author mstepnie
**/

package vrl.dnmj.plugin;

//import eu.mihosoft.vrl.annotation.ComponentInfo;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

//import edu.gcsc.vrl.ug.api.I_ApproximationSpace;
//import edu.gcsc.vrl.ug.api.I_UserNumber;
//import edu.gcsc.vrl.ug.api.I_UserVector;
import edu.gcsc.vrl.ug.api.*;

//import edu.gcsc.vrl.userdata.*;
import edu.gcsc.vrl.userdata.UserDataTuple;
import edu.gcsc.vrl.userdata.FunctionDefinition;
import edu.gcsc.vrl.userdata.UserDependentSubsetModel;

import eu.mihosoft.vrl.annotation.ComponentInfo;
import eu.mihosoft.vrl.annotation.MethodInfo;
import eu.mihosoft.vrl.annotation.OutputInfo;
import eu.mihosoft.vrl.annotation.ParamInfo;
import eu.mihosoft.vrl.annotation.ParamGroupInfo;


@ComponentInfo(name="VRL-DNMJ-Plugin", category="VRL-DNMJ-Plugin", description="Problem setup component")

public class DNMJ implements java.io.Serializable{
    private static final long serialVersionUID = 1L;

    @MethodInfo(valueStyle="multi-out", interactive = false)
    @OutputInfo
    (
        style="multi-out",
        elemStyles = {"silent", "silent", "silent", "silent", "silent", "silent", "silent"},
        elemNames = {"Domain Disc", "Approximation Space", "Initial Solution", "Stimulation start", "Stimulation end", "Simulation dt", "VGCC dt"},
        //elemNames = {"Domain Disc", "Approximation Space", "Initial Solution"},
        elemTypes = {I_DomainDiscretization.class, I_ApproximationSpace.class, UserDataTuple[].class, double.class, double.class, double.class, double.class}
        //elemTypes = {I_DomainDiscretization.class, I_ApproximationSpace.class, UserDataTuple[].class}
    )
    
    public Object[] invoke
    (
    //  DOMAIN
        @ParamGroupInfo(group="DOMAIN")
        @ParamInfo(name="Grid", style="ugx-load-dialog", options="ugx_tag=\"gridFile\"")
        java.io.File gridFile,
            
        @ParamGroupInfo(group="DOMAIN|false;Refiner")
        @ParamInfo(name="Refinements", style="default", options="value=0")
        Integer num_refs,
        
        @ParamGroupInfo(group="DOMAIN|false;Refiner")
        @ParamInfo(name="Pre-refinements", style="default", options="value=0")
        Integer num_prerefs,
        
    //  SIMULATION PARAMS
        @ParamGroupInfo(group="SIMULATION PARAMS|false;Stimulation protocol")
        @ParamInfo(name="Stimulation start [s]", style="default", options="value=0.0")
        double stimTimeBegin,
        
        @ParamGroupInfo(group="SIMULATION PARAMS|false;Stimulation protocol")
        @ParamInfo(name="Stimulation end [s]", style="default", options="value=0.01")
        double stimTimeEnd,
        
        @ParamGroupInfo(group="SIMULATION PARAMS|false;Stimulation protocol")
        @ParamInfo(name="Stimulation frequency [Hz]", style="default", options="value=40")
        double stimFrequency,
        
        @ParamGroupInfo(group="SIMULATION PARAMS|false;Stimulation protocol")
        @ParamInfo(name="Simulation time step size dt [s]", style="default", options="value=0.001")
        double dt,
        
        @ParamGroupInfo(group="SIMULATION PARAMS|false;Stimulation protocol")
        @ParamInfo(name="VGCC time step size VGCC_dt [s]", style="default", options="value=0.0001")
        double VGCC_dt,                
        
    //  PROBLEM DEFINITION
        @ParamGroupInfo(group="PROBLEM DEFINITION|true;Functions|true")
        @ParamInfo(name="Involved ion species", style="array", options="ugx_tag=\"gridFile\"; fct_tag=\"fctDef\"; minArraySize=1;")
        FunctionDefinition[] problemDefinition,
        
        @ParamGroupInfo(group="PROBLEM DEFINITION|true; Diffusion|false")
        @ParamInfo(name="Involved ion species", style="array", options="ugx_tag=\"gridFile\"; fct_tag=\"fctDef\"; minArraySize=1; type=\"S1|mnn:function & subset, diffusion [um^2/s], reaction rate, reaction term\"")
        UserDataTuple[] diffusionData,
        
        @ParamGroupInfo(group="PROBLEM DEFINITION|true; Buffering|false")
        @ParamInfo(name="Involved ion species [uM]", style="array", options="ugx_tag=\"gridFile\"; fct_tag=\"fctDef\"; minArraySize=0; type=\"S2|nnn:buffering substance, buffered substance, total buffer [mM], association rate [1/(mM*s)], dissociation rate [1/s]\"")
        UserDataTuple[] bufferingData,
        
        @ParamGroupInfo(group="PROBLEM DEFINITION|true; Plasma membrane|false; PMCA")
        @ParamInfo(name="", style="default", options="ugx_tag=\"gridFile\"; fct_tag=\"fctDef\"; type=\"S2|n:cytosolic calcium, extracellular calcium, density [1/um^2]\"")
        UserDataTuple pmcaData,
        
        @ParamGroupInfo(group="PROBLEM DEFINITION|true; Plasma membrane|false; NCX")
        @ParamInfo(name="", style="default", options="ugx_tag=\"gridFile\"; fct_tag=\"fctDef\"; type=\"S2|n:cytosolic calcium, extracellular calcium, density [1/um^2]\"")
        UserDataTuple ncxData,      
        
        @ParamGroupInfo(group="PROBLEM DEFINITION|true; Plasma membrane|false; VGCC")
        @ParamInfo(name="", style="default", options="ugx_tag=\"gridFile\"; fct_tag=\"fctDef\"; type=\"S2|n:extracellular calcium, cytosolic calcium, density [1/um^2]\"")
        UserDataTuple vgccData,
        
        @ParamGroupInfo(group="PROBLEM DEFINITION|true; Plasma membrane|false; VGCC")
        @ParamInfo(name="Channel type", style="selection", options="value=[\"L\",\"N\",\"T\"]")
        String vgccChannelType,        
        
        @ParamGroupInfo(group="PROBLEM DEFINITION|true; Start Value|false")
        @ParamInfo(name="Start Value [mM]", style="array", options="ugx_tag=\"gridFile\"; fct_tag=\"fctDef\"; minArraySize=1; type=\"S1|n:function & subset, start value\"")
        UserDataTuple[] startValue
    )            
    {
    //  get selected geometry and its dim
        String gridFileName = gridFile.getAbsolutePath();
        UGXFileInfo ugxFI = new UGXFileInfo();

        
    //  parse ugx file for world dimension
        if (ugxFI.parse_file(gridFileName) == false)
            throw new RuntimeException("Unable to parse ugx-File: " + gridFileName);
        
        if (ugxFI.const__num_grids() != 1)
            throw new RuntimeException("UGX file must contain exactly one grid.");
        
        int dim = ugxFI.const__grid_world_dimension(0);

        
    //  Init UG for dimension and algebra
        F_InitUG.invoke(dim, new AlgebraType("CPU", 1));
        
        
    //  create, load, refine and distribute domain
        /*
        System.out.println("Create, refine and distribute domain");
        I_Domain dom = new Domain();
        
        System.out.print("Loading domain " + gridFileName + " ... ");
        F_LoadDomain.invoke(dom, gridFileName);
	System.out.print("done.\n");
        
        I_IRefiner refiner = null;
	if (num_refs > 0)
        {
            System.out.print("Refining(" + num_refs + "): ");
            refiner = F_GlobalDomainRefiner.invoke(dom);
	
            // performing post refs
            for (int i=0; i<num_refs; i++)
            {
                refiner.refine();
                System.out.print(i + " ");
            }
            
            System.out.print("done.\n");
        }
        */
        
        
        // create, load, refine and distribute domain
        System.out.println("Create, refine and distribute domain");
        String[] neededSubsets = {};
        String distributionMethod = "metisReweigh";
        I_InterSubsetPartitionWeighting weightingFct = new InterSubsetPartitionWeighting();
        weightingFct.set_default_weights(1,1);
        weightingFct.set_inter_subset_weight(0, 1, 1000);
        I_Domain dom = createAndDistributeDomain(gridFileName, num_refs, num_prerefs, neededSubsets, distributionMethod, true, -1, -1, weightingFct);
        
   
    //  create approximation space
        System.out.println("Create approximation space");
        I_ApproximationSpace approxSpace = new ApproximationSpace(dom);
        
    //  defining approximation space according to function definition
        for (FunctionDefinition fd: problemDefinition)
        {
            String subsets = "";
            if (fd.getFctData().subsetList.isEmpty())
                throw new RuntimeException("No subset definition for function '"+fd.getFctData().fctName+"'!");
            
            for (String s: fd.getFctData().subsetList) subsets = subsets + ", " + s;
            subsets = subsets.substring(2);
            approxSpace.add_fct(fd.getFctData().fctName, "Lagrange", 1, subsets);
        }

        approxSpace.init_levels();
        approxSpace.const__print_layout_statistic();
        approxSpace.const__print_statistic();
        
        
    //////////////////////////
    //  Discretization setup 
    //////////////////////////
        
        I_DomainDiscretization domainDisc = new DomainDiscretization(approxSpace);
        
        
    ////////////////////////////////////////////////////
    //  Diffusion processes
        I_ConvectionDiffusionFV1[] diffusionDisc = new ConvectionDiffusionFV1[diffusionData.length];
        for (int i = 0; i < diffusionData.length; i++)
        { 
        //  get selected function and selected subsets
            UserDependentSubsetModel.FSDataType fctSsSel = (UserDependentSubsetModel.FSDataType) diffusionData[i].getData(0);
            if (fctSsSel.getSelFct().length != 1) throw new RuntimeException("Diffusion process "+i+" needs exactly one function, but has "+fctSsSel.getSelFct().length+".");
            String fct = fctSsSel.getSelFct()[0];
            String[] ss = fctSsSel.getSelSs();
            String ssString = "";
            if (ss.length == 0) throw new RuntimeException("No subset definition for function '"+fct+"' in diffusion process "+i+"!");
            for (String s: ss) ssString = ssString + ", " + s;
            ssString = ssString.substring(2);
            
        //  create elemDisc
            diffusionDisc[i] = new ConvectionDiffusionFV1();
            diffusionDisc[i].constructor(fct, ssString);

        //  upwinding not needed, no convection
            I_IConvectionShapes upwind = new NoUpwind();
            diffusionDisc[i].set_upwind(upwind);
            
        //  get parameters for diffusion / reaction
            I_CplUserMatrix diffTensor = (I_CplUserMatrix) diffusionData[i].getMatrixData(1);
            I_CplUserNumber reactionRate = (I_CplUserNumber) diffusionData[i].getNumberData(2);
            I_CplUserNumber reactionTerm = (I_CplUserNumber) diffusionData[i].getNumberData(3);
            
            diffusionDisc[i].set_diffusion(diffTensor);
            diffusionDisc[i].set_reaction_rate(reactionRate);
            diffusionDisc[i].set_reaction(reactionTerm);
            
        //  add to domain discretization
            domainDisc.add(diffusionDisc[i]);
        }      
        
        
    ////////////////////////////////////////////////////
    //  Buffering processes
        I_BufferFV1[] bufferingDisc = new BufferFV1[bufferingData.length];
        for (int i = 0; i < bufferingData.length; i++)
        { 
        //  get selected function and selected subsets
            UserDependentSubsetModel.FSDataType fctSsSel = (UserDependentSubsetModel.FSDataType) bufferingData[i].getData(0);
            String[] fcts = fctSsSel.getSelFct();
            if (fcts.length != 2) throw new RuntimeException("Buffering reaction "+i+" needs exactly two functions, but has "+fcts.length+".");
            String[] ss = fctSsSel.getSelSs();
            String ssString = "";
            if (ss.length == 0) throw new RuntimeException("No subset definition in buffering reaction "+i+"!");
            for (String s: ss) ssString = ssString + ", " + s;
            ssString = ssString.substring(2);
            
        //  create elemDisc
            bufferingDisc[i] = new BufferFV1();
            bufferingDisc[i].constructor(ssString);
            
            bufferingDisc[i].add_reaction(fcts[0], fcts[1],
                            (I_CplUserNumber) bufferingData[i].getNumberData(1),
                            (I_CplUserNumber) bufferingData[i].getNumberData(2),
                            (I_CplUserNumber) bufferingData[i].getNumberData(3));
            
        //  add to domain discretization
            domainDisc.add(bufferingDisc[i]);
        }
        
        
    ////////////////////////////////////////////////////
    //  PM membrane transport mechanisms
        
    //  VGCC
        String[] vgccSelFcts = ((UserDependentSubsetModel.FSDataType) vgccData.getData(0)).getSelFct();
        if (vgccSelFcts.length != 2) throw new RuntimeException("Voltage-gated channel mechanism needs exactly two functions, but has "+vgccSelFcts.length+".");
        String vgccFcts = vgccSelFcts[0] + ", " + vgccSelFcts[1];
        
        String[] vgccSelSs = ((UserDependentSubsetModel.FSDataType) vgccData.getData(0)).getSelSs();
        String vgccSsString = "";
        if (vgccSelSs.length == 0) throw new RuntimeException("No subset definition in voltage-gated channel definition!");
        for (String s: vgccSelSs) vgccSsString = vgccSsString + ", " + s;
        vgccSsString = vgccSsString.substring(2);
        
        I_CplUserNumber vgccDensityFct = (I_CplUserNumber) vgccData.getNumberData(1);
                       
        //  Channel object setup                
            double basicVoltage         = -65.0;        
            int VGCC_numTimeSteps       = (int)(stimTimeEnd / VGCC_dt);         
            int AP_duration_in_ms 	= (int)(1000 / stimFrequency);                          
        
            I_VoltageGatedChannels CalciumChannels = new VoltageGatedChannels();     

            if ("N".equals(vgccChannelType)) CalciumChannels.install_ca_N_gates_for_CFP_model(2e-12);
            else if ("L".equals(vgccChannelType)) CalciumChannels.install_ca_L_gates_for_CFP_model(3e-9);
            else if ("T".equals(vgccChannelType)) CalciumChannels.install_ca_T_gates_for_CFP_model(3e-9);

            CalciumChannels.calc_initial_gating(-65.0);
            CalciumChannels.calc_current(0.0, 0.0001, -65.0, 5e-5, 1.5);               

            for(int i = 1; i <= VGCC_numTimeSteps; i++)
            {
                double time = i*VGCC_dt;
                double V = F_Membrane_potential.invoke(time, stimTimeBegin, stimTimeEnd, AP_duration_in_ms, basicVoltage);
                CalciumChannels.update_gating(VGCC_dt, V);
            }

        I_FV1VGCCBoundary vgccDisc = new FV1VGCCBoundary();
        vgccDisc.constructor(vgccFcts, vgccSsString, CalciumChannels, VGCC_dt, stimTimeBegin, stimTimeEnd, (double)(AP_duration_in_ms), basicVoltage);
        vgccDisc.set_CFP_model();
        vgccDisc.set_density_function(vgccDensityFct);
        
        domainDisc.add(vgccDisc);
           
        
    //  PMCA
        String[] pmcaSelFcts = ((UserDependentSubsetModel.FSDataType) pmcaData.getData(0)).getSelFct();
        if (pmcaSelFcts.length != 2) throw new RuntimeException("PMCA pump mechanism needs exactly two functions, but has "+pmcaSelFcts.length+".");
        String pmcaFcts = pmcaSelFcts[0] + ", " + pmcaSelFcts[1];        
        
        String[] pmcaSelSs = ((UserDependentSubsetModel.FSDataType) pmcaData.getData(0)).getSelSs();
        String pmcaSsString = "";
        if (pmcaSelSs.length == 0) throw new RuntimeException("No subset definition in PMCA pump definition!");
        for (String s: pmcaSelSs) pmcaSsString = pmcaSsString + ", " + s;
        pmcaSsString = pmcaSsString.substring(2);
        
        I_CplUserNumber pmcaDensityFct = (I_CplUserNumber) pmcaData.getNumberData(1);
        
        //  PMCA object setup
            String PMCA_isoform = "h2b";
            double f_stim = 0.0;
            I_PMCA PMCA = new PMCA();
            PMCA.constructor(PMCA_isoform);

        I_FV1PMCABoundary pmcaDisc = new FV1PMCABoundary();
        pmcaDisc.constructor(pmcaFcts, pmcaSsString, PMCA, f_stim);
        pmcaDisc.set_density_function(pmcaDensityFct);
        
        domainDisc.add(pmcaDisc);
        
        
    //  NCX
        String[] ncxSelFcts = ((UserDependentSubsetModel.FSDataType) ncxData.getData(0)).getSelFct();
        if (ncxSelFcts.length != 2) throw new RuntimeException("NCX mechanism needs exactly two functions, but has "+ncxSelFcts.length+".");
        String ncxFcts = ncxSelFcts[0] + ", " + ncxSelFcts[1];        
        
        String[] ncxSelSs = ((UserDependentSubsetModel.FSDataType) ncxData.getData(0)).getSelSs();
        String ncxSsString = "";
        if (ncxSelSs.length == 0) throw new RuntimeException("No subset definition in NCX definition!");
        for (String s: ncxSelSs) ncxSsString = ncxSsString + ", " + s;
        ncxSsString = ncxSsString.substring(2);
        
        I_CplUserNumber ncxDensityFct = (I_CplUserNumber) ncxData.getNumberData(1);        
        
        //  NCX object setup
            I_NCX NCX = new NCX();
               
        I_FV1NCXBoundary ncxDisc = new FV1NCXBoundary();
        ncxDisc.constructor(ncxFcts, ncxSsString, NCX);
        ncxDisc.set_density_function(ncxDensityFct);
        
        domainDisc.add(ncxDisc);
        
               
    ////////////////////////////////////////////////////
    /// start value

    //  check that every function has been initialized on each of its subsets
        for (FunctionDefinition fd: problemDefinition)
        {
            // construct list of all subset lists in initial solution definition
            // for the given function
            List<List<String>> dssll = new ArrayList<List<String>>();
            for (UserDataTuple udt: startValue)
            {
                if (((UserDependentSubsetModel.FSDataType) udt.getData(0)).getSelFct().length != 1)
                {
                    throw new RuntimeException("Start value definition needs exactly one function at a time, but has "
                        + ((UserDependentSubsetModel.FSDataType) udt.getData(0)).getSelFct().length+".");
                }
                String fct = ((UserDependentSubsetModel.FSDataType) udt.getData(0)).getSelFct()[0];

                if (fct.equals(fd.getFctData().fctName))
                    dssll.add(fd.getFctData().subsetList);
            }

            // check that each item on the definition subset list
            // is somewhere in the initial solution subset list list
            List<String> dssl = fd.getFctData().subsetList;
            for (String dss: dssl)
            {
                boolean found = false;
                outerSearchLoop:
                for (List<String> issl: dssll)
                {
                    for (String iss: issl)
                        if (iss.equals(dss))
                        {
                            found = true;
                            break outerSearchLoop;
                        }
                }
                if (!found)
                {
                    throw new RuntimeException("Initialization value for function '"
                            +fd.getFctData().fctName+"' on subset '"+dss+"' is missing!");
                }
            }
        }
        
        return new Object[]{domainDisc, approxSpace, startValue, stimTimeBegin, stimTimeEnd, dt, VGCC_dt};
    }
    
    
    // The following functions are mere java equivalents of lua script functions
    // defined in ug_util.lua and domain_distribution_util.lua.
    
    private I_Domain createAndDistributeDomain(String gridName, int numRefs, int numPreRefs, String[] neededSubsets, String distributionMethod, boolean verticalInterfaces, int numTargetProcs, int distributionLevel, I_PartitionWeighting wFct)
    {
        // create Instance of a Domain
	I_Domain dom = new Domain();
	
	// load domain
        System.out.print("Loading domain " + gridName + " ... ");
        F_LoadDomain.invoke(dom, gridName);
	System.out.print("done.\n");
        
	// create Refiner
	if (numPreRefs > numRefs)
            throw new RuntimeException("numPreRefs must be smaller than numRefs.");
	
	if (numPreRefs > numRefs)  numPreRefs = numRefs;
	
	// Create a refiner instance. This is a factory method
	// which automatically creates a parallel refiner if required.
	I_IRefiner refiner = null;
	if (numRefs > 0)
        {
            System.out.print("Refining(" + numRefs + "): ");
            refiner = F_GlobalDomainRefiner.invoke(dom);
	
            // performing pre-refines
            for (int i=0; i<numPreRefs; i++)
            {
                refiner.refine();
                System.out.print(i + " ");
            }
            
            System.out.print("done.\n");
        }
	
        // distribute the domain to all involved processes
	if (!distributeDomain(dom, distributionMethod, verticalInterfaces, numTargetProcs, distributionLevel, wFct))
            throw new RuntimeException("Error while Distributing Grid.");
	
	// perform post-refine
	for (int i=numPreRefs; i<numRefs; i++) refiner.refine();
	
	// loop all subsets and search for them
        // in the SubsetHandler of the domain
	if (neededSubsets != null)
        {
            if (!checkSubsets(dom, neededSubsets))
                throw new RuntimeException("Something wrong with required subsets. Aborting.");
        }
	
	// return the created domain
	return dom;
    }
    
    
    private boolean distributeDomain(I_Domain dom, String partitioningMethod, boolean verticalInterfaces, int numTargetProcs, int distributionLevel, I_PartitionWeighting wFct)
    {
        if (F_NumProcs.invoke() == 1) return true;
	
	I_PartitionMap partitionMap = new PartitionMap();
        
	if (partitioningMethod == null) partitioningMethod = "bisection";
	
	if (numTargetProcs <= 0) numTargetProcs = F_NumProcs.invoke();
        
	if (distributionLevel < 0)
        {
		distributionLevel = dom.grid().const__num_levels() - 1;
		if (distributionLevel < 0) distributionLevel = 0;
        }
	
	if (dom.const__domain_info().const__num_elements_on_level(distributionLevel) < numTargetProcs)
        {
            System.out.println("\nWARNING in DistributeDomain:");
            System.out.println("    There are less elements on distributionLevel than there are target processes.");
            System.out.println("    If ug hangs during parallel execution, consider increasing numPreRefs to avoid this!");
            System.out.println("    num elements on level " + distributionLevel + ": "
                + dom.const__domain_info().const__num_elements_on_level(distributionLevel));
            System.out.println("    num target processes: " + numTargetProcs);
            System.out.println("");
        }
	
        /*
	if ("bisection".equals(partitioningMethod))
        {
            if (distributionLevel < dom.grid().const__num_levels() - 1)
            {
                System.out.println("WARNING in util.DistributeDomain: 'bisection' can currently "
                    + "only be performed on the top level. Sorry...");
            }
            partitionMapBisection(dom, partitionMap, numTargetProcs);
        }
        */        
        //else if ("metis".equals(partitioningMethod))
        if ("metis".equals(partitioningMethod))
        {
            partitionMapMetis(dom, partitionMap, numTargetProcs, distributionLevel);
        }
        else if ("metisReweigh".equals(partitioningMethod))
        {
            if (wFct != null)
                partitionMapMetisReweigh(dom, partitionMap, numTargetProcs, distributionLevel, wFct);
            else 
            {
                System.out.println("ERROR in CalciumDynamics::distributeDomain: "
                        + "requested partitionMethod \"metisReweigh\", but no weightingFct given.");
                return false;
            }
        }
        else
        {
            System.out.println("ERROR in util.DistributeDomain: Unknown partitioning method.\n"
                + "  Valid partitioning methods are: 'bisection', 'metis' and 'metisReweigh'");
            return false;
        }
	
	boolean success = F_DistributeDomain.invoke(dom, partitionMap, verticalInterfaces);
	
	return success;
    }   
    
    
    /*
    private void partitionMapBisection(I_Domain dom, I_PartitionMap partitionMapOut, int numProcs)
    {
        if (partitionMapOut.num_target_procs() != numProcs)
        {
            partitionMapOut.clear();
            partitionMapOut.add_target_procs(0, numProcs);
        }
        
        I_ProcessHierarchy procH = new ProcessHierarchy();
	if (dom.grid().const__num_levels() > 0)
            procH.add_hierarchy_level(dom.grid().const__num_levels() - 1, numProcs);
	else
            procH.add_hierarchy_level(0, numProcs);
        
        if (dom.const__domain_info().const__element_type() == dom.const__get_dim() - 2)
        {
            I_HyperManifoldPartitioner_DynamicBisection partitioner = new HyperManifoldPartitioner_DynamicBisection(dom);
            partitioner.enable_clustered_siblings(false);
            partitioner.set_verbose(false);
            partitioner.enable_static_partitioning(true);
            partitioner.set_subset_handler(partitionMapOut.get_partition_handler());
            partitioner.set_next_process_hierarchy(procH);
            partitioner.partition(0, 0);
        }
        else if (dom.const__domain_info().const__element_type() == dom.const__get_dim() - 1)
        {
            I_ManifoldPartitioner_DynamicBisection partitioner = new ManifoldPartitioner_DynamicBisection(dom);
            partitioner.enable_clustered_siblings(false);
            partitioner.set_verbose(false);
            partitioner.enable_static_partitioning(true);
            partitioner.set_subset_handler(partitionMapOut.get_partition_handler());
            partitioner.set_next_process_hierarchy(procH);
            partitioner.partition(0, 0);
        }
        else if (dom.const__domain_info().const__element_type() == dom.const__get_dim())
        {
            Partitioner_DynamicBisection partitioner = new Partitioner_DynamicBisection(dom);
            partitioner.enable_clustered_siblings(false);
            partitioner.set_verbose(false);
            partitioner.enable_static_partitioning(true);
            partitioner.set_subset_handler(partitionMapOut.get_partition_handler());
            partitioner.set_next_process_hierarchy(procH);
            partitioner.partition(0, 0);
        }
    }
    */
    
    private void partitionMapMetis(I_Domain dom, I_PartitionMap partitionMapOut, int numProcs, int baseLevel)
    {
        if (partitionMapOut.num_target_procs() != numProcs)
        {
            partitionMapOut.clear();
            partitionMapOut.add_target_procs(0, numProcs);
        }
        F_PartitionDomain_MetisKWay.invoke(dom, partitionMapOut, numProcs, baseLevel, 1, 1);
    }
    
    private void partitionMapMetisReweigh(I_Domain dom, I_PartitionMap partitionMapOut, int numProcs, int baseLevel, I_PartitionWeighting wFct)
    {
        if (partitionMapOut.num_target_procs() != numProcs)
        {
            partitionMapOut.clear();
            partitionMapOut.add_target_procs(0, numProcs);
        }
	F_PartitionDomain_MetisKWay.invoke(dom, partitionMapOut, numProcs, baseLevel, wFct);
    }
    
    private boolean checkSubsets(I_Domain dom, String[] neededSubsets)
    {
        I_MGSubsetHandler sh = dom.subset_handler();
	for (String s: neededSubsets)
        {
            if (sh.const__get_subset_index(s) == -1)
            {
                System.out.print("Domain does not contain subset '" + s + "'.");
                return false;
            }
        }
	
	return true;
    }
    
    
    private void errorExit(String s)
    {
        
        eu.mihosoft.vrl.system.VMessage.exception("Setup Error in DNMJ: ", s);
    }
}
