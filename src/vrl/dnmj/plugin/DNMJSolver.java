/**
 * VRL-DNMJ-Plugin
 * DNMJSolver.java
 * 
 * @date 2014-06-18
 * @author mstepnie
**/

package vrl.dnmj.plugin;

import edu.gcsc.vrl.ug.api.AssembledOperator;
import edu.gcsc.vrl.ug.api.BiCGStab;
import edu.gcsc.vrl.ug.api.CG;
import edu.gcsc.vrl.ug.api.CompositeConvCheck;
import edu.gcsc.vrl.ug.api.ConvCheck;
import edu.gcsc.vrl.ug.api.F_ComputeVolume;
import edu.gcsc.vrl.ug.api.F_Interpolate;
import edu.gcsc.vrl.ug.api.F_Print;
import edu.gcsc.vrl.ug.api.F_TakeMeasurement;
import edu.gcsc.vrl.ug.api.F_VecScaleAssign;
import edu.gcsc.vrl.ug.api.GaussSeidel;
import edu.gcsc.vrl.ug.api.GeometricMultiGrid;
import edu.gcsc.vrl.ug.api.GridFunction;
import edu.gcsc.vrl.ug.api.ILU;
import edu.gcsc.vrl.ug.api.I_ApproximationSpace;
import edu.gcsc.vrl.ug.api.I_AssembledOperator;
import edu.gcsc.vrl.ug.api.I_CompositeConvCheck;
import edu.gcsc.vrl.ug.api.I_ConvCheck;
import edu.gcsc.vrl.ug.api.I_CplUserNumber;
import edu.gcsc.vrl.ug.api.I_DomainDiscretization;
//import edu.gcsc.vrl.ug.api.I_FV1BorgGrahamWithVM2UG;
import edu.gcsc.vrl.ug.api.I_GridFunction;
import edu.gcsc.vrl.ug.api.I_ILinearIterator;
import edu.gcsc.vrl.ug.api.I_ILinearOperatorInverse;
import edu.gcsc.vrl.ug.api.I_IPreconditionedLinearOperatorInverse;
import edu.gcsc.vrl.ug.api.I_MGSubsetHandler;
import edu.gcsc.vrl.ug.api.I_NewtonSolver;
import edu.gcsc.vrl.ug.api.I_SolutionTimeSeries;
//import edu.gcsc.vrl.ug.api.I_Test;
import edu.gcsc.vrl.ug.api.I_ThetaTimeStep;
import edu.gcsc.vrl.ug.api.I_VTKOutput;
import edu.gcsc.vrl.ug.api.I_VTest;
import edu.gcsc.vrl.ug.api.I_Vector;
import edu.gcsc.vrl.ug.api.Jacobi;
import edu.gcsc.vrl.ug.api.LU;
import edu.gcsc.vrl.ug.api.LinearSolver;
import edu.gcsc.vrl.ug.api.NewtonSolver;
import edu.gcsc.vrl.ug.api.SolutionTimeSeries;
//import edu.gcsc.vrl.ug.api.Test;
import edu.gcsc.vrl.ug.api.ThetaTimeStep;
import edu.gcsc.vrl.ug.api.VTKOutput;
import edu.gcsc.vrl.userdata.UserDataTuple;
import edu.gcsc.vrl.userdata.UserDependentSubsetModel;
import eu.mihosoft.vrl.annotation.ComponentInfo;
import eu.mihosoft.vrl.annotation.MethodInfo;
import eu.mihosoft.vrl.annotation.OutputInfo;
import eu.mihosoft.vrl.annotation.ParamGroupInfo;
import eu.mihosoft.vrl.annotation.ParamInfo;
import eu.mihosoft.vrl.math.Trajectory;
import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author mstepnie, adapted from mbreit & A. Vogel
 */
@ComponentInfo(name="DNMJSolver", category="VRL-DNMJ-Plugin")
public class DNMJSolver implements Serializable
{
    private static final long serialVersionUID = 1L;
    private boolean stopSolver;
    private Trajectory[] vEvalTrajectory;
    
    @MethodInfo(valueStyle="multi-out", interactive = false)
    @OutputInfo
    (
      style="multi-out",
      elemNames = {"PVD File"},
      elemTypes = {File.class}
    )
    public Object[] invoke
    (
    //  PROBLEM SETUP
        @ParamGroupInfo(group="PROBLEM SETUP|false")
        @ParamInfo(name="Domain Disc", style="default")
        I_DomainDiscretization domainDisc,            
            
        @ParamGroupInfo(group="PROBLEM SETUP|false")
        @ParamInfo(name="Approximation Space", style="default")
        I_ApproximationSpace approxSpace,
        
        @ParamGroupInfo(group="PROBLEM SETUP|false")
        @ParamInfo(name="Initial Solution", style="default")
        UserDataTuple[] initValues,        
        
    //  TIME STEPPING
        @ParamGroupInfo(group="TIME STEPPING|false")
        @ParamInfo(name="Stimulation start [s]", style="default")
        double stimTimeBegin,
        
        @ParamGroupInfo(group="TIME STEPPING|false")
        @ParamInfo(name="Stimulation end [s]", style="default")
        double stimTimeEnd,
        
        @ParamGroupInfo(group="TIME STEPPING|false")
        @ParamInfo(name="Simulation time step size dt [s]", style="default")
        double user_specified_dt,
                
        @ParamGroupInfo(group="TIME STEPPING|false")
        @ParamInfo(name="VGCC time step size VGCC_dt [s]", style="default")
        double VGCC_dt,
        
        @ParamGroupInfo(group="TIME STEPPING|false|Time Solver Parameters")
        @ParamInfo(name="Time Solver", style="selection", options="value=[\"Implicit Euler\", \"Crank-Nicolson\", \"Explicit Euler\"]")
        String timeSolverName,
        
        @ParamGroupInfo(group="TIME STEPPING|false")
        @ParamInfo(name="Simulation start Time [s]", style="default", options="value=0.00D")
        double timeStart,
        
        @ParamGroupInfo(group="TIME STEPPING|false")
        @ParamInfo(name="Simulation end Time [s]", style="default", options="value=0.1D")
        double timeEnd,
               
    //  NONLINEAR SOLVER
        @ParamGroupInfo(group="NONLINEAR SOLVER SETUP|false")
        @ParamInfo(name="Maximal Number Iterations", style="default", options="value=30")
        int maxNumIterNonLinear,
        
        @ParamGroupInfo(group="NONLINEAR SOLVER SETUP|false")
        @ParamInfo(name="Minimal Defect Reduction", style="default", options="value=1e-7D")
        double minRedNonLinear,
        
        @ParamGroupInfo(group="NONLINEAR SOLVER SETUP|false")
        @ParamInfo(name="Minimal Residual Norm", style="default", options="value=1e-15D")
        double absTolNonLinear,

    //  linear solver
        @ParamGroupInfo(group="NONLINEAR SOLVER SETUP|false; Linear Solver Setup|false|Solver Parameters")
        @ParamInfo(name="Linear Solver", style="selection", options="value=[\"Bi-CGSTAB\", \"Linear Solver\", \"CG\", \"LU\"]")
        String solverName,
        
        @ParamGroupInfo(group="NONLINEAR SOLVER SETUP|false; Linear Solver Setup|false; Preconditioning|false")
        @ParamInfo(name="Preconditioner", style="selection", options="value=[\"GMG\", \"Jacobi\", \"Gauss-Seidel\", \"ILU\", \"None\"]")
        String precondName,
        
        @ParamGroupInfo(group="NONLINEAR SOLVER SETUP|false; Linear Solver Setup; Preconditioning|false; GMG options|false|GMG")
        @ParamInfo(name="Base Level", style="default", options="value=0")
        int gmgBaseLevel,
        
        @ParamGroupInfo(group="NONLINEAR SOLVER SETUP|false; Linear Solver Setup|false; Preconditioning|false; GMG options|false|GMG")
        @ParamInfo(name="Base Solver", style="selection", options="value=[\"Bi-CGStab\", \"Jacobi\", \"Gauss-Seidel\", \"ILU\", \"LU\"]")
        String gmgBaseSolver,
        
        @ParamGroupInfo(group="NONLINEAR SOLVER SETUP|false; Linear Solver Setup|false; Preconditioning|false; GMG options|false")
        @ParamInfo(name="Maximal Number Iterations", style="default", options="value=1000")
        int maxNumIterBase,
        
        @ParamGroupInfo(group="NONLINEAR SOLVER SETUP|false; Linear Solver Setup|false; Preconditioning|false; GMG options|false")
        @ParamInfo(name="Minimal Defect Reduction", style="default", options="value=1E-10")
        double minRedBase,
        
        @ParamGroupInfo(group="NONLINEAR SOLVER SETUP|false; Linear Solver Setup|false; Preconditioning|false; GMG options|false")
        @ParamInfo(name="Minimal Residual Norm", style="default", options="value=1E-10")
        double absTolBase,
        
        @ParamGroupInfo(group="NONLINEAR SOLVER SETUP|false; Linear Solver Setup|false; Preconditioning|false; GMG options|false")
        @ParamInfo(name="Smoother", style="selection", options="value=[\"ILU\", \"Jacobi\", \"Gauss-Seidel\"]")
        String gmgSmoother,
        
        @ParamGroupInfo(group="NONLINEAR SOLVER SETUP|false; Linear Solver Setup|false; Preconditioning|false; GMG options|false|GMG")
        @ParamInfo(name="Cycle Type", style="selection", options="value=[\"V\", \"W\"]")
        String gmgCycleType,
        
        @ParamGroupInfo(group="NONLINEAR SOLVER SETUP|false; Linear Solver Setup|false; Preconditioning|false; GMG options|false|GMG")
        @ParamInfo(name="# Presmoothing", style="default", options="value=5")
        int gmgNumPreSmooth,
        
        @ParamGroupInfo(group="NONLINEAR SOLVER SETUP|false; Linear Solver Setup|false; Preconditioning|false; GMG options|false|GMG")
        @ParamInfo(name="# Postsmoothing", style="default", options="value=5")
        int gmgnumPostSmooth,
        
        @ParamGroupInfo(group="NONLINEAR SOLVER SETUP|false; Linear Solver Setup|false")
        @ParamInfo(name="Maximal Number Iterations", style="default", options="value=200")
        int maxNumIterLinear,
        
        @ParamGroupInfo(group="NONLINEAR SOLVER SETUP|false; Linear Solver Setup|false")
        @ParamInfo(name="Minimal Defect Reduction", style="default", options="value=1e-6")
        double minRedLinear,
        
        @ParamGroupInfo(group="NONLINEAR SOLVER SETUP|false; Linear Solver Setup|false")
        @ParamInfo(name="Minimal Residual Norm", style="default", options="value=1e-22")
        double absTolLinear,
        
        // OUTPUT
        @ParamGroupInfo(group="OUTPUT|false")
        @ParamInfo(name="Output path", style="save-folder-dialog")
        String outputPath,
        
        // measurements
        @ParamGroupInfo(group="OUTPUT|false; Measurements|false")
        @ParamInfo(name="Integration Subset", style="array", options="fct_tag=\"fctDef\"; minArraySize=0; type=\"S1:function & subset\"")
        UserDataTuple[] meas,
        
        // plotting
        @ParamGroupInfo(group="OUTPUT|false; VTK|false")
        @ParamInfo(name="do plot")
        boolean generateVTKoutput,
        
        @ParamGroupInfo(group="OUTPUT|false; VTK|false")
        @ParamInfo(name="Plotting step", style="default", options="value=0.001")
        double plotStep
    )
    {
        //if (bEraseOldFiles) eraseAllFilesInFolder(fileOut, "vtu");
          
    //  set abortion flag to false initially (can be changed using stopSolver-Method)
        stopSolver = false;
        
        
        //////////////
        // Operator //
        //////////////
        
    //  create time discretization
        I_ThetaTimeStep timeDisc = new ThetaTimeStep(domainDisc);
        if      ("Implicit Euler".equals(timeSolverName)) timeDisc.set_theta(1.0);
        else if ("Crank-Nicolson".equals(timeSolverName)) timeDisc.set_theta(0.5);
        else if ("Explicit Euler".equals(timeSolverName)) timeDisc.set_theta(0.0);
        else errorExit("Cannot find time solver: " + timeSolverName);

    //  create operator from discretization
        I_AssembledOperator op = new AssembledOperator();
        op.set_discretization(timeDisc);
        op.init();
        
        //I_MGSubsetHandler sh = dom.subset_handler();

        
        /////////////
        // Algebra //
        /////////////

    //  create Convergence Check for linear solver
        I_ConvCheck convCheckLinear = new ConvCheck();
        convCheckLinear.set_maximum_steps(maxNumIterLinear);
        convCheckLinear.set_minimum_defect(absTolLinear);
        convCheckLinear.set_reduction(minRedLinear);

    //  create solver (with preconditioner)
        I_ILinearOperatorInverse solver = null;
        if      ("Bi-CGSTAB".equals(solverName)) solver = new BiCGStab();
        else if ("CG".equals(solverName)) solver = new CG();
        else if ("Linear Solver".equals(solverName)) solver = new LinearSolver();
        else if ("LU".equals(solverName)) solver = new LU();
        else errorExit("Cannot find solver: " + solverName);
        solver.set_convergence_check(convCheckLinear);

    //  create preconditioner (if suitable)
        if (!"None".equals(precondName) && !"LU".equals(solverName))
        {
            I_ILinearIterator precond = null;
 
            if      ("Jacobi".equals(precondName)) precond = new Jacobi();
            else if ("Gauss-Seidel".equals(precondName)) precond = new GaussSeidel();
            else if ("ILU".equals(precondName)) precond = new ILU();
            else if ("GMG".equals(precondName))
            {
            //  create GMG according to settings
                precond = new GeometricMultiGrid(approxSpace);
                ((GeometricMultiGrid) precond).set_discretization(timeDisc);
                ((GeometricMultiGrid) precond).set_base_level(gmgBaseLevel);
                
            //  conv check for base
                I_ConvCheck convCheckBase = new ConvCheck();
                convCheckBase.set_maximum_steps(maxNumIterBase);
                convCheckBase.set_minimum_defect(absTolBase);
                convCheckBase.set_reduction(minRedBase);
                
            //  base solver
                I_ILinearOperatorInverse baseSolver;
                if ("LU".equals(gmgBaseSolver)) baseSolver = new LU();
                else if ("Bi-CGStab".equals(gmgBaseSolver)) {baseSolver = new BiCGStab(); ((I_IPreconditionedLinearOperatorInverse)baseSolver).set_preconditioner(new ILU());}
                else
                {
                    baseSolver = new LinearSolver();
                    if      ("Jacobi".equals(gmgBaseSolver))
                        ((I_IPreconditionedLinearOperatorInverse)baseSolver).set_preconditioner(new Jacobi());
                    else if ("Gauss-Seidel".equals(gmgBaseSolver))
                        ((I_IPreconditionedLinearOperatorInverse)baseSolver).set_preconditioner(new GaussSeidel());
                    else if ("ILU".equals(gmgBaseSolver))
                        ((I_IPreconditionedLinearOperatorInverse)baseSolver).set_preconditioner(new ILU());
                    else errorExit("Cannot find preconditioner '" + gmgBaseSolver + "' for GMG.");
                }
                baseSolver.set_convergence_check(convCheckBase);
                ((GeometricMultiGrid) precond).set_base_solver(baseSolver);
                
            //  smoother
                if      ("Jacobi".equals(gmgSmoother))
                    ((GeometricMultiGrid) precond).set_smoother(new Jacobi());
                else if ("Gauss-Seidel".equals(gmgSmoother))
                    ((GeometricMultiGrid) precond).set_smoother(new GaussSeidel());
                else if ("ILU".equals(gmgSmoother))
                    ((GeometricMultiGrid) precond).set_smoother(new ILU());
                else errorExit("Cannot find smoother '" + gmgSmoother + "' for GMG.");
                
            //  cycle
                if      ("V".equals(gmgCycleType)) ((GeometricMultiGrid) precond).set_cycle_type(1);
                else if ("W".equals(gmgCycleType)) ((GeometricMultiGrid) precond).set_cycle_type(2);
                else errorExit("Cannot find cycle type '" + gmgCycleType + "' for GMG.");
                    
                ((GeometricMultiGrid) precond).set_num_presmooth(gmgNumPreSmooth);
                ((GeometricMultiGrid) precond).set_num_postsmooth(gmgnumPostSmooth);
            }
            else errorExit("Cannot find preconditioner: "+precondName);
 
            try {((I_IPreconditionedLinearOperatorInverse)solver).set_preconditioner(precond);}
            catch (Exception e)
            {
                // in case of solver is LU, there maybe a CastExeption
                // No, this cannot happen because of the 
                // if (!"LU".equals(solverName))' statement
                System.err.println(getClass().getSimpleName() +": solver is LU therefore CastException !?");
            }
        }
        approxSpace.domain().const__get_dim();
        I_CompositeConvCheck convCheckNewton = new CompositeConvCheck(approxSpace, maxNumIterNonLinear, absTolNonLinear, minRedNonLinear);
        convCheckNewton.set_verbose(true);
        convCheckNewton.set_time_measurement(true);

    //  create Newton solver
        I_NewtonSolver newtonSolver = new NewtonSolver();
        newtonSolver.set_linear_solver(solver);
        newtonSolver.set_convergence_check(convCheckNewton);
        newtonSolver.init(op);
      
        
        /////////
        // I/O //
        /////////
        
    //  append / to output path
        outputPath = outputPath + "/";
        
    //  VTK output
        I_VTKOutput vtkOut = null;
        if (generateVTKoutput) vtkOut = new VTKOutput();       
        
        
        ////////////////
        // Simulation //
        ////////////////
        
    //  start
        double time = timeStart;
        int step = 0;
        
    //  initialize solution
    //  get grid function
        F_Print.invoke("Initializing solution.\n");
        I_GridFunction u = new GridFunction(approxSpace);
        
    //  interpolate start values
        int cntUDT = 0;
        for (UserDataTuple udt: initValues)
        {
        //  get function to interpolate for
            String[] selFct = ((UserDependentSubsetModel.FSDataType) udt.getData(0)).getSelFct();
            if (selFct.length != 1) throw new RuntimeException("Start value definition needs exactly one function at a time, but has "+selFct.length+".");
            String fct = selFct[0];

        //  get subsets to interpolate for
            String[] selSs = ((UserDependentSubsetModel.FSDataType) udt.getData(0)).getSelSs();
            String ssString = "";
            if (selSs.length == 0) throw new RuntimeException("No subset selection in start value definition "+cntUDT+".");
            for (String s: selSs) ssString = ssString + ", " + s;
            ssString = ssString.substring(2);
            
        //  get start value
            I_CplUserNumber value = (I_CplUserNumber) udt.getNumberData(1);
            
        //  interpolate grid function for time
            F_Interpolate.invoke(value, u, fct, ssString, time);
            
            cntUDT++;
        }                
        
    //  write initial solution to vtk file
        if (generateVTKoutput){       
            //vtkOut.print(outputPath + "vtk/result", u, (int) Math.floor(time/plotStep+0.5), time);        
            vtkOut.print(outputPath + "vtk/result", u, step, time);       
        }        

    //  prepare measurements
        List<String> measFct = new ArrayList<String>();
        List<String> measSs = new ArrayList<String>();
        cntUDT = 0;
        for (UserDataTuple udt: meas)
        {
        //  get function to interpolate for
            String[] selFct = ((UserDependentSubsetModel.FSDataType) udt.getData(0)).getSelFct();
            if (selFct.length != 1) throw new RuntimeException("Measurement definition needs exactly one function at a time, but has "+selFct.length+".");
            measFct.add(selFct[0]);

        //  get subsets to interpolate for
            String[] selSs = ((UserDependentSubsetModel.FSDataType) udt.getData(0)).getSelSs();
            String ssString = "";
            if (selSs.length == 0) throw new RuntimeException("No subset selection in measurement definition "+cntUDT+".");
            for (String s: selSs) ssString = ssString + ", " + s;
            measSs.add(ssString.substring(2));
            
            cntUDT++;
        }
        
    //  take first measurement
        for (int i=0; i<measFct.size(); i++)
        {
            F_TakeMeasurement.invoke(u, approxSpace, time, measSs.get(i),
                                     measFct.get(i), outputPath + "meas/data");
        }
        
    //  create new grid function for old value
        I_GridFunction uOld = u.clone();
        
    //  store grid function in vector of old solutions
        I_SolutionTimeSeries solTimeSeries = new SolutionTimeSeries();
        solTimeSeries.push(uOld, time);
        
        I_MGSubsetHandler sh = approxSpace.domain().subset_handler();
        
    //  display volumes/areas of subsets
        String allSubsets = "";
        for (int i=0; i<sh.const__num_subsets(); i++)
            allSubsets = allSubsets + ", " + sh.const__get_subset_name(i);
        if (sh.const__num_subsets()>0) allSubsets = allSubsets.substring(2);
        
        F_ComputeVolume.invoke(approxSpace, allSubsets);
        
        
    //  computations for time stepping
    //  choose length of time step at the beginning
    //  if not timeStepStart = 2^(-n)*timeStep, take nearest lower number of that form
        double timeStepStart = user_specified_dt;
        double maxStepSize = timeStepStart;
        double minStepSize = VGCC_dt;
        
        int startLv =  (int) Math.ceil(log2(maxStepSize/timeStepStart));
        double timeStepStartNew = maxStepSize / Math.pow(2, startLv);
        timeStepStart = timeStepStartNew;
        double dt = timeStepStart;
        
    //  same for minStepSize
        int LowLv =  (int) Math.ceil(log2(maxStepSize/minStepSize));
        double minStepSizeNew = maxStepSize / Math.pow(2, LowLv);
        minStepSize = minStepSizeNew;
        
        int checkbackInterval = 10;
        int lv = startLv;
        double levelUpDelay = 0.0; //caEntryDuration + (nSpikes - 1) * 1.0/freq;
        int[] checkbackCounter = new int[LowLv+1];
        for (int i=0; i<checkbackCounter.length; i++) checkbackCounter[i] = 0;
        
        
    //  begin simulation loop
        while (time < timeEnd)
        {
            F_Print.invoke("++++++ POINT IN TIME  " + Math.floor((time+dt)/dt+0.5)*dt + "s  BEGIN ++++++");

        //  setup time disc for old solutions and timestep            
            timeDisc.prepare_step(solTimeSeries, dt);            

        //  prepare Newton solver
            if (!newtonSolver.prepare(u))
            {
                F_Print.invoke("Newton solver failed at point in time "
                        + Math.floor((time+dt)/dt+0.5)*dt + ".");
                errorExit("Newton solver failed at point in time "
                        + Math.floor((time+dt)/dt+0.5)*dt + ".");
            }                                    
            
        //  apply Newton solver
            if (!newtonSolver.apply(u))
            {                                
            //  in case of failure:
                F_Print.invoke("Newton solver failed at point in time "
                               + Math.floor((time+dt)/dt+0.5)*dt
                               + " with time step " + dt + ".");

                dt = dt/2;
                lv = lv++;
                F_VecScaleAssign.invoke(u, 1.0, solTimeSeries.latest());

            //  half time step and try again unless time step below minimum
                if (dt < minStepSize)
                {
                    F_Print.invoke("Time step below minimum. Aborting. Failed at point in time "
                            + Math.floor((time+dt)/dt+0.5)*dt + ".");
                    time = timeEnd;
                }
                else
                {    
                    F_Print.invoke("Trying with half the time step...");
                    checkbackCounter[lv] = 0;
                }
            }
            else
            {
            //  update new time
                time = solTimeSeries.const__time(0) + dt;

            //  update checkback counter and if applicable, reset dt
                checkbackCounter[lv]++;
                while (checkbackCounter[lv] % (2*checkbackInterval) == 0 && lv > 0 && (time >= levelUpDelay || lv > startLv))
                {
                    F_Print.invoke("Doubling time due to continuing convergence; now: " + 2*dt);
                    dt = 2*dt;
                    lv--;
                    checkbackCounter[lv] = checkbackCounter[lv] + checkbackCounter[lv+1] / 2;
                    checkbackCounter[lv+1] = 0;
                }

            //  plot solution every plotStep seconds
                if (generateVTKoutput)
                {
                    if (Math.abs(time/plotStep - Math.floor(time/plotStep+0.5)) < 1e-5)
                        vtkOut.print(outputPath + "vtk/result", u, (int) Math.floor(time/plotStep+0.5), time);
                }

            //  take measurement in bouton every timeStep seconds
                for (int i=0; i<measFct.size(); i++)
                {
                    F_TakeMeasurement.invoke(u, approxSpace, time, measSs.get(i),
                                             measFct.get(i), outputPath + "meas/data");
                }

            //  export solution of ca on mem_er
                //F_ExportSolution.invoke(u, approxSpace, time, "mem_cyt", "ca_cyt", outputPath + "sol/sol");

                // get oldest solution
                I_Vector oldestSol = solTimeSeries.oldest();

                // copy values into oldest solution (we reuse the memory here)
                F_VecScaleAssign.invoke(oldestSol, 1.0, u);

                // push oldest solutions with new values to front, oldest sol pointer is popped from end
                solTimeSeries.push_discard_oldest(oldestSol, time);

                F_Print.invoke("++++++ POINT IN TIME  " + Math.floor(time/dt+0.5)*dt + "s  END ++++++++");
                if (stopSolver)
                    break;
            }
        }
        
        // end timeseries, produce gathering file
        if (generateVTKoutput) vtkOut.write_time_pvd(outputPath + "vtk/result", u);
        
        if (generateVTKoutput) return new Object[]{new File(outputPath + "vtk/result.pvd")};
        return new Object[]{null};
    }

    @MethodInfo(name="stopSolver", buttonText="Stop time stepping", hideCloseIcon=true)   
    public void stopSolver()
    {
        stopSolver=true;
    }
   
    private void errorExit(String s)
    {
        eu.mihosoft.vrl.system.VMessage.exception("Setup Error in KineticSolver: ", s);
    }

    @MethodInfo(name="", valueName="Trajectories", valueStyle="default", valueOptions="", interactive = false)
    public Trajectory[] getTrajectories()
    {
       return vEvalTrajectory;
    }
    
    private double log2(double x)
    {
	return Math.log(x)/Math.log(2.0);
    }
    

}
