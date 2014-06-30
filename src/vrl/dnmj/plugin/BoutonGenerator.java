/**
 * VRL-DNMJ-Plugin
 * BoutonGenerator.java
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


@ComponentInfo(name="BoutonGenerator", category="VRL-DNMJ-Plugin", description="Problem setup component")


public class BoutonGenerator implements java.io.Serializable{
    private static final long serialVersionUID = 1L;

    
    public void BuildBouton(@ParamInfo(name="file", style="save-dialog", options="tag=\"TheFile\"") java.io.File file,
                            @ParamInfo(name="radius [um]", options="value=1.0") double radius,
                            @ParamInfo(name="numRefinements", options="value=3") int numRefinements,
                            @ParamInfo(name="numReleaseSites", options="value=12") int numReleaseSites,
                            @ParamInfo(name="TbarHeight [um]" , options="value=0.07") double TbarHeight,
                            @ParamInfo(name="TbarLegRadius [um]", options="value=0.02") double TbarLegRadius,
                            @ParamInfo(name="TbarTopRadius [um]", options="value=0.15") double TbarTopRadius,
                            @ParamInfo(name="TbarTopHeight [um]", options="value=0.015") double TbarTopHeight) {
                              
        
        String fileName = file.getAbsoluteFile().getAbsolutePath();
        
        edu.gcsc.vrl.ug.api.F_BuildBouton.invoke(radius, 
                                                  numRefinements, 
                                                  numReleaseSites, 
                                                  TbarHeight, 
                                                  TbarLegRadius, 
                                                  TbarTopRadius, 
                                                  TbarTopHeight,
                                                  fileName);
        
        //return "Congratulations! "
        //        + "Your new project works!" + " " + a;
    }
    

    
}
