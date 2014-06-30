/**
 * VRL-DNMJ-Plugin
 * DNMJPluginConfigurator.java
 * 
 * @date 2014-06-18
 * @author mstepnie
**/

package vrl.dnmj.plugin;

import eu.mihosoft.vrl.system.InitPluginAPI;
import eu.mihosoft.vrl.system.PluginAPI;
import eu.mihosoft.vrl.system.PluginDependency;
import eu.mihosoft.vrl.system.PluginIdentifier;
import eu.mihosoft.vrl.system.VPluginAPI;
import eu.mihosoft.vrl.system.VPluginConfigurator;

public class DNMJPluginConfigurator extends VPluginConfigurator{

public DNMJPluginConfigurator() {
    //specify the plugin name and version
   setIdentifier(new PluginIdentifier("VRL-DNMJ-Plugin", "0.1"));

   // describe the plugin
   setDescription("Plugin Description");

   // copyright info
   setCopyrightInfo("VRL-DNMJ-Plugin",
           "(c) Martin Stepniewski",
           "www.you.com", "License Name", "License Text...");
   
       addDependency(new PluginDependency("VRL", "0.4.2", "0.4.2"));
       addDependency(new PluginDependency("VRL-UG", "0.2", "0.2"));
       addDependency(new PluginDependency("VRL-UG-API", "0.2", "0.2"));
       addDependency(new PluginDependency("VRL-UserData", "0.2", "0.2"));
}

@Override
public void register(PluginAPI api) {

   // register plugin with canvas
   if (api instanceof VPluginAPI) {
       VPluginAPI vapi = (VPluginAPI) api;
       vapi.addComponent(DNMJ.class);
       vapi.addComponent(DNMJSolver.class);
       vapi.addComponent(BoutonGenerator.class);       
       vapi.addTypeRepresentation(SilentDoubleType.class);
       }
}

@Override
public void unregister(PluginAPI api) {
   // nothing to unregister
}

@Override
public void init(InitPluginAPI iApi) {
   // nothing to init
  }
}