import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.*;

// The COMSOL reference of AdaptiveTimeStepper: the model the step sequences in
// test/reference/comsol_steps.hpp were recorded from.
//
// One Global Equation, `u + u' = 0`, `u(0) = 1`, over [0, 10] s, run once per
// configuration with the solver's own settings read back from its node. The
// steps are written as the solver's stored times, one column per step, in a
// Data export; the solver log is written next to them.
//
// COMSOL's Java compiler targets an old source level: no diamond, no
// try-with-resources, no multi-catch.
public class AtsReference {

    static final String DIR = "reference";

    static Model build() {
        Model model = ModelUtil.create("Model");
        String comp = "comp1";
        model.component().create(comp, true);
        model.component(comp).geom().create("geom1", 2);
        model.component(comp).geom("geom1").create("sq1", "Square");
        model.component(comp).geom("geom1").feature("sq1").set("size", 1.0);
        model.component(comp).geom("geom1").run();
        model.component(comp).physics().create("ge", "GlobalEquations", "geom1");
        model.component(comp).physics("ge").feature("ge1").set("equation", "u+ut");
        model.component(comp).physics("ge").feature("ge1").set("name", "u");
        model.component(comp).physics("ge").feature("ge1").set("initialValueU", "1");
        model.component(comp).mesh().create("mesh1");
        model.component(comp).mesh("mesh1").autoMeshSize(9);
        model.component(comp).mesh("mesh1").run();
        model.study().create("std1");
        model.study("std1").create("time", "Transient");
        model.study("std1").feature("time").set("tlist", "range(0,1,10)");
        return model;
    }

    static void configure(Model model, Map props) throws Exception {
        for (Iterator it = props.entrySet().iterator(); it.hasNext();) {
            Map.Entry e = (Map.Entry) it.next();
            if (e.getValue() instanceof Integer)
                model.sol("sol1").feature("t1").set((String) e.getKey(),
                    ((Integer) e.getValue()).intValue());
            else if (e.getValue() instanceof Double)
                model.sol("sol1").feature("t1").set((String) e.getKey(),
                    ((Double) e.getValue()).doubleValue());
            else
                model.sol("sol1").feature("t1").set((String) e.getKey(),
                    (String) e.getValue());
        }
    }

    static Map cfg(Object[] kv) {
        Map m = new LinkedHashMap();
        for (int i = 0; i < kv.length; i += 2) m.put((String) kv[i], kv[i + 1]);
        return m;
    }

    static void export(Model model, String path) {
        model.result().export().create("data1", "Data");
        model.result().export("data1").set("data", "dset1");
        model.result().export("data1").set("filename", path);
        model.result().export("data1").set("expr", new String[] {"u"});
        model.result().export("data1").run();
    }

    static void write(String path, String text) throws IOException {
        FileWriter writer = new FileWriter(path);
        try {
            writer.write(text);
        } finally {
            writer.close();
        }
    }

    static void run(String name, Map props, boolean keepLog) {
        try {
            Model model = build();
            model.study("std1").run();
            configure(model, props);
            if (keepLog) model.sol("sol1").feature("t1").set("keeplog", "on");
            model.sol("sol1").runAll();

            String file = new File(DIR, name + ".txt").getAbsolutePath();
            export(model, file);
            System.out.println("CASE " + name + " written");

            if (keepLog) {
                String log = model.sol("sol1").feature("t1").getString("message");
                write(new File(DIR, "solver-log.txt").getAbsolutePath(), log);
                System.out.println("LOG written");
            }
            ModelUtil.clear();
        } catch (Throwable error) {
            System.out.println("CASE " + name + " FAILED: " + error);
            error.printStackTrace();
        }
    }

    public static void main(String[] args) throws Exception {
        new File(DIR).mkdirs();

        // The default run, with the solver log: the step sizes and the order of
        // every step it prints.
        run("default", cfg(new Object[] {"tout", "tsteps"}), true);
        run("max_order_1", cfg(new Object[] {"tout", "tsteps", "maxorder", new Integer(1)}), false);
        run("max_order_2", cfg(new Object[] {"tout", "tsteps", "maxorder", new Integer(2)}), false);
        run("max_order_3", cfg(new Object[] {"tout", "tsteps", "maxorder", new Integer(3)}), false);
        run("max_order_4", cfg(new Object[] {"tout", "tsteps", "maxorder", new Integer(4)}), false);
        run("rtol_1e3", cfg(new Object[] {"tout", "tsteps", "rtol", new Double(1e-3)}), false);
        run("rtol_1e5", cfg(new Object[] {"tout", "tsteps", "rtol", new Double(1e-5)}), false);
        run("strict", cfg(new Object[] {"tout", "tsteps", "tstepsbdf", "strict"}), false);
        run("intermediate", cfg(new Object[] {"tout", "tsteps", "tstepsbdf", "intermediate"}), false);
        run("max_step_5e3", cfg(new Object[] {"tout", "tsteps",
            "maxstepconstraintbdf", "const", "maxstepbdf", new Double(0.05)}), false);

        write(new File(DIR, "properties.txt").getAbsolutePath(),
            "The solver node's own values, read back after a run.\n");
        System.out.println("ALL_DONE");
    }
}
