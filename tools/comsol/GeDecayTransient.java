import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;

// Full-precision solver trajectories for u'+u=0, u(0)=1. Unlike Data-export
// labels, getPVals()/getU()/getUDot() preserve the solver's actual doubles.
// All nontrivial API calls were exercised by TimeSteppingDemonstration.
public class GeDecayTransient {
    static String directory;
    static PrintWriter manifest;

    static Model build(String times) {
        Model m = ModelUtil.create("Model");
        m.component().create("comp1", true);
        m.component("comp1").geom().create("geom1", 2);
        m.component("comp1").geom("geom1").create("sq1", "Square");
        m.component("comp1").geom("geom1").run();
        m.component("comp1").physics().create("ge", "GlobalEquations", "geom1");
        m.component("comp1").physics("ge").feature("ge1").set("equation", "u+ut");
        m.component("comp1").physics("ge").feature("ge1").set("name", "u");
        m.component("comp1").physics("ge").feature("ge1").set("initialValueU", "1");
        m.component("comp1").mesh().create("mesh1");
        m.component("comp1").mesh("mesh1").autoMeshSize(9);
        m.component("comp1").mesh("mesh1").run();
        m.study().create("std1");
        m.study("std1").create("time", "Transient");
        m.study("std1").feature("time").set("tlist", times);
        m.study("std1").createAutoSequences("all");
        return m;
    }

    static void run(String name, int order, double tolerance, String mode,
                    String maxStepKind, double maxStep, boolean endInterpolation,
                    boolean exact, String times, String mph) throws Exception {
        Model m = build(times);
        SolverFeature f = m.sol("sol1").feature("t1");
        f.set("tout", "tsteps");
        f.set("keeplog", "on");
        f.set("maxorder", order);
        f.set("rtol", tolerance);
        f.set("tstepsbdf", mode);
        String store = name.equals("store_closest") ? "tstepsclosest"
            : name.equals("store_interpolate") ? "tlist" : "tsteps";
        f.set("tout", store);
        String manualStep = name.equals("manual_large") ? "0.2" : "0.1";
        if (mode.equals("manual"))
            f.set("timestepbdf", manualStep);
        f.set("endtimeinterpolation", endInterpolation ? "on" : "off");
        if (!maxStepKind.equals("auto")) {
            f.set("maxstepconstraintbdf", maxStepKind);
            f.set("maxstepbdf", maxStep);
        }
        if (exact) {
            f.feature("fc1").set("jtech", "onevery");
            f.feature("fc1").set("ntolfact", 1e-8);
        }
        m.sol("sol1").runAll();
        double[] t = m.sol("sol1").getPVals();
        PrintWriter data = new PrintWriter(new File(directory, name + ".csv"));
        data.println("time,value,derivative");
        for (int i = 0; i < t.length; ++i) {
            data.println(Double.toString(t[i]) + ","
                + Double.toString(m.sol("sol1").getU(i + 1)[0]) + ","
                + Double.toString(m.sol("sol1").getUDot(i + 1)[0]));
        }
        data.close();
        PrintWriter log = new PrintWriter(new File(directory, name + ".log"));
        log.print(f.getString("message"));
        log.close();
        manifest.println(name + "," + order + "," + tolerance + "," + mode
            + "," + maxStepKind + "," + maxStep + "," + endInterpolation + "," + exact
            + "," + store + "," + (mode.equals("manual") ? manualStep : "0")
            + ",\"" + times + "\"");
        manifest.flush();
        if (mph != null) m.save(mph);
        ModelUtil.clear();
        System.out.println("CASE_OK " + name + " samples=" + t.length);
    }

    public static void main(String[] args) throws Exception {
        directory = args != null && args.length > 1 ? args[1] : "trajectories";
        new File(directory).mkdirs();
        manifest = new PrintWriter(new File(directory, "manifest.csv"));
        manifest.println("name,max_order,tolerance,mode,max_step_kind,max_step_seconds,end_interpolation,exact,store,manual_step,times");
        String times = "0 1 2 3 4 5 6 7 8 9 10";
        run("default", 5, 0.01, "free", "auto", 0.0, true, false, times, null);
        run("max_order_1", 1, 0.01, "free", "auto", 0.0, true, false, times, null);
        run("max_order_2", 2, 0.01, "free", "auto", 0.0, true, false, times, null);
        run("max_order_3", 3, 0.01, "free", "auto", 0.0, true, false, times, null);
        run("max_order_4", 4, 0.01, "free", "auto", 0.0, true, false, times, null);
        run("rtol_1e3", 5, 0.001, "free", "auto", 0.0, true, false, times, null);
        run("rtol_1e5", 5, 0.00001, "free", "auto", 0.0, true, false, times, null);
        run("strict", 5, 0.01, "strict", "auto", 0.0, true, false, times, null);
        run("intermediate", 5, 0.01, "intermediate", "auto", 0.0, true, false, times, null);
        run("max_step_5e3", 5, 0.01, "free", "const", 0.05, true, false, times, null);
        run("end_landing", 5, 0.01, "free", "auto", 0.0, false, false, times, null);
        run("exact_default", 5, 0.01, "free", "auto", 0.0, true, true, times, null);
        run("exact_tight", 5, 0.00001, "free", "auto", 0.0, true, true, times, null);
        run("exact_strict", 5, 0.01, "strict", "auto", 0.0, false, true, times, null);
        String irregular = "0 0.001 0.01 0.05 0.2 0.3 0.7 1 1.1 1.4 2 3 5 8 10";
        run("irregular_strict", 5, 0.01, "strict", "auto", 0.0, true, false, irregular, null);
        run("irregular_intermediate", 5, 0.01, "intermediate", "auto", 0.0, true, false, irregular, null);
        run("manual_default", 5, 0.01, "manual", "auto", 0.0, true, false, times, null);
        run("manual_exact", 5, 0.01, "manual", "auto", 0.0, true, true, times, null);
        run("manual_large", 5, 0.01, "manual", "auto", 0.0, true, false, times, null);
        run("store_closest", 5, 0.01, "free", "auto", 0.0, true, false, times, null);
        run("store_interpolate", 5, 0.01, "free", "auto", 0.0, true, false, times,
            args != null && args.length > 0 ? args[0] : "GeDecayTransient.mph");
        manifest.close();
        System.out.println("ALL_CASES_OK");
    }
}
