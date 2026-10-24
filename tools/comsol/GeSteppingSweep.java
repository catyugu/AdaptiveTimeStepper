import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.*;

// Stepping-setting sweep for the same Global Equations model as
// GeDecayTransient.java (u'+u=0, u(0)=1). Every case writes name.csv
// (time,value,derivative), name.log (solver log) and one manifest.csv row.
//
// The property names below were read back from SolverFeature.properties() of
// sol1/t1 on COMSOL 6.2.0.290 instead of guessed:
//   maxstepconstraintbdf = auto | const | expr
//   maxstepbdf            constant maximum step value
//   maxstepexpressionbdf  maximum step expression, evaluated while solving
//   initialstepbdfactive  on | off, the "Initial step" check box
//   initialstepbdf        specified initial step value
//   tstepsbdf             free | manual, "Steps taken by solver"
//   timestepbdf           manual time step
//   bdforder              manual BDF order, 1-5
//   initialstepfractionbdf-<order>    manual startup initial step fraction
//   initialstepgrowthratebdf-<order>  manual startup initial step growth rate
//   tstepsstore           "Store every Nth step"
//   tout                  tsteps | tstepsclosest | tlist
//   starttime, tlist, rtol, maxorder, keeplog
//   sol1/t1/aDef logsampling = 0 prints every accepted step; logging only
//
// comsolcompile emits a single class file, so the case table is a string list
// and all case state lives in static fields.
public class GeSteppingSweep {
    static final String[] SPECS = {
        // Maximum step constraint: automatic on two spans, constant, expression.
        "auto_span10",
        "auto_span100 end=100",
        "maxstep_const_005 maxstep=const:0.05",
        "maxstep_expr_sin maxstep=expr:0.05*(1+0.5*sin(t))",
        // Initial step: auto_span10 is the automatic reference of this group.
        "initstep_002 initstep=specified:0.02",
        "initstep_01 initstep=specified:0.1",
        // Start time: separates 0.1% of the end time from 0.1% of the span.
        "t0_2_span10 t0=2 end=12",
        // Store every Nth step. store_every1 repeats auto_span10 on purpose:
        // two independent solves of one configuration are a determinism control.
        "store_every1",
        "store_every3 every=3",
        // Manual BDF startup at the default fraction and growth rate.
        "manual_order1 manual=1",
        "manual_order2 manual=2",
        "manual_order3 manual=3",
        "manual_order5 manual=5",
        "manual_order2_frac1e-4 manual=2 fraction=1e-4",
        "manual_order2_frac1e-2 manual=2 fraction=1e-2",
        // Growth rate 1.5 is the order-2 default, so that row is a control;
        // order 3 (default fraction 0.1) is where the growth rate is observable.
        "manual_order2_growth1.5 manual=2 growth=1.5",
        "manual_order2_growth3 manual=2 growth=3",
        "manual_order3_growth1.5 manual=3 growth=1.5",
        "manual_order3_growth3 manual=3 growth=3",
        // Loose tolerance: the error test stops binding, so the automatic
        // maximum step and the automatic initial step become observable.
        "rtol1_span10 rtol=1",
        "rtol1_span100 rtol=1 end=100",
        "rtol1_t0_2 rtol=1 t0=2 end=12",
        // Control: identical to auto_span10 but with COMSOL's default log
        // sampling, which prints only a subset of the accepted steps.
        "log_default_sampling log=default",
    };

    static String directory;
    static PrintWriter manifest;

    static String name;
    static double t0;
    static double end;
    static String mode;
    static int maxOrder;
    static double rtol;
    static String maxStepKind;
    static String maxStepValue;
    static String initialStepKind;
    static double initialStepValue;
    static String storeMode;
    static int storeEvery;
    static int manualOrder;
    static double fraction;
    static double growth;
    static boolean fullLog;

    static void reset() {
        name = "case";
        t0 = 0.0;
        end = 10.0;
        mode = "free";
        maxOrder = 5;
        rtol = 0.01;
        maxStepKind = "auto";
        maxStepValue = "";
        initialStepKind = "auto";
        initialStepValue = 0.0;
        storeMode = "tsteps";
        storeEvery = 1;
        manualOrder = 0;
        fraction = 0.0;
        growth = 0.0;
        fullLog = true;
    }

    // "name key=value ..."; every omitted key keeps the default above, so a row
    // shows only what differs from a plain transient solve.
    static void apply(String spec) {
        reset();
        String[] tokens = spec.trim().split("\\s+");
        name = tokens[0];
        for (int i = 1; i < tokens.length; ++i) {
            int split = tokens[i].indexOf('=');
            if (split < 0) throw new IllegalArgumentException("bad token: " + tokens[i]);
            String key = tokens[i].substring(0, split);
            String value = tokens[i].substring(split + 1);
            if (key.equals("t0")) {
                t0 = Double.parseDouble(value);
            } else if (key.equals("end")) {
                end = Double.parseDouble(value);
            } else if (key.equals("mode")) {
                mode = value;
            } else if (key.equals("order")) {
                maxOrder = Integer.parseInt(value);
            } else if (key.equals("rtol")) {
                rtol = Double.parseDouble(value);
            } else if (key.equals("maxstep")) {
                if (!value.startsWith("const:") && !value.startsWith("expr:")) {
                    throw new IllegalArgumentException("maxstep needs const:<v> or expr:<v>: " + name);
                }
                maxStepKind = value.substring(0, value.indexOf(':'));
                maxStepValue = value.substring(value.indexOf(':') + 1);
            } else if (key.equals("initstep")) {
                if (!value.startsWith("specified:")) {
                    throw new IllegalArgumentException("initstep needs specified:<v>: " + name);
                }
                initialStepKind = "specified";
                initialStepValue = Double.parseDouble(value.substring(value.indexOf(':') + 1));
            } else if (key.equals("every")) {
                storeEvery = Integer.parseInt(value);
            } else if (key.equals("manual")) {
                mode = "manual";
                manualOrder = Integer.parseInt(value);
            } else if (key.equals("fraction")) {
                fraction = Double.parseDouble(value);
            } else if (key.equals("growth")) {
                growth = Double.parseDouble(value);
            } else if (key.equals("log")) {
                fullLog = value.equals("full");
            } else {
                throw new IllegalArgumentException("unknown key: " + key);
            }
        }
        if (mode.equals("manual") != (manualOrder > 0)) {
            throw new IllegalArgumentException("manual needs manual=<order>: " + name);
        }
        if (maxStepKind.equals("auto") != maxStepValue.isEmpty()) {
            throw new IllegalArgumentException("maxstep needs a kind and a value: " + name);
        }
        if (manualOrder > 0 && maxOrder < manualOrder) {
            throw new IllegalArgumentException("maxorder must cover bdforder: " + name);
        }
    }

    static String num(double value) {
        if (value == Math.rint(value)) return Long.toString((long) value);
        return Double.toString(value);
    }

    static String csv(String value) {
        return "\"" + value.replace("\"", "\"\"") + "\"";
    }

    static Model build() {
        double step = (end - t0) / 10.0;
        String times = "range(" + num(t0) + "," + num(step) + "," + num(end) + ")";
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

    static void run() throws Exception {
        Model m = build();
        SolverFeature f = m.sol("sol1").feature("t1");
        f.set("tout", storeMode);
        f.set("tstepsstore", Integer.toString(storeEvery));
        f.set("keeplog", "on");
        f.set("maxorder", Integer.toString(maxOrder));
        f.set("rtol", Double.toString(rtol));
        f.set("tstepsbdf", mode);
        if (t0 != 0.0) f.set("starttime", num(t0));
        if (!maxStepKind.equals("auto")) {
            f.set("maxstepconstraintbdf", maxStepKind);
            if (maxStepKind.equals("const")) {
                f.set("maxstepbdf", maxStepValue);
            } else {
                f.set("maxstepexpressionbdf", maxStepValue);
            }
        }
        if (initialStepKind.equals("specified")) {
            f.set("initialstepbdfactive", "on");
            f.set("initialstepbdf", Double.toString(initialStepValue));
        }
        double usedFraction = 0.0;
        double usedGrowth = 0.0;
        if (manualOrder > 0) {
            String order = Integer.toString(manualOrder);
            f.set("timestepbdf", "0.1");
            f.set("bdforder", order);
            if (fraction != 0.0) {
                f.set("initialstepfractionbdf-" + order, Double.toString(fraction));
            }
            if (growth != 0.0) {
                f.set("initialstepgrowthratebdf-" + order, Double.toString(growth));
            }
            usedFraction = Double.parseDouble(f.getString("initialstepfractionbdf-" + order));
            usedGrowth = Double.parseDouble(f.getString("initialstepgrowthratebdf-" + order));
        }
        if (fullLog) f.feature("aDef").set("logsampling", "0");
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
        PrintWriter logFile = new PrintWriter(new File(directory, name + ".log"));
        logFile.print(f.getString("message"));
        logFile.close();
        // The trailing columns are read back from the solver, so a property that
        // does not exist or was rejected cannot pass silently.
        String[] row = {
            name,
            num(t0),
            num(end),
            mode,
            Integer.toString(maxOrder),
            Double.toString(rtol),
            maxStepKind,
            csv(maxStepValue),
            initialStepKind,
            Double.toString(initialStepValue),
            storeMode,
            Integer.toString(storeEvery),
            Integer.toString(manualOrder),
            Double.toString(usedFraction),
            Double.toString(usedGrowth),
            fullLog ? "full" : "default",
            f.getString("maxstepconstraintbdf"),
            f.getString("tout"),
            f.getString("initialstepbdfactive"),
            f.getString("tstepsbdf"),
            f.getString("starttime"),
        };
        manifest.println(String.join(",", row));
        manifest.flush();
        ModelUtil.clear();
        System.out.println("CASE_OK " + name + " samples=" + t.length
            + " first=" + t[0] + " last=" + t[t.length - 1]);
    }

    // Prints the stepping properties and their values once per sweep, so the
    // names the cases set can be checked against the solver instead of trusted.
    static void reportProperties() {
        Model m = build();
        SolverFeature f = m.sol("sol1").feature("t1");
        String[] properties = f.properties();
        Arrays.sort(properties);
        for (String property : properties) {
            if (!property.matches("(maxstep|initialstep|tsteps|timestep|bdf|starttime"
                + "|rtol|maxorder|minorder|tout|keeplog|tstepping).*")) continue;
            String value;
            try {
                value = f.getString(property);
            } catch (Exception e) {
                value = "<not a string>";
            }
            System.out.println("PROP " + property + " = " + value);
        }
        ModelUtil.clear();
    }

    public static void main(String[] args) throws Exception {
        directory = args != null && args.length > 1 ? args[1] : "stepping";
        new File(directory).mkdirs();
        reportProperties();
        manifest = new PrintWriter(new File(directory, "manifest.csv"));
        manifest.println("name,t0,end,mode,max_order,rtol,max_step_kind,"
            + "max_step_value,initial_step_kind,initial_step_value,store_mode,"
            + "store_every,manual_order,initial_fraction,growth_rate,log_sampling,"
            + "max_step_kind_readback,store_mode_readback,initial_step_active_readback,"
            + "step_mode_readback,start_time_readback");
        for (String spec : SPECS) {
            apply(spec);
            run();
        }
        manifest.close();
        System.out.println("ALL_CASES_OK");
    }
}
