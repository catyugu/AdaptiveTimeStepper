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
        // The controller: the automatic maximum step on two spans, the start
        // time that separates 0.1% of the end from 0.1% of the span, the
        // constant maximum step, and the expression one.
        "auto_span10",
        "auto_span100 end=100",
        "t0_2_span10 t0=2 end=12",
        "maxstep_const_005 maxstep=const:0.05",
        "maxstep_expr_sin maxstep=expr:0.05*(1+0.5*sin(t))",
        // The initial step: the automatic rule is auto_span10's, this is the
        // one the run names.
        "initstep_002 initstep=specified:0.02",
        // The tolerance: loose enough that the automatic maximum step and the
        // automatic initial step are what the steps show, and tight enough that
        // the order selection runs to five.
        "rtol1_span10 rtol=1",
        "rtol_1e5 rtol=1e-5",
        "max_order_3 order=3",
        // "Store every Nth step", whose modulus is anchored on the initial
        // state; the run rejects a step, so the phase is measured too.
        "store_every3 every=3",
        // The manual step: the per-order startup curve, and the fraction a run
        // may override.
        "manual_order3 manual=3",
        "manual_order2_frac1e-2 manual=2 fraction=1e-2",
        // The placement: a step at every output time, and a step inside every
        // output interval, on a uniform and on an uneven schedule.
        "strict mode=strict",
        "intermediate mode=intermediate",
        "irregular_strict mode=strict tlist=0 0.001 0.01 0.05 0.2 0.3 0.7 1 1.1 1.4 2 3 5 8 10",
        // The storage: the closest step, the values at the output times, and
        // the run that lands on the end instead of interpolating it.
        "store_closest tout=closest",
        "store_interpolate tout=tlist",
        "end_landing endinterp=off",
        // The solver's own Newton made exact, which the replay reaches with an
        // exactly solved step.
        "exact_default exact=on",
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
    static String tlist;
    static int storeEvery;
    static int manualOrder;
    static double fraction;
    static double growth;
    static boolean endInterpolation;
    static boolean exact;
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
        tlist = "";
        storeEvery = 1;
        manualOrder = 0;
        fraction = 0.0;
        growth = 0.0;
        endInterpolation = true;
        exact = false;
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
            } else if (key.equals("tout")) {
                if (value.equals("closest")) {
                    storeMode = "tstepsclosest";
                } else if (value.equals("tlist")) {
                    storeMode = "tlist";
                } else {
                    throw new IllegalArgumentException("unknown tout: " + name);
                }
            } else if (key.equals("tlist")) {
                // A schedule is a list, so this key takes the rest of the line.
                StringBuilder schedule = new StringBuilder(value);
                for (int j = i + 1; j < tokens.length; ++j) {
                    schedule.append(' ').append(tokens[j]);
                }
                tlist = schedule.toString();
                i = tokens.length;
            } else if (key.equals("endinterp")) {
                endInterpolation = !value.equals("off");
            } else if (key.equals("exact")) {
                exact = value.equals("on");
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
        // The study's schedule: ten equal intervals of the span, unless the case
        // names its own.
        String times = tlist.isEmpty()
            ? "range(" + num(t0) + "," + num(step) + "," + num(end) + ")"
            : tlist;
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
        f.set("endtimeinterpolation", endInterpolation ? "on" : "off");
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
        if (exact) {
            f.feature("fc1").set("jtech", "onevery");
            f.feature("fc1").set("ntolfact", "1e-8");
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
            endInterpolation ? "on" : "off",
            exact ? "on" : "off",
            csv(tlist),
            f.getString("maxstepconstraintbdf"),
            f.getString("tout"),
            f.getString("initialstepbdfactive"),
            f.getString("tstepsbdf"),
            f.getString("starttime"),
            f.getString("endtimeinterpolation"),
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
        reset();
        reportProperties();
        manifest = new PrintWriter(new File(directory, "manifest.csv"));
        manifest.println("name,t0,end,mode,max_order,rtol,max_step_kind,"
            + "max_step_value,initial_step_kind,initial_step_value,store_mode,"
            + "store_every,manual_order,initial_fraction,growth_rate,log_sampling,"
            + "end_interpolation,exact,tlist,"
            + "max_step_kind_readback,store_mode_readback,initial_step_active_readback,"
            + "step_mode_readback,start_time_readback,end_interpolation_readback");
        for (String spec : SPECS) {
            apply(spec);
            run();
        }
        manifest.close();
        System.out.println("ALL_CASES_OK");
    }
}
