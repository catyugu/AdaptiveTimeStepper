import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.*;

// COMSOL's "Consistent initialization" sweep: one model integrated with the
// setting on and off, so the library's switch has a recorded trajectory on each
// side.
//
// The property below was read back from SolverFeature.properties() of sol1/t1 on
// COMSOL 6.2.0.290 instead of guessed:
//   consistent = bweuler | off    "Consistent initialization"
// `bweuler` is the default; the domain printed by main() shows which values the
// solver accepts. The other properties are the ones GeSteppingSweep reads back.
//
// Every case rebuilds the Fully Coupled node's Jacobian at every iteration
// (jtech=onevery, ntolfact=1e-8), so the record is the time integration's alone.
//
// Two models, each one Global Equation:
//   first   u'+u=0,           u(0)=1
//   second  qtt+0.1*qt+q=0,   q(0)=1, q'(0)=100
// and a case writes name.csv (time, value, derivative), name.log (the solver log)
// and one manifest.csv row.
public class GeConsistentInit {
    static final String[] SPECS = {
        // The setting itself, on the first-order model: the two backward-Euler
        // probes reconcile the initial value with the equation, and the
        // derivative they measure caps the first step.
        "on",
        "off cons=off",
        // The setting off, with the first step a case names instead of the
        // automatic one.
        "off_initial_step cons=off initstep=specified:0.02",
        // A manual run, whose startup curve replaces the automatic first step. Its
        // accumulated clock stops a few ulps short of the end, which is the margin
        // a run's end time has to allow.
        "off_manual cons=off manual=2 step=0.1",
        // The second-order form: its two initial values are reconciled by the
        // probes when the setting is on, and kept as they are when it is off.
        "second_on second=on",
        "second_off second=on cons=off",
    };

    static String directory;
    static PrintWriter manifest;

    static String name;
    static boolean secondOrder;
    static String consistent;
    static double tolerance;
    static int maxOrder;
    static String initialStepKind;
    static double initialStepValue;
    static int manualOrder;
    static double manualStep;
    static double t0;
    static double end;

    static void reset() {
        name = "case";
        secondOrder = false;
        consistent = "bweuler";
        tolerance = 0.01;
        maxOrder = 5;
        initialStepKind = "auto";
        initialStepValue = 0.0;
        manualOrder = 0;
        manualStep = 0.0;
        t0 = 0.0;
        end = 10.0;
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
            if (key.equals("second")) {
                secondOrder = value.equals("on");
            } else if (key.equals("cons")) {
                consistent = value;
            } else if (key.equals("rtol")) {
                tolerance = Double.parseDouble(value);
            } else if (key.equals("order")) {
                maxOrder = Integer.parseInt(value);
            } else if (key.equals("initstep")) {
                if (!value.startsWith("specified:")) {
                    throw new IllegalArgumentException("initstep needs specified:<v>: " + name);
                }
                initialStepKind = "specified";
                initialStepValue = Double.parseDouble(value.substring(value.indexOf(':') + 1));
            } else if (key.equals("manual")) {
                manualOrder = Integer.parseInt(value);
            } else if (key.equals("step")) {
                manualStep = Double.parseDouble(value);
            } else if (key.equals("t0")) {
                t0 = Double.parseDouble(value);
            } else if (key.equals("end")) {
                end = Double.parseDouble(value);
            } else {
                throw new IllegalArgumentException("unknown key: " + key);
            }
        }
        if (manualOrder > 0 && manualStep <= 0.0) manualStep = 0.1;
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
        Model m = ModelUtil.create("Model");
        m.component().create("comp1", true);
        m.component("comp1").geom().create("geom1", 1);
        m.component("comp1").geom("geom1").create("i1", "Interval");
        m.component("comp1").geom("geom1").run();
        m.component("comp1").mesh().create("mesh1");
        m.component("comp1").mesh("mesh1").run();
        m.component("comp1").physics().create("ge", "GlobalEquations", "geom1");
        if (secondOrder) {
            m.component("comp1").physics("ge").feature("ge1")
                .set("equation", "qtt+0.1*qt+q");
            m.component("comp1").physics("ge").feature("ge1").set("name", "q");
            m.component("comp1").physics("ge").feature("ge1")
                .set("initialValueU", "1");
            m.component("comp1").physics("ge").feature("ge1")
                .set("initialValueUt", "100");
        } else {
            m.component("comp1").physics("ge").feature("ge1").set("equation", "u+ut");
            m.component("comp1").physics("ge").feature("ge1").set("name", "u");
            m.component("comp1").physics("ge").feature("ge1").set("initialValueU", "1");
        }
        m.study().create("std1");
        m.study("std1").create("time", "Transient");
        m.study("std1").feature("time").set("tlist",
            "range(" + num(t0) + "," + num(step) + "," + num(end) + ")");
        m.study("std1").createAutoSequences("all");
        return m;
    }

    static void run() throws Exception {
        Model m = build();
        SolverFeature f = m.sol("sol1").feature("t1");
        f.set("tout", "tsteps");
        f.set("keeplog", "on");
        f.set("maxorder", Integer.toString(maxOrder));
        f.set("rtol", Double.toString(tolerance));
        f.set("consistent", consistent);
        if (t0 != 0.0) f.set("starttime", num(t0));
        if (initialStepKind.equals("specified")) {
            f.set("initialstepbdfactive", "on");
            f.set("initialstepbdf", Double.toString(initialStepValue));
        }
        if (manualOrder > 0) {
            f.set("tstepsbdf", "manual");
            f.set("timestepbdf", Double.toString(manualStep));
            f.set("bdforder", Integer.toString(manualOrder));
        }
        f.feature("fc1").set("jtech", "onevery");
        f.feature("fc1").set("ntolfact", "1e-8");
        f.feature("aDef").set("logsampling", "0");
        // Read back before the solve: a property that does not exist or was
        // rejected cannot then pass silently through the manifest.
        String consistentReadback = f.getString("consistent");
        String startReadback = f.getString("starttime");
        String initialStepReadback = f.getString("initialstepbdfactive");
        String stepModeReadback = f.getString("tstepsbdf");
        String jtechReadback = f.feature("fc1").getString("jtech");
        m.sol("sol1").runAll();
        double[] t = m.sol("sol1").getPVals();
        if (m.sol("sol1").getU(1).length != 1) {
            throw new IllegalStateException("expected one solved dof: " + name);
        }
        // The first sample, kept before the model is cleared: the record's starting
        // state and derivative are what the setting decides.
        double firstValue = m.sol("sol1").getU(1)[0];
        double firstDerivative = m.sol("sol1").getUDot(1)[0];
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
        String[] row = {
            name,
            secondOrder ? "second" : "first",
            consistent,
            Double.toString(tolerance),
            Integer.toString(maxOrder),
            initialStepKind,
            Double.toString(initialStepValue),
            Integer.toString(manualOrder),
            Double.toString(manualStep),
            num(t0),
            num(end),
            Integer.toString(t.length),
            consistentReadback,
            initialStepReadback,
            startReadback,
            stepModeReadback,
            jtechReadback,
        };
        manifest.println(String.join(",", row));
        manifest.flush();
        ModelUtil.clear();
        System.out.println("CASE_OK " + name + " consistent=" + consistentReadback
            + " samples=" + t.length + " first=(" + t[0] + ","
            + firstValue + "," + firstDerivative + ")");
    }

    // Prints the properties and the values the solver accepts, once per sweep, so
    // the name and the domain of the setting are measured rather than trusted.
    static void reportSetting() {
        Model probe = build();
        SolverFeature f = probe.sol("sol1").feature("t1");
        for (String property : f.properties()) {
            if (property.toLowerCase().contains("cons")) {
                System.out.println("PROP " + property + " = " + f.getString(property));
            }
        }
        ModelUtil.clear();
        for (String value : new String[] { "bweuler", "off", "none" }) {
            Model trial = build();
            SolverFeature g = trial.sol("sol1").feature("t1");
            String read;
            try {
                g.set("consistent", value);
                read = g.getString("consistent");
            } catch (Exception e) {
                read = "rejected";
            }
            System.out.println("DOMAIN consistent=" + value + " -> " + read);
            ModelUtil.clear();
        }
    }

    public static void main(String[] args) throws Exception {
        directory = args != null && args.length > 1 ? args[1] : "consistent";
        new File(directory).mkdirs();
        reset();
        reportSetting();
        manifest = new PrintWriter(new File(directory, "manifest.csv"));
        manifest.println("name,equation,consistent,tolerance,max_order,"
            + "initial_step_kind,initial_step_value,manual_order,manual_step,"
            + "t0,end,samples,consistent_readback,initial_step_active_readback,"
            + "start_time_readback,step_mode_readback,jtech");
        for (String spec : SPECS) {
            apply(spec);
            run();
        }
        manifest.close();
        System.out.println("ALL_CASES_OK");
    }
}
